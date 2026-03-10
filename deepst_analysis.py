"""
DeepST 完整分析流程
单样本空间域识别

适用于 macOS Intel, 自动 CPU/GPU 兼容
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
import torch
from torch_geometric.data import Data
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from tqdm import tqdm
import warnings

# 导入自定义模块
from deepst_core import (
    DeviceManager,
    SpatialGraphConstructor,
    DeepSTModel,
    DeepSTTrainer,
    SpatialDomainIdentifier,
    preprocess_data
)

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0

# 设置随机种子
def set_seed(seed=42):
    """设置随机种子以保证可重复性"""
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

set_seed(42)


class DeepSTAnalyzer:
    """DeepST 分析器 - 完整流程封装"""
    
    def __init__(
        self,
        n_domains: int = 7,
        hidden_dim: int = 512,
        latent_dim: int = 128,
        use_gpu: bool = True,
        output_dir: str = './results'
    ):
        """
        参数:
            n_domains: 空间域数量
            hidden_dim: 隐藏层维度
            latent_dim: 潜在表示维度
            use_gpu: 是否使用GPU (如果可用)
            output_dir: 输出目录
        """
        self.n_domains = n_domains
        self.hidden_dim = hidden_dim
        self.latent_dim = latent_dim
        self.output_dir = output_dir
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(f"{output_dir}/figures", exist_ok=True)
        os.makedirs(f"{output_dir}/embeddings", exist_ok=True)
        
        # 设备管理
        self.device_manager = DeviceManager(prefer_gpu=use_gpu)
        
        # 初始化组件
        self.graph_constructor = SpatialGraphConstructor(
            n_neighbors=6,
            distance_type='KDTree'
        )
        
        self.domain_identifier = SpatialDomainIdentifier(
            n_domains=n_domains,
            use_spatial_refine=True,
            refine_iterations=10
        )
        
        self.adata = None
        self.model = None
        self.trainer = None
    
    def load_data(
        self,
        file_path: str,
        platform: str = 'visium'
    ) -> ad.AnnData:
        """
        加载空间转录组学数据
        
        参数:
            file_path: 数据文件路径 (.h5ad 或目录)
            platform: 平台类型 ('visium', 'slideseq', 'stereo', etc.)
            
        返回:
            adata: AnnData对象
        """
        print(f"\n{'='*60}")
        print(f"加载数据: {file_path}")
        print(f"{'='*60}")
        
        if file_path.endswith('.h5ad'):
            adata = sc.read_h5ad(file_path)
        elif platform.lower() == 'visium':
            adata = sc.read_visium(file_path)
        else:
            raise ValueError(f"不支持的平台: {platform}")
        
        # 基本质控
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.filter_cells(adata, min_genes=200)
        
        # 归一化
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # 识别高变基因
        sc.pp.highly_variable_genes(adata, n_top_genes=3000)
        
        print(f"✓ 数据加载完成:")
        print(f"  - Spots: {adata.n_obs}")
        print(f"  - Genes: {adata.n_vars}")
        print(f"  - 高变基因: {adata.var['highly_variable'].sum()}")
        
        self.adata = adata
        return adata
    
    def prepare_data(
        self,
        use_highly_variable: bool = True,
        pca_components: int = 200,
        spatial_smooth: bool = True,
        smooth_alpha: float = 0.3
    ):
        """
        准备训练数据
        
        参数:
            use_highly_variable: 是否只使用高变基因
            pca_components: PCA降维维数
            spatial_smooth: 是否进行空间平滑
            smooth_alpha: 平滑强度
        """
        print(f"\n{'='*60}")
        print("数据预处理")
        print(f"{'='*60}")
        
        # 提取表达矩阵
        if use_highly_variable and 'highly_variable' in self.adata.var:
            gene_expression = self.adata[:, self.adata.var['highly_variable']].X
        else:
            gene_expression = self.adata.X
        
        # 转为稠密矩阵
        if hasattr(gene_expression, 'toarray'):
            gene_expression = gene_expression.toarray()
        
        # 空间坐标
        spatial_coords = self.adata.obsm['spatial']
        
        # 空间平滑
        if spatial_smooth:
            print("执行空间平滑...")
            gene_expression = self.graph_constructor.spatial_smoothing(
                gene_expression,
                spatial_coords,
                alpha=smooth_alpha
            )
        
        # PCA降维
        print("PCA降维...")
        from sklearn.decomposition import PCA
        pca = PCA(n_components=min(pca_components, gene_expression.shape[1]))
        reduced_expression = pca.fit_transform(gene_expression)
        
        variance_explained = pca.explained_variance_ratio_.sum()
        print(f"✓ PCA完成: {pca_components}维, 解释方差={variance_explained:.2%}")
        
        # 构建空间图
        print("构建空间邻接图...")
        edge_index, edge_weight = self.graph_constructor.build_graph(spatial_coords)
        
        # 创建 PyTorch Geometric Data
        self.graph_data = Data(
            x=torch.FloatTensor(reduced_expression),
            edge_index=torch.LongTensor(edge_index),
            edge_attr=torch.FloatTensor(edge_weight),
            pos=torch.FloatTensor(spatial_coords)
        )
        
        # 保存到 adata
        self.adata.obsm['X_pca'] = reduced_expression
        
        print(f"✓ 数据准备完成")
    
    def build_model(self):
        """构建 DeepST 模型"""
        print(f"\n{'='*60}")
        print("构建模型")
        print(f"{'='*60}")
        
        n_genes = self.graph_data.x.shape[1]
        n_spots = self.graph_data.x.shape[0]
        
        self.model = DeepSTModel(
            n_genes=n_genes,
            n_spots=n_spots,
            hidden_dim=self.hidden_dim,
            latent_dim=self.latent_dim,
            use_batch_correction=False
        )
        
        # 计算参数量
        n_params = sum(p.numel() for p in self.model.parameters())
        print(f"✓ 模型构建完成")
        print(f"  - 参数量: {n_params:,}")
        print(f"  - 输入维度: {n_genes}")
        print(f"  - 潜在维度: {self.latent_dim}")
        
        # 创建训练器
        self.trainer = DeepSTTrainer(
            model=self.model,
            device_manager=self.device_manager,
            learning_rate=1e-3
        )
    
    def train(
        self,
        n_epochs: int = 500,
        print_every: int = 50
    ):
        """
        训练模型
        
        参数:
            n_epochs: 训练轮数
            print_every: 打印间隔
        """
        print(f"\n{'='*60}")
        print(f"开始训练 ({n_epochs} epochs)")
        print(f"{'='*60}")
        
        history = {
            'total_loss': [],
            'recon_loss': []
        }
        
        progress_bar = tqdm(range(n_epochs), desc="训练进度")
        
        for epoch in progress_bar:
            # 训练一个epoch
            losses = self.trainer.train_epoch(
                self.graph_data,
                lambda_recon=1.0
            )
            
            history['total_loss'].append(losses['total_loss'])
            history['recon_loss'].append(losses['recon_loss'])
            
            # 更新进度条
            progress_bar.set_postfix({
                'Loss': f"{losses['total_loss']:.4f}",
                'Recon': f"{losses['recon_loss']:.4f}"
            })
            
            # 定期打印
            if (epoch + 1) % print_every == 0:
                print(f"\nEpoch {epoch+1}/{n_epochs} - "
                      f"Loss: {losses['total_loss']:.4f}, "
                      f"Recon: {losses['recon_loss']:.4f}")
            
            # 学习率调度
            self.trainer.scheduler.step()
        
        print("\n✓ 训练完成")
        
        # 绘制训练曲线
        self._plot_training_history(history)
        
        return history
    
    def get_embeddings(self) -> np.ndarray:
        """获取嵌入表示"""
        print("\n获取嵌入表示...")
        embeddings = self.trainer.get_embeddings(self.graph_data)
        self.adata.obsm['DeepST_embed'] = embeddings
        
        # 保存
        np.save(f"{self.output_dir}/embeddings/deepst_embeddings.npy", embeddings)
        print(f"✓ 嵌入保存至: {self.output_dir}/embeddings/")
        
        return embeddings
    
    def identify_domains(self) -> np.ndarray:
        """识别空间域"""
        print(f"\n{'='*60}")
        print("识别空间域")
        print(f"{'='*60}")
        
        embeddings = self.adata.obsm['DeepST_embed']
        spatial_coords = self.adata.obsm['spatial']
        
        domain_labels = self.domain_identifier.identify_domains(
            embeddings,
            spatial_coords
        )
        
        self.adata.obs['DeepST_domain'] = domain_labels.astype(str)
        
        print(f"✓ 识别到 {len(np.unique(domain_labels))} 个空间域")
        print("域分布:")
        for i in range(self.n_domains):
            count = (domain_labels == i).sum()
            print(f"  域 {i}: {count} spots ({count/len(domain_labels)*100:.1f}%)")
        
        return domain_labels
    
    def visualize_results(self):
        """可视化结果"""
        print(f"\n{'='*60}")
        print("生成可视化")
        print(f"{'='*60}")
        
        # 1. 空间域可视化
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        
        sc.pl.spatial(
            self.adata,
            color=['DeepST_domain'],
            spot_size=150,
            frameon=False,
            title='DeepST Spatial Domains',
            ax=ax,
            show=False
        )
        
        plt.tight_layout()
        plt.savefig(
            f"{self.output_dir}/figures/spatial_domains.png",
            dpi=300,
            bbox_inches='tight'
        )
        print(f"✓ 保存: spatial_domains.png")
        plt.close()
        
        # 2. UMAP可视化
        print("计算UMAP...")
        sc.pp.neighbors(self.adata, use_rep='DeepST_embed', n_neighbors=15)
        sc.tl.umap(self.adata)
        
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        sc.pl.umap(
            self.adata,
            color='DeepST_domain',
            title='UMAP - DeepST Domains',
            ax=axes[0],
            show=False,
            frameon=False
        )
        
        # 如果有真实标签
        if 'ground_truth' in self.adata.obs:
            sc.pl.umap(
                self.adata,
                color='ground_truth',
                title='UMAP - Ground Truth',
                ax=axes[1],
                show=False,
                frameon=False
            )
        
        plt.tight_layout()
        plt.savefig(
            f"{self.output_dir}/figures/umap_visualization.png",
            dpi=300,
            bbox_inches='tight'
        )
        print(f"✓ 保存: umap_visualization.png")
        plt.close()
        
        # 3. 嵌入热图
        self._plot_embedding_heatmap()
    
    def evaluate_performance(
        self,
        ground_truth_key: str = 'ground_truth'
    ) -> Dict:
        """
        评估性能 (如果有真实标签)
        
        参数:
            ground_truth_key: 真实标签的key
            
        返回:
            metrics: 评估指标字典
        """
        if ground_truth_key not in self.adata.obs:
            print(f"警告: 未找到真实标签 '{ground_truth_key}'")
            return {}
        
        print(f"\n{'='*60}")
        print("性能评估")
        print(f"{'='*60}")
        
        pred_labels = self.adata.obs['DeepST_domain'].astype(int)
        true_labels = self.adata.obs[ground_truth_key].astype(int)
        
        # 计算指标
        ari = adjusted_rand_score(true_labels, pred_labels)
        nmi = normalized_mutual_info_score(true_labels, pred_labels)
        
        metrics = {
            'ARI': ari,
            'NMI': nmi
        }
        
        print(f"✓ ARI: {ari:.4f}")
        print(f"✓ NMI: {nmi:.4f}")
        
        # 保存指标
        metrics_df = pd.DataFrame([metrics])
        metrics_df.to_csv(f"{self.output_dir}/metrics.csv", index=False)
        print(f"✓ 指标保存至: {self.output_dir}/metrics.csv")
        
        return metrics
    
    def save_results(self):
        """保存结果"""
        print(f"\n{'='*60}")
        print("保存结果")
        print(f"{'='*60}")
        
        # 保存 AnnData
        output_path = f"{self.output_dir}/deepst_results.h5ad"
        self.adata.write_h5ad(output_path)
        print(f"✓ AnnData 保存至: {output_path}")
        
        # 保存域标签
        domain_df = pd.DataFrame({
            'spot_id': self.adata.obs_names,
            'domain': self.adata.obs['DeepST_domain']
        })
        domain_df.to_csv(f"{self.output_dir}/domain_labels.csv", index=False)
        print(f"✓ 域标签保存至: {self.output_dir}/domain_labels.csv")
        
        # 保存模型
        torch.save(
            self.model.state_dict(),
            f"{self.output_dir}/deepst_model.pth"
        )
        print(f"✓ 模型保存至: {self.output_dir}/deepst_model.pth")
    
    def _plot_training_history(self, history):
        """绘制训练历史"""
        fig, axes = plt.subplots(1, 2, figsize=(14, 4))
        
        # 总损失
        axes[0].plot(history['total_loss'], label='Total Loss')
        axes[0].set_xlabel('Epoch')
        axes[0].set_ylabel('Loss')
        axes[0].set_title('Training Loss')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # 重构损失
        axes[1].plot(history['recon_loss'], label='Reconstruction Loss', color='orange')
        axes[1].set_xlabel('Epoch')
        axes[1].set_ylabel('Loss')
        axes[1].set_title('Reconstruction Loss')
        axes[1].legend()
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(
            f"{self.output_dir}/figures/training_history.png",
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()
    
    def _plot_embedding_heatmap(self):
        """绘制嵌入热图"""
        embeddings = self.adata.obsm['DeepST_embed']
        domain_labels = self.adata.obs['DeepST_domain'].astype(int).values
        
        # 按域排序
        sorted_indices = np.argsort(domain_labels)
        sorted_embeddings = embeddings[sorted_indices]
        sorted_labels = domain_labels[sorted_indices]
        
        # 绘制热图
        fig, ax = plt.subplots(figsize=(12, 8))
        
        im = ax.imshow(
            sorted_embeddings.T,
            aspect='auto',
            cmap='viridis',
            interpolation='nearest'
        )
        
        # 添加域边界
        unique_domains = np.unique(sorted_labels)
        boundaries = []
        for domain in unique_domains[:-1]:
            boundary = np.where(sorted_labels == domain)[0][-1]
            boundaries.append(boundary)
            ax.axvline(boundary, color='red', linestyle='--', linewidth=1, alpha=0.5)
        
        ax.set_xlabel('Spots (sorted by domain)')
        ax.set_ylabel('Embedding Dimensions')
        ax.set_title('DeepST Embeddings Heatmap')
        
        plt.colorbar(im, ax=ax, label='Embedding Value')
        plt.tight_layout()
        plt.savefig(
            f"{self.output_dir}/figures/embedding_heatmap.png",
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()


def main():
    """主函数 - 示例流程"""
    
    print("\n" + "="*60)
    print("DeepST 空间转录组学分析")
    print("macOS Intel 兼容版本")
    print("="*60 + "\n")
    
    # 配置参数
    config = {
        'data_path': './data/sample.h5ad',  # 修改为您的数据路径
        'n_domains': 7,
        'hidden_dim': 512,
        'latent_dim': 128,
        'n_epochs': 500,
        'use_gpu': True,
        'output_dir': './results'
    }
    
    # 创建分析器
    analyzer = DeepSTAnalyzer(
        n_domains=config['n_domains'],
        hidden_dim=config['hidden_dim'],
        latent_dim=config['latent_dim'],
        use_gpu=config['use_gpu'],
        output_dir=config['output_dir']
    )
    
    # 1. 加载数据
    # analyzer.load_data(config['data_path'])
    
    # 如果没有数据，创建模拟数据用于测试
    print("生成模拟数据用于演示...")
    from sklearn.datasets import make_blobs
    
    # 模拟空间坐标和基因表达
    n_spots = 1000
    n_genes = 2000
    
    spatial_coords, true_labels = make_blobs(
        n_samples=n_spots,
        n_features=2,
        centers=7,
        cluster_std=5.0,
        random_state=42
    )
    
    gene_expression = np.random.negative_binomial(5, 0.3, (n_spots, n_genes))
    
    # 创建 AnnData
    adata = ad.AnnData(X=gene_expression)
    adata.obsm['spatial'] = spatial_coords
    adata.obs['ground_truth'] = true_labels.astype(str)
    
    analyzer.adata = adata
    
    # 2. 数据预处理
    analyzer.prepare_data(
        use_highly_variable=False,
        pca_components=200,
        spatial_smooth=True
    )
    
    # 3. 构建模型
    analyzer.build_model()
    
    # 4. 训练
    analyzer.train(n_epochs=config['n_epochs'])
    
    # 5. 获取嵌入
    analyzer.get_embeddings()
    
    # 6. 识别空间域
    analyzer.identify_domains()
    
    # 7. 可视化
    analyzer.visualize_results()
    
    # 8. 评估 (如果有真实标签)
    analyzer.evaluate_performance('ground_truth')
    
    # 9. 保存结果
    analyzer.save_results()
    
    print(f"\n{'='*60}")
    print("分析完成!")
    print(f"结果保存在: {config['output_dir']}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
