"""
DeepST 多样本整合分析
支持批次效应校正和跨平台数据整合

适用于 AIDD 多组学数据分析
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
from typing import List, Dict
from tqdm import tqdm

from deepst_core import (
    DeviceManager,
    SpatialGraphConstructor,
    DeepSTModel,
    DeepSTTrainer,
    SpatialDomainIdentifier
)


class MultiSampleIntegrator:
    """多样本整合器 - 批次效应校正"""
    
    def __init__(
        self,
        n_domains: int = 7,
        hidden_dim: int = 512,
        latent_dim: int = 128,
        use_gpu: bool = True,
        output_dir: str = './integration_results'
    ):
        self.n_domains = n_domains
        self.hidden_dim = hidden_dim
        self.latent_dim = latent_dim
        self.output_dir = output_dir
        
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(f"{output_dir}/figures", exist_ok=True)
        
        self.device_manager = DeviceManager(prefer_gpu=use_gpu)
        self.graph_constructor = SpatialGraphConstructor(n_neighbors=6)
        self.domain_identifier = SpatialDomainIdentifier(n_domains=n_domains)
        
        self.sample_adatas = []
        self.combined_adata = None
        self.model = None
    
    def load_samples(
        self,
        data_paths: List[str],
        sample_names: List[str],
        platform: str = 'visium'
    ):
        """
        加载多个样本
        
        参数:
            data_paths: 数据文件路径列表
            sample_names: 样本名称列表
            platform: 平台类型
        """
        print(f"\n{'='*60}")
        print(f"加载 {len(data_paths)} 个样本")
        print(f"{'='*60}")
        
        for path, name in zip(data_paths, sample_names):
            print(f"\n加载样本: {name}")
            
            if path.endswith('.h5ad'):
                adata = sc.read_h5ad(path)
            elif platform.lower() == 'visium':
                adata = sc.read_visium(path)
            else:
                raise ValueError(f"不支持的平台: {platform}")
            
            # 添加批次标签
            adata.obs['batch'] = name
            adata.obs['batch_id'] = len(self.sample_adatas)
            
            # 基本预处理
            sc.pp.filter_genes(adata, min_cells=10)
            sc.pp.filter_cells(adata, min_genes=200)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            self.sample_adatas.append(adata)
            
            print(f"  ✓ {adata.n_obs} spots, {adata.n_vars} genes")
        
        print(f"\n✓ 总计加载 {len(self.sample_adatas)} 个样本")
    
    def prepare_integration(
        self,
        n_top_genes: int = 3000,
        pca_components: int = 200
    ):
        """
        准备整合数据
        
        参数:
            n_top_genes: 高变基因数量
            pca_components: PCA维度
        """
        print(f"\n{'='*60}")
        print("准备整合数据")
        print(f"{'='*60}")
        
        # 找共同的高变基因
        print("识别共同高变基因...")
        for adata in self.sample_adatas:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                flavor='seurat_v3'
            )
        
        # 取交集
        hvg_lists = [set(adata.var_names[adata.var.highly_variable]) 
                     for adata in self.sample_adatas]
        common_hvgs = list(set.intersection(*hvg_lists))
        
        print(f"✓ 找到 {len(common_hvgs)} 个共同高变基因")
        
        # 合并样本
        print("合并样本...")
        combined = ad.concat(
            self.sample_adatas,
            join='inner',
            label='batch_name',
            keys=[adata.obs['batch'][0] for adata in self.sample_adatas]
        )
        
        # 只保留共同高变基因
        combined = combined[:, common_hvgs].copy()
        
        # PCA降维
        print("PCA降维...")
        from sklearn.decomposition import PCA
        
        expression = combined.X
        if hasattr(expression, 'toarray'):
            expression = expression.toarray()
        
        pca = PCA(n_components=pca_components)
        reduced = pca.fit_transform(expression)
        combined.obsm['X_pca'] = reduced
        
        print(f"✓ PCA: {pca_components}维, 解释方差={pca.explained_variance_ratio_.sum():.2%}")
        
        self.combined_adata = combined
        
        # 构建图
        self._build_combined_graph()
    
    def _build_combined_graph(self):
        """构建整合后的空间图"""
        print("\n构建整合空间图...")
        
        all_edges = []
        all_weights = []
        offset = 0
        
        # 为每个批次单独构建图
        for batch_name in self.combined_adata.obs['batch_name'].unique():
            batch_mask = self.combined_adata.obs['batch_name'] == batch_name
            batch_coords = self.combined_adata[batch_mask].obsm['spatial']
            
            edge_index, edge_weight = self.graph_constructor.build_graph(batch_coords)
            
            # 调整索引
            edge_index = edge_index + offset
            
            all_edges.append(edge_index)
            all_weights.append(edge_weight)
            
            offset += batch_mask.sum()
        
        # 合并所有边
        combined_edges = np.concatenate(all_edges, axis=1)
        combined_weights = np.concatenate(all_weights)
        
        # 创建 PyG Data
        self.graph_data = Data(
            x=torch.FloatTensor(self.combined_adata.obsm['X_pca']),
            edge_index=torch.LongTensor(combined_edges),
            edge_attr=torch.FloatTensor(combined_weights),
            pos=torch.FloatTensor(self.combined_adata.obsm['spatial'])
        )
        
        print(f"✓ 整合图: {combined_edges.shape[1]} 条边")
    
    def build_model(self):
        """构建带域对抗的模型"""
        print(f"\n{'='*60}")
        print("构建整合模型")
        print(f"{'='*60}")
        
        n_genes = self.graph_data.x.shape[1]
        n_spots = self.graph_data.x.shape[0]
        n_batches = len(self.sample_adatas)
        
        self.model = DeepSTModel(
            n_genes=n_genes,
            n_spots=n_spots,
            hidden_dim=self.hidden_dim,
            latent_dim=self.latent_dim,
            n_domains=n_batches,  # 域数量 = 批次数量
            use_batch_correction=True
        )
        
        n_params = sum(p.numel() for p in self.model.parameters())
        print(f"✓ 整合模型构建完成")
        print(f"  - 参数量: {n_params:,}")
        print(f"  - 批次数: {n_batches}")
        print(f"  - 使用域对抗: 是")
        
        self.trainer = DeepSTTrainer(
            model=self.model,
            device_manager=self.device_manager,
            learning_rate=1e-3
        )
    
    def train(
        self,
        n_epochs: int = 500,
        lambda_domain: float = 0.1
    ):
        """
        训练整合模型
        
        参数:
            n_epochs: 训练轮数
            lambda_domain: 域对抗损失权重
        """
        print(f"\n{'='*60}")
        print(f"开始整合训练 ({n_epochs} epochs)")
        print(f"{'='*60}")
        
        # 准备批次标签
        batch_labels = torch.LongTensor(
            self.combined_adata.obs['batch_id'].values
        )
        
        history = {
            'total_loss': [],
            'recon_loss': [],
            'domain_loss': []
        }
        
        progress_bar = tqdm(range(n_epochs), desc="训练进度")
        
        for epoch in progress_bar:
            # 域对抗的 alpha 参数 (逐渐增加)
            alpha = min(1.0, epoch / (n_epochs * 0.5))
            
            losses = self.trainer.train_epoch(
                self.graph_data,
                domain_labels=batch_labels,
                alpha=alpha,
                lambda_recon=1.0,
                lambda_domain=lambda_domain
            )
            
            history['total_loss'].append(losses['total_loss'])
            history['recon_loss'].append(losses['recon_loss'])
            history['domain_loss'].append(losses['domain_loss'])
            
            progress_bar.set_postfix({
                'Loss': f"{losses['total_loss']:.4f}",
                'Recon': f"{losses['recon_loss']:.4f}",
                'Domain': f"{losses['domain_loss']:.4f}"
            })
            
            self.trainer.scheduler.step()
        
        print("\n✓ 整合训练完成")
        self._plot_integration_history(history)
        
        return history
    
    def get_integrated_embeddings(self):
        """获取整合后的嵌入"""
        print("\n获取整合嵌入...")
        
        embeddings = self.trainer.get_embeddings(self.graph_data)
        self.combined_adata.obsm['integrated_embed'] = embeddings
        
        np.save(
            f"{self.output_dir}/integrated_embeddings.npy",
            embeddings
        )
        
        print(f"✓ 整合嵌入保存")
        return embeddings
    
    def identify_integrated_domains(self):
        """识别整合后的空间域"""
        print(f"\n{'='*60}")
        print("识别整合空间域")
        print(f"{'='*60}")
        
        embeddings = self.combined_adata.obsm['integrated_embed']
        spatial_coords = self.combined_adata.obsm['spatial']
        
        domain_labels = self.domain_identifier.identify_domains(
            embeddings,
            spatial_coords
        )
        
        self.combined_adata.obs['integrated_domain'] = domain_labels.astype(str)
        
        print(f"✓ 识别到 {len(np.unique(domain_labels))} 个整合空间域")
        
        # 统计每个批次的域分布
        print("\n批次-域分布:")
        for batch in self.combined_adata.obs['batch_name'].unique():
            batch_mask = self.combined_adata.obs['batch_name'] == batch
            batch_domains = domain_labels[batch_mask]
            unique, counts = np.unique(batch_domains, return_counts=True)
            print(f"\n{batch}:")
            for d, c in zip(unique, counts):
                print(f"  域 {d}: {c} spots")
        
        return domain_labels
    
    def visualize_integration(self):
        """可视化整合结果"""
        print(f"\n{'='*60}")
        print("生成整合可视化")
        print(f"{'='*60}")
        
        # 1. UMAP - 批次和域
        sc.pp.neighbors(self.combined_adata, use_rep='integrated_embed', n_neighbors=15)
        sc.tl.umap(self.combined_adata)
        
        fig, axes = plt.subplots(1, 3, figsize=(20, 5))
        
        sc.pl.umap(
            self.combined_adata,
            color='batch_name',
            title='UMAP - Batches',
            ax=axes[0],
            show=False,
            frameon=False
        )
        
        sc.pl.umap(
            self.combined_adata,
            color='integrated_domain',
            title='UMAP - Integrated Domains',
            ax=axes[1],
            show=False,
            frameon=False
        )
        
        # 批次混合度
        from sklearn.metrics import silhouette_score
        
        batch_ids = self.combined_adata.obs['batch_id'].values
        embeddings = self.combined_adata.obsm['integrated_embed']
        
        sil_score = silhouette_score(embeddings, batch_ids)
        
        axes[2].text(
            0.5, 0.5,
            f"批次混合度\n(Silhouette Score)\n\n{sil_score:.4f}\n\n"
            f"(越接近0越好)",
            ha='center',
            va='center',
            fontsize=16,
            transform=axes[2].transAxes
        )
        axes[2].axis('off')
        
        plt.tight_layout()
        plt.savefig(
            f"{self.output_dir}/figures/integration_umap.png",
            dpi=300,
            bbox_inches='tight'
        )
        print("✓ 保存: integration_umap.png")
        plt.close()
        
        # 2. 每个样本的空间域
        n_samples = len(self.sample_adatas)
        fig, axes = plt.subplots(1, n_samples, figsize=(6*n_samples, 6))
        
        if n_samples == 1:
            axes = [axes]
        
        for i, batch_name in enumerate(self.combined_adata.obs['batch_name'].unique()):
            batch_data = self.combined_adata[
                self.combined_adata.obs['batch_name'] == batch_name
            ]
            
            sc.pl.spatial(
                batch_data,
                color='integrated_domain',
                spot_size=150,
                title=f'{batch_name} - Integrated Domains',
                ax=axes[i],
                show=False,
                frameon=False
            )
        
        plt.tight_layout()
        plt.savefig(
            f"{self.output_dir}/figures/spatial_domains_all_samples.png",
            dpi=300,
            bbox_inches='tight'
        )
        print("✓ 保存: spatial_domains_all_samples.png")
        plt.close()
    
    def save_integration_results(self):
        """保存整合结果"""
        print(f"\n{'='*60}")
        print("保存整合结果")
        print(f"{'='*60}")
        
        # 保存整合后的 AnnData
        self.combined_adata.write_h5ad(
            f"{self.output_dir}/integrated_results.h5ad"
        )
        print(f"✓ 整合数据: integrated_results.h5ad")
        
        # 保存域标签
        domain_df = pd.DataFrame({
            'spot_id': self.combined_adata.obs_names,
            'batch': self.combined_adata.obs['batch_name'],
            'domain': self.combined_adata.obs['integrated_domain']
        })
        domain_df.to_csv(
            f"{self.output_dir}/integrated_domains.csv",
            index=False
        )
        print(f"✓ 域标签: integrated_domains.csv")
        
        # 保存模型
        torch.save(
            self.model.state_dict(),
            f"{self.output_dir}/integration_model.pth"
        )
        print(f"✓ 模型: integration_model.pth")
    
    def _plot_integration_history(self, history):
        """绘制整合训练历史"""
        fig, axes = plt.subplots(1, 3, figsize=(18, 4))
        
        axes[0].plot(history['total_loss'])
        axes[0].set_title('Total Loss')
        axes[0].set_xlabel('Epoch')
        axes[0].set_ylabel('Loss')
        axes[0].grid(True, alpha=0.3)
        
        axes[1].plot(history['recon_loss'], color='orange')
        axes[1].set_title('Reconstruction Loss')
        axes[1].set_xlabel('Epoch')
        axes[1].set_ylabel('Loss')
        axes[1].grid(True, alpha=0.3)
        
        axes[2].plot(history['domain_loss'], color='green')
        axes[2].set_title('Domain Adversarial Loss')
        axes[2].set_xlabel('Epoch')
        axes[2].set_ylabel('Loss')
        axes[2].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(
            f"{self.output_dir}/figures/integration_training.png",
            dpi=300,
            bbox_inches='tight'
        )
        plt.close()


def example_integration():
    """示例：多样本整合分析"""
    
    print("\n" + "="*60)
    print("DeepST 多样本整合分析示例")
    print("="*60 + "\n")
    
    # 配置
    sample_paths = [
        './data/sample1.h5ad',
        './data/sample2.h5ad',
        './data/sample3.h5ad'
    ]
    
    sample_names = ['Sample_1', 'Sample_2', 'Sample_3']
    
    # 创建整合器
    integrator = MultiSampleIntegrator(
        n_domains=7,
        hidden_dim=512,
        latent_dim=128,
        use_gpu=True,
        output_dir='./integration_results'
    )
    
    # 1. 加载样本
    # integrator.load_samples(sample_paths, sample_names)
    
    # 示例：生成模拟数据
    print("生成模拟多样本数据...")
    from sklearn.datasets import make_blobs
    import anndata as ad
    
    for i in range(3):
        spatial_coords, _ = make_blobs(
            n_samples=500,
            n_features=2,
            centers=7,
            cluster_std=5.0,
            random_state=42 + i
        )
        
        gene_expression = np.random.negative_binomial(5, 0.3, (500, 2000))
        
        adata = ad.AnnData(X=gene_expression)
        adata.obsm['spatial'] = spatial_coords
        adata.obs['batch'] = f'Sample_{i+1}'
        adata.obs['batch_id'] = i
        
        integrator.sample_adatas.append(adata)
    
    # 2. 准备整合
    integrator.prepare_integration()
    
    # 3. 构建模型
    integrator.build_model()
    
    # 4. 训练
    integrator.train(n_epochs=300)
    
    # 5. 获取嵌入
    integrator.get_integrated_embeddings()
    
    # 6. 识别域
    integrator.identify_integrated_domains()
    
    # 7. 可视化
    integrator.visualize_integration()
    
    # 8. 保存
    integrator.save_integration_results()
    
    print(f"\n{'='*60}")
    print("整合分析完成!")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    example_integration()
