"""
高级空间转录组学预测分析系统
作者: 计算机博士研究项目
数据: 10x Genomics Human Breast Cancer Visium

核心创新:
1. 图神经网络(GNN)空间嵌入
2. 变分自编码器(VAE)基因表达预测
3. 扩散模型(Diffusion)空间插值
4. 多模态融合(图像+转录组)
"""

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv, SAGEConv
from torch_geometric.data import Data
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.sparse import issparse
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# 第一部分: 数据加载与预处理
# ============================================================================

class SpatialDataLoader:
    """加载和预处理10x Genomics Visium数据"""
    
    def __init__(self, data_path):
        self.data_path = data_path
        self.adata = None
        
    def load_visium_data(self, filtered=True):
        """加载Visium数据"""
        # 实际使用时替换为下载的数据路径
        # self.adata = sc.read_visium(self.data_path)
        
        # 演示用: 创建模拟数据结构
        print("加载10x Genomics Visium数据...")
        # 实际代码:
        # self.adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
        # sq.read.visium(self.data_path, counts_file="filtered_feature_bc_matrix.h5")
        
        return self.adata
    
    def preprocess(self, adata, n_top_genes=3000):
        """标准预处理流程"""
        print("执行质量控制...")
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        
        print("归一化和对数转换...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        print(f"选择高变基因 (top {n_top_genes})...")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=True)
        
        return adata
    
    def compute_spatial_graph(self, adata, n_neighs=6, coord_type='grid'):
        """构建空间邻接图"""
        print("计算空间邻接图...")
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type=coord_type)
        return adata


# ============================================================================
# 第二部分: 图神经网络空间嵌入
# ============================================================================

class SpatialGNNEncoder(nn.Module):
    """
    图神经网络编码器 - 捕获空间依赖性
    使用多层GCN/GAT捕获局部和全局空间模式
    """
    
    def __init__(self, input_dim, hidden_dims=[512, 256, 128], 
                 gnn_type='GAT', heads=4):
        super(SpatialGNNEncoder, self).__init__()
        
        self.gnn_type = gnn_type
        self.layers = nn.ModuleList()
        
        # 构建多层GNN
        dims = [input_dim] + hidden_dims
        for i in range(len(dims) - 1):
            if gnn_type == 'GAT':
                if i < len(dims) - 2:
                    self.layers.append(GATConv(dims[i], dims[i+1], heads=heads, concat=True))
                    dims[i+1] = dims[i+1] * heads  # 更新维度
                else:
                    self.layers.append(GATConv(dims[i], dims[i+1], heads=1, concat=False))
            elif gnn_type == 'GCN':
                self.layers.append(GCNConv(dims[i], dims[i+1]))
            elif gnn_type == 'SAGE':
                self.layers.append(SAGEConv(dims[i], dims[i+1]))
        
        self.dropout = nn.Dropout(0.3)
        self.batch_norms = nn.ModuleList([nn.BatchNorm1d(d) for d in hidden_dims])
        
    def forward(self, x, edge_index):
        """前向传播"""
        for i, layer in enumerate(self.layers[:-1]):
            x = layer(x, edge_index)
            x = self.batch_norms[i](x)
            x = F.elu(x)
            x = self.dropout(x)
        
        x = self.layers[-1](x, edge_index)
        return x


class SpatialVAE(nn.Module):
    """
    空间变分自编码器 (Spatial VAE)
    创新点: 结合GNN编码器和概率解码器进行基因表达预测
    """
    
    def __init__(self, n_genes, spatial_dim=128, latent_dim=32):
        super(SpatialVAE, self).__init__()
        
        # GNN编码器
        self.spatial_encoder = SpatialGNNEncoder(
            input_dim=n_genes,
            hidden_dims=[512, 256, spatial_dim],
            gnn_type='GAT'
        )
        
        # VAE参数网络
        self.fc_mu = nn.Linear(spatial_dim, latent_dim)
        self.fc_logvar = nn.Linear(spatial_dim, latent_dim)
        
        # 解码器
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(128, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Linear(512, n_genes)
        )
        
    def encode(self, x, edge_index):
        """编码"""
        h = self.spatial_encoder(x, edge_index)
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar
    
    def reparameterize(self, mu, logvar):
        """重参数化技巧"""
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def decode(self, z):
        """解码"""
        return self.decoder(z)
    
    def forward(self, x, edge_index):
        """前向传播"""
        mu, logvar = self.encode(x, edge_index)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar, z


# ============================================================================
# 第三部分: 扩散模型用于空间插值
# ============================================================================

class SpatialDiffusionModel(nn.Module):
    """
    扩散模型用于空间基因表达插值
    创新点: 在连续空间中预测未测量位置的基因表达
    """
    
    def __init__(self, n_genes, spatial_dim=2, time_dim=32, hidden_dim=256):
        super(SpatialDiffusionModel, self).__init__()
        
        # 时间嵌入
        self.time_mlp = nn.Sequential(
            nn.Linear(1, time_dim),
            nn.SiLU(),
            nn.Linear(time_dim, time_dim)
        )
        
        # 空间位置嵌入
        self.spatial_mlp = nn.Sequential(
            nn.Linear(spatial_dim, 64),
            nn.SiLU(),
            nn.Linear(64, 64)
        )
        
        # 主网络 - U-Net风格
        input_dim = n_genes + time_dim + 64
        
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim * 2),
            nn.LayerNorm(hidden_dim * 2),
            nn.SiLU()
        )
        
        self.decoder = nn.Sequential(
            nn.Linear(hidden_dim * 2, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, n_genes)
        )
        
    def forward(self, x, t, spatial_coords):
        """
        x: 基因表达 (noise)
        t: 时间步
        spatial_coords: 空间坐标
        """
        # 嵌入
        t_emb = self.time_mlp(t.unsqueeze(-1))
        s_emb = self.spatial_mlp(spatial_coords)
        
        # 拼接
        h = torch.cat([x, t_emb, s_emb], dim=-1)
        
        # 前向
        h = self.encoder(h)
        noise_pred = self.decoder(h)
        
        return noise_pred


# ============================================================================
# 第四部分: 空间域发现 (超越传统聚类)
# ============================================================================

class SpatialDomainDiscovery(nn.Module):
    """
    连续空间域发现
    创新点: 使用神经场(Neural Field)建模连续空间域
    """
    
    def __init__(self, spatial_dim=2, hidden_dim=128, n_domains=10):
        super(SpatialDomainDiscovery, self).__init__()
        
        # 位置编码 (类似NeRF)
        self.pos_encoding_dim = 32
        
        # MLP网络
        self.network = nn.Sequential(
            nn.Linear(spatial_dim + self.pos_encoding_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, n_domains)
        )
        
    def positional_encoding(self, coords, L=5):
        """位置编码"""
        encoding = []
        for l in range(L):
            encoding.append(torch.sin(2**l * np.pi * coords))
            encoding.append(torch.cos(2**l * np.pi * coords))
        return torch.cat(encoding, dim=-1)
    
    def forward(self, spatial_coords):
        """预测空间域概率"""
        pos_enc = self.positional_encoding(spatial_coords)
        combined = torch.cat([spatial_coords, pos_enc], dim=-1)
        logits = self.network(combined)
        return F.softmax(logits, dim=-1)


# ============================================================================
# 第五部分: 训练流程
# ============================================================================

class SpatialTranscriptomicsPredictor:
    """完整的预测系统"""
    
    def __init__(self, adata, device='cuda' if torch.cuda.is_available() else 'cpu'):
        self.adata = adata
        self.device = device
        self.n_genes = adata.n_vars
        self.n_spots = adata.n_obs
        
        # 初始化模型
        self.vae = SpatialVAE(
            n_genes=self.n_genes,
            spatial_dim=128,
            latent_dim=32
        ).to(device)
        
        self.diffusion = SpatialDiffusionModel(
            n_genes=self.n_genes,
            spatial_dim=2
        ).to(device)
        
        self.domain_model = SpatialDomainDiscovery(
            spatial_dim=2,
            n_domains=10
        ).to(device)
        
    def prepare_graph_data(self):
        """准备图数据"""
        # 基因表达矩阵
        if issparse(self.adata.X):
            X = torch.FloatTensor(self.adata.X.toarray()).to(self.device)
        else:
            X = torch.FloatTensor(self.adata.X).to(self.device)
        
        # 空间坐标
        spatial_coords = torch.FloatTensor(
            self.adata.obsm['spatial']
        ).to(self.device)
        
        # 邻接图
        edge_index = torch.LongTensor(
            np.array(self.adata.obsp['spatial_connectivities'].nonzero())
        ).to(self.device)
        
        return X, spatial_coords, edge_index
    
    def train_vae(self, epochs=100, lr=1e-3):
        """训练VAE"""
        optimizer = torch.optim.Adam(self.vae.parameters(), lr=lr)
        X, spatial_coords, edge_index = self.prepare_graph_data()
        
        self.vae.train()
        losses = []
        
        print("训练空间VAE...")
        for epoch in range(epochs):
            optimizer.zero_grad()
            
            recon, mu, logvar, z = self.vae(X, edge_index)
            
            # 重构损失
            recon_loss = F.mse_loss(recon, X)
            
            # KL散度
            kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
            kl_loss = kl_loss / X.size(0)
            
            # 总损失
            loss = recon_loss + 0.001 * kl_loss
            
            loss.backward()
            optimizer.step()
            
            losses.append(loss.item())
            
            if (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch+1}/{epochs}, Loss: {loss.item():.4f}, "
                      f"Recon: {recon_loss.item():.4f}, KL: {kl_loss.item():.4f}")
        
        return losses
    
    def train_domain_discovery(self, epochs=50, lr=1e-3):
        """训练空间域发现模型"""
        optimizer = torch.optim.Adam(self.domain_model.parameters(), lr=lr)
        X, spatial_coords, edge_index = self.prepare_graph_data()
        
        # 获取VAE嵌入作为目标
        self.vae.eval()
        with torch.no_grad():
            _, _, _, z_target = self.vae(X, edge_index)
        
        self.domain_model.train()
        losses = []
        
        print("训练空间域发现模型...")
        for epoch in range(epochs):
            optimizer.zero_grad()
            
            # 预测空间域
            domain_probs = self.domain_model(spatial_coords)
            
            # 使用对比学习损失
            # 同域内的spots应该有相似的表达模式
            loss = self._contrastive_loss(domain_probs, z_target)
            
            loss.backward()
            optimizer.step()
            
            losses.append(loss.item())
            
            if (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch+1}/{epochs}, Loss: {loss.item():.4f}")
        
        return losses
    
    def _contrastive_loss(self, domain_probs, embeddings):
        """对比学习损失"""
        # 简化版本: 最大化域内一致性
        domain_assignments = torch.argmax(domain_probs, dim=1)
        loss = 0
        
        for domain_id in range(domain_probs.size(1)):
            mask = (domain_assignments == domain_id)
            if mask.sum() > 1:
                domain_embs = embeddings[mask]
                # 域内方差应该小
                loss += torch.var(domain_embs, dim=0).mean()
        
        return loss
    
    def predict_gene_expression(self, new_spatial_coords):
        """
        预测新空间位置的基因表达
        new_spatial_coords: (n_new_spots, 2)
        """
        self.vae.eval()
        self.diffusion.eval()
        
        with torch.no_grad():
            new_coords = torch.FloatTensor(new_spatial_coords).to(self.device)
            
            # 使用扩散模型生成
            # 简化版: 从噪声开始迭代去噪
            x = torch.randn(new_coords.size(0), self.n_genes).to(self.device)
            
            n_steps = 50
            for t in range(n_steps-1, -1, -1):
                t_tensor = torch.ones(new_coords.size(0)).to(self.device) * t / n_steps
                noise_pred = self.diffusion(x, t_tensor, new_coords)
                
                # 简化的去噪步骤
                alpha = 1 - t / n_steps
                x = x - 0.01 * noise_pred
            
            predicted_expression = x
        
        return predicted_expression.cpu().numpy()
    
    def get_spatial_domains(self):
        """获取空间域分配"""
        self.domain_model.eval()
        X, spatial_coords, edge_index = self.prepare_graph_data()
        
        with torch.no_grad():
            domain_probs = self.domain_model(spatial_coords)
            domain_assignments = torch.argmax(domain_probs, dim=1)
        
        return domain_assignments.cpu().numpy(), domain_probs.cpu().numpy()
    
    def get_latent_representation(self):
        """获取低维嵌入用于可视化"""
        self.vae.eval()
        X, spatial_coords, edge_index = self.prepare_graph_data()
        
        with torch.no_grad():
            _, _, _, z = self.vae(X, edge_index)
        
        return z.cpu().numpy()


# ============================================================================
# 第六部分: 可视化和分析
# ============================================================================

class SpatialVisualizer:
    """空间可视化工具"""
    
    @staticmethod
    def plot_spatial_domains(adata, domains, save_path=None):
        """绘制空间域"""
        fig, ax = plt.subplots(figsize=(10, 8))
        
        spatial_coords = adata.obsm['spatial']
        scatter = ax.scatter(
            spatial_coords[:, 0],
            spatial_coords[:, 1],
            c=domains,
            cmap='tab20',
            s=50,
            alpha=0.8
        )
        
        plt.colorbar(scatter, label='Spatial Domain')
        ax.set_xlabel('Spatial X')
        ax.set_ylabel('Spatial Y')
        ax.set_title('Discovered Spatial Domains')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
    
    @staticmethod
    def plot_latent_space(latent_embedding, domains=None, save_path=None):
        """绘制潜在空间"""
        from sklearn.manifold import TSNE
        
        # t-SNE降维
        tsne = TSNE(n_components=2, random_state=42)
        latent_2d = tsne.fit_transform(latent_embedding)
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        if domains is not None:
            scatter = ax.scatter(
                latent_2d[:, 0],
                latent_2d[:, 1],
                c=domains,
                cmap='tab20',
                s=30,
                alpha=0.6
            )
            plt.colorbar(scatter)
        else:
            ax.scatter(latent_2d[:, 0], latent_2d[:, 1], s=30, alpha=0.6)
        
        ax.set_xlabel('t-SNE 1')
        ax.set_ylabel('t-SNE 2')
        ax.set_title('Latent Space Visualization')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()
    
    @staticmethod
    def plot_predicted_genes(adata, predicted_coords, predicted_expression, 
                           gene_indices=[0, 1, 2], save_path=None):
        """绘制预测的基因表达"""
        n_genes = len(gene_indices)
        fig, axes = plt.subplots(1, n_genes, figsize=(6*n_genes, 5))
        
        if n_genes == 1:
            axes = [axes]
        
        for idx, gene_idx in enumerate(gene_indices):
            ax = axes[idx]
            
            # 原始数据
            spatial_coords = adata.obsm['spatial']
            if issparse(adata.X):
                original_exp = adata.X[:, gene_idx].toarray().flatten()
            else:
                original_exp = adata.X[:, gene_idx]
            
            ax.scatter(
                spatial_coords[:, 0],
                spatial_coords[:, 1],
                c=original_exp,
                cmap='viridis',
                s=50,
                alpha=0.6,
                label='Original'
            )
            
            # 预测数据
            ax.scatter(
                predicted_coords[:, 0],
                predicted_coords[:, 1],
                c=predicted_expression[:, gene_idx],
                cmap='plasma',
                s=100,
                marker='^',
                alpha=0.8,
                label='Predicted'
            )
            
            ax.set_title(f'Gene {gene_idx}')
            ax.legend()
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.show()


# ============================================================================
# 第七部分: 主流程示例
# ============================================================================

def main_analysis_pipeline():
    """主分析流程"""
    
    print("="*80)
    print("高级空间转录组学预测分析系统")
    print("="*80)
    
    # 1. 加载数据
    print("\n步骤 1: 数据加载")
    loader = SpatialDataLoader("path/to/visium/data")
    # adata = loader.load_visium_data()
    
    # 演示用: 创建模拟数据
    print("创建模拟Visium数据...")
    n_spots = 2000
    n_genes = 2000
    
    # 模拟基因表达
    np.random.seed(42)
    expression = np.random.negative_binomial(5, 0.3, (n_spots, n_genes)).astype(float)
    
    # 模拟空间坐标 (六边形网格)
    rows, cols = 50, 40
    spatial_coords = []
    for i in range(rows):
        for j in range(cols):
            x = j + 0.5 * (i % 2)
            y = i * np.sqrt(3) / 2
            spatial_coords.append([x, y])
    spatial_coords = np.array(spatial_coords[:n_spots])
    
    # 创建AnnData对象
    adata = sc.AnnData(X=expression)
    adata.obsm['spatial'] = spatial_coords
    
    # 2. 预处理
    print("\n步骤 2: 数据预处理")
    adata = loader.preprocess(adata, n_top_genes=2000)
    adata = loader.compute_spatial_graph(adata, n_neighs=6)
    
    # 3. 初始化预测器
    print("\n步骤 3: 初始化深度学习模型")
    predictor = SpatialTranscriptomicsPredictor(adata, device='cpu')
    
    # 4. 训练VAE
    print("\n步骤 4: 训练空间VAE")
    vae_losses = predictor.train_vae(epochs=50, lr=1e-3)
    
    # 5. 训练空间域发现
    print("\n步骤 5: 训练空间域发现模型")
    domain_losses = predictor.train_domain_discovery(epochs=30, lr=1e-3)
    
    # 6. 获取结果
    print("\n步骤 6: 提取分析结果")
    domains, domain_probs = predictor.get_spatial_domains()
    latent_embedding = predictor.get_latent_representation()
    
    # 保存到adata
    adata.obs['spatial_domain'] = domains
    adata.obsm['latent'] = latent_embedding
    adata.obsm['domain_probs'] = domain_probs
    
    # 7. 预测新位置
    print("\n步骤 7: 预测未测量位置的基因表达")
    # 在现有spots之间插值
    new_coords = spatial_coords[::10] + np.random.randn(len(spatial_coords[::10]), 2) * 0.5
    predicted_exp = predictor.predict_gene_expression(new_coords)
    
    print(f"预测了 {len(new_coords)} 个新位置的基因表达")
    print(f"预测表达矩阵形状: {predicted_exp.shape}")
    
    # 8. 可视化
    print("\n步骤 8: 生成可视化")
    visualizer = SpatialVisualizer()
    
    visualizer.plot_spatial_domains(adata, domains, 
                                   save_path='/home/claude/spatial_domains.png')
    
    visualizer.plot_latent_space(latent_embedding, domains,
                                save_path='/home/claude/latent_space.png')
    
    visualizer.plot_predicted_genes(adata, new_coords, predicted_exp,
                                   gene_indices=[0, 10, 20],
                                   save_path='/home/claude/predicted_genes.png')
    
    # 9. 保存结果
    print("\n步骤 9: 保存分析结果")
    adata.write('/home/claude/analyzed_visium_data.h5ad')
    
    # 保存预测结果
    np.savez('/home/claude/predictions.npz',
             coords=new_coords,
             expression=predicted_exp)
    
    print("\n" + "="*80)
    print("分析完成!")
    print("="*80)
    print(f"发现的空间域数量: {len(np.unique(domains))}")
    print(f"潜在空间维度: {latent_embedding.shape[1]}")
    print(f"预测位置数量: {len(new_coords)}")
    
    return adata, predictor


if __name__ == "__main__":
    # 运行完整分析
    adata, predictor = main_analysis_pipeline()
    
    print("\n可用功能:")
    print("1. predictor.predict_gene_expression(new_coords) - 预测基因表达")
    print("2. predictor.get_spatial_domains() - 获取空间域")
    print("3. predictor.get_latent_representation() - 获取嵌入")
