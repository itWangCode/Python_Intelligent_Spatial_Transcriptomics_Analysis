"""
DeepST 核心实现 - macOS Intel 兼容版本
支持 CPU 和 GPU 自动检测
适用于博士研究和 AIDD 应用

作者: 扩展版本用于研究
版本: 2.0
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv
from torch_geometric.data import Data
import numpy as np
from sklearn.neighbors import KDTree, BallTree
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from typing import Optional, Tuple, Dict, List
import warnings

warnings.filterwarnings('ignore')


class DeviceManager:
    """设备管理器 - 自动检测和管理 CPU/GPU"""
    
    def __init__(self, prefer_gpu: bool = True):
        self.prefer_gpu = prefer_gpu
        self.device = self._get_device()
        
    def _get_device(self) -> torch.device:
        """自动检测最佳设备"""
        if self.prefer_gpu:
            # 检查 CUDA (NVIDIA GPU)
            if torch.cuda.is_available():
                device = torch.device('cuda')
                print(f"✓ 使用 CUDA GPU: {torch.cuda.get_device_name(0)}")
                return device
            
            # 检查 MPS (Apple Silicon GPU) - macOS Intel 不支持
            if hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
                device = torch.device('mps')
                print("✓ 使用 Apple MPS GPU")
                return device
        
        # 默认使用 CPU
        device = torch.device('cpu')
        print(f"✓ 使用 CPU ({torch.get_num_threads()} 线程)")
        return device
    
    def to_device(self, tensor):
        """将张量移动到设备"""
        return tensor.to(self.device)
    
    def __repr__(self):
        return f"DeviceManager(device={self.device})"


class SpatialGraphConstructor:
    """空间图构建器 - 构建spot之间的空间邻接关系"""
    
    def __init__(
        self,
        n_neighbors: int = 6,
        distance_type: str = 'KDTree',
        radius: Optional[float] = None
    ):
        """
        参数:
            n_neighbors: 每个spot的邻居数量
            distance_type: 距离计算方法 ('KDTree' 或 'BallTree')
            radius: 半径阈值 (可选)
        """
        self.n_neighbors = n_neighbors
        self.distance_type = distance_type
        self.radius = radius
        
    def build_graph(
        self,
        spatial_coords: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        构建空间邻接图
        
        参数:
            spatial_coords: 空间坐标矩阵 [n_spots, 2]
            
        返回:
            edge_index: 边索引 [2, n_edges]
            edge_weight: 边权重 [n_edges]
        """
        n_spots = spatial_coords.shape[0]
        
        # 选择距离计算方法
        if self.distance_type == 'KDTree':
            tree = KDTree(spatial_coords, leaf_size=40)
        else:
            tree = BallTree(spatial_coords, leaf_size=40)
        
        # 查询 k 近邻
        distances, indices = tree.query(
            spatial_coords,
            k=self.n_neighbors + 1  # +1 包含自己
        )
        
        # 构建边
        edge_list = []
        edge_weights = []
        
        for i in range(n_spots):
            for j, dist in zip(indices[i][1:], distances[i][1:]):  # 跳过自己
                if self.radius is None or dist <= self.radius:
                    edge_list.append([i, j])
                    # 高斯核权重
                    weight = np.exp(-dist**2 / (2 * np.mean(distances)**2))
                    edge_weights.append(weight)
        
        edge_index = np.array(edge_list).T
        edge_weight = np.array(edge_weights)
        
        print(f"✓ 构建空间图: {n_spots} 个节点, {len(edge_weights)} 条边")
        
        return edge_index, edge_weight
    
    def spatial_smoothing(
        self,
        gene_expression: np.ndarray,
        spatial_coords: np.ndarray,
        alpha: float = 0.5
    ) -> np.ndarray:
        """
        基于空间邻域的基因表达平滑
        
        参数:
            gene_expression: 基因表达矩阵 [n_spots, n_genes]
            spatial_coords: 空间坐标 [n_spots, 2]
            alpha: 平滑强度 (0-1)
            
        返回:
            smoothed_expression: 平滑后的表达矩阵
        """
        if self.distance_type == 'KDTree':
            tree = KDTree(spatial_coords)
        else:
            tree = BallTree(spatial_coords)
        
        distances, indices = tree.query(spatial_coords, k=self.n_neighbors + 1)
        
        smoothed = gene_expression.copy()
        
        for i in range(len(gene_expression)):
            # 计算邻居权重
            neighbor_dist = distances[i][1:]
            weights = np.exp(-neighbor_dist**2 / (2 * np.mean(distances)**2))
            weights = weights / weights.sum()
            
            # 加权平均
            neighbor_expr = gene_expression[indices[i][1:]]
            smoothed[i] = (1 - alpha) * gene_expression[i] + alpha * (weights[:, None] * neighbor_expr).sum(axis=0)
        
        print(f"✓ 空间平滑完成 (alpha={alpha})")
        return smoothed


class GNNEncoder(nn.Module):
    """图神经网络编码器"""
    
    def __init__(
        self,
        in_channels: int,
        hidden_channels: int = 512,
        out_channels: int = 128,
        num_layers: int = 2,
        dropout: float = 0.1,
        gnn_type: str = 'GCN'
    ):
        super().__init__()
        
        self.num_layers = num_layers
        self.dropout = dropout
        
        # 选择 GNN 类型
        GNN = GCNConv if gnn_type == 'GCN' else GATConv
        
        self.convs = nn.ModuleList()
        self.convs.append(GNN(in_channels, hidden_channels))
        
        for _ in range(num_layers - 2):
            self.convs.append(GNN(hidden_channels, hidden_channels))
        
        self.convs.append(GNN(hidden_channels, out_channels))
        
        self.batch_norms = nn.ModuleList([
            nn.BatchNorm1d(hidden_channels) for _ in range(num_layers - 1)
        ])
        
    def forward(self, x, edge_index, edge_weight=None):
        for i, conv in enumerate(self.convs[:-1]):
            x = conv(x, edge_index, edge_weight=edge_weight)
            x = self.batch_norms[i](x)
            x = F.relu(x)
            x = F.dropout(x, p=self.dropout, training=self.training)
        
        x = self.convs[-1](x, edge_index, edge_weight=edge_weight)
        return x


class DenoisingAutoencoder(nn.Module):
    """去噪自编码器"""
    
    def __init__(
        self,
        in_features: int,
        hidden_dims: List[int] = [512, 256, 128],
        latent_dim: int = 128,
        dropout: float = 0.1
    ):
        super().__init__()
        
        # 编码器
        encoder_layers = []
        prev_dim = in_features
        for hidden_dim in hidden_dims:
            encoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ])
            prev_dim = hidden_dim
        encoder_layers.append(nn.Linear(prev_dim, latent_dim))
        self.encoder = nn.Sequential(*encoder_layers)
        
        # 解码器
        decoder_layers = []
        prev_dim = latent_dim
        for hidden_dim in reversed(hidden_dims):
            decoder_layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout)
            ])
            prev_dim = hidden_dim
        decoder_layers.append(nn.Linear(prev_dim, in_features))
        self.decoder = nn.Sequential(*decoder_layers)
    
    def forward(self, x, noise_factor: float = 0.1):
        # 添加噪声
        if self.training:
            noise = torch.randn_like(x) * noise_factor
            x_noisy = x + noise
        else:
            x_noisy = x
        
        # 编码-解码
        latent = self.encoder(x_noisy)
        reconstruction = self.decoder(latent)
        
        return latent, reconstruction


class DomainAdversarialNetwork(nn.Module):
    """域对抗网络 - 用于批次效应校正"""
    
    def __init__(
        self,
        in_features: int,
        n_domains: int,
        hidden_dim: int = 128
    ):
        super().__init__()
        
        self.domain_classifier = nn.Sequential(
            nn.Linear(in_features, hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Linear(hidden_dim // 2, n_domains)
        )
    
    def forward(self, x, alpha: float = 1.0):
        # 梯度反转层
        x = GradientReversal.apply(x, alpha)
        domain_pred = self.domain_classifier(x)
        return domain_pred


class GradientReversal(torch.autograd.Function):
    """梯度反转层"""
    
    @staticmethod
    def forward(ctx, x, alpha):
        ctx.alpha = alpha
        return x.view_as(x)
    
    @staticmethod
    def backward(ctx, grad_output):
        return grad_output.neg() * ctx.alpha, None


class DeepSTModel(nn.Module):
    """DeepST 主模型"""
    
    def __init__(
        self,
        n_genes: int,
        n_spots: int,
        hidden_dim: int = 512,
        latent_dim: int = 128,
        n_domains: Optional[int] = None,
        use_batch_correction: bool = False,
        gnn_type: str = 'GCN'
    ):
        super().__init__()
        
        self.n_genes = n_genes
        self.n_spots = n_spots
        self.latent_dim = latent_dim
        self.use_batch_correction = use_batch_correction
        
        # GNN 编码器
        self.gnn_encoder = GNNEncoder(
            in_channels=n_genes,
            hidden_channels=hidden_dim,
            out_channels=latent_dim,
            gnn_type=gnn_type
        )
        
        # 去噪自编码器
        self.dae = DenoisingAutoencoder(
            in_features=n_genes,
            hidden_dims=[hidden_dim, hidden_dim // 2, hidden_dim // 4],
            latent_dim=latent_dim
        )
        
        # 域对抗网络 (可选)
        if use_batch_correction and n_domains is not None:
            self.dan = DomainAdversarialNetwork(
                in_features=latent_dim,
                n_domains=n_domains
            )
        else:
            self.dan = None
        
        # 融合层
        self.fusion = nn.Sequential(
            nn.Linear(latent_dim * 2, latent_dim),
            nn.BatchNorm1d(latent_dim),
            nn.ReLU()
        )
    
    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        edge_weight: Optional[torch.Tensor] = None,
        alpha: float = 1.0
    ):
        # GNN 分支
        gnn_latent = self.gnn_encoder(x, edge_index, edge_weight)
        
        # DAE 分支
        dae_latent, reconstruction = self.dae(x)
        
        # 融合两个分支
        combined_latent = torch.cat([gnn_latent, dae_latent], dim=1)
        fused_latent = self.fusion(combined_latent)
        
        # 域对抗 (如果启用)
        domain_pred = None
        if self.dan is not None and self.training:
            domain_pred = self.dan(fused_latent, alpha)
        
        return fused_latent, reconstruction, domain_pred


class DeepSTTrainer:
    """DeepST 训练器"""
    
    def __init__(
        self,
        model: DeepSTModel,
        device_manager: DeviceManager,
        learning_rate: float = 1e-3,
        weight_decay: float = 1e-4
    ):
        self.model = model.to(device_manager.device)
        self.device = device_manager.device
        self.device_manager = device_manager
        
        self.optimizer = torch.optim.Adam(
            model.parameters(),
            lr=learning_rate,
            weight_decay=weight_decay
        )
        
        self.scheduler = torch.optim.lr_scheduler.StepLR(
            self.optimizer,
            step_size=100,
            gamma=0.9
        )
    
    def train_epoch(
        self,
        data: Data,
        domain_labels: Optional[torch.Tensor] = None,
        alpha: float = 1.0,
        lambda_recon: float = 1.0,
        lambda_domain: float = 0.1
    ) -> Dict[str, float]:
        """训练一个 epoch"""
        self.model.train()
        
        # 准备数据
        x = self.device_manager.to_device(data.x)
        edge_index = self.device_manager.to_device(data.edge_index)
        edge_weight = self.device_manager.to_device(data.edge_attr) if data.edge_attr is not None else None
        
        # 前向传播
        latent, reconstruction, domain_pred = self.model(
            x, edge_index, edge_weight, alpha
        )
        
        # 重构损失
        recon_loss = F.mse_loss(reconstruction, x)
        
        # 总损失
        total_loss = lambda_recon * recon_loss
        
        # 域对抗损失 (如果启用)
        domain_loss = torch.tensor(0.0).to(self.device)
        if domain_pred is not None and domain_labels is not None:
            domain_labels = self.device_manager.to_device(domain_labels)
            domain_loss = F.cross_entropy(domain_pred, domain_labels)
            total_loss = total_loss + lambda_domain * domain_loss
        
        # 反向传播
        self.optimizer.zero_grad()
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
        self.optimizer.step()
        
        return {
            'total_loss': total_loss.item(),
            'recon_loss': recon_loss.item(),
            'domain_loss': domain_loss.item() if isinstance(domain_loss, torch.Tensor) else 0.0
        }
    
    def get_embeddings(self, data: Data) -> np.ndarray:
        """获取嵌入表示"""
        self.model.eval()
        
        with torch.no_grad():
            x = self.device_manager.to_device(data.x)
            edge_index = self.device_manager.to_device(data.edge_index)
            edge_weight = self.device_manager.to_device(data.edge_attr) if data.edge_attr is not None else None
            
            latent, _, _ = self.model(x, edge_index, edge_weight)
            
        return latent.cpu().numpy()


class SpatialDomainIdentifier:
    """空间域识别器"""
    
    def __init__(
        self,
        n_domains: int,
        use_spatial_refine: bool = True,
        refine_iterations: int = 10
    ):
        self.n_domains = n_domains
        self.use_spatial_refine = use_spatial_refine
        self.refine_iterations = refine_iterations
    
    def identify_domains(
        self,
        embeddings: np.ndarray,
        spatial_coords: np.ndarray
    ) -> np.ndarray:
        """
        识别空间域
        
        参数:
            embeddings: 低维嵌入 [n_spots, latent_dim]
            spatial_coords: 空间坐标 [n_spots, 2]
            
        返回:
            domain_labels: 域标签 [n_spots]
        """
        # 初始聚类
        kmeans = KMeans(n_clusters=self.n_domains, random_state=42, n_init=20)
        initial_labels = kmeans.fit_predict(embeddings)
        
        if not self.use_spatial_refine:
            return initial_labels
        
        # 空间精炼
        refined_labels = self._spatial_refine(
            initial_labels,
            spatial_coords
        )
        
        return refined_labels
    
    def _spatial_refine(
        self,
        labels: np.ndarray,
        spatial_coords: np.ndarray
    ) -> np.ndarray:
        """基于空间邻域的标签精炼"""
        refined = labels.copy()
        
        # 构建 KNN 图
        tree = KDTree(spatial_coords)
        
        for iteration in range(self.refine_iterations):
            changed = 0
            
            for i in range(len(labels)):
                # 找到邻居
                _, indices = tree.query([spatial_coords[i]], k=7)
                neighbors = indices[0][1:]  # 排除自己
                
                # 多数投票
                neighbor_labels = refined[neighbors]
                unique, counts = np.unique(neighbor_labels, return_counts=True)
                majority_label = unique[np.argmax(counts)]
                
                if refined[i] != majority_label:
                    refined[i] = majority_label
                    changed += 1
            
            if changed == 0:
                break
        
        print(f"✓ 空间精炼完成 ({iteration + 1} 次迭代)")
        return refined


def preprocess_data(
    gene_expression: np.ndarray,
    n_top_genes: int = 3000,
    pca_components: int = 200
) -> np.ndarray:
    """
    数据预处理
    
    参数:
        gene_expression: 基因表达矩阵 [n_spots, n_genes]
        n_top_genes: 保留的高变基因数量
        pca_components: PCA降维维数
        
    返回:
        processed_data: 预处理后的数据
    """
    # 对数变换
    expression = np.log1p(gene_expression)
    
    # 标准化
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler()
    expression = scaler.fit_transform(expression)
    
    # 高变基因选择
    if expression.shape[1] > n_top_genes:
        variances = np.var(expression, axis=0)
        top_indices = np.argsort(variances)[-n_top_genes:]
        expression = expression[:, top_indices]
        print(f"✓ 选择前 {n_top_genes} 个高变基因")
    
    # PCA 降维
    if pca_components < expression.shape[1]:
        pca = PCA(n_components=pca_components, random_state=42)
        expression = pca.fit_transform(expression)
        print(f"✓ PCA 降维至 {pca_components} 维 (解释方差: {pca.explained_variance_ratio_.sum():.2%})")
    
    return expression


if __name__ == "__main__":
    print("DeepST 核心模块加载成功!")
    print("\n可用组件:")
    print("  - DeviceManager: CPU/GPU 自动管理")
    print("  - SpatialGraphConstructor: 空间图构建")
    print("  - DeepSTModel: 主模型")
    print("  - DeepSTTrainer: 训练器")
    print("  - SpatialDomainIdentifier: 空间域识别")
    
    # 测试设备检测
    dm = DeviceManager(prefer_gpu=True)
    print(f"\n当前设备: {dm.device}")
