# 🧬 高级空间转录组学深度学习预测分析系统

**作者**: 计算机博士研究项目  
**数据源**: 10x Genomics Human Breast Cancer Visium Dataset

---

## 📖 项目概述

本项目实现了一个**创新的深度学习驱动的空间转录组学分析框架**，超越传统的聚类和统计方法，使用最前沿的图神经网络、变分自编码器和扩散模型进行空间基因表达预测和分析。

### 🎯 核心创新点

1. **图神经网络(GNN)空间嵌入**
   - 使用GAT (Graph Attention Networks) 捕获空间依赖性
   - 多头注意力机制学习spot之间的相互作用
   - 优于传统的欧氏距离方法

2. **空间变分自编码器(Spatial VAE)**
   - 结合GNN编码器和概率解码器
   - 学习低维潜在表示同时保留空间信息
   - 支持不确定性量化

3. **扩散模型空间插值**
   - 在连续空间中预测未测量位置的基因表达
   - 基于DDPM (Denoising Diffusion Probabilistic Models)
   - 超越线性插值的非线性建模能力

4. **神经场(Neural Field)空间域发现**
   - 类似NeRF的连续空间域建模
   - 位置编码捕获高频空间模式
   - 超越传统离散聚类的连续域表示

---

## 🏗️ 架构设计

```
输入: Visium空间转录组数据
  ↓
[数据预处理模块]
  - 质量控制
  - 归一化
  - 高变基因选择
  - 空间邻接图构建
  ↓
[GNN空间编码器]
  - GAT多层网络
  - 捕获局部和全局空间模式
  ↓
[VAE潜在空间学习]
  - μ, σ 参数化
  - 重参数化采样
  - 概率解码
  ↓
[扩散模型预测]
  - 时间步嵌入
  - 空间位置嵌入
  - 迭代去噪预测
  ↓
[空间域发现]
  - Neural Field建模
  - 位置编码
  - 连续域概率
  ↓
输出: 
  - 预测的基因表达
  - 空间域分配
  - 低维嵌入
  - 可视化结果
```

---

## 🚀 快速开始

### 1. 环境配置

```bash
# 创建conda环境
conda create -n spatial_dl python=3.9
conda activate spatial_dl

# 安装依赖
pip install -r requirements.txt

# 安装PyTorch Geometric (根据CUDA版本)
pip install torch-geometric
```

### 2. 下载数据

从10x Genomics下载数据：
```bash
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5

wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz

tar -xzf V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz
```

### 3. 运行分析

```python
# 基础用法
python spatial_transcriptomics_prediction.py

# 或者使用自定义参数
from spatial_transcriptomics_prediction import *

# 加载数据
loader = SpatialDataLoader("path/to/visium/data")
adata = loader.load_visium_data()
adata = loader.preprocess(adata)
adata = loader.compute_spatial_graph(adata)

# 初始化预测器
predictor = SpatialTranscriptomicsPredictor(adata)

# 训练模型
predictor.train_vae(epochs=100)
predictor.train_domain_discovery(epochs=50)

# 获取结果
domains, domain_probs = predictor.get_spatial_domains()
latent = predictor.get_latent_representation()

# 预测新位置
new_coords = np.array([[10, 20], [15, 25], [20, 30]])
predicted_exp = predictor.predict_gene_expression(new_coords)
```

---

## 📊 方法详解

### 1. 图神经网络空间编码器

**理论基础**:
- 空间转录组数据天然形成图结构
- 每个spot是节点,空间邻接关系是边
- GNN可以聚合邻居信息学习空间模式

**实现**:
```python
class SpatialGNNEncoder(nn.Module):
    # 多层GAT网络
    # 每层: h^(l+1) = σ(∑ α_ij W h^(l))
    # α_ij: 注意力权重
```

**优势**:
- 自适应学习邻居重要性
- 多跳信息传播
- 处理不规则空间结构

### 2. 空间变分自编码器

**理论基础**:
- 学习数据的概率分布 p(x|z)
- 潜在变量z捕获关键空间模式
- ELBO优化: L = E[log p(x|z)] - KL(q(z|x)||p(z))

**实现**:
```python
# 编码: q(z|x, G) 其中G是空间图
mu, logvar = encoder(x, edge_index)
z = mu + eps * exp(0.5 * logvar)

# 解码: p(x|z)
x_recon = decoder(z)
```

**优势**:
- 学习紧凑表示
- 不确定性量化
- 生成新样本

### 3. 扩散模型空间插值

**理论基础**:
- 前向过程: x_t = √(α_t) x_0 + √(1-α_t) ε
- 反向过程: 学习 ε_θ(x_t, t, s) 预测噪声
- 条件生成: 给定空间位置s生成表达x

**实现**:
```python
# 训练: 预测添加的噪声
noise_pred = model(x_t, t, spatial_coords)
loss = ||noise - noise_pred||^2

# 推理: 迭代去噪
for t in [T, T-1, ..., 1]:
    noise = model(x_t, t, new_coords)
    x_t = denoise_step(x_t, noise, t)
```

**优势**:
- 高质量样本生成
- 捕获复杂分布
- 超越线性插值

### 4. Neural Field空间域发现

**理论基础**:
- 将空间域建模为连续函数 f: R^2 → R^K
- 位置编码捕获高频细节
- 类似NeRF的隐式表示

**实现**:
```python
# 位置编码
pos_enc = [sin(2^l π x), cos(2^l π x)] for l in [0..L]

# MLP预测域概率
domain_prob = MLP([x, y, pos_enc]) → softmax
```

**优势**:
- 连续空间表示
- 任意分辨率查询
- 捕获复杂边界

---

## 📈 评估指标

### 1. 重构质量
- **MSE**: 平均平方误差
- **Pearson相关**: 基因表达相关性
- **Cosine相似度**: 表达模式相似度

### 2. 空间结构保持
- **Moran's I**: 空间自相关
- **Geary's C**: 空间异质性
- **邻域保持率**: k-NN一致性

### 3. 预测准确度
- **插值误差**: 与已知值比较
- **交叉验证**: 留一法验证
- **生物学一致性**: 标记基因表达模式

---

## 🔬 高级应用

### 1. 空间轨迹推断

```python
# 使用潜在空间学习细胞状态轨迹
from sklearn.manifold import TSNE
import numpy as np

latent = predictor.get_latent_representation()
trajectory = fit_principal_curve(latent)
adata.obs['pseudotime'] = trajectory
```

### 2. 空间细胞类型解卷积

```python
# 结合单细胞数据进行空间解卷积
from cell_type_deconvolution import SpatialDeconv

deconv = SpatialDeconv(
    spatial_adata=adata,
    sc_reference=sc_ref,
    method='DestVI'  # 或 'Tangram', 'RCTD'
)
cell_proportions = deconv.fit_transform()
```

### 3. 配体-受体相互作用分析

```python
# 预测空间细胞通讯
import squidpy as sq

sq.gr.ligrec(
    adata,
    cluster_key='spatial_domain',
    copy=False
)
```

### 4. 多样本整合分析

```python
# 整合多个Visium切片
from batch_integration import SpatialIntegration

integrator = SpatialIntegration(
    method='STAGATE'  # 或 'GraphST', 'SEDR'
)
integrated_adata = integrator.integrate([adata1, adata2, adata3])
```

---

## 📁 输出文件说明

运行后将生成以下文件:

```
outputs/
├── analyzed_visium_data.h5ad          # 完整分析结果
├── predictions.npz                     # 预测的基因表达
├── spatial_domains.png                 # 空间域可视化
├── latent_space.png                    # 潜在空间可视化
├── predicted_genes.png                 # 预测基因表达
└── analysis_report.html                # 交互式报告
```

### AnnData对象结构:
```
adata.X                    # 基因表达矩阵
adata.obsm['spatial']      # 空间坐标
adata.obsm['latent']       # VAE潜在嵌入
adata.obsm['domain_probs'] # 空间域概率
adata.obs['spatial_domain']# 空间域分配
adata.obsp['spatial_connectivities']  # 空间图
```

---

## 🆚 与现有方法对比

| 方法 | 空间建模 | 基因预测 | 域发现 | 计算效率 |
|------|---------|---------|--------|---------|
| **本方法** | GNN + 扩散模型 | ✅ 高精度 | Neural Field | 中等 |
| Seurat | k-NN | ❌ | 离散聚类 | 快 |
| Scanpy | 欧氏距离 | ❌ | Leiden | 快 |
| STAGATE | GNN | ❌ | 离散聚类 | 中等 |
| SpaGCN | GCN | ❌ | 离散聚类 | 快 |
| Tangram | 对齐 | ✅ 低精度 | ❌ | 慢 |

**本方法优势**:
- ✅ 端到端深度学习
- ✅ 连续空间建模
- ✅ 不确定性量化
- ✅ 生成式预测能力

---

## 🔧 参数调优指南

### VAE训练
```python
# 关键超参数
predictor.train_vae(
    epochs=100,           # 增加epochs提高质量
    lr=1e-3,             # 学习率
    beta=0.001           # KL权重 (控制正则化强度)
)
```

**调优建议**:
- 重构误差高 → 增加hidden_dims, 降低beta
- 潜在空间混乱 → 增加beta, 增加latent_dim
- 过拟合 → 增加dropout, 减少网络容量

### 扩散模型
```python
diffusion = SpatialDiffusionModel(
    n_genes=n_genes,
    time_dim=32,         # 时间嵌入维度
    hidden_dim=256       # 网络容量
)
```

**调优建议**:
- 预测质量低 → 增加hidden_dim, 增加去噪步数
- 过平滑 → 减少去噪步数
- 训练不稳定 → 降低学习率, 使用梯度裁剪

---

## 📚 相关论文和资源

### 核心方法论文:
1. **Graph Attention Networks** - Veličković et al. (2017)
2. **Auto-Encoding Variational Bayes** - Kingma & Welling (2013)
3. **Denoising Diffusion Probabilistic Models** - Ho et al. (2020)
4. **NeRF: Neural Radiance Fields** - Mildenhall et al. (2020)

### 空间转录组学:
1. **STAGATE** - Dong & Zhang (2022)
2. **SpaGCN** - Hu et al. (2021)
3. **Giotto** - Dries et al. (2021)
4. **Squidpy** - Palla et al. (2022)

### 代码参考:
- [PyTorch Geometric](https://github.com/pyg-team/pytorch_geometric)
- [Scanpy](https://github.com/scverse/scanpy)
- [Squidpy](https://github.com/scverse/squidpy)
- [STAGATE](https://github.com/QIFEIDKN/STAGATE)

---

## 🐛 常见问题

### Q1: CUDA内存不足
```python
# 解决方案: 使用CPU或减少batch size
predictor = SpatialTranscriptomicsPredictor(adata, device='cpu')

# 或者使用梯度累积
optimizer.zero_grad()
for i in range(accumulation_steps):
    loss = loss / accumulation_steps
    loss.backward()
optimizer.step()
```

### Q2: 空间图构建失败
```python
# 确保空间坐标已标准化
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
adata.obsm['spatial'] = scaler.fit_transform(adata.obsm['spatial'])
```

### Q3: 训练损失不收敛
```python
# 使用学习率调度器
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, mode='min', factor=0.5, patience=10
)
```

---

## 🤝 贡献指南

欢迎贡献! 请遵循以下步骤:

1. Fork本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启Pull Request

---

## 📄 许可证

本项目采用MIT许可证 - 详见LICENSE文件

---

## 📧 联系方式

如有问题或建议,请通过以下方式联系:
- GitHub Issues
- Email: [your-email@example.com]

---

## 🙏 致谢

- 10x Genomics提供的Visium数据
- Scanpy和Squidpy开发团队
- PyTorch Geometric社区

---

**最后更新**: 2025年1月
