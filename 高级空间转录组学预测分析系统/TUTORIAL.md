# 🚀 快速开始教程

## 5分钟快速上手

### 步骤 1: 安装依赖

```bash
# 创建环境
conda create -n spatial_ai python=3.9
conda activate spatial_ai

# 安装依赖
pip install -r requirements.txt

# 安装PyTorch Geometric (根据你的CUDA版本)
# CUDA 11.8
pip install torch-geometric torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.0.0+cu118.html

# CPU版本
pip install torch-geometric torch-scatter torch-sparse
```

### 步骤 2: 下载数据

```bash
# 方法1: 使用我们的下载脚本 (推荐)
python download_visium_data.py --dataset breast_cancer

# 方法2: 手动下载
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Breast_Cancer_Block_A_Section_1/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz
tar -xzf V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz
```

### 步骤 3: 运行分析

```bash
# 基础分析 (使用默认参数)
python spatial_transcriptomics_prediction.py

# 高级分析
python advanced_biological_analysis.py
```

### 步骤 4: 查看结果

分析完成后,检查以下文件:
```
outputs/
├── analyzed_visium_data.h5ad          # 完整分析结果
├── predictions.npz                     # 预测的基因表达
├── spatial_domains.png                 # 空间域可视化
├── latent_space.png                    # 潜在空间可视化
└── predicted_genes.png                 # 预测基因表达
```

---

## 🎓 进阶教程

### 自定义分析流程

```python
import scanpy as sc
from spatial_transcriptomics_prediction import *

# 1. 加载你自己的数据
adata = sc.read_visium("path/to/your/data")

# 2. 预处理
loader = SpatialDataLoader(".")
adata = loader.preprocess(adata, n_top_genes=3000)
adata = loader.compute_spatial_graph(adata, n_neighs=6)

# 3. 初始化模型
predictor = SpatialTranscriptomicsPredictor(
    adata,
    device='cuda'  # 使用GPU加速
)

# 4. 训练 (调整超参数)
predictor.train_vae(
    epochs=200,      # 更多epochs
    lr=1e-3          # 学习率
)

predictor.train_domain_discovery(
    epochs=100,
    lr=5e-4
)

# 5. 获取结果
domains, domain_probs = predictor.get_spatial_domains()
latent = predictor.get_latent_representation()

# 6. 预测新位置
# 创建规则网格
x = np.linspace(adata.obsm['spatial'][:, 0].min(), 
                adata.obsm['spatial'][:, 0].max(), 50)
y = np.linspace(adata.obsm['spatial'][:, 1].min(), 
                adata.obsm['spatial'][:, 1].max(), 50)
xx, yy = np.meshgrid(x, y)
new_coords = np.c_[xx.ravel(), yy.ravel()]

predicted_exp = predictor.predict_gene_expression(new_coords)

# 7. 可视化
visualizer = SpatialVisualizer()
visualizer.plot_spatial_domains(adata, domains)
visualizer.plot_predicted_genes(adata, new_coords, predicted_exp)
```

### 批处理多个样本

```python
import glob

# 批量处理所有Visium数据
data_paths = glob.glob("data/*/")

results = {}
for path in data_paths:
    sample_name = path.split('/')[-2]
    print(f"处理样本: {sample_name}")
    
    # 加载
    adata = sc.read_visium(path)
    
    # 分析
    predictor = SpatialTranscriptomicsPredictor(adata)
    predictor.train_vae(epochs=100)
    
    # 保存
    results[sample_name] = {
        'adata': adata,
        'predictor': predictor
    }
```

### 整合单细胞数据

```python
# 使用单细胞参考进行空间注释
import scanpy as sc

# 加载单细胞参考
sc_ref = sc.read_h5ad("reference_scRNA.h5ad")

# 使用Tangram进行空间映射
import tangram as tg

ad_map = tg.map_cells_to_space(
    adata_sc=sc_ref,
    adata_sp=adata,
    mode="cells",
    device="cuda"
)

# 将细胞类型映射到空间
tg.project_cell_annotations(ad_map, adata, annotation="cell_type")

# 现在可以可视化细胞类型的空间分布
sc.pl.spatial(adata, color="cell_type", spot_size=50)
```

---

## 🔧 常见问题解决

### Q1: 内存不足

```python
# 解决方案1: 使用CPU
predictor = SpatialTranscriptomicsPredictor(adata, device='cpu')

# 解决方案2: 减少基因数
adata = adata[:, adata.var['highly_variable']].copy()

# 解决方案3: 分批处理
batch_size = 500
for i in range(0, adata.n_obs, batch_size):
    batch = adata[i:i+batch_size]
    # 处理batch
```

### Q2: 训练不收敛

```python
# 尝试不同的学习率
for lr in [1e-2, 1e-3, 1e-4]:
    predictor = SpatialTranscriptomicsPredictor(adata)
    losses = predictor.train_vae(epochs=50, lr=lr)
    print(f"LR={lr}, Final Loss={losses[-1]}")

# 使用学习率调度
import torch.optim as optim

optimizer = optim.Adam(predictor.vae.parameters(), lr=1e-3)
scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min')
```

### Q3: 空间图构建失败

```python
# 检查坐标范围
print("坐标统计:")
print(adata.obsm['spatial'].min(axis=0))
print(adata.obsm['spatial'].max(axis=0))

# 标准化坐标
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
adata.obsm['spatial'] = scaler.fit_transform(adata.obsm['spatial'])

# 重新构建图
import squidpy as sq
sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type='generic')
```

---

## 📊 性能优化

### GPU加速

```python
# 确认GPU可用
import torch
print(f"CUDA available: {torch.cuda.is_available()}")
print(f"GPU device: {torch.cuda.get_device_name(0)}")

# 使用混合精度训练
from torch.cuda.amp import autocast, GradScaler

scaler = GradScaler()
with autocast():
    recon, mu, logvar, z = predictor.vae(X, edge_index)
    loss = criterion(recon, X, mu, logvar)
    
scaler.scale(loss).backward()
scaler.step(optimizer)
scaler.update()
```

### 并行处理

```python
# 使用多进程加速数据预处理
from multiprocessing import Pool

def process_sample(path):
    adata = sc.read_visium(path)
    # 预处理
    return adata

with Pool(4) as p:
    results = p.map(process_sample, data_paths)
```

---

## 📈 可视化技巧

### 交互式可视化

```python
import plotly.graph_objects as go

# 3D可视化潜在空间
fig = go.Figure(data=[go.Scatter3d(
    x=latent[:, 0],
    y=latent[:, 1],
    z=latent[:, 2],
    mode='markers',
    marker=dict(
        size=3,
        color=domains,
        colorscale='Viridis',
    )
)])

fig.update_layout(title='3D Latent Space')
fig.show()
```

### 发布质量图表

```python
import matplotlib.pyplot as plt
import seaborn as sns

# 设置风格
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 12

# 创建多面板图
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# 各种可视化...

plt.tight_layout()
plt.savefig('publication_figure.pdf', dpi=300, bbox_inches='tight')
```

---

## 🎯 实际应用案例

### 案例1: 肿瘤微环境分析

```python
# 识别肿瘤区域和微环境
domains, _ = predictor.get_spatial_domains()

# 标注区域
domain_names = {
    0: 'Tumor core',
    1: 'Tumor edge',
    2: 'Stroma',
    3: 'Immune infiltrate',
    4: 'Normal tissue'
}

adata.obs['region'] = [domain_names.get(d, 'Unknown') for d in domains]

# 分析肿瘤-免疫相互作用
from advanced_biological_analysis import AdvancedBiologicalAnalysis

analyzer = AdvancedBiologicalAnalysis(adata)
lr_results = analyzer.cell_cell_communication()
```

### 案例2: 发育轨迹重建

```python
# 推断发育伪时间
pseudotime = analyzer.spatial_trajectory_inference(start_domain=0)

# 识别轨迹相关基因
from scipy.stats import spearmanr

trajectory_genes = []
for gene in adata.var_names:
    gene_exp = adata[:, gene].X.toarray().flatten()
    corr, pval = spearmanr(pseudotime, gene_exp)
    if abs(corr) > 0.5 and pval < 0.01:
        trajectory_genes.append((gene, corr, pval))

# 按相关性排序
trajectory_genes.sort(key=lambda x: abs(x[1]), reverse=True)
```

### 案例3: 药物响应预测

```python
# 预测治疗后的空间表达模式
# (需要配对的治疗前后样本)

# 训练条件生成模型
class ConditionalSpatialVAE(SpatialVAE):
    def __init__(self, n_genes, n_conditions=2):
        super().__init__(n_genes)
        self.condition_emb = nn.Embedding(n_conditions, 64)
    
    def forward(self, x, edge_index, condition):
        # 加入条件信息
        cond_emb = self.condition_emb(condition)
        # ... 其余实现
```

---

## 💡 最佳实践

1. **始终保存原始数据**
   ```python
   adata.layers['counts'] = adata.X.copy()
   ```

2. **记录分析参数**
   ```python
   adata.uns['analysis_params'] = {
       'vae_epochs': 100,
       'vae_lr': 1e-3,
       'n_domains': 10,
       'date': '2025-01-30'
   }
   ```

3. **验证结果**
   ```python
   # 检查重构质量
   from sklearn.metrics import r2_score
   
   with torch.no_grad():
       recon, _, _, _ = predictor.vae(X, edge_index)
   
   r2 = r2_score(X.cpu().numpy(), recon.cpu().numpy())
   print(f"Reconstruction R²: {r2:.3f}")
   ```

4. **可重复性**
   ```python
   # 设置随机种子
   import random
   random.seed(42)
   np.random.seed(42)
   torch.manual_seed(42)
   if torch.cuda.is_available():
       torch.cuda.manual_seed_all(42)
   ```

---

**需要帮助?** 
- 查看完整文档: README.md
- 查看示例: examples/
- 提交Issue: GitHub Issues
