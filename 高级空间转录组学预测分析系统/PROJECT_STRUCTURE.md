# 📁 项目结构说明

```
spatial-transcriptomics-ai/
│
├── README.md                                    # 项目主文档
├── TUTORIAL.md                                  # 快速开始教程
├── requirements.txt                             # Python依赖
├── config.ini                                   # 配置文件
│
├── spatial_transcriptomics_prediction.py        # 🔥 主分析脚本
├── advanced_biological_analysis.py              # 🔬 高级下游分析
├── download_visium_data.py                      # 📥 数据下载工具
│
├── models/                                      # 深度学习模型
│   ├── __init__.py
│   ├── spatial_vae.py                          # 空间VAE模型
│   ├── diffusion_model.py                      # 扩散模型
│   ├── domain_discovery.py                     # 域发现模型
│   └── gnn_encoder.py                          # GNN编码器
│
├── utils/                                       # 工具函数
│   ├── __init__.py
│   ├── data_loader.py                          # 数据加载
│   ├── preprocessing.py                        # 预处理
│   ├── visualization.py                        # 可视化
│   └── metrics.py                              # 评估指标
│
├── analysis/                                    # 分析模块
│   ├── __init__.py
│   ├── differential_expression.py              # 差异表达
│   ├── pathway_enrichment.py                   # 通路富集
│   ├── spatial_trajectory.py                   # 空间轨迹
│   └── cell_communication.py                   # 细胞通讯
│
├── examples/                                    # 示例脚本
│   ├── basic_analysis.py                       # 基础分析示例
│   ├── custom_model.py                         # 自定义模型示例
│   ├── batch_processing.py                     # 批处理示例
│   └── integration_with_scrna.py               # 单细胞整合示例
│
├── notebooks/                                   # Jupyter笔记本
│   ├── 01_data_exploration.ipynb               # 数据探索
│   ├── 02_model_training.ipynb                 # 模型训练
│   ├── 03_prediction_analysis.ipynb            # 预测分析
│   └── 04_visualization.ipynb                  # 可视化
│
├── data/                                        # 数据目录
│   ├── breast_cancer/                          # 乳腺癌数据
│   │   ├── filtered_feature_bc_matrix.h5
│   │   ├── spatial/
│   │   └── processed_data.h5ad
│   └── mouse_brain/                            # 小鼠脑数据
│
├── outputs/                                     # 输出目录
│   ├── analyzed_visium_data.h5ad               # 分析结果
│   ├── predictions.npz                         # 预测数据
│   ├── figures/                                # 图片
│   │   ├── spatial_domains.png
│   │   ├── latent_space.png
│   │   └── predicted_genes.png
│   ├── tables/                                 # 表格
│   │   ├── DE_domain_0.csv
│   │   └── SVGs.csv
│   └── reports/                                # 报告
│       └── analysis_report.html
│
├── tests/                                       # 测试
│   ├── test_models.py
│   ├── test_preprocessing.py
│   └── test_analysis.py
│
└── docs/                                        # 文档
    ├── API.md                                   # API文档
    ├── METHODS.md                               # 方法说明
    └── FAQ.md                                   # 常见问题
```

---

## 📋 文件说明

### 核心脚本

#### `spatial_transcriptomics_prediction.py`
**主要功能分析脚本,包含:**
- `SpatialDataLoader`: 数据加载和预处理
- `SpatialGNNEncoder`: 图神经网络编码器
- `SpatialVAE`: 空间变分自编码器
- `SpatialDiffusionModel`: 扩散模型用于插值
- `SpatialDomainDiscovery`: 神经场域发现
- `SpatialTranscriptomicsPredictor`: 完整预测系统
- `SpatialVisualizer`: 可视化工具

**使用方法:**
```bash
python spatial_transcriptomics_prediction.py
```

#### `advanced_biological_analysis.py`
**高级下游生物学分析,包含:**
- `AdvancedBiologicalAnalysis`: 分析工具集
  - 差异表达分析
  - 通路富集分析
  - 空间轨迹推断
  - 细胞-细胞通讯
  - 空间自相关分析
  - 生态位聚类
  - 基因共表达网络

**使用方法:**
```bash
python advanced_biological_analysis.py
```

#### `download_visium_data.py`
**自动下载10x Genomics数据**

**使用方法:**
```bash
# 下载乳腺癌数据
python download_visium_data.py --dataset breast_cancer

# 下载所有文件
python download_visium_data.py --dataset breast_cancer --download-all

# 只加载数据不下载
python download_visium_data.py --dataset breast_cancer --skip-download
```

---

## 🔧 配置文件

### `config.ini`
包含所有可调参数,分为以下部分:
- `[data]`: 数据相关
- `[preprocessing]`: 预处理参数
- `[spatial_graph]`: 空间图构建
- `[model_vae]`: VAE模型
- `[model_diffusion]`: 扩散模型
- `[model_domain]`: 域发现模型
- `[training_*]`: 训练参数
- `[prediction]`: 预测参数
- `[analysis]`: 分析参数
- `[visualization]`: 可视化参数
- `[output]`: 输出设置
- `[advanced]`: 高级选项

**使用配置文件:**
```python
import configparser

config = configparser.ConfigParser()
config.read('config.ini')

# 读取参数
n_top_genes = config.getint('data', 'n_top_genes')
vae_epochs = config.getint('training_vae', 'epochs')
```

---

## 📊 输出文件

### AnnData对象结构 (`.h5ad`)

```python
adata
├── .X                          # 基因表达矩阵 (归一化后)
├── .layers
│   └── ['counts']             # 原始计数
├── .obs                        # Spot注释
│   ├── ['spatial_domain']     # 空间域分配
│   ├── ['niche']              # 生态位
│   ├── ['dpt_pseudotime']     # 伪时间
│   └── QC指标...
├── .var                        # 基因注释
│   ├── ['highly_variable']    # 高变基因标记
│   └── 统计信息...
├── .obsm                       # 多维注释
│   ├── ['spatial']            # 空间坐标 (n_spots, 2)
│   ├── ['latent']             # VAE嵌入 (n_spots, latent_dim)
│   ├── ['domain_probs']       # 域概率 (n_spots, n_domains)
│   └── ['X_diffmap']          # 扩散图
├── .obsp                       # Spot间关系
│   └── ['spatial_connectivities']  # 空间邻接矩阵
├── .uns                        # 非结构化数据
│   ├── ['spatial']            # 组织图像
│   ├── ['rank_genes_groups']  # 差异表达结果
│   ├── ['ligrec']             # 配体-受体分析
│   └── ['analysis_params']    # 分析参数
└── .varm / .varp              # 基因间关系
```

### 预测结果 (`predictions.npz`)

```python
predictions = np.load('predictions.npz')

predictions['coords']          # 预测位置坐标 (n_new, 2)
predictions['expression']      # 预测基因表达 (n_new, n_genes)
predictions['uncertainty']     # 预测不确定性 (可选)
```

### 差异表达结果 (`DE_*.csv`)

| gene | score | pval | pval_adj | logfoldchange |
|------|-------|------|----------|---------------|
| GENE1 | 5.2 | 1e-10 | 1e-8 | 2.3 |
| GENE2 | 4.8 | 1e-9 | 1e-7 | -1.9 |

### 空间可变基因 (`SVGs.csv`)

| gene | Morans_I | pvalue |
|------|----------|--------|
| SVG1 | 0.85 | 1e-15 |
| SVG2 | 0.78 | 1e-12 |

---

## 🔬 扩展开发

### 添加新模型

1. 在 `models/` 下创建新文件
2. 继承基类并实现必要方法
3. 在主脚本中导入使用

示例:
```python
# models/custom_model.py
from spatial_transcriptomics_prediction import SpatialVAE

class CustomSpatialModel(SpatialVAE):
    def __init__(self, n_genes, **kwargs):
        super().__init__(n_genes, **kwargs)
        # 添加自定义层
        self.custom_layer = nn.Linear(32, 64)
    
    def forward(self, x, edge_index):
        # 自定义前向传播
        pass
```

### 添加新分析

1. 在 `analysis/` 下创建新模块
2. 实现分析函数
3. 在高级分析脚本中集成

示例:
```python
# analysis/custom_analysis.py
def my_spatial_analysis(adata):
    """自定义空间分析"""
    # 实现分析逻辑
    return results
```

---

## 📚 相关资源

### 数据集
- [10x Genomics公开数据集](https://www.10xgenomics.com/resources/datasets)
- [Spatial Research](https://www.spatialresearch.org/)

### 工具
- [Scanpy](https://scanpy.readthedocs.io/)
- [Squidpy](https://squidpy.readthedocs.io/)
- [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/)

### 论文
- STAGATE (Nature Communications, 2022)
- SpaGCN (Nature Methods, 2021)
- Spatial transcriptomics review (Nature Reviews Genetics, 2023)

---

## 🤝 贡献

欢迎贡献! 请确保:
1. 代码符合PEP 8规范
2. 添加必要的文档字符串
3. 提供单元测试
4. 更新相关文档

---

## 📝 版本历史

### v1.0.0 (2025-01-30)
- 初始版本发布
- 实现GNN+VAE+扩散模型框架
- 完整的下游分析工具
- 详细文档和教程

---

## 📄 许可证

MIT License - 详见LICENSE文件
