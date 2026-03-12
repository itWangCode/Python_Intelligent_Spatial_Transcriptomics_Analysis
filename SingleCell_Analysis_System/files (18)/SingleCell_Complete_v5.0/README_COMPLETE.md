# 单细胞RNA-seq完整分析系统 v5.0 - 终极完整版

## 🎯 这是什么？

这是一个**终极完整的单细胞RNA-seq分析系统**，包含Module 1（核心分析）和Module 2（高级分析）的所有功能。

## ✨ 完整功能列表

### Module 1: 核心分析（9个步骤）
1. ✅ 10X数据读取和合并
2. ✅ QC质控（4个图表）
3. ✅ 归一化和降维（4个图表）
4. ✅ Marker基因识别（3个图表）
5. ✅ 自动细胞注释（SingleR）
6. ✅ 手动细胞注释
7. ✅ 差异表达分析（3个图表）
8. ✅ 细胞比例分析（4个图表）
9. ✅ 结果保存和导出

### Module 2: 高级分析（5个模块）
1. ✅ **CellChat** - 细胞间通讯分析
   - 配体-受体相互作用
   - 信号通路网络
   - 和弦图、层级图
   - 气泡图

2. ✅ **Virtual Knockout** - 虚拟基因敲除
   - scTenifoldKnk分析
   - 敲除效应预测
   - 火山图、柱状图、饼图

3. ✅ **GO Enrichment** - GO富集分析
   - BP/CC/MF三大类
   - 柱状图、气泡图
   - GO网络图
   - Circos图

4. ✅ **Gene Visualization** - 基因可视化
   - FeaturePlot (散点图)
   - ViolinPlot (小提琴图)
   - DotPlot (点图)
   - RidgePlot (岭线图)
   - Heatmap (热图)

5. ✅ **Advanced Plots** - 高级综合图表
   - 细胞类型marker气泡图
   - 细胞比例堆叠图
   - UMAP密度图
   - 细胞周期分析
   - QC综合图
   - 细胞类型饼图

### 额外功能
- ✅ **jjVolcano** - 多聚类火山图
- ✅ **Gene Cluster Heatmap** - 基因聚类富集热图

## 📦 文件清单

### 核心模块（必需）
```
Module1_Core_Analysis_ULTIMATE.R              (24KB) - Module 1核心分析
Module1_DiffAnalysis_Proportions.R            (15KB) - 差异和比例分析
Module2_Advanced_Analysis_ULTIMATE_Part1.R    (28KB) - CellChat + 敲除 + GO
Module2_Advanced_Analysis_ULTIMATE_Part2.R    (16KB) - 基因可视化 + 高级图
Advanced_Visualizations.R                     (17KB) - jjVolcano + Cluster Heatmap
Run_COMPLETE_Analysis_with_Module2.R          (15KB) - 主运行脚本
```

### 文档（强烈推荐）
```
README_COMPLETE.md          - 完整说明（本文件）
USAGE_GUIDE.md              - 详细使用指南
QUICK_REFERENCE.txt         - 快速参考卡
```

**总计**: 8个文件，约115KB

## 🚀 快速开始（3步）

### 步骤1: 准备数据
```bash
mkdir SingleCell_Analysis
cd SingleCell_Analysis
mkdir data
# 将10X数据放入data/目录
```

### 步骤2: 复制所有文件到工作目录
```bash
# 复制所有.R文件
```

### 步骤3: 运行完整分析
```r
setwd("/path/to/SingleCell_Analysis")
source("Run_COMPLETE_Analysis_with_Module2.R")
```

## ⚙️ 参数配置

在`Run_COMPLETE_Analysis_with_Module2.R`中修改：

```r
params <- list(
  # ===== Module 1 参数 =====
  min_features = 200,
  max_features = 6000,
  mt_cutoff = 20,
  resolution = 0.8,
  
  # ===== Module 2 参数 =====
  # CellChat
  run_cellchat = TRUE,
  cellchat_db = "All",  # 或 "Secreted Signaling"
  
  # 虚拟敲除
  run_knockout = TRUE,
  knockout_target_gene = "VCAN",  # ← 修改为您的目标基因
  knockout_max_cells = 5000,
  
  # GO富集
  run_go_enrichment = TRUE,
  
  # 基因可视化
  run_gene_visualization = TRUE,
  gene_of_interest = "CD3D",  # ← 修改为感兴趣的基因
  
  # 高级图表
  run_advanced_plots = TRUE
)
```

## 📊 输出文件结构

```
SingleCell_Analysis/
├── Module1_Results/                    # 核心分析
│   ├── 01_QC/                         # 4个QC图
│   ├── 02_Basic_Analysis/             # 4个降维图
│   ├── 03_Clustering/                 # 2个UMAP图
│   ├── 04_Marker_Genes/               # 3个Marker图
│   ├── 05_Heatmaps/                   # 2个热图
│   ├── 06_Cell_Annotation/            # 3个注释图
│   ├── 07_Differential_Expression/    # 3个差异图
│   └── 08_Cell_Proportions/           # 4个比例图
│
├── Module2_Results/                    # 高级分析
│   ├── 01_CellChat/                   # 细胞通讯
│   │   ├── Network_Count.pdf
│   │   ├── Network_Weight.pdf
│   │   ├── Pathway_*.pdf
│   │   └── Communication_Network.txt
│   ├── 02_Virtual_Knockout/           # 虚拟敲除
│   │   ├── Volcano_Plot.pdf
│   │   ├── Barplot.pdf
│   │   ├── Pie_Chart.pdf
│   │   └── Significant_Diff_Genes.txt
│   ├── 03_GO_Enrichment/              # GO富集
│   │   ├── GO_Barplot.pdf
│   │   ├── GO_DotPlot.pdf
│   │   ├── GO_Network.pdf
│   │   ├── GO_Circos_Plot.pdf
│   │   └── GO_Enrichment_Results.txt
│   ├── 04_Gene_Visualization/         # 基因可视化
│   │   ├── Gene_FeaturePlot.pdf
│   │   ├── Gene_ViolinPlot.pdf
│   │   ├── Gene_DotPlot.pdf
│   │   ├── Gene_RidgePlot.pdf
│   │   └── Gene_Heatmap.pdf
│   └── 05_Advanced_Plots/             # 高级图表
│       ├── CellType_Markers_DotPlot.pdf
│       ├── Cell_Proportions_Stacked.pdf
│       ├── UMAP_Density.pdf
│       ├── Cell_Cycle_UMAP.pdf
│       ├── QC_Metrics_Combined.pdf
│       └── CellType_PieChart.pdf
│
├── Advanced_Visualizations/            # 额外可视化
│   ├── jjVolcano_Plot.pdf
│   ├── Gene_Cluster_Enrichment_Heatmap.pdf
│   └── Gene_Cluster_Heatmap_Simple.pdf
│
├── RData_Storage/                      # 数据存储
│   ├── Seurat.RData
│   ├── Seurat_Object.rds
│   └── Normalized_Expression_Matrix.csv.gz
│
└── Statistical_Reports/                # 统计报告
    └── Complete_Analysis_Report.txt
```

**总计**: 50+ PDF图表 + 15+ 数据文件

## ⏱️ 预期运行时间

| 数据规模 | Module 1 | Module 2 | 总计 |
|---------|---------|----------|------|
| <10K细胞 | 15-30分钟 | 20-40分钟 | 35-70分钟 |
| 10-50K细胞 | 30-90分钟 | 40-120分钟 | 70-210分钟 |
| >50K细胞 | 1.5-3小时 | 1-2小时 | 2.5-5小时 |

## 🔧 控制Module 2运行

如果您只想运行部分Module 2分析，可以修改参数：

```r
params <- list(
  # ... Module 1参数 ...
  
  # 关闭不需要的Module 2分析
  run_cellchat = FALSE,           # 关闭CellChat
  run_knockout = FALSE,           # 关闭虚拟敲除
  run_go_enrichment = TRUE,       # 保留GO富集
  run_gene_visualization = TRUE,  # 保留基因可视化
  run_advanced_plots = TRUE       # 保留高级图表
)
```

## 📚 Module 2 详细说明

### 2.1 CellChat - 细胞通讯分析
**功能**: 分析细胞类型间的配体-受体相互作用

**需要的包**:
```r
remotes::install_github("sqjin/CellChat")
```

**输出**:
- 交互数量和强度网络图
- 信号通路层级图和和弦图
- 配体-受体对气泡图
- 通讯网络数据表

**使用场景**: 研究细胞间如何通过信号分子进行通讯

---

### 2.2 Virtual Knockout - 虚拟基因敲除
**功能**: 预测敲除某个基因对其他基因的影响

**需要的包**:
```r
remotes::install_github("cailab-tamu/scTenifoldKnk")
```

**参数**:
```r
knockout_target_gene = "VCAN"    # 修改为您的目标基因
knockout_max_cells = 5000        # 限制细胞数以节省时间
```

**输出**:
- 差异基因火山图
- Top基因柱状图
- 显著基因比例饼图
- 差异基因列表

**使用场景**: 
- 研究关键基因的功能
- 预测基因敲除的下游效应
- 验证基因功能假说

---

### 2.3 GO富集分析
**功能**: 识别差异基因富集的生物学过程

**需要的包**: (已包含在Bioconductor)
```r
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
```

**输出**:
- GO柱状图
- GO气泡图
- GO关系网络图
- GO Circos图（环形可视化）

**使用场景**:
- 解释差异基因的生物学意义
- 发现富集的功能通路
- 文献阅读和假说生成

---

### 2.4 基因可视化
**功能**: 从多个角度展示单个基因的表达模式

**参数**:
```r
gene_of_interest = "CD3D"  # 修改为您感兴趣的基因
```

**输出**:
- FeaturePlot (UMAP上的表达)
- ViolinPlot (各细胞类型表达分布)
- DotPlot (跨细胞类型的表达模式)
- RidgePlot (表达密度图)
- Heatmap (表达热图)

**使用场景**:
- 验证marker基因
- 研究特定基因在不同细胞类型中的表达
- 准备发表图表

---

### 2.5 高级综合图表
**功能**: 生成publication-ready的综合图表

**输出**:
- 细胞类型marker基因气泡图（多基因）
- 细胞比例堆叠图（组间比较）
- UMAP密度图
- 细胞周期分析图
- QC指标综合图
- 细胞类型饼图

**使用场景**:
- 准备发表图表
- 展示数据质量
- 呈现细胞组成

## 💡 常见问题

### Q1: CellChat安装失败？
```r
# 方法1: 从GitHub安装
remotes::install_github("sqjin/CellChat")

# 方法2: 如果GitHub访问困难，设置run_cellchat = FALSE
```

### Q2: 虚拟敲除分析太慢？
```r
# 减少细胞数
knockout_max_cells = 2000

# 或减少网络数
knockout_nc_nNet = 5
```

### Q3: GO富集没有结果？
可能原因：
- 差异基因太少
- 基因ID转换失败
- 阈值太严格

解决：检查差异基因列表，确保有足够的显著基因

### Q4: 只想运行Module 1？
```r
# 使用原来的运行脚本
source("Run_Complete_Analysis.R")
```

### Q5: 如何跳过某个Module 2分析？
在参数中设置对应的开关为`FALSE`：
```r
params <- list(
  run_cellchat = FALSE,      # 跳过CellChat
  run_knockout = FALSE,      # 跳过虚拟敲除
  # ... 其他参数保持TRUE
)
```

## 📖 详细文档

1. **USAGE_GUIDE.md** - 完整使用教程
2. **QUICK_REFERENCE.txt** - 快速参考卡
3. **本文件** - 完整功能说明

## 🎯 使用建议

### 首次使用
1. 先运行小数据集测试
2. 只启用Module 1
3. 熟悉后再启用Module 2

### 日常使用
1. 根据需求选择Module 2功能
2. 调整参数优化运行时间
3. 利用断点续传功能

### 发表准备
1. Module 1: 基础QC和注释图
2. Module 2: 深入分析图表
3. 所有图表都是600 DPI，可直接用于发表

## 📊 系统要求

**最低配置**:
- R ≥ 4.0.0
- 内存 ≥ 16 GB (Module 2需要更多内存)
- 硬盘 ≥ 20 GB

**推荐配置**:
- R ≥ 4.3.0
- 内存 ≥ 32 GB
- 硬盘 ≥ 100 GB
- 多核CPU

## ✅ 已修复的所有问题

1. ✅ GetAssayData错误 (slot → layer) - 所有文件已修复
2. ✅ 内存限制 - 批处理和采样
3. ✅ 缺失图表 - 新增50+图表
4. ✅ SCI配色 - 柔和专业配色
5. ✅ 高级功能 - Module 2完整集成
6. ✅ 错误处理 - 优雅降级
7. ✅ 断点续传 - 自动保存恢复

## 🎉 核心特性

- ✅ **超级完整**: Module 1 + Module 2 全功能
- ✅ **50+ 图表**: 所有分析的可视化
- ✅ **Publication-ready**: 600 DPI，SCI配色
- ✅ **灵活控制**: 可选择启用的分析
- ✅ **智能优化**: 内存管理，批处理
- ✅ **详细文档**: 完整使用指南

---

**版本**: v5.0 终极完整版  
**发布日期**: 2024-12-24  
**包含**: Module 1 + Module 2 + 所有高级功能

祝您分析顺利！🎊
