# ============================================================================
# 单细胞RNA-seq分析系统 v5.0 - 完整使用指南
# ============================================================================

## 📋 目录

1. [系统概述](#1-系统概述)
2. [新增功能](#2-新增功能)
3. [快速开始](#3-快速开始)
4. [详细使用说明](#4-详细使用说明)
5. [输出文件说明](#5-输出文件说明)
6. [常见问题](#6-常见问题)
7. [高级用法](#7-高级用法)

---

## 1. 系统概述

### 1.1 版本信息
- **版本**: v5.0 终极完整版
- **发布日期**: 2024-12-24
- **主要特性**:
  - ✅ 修复所有 GetAssayData() 错误（slot → layer）
  - ✅ 内存优化（批处理、细胞数限制）
  - ✅ 完整断点续传功能
  - ✅ 20+ SCI级别可视化图表
  - ✅ jjVolcano火山图
  - ✅ Gene Cluster Enrichment Heatmap
  - ✅ 自动化错误处理

### 1.2 系统要求
- **R版本**: ≥ 4.0.0
- **内存**: 推荐 ≥ 16 GB
- **操作系统**: Windows / macOS / Linux
- **必需R包**:
  - Seurat (≥ 5.0.0)
  - dplyr, ggplot2, patchwork
  - SingleR, celldex
  - scRNAtoolVis, ClusterGVis

---

## 2. 新增功能

### 2.1 错误修复
```r
# ❌ 旧代码（已废弃）
expr_data <- GetAssayData(pbmc, assay = "RNA", slot = "scale.data")

# ✅ 新代码（正确）
expr_data <- GetAssayData(pbmc, assay = "RNA", layer = "scale.data")
```

### 2.2 内存优化
- **批处理**: Marker分析每次处理5个聚类
- **细胞限制**: max.cells.per.ident = 500
- **基因限制**: 热图最多显示50个基因
- **采样导出**: 表达矩阵限制5000细胞

### 2.3 新增可视化

#### 质控阶段（4个图）
1. `01_QC_ViolinPlot.pdf` - QC指标小提琴图
2. `02_QC_FeatureScatter.pdf` - 特征散点图
3. `03_QC_Feature_Density.pdf` - 特征密度图
4. `04_QC_Statistics_Comparison.pdf` - 质控前后对比

#### 降维阶段（4个图）
1. `01_VariableFeatures.pdf` - 高变基因
2. `02_PCA_Plot.pdf` - PCA可视化
3. `04_ElbowPlot.pdf` - PC选择图
4. `01_UMAP_Clusters.pdf/png` - UMAP聚类

#### Marker分析（3个图）
1. `01_Markers_Heatmap.pdf` - Marker热图
2. `02_Markers_DotPlot.pdf` - Dot plot
3. `03_Top_Markers_FeaturePlot.pdf` - 特征图

#### 差异分析（3个图）
1. `01_Volcano_Plot_Enhanced.pdf` - 增强版火山图
2. `02_MA_Plot.pdf` - MA图
3. `03_Top_DiffGenes_FeaturePlot.pdf` - 差异基因特征图

#### 细胞比例（4个图）
1. `01_Proportions_Stacked.pdf` - 堆叠柱状图
2. `02_Proportions_Grouped.pdf` - 分组柱状图
3. `03_Proportions_Bubble.pdf` - 气泡图
4. `04_Pie_Chart_*.pdf` - 饼图

#### 高级可视化（3个图）
1. `jjVolcano_Plot.pdf/png` - jjVolcano火山图
2. `Gene_Cluster_Enrichment_Heatmap.pdf` - 基因聚类富集热图
3. `Gene_Cluster_Heatmap_Simple.pdf` - 简版热图

---

## 3. 快速开始

### 3.1 准备工作

**步骤1**: 创建工作目录
```bash
mkdir SingleCell_Analysis
cd SingleCell_Analysis
```

**步骤2**: 准备数据
```bash
# 创建data目录
mkdir data

# 将10X数据放入data目录
# 目录结构应该是:
# data/
#   ├── sample1/
#   │   ├── barcodes.tsv.gz
#   │   ├── features.tsv.gz
#   │   └── matrix.mtx.gz
#   ├── sample2/
#   │   └── ...
```

**步骤3**: 复制分析脚本
```bash
# 将以下文件复制到工作目录:
# - Module1_Core_Analysis_ULTIMATE.R
# - Module1_DiffAnalysis_Proportions.R
# - Advanced_Visualizations.R
# - Run_Complete_Analysis.R
```

### 3.2 一键运行

**在R中运行**:
```r
# 设置工作目录
setwd("/path/to/SingleCell_Analysis")

# 运行完整分析
source("Run_Complete_Analysis.R")
```

**或者在命令行运行**:
```bash
cd /path/to/SingleCell_Analysis
Rscript Run_Complete_Analysis.R
```

### 3.3 预期运行时间
- 小数据集（<10,000细胞）: 15-30分钟
- 中等数据集（10,000-50,000细胞）: 30-90分钟
- 大数据集（>50,000细胞）: 1.5-3小时

---

## 4. 详细使用说明

### 4.1 修改参数

编辑 `Run_Complete_Analysis.R` 中的参数部分：

```r
params <- list(
  # QC参数
  min_cells = 3,           # 基因至少在多少细胞中表达
  min_features = 200,      # 细胞最少基因数
  max_features = 6000,     # 细胞最多基因数
  mt_cutoff = 20,          # 线粒体比例阈值(%)
  
  # 降维参数
  n_features = 2000,       # 高变基因数
  n_pcs = 30,              # 主成分数
  
  # 聚类参数
  resolution = 0.8         # 聚类分辨率(0.4-1.2)
)
```

**参数调整建议**:
- `resolution`: 
  - 更小(0.4-0.6)→ 更少聚类
  - 更大(1.0-1.5)→ 更多聚类
- `mt_cutoff`: 
  - 正常组织: 15-20%
  - 癌症/受损组织: 20-30%

### 4.2 修改细胞注释

编辑注释映射：

```r
annotation_map <- c(
  "0" = "CD4+ T cells",
  "1" = "CD8+ T cells", 
  "2" = "B cells",
  "3" = "NK cells",
  "4" = "Monocytes",
  "5" = "Macrophages",
  "6" = "Dendritic cells",
  "7" = "Neutrophils"
  # 添加更多聚类...
)
```

**如何确定注释**:
1. 查看自动注释结果: `01_Auto_Annotation_UMAP.pdf`
2. 查看Marker基因: `All_Markers.csv`
3. 参考文献中的Marker基因
4. 使用在线工具: CellMarker, PanglaoDB

### 4.3 断点续传

系统自动保存断点：
```r
# 如果分析中断，再次运行会自动从断点恢复
source("Run_Complete_Analysis.R")

# 系统会检测以下断点文件:
# - checkpoint_read_data.rds
# - checkpoint_seurat_qc.rds
# - checkpoint_normalized.rds
# - checkpoint_markers.rds
# - checkpoint_annotation.rds
# - checkpoint_diff_expr.rds
```

**清除断点重新分析**:
```r
# 删除所有断点文件
file.remove(list.files(pattern = "^checkpoint_.*\\.rds$"))
```

### 4.4 仅运行特定模块

如果只想运行某个模块：

```r
# 例如：只运行高级可视化
source("Advanced_Visualizations.R")

work_dir <- "/path/to/your/results"
data_dir <- file.path(work_dir, "Module1_Results/07_Differential_Expression")

advanced_results <- run_advanced_visualizations(
  work_dir = work_dir,
  data_dir = data_dir
)
```

---

## 5. 输出文件说明

### 5.1 目录结构

```
SingleCell_Analysis/
├── data/                                    # 输入数据
├── Module1_Results/                         # 核心分析结果
│   ├── 01_QC/                              # 质控图表
│   ├── 02_Basic_Analysis/                  # 基础分析
│   ├── 03_Clustering/                      # 聚类结果
│   ├── 04_Marker_Genes/                    # Marker基因
│   ├── 05_Heatmaps/                        # 热图
│   ├── 06_Cell_Annotation/                 # 细胞注释
│   ├── 07_Differential_Expression/         # 差异分析
│   └── 08_Cell_Proportions/                # 细胞比例
├── Advanced_Visualizations/                 # 高级可视化
│   ├── jjVolcano_Plot.pdf
│   ├── Gene_Cluster_Enrichment_Heatmap.pdf
│   └── ...
├── RData_Storage/                          # RData文件
│   ├── Seurat.RData
│   ├── Seurat_Object.rds
│   └── Normalized_Expression_Matrix.csv.gz
├── Statistical_Reports/                     # 统计报告
│   ├── Analysis_Summary_Report.txt
│   ├── Cell_Metadata.csv
│   └── Final_Analysis_Summary.csv
└── checkpoint_*.rds                        # 断点文件
```

### 5.2 重要文件说明

#### 核心数据文件
- **Seurat.RData**: 完整的Seurat对象（用于进一步分析）
- **Seurat_Object.rds**: Seurat对象的RDS格式
- **All_Markers.csv**: 所有Marker基因列表
- **Differential_Genes.csv**: 差异表达基因
- **Cell_Metadata.csv**: 细胞元数据

#### 关键可视化文件
- **01_UMAP_Clusters.png**: 最终聚类图（推荐用于发表）
- **03_Manual_Annotation_UMAP.png**: 注释后的UMAP（推荐用于发表）
- **01_Volcano_Plot_Enhanced.pdf**: 火山图（推荐用于发表）
- **jjVolcano_Plot.pdf**: 多聚类火山图
- **Gene_Cluster_Enrichment_Heatmap.pdf**: 基因聚类热图

---

## 6. 常见问题

### Q1: 安装包时出错怎么办？

```r
# 方法1: 更换CRAN镜像
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# 方法2: 使用清华镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# 方法3: 手动安装
install.packages("package_name", dependencies = TRUE)
```

### Q2: 内存不足怎么办？

```r
# 1. 增加内存限制（Windows）
memory.limit(size = 32000)  # 32 GB

# 2. 减少细胞数进行测试
pbmc <- subset(pbmc, cells = sample(colnames(pbmc), 5000))

# 3. 修改参数减少内存使用
params$n_features <- 1000   # 减少高变基因数
```

### Q3: GetAssayData错误？

**已修复！** v5.0已全部修复，确保使用最新代码。

### Q4: 没有生成某些图表？

检查以下内容：
```r
# 1. 查看是否有错误信息
cat("查看上方的⚠警告信息\n")

# 2. 检查是否有必要的分组信息
# 某些分析需要多个样本/组别
table(pbmc$orig.ident)

# 3. 手动运行单个函数进行调试
analyze_proportions(pbmc, output_dirs)
```

### Q5: 如何修改图表颜色？

编辑模块文件中的配色方案：

```r
# 在 Module1_Core_Analysis_ULTIMATE.R 中找到:
colors_sci <- colorRampPalette(c("#E8F4F8", "#A8D8E8", "#7BA7BC", 
                                 "#5B8FA4", "#D67280"))

# 修改为您喜欢的颜色:
my_colors <- c("#颜色1", "#颜色2", "#颜色3")
```

### Q6: jjVolcano或ClusterGVis安装失败？

```r
# jjVolcano需要scRNAtoolVis包
install.packages("scRNAtoolVis")

# ClusterGVis需要从GitHub安装
install.packages("remotes")
remotes::install_github("junjunlab/ClusterGVis")

# 如果GitHub访问困难，跳过高级可视化:
# 注释掉 Run_Complete_Analysis.R 中的这一行:
# advanced_results <- run_advanced_visualizations(...)
```

---

## 7. 高级用法

### 7.1 批量处理多个数据集

```r
# 创建批处理脚本
datasets <- c("dataset1", "dataset2", "dataset3")

for (dataset in datasets) {
  cat(sprintf("\n处理: %s\n", dataset))
  
  # 设置路径
  work_dir <- file.path("/path/to/results", dataset)
  data_dir <- file.path("/path/to/data", dataset)
  
  # 运行分析
  setwd(work_dir)
  source("Run_Complete_Analysis.R")
}
```

### 7.2 自定义可视化

```r
# 加载已保存的Seurat对象
pbmc <- readRDS("RData_Storage/Seurat_Object.rds")

# 自定义UMAP图
library(ggplot2)
p <- DimPlot(pbmc, reduction = "umap", 
             group.by = "cell_type",
             cols = c("red", "blue", "green")) +
  labs(title = "My Custom UMAP") +
  theme_minimal()

ggsave("My_Custom_UMAP.pdf", p, width = 10, height = 8)
```

### 7.3 导出到Scanpy (Python)

```r
# 1. 保存为h5ad格式
library(SeuratDisk)
SaveH5Seurat(pbmc, filename = "pbmc.h5Seurat")
Convert("pbmc.h5Seurat", dest = "h5ad")

# 2. 在Python中加载
# import scanpy as sc
# adata = sc.read_h5ad("pbmc.h5ad")
```

### 7.4 整合多个样本

```r
# 如果有多个样本需要整合
pbmc.list <- SplitObject(pbmc, split.by = "orig.ident")

# 使用Harmony整合
library(harmony)
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident")
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:30)
```

---

## 8. 引用和致谢

### 主要工具
- **Seurat**: Hao et al., Cell (2021)
- **SingleR**: Aran et al., Nat Immunol (2019)
- **ClusterGVis**: https://github.com/junjunlab/ClusterGVis
- **scRNAtoolVis**: https://github.com/junjunlab/scRNAtoolVis

### 支持
- 如有问题，请查看文档或联系开发者
- GitHub: [您的项目链接]
- Email: [您的邮箱]

---

## 9. 更新日志

### v5.0 (2024-12-24)
- ✅ 修复所有GetAssayData错误
- ✅ 新增内存优化
- ✅ 新增20+可视化图表
- ✅ 集成jjVolcano和Gene Cluster Heatmap
- ✅ 完整断点续传
- ✅ 自动错误处理

### v4.0 (之前版本)
- 基础分析流程
- Marker识别
- 细胞注释

---

## 10. 快速命令参考

```r
# 运行完整分析
source("Run_Complete_Analysis.R")

# 只运行高级可视化
source("Advanced_Visualizations.R")
run_advanced_visualizations(work_dir, data_dir)

# 加载保存的对象
pbmc <- readRDS("RData_Storage/Seurat_Object.rds")

# 查看细胞类型
table(pbmc$cell_type)

# 查看聚类
table(pbmc$seurat_clusters)

# 清除断点
file.remove(list.files(pattern = "^checkpoint_"))

# 生成自定义图
DimPlot(pbmc, group.by = "cell_type")
FeaturePlot(pbmc, features = c("CD3D", "CD19"))
```

---

**祝分析顺利！** 🎉
