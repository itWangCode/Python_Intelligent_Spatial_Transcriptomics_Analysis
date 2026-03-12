# 单细胞RNA测序完整分析系统 v4.0 使用说明

## 系统概述

这是一个完整的单细胞RNA测序数据分析系统，整合了从数据预处理到高级分析的全部功能。系统采用模块化设计，支持断点续传，可以灵活运行不同的分析模块。

### 核心功能

#### 模块1：数据预处理和基础分析
- 10X Genomics数据读取和合并
- 质量控制（QC）和过滤
- 数据归一化和标准化
- 降维分析（PCA、UMAP）
- 细胞聚类
- Marker基因识别
- 自动细胞类型注释（SingleR）
- 手动细胞类型注释
- 组间差异表达分析
- 细胞比例分析

#### 模块2：高级分析和可视化
- **CellChat** 细胞通讯分析
- **scTenifoldKnk** 虚拟基因敲除
- **GO富集分析**（包括圈图可视化）
- **基因表达可视化**
- **高级图表**（热图、气泡图等）

### 系统特点

**断点续传**：支持从中断处继续运行，节省时间
**容错处理**：自动重试失败步骤，提高稳定性
**模块化设计**：可独立运行任何模块
**自动化包管理**：自动检测和安装缺失的R包
**进度追踪**：实时显示分析进度
**全面兼容**：支持多种数据格式和版本

## 文件结构

```
SingleCell_Analysis_System/
├── SingleCell_Complete_Analysis_System.R    # 主控制脚本
├── Module1_Core_Analysis.R                  # 模块1：基础分析
├── Module2_Advanced_Analysis.R              # 模块2：高级分析
├── Run_Complete_Analysis.R                  # 执行脚本（推荐使用）
├── README.md                                # 本文档
└── 数据文件/
    ├── filtered_single_cell_data.rds        # （可选）预处理的RDS文件
    ├── cell_markers_bcbm.txt               # （可选）细胞标记基因列表
    └── 样本目录/                            # 10X Genomics原始数据
        ├── Sample1/
        │   ├── barcodes.tsv.gz
        │   ├── features.tsv.gz
        │   └── matrix.mtx.gz
        └── Sample2/
            └── ...
```

## 快速开始

### 步骤1：准备数据

**选项A：如果有10X Genomics原始数据**
```
工作目录/
├── Sample1/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── Sample2/
│   └── ...
└── Run_Complete_Analysis.R
```

**选项B：如果有预处理的RDS文件**
```
工作目录/
├── filtered_single_cell_data.rds
└── Run_Complete_Analysis.R
```

### 步骤2：配置参数

编辑 `Run_Complete_Analysis.R` 文件，修改以下参数：

```r
# 工作目录
WORK_DIR <- "/Users/wangyang/Desktop/SingleCell_Analysis_System" # 这个是我的地址❤️

# 数据文件
DATA_FILE <- "filtered_single_cell_data.rds"  # 如果有RDS文件
# 如果没有RDS，留空或注释掉，系统会自动读取10X数据

# 运行选项
RUN_MODULE1 <- TRUE   # 是否运行基础分析
RUN_MODULE2 <- TRUE   # 是否运行高级分析
RESET_PROGRESS <- FALSE  # 是否从头开始

# 细胞类型注释（在首次运行后根据结果修改）
CLUSTER_ANNOTATION <- c(
  "0" = "T cells",
  "1" = "T cells",
  "2" = "Monocytes",
  "3" = "NK cells",
  "4" = "B cells"
  # ... 根据实际聚类数量添加
)

# 高级分析参数
KNOCKOUT_TARGET <- "VCAN"      # 虚拟敲除的目标基因
GENE_OF_INTEREST <- "VCAN"     # 可视化的目标基因
```

### 步骤3：运行分析

**在RStudio中**：
1. 打开 `Run_Complete_Analysis.R`
2. 点击 "Source" 按钮
3. 等待分析完成

**在命令行中**：
```bash
source("Run_Complete_Analysis.R")
```

## 输出结果

### 模块1输出

#### 1. Module1_Basic_Analysis/
- Seurat对象和基础分析结果

#### 2. Module1_QC_Reports/
- `QC_Statistics.csv` - 质控统计
- `QC_Violin_Plots.pdf` - 质控小提琴图
- `QC_Scatter_Plots.pdf` - 质控散点图
- `Variable_Features.pdf` - 高变基因图

#### 3. Module1_Clustering/
- `UMAP_Basic.pdf` - 基础UMAP图
- `UMAP_Split_by_Group.pdf` - 分组UMAP图
- `PCA_Analysis.pdf` - PCA分析图
- `Cluster_Assignment.txt` - 聚类分配表

#### 4. Module1_Marker_Genes/
- `All_Markers.csv` - 所有marker基因
- `Significant_Markers.csv` - 显著marker基因
- `Top10_Markers_Per_Cluster.csv` - 每个聚类的Top10 marker
- `Top_Markers_Violin.pdf` - Marker基因小提琴图
- `Top_Markers_Feature.pdf` - Marker基因特征图

#### 5. Module1_Cell_Annotation/
- `SingleR_Annotation.csv` - 自动注释结果
- `Final_Cell_Type_Mapping.csv` - 最终细胞类型映射
- `Annotated_UMAP.pdf` - 注释后的UMAP图
- `Cluster_CellType_Table.csv` - 聚类-细胞类型对应表

#### 6. Module1_Heatmaps/
- `Marker_Heatmap.pdf` - Marker基因热图

#### 7. Differential_Analysis/
- `All_Differential_Genes.csv` - 所有差异基因
- `Significant_Differential_Genes.csv` - 显著差异基因
- `Volcano_Plot.pdf` - 火山图

#### 8. Cell_Proportions/
- `Cell_Proportions.csv` - 细胞比例数据
- `Cell_Proportions.pdf` - 细胞比例柱状图
- `Chi_Square_Test.csv` - 卡方检验结果

### 模块2输出

#### 9. Module2_CellChat_Analysis/
- `Communication_Network.txt` - 通讯网络数据
- `Network_Count.pdf` - 交互数量网络图
- `Network_Weight.pdf` - 交互强度网络图
- `Individual_Cell_Networks.pdf` - 单细胞类型网络图
- `Network_[CellType].pdf` - 各细胞类型独立网络图
- `Bubble_Plot.pdf` - 气泡图
- `Pathway_[PathwayName].pdf` - 信号通路图
- `CellChat_Object.rds` - CellChat对象

#### 10. Module2_Virtual_Knockout/
- `Significant_Diff_Genes.txt` - 敲除后显著差异基因
- `All_Diff_Genes.txt` - 所有差异基因
- `Barplot.pdf` - 柱状图
- `Volcano_Plot.pdf` - 火山图
- `Pie_Chart.pdf` - 饼图
- `scTenifoldKnk_Result.rds` - 分析结果对象

#### 11. Module2_GO_Enrichment/
- `GO_All_Results.txt` - 所有GO富集结果
- `GO_Significant_Results.txt` - 显著GO term
- `GO_Barplot.pdf` - GO柱状图
- `GO_Bubble.pdf` - GO气泡图
- `GO_Circlize.pdf` - GO圈图

#### 12. Module2_Gene_Visualization/
- `Gene_FeaturePlot.pdf` - 基因表达散点图
- `Gene_ViolinPlot.pdf` - 基因表达小提琴图
- `Gene_DotPlot.pdf` - 基因表达点图
- `Gene_FeaturePlot_Split.pdf` - 分组基因表达图
- `Gene_RidgePlot.pdf` - 基因表达岭线图

#### 13. Module2_Advanced_Plots/
- `CellType_Markers_DotPlot.pdf` - 细胞类型标记基因气泡图
- `Cell_Proportion_Stacked.pdf` - 细胞比例堆叠图
- `UMAP_Enhanced.pdf` - 美化版UMAP图

#### 14. RData_Storage/
- `Seurat.RData` - **核心Seurat对象**（包含pbmc和cellAnn）
- `Seurat_Object.rds` - Seurat对象RDS格式
- `Expression_Matrix.rds` - 表达矩阵
- `Metadata.csv` - 元数据

#### 15. Statistical_Reports/
- `Module1_Summary.csv` - 模块1统计摘要
- 其他统计报告

## 高级使用

### 断点续传

如果分析中断，只需重新运行脚本：

```r
# 不要重置进度
RESET_PROGRESS <- FALSE

# 重新运行
Rscript Run_Complete_Analysis.R
```

系统会自动跳过已完成的步骤，从中断处继续。

### 只运行特定模块

**只运行模块1（基础分析）**：
```r
RUN_MODULE1 <- TRUE
RUN_MODULE2 <- FALSE
```

**只运行模块2（高级分析）**：
```r
RUN_MODULE1 <- FALSE
RUN_MODULE2 <- TRUE
```

注意：模块2需要模块1的输出文件（Seurat.RData）。

### 重新开始分析

```r
RESET_PROGRESS <- TRUE
```

这会删除所有进度记录，从头开始分析。

### 自定义细胞类型注释

1. 首次运行后，查看以下文件：
   - `Module1_Cell_Annotation/SingleR_Annotation.csv` - 自动注释建议
   - `Module1_Marker_Genes/Top10_Markers_Per_Cluster.csv` - Marker基因

2. 根据结果修改 `CLUSTER_ANNOTATION`：
```r
CLUSTER_ANNOTATION <- c(
  "0" = "你的细胞类型1",
  "1" = "你的细胞类型2",
  "2" = "你的细胞类型3"
  # ...
)
```

3. 重新运行，系统会自动应用新的注释。

## 功能详解

### CellChat 细胞通讯分析

分析细胞类型之间的配体-受体相互作用。

**可调参数**：
```r
CELLCHAT_DB <- "Secreted Signaling"  # 或 "All"
CELLCHAT_MIN_CELLS <- 10
```

**输出解读**：
- **Network图**：显示细胞类型之间的通讯强度
- **Bubble图**：显示特定配体-受体对的表达
- **通路图**：显示信号通路的激活情况

### scTenifoldKnk 虚拟基因敲除

模拟特定基因被敲除后，基因调控网络的变化。

**可调参数**：
```r
KNOCKOUT_TARGET <- "VCAN"           # 目标基因
KNOCKOUT_NC_NNET <- 10              # 子网络数量
KNOCKOUT_NC_NCELLS <- 500           # 每个网络的细胞数
```

**输出解读**：
- **火山图**：显示敲除后的差异基因
- **柱状图**：Top差异基因
- **饼图**：显著基因比例

### GO富集分析

对差异基因进行功能富集分析。

**可调参数**：
```r
GO_PVALUE_FILTER <- 0.05
GO_ADJ_PVAL_FILTER <- 1
```

**输出解读**：
- **柱状图/气泡图**：显著富集的GO term
- **圈图**：GO term之间的关系和富集程度

## 依赖包列表

### 核心包
- Seurat
- limma
- SingleR
- celldex

### 数据处理
- dplyr
- magrittr
- tidyr
- stringr

### 可视化
- ggplot2
- ggpubr
- RColorBrewer
- viridis
- patchwork
- ggrepel

### 高级分析
- CellChat
- scTenifoldKnk
- clusterProfiler
- org.Hs.eg.db
- ComplexHeatmap
- circlize

所有包会在运行时自动检测和安装。

## 常见问题

### Q1: 如何确定聚类数量合适？

A: 查看以下文件：
- `Module1_Clustering/PCA_Analysis.pdf` 中的Elbow图
- 调整 `cluster_resolution` 参数（默认0.4）
- 较大的resolution会产生更多聚类

### Q2: 自动注释结果不准确怎么办？

A: 使用手动注释：
1. 查看 `Module1_Marker_Genes/Top10_Markers_Per_Cluster.csv`
2. 根据已知marker基因判断细胞类型
3. 修改 `CLUSTER_ANNOTATION` 参数
4. 重新运行

### Q3: 模块2某些分析失败了？

A: 可能原因：
- 缺少必要的R包
- 目标基因不存在
- 数据质量问题

解决方法：
- 检查错误信息
- 确认基因名称拼写正确
- 尝试安装缺失的包

### Q4: 如何处理大数据集？

A: 建议：
- 增加 `knockout_nc_nCells` 参数
- 减少细胞数量（下采样）
- 使用高性能计算集群

### Q5: 分析需要多长时间？

A: 取决于：
- 数据大小（细胞数、基因数）
- 计算资源
- 运行的模块

典型：
- 10,000细胞：模块1 ~30分钟，模块2 ~1-2小时
- 50,000细胞：模块1 ~1-2小时，模块2 ~3-5小时

## 引用

如果使用本系统发表文章，请引用以下软件：

- **Seurat**: Hao et al. (2021) Cell
- **SingleR**: Aran et al. (2019) Nature Immunology
- **CellChat**: Jin et al. (2021) Nature Communications
- **scTenifoldKnk**: Osorio et al. (2022) Patterns
- **clusterProfiler**: Yu et al. (2012) OMICS
