# 快速参考指南 - 单细胞分析系统

## 🚀 一分钟上手

```r
# 1. 设置工作目录
setwd("/path/to/your/data")

# 2. 运行分析
source("Run_Complete_Analysis.R")
```

## 📝 关键配置（在Run_Complete_Analysis.R中修改）

```r
# ===== 必改参数 =====
WORK_DIR <- "/your/work/directory"                    # 工作目录
DATA_FILE <- "filtered_single_cell_data.rds"         # 数据文件

# ===== 运行选项 =====
RUN_MODULE1 <- TRUE      # 基础分析
RUN_MODULE2 <- TRUE      # 高级分析  
RESET_PROGRESS <- FALSE  # 断点续传 vs 重新开始

# ===== 细胞注释 =====
CLUSTER_ANNOTATION <- c(
  "0" = "T cells",
  "1" = "B cells",
  "2" = "Monocytes"
  # 根据聚类数量添加更多
)

# ===== 高级分析参数 =====
KNOCKOUT_TARGET <- "VCAN"     # 虚拟敲除的基因
GENE_OF_INTEREST <- "VCAN"    # 可视化的基因
```

## 📊 核心输出文件

| 文件位置 | 说明 |
|---------|------|
| `RData_Storage/Seurat.RData` | **最重要** - Seurat对象 |
| `Module1_Cell_Annotation/SingleR_Annotation.csv` | 自动注释参考 |
| `Module1_Marker_Genes/Top10_Markers_Per_Cluster.csv` | Marker基因 |
| `Differential_Analysis/Significant_Differential_Genes.csv` | 差异基因 |
| `Module2_CellChat_Analysis/` | 细胞通讯结果 |
| `Module2_Virtual_Knockout/` | 虚拟敲除结果 |

## 🔄 常用工作流

### 工作流1：完整新分析
```r
RESET_PROGRESS <- TRUE
RUN_MODULE1 <- TRUE
RUN_MODULE2 <- TRUE
source("Run_Complete_Analysis.R")
```

### 工作流2：只做基础分析
```r
RUN_MODULE1 <- TRUE
RUN_MODULE2 <- FALSE
source("Run_Complete_Analysis.R")
```

### 工作流3：更新细胞注释后重新运行
```r
# 1. 修改 CLUSTER_ANNOTATION
# 2. 运行
RESET_PROGRESS <- FALSE  # 保留其他步骤
source("Run_Complete_Analysis.R")
```

### 工作流4：只运行高级分析
```r
RUN_MODULE1 <- FALSE
RUN_MODULE2 <- TRUE
source("Run_Complete_Analysis.R")
```

### 工作流5：更换目标基因后重新分析
```r
# 1. 修改参数
KNOCKOUT_TARGET <- "NEW_GENE"
GENE_OF_INTEREST <- "NEW_GENE"

# 2. 重置模块2相关步骤
# 删除 complete_analysis_checkpoint.rds 中模块2的记录
# 或设置 RESET_PROGRESS <- TRUE

# 3. 运行
source("Run_Complete_Analysis.R")
```

## 🛠️ 常见操作

### 查看进度
```r
# R控制台中
tracker <- ProgressTracker$new()
tracker$print_progress()
```

### 重置特定模块
```r
# 删除检查点文件
file.remove("complete_analysis_checkpoint.rds")

# 然后重新运行
source("Run_Complete_Analysis.R")
```

### 加载已有结果
```r
# 加载Seurat对象
load("RData_Storage/Seurat.RData")

# pbmc 是Seurat对象
# cellAnn 是细胞注释

# 查看基本信息
print(pbmc)
table(pbmc$cell_type)
```

### 自定义单个分析
```r
# 只运行CellChat
source("Module2_Advanced_Analysis.R")
load("RData_Storage/Seurat.RData")
seurat_obj <- pbmc

# 设置参数
params <- list(
  cellchat_db = "Secreted Signaling",
  cellchat_min_cells = 10
)

output_dirs <- list(
  module2_cellchat = "Module2_CellChat_Analysis"
)

# 运行
result <- run_cellchat_analysis(seurat_obj, params, output_dirs)
```

## 🎯 参数调优建议

### 聚类分辨率
```r
# 数据预处理参数（在setup_complete_workspace中）
cluster_resolution = 0.4   # 默认
# 0.2-0.4: 粗聚类（大类）
# 0.5-0.8: 中等
# 0.9-1.5: 细聚类（亚型）
```

### 差异分析阈值
```r
logFCfilter = 1           # log2 fold change阈值
adjPvalFilter = 0.05      # 校正p值阈值
```

### CellChat设置
```r
cellchat_db = "Secreted Signaling"  # 分泌型信号
# 或
cellchat_db = "All"                 # 所有信号

cellchat_min_cells = 10    # 最小细胞数
```

### 虚拟敲除设置
```r
knockout_nc_nNet = 10       # 子网络数（10-50）
knockout_nc_nCells = 500    # 每网络细胞数（500-2000）
# 增大参数可提高准确性但增加计算时间
```

## 📦 快速安装依赖

```r
# 核心Bioconductor包
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "Seurat", "limma", "SingleR", "celldex",
  "org.Hs.eg.db", "clusterProfiler", 
  "ComplexHeatmap", "enrichplot"
))

# CRAN包
install.packages(c(
  "dplyr", "ggplot2", "patchwork", 
  "ggrepel", "R6", "progress"
))

# GitHub包
devtools::install_github("sqjin/CellChat")
devtools::install_github("immunogenomics/presto")
remotes::install_github("cailab-tamu/scTenifoldKnk")
```

## 🔍 快速诊断

### 问题：分析失败
```r
# 检查1：查看错误信息
# 终端会显示详细错误

# 检查2：验证数据
file.exists(DATA_FILE)
load("RData_Storage/Seurat.RData")
print(pbmc)

# 检查3：包是否安装
library(Seurat)
library(CellChat)
library(scTenifoldKnk)
```

### 问题：内存不足
```r
# 解决方案1：增加R内存
# 在启动R前设置：
# R --max-mem-size=32G

# 解决方案2：下采样
pbmc_sub <- subset(pbmc, downsample = 10000)

# 解决方案3：减少参数
knockout_nc_nCells = 300  # 降低
```

### 问题：找不到基因
```r
# 检查基因名是否存在
"VCAN" %in% rownames(pbmc)

# 查看所有基因名
head(rownames(pbmc), 20)

# 查找相似基因
grep("VCA", rownames(pbmc), value = TRUE)
```

## 📊 快速数据探索

```r
# 加载数据
load("RData_Storage/Seurat.RData")

# 基本统计
ncol(pbmc)                          # 细胞数
nrow(pbmc)                          # 基因数
table(pbmc$cell_type)               # 细胞类型分布
table(pbmc$Type)                    # 样本分组

# 快速可视化
DimPlot(pbmc, label = TRUE)         # UMAP图
VlnPlot(pbmc, features = "CD3D")    # 小提琴图
FeaturePlot(pbmc, features = "CD3D") # 特征图

# 导出特定数据
# 表达矩阵
expr <- GetAssayData(pbmc, slot = "data")
write.csv(as.matrix(expr[1:100, 1:100]), "sample_expression.csv")

# 元数据
write.csv(pbmc@meta.data, "metadata.csv")
```

## 💡 效率提示

1. **首次运行**：先用小数据集测试
2. **参数调整**：逐步优化，不要一次改太多
3. **保存关键点**：定期备份Seurat.RData
4. **并行计算**：可用的地方启用多核
5. **监控资源**：注意内存和CPU使用

## 🔗 相关资源

- [Seurat教程](https://satijalab.org/seurat/)
- [CellChat文档](https://github.com/sqjin/CellChat)
- [scTenifoldKnk文档](https://github.com/cailab-tamu/scTenifoldKnk)
- [clusterProfiler手册](https://yulab-smu.top/biomedical-knowledge-mining-book/)

---

**遇到问题？查看完整的 README.md 文档！**
