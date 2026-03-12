# tcga_complete_analysis.R 使用说明

## ✨ 主要特性

### 1. 🔄 断点续传功能
- 自动记录每个步骤的完成状态
- 重新运行时自动跳过已完成的步骤
- 进度保存在 `logs/analysis_progress.rds`
- 可以手动删除该文件以重新运行所有步骤

### 2. 🛡️ 错误处理机制
- 每个步骤都有错误捕获
- 失败的步骤会记录详细错误信息
- 支持降级策略（DESeq2失败→limma）
- 所有日志保存在 `logs/` 目录

### 3. 🎨 优化的可视化
- 热图：Top差异表达基因
- 火山图：SCI期刊风格，自动标注top基因
- PCA图：带置信椭圆
- MA图：展示表达量vs变化倍数
- 生存曲线：基于基因表达分组

### 4. 📊 自动化流程
- 自动识别肿瘤和正常样本
- 自动下载临床数据
- 自动生成统计报告
- 完整的目录结构

---

## 🚀 快速开始

### 第一次运行

```r
# 1. 修改工作目录（脚本第17行）
setwd("/Users/wangyang/Desktop/BCBM_TCGA_download")

# 2. 运行脚本
source("tcga_complete_analysis.R")

# 就这么简单！脚本会自动完成所有步骤
```

### 继续中断的分析

如果分析被中断（网络问题、电脑关机等），只需再次运行：

```r
source("tcga_complete_analysis.R")
```

脚本会自动从上次中断的地方继续！

### 强制重新运行

```r
# 删除进度文件
file.remove("logs/analysis_progress.rds")

# 重新运行
source("tcga_complete_analysis.R")
```

---

## 📁 输出文件结构

```
BCBM_TCGA_download/
├── results/
│   ├── figures/              # 所有图表
│   │   ├── heatmap_DEGs.pdf
│   │   ├── volcano_plot.pdf
│   │   ├── PCA_plot.pdf
│   │   ├── MA_plot.pdf
│   │   └── survival_plot.pdf
│   │
│   ├── tables/               # 所有数据表格
│   │   ├── differential_expression_all.csv
│   │   ├── differential_expression_significant.csv
│   │   ├── DEG_gene_list.txt
│   │   ├── clinical_data_complete.csv
│   │   ├── sample_info.csv
│   │   ├── gene_info.csv
│   │   └── clinical_summary.txt
│   │
│   ├── data/                 # 中间数据文件
│   │   ├── TCGA_RNAseq_raw.rds
│   │   ├── expr_filtered.rds
│   │   ├── grouped_data.rds
│   │   ├── deseq2_results.rds
│   │   └── pca_results.rds
│   │
│   └── analysis_report.txt   # 分析总结报告
│
├── logs/                     # 日志文件
│   ├── analysis_progress.rds # 进度记录
│   ├── analysis_YYYYMMDD_HHMMSS.log
│   └── analysis_parameters.rds
│
└── GDCdata/                  # 原始下载数据（自动创建）
```

---

## 🔧 自定义参数

### 方法1: 直接修改脚本

在脚本的第19-28行修改参数：

```r
PARAMS <- list(
  project = "TCGA-BRCA",        # 改为其他癌症类型
  threshold_logFC = 1,          # log2FC阈值
  threshold_adjP = 0.05,        # p值阈值
  max_display_genes = 50,       # 热图基因数
  min_sample_count = 10,        # 过滤参数
  min_count = 10
)
```

### 方法2: 使用配置文件（推荐）

```r
# 1. 复制配置模板
file.copy("config_template.R", "config.R")

# 2. 编辑 config.R 修改参数

# 3. 在脚本开头加载配置
source("config.R")
```

---

## 📊 分析步骤详解

### Step 1: 下载RNA-Seq数据
- 从GDC下载STAR-Counts数据
- 自动保存到本地
- **时间**: 1-3小时（取决于网速）

### Step 2: 准备数据
- 整理成SummarizedExperiment对象
- 保存为RDS文件
- **时间**: 5-15分钟

### Step 3: 数据预处理
- 提取表达矩阵
- 处理重复基因名
- 过滤低表达基因
- **时间**: 2-5分钟

### Step 4: 自动识别分组
- 从TCGA barcode识别样本类型
  - `01` = 肿瘤样本
  - `11` = 正常样本
- 生成样本信息表
- **时间**: <1分钟

### Step 5: 差异表达分析
- **优先使用**: DESeq2（更准确）
- **降级方案**: limma（更快速）
- 自动处理错误
- **时间**: 10-30分钟

### Step 6-10: 生成可视化
- 热图
- 火山图
- PCA图
- MA图
- **时间**: 每个1-2分钟

### Step 11: 生存分析
- 基于基因表达分组
- 自动选择top差异基因
- **时间**: 1-2分钟

### Step 12: 临床数据
- 下载完整临床信息
- 生成统计摘要
- **时间**: 1-5分钟

---

## 🎨 可视化说明

### 热图 (heatmap_DEGs.pdf)
- **内容**: Top差异表达基因
- **配色**: 红蓝配色（红=高表达，蓝=低表达）
- **注释**: 样本分组（肿瘤/正常）
- **用途**: 展示基因表达模式

### 火山图 (volcano_plot.pdf)
- **内容**: 所有基因的FC vs p值
- **配色**: 
  - 红色 = 上调基因
  - 蓝色 = 下调基因
  - 灰色 = 不显著
- **标注**: Top 20显著差异基因
- **用途**: 快速识别关键基因

### PCA图 (PCA_plot.pdf)
- **内容**: 样本在主成分空间的分布
- **特点**: 
  - 椭圆显示组内分布
  - 显示方差解释度
- **用途**: 评估样本分组是否明显

### MA图 (MA_plot.pdf)
- **内容**: 平均表达量 vs 变化倍数
- **用途**: 识别表达量偏倚

### 生存曲线 (survival_plot.pdf)
- **内容**: 高/低表达组的生存差异
- **特点**: 
  - 包含p值
  - 风险表
  - 置信区间
- **用途**: 评估基因预后价值

---

## 💡 常见问题

### Q1: 下载速度很慢
**A**: 
```r
# 修改下载参数（脚本第130行）
GDCdownload(
    query = query_rna,
    method = "api",
    files.per.chunk = 5  # 减少批量大小
)
```

### Q2: 内存不足
**A**:
```r
# 在脚本开头添加
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB
gc()  # 强制垃圾回收
```

### Q3: DESeq2分析失败
**A**: 脚本会自动切换到limma，无需手动操作

### Q4: 想分析其他癌症类型
**A**:
```r
# 修改PARAMS（第21行）
PARAMS <- list(
  project = "TCGA-LUAD",  # 改为肺腺癌
  ...
)
```

### Q5: 想重新运行某个步骤
**A**:
```r
# 方法1: 删除该步骤的输出文件
file.remove("results/figures/volcano_plot.pdf")

# 方法2: 编辑进度文件
progress <- readRDS("logs/analysis_progress.rds")
progress$step8_volcano <- FALSE
saveRDS(progress, "logs/analysis_progress.rds")

# 方法3: 删除整个进度文件（重新运行所有）
file.remove("logs/analysis_progress.rds")
```

### Q6: 如何查看详细日志
**A**:
```r
# 查看最新日志
log_files <- list.files("logs", pattern = "analysis_.*\\.log", full.names = TRUE)
latest_log <- log_files[which.max(file.info(log_files)$mtime)]
readLines(latest_log)
```

---

## 🔬 进阶使用

### 1. 批量分析多个项目

```r
projects <- c("TCGA-BRCA", "TCGA-LUAD", "TCGA-COAD")

for (proj in projects) {
  # 创建项目目录
  proj_dir <- paste0("results/", proj)
  dir.create(proj_dir, recursive = TRUE)
  
  # 修改工作目录和项目
  setwd(proj_dir)
  PARAMS$project <- proj
  
  # 运行分析
  source("../../tcga_complete_analysis.R")
}
```

### 2. 自定义可视化主题

在脚本中搜索 `theme_classic` 并修改：

```r
# 改为
theme_minimal() +
theme(
  plot.title = element_text(family = "Arial", face = "bold"),
  ...
)
```

### 3. 导出特定基因的表达数据

```r
# 在分析完成后
genes_of_interest <- c("TP53", "BRCA1", "ERBB2")

expr_subset <- expr_valid[genes_of_interest, ]
write.csv(expr_subset, "results/tables/selected_genes_expression.csv")
```

### 4. 批次效应校正

```r
# 如果怀疑有批次效应，在step5后添加：
library(sva)

# 找出批次信息
batch <- sample_info_valid$plate  # 假设plate是批次

# 使用ComBat校正
expr_combat <- ComBat(dat = expr_valid, batch = batch)
```

---

## 📈 性能优化

### 并行计算

```r
# 在脚本开头添加
library(BiocParallel)
register(MulticoreParam(workers = 4))  # 使用4核

# DESeq2会自动使用并行
dds <- DESeq(dds, parallel = TRUE)
```

### 减少内存使用

```r
# 定期清理
rm(large_object)
gc()

# 使用更紧凑的数据格式
expr_sparse <- as(expr_matrix, "dgCMatrix")
```

---

## 📚 输出文件说明

### 主要结果文件

1. **differential_expression_significant.csv**
   - 显著差异表达基因列表
   - 包含logFC、p值、调控方向等

2. **DEG_gene_list.txt**
   - 纯基因名列表
   - 可直接用于富集分析（DAVID、KEGG等）

3. **sample_info.csv**
   - 样本详细信息
   - 包含分组、临床特征等

4. **clinical_data_complete.csv**
   - 完整的临床数据
   - 可用于进一步的临床相关性分析

---

## 🎯 下一步分析建议

### 1. 功能富集分析

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# 读取差异基因
deg_genes <- read.table("results/tables/DEG_gene_list.txt")$V1

# GO富集
ego <- enrichGO(gene = deg_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# KEGG富集
kegg <- enrichKEGG(gene = bitr(deg_genes, 
                               fromType = "SYMBOL",
                               toType = "ENTREZID",
                               OrgDb = org.Hs.eg.db)$ENTREZID,
                   organism = "hsa")
```

### 2. 蛋白互作网络

```r
library(STRINGdb)

string_db <- STRINGdb$new(version = "11", species = 9606)
mapped <- string_db$map(deg_genes, "gene")
hits <- mapped$STRING_id[!is.na(mapped$STRING_id)]

# 绘制网络
string_db$plot_network(hits)
```

### 3. 基因集富集分析(GSEA)

```r
library(fgsea)

# 准备基因排序列表
res_df <- read.csv("results/tables/differential_expression_all.csv")
gene_list <- res_df$log2FoldChange
names(gene_list) <- res_df$gene
gene_list <- sort(gene_list, decreasing = TRUE)

# 运行GSEA
fgsea_results <- fgsea(pathways = hallmark_pathways,
                       stats = gene_list,
                       minSize = 15,
                       maxSize = 500)
```

---

## 🆘 获取帮助

1. **查看日志**: `logs/analysis_YYYYMMDD_HHMMSS.log`
2. **检查进度**: `readRDS("logs/analysis_progress.rds")`
3. **重新运行**: 删除进度文件
4. **Github Issues**: 提交问题到项目仓库

---

**祝分析顺利！** 🎉
