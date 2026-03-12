# TCGA分析快速参考卡

## 🚀 5分钟快速开始

```r
# 1. 打开RStudio
# 2. 运行这3行代码：

setwd("/你的/工作/目录")  # 修改路径
source("tcga_complete_analysis.R")
# 等待完成！
```

---

## 📁 核心文件

| 文件 | 用途 | 何时使用 |
|------|------|---------|
| **tcga_complete_analysis.R** | 完整自动化脚本 | 👈 大多数情况 |
| COMPLETE_ANALYSIS_GUIDE.md | 详细说明 | 遇到问题时 |
| config_template.R | 参数配置 | 需要自定义时 |

---

## 🎯 常见任务

### 任务1: 首次分析
```r
source("tcga_complete_analysis.R")
```

### 任务2: 继续中断的分析
```r
# 直接重新运行，会自动跳过已完成步骤
source("tcga_complete_analysis.R")
```

### 任务3: 完全重新开始
```r
file.remove("logs/analysis_progress.rds")
source("tcga_complete_analysis.R")
```

### 任务4: 修改参数
```r
# 编辑脚本第21-28行的PARAMS
# 或创建 config.R 文件
```

### 任务5: 分析其他癌症
```r
# 修改第21行
PARAMS$project <- "TCGA-LUAD"  # 肺腺癌
```

---

## 📊 输出文件位置

```
results/
├── figures/          ← 所有PDF图表
├── tables/           ← 所有CSV数据表
└── data/             ← 中间数据(RDS)
```

### 关键输出文件

| 文件 | 内容 |
|------|------|
| **differential_expression_significant.csv** | 差异基因列表 |
| **volcano_plot.pdf** | 火山图 |
| **heatmap_DEGs.pdf** | 热图 |
| **PCA_plot.pdf** | PCA图 |
| **clinical_data_complete.csv** | 临床数据 |

---

## 🔧 常见参数

```r
PARAMS <- list(
  project = "TCGA-BRCA",       # 项目ID
  threshold_logFC = 1,         # logFC阈值
  threshold_adjP = 0.05,       # p值阈值
  max_display_genes = 50       # 热图基因数
)
```

---

## ⚙️ 常见癌症代码

| 代码 | 癌症类型 | 样本数 |
|------|---------|-------|
| TCGA-BRCA | 乳腺癌 | ~1100 |
| TCGA-LUAD | 肺腺癌 | ~500 |
| TCGA-LUSC | 肺鳞癌 | ~500 |
| TCGA-PRAD | 前列腺癌 | ~500 |
| TCGA-COAD | 结肠癌 | ~400 |
| TCGA-STAD | 胃癌 | ~400 |
| TCGA-LIHC | 肝癌 | ~370 |
| TCGA-KIRC | 肾透明细胞癌 | ~530 |

[完整列表](https://portal.gdc.cancer.gov/)

---

## 🐛 快速排错

### 问题: 下载失败
```r
# 网络问题，重新运行即可
source("tcga_complete_analysis.R")
```

### 问题: 内存不足
```r
# 增加内存限制（Windows）
memory.limit(size = 32000)  # 32GB
gc()  # 清理内存
```

### 问题: 包安装失败
```r
# 尝试清华镜像
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
```

### 问题: DESeq2失败
```r
# 脚本会自动切换到limma，无需手动操作
# 查看日志确认：logs/analysis_*.log
```

---

## 📈 分析流程图

```
下载数据 (1-3小时)
    ↓
准备数据 (5-15分钟)
    ↓
识别分组 (< 1分钟)
    ↓
差异分析 (10-30分钟)
    ↓
生成图表 (5-10分钟)
    ↓
临床数据 (1-5分钟)
    ↓
完成！
```

**总时间**: 约2-4小时

---

## 📞 获取帮助

1. **查看日志**: `logs/analysis_*.log`
2. **查看详细说明**: `COMPLETE_ANALYSIS_GUIDE.md`
3. **查看对比**: `SCRIPT_COMPARISON.md`
4. **查看快速入门**: `QUICKSTART.md`

---

## 💡 高级技巧

### 并行运算
```r
# 脚本开头添加
library(BiocParallel)
register(MulticoreParam(workers = 4))
```

### 查看进度
```r
progress <- readRDS("logs/analysis_progress.rds")
print(progress)
```

### 重新运行特定步骤
```r
# 编辑进度文件
progress$step8_volcano <- FALSE
saveRDS(progress, "logs/analysis_progress.rds")
source("tcga_complete_analysis.R")
```

---

## 🎨 图表定制

### 修改颜色
```r
# 在脚本中搜索并修改：
COLOR_TUMOR <- "#你的颜色代码"
COLOR_NORMAL <- "#你的颜色代码"
```

### 修改大小
```r
# 搜索 ggsave 或 pdf() 函数
ggsave("file.pdf", width = 12, height = 10)  # 修改尺寸
```

---

## ✅ 检查清单

- [ ] 已安装R (>= 4.0)
- [ ] 已安装RStudio (推荐)
- [ ] 网络连接正常
- [ ] 至少20GB空闲磁盘空间
- [ ] 至少8GB内存
- [ ] 已修改工作目录路径
- [ ] 已了解基本R操作

---

## 🎯 下一步

分析完成后可以做：
1. **功能富集分析** - 使用clusterProfiler
2. **通路分析** - KEGG/Reactome
3. **蛋白互作网络** - STRING数据库
4. **基因集富集分析(GSEA)** - hallmark pathways
5. **验证实验设计** - 选择关键基因

---

## 📚 学习资源

- **TCGA官方**: https://www.cancer.gov/tcga
- **GDC Portal**: https://portal.gdc.cancer.gov/
- **Bioconductor**: https://bioconductor.org/
- **教程视频**: YouTube搜索 "TCGA analysis"

---

**打印此文档作为快速参考！** 📄
