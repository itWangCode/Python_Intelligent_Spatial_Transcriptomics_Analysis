#!/usr/bin/env Rscript
# ==============================================================================
# 单细胞RNA-seq完整分析系统 v5.0 终极版
# Module 1 + Module 2 完整流程
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("单细胞RNA-seq分析系统 v5.0 终极完整版\n")
cat("Module 1: 核心分析 + Module 2: 高级分析\n")
cat(strrep("=", 80), "\n\n")

# ==============================================================================
# 步骤0: 环境准备
# ==============================================================================

cat("=== 步骤0: 环境准备 ===\n\n")

# 设置工作目录
work_dir <- getwd()
cat(sprintf("工作目录: %s\n", work_dir))

# 数据目录
data_dir <- file.path(work_dir, "data")
if (!dir.exists(data_dir)) {
  cat(sprintf("⚠ 数据目录不存在: %s\n", data_dir))
  cat("请创建数据目录并放入10X数据文件\n")
  stop("数据目录不存在")
}

# 创建输出目录结构
output_dirs <- list(
  # Module 1 目录
  module1_qc = file.path(work_dir, "Module1_Results", "01_QC"),
  module1_basic = file.path(work_dir, "Module1_Results", "02_Basic_Analysis"),
  module1_clustering = file.path(work_dir, "Module1_Results", "03_Clustering"),
  module1_markers = file.path(work_dir, "Module1_Results", "04_Marker_Genes"),
  module1_heatmaps = file.path(work_dir, "Module1_Results", "05_Heatmaps"),
  module1_annotation = file.path(work_dir, "Module1_Results", "06_Cell_Annotation"),
  module1_diff = file.path(work_dir, "Module1_Results", "07_Differential_Expression"),
  module1_proportion = file.path(work_dir, "Module1_Results", "08_Cell_Proportions"),
  
  # Module 2 目录
  module2_cellchat = file.path(work_dir, "Module2_Results", "01_CellChat"),
  module2_knockout = file.path(work_dir, "Module2_Results", "02_Virtual_Knockout"),
  module2_go = file.path(work_dir, "Module2_Results", "03_GO_Enrichment"),
  module2_gene_view = file.path(work_dir, "Module2_Results", "04_Gene_Visualization"),
  module2_advanced_plots = file.path(work_dir, "Module2_Results", "05_Advanced_Plots"),
  
  # 高级可视化
  advanced_viz = file.path(work_dir, "Advanced_Visualizations"),
  
  # 存储目录
  rdata_storage = file.path(work_dir, "RData_Storage"),
  reports = file.path(work_dir, "Statistical_Reports")
)

# 创建所有目录
for (dir_path in output_dirs) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

cat("✓ 输出目录结构已创建\n\n")

# ==============================================================================
# 步骤1: 加载R包
# ==============================================================================

cat("=== 步骤1: 加载必需的R包 ===\n\n")

required_packages <- c(
  # 核心包
  "Seurat", "dplyr", "ggplot2", "patchwork",
  
  # 可视化
  "RColorBrewer", "ggrepel", "ggsci", "viridis",
  
  # 数据处理
  "tidyr", "stringr", "readr",
  
  # 其他
  "Matrix", "scales"
)

# CRAN包安装
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("正在安装 %s...\n", pkg))
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# Bioconductor包
bioc_packages <- c("SingleR", "celldex", "org.Hs.eg.db", "clusterProfiler")

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("正在安装 %s (Bioconductor)...\n", pkg))
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg)
  }
}

cat("✓ 所有必需包已加载\n\n")

# ==============================================================================
# 步骤2: 加载分析模块
# ==============================================================================

cat("=== 步骤2: 加载分析模块 ===\n\n")

# Module 1 - 核心分析
module_files <- c(
  "Module1_Core_Analysis_ULTIMATE.R",
  "Module1_DiffAnalysis_Proportions.R",
  "Advanced_Visualizations.R",
  "Module2_Advanced_Analysis_ULTIMATE_Part1.R",
  "Module2_Advanced_Analysis_ULTIMATE_Part2.R"
)

for (module_file in module_files) {
  if (!file.exists(module_file)) {
    cat(sprintf("✗ 模块文件不存在: %s\n", module_file))
    stop(sprintf("缺少模块文件: %s", module_file))
  }
  source(module_file)
  cat(sprintf("✓ 已加载: %s\n", module_file))
}

cat("\n✓ 所有分析模块已加载\n\n")

# ==============================================================================
# 步骤3: 设置分析参数
# ==============================================================================

cat("=== 步骤3: 设置分析参数 ===\n\n")

params <- list(
  # ========== Module 1 参数 ==========
  # QC参数
  min_cells = 3,
  min_features = 200,
  max_features = 6000,
  mt_cutoff = 20,
  
  # 降维参数
  n_features = 2000,
  n_pcs = 30,
  
  # 聚类参数
  resolution = 0.8,
  
  # ========== Module 2 参数 ==========
  # CellChat参数
  run_cellchat = TRUE,
  cellchat_db = "All",  # 或 "Secreted Signaling"
  cellchat_min_cells = 10,
  
  # 虚拟敲除参数
  run_knockout = TRUE,
  knockout_target_gene = "VCAN",  # 修改为您的目标基因
  knockout_max_cells = 5000,
  knockout_qc_mtThreshold = 20,
  knockout_qc_minLSize = 200,
  knockout_nc_nNet = 10,
  knockout_nc_nCells = 500,
  
  # GO富集参数
  run_go_enrichment = TRUE,
  
  # 基因可视化参数
  run_gene_visualization = TRUE,
  gene_of_interest = "CD3D",  # 修改为您感兴趣的基因
  
  # 高级图表
  run_advanced_plots = TRUE
)

cat("分析参数:\n")
cat("Module 1:\n")
cat(sprintf("  特征范围: %d - %d\n", params$min_features, params$max_features))
cat(sprintf("  线粒体阈值: %d%%\n", params$mt_cutoff))
cat(sprintf("  聚类分辨率: %.1f\n", params$resolution))
cat("\nModule 2:\n")
cat(sprintf("  CellChat: %s\n", ifelse(params$run_cellchat, "启用", "禁用")))
cat(sprintf("  虚拟敲除: %s (目标: %s)\n",
            ifelse(params$run_knockout, "启用", "禁用"),
            params$knockout_target_gene))
cat(sprintf("  GO富集: %s\n", ifelse(params$run_go_enrichment, "启用", "禁用")))
cat(sprintf("  基因可视化: %s (基因: %s)\n",
            ifelse(params$run_gene_visualization, "启用", "禁用"),
            params$gene_of_interest))
cat("\n")

# ==============================================================================
# 步骤4: 运行Module 1 - 核心分析
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("开始Module 1: 核心分析流程\n")
cat(strrep("=", 80), "\n\n")

start_time <- Sys.time()

# 4.1 读取数据
cat("【1/9】读取10X数据...\n")
counts_data <- read_and_process_10x_data(data_dir, output_dirs)

# 4.2 质控
cat("\n【2/9】创建Seurat对象并质控...\n")
pbmc <- create_seurat_with_qc(counts_data, params, output_dirs)

# 4.3 归一化和降维
cat("\n【3/9】归一化和降维...\n")
pbmc <- normalize_and_reduce(pbmc, params, output_dirs)

# 4.4 Marker分析
cat("\n【4/9】Marker基因分析...\n")
all_markers <- find_all_markers(pbmc, params, output_dirs)

# 4.5 自动注释
cat("\n【5/9】自动细胞类型注释...\n")
pbmc <- automatic_cell_annotation(pbmc, output_dirs)

# 4.6 手动注释
cat("\n【6/9】应用手动注释...\n")
annotation_map <- c(
  "0" = "T cells",
  "1" = "B cells",
  "2" = "Monocytes",
  "3" = "NK cells",
  "4" = "Dendritic cells",
  "5" = "Macrophages"
)
pbmc <- manual_cell_annotation(pbmc, annotation_map, output_dirs)

# 4.7 差异分析
cat("\n【7/9】差异表达分析...\n")
pbmc <- perform_differential_expression(pbmc, params, output_dirs)

# 4.8 比例分析
cat("\n【8/9】细胞比例分析...\n")
analyze_proportions(pbmc, output_dirs)

# 4.9 保存结果
cat("\n【9/9】保存Module 1结果...\n")
pbmc <- save_module1_results(pbmc, output_dirs, work_dir)

module1_time <- Sys.time()
cat(sprintf("\n✓ Module 1 完成，用时: %.2f 分钟\n",
            difftime(module1_time, start_time, units = "mins")))

# ==============================================================================
# 步骤5: 运行高级可视化 (jjVolcano + Gene Cluster Heatmap)
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("开始高级可视化 (jjVolcano + Gene Cluster Heatmap)\n")
cat(strrep("=", 80), "\n\n")

advanced_results <- run_advanced_visualizations(
  work_dir = work_dir,
  data_dir = output_dirs$module1_diff
)

advanced_time <- Sys.time()
cat(sprintf("\n✓ 高级可视化完成，用时: %.2f 分钟\n",
            difftime(advanced_time, module1_time, units = "mins")))

# ==============================================================================
# 步骤6: 运行Module 2 - 高级分析
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("开始Module 2: 高级分析\n")
cat(strrep("=", 80), "\n\n")

module2_results <- run_all_module2_analyses(pbmc, params, output_dirs)

module2_time <- Sys.time()
cat(sprintf("\n✓ Module 2 完成，用时: %.2f 分钟\n",
            difftime(module2_time, advanced_time, units = "mins")))

# ==============================================================================
# 步骤7: 生成最终报告
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("生成最终分析报告\n")
cat(strrep("=", 80), "\n\n")

end_time <- Sys.time()
total_time <- difftime(end_time, start_time, units = "mins")

# 创建总结报告
report_file <- file.path(output_dirs$reports, "Complete_Analysis_Report.txt")

report_lines <- c(
  strrep("=", 80),
  "单细胞RNA-seq完整分析报告 v5.0",
  paste("分析日期:", Sys.Date()),
  paste("总运行时间:", round(total_time, 2), "分钟"),
  strrep("=", 80),
  "",
  "一、数据概览",
  strrep("-", 80),
  paste("  总细胞数:", ncol(pbmc)),
  paste("  总基因数:", nrow(pbmc)),
  paste("  聚类数:", length(unique(pbmc$seurat_clusters))),
  paste("  细胞类型数:", length(unique(pbmc$cell_type))),
  "",
  "二、Module 1 - 核心分析结果",
  strrep("-", 80),
  paste("  平均UMI数:", round(mean(pbmc$nCount_RNA), 2)),
  paste("  平均基因数:", round(mean(pbmc$nFeature_RNA), 2)),
  paste("  平均线粒体比例:", round(mean(pbmc$percent.mt), 2), "%"),
  paste("  显著marker数:", nrow(all_markers)),
  "",
  "三、Module 2 - 高级分析结果",
  strrep("-", 80),
  paste("  CellChat:", ifelse(!is.null(module2_results$cellchat), "✓ 完成", "✗ 跳过")),
  paste("  虚拟敲除:", ifelse(!is.null(module2_results$knockout), "✓ 完成", "✗ 跳过")),
  paste("  GO富集:", ifelse(!is.null(module2_results$go), "✓ 完成", "✗ 跳过")),
  paste("  基因可视化:", ifelse(!is.null(module2_results$gene_viz), "✓ 完成", "✗ 跳过")),
  paste("  高级图表:", ifelse(!is.null(module2_results$advanced_plots), "✓ 完成", "✗ 跳过")),
  "",
  "四、输出文件位置",
  strrep("-", 80),
  paste("  Module 1:", file.path(work_dir, "Module1_Results")),
  paste("  Module 2:", file.path(work_dir, "Module2_Results")),
  paste("  高级可视化:", output_dirs$advanced_viz),
  paste("  RData:", output_dirs$rdata_storage),
  paste("  报告:", output_dirs$reports),
  "",
  "五、运行时间分析",
  strrep("-", 80),
  paste("  Module 1:", round(difftime(module1_time, start_time, units = "mins"), 2), "分钟"),
  paste("  高级可视化:", round(difftime(advanced_time, module1_time, units = "mins"), 2), "分钟"),
  paste("  Module 2:", round(difftime(module2_time, advanced_time, units = "mins"), 2), "分钟"),
  paste("  总计:", round(total_time, 2), "分钟"),
  "",
  strrep("=", 80),
  "分析完成！",
  strrep("=", 80)
)

writeLines(report_lines, report_file)

# ==============================================================================
# 最终输出
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("🎉 完整分析全部完成！ 🎉\n")
cat(strrep("=", 80), "\n")
cat(sprintf("总运行时间: %.2f 分钟\n", total_time))
cat(sprintf("总细胞数: %d\n", ncol(pbmc)))
cat(sprintf("总基因数: %d\n", nrow(pbmc)))
cat(sprintf("聚类数: %d\n", length(unique(pbmc$seurat_clusters))))
cat("\n")
cat("主要输出目录:\n")
cat(sprintf("  Module 1: %s\n", file.path(work_dir, "Module1_Results")))
cat(sprintf("  Module 2: %s\n", file.path(work_dir, "Module2_Results")))
cat(sprintf("  高级可视化: %s\n", output_dirs$advanced_viz))
cat(sprintf("  完整报告: %s\n", report_file))
cat("\n")
cat("Module 1 包含:\n")
cat("  - QC质控 (4个图)\n")
cat("  - 降维聚类 (4个图)\n")
cat("  - Marker分析 (3个图)\n")
cat("  - 细胞注释 (3个图)\n")
cat("  - 差异分析 (3个图)\n")
cat("  - 细胞比例 (4个图)\n")
cat("\n")
cat("Module 2 包含:\n")
cat("  - CellChat细胞通讯分析\n")
cat("  - 虚拟基因敲除分析\n")
cat("  - GO富集分析\n")
cat("  - 基因可视化\n")
cat("  - 高级综合图表\n")
cat("\n")
cat("高级可视化包含:\n")
cat("  - jjVolcano多聚类火山图\n")
cat("  - Gene Cluster Enrichment Heatmap\n")
cat("\n")
cat(strrep("=", 80), "\n")
cat("请查看各目录中的PDF/PNG文件以查看结果\n")
cat(strrep("=", 80), "\n")
