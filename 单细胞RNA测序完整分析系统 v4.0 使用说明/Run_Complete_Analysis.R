#!/usr/bin/env Rscript
# ==============================================================================
# 单细胞RNA测序完整分析系统 - 最终执行脚本
# 版本: 4.0 Ultimate Edition
# ==============================================================================

# ==============================================================================
# 源代码加载
# ==============================================================================

# 获取脚本所在目录
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0 || script_dir == "") {
  script_dir <- getwd()
}

# 加载主控制脚本
if (file.exists(file.path(script_dir, "SingleCell_Complete_Analysis_System.R"))) {
  source(file.path(script_dir, "SingleCell_Complete_Analysis_System.R"))
} else if (file.exists("SingleCell_Complete_Analysis_System.R")) {
  source("SingleCell_Complete_Analysis_System.R")
} else {
  stop("找不到主控制脚本文件")
}

# 加载模块1脚本
if (file.exists(file.path(script_dir, "Module1_Core_Analysis.R"))) {
  source(file.path(script_dir, "Module1_Core_Analysis.R"))
} else if (file.exists("Module1_Core_Analysis.R")) {
  source("Module1_Core_Analysis.R")
} else {
  warning("找不到模块1脚本，将使用默认功能")
}

# 加载模块2脚本
if (file.exists(file.path(script_dir, "Module2_Advanced_Analysis.R"))) {
  source(file.path(script_dir, "Module2_Advanced_Analysis.R"))
} else if (file.exists("Module2_Advanced_Analysis.R")) {
  source("Module2_Advanced_Analysis.R")
} else {
  warning("找不到模块2脚本，将使用默认功能")
}

# ==============================================================================
# 全局参数配置
# ==============================================================================

# 工作目录设置（请根据实际情况修改）
WORK_DIR <- getwd()  # 默认使用当前目录
# WORK_DIR <- "/path/to/your/work/directory"  # 或指定具体路径

# 数据文件设置
DATA_FILE <- "filtered_single_cell_data.rds"  # 如果有预处理的RDS文件
# 如果没有RDS文件，将自动从子目录读取10X数据

# 细胞标记基因文件（可选）
MARKER_FILE <- "cell_markers_bcbm.txt"

# 是否重置进度（TRUE = 从头开始，FALSE = 断点续传）
RESET_PROGRESS <- FALSE

# 运行哪些模块
RUN_MODULE1 <- TRUE   # 基础分析
RUN_MODULE2 <- TRUE   # 高级分析

# ==============================================================================
# 手动细胞类型注释（请根据实际聚类结果修改）
# ==============================================================================

# 方式1：根据SingleR自动注释的结果来定义
# 运行模块1后，查看 Module1_Cell_Annotation/SingleR_Annotation.csv
# 然后根据聚类-细胞类型对应关系来定义

# 方式2：根据marker基因手动定义
# 查看 Module1_Marker_Genes/Top10_Markers_Per_Cluster.csv
# 根据已知marker基因判断细胞类型

CLUSTER_ANNOTATION <- c(
  "0" = "T cells",
  "1" = "T cells",
  "2" = "Monocytes",
  "3" = "NK cells",
  "4" = "B cells",
  "5" = "T cells",
  "6" = "Dendritic cells",
  "7" = "Plasma cells"
)

# 注意：
# 1. 聚类编号（"0", "1", "2"...）需要与实际数据匹配
# 2. 细胞类型名称可以自定义
# 3. 如果不确定，可以先设置为NULL，使用自动注释

# ==============================================================================
# 高级分析参数
# ==============================================================================

# CellChat参数
CELLCHAT_DB <- "Secreted Signaling"  # 或 "All"
CELLCHAT_MIN_CELLS <- 10

# 虚拟敲除参数
KNOCKOUT_TARGET <- "VCAN"  # 目标基因
KNOCKOUT_NC_NNET <- 10
KNOCKOUT_NC_NCELLS <- 500
KNOCKOUT_QC_MTTHRESHOLD <- 0.1
KNOCKOUT_QC_MINLSIZE <- 1000

# GO富集参数
GO_PVALUE_FILTER <- 0.05
GO_ADJ_PVAL_FILTER <- 1

# 基因可视化参数
GENE_OF_INTEREST <- "VCAN"

# ==============================================================================
# 主执行流程
# ==============================================================================

main_execution <- function() {
  
  cat("\n")
  cat("================================================================================\n")
  cat("       单细胞RNA测序完整分析系统 v4.0\n")
  cat("================================================================================\n\n")
  
  # 打印配置信息
  cat("配置信息:\n")
  cat(sprintf("  工作目录: %s\n", WORK_DIR))
  cat(sprintf("  数据文件: %s\n", DATA_FILE))
  cat(sprintf("  运行模块1: %s\n", RUN_MODULE1))
  cat(sprintf("  运行模块2: %s\n", RUN_MODULE2))
  cat(sprintf("  重置进度: %s\n", RESET_PROGRESS))
  cat("\n")
  
  # 创建进度追踪器
  tracker <- ProgressTracker$new()
  
  if (RESET_PROGRESS) {
    tracker$reset()
    cat("✓ 进度已重置，将从头开始分析\n\n")
  } else {
    tracker$print_progress()
  }
  
  # 加载所有包
  cat("正在加载R包...\n")
  safe_execute(load_all_packages, "load_packages", tracker)
  
  # 设置工作环境
  cat("\n正在设置工作环境...\n")
  workspace <- safe_execute(setup_complete_workspace, "setup_workspace",
                           tracker, workDir = WORK_DIR)
  
  if (is.null(workspace)) {
    stop("❌ 工作环境设置失败")
  }
  
  workDir <- workspace$workDir
  output_dirs <- workspace$output_dirs
  params <- workspace$params
  
  # 更新参数
  params$cellchat_db <- CELLCHAT_DB
  params$cellchat_min_cells <- CELLCHAT_MIN_CELLS
  params$knockout_target <- KNOCKOUT_TARGET
  params$knockout_nc_nNet <- KNOCKOUT_NC_NNET
  params$knockout_nc_nCells <- KNOCKOUT_NC_NCELLS
  params$knockout_qc_mtThreshold <- KNOCKOUT_QC_MTTHRESHOLD
  params$knockout_qc_minLSize <- KNOCKOUT_QC_MINLSIZE
  params$go_pvalue_filter <- GO_PVALUE_FILTER
  params$go_adj_pval_filter <- GO_ADJ_PVAL_FILTER
  params$gene_of_interest <- GENE_OF_INTEREST
  
  # 结果容器
  results <- list()
  
  # ============================================================================
  # 模块1：基础分析
  # ============================================================================
  
  if (RUN_MODULE1) {
    cat("\n")
    cat("================================================================================\n")
    cat("                     模块1：数据预处理和基础分析\n")
    cat("================================================================================\n\n")
    
    results$module1 <- run_module1_analysis(
      workDir = workDir,
      data_file = DATA_FILE,
      marker_file = MARKER_FILE,
      params = params,
      output_dirs = output_dirs,
      tracker = tracker,
      cluster_annotation = CLUSTER_ANNOTATION
    )
    
    cat("\n")
    cat("================================================================================\n")
    cat("                     ✓ 模块1完成\n")
    cat("================================================================================\n\n")
  }
  
  # ============================================================================
  # 模块2：高级分析
  # ============================================================================
  
  if (RUN_MODULE2) {
    cat("\n")
    cat("================================================================================\n")
    cat("                     模块2：高级分析和可视化\n")
    cat("================================================================================\n\n")
    
    # 检查Seurat.RData是否存在
    seurat_file <- file.path(output_dirs$rdata_storage, "Seurat.RData")
    
    if (!file.exists(seurat_file)) {
      cat("❌ 未找到Seurat.RData文件\n")
      cat("   请先运行模块1或确保文件存在于: %s\n", seurat_file)
    } else {
      results$module2 <- run_module2_analysis(
        workDir = workDir,
        params = params,
        output_dirs = output_dirs,
        tracker = tracker
      )
      
      cat("\n")
      cat("================================================================================\n")
      cat("                     ✓ 模块2完成\n")
      cat("================================================================================\n\n")
    }
  }
  
  # ============================================================================
  # 最终总结
  # ============================================================================
  
  cat("\n")
  cat("================================================================================\n")
  cat("                     🎉 完整分析流程已完成！ 🎉\n")
  cat("================================================================================\n\n")
  
  tracker$print_progress()
  
  cat("结果文件位置:\n")
  cat(sprintf("  工作目录: %s\n", workDir))
  cat("\n模块1结果:\n")
  cat(sprintf("  - 基础分析: %s\n", output_dirs$module1_basic))
  cat(sprintf("  - 聚类结果: %s\n", output_dirs$module1_clustering))
  cat(sprintf("  - Marker基因: %s\n", output_dirs$module1_markers))
  cat(sprintf("  - 细胞注释: %s\n", output_dirs$module1_annotation))
  
  if (RUN_MODULE2) {
    cat("\n模块2结果:\n")
    cat(sprintf("  - CellChat: %s\n", output_dirs$module2_cellchat))
    cat(sprintf("  - 虚拟敲除: %s\n", output_dirs$module2_knockout))
    cat(sprintf("  - GO富集: %s\n", output_dirs$module2_enrichment))
    cat(sprintf("  - 基因可视化: %s\n", output_dirs$module2_gene_view))
  }
  
  cat("\n重要文件:\n")
  cat(sprintf("  - Seurat对象: %s/Seurat.RData\n", output_dirs$rdata_storage))
  cat(sprintf("  - 表达矩阵: %s/Expression_Matrix.rds\n", output_dirs$rdata_storage))
  cat(sprintf("  - 元数据: %s/Metadata.csv\n", output_dirs$rdata_storage))
  
  cat("\n")
  cat("================================================================================\n\n")
  
  return(list(
    results = results,
    tracker = tracker,
    workspace = workspace
  ))
}

# ==============================================================================
# 执行分析
# ==============================================================================

if (!interactive()) {
  # 捕获可能的错误
  result <- tryCatch({
    main_execution()
  }, error = function(e) {
    cat("\n")
    cat("================================================================================\n")
    cat("                     ❌ 分析过程中出现错误\n")
    cat("================================================================================\n\n")
    cat(sprintf("错误信息: %s\n", e$message))
    cat("\n")
    cat("请检查:\n")
    cat("  1. 数据文件是否存在且格式正确\n")
    cat("  2. 必要的R包是否已安装\n")
    cat("  3. 工作目录是否有写入权限\n")
    cat("  4. 参数设置是否正确\n")
    cat("\n")
    cat("如需帮助，请查看完整的错误堆栈:\n")
    print(e)
    cat("\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    cat("✓ 分析成功完成！\n")
    cat("✓ 所有结果已保存\n")
    cat("✓ 可以查看各个目录下的图表和数据文件\n\n")
  }
}
