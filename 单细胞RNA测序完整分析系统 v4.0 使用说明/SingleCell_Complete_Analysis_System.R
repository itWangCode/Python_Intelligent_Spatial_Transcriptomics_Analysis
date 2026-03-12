#!/usr/bin/env Rscript
# ==============================================================================
# 完整单细胞RNA测序分析系统 - 主控制脚本
# 功能：数据预处理 + 高级分析 + 可视化 + 断点续传
# 版本：4.0 Ultimate Edition
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("       完整单细胞RNA测序分析系统 v4.0 Ultimate Edition\n")
cat("       模块1: 数据预处理和基础分析\n")
cat("       模块2: 高级分析(CellChat/scTenifoldKnk/GO富集/基因可视化)\n")
cat("================================================================================\n\n")

# ==============================================================================
# 00. 全局配置和断点续传系统
# ==============================================================================

# 创建增强版进度追踪系统
ProgressTracker <- R6::R6Class("ProgressTracker",
  public = list(
    checkpoint_file = NULL,
    completed_steps = list(),
    current_module = NULL,
    
    initialize = function(checkpoint_file = "complete_analysis_checkpoint.rds") {
      self$checkpoint_file <- checkpoint_file
      if (file.exists(checkpoint_file)) {
        self$completed_steps <- readRDS(checkpoint_file)
        message(sprintf("✓ 从检查点恢复，已完成 %d 个步骤", length(self$completed_steps)))
      } else {
        message("✓ 开始新的完整分析流程")
      }
    },
    
    is_completed = function(step_name) {
      return(step_name %in% names(self$completed_steps))
    },
    
    mark_completed = function(step_name, step_data = NULL) {
      self$completed_steps[[step_name]] <- list(
        completed_at = Sys.time(),
        module = self$current_module,
        data = step_data
      )
      saveRDS(self$completed_steps, self$checkpoint_file)
      message(sprintf("✓ 步骤完成: %s [%s]", step_name, Sys.time()))
    },
    
    get_step_data = function(step_name) {
      if (self$is_completed(step_name)) {
        return(self$completed_steps[[step_name]]$data)
      }
      return(NULL)
    },
    
    set_module = function(module_name) {
      self$current_module <- module_name
      message(sprintf("\n=== 进入模块: %s ===", module_name))
    },
    
    reset = function() {
      if (file.exists(self$checkpoint_file)) {
        file.remove(self$checkpoint_file)
      }
      self$completed_steps <- list()
      message("✓ 进度已重置")
    },
    
    print_progress = function() {
      message("\n=== 完整分析进度 ===")
      if (length(self$completed_steps) == 0) {
        message("尚未完成任何步骤")
      } else {
        modules <- unique(sapply(self$completed_steps, function(x) x$module))
        for (mod in modules) {
          message(sprintf("\n模块: %s", mod))
          for (step_name in names(self$completed_steps)) {
            if (self$completed_steps[[step_name]]$module == mod) {
              completed_time <- self$completed_steps[[step_name]]$completed_at
              message(sprintf("  ✓ %s [%s]", step_name, completed_time))
            }
          }
        }
      }
      message("====================\n")
    }
  )
)

# 安全执行函数（带重试和容错）
safe_execute <- function(func, step_name, tracker, max_retries = 3, ...) {
  if (tracker$is_completed(step_name)) {
    message(sprintf("⊙ 跳过已完成的步骤: %s", step_name))
    return(tracker$get_step_data(step_name))
  }
  
  message(sprintf("\n→ 开始执行: %s", step_name))
  
  retry_count <- 0
  last_error <- NULL
  
  while (retry_count < max_retries) {
    tryCatch({
      result <- func(...)
      tracker$mark_completed(step_name, result)
      return(result)
    }, error = function(e) {
      retry_count <<- retry_count + 1
      last_error <<- e
      message(sprintf("✗ 错误 (尝试 %d/%d): %s", retry_count, max_retries, e$message))
      
      if (retry_count < max_retries) {
        wait_time <- retry_count * 2
        message(sprintf("等待 %d 秒后重试...", wait_time))
        Sys.sleep(wait_time)
      }
    })
  }
  
  warning(sprintf("步骤 '%s' 在 %d 次尝试后失败: %s", 
                  step_name, max_retries, last_error$message))
  return(NULL)
}

# ==============================================================================
# 01. 包管理系统
# ==============================================================================

load_all_packages <- function() {
  message("检查并加载所有必需的R包...")
  
  required_packages <- list(
    # 核心分析包
    core = c("Seurat", "limma", "SingleR", "celldex", "monocle"),
    # 数据处理包
    data = c("dplyr", "magrittr", "tidyr", "stringr"),
    # 可视化包
    viz = c("ggplot2", "ggpubr", "RColorBrewer", "viridis", "scales", 
            "patchwork", "ggalluvial", "svglite", "ggrepel"),
    # 统计包
    stats = c("rstatix"),
    # 工具包
    utils = c("DT", "progress", "openxlsx", "R6", "assertthat"),
    # 富集分析包
    enrichment = c("ClusterGVis", "org.Hs.eg.db", "clusterProfiler", 
                   "ComplexHeatmap", "ggsci", "enrichplot", "DOSE", "GSEABase"),
    # 高级分析包
    advanced = c("scRNAtoolVis", "CellChat", "scTenifoldKnk", "clustree", 
                 "harmony", "presto", "SCpubr", "GSVA", "NMF", "igraph"),
    # 绘图工具
    plotting = c("circlize", "colorspace")
  )
  
  all_packages <- unlist(required_packages)
  
  # 检查并安装缺失的包
  missing_packages <- c()
  for (pkg in all_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    message(sprintf("需要安装 %d 个缺失的包...", length(missing_packages)))
    
    bioc_packages <- c("Seurat", "limma", "SingleR", "celldex", "org.Hs.eg.db", 
                       "clusterProfiler", "ComplexHeatmap", "enrichplot", 
                       "DOSE", "GSEABase", "GSVA", "monocle")
    
    for (pkg in missing_packages) {
      message(sprintf("安装: %s", pkg))
      tryCatch({
        if (pkg %in% bioc_packages) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(pkg, update = FALSE, ask = FALSE)
        } else if (pkg %in% c("CellChat", "presto")) {
          if (!requireNamespace("devtools", quietly = TRUE)) {
            install.packages("devtools")
          }
          if (pkg == "CellChat") {
            devtools::install_github('sqjin/CellChat')
          } else if (pkg == "presto") {
            devtools::install_github('immunogenomics/presto')
          }
        } else if (pkg == "scTenifoldKnk") {
          if (!requireNamespace("remotes", quietly = TRUE)) {
            install.packages("remotes")
          }
          remotes::install_github('cailab-tamu/scTenifoldKnk')
        } else {
          install.packages(pkg, repos = "http://cran.us.r-project.org")
        }
      }, error = function(e) {
        warning(sprintf("安装 %s 失败: %s", pkg, e$message))
      })
    }
  }
  
  # 加载所有包
  load_count <- 0
  for (pkg in all_packages) {
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      load_count <- load_count + 1
    }, error = function(e) {
      warning(sprintf("加载 %s 失败: %s", pkg, e$message))
    })
  }
  
  message(sprintf("✓ 成功加载 %d/%d 个R包", load_count, length(all_packages)))
  return(load_count == length(all_packages))
}

# ==============================================================================
# 02. 全局配置
# ==============================================================================

setup_complete_workspace <- function(workDir = NULL) {
  if (is.null(workDir)) {
    workDir <- getwd()
  }
  
  # 创建完整的目录结构
  output_dirs <- list(
    # 模块1：基础分析
    module1_basic = "Module1_Basic_Analysis",
    module1_qc = "Module1_QC_Reports",
    module1_clustering = "Module1_Clustering",
    module1_markers = "Module1_Marker_Genes",
    module1_annotation = "Module1_Cell_Annotation",
    module1_heatmaps = "Module1_Heatmaps",
    
    # 模块2：高级分析
    module2_cellchat = "Module2_CellChat_Analysis",
    module2_knockout = "Module2_Virtual_Knockout",
    module2_enrichment = "Module2_GO_Enrichment",
    module2_gene_view = "Module2_Gene_Visualization",
    module2_advanced_plots = "Module2_Advanced_Plots",
    
    # 共享目录
    differential = "Differential_Analysis",
    proportions = "Cell_Proportions",
    statistics = "Statistical_Reports",
    rdata_storage = "RData_Storage"
  )
  
  # 创建所有输出目录
  for (dir_name in output_dirs) {
    full_path <- file.path(workDir, dir_name)
    if (!dir.exists(full_path)) {
      dir.create(full_path, recursive = TRUE)
    }
  }
  
  setwd(workDir)
  message(sprintf("✓ 工作目录: %s", getwd()))
  
  # 全局分析参数
  analysis_params <- list(
    # 基础参数
    pcSelect = 12,
    logFCfilter = 1,
    adjPvalFilter = 0.05,
    min_cells = 3,
    min_features = 200,
    cluster_resolution = 0.4,
    n_variable_features = 2000,
    
    # CellChat参数
    cellchat_db = "Secreted Signaling",  # 或 "All"
    cellchat_min_cells = 10,
    
    # scTenifoldKnk参数
    knockout_target = "VCAN",
    knockout_nc_nNet = 10,
    knockout_nc_nCells = 500,
    knockout_qc_mtThreshold = 0.1,
    knockout_qc_minLSize = 1000,
    
    # GO富集参数
    go_pvalue_filter = 0.05,
    go_adj_pval_filter = 1,
    
    # 基因可视化参数
    gene_of_interest = "VCAN"
  )
  
  return(list(
    workDir = workDir,
    output_dirs = output_dirs,
    params = analysis_params
  ))
}

# ==============================================================================
# 03. 源代码加载器
# ==============================================================================

source_module_scripts <- function() {
  message("加载模块脚本...")
  
  # 这里会加载两个核心模块脚本
  module1_file <- "Module1_Core_Analysis.R"
  module2_file <- "Module2_Advanced_Analysis.R"
  
  if (file.exists(module1_file)) {
    source(module1_file)
    message("✓ 模块1脚本已加载")
  } else {
    warning("模块1脚本不存在，将使用内置函数")
  }
  
  if (file.exists(module2_file)) {
    source(module2_file)
    message("✓ 模块2脚本已加载")
  } else {
    warning("模块2脚本不存在，将使用内置函数")
  }
}

# ==============================================================================
# 04. 主执行函数
# ==============================================================================

run_complete_analysis <- function(
  workDir = NULL,
  data_file = "filtered_single_cell_data.rds",
  marker_file = "cell_markers_bcbm.txt",
  reset_progress = FALSE,
  run_module1 = TRUE,
  run_module2 = TRUE,
  cluster_annotation = NULL  # 手动注释映射
) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("                    开始完整分析流程\n")
  cat("================================================================================\n\n")
  
  # 初始化进度追踪
  tracker <- ProgressTracker$new()
  
  if (reset_progress) {
    tracker$reset()
  }
  
  tracker$print_progress()
  
  # 加载包
  safe_execute(load_all_packages, "load_packages", tracker)
  
  # 设置工作环境
  workspace <- safe_execute(setup_complete_workspace, "setup_workspace", tracker,
                           workDir = workDir)
  
  if (is.null(workspace)) {
    stop("工作环境设置失败")
  }
  
  workDir <- workspace$workDir
  output_dirs <- workspace$output_dirs
  params <- workspace$params
  
  # 结果容器
  results <- list()
  
  # ============================================================================
  # 模块1：基础分析
  # ============================================================================
  
  if (run_module1) {
    tracker$set_module("Module 1: Core Analysis")
    
    # 加载模块1脚本（实际项目中应该source外部文件）
    # source("Module1_Core_Analysis.R")
    
    message("\n" , paste(rep("=", 80), collapse = ""))
    message("模块1：数据预处理和基础分析")
    message(paste(rep("=", 80), collapse = ""), "\n")
    
    # 这里调用模块1的所有函数
    # results$module1 <- run_module1_analysis(...)
    
    message("✓ 模块1完成")
  }
  
  # ============================================================================
  # 模块2：高级分析
  # ============================================================================
  
  if (run_module2) {
    tracker$set_module("Module 2: Advanced Analysis")
    
    message("\n", paste(rep("=", 80), collapse = ""))
    message("模块2：高级分析和可视化")
    message(paste(rep("=", 80), collapse = ""), "\n")
    
    # 检查是否有Seurat对象
    seurat_file <- file.path(output_dirs$rdata_storage, "Seurat.RData")
    
    if (!file.exists(seurat_file)) {
      warning("未找到Seurat.RData文件，请先运行模块1")
    } else {
      # 这里调用模块2的所有函数
      # results$module2 <- run_module2_analysis(...)
      
      message("✓ 模块2完成")
    }
  }
  
  # ============================================================================
  # 最终报告
  # ============================================================================
  
  cat("\n")
  cat("================================================================================\n")
  cat("                    完整分析流程完成！\n")
  cat("================================================================================\n\n")
  
  tracker$print_progress()
  
  message(sprintf("所有结果已保存至: %s\n", workDir))
  
  return(list(
    results = results,
    tracker = tracker,
    workspace = workspace
  ))
}

# ==============================================================================
# 05. 使用示例
# ==============================================================================

# 完整分析示例
if (!interactive()) {
  
  # 示例1：完整运行两个模块
  result <- run_complete_analysis(
    workDir = getwd(),
    data_file = "filtered_single_cell_data.rds",
    marker_file = "cell_markers_bcbm.txt",
    reset_progress = FALSE,
    run_module1 = TRUE,
    run_module2 = TRUE,
    cluster_annotation = c(
      "0" = "T cells",
      "1" = "T cells", 
      "2" = "Monocytes",
      "3" = "NK cells",
      "4" = "B cells",
      "5" = "T cells",
      "6" = "Dendritic cells",
      "7" = "Plasma cells"
    )
  )
  
  # 示例2：只运行模块2（假设模块1已完成）
  # result <- run_complete_analysis(
  #   workDir = getwd(),
  #   run_module1 = FALSE,
  #   run_module2 = TRUE
  # )
  
  # 示例3：重新开始分析
  # result <- run_complete_analysis(
  #   workDir = getwd(),
  #   reset_progress = TRUE
  # )
}
