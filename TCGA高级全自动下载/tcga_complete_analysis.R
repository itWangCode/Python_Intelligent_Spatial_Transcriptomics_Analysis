#!/usr/bin/env Rscript
# ============================================================================
# TCGA完整自动化分析流程
# 功能：断点续传 + 错误处理 + 自动分组识别 + 精美可视化
# 作者: itwangyang
# 日期: 2025-11-27
# ============================================================================

# ========================== 全局配置 ==========================
cat("="*80, "\n")
cat("TCGA完整自动化分析流程\n")
cat("="*80, "\n\n")

# 工作目录设置
WORK_DIR <- "/Users/wangyang/Desktop/BCBM_TCGA_download"
setwd(WORK_DIR)

# 创建输出目录结构
dir.create("results", showWarnings = FALSE)
dir.create("results/figures", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
dir.create("results/data", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# 进度记录文件
PROGRESS_FILE <- "logs/analysis_progress.rds"

# 分析参数
PARAMS <- list(
  project = "TCGA-BRCA",
  threshold_logFC = 1,
  threshold_adjP = 0.05,
  max_display_genes = 50,
  min_sample_count = 10,  # 过滤基因：至少在N个样本中表达
  min_count = 10           # 过滤基因：最小count值
)

# 保存参数
saveRDS(PARAMS, "logs/analysis_parameters.rds")

# ========================== 日志和进度管理 ==========================

# 初始化日志
LOG_FILE <- paste0("logs/analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] [%s] %s\n", timestamp, level, msg)
  cat(log_line)
  cat(log_line, file = LOG_FILE, append = TRUE)
}

# 加载或初始化进度
load_progress <- function() {
  if (file.exists(PROGRESS_FILE)) {
    progress <- readRDS(PROGRESS_FILE)
    log_message("加载已有进度记录")
    return(progress)
  } else {
    log_message("初始化新的进度记录")
    return(list(
      step1_download = FALSE,
      step2_prepare = FALSE,
      step3_preprocess = FALSE,
      step4_grouping = FALSE,
      step5_differential = FALSE,
      step6_save_results = FALSE,
      step7_heatmap = FALSE,
      step8_volcano = FALSE,
      step9_pca = FALSE,
      step10_ma = FALSE,
      step11_survival = FALSE,
      step12_clinical = FALSE
    ))
  }
}

# 保存进度
save_progress <- function(progress) {
  saveRDS(progress, PROGRESS_FILE)
}

# 标记步骤完成
mark_step_complete <- function(progress, step_name) {
  progress[[step_name]] <- TRUE
  save_progress(progress)
  log_message(paste("步骤完成:", step_name), "SUCCESS")
  return(progress)
}

# 错误处理包装器
safe_execute <- function(expr, step_name, progress, skip_if_done = TRUE) {
  # 如果已完成且设置跳过，则跳过
  if (skip_if_done && progress[[step_name]]) {
    log_message(paste("跳过已完成的步骤:", step_name), "SKIP")
    return(list(success = TRUE, result = NULL, skipped = TRUE))
  }
  
  log_message(paste("开始执行:", step_name), "START")
  start_time <- Sys.time()
  
  result <- tryCatch({
    res <- eval(expr)
    elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
    log_message(paste("步骤成功完成, 耗时:", elapsed, "秒"), "SUCCESS")
    progress <- mark_step_complete(progress, step_name)
    list(success = TRUE, result = res, elapsed = elapsed, skipped = FALSE)
  }, error = function(e) {
    elapsed <- round(difftime(Sys.time(), start_time, units = "secs"), 2)
    log_message(paste("步骤执行失败:", e$message), "ERROR")
    log_message(paste("错误详情:", paste(capture.output(traceback()), collapse = "\n")), "ERROR")
    list(success = FALSE, error = e$message, elapsed = elapsed, skipped = FALSE)
  })
  
  return(result)
}

# ========================== 包管理 ==========================

log_message("检查和安装必要的R包...")

# 检查并安装CRAN包
install_cran_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    log_message(paste("安装CRAN包:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

# 检查并安装Bioconductor包
install_bioc_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    log_message(paste("安装Bioconductor包:", pkg))
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg, update = FALSE)
  }
}

# 必需的包列表
cran_packages <- c("ggplot2", "dplyr", "pheatmap", "RColorBrewer", 
                   "ggrepel", "gridExtra", "scales", "tidyr")
bioc_packages <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", 
                   "limma", "edgeR", "survival", "survminer")

# 安装包
for (pkg in cran_packages) {
  install_cran_package(pkg)
}

for (pkg in bioc_packages) {
  install_bioc_package(pkg)
}

# 加载包
log_message("加载R包...")
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(survival)
  library(survminer)
  library(gridExtra)
  library(scales)
  library(tidyr)
})

# ========================== 加载进度 ==========================
progress <- load_progress()

# ========================== Step 1: 下载数据 ==========================
step1_result <- safe_execute({
  query_rna <- GDCquery(
    project = PARAMS$project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  log_message(paste("找到", nrow(getResults(query_rna)), "个RNA-Seq文件"))
  
  # 下载数据
  GDCdownload(
    query = query_rna,
    method = "api",
    files.per.chunk = 10
  )
  
  # 保存查询对象
  saveRDS(query_rna, "results/data/query_rna.rds")
  
  query_rna
}, "step1_download", progress)

# 检查是否成功
if (!step1_result$success) {
  log_message("数据下载失败，程序终止", "ERROR")
  quit(status = 1)
}

# 如果是跳过的，加载已有数据
if (step1_result$skipped) {
  query_rna <- readRDS("results/data/query_rna.rds")
} else {
  query_rna <- step1_result$result
}

# ========================== Step 2: 准备数据 ==========================
step2_result <- safe_execute({
  data_rna <- GDCprepare(query_rna)
  
  # 保存原始数据
  saveRDS(data_rna, "results/data/TCGA_RNAseq_raw.rds")
  
  log_message(paste("数据维度:", nrow(data_rna), "基因 x", ncol(data_rna), "样本"))
  
  data_rna
}, "step2_prepare", progress)

if (!step2_result$success) {
  log_message("数据准备失败，程序终止", "ERROR")
  quit(status = 1)
}

if (step2_result$skipped) {
  data_rna <- readRDS("results/data/TCGA_RNAseq_raw.rds")
} else {
  data_rna <- step2_result$result
}

# ========================== Step 3: 数据预处理 ==========================
step3_result <- safe_execute({
  # 提取表达矩阵
  expr_matrix <- assay(data_rna, "unstranded")
  
  # 使用gene_name作为行名（如果重复则使用gene_id）
  gene_names <- rowData(data_rna)$gene_name
  gene_ids <- rowData(data_rna)$gene_id
  
  # 处理重复基因名
  duplicated_genes <- duplicated(gene_names) | is.na(gene_names)
  gene_names[duplicated_genes] <- gene_ids[duplicated_genes]
  rownames(expr_matrix) <- gene_names
  
  # 保存基因信息
  gene_info <- as.data.frame(rowData(data_rna))
  write.csv(gene_info, "results/tables/gene_info.csv", row.names = FALSE)
  
  # 过滤低表达基因
  min_samples <- PARAMS$min_sample_count
  min_count <- PARAMS$min_count
  keep_genes <- rowSums(expr_matrix >= min_count) >= min_samples
  
  expr_filtered <- expr_matrix[keep_genes, ]
  
  log_message(paste("过滤前:", nrow(expr_matrix), "基因"))
  log_message(paste("过滤后:", nrow(expr_filtered), "基因"))
  
  # 保存过滤后的数据
  saveRDS(expr_filtered, "results/data/expr_filtered.rds")
  write.csv(expr_filtered, "results/tables/expression_matrix_filtered.csv")
  
  list(
    expr_raw = expr_matrix,
    expr_filtered = expr_filtered,
    gene_info = gene_info
  )
}, "step3_preprocess", progress)

if (!step3_result$success) {
  log_message("数据预处理失败，程序终止", "ERROR")
  quit(status = 1)
}

if (step3_result$skipped) {
  expr_filtered <- readRDS("results/data/expr_filtered.rds")
  gene_info <- read.csv("results/tables/gene_info.csv")
} else {
  expr_filtered <- step3_result$result$expr_filtered
  gene_info <- step3_result$result$gene_info
}

# ========================== Step 4: 自动识别分组 ==========================
step4_result <- safe_execute({
  # 获取临床数据
  clinical_data <- as.data.frame(colData(data_rna))
  
  # 从样本barcode自动识别分组
  sample_barcodes <- colnames(expr_filtered)
  
  # TCGA barcode格式：TCGA-XX-XXXX-XXA-XXX
  # 第4段的前两位数字表示样本类型：
  # 01 = Primary Solid Tumor
  # 11 = Solid Tissue Normal
  
  sample_types <- sapply(strsplit(sample_barcodes, "-"), function(x) {
    if (length(x) >= 4) {
      substr(x[4], 1, 2)
    } else {
      NA
    }
  })
  
  # 转换为分组标签
  groups <- ifelse(sample_types == "01", "Tumor",
                   ifelse(sample_types == "11", "Normal", NA))
  
  # 统计分组
  group_table <- table(groups, useNA = "always")
  log_message("样本分组统计:")
  log_message(paste(capture.output(print(group_table)), collapse = "\n"))
  
  # 创建样本信息表
  sample_info <- data.frame(
    barcode = sample_barcodes,
    sample_type_code = sample_types,
    group = groups,
    stringsAsFactors = FALSE
  )
  
  # 添加临床信息（如果有）
  if (nrow(clinical_data) > 0) {
    clinical_subset <- clinical_data[match(sample_barcodes, clinical_data$barcode), ]
    sample_info <- cbind(sample_info, clinical_subset)
  }
  
  # 保存样本信息
  write.csv(sample_info, "results/tables/sample_info.csv", row.names = FALSE)
  
  # 过滤掉未识别的样本
  valid_samples <- !is.na(groups)
  expr_valid <- expr_filtered[, valid_samples]
  sample_info_valid <- sample_info[valid_samples, ]
  
  # 分离肿瘤和正常样本
  tumor_idx <- sample_info_valid$group == "Tumor"
  normal_idx <- sample_info_valid$group == "Normal"
  
  data_tumor <- expr_valid[, tumor_idx]
  data_normal <- expr_valid[, normal_idx]
  
  log_message(paste("有效肿瘤样本:", sum(tumor_idx), "个"))
  log_message(paste("有效正常样本:", sum(normal_idx), "个"))
  
  # 保存分组数据
  saveRDS(list(
    expr_valid = expr_valid,
    sample_info = sample_info_valid,
    data_tumor = data_tumor,
    data_normal = data_normal
  ), "results/data/grouped_data.rds")
  
  list(
    expr_valid = expr_valid,
    sample_info = sample_info_valid,
    data_tumor = data_tumor,
    data_normal = data_normal
  )
}, "step4_grouping", progress)

if (!step4_result$success) {
  log_message("分组识别失败，程序终止", "ERROR")
  quit(status = 1)
}

if (step4_result$skipped) {
  grouped_data <- readRDS("results/data/grouped_data.rds")
  expr_valid <- grouped_data$expr_valid
  sample_info_valid <- grouped_data$sample_info
  data_tumor <- grouped_data$data_tumor
  data_normal <- grouped_data$data_normal
} else {
  expr_valid <- step4_result$result$expr_valid
  sample_info_valid <- step4_result$result$sample_info
  data_tumor <- step4_result$result$data_tumor
  data_normal <- step4_result$result$data_normal
}

# ========================== Step 5: 差异表达分析（DESeq2） ==========================
step5_result <- safe_execute({
  # 创建DESeq2数据集
  dds <- DESeqDataSetFromMatrix(
    countData = expr_valid,
    colData = sample_info_valid,
    design = ~ group
  )
  
  # 运行DESeq2
  dds <- DESeq(dds)
  
  # 获取结果
  res <- results(dds, contrast = c("group", "Tumor", "Normal"))
  
  # 转换为数据框
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # 添加调控方向
  res_df$regulation <- ifelse(res_df$log2FoldChange > 0, "Up", 
                               ifelse(res_df$log2FoldChange < 0, "Down", "No change"))
  
  # 添加显著性标签
  res_df$significant <- ifelse(
    abs(res_df$log2FoldChange) > PARAMS$threshold_logFC & 
      res_df$padj < PARAMS$threshold_adjP,
    res_df$regulation,
    "Not Significant"
  )
  
  # 筛选显著差异基因
  sig_genes <- res_df[res_df$significant != "Not Significant", ]
  sig_genes <- sig_genes[order(sig_genes$padj), ]
  
  log_message(paste("差异表达基因总数:", nrow(sig_genes)))
  log_message(paste("上调基因:", sum(sig_genes$regulation == "Up")))
  log_message(paste("下调基因:", sum(sig_genes$regulation == "Down")))
  
  # 保存结果
  saveRDS(list(dds = dds, results = res_df, sig_genes = sig_genes), 
          "results/data/deseq2_results.rds")
  
  list(dds = dds, results = res_df, sig_genes = sig_genes)
}, "step5_differential", progress)

if (!step5_result$success) {
  log_message("差异表达分析失败", "WARNING")
  # 尝试使用limma替代
  log_message("尝试使用limma进行差异分析...", "INFO")
  
  step5_result <- safe_execute({
    # log2转换（添加伪计数）
    expr_log <- log2(expr_valid + 1)
    
    # 归一化
    expr_norm <- normalizeBetweenArrays(expr_log)
    
    # 设计矩阵
    design <- model.matrix(~ 0 + sample_info_valid$group)
    colnames(design) <- c("Normal", "Tumor")
    
    # 拟合模型
    fit <- lmFit(expr_norm, design)
    
    # 对比矩阵
    contrast_matrix <- makeContrasts(Tumor - Normal, levels = design)
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # 获取结果
    res_df <- topTable(fit2, number = Inf, adjust.method = "fdr")
    res_df$gene <- rownames(res_df)
    
    # 添加调控方向和显著性
    res_df$regulation <- ifelse(res_df$logFC > 0, "Up", 
                                 ifelse(res_df$logFC < 0, "Down", "No change"))
    res_df$significant <- ifelse(
      abs(res_df$logFC) > PARAMS$threshold_logFC & 
        res_df$adj.P.Val < PARAMS$threshold_adjP,
      res_df$regulation,
      "Not Significant"
    )
    
    # 重命名列以匹配DESeq2格式
    colnames(res_df)[colnames(res_df) == "logFC"] <- "log2FoldChange"
    colnames(res_df)[colnames(res_df) == "adj.P.Val"] <- "padj"
    colnames(res_df)[colnames(res_df) == "P.Value"] <- "pvalue"
    
    sig_genes <- res_df[res_df$significant != "Not Significant", ]
    sig_genes <- sig_genes[order(sig_genes$padj), ]
    
    log_message(paste("差异表达基因总数:", nrow(sig_genes)))
    log_message(paste("上调基因:", sum(sig_genes$regulation == "Up")))
    log_message(paste("下调基因:", sum(sig_genes$regulation == "Down")))
    
    saveRDS(list(results = res_df, sig_genes = sig_genes), 
            "results/data/limma_results.rds")
    
    list(results = res_df, sig_genes = sig_genes, method = "limma")
  }, "step5_differential", progress, skip_if_done = FALSE)
}

if (step5_result$skipped) {
  if (file.exists("results/data/deseq2_results.rds")) {
    deseq_results <- readRDS("results/data/deseq2_results.rds")
    res_df <- deseq_results$results
    sig_genes <- deseq_results$sig_genes
  } else {
    limma_results <- readRDS("results/data/limma_results.rds")
    res_df <- limma_results$results
    sig_genes <- limma_results$sig_genes
  }
} else {
  res_df <- step5_result$result$results
  sig_genes <- step5_result$result$sig_genes
}

# ========================== Step 6: 保存差异分析结果表格 ==========================
step6_result <- safe_execute({
  write.csv(res_df, "results/tables/differential_expression_all.csv", row.names = FALSE)
  write.csv(sig_genes, "results/tables/differential_expression_significant.csv", row.names = FALSE)
  
  # 保存差异基因列表
  write.table(sig_genes$gene, "results/tables/DEG_gene_list.txt", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  TRUE
}, "step6_save_results", progress)

# ========================== Step 7: 热图（优化版） ==========================
step7_result <- safe_execute({
  # 获取归一化后的表达数据用于可视化
  if (exists("dds")) {
    # 使用DESeq2的rlog转换
    rld <- rlog(dds, blind = FALSE)
    expr_for_heatmap <- assay(rld)
  } else {
    # 使用log2转换
    expr_for_heatmap <- log2(expr_valid + 1)
  }
  
  # 选择展示的基因
  n_genes <- min(PARAMS$max_display_genes * 2, nrow(sig_genes))
  
  # 按logFC排序，选择top上调和下调基因
  sig_genes_ordered <- sig_genes[order(sig_genes$log2FoldChange, decreasing = TRUE), ]
  
  top_up <- head(sig_genes_ordered, PARAMS$max_display_genes)
  top_down <- tail(sig_genes_ordered, PARAMS$max_display_genes)
  
  selected_genes <- c(top_up$gene, top_down$gene)
  selected_genes <- selected_genes[selected_genes %in% rownames(expr_for_heatmap)]
  
  heatmap_data <- expr_for_heatmap[selected_genes, ]
  
  # 准备注释
  annotation_col <- data.frame(
    Group = sample_info_valid$group,
    row.names = colnames(heatmap_data)
  )
  
  # 颜色设置
  ann_colors <- list(
    Group = c(Tumor = "#E74C3C", Normal = "#3498DB")
  )
  
  color_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  # 绘制热图
  pdf("results/figures/heatmap_DEGs.pdf", width = 12, height = 10)
  pheatmap(
    heatmap_data,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = color_palette,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    scale = "row",
    fontsize = 10,
    fontsize_row = 8,
    border_color = NA,
    main = paste("Top", length(selected_genes), "Differentially Expressed Genes")
  )
  dev.off()
  
  log_message("热图已保存")
  
  TRUE
}, "step7_heatmap", progress)

# ========================== Step 8: 火山图（SCI风格） ==========================
step8_result <- safe_execute({
  # 准备数据
  volcano_data <- res_df %>%
    filter(!is.na(padj) & !is.na(log2FoldChange))
  
  # 选择top基因用于标注
  top_n <- 20
  top_genes <- volcano_data %>%
    filter(significant != "Not Significant") %>%
    arrange(padj) %>%
    head(top_n)
  
  # 设置颜色
  colors <- c(
    "Up" = "#E74C3C",
    "Down" = "#3498DB",
    "Not Significant" = "#95A5A6"
  )
  
  # 绘制火山图
  p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), size = 2, alpha = 0.6) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-PARAMS$threshold_logFC, PARAMS$threshold_logFC), 
               linetype = "dashed", color = "gray40", size = 0.7) +
    geom_hline(yintercept = -log10(PARAMS$threshold_adjP), 
               linetype = "dashed", color = "gray40", size = 0.7) +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3.5,
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "gray50"
    ) +
    labs(
      title = "Volcano Plot of Differential Gene Expression",
      subtitle = paste0("Upregulated: ", sum(sig_genes$regulation == "Up"),
                       " | Downregulated: ", sum(sig_genes$regulation == "Down")),
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~Adjusted~P~value),
      color = "Regulation"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "gray90", linetype = "dotted")
    )
  
  ggsave("results/figures/volcano_plot.pdf", p, width = 10, height = 8)
  
  log_message("火山图已保存")
  
  TRUE
}, "step8_volcano", progress)

# ========================== Step 9: PCA分析（优化版） ==========================
step9_result <- safe_execute({
  # 准备PCA输入数据
  if (exists("dds")) {
    # 使用variance stabilizing transformation
    vsd <- vst(dds, blind = FALSE)
    pca_data <- assay(vsd)
  } else {
    # 使用log2转换
    pca_data <- log2(expr_valid + 1)
  }
  
  # 转置（样本为行，基因为列）
  pca_input <- t(pca_data)
  
  # 运行PCA
  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
  
  # 提取前两个主成分
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Group = sample_info_valid$group,
    Sample = rownames(pca_result$x)
  )
  
  # 计算方差解释度
  var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)
  
  # 绘制PCA图
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Group)) +
    stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, 
                 color = NA, show.legend = FALSE) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c(Tumor = "#E74C3C", Normal = "#3498DB")) +
    scale_fill_manual(values = c(Tumor = "#E74C3C", Normal = "#3498DB")) +
    scale_shape_manual(values = c(Tumor = 16, Normal = 17)) +
    labs(
      title = "Principal Component Analysis",
      subtitle = paste(PARAMS$project, "RNA-Seq Data"),
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)")
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
      axis.title = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "gray90", linetype = "dotted")
    )
  
  ggsave("results/figures/PCA_plot.pdf", p, width = 10, height = 8)
  
  # 保存PCA结果
  saveRDS(list(pca_result = pca_result, pca_df = pca_df), 
          "results/data/pca_results.rds")
  
  log_message("PCA分析已完成")
  
  TRUE
}, "step9_pca", progress)

# ========================== Step 10: MA图 ==========================
step10_result <- safe_execute({
  ma_data <- res_df %>%
    filter(!is.na(baseMean) & !is.na(log2FoldChange))
  
  p <- ggplot(ma_data, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
    geom_point(aes(color = significant), size = 1.5, alpha = 0.5) +
    scale_color_manual(values = c(Up = "#E74C3C", Down = "#3498DB", 
                                   "Not Significant" = "#95A5A6")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    geom_hline(yintercept = c(-PARAMS$threshold_logFC, PARAMS$threshold_logFC), 
               linetype = "dashed", color = "gray40") +
    labs(
      title = "MA Plot",
      subtitle = "Mean Expression vs Log Fold Change",
      x = expression(Log[10]~Mean~Expression),
      y = expression(Log[2]~Fold~Change),
      color = "Regulation"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
      axis.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  ggsave("results/figures/MA_plot.pdf", p, width = 10, height = 7)
  
  log_message("MA图已保存")
  
  TRUE
}, "step10_ma", progress)

# ========================== Step 11: 生存分析 ==========================
step11_result <- safe_execute({
  # 准备生存数据
  survival_data <- sample_info_valid %>%
    filter(group == "Tumor") %>%
    mutate(
      OS_time = as.numeric(days_to_death),
      OS_time = ifelse(is.na(OS_time), as.numeric(days_to_last_follow_up), OS_time),
      OS_status = ifelse(vital_status == "Dead", 1, 0)
    ) %>%
    filter(!is.na(OS_time) & OS_time > 0)
  
  if (nrow(survival_data) < 10) {
    log_message("生存数据样本数不足，跳过生存分析", "WARNING")
    return(FALSE)
  }
  
  # 选择一个感兴趣的基因（例如top差异基因）
  gene_of_interest <- sig_genes$gene[1]
  
  # 获取该基因的表达值
  if (exists("dds")) {
    expr_values <- counts(dds, normalized = TRUE)[gene_of_interest, survival_data$barcode]
  } else {
    expr_values <- expr_valid[gene_of_interest, survival_data$barcode]
  }
  
  # 按中位数分组
  survival_data$expression_group <- ifelse(
    expr_values > median(expr_values, na.rm = TRUE),
    "High", "Low"
  )
  
  # 创建生存对象
  surv_obj <- Surv(time = survival_data$OS_time / 365, 
                   event = survival_data$OS_status)
  
  # 拟合生存曲线
  fit <- survfit(surv_obj ~ expression_group, data = survival_data)
  
  # 绘制生存曲线
  p <- ggsurvplot(
    fit,
    data = survival_data,
    pval = TRUE,
    risk.table = TRUE,
    conf.int = TRUE,
    palette = c("#3498DB", "#E74C3C"),
    title = paste("Survival Analysis:", gene_of_interest),
    xlab = "Time (Years)",
    ylab = "Survival Probability",
    legend.title = "Expression",
    legend.labs = c("High", "Low"),
    risk.table.height = 0.25
  )
  
  pdf("results/figures/survival_plot.pdf", width = 10, height = 10)
  print(p)
  dev.off()
  
  log_message("生存分析已完成")
  
  TRUE
}, "step11_survival", progress)

# ========================== Step 12: 下载和处理临床数据 ==========================
step12_result <- safe_execute({
  # 查询临床数据
  clinical <- GDCquery_clinic(project = PARAMS$project, type = "clinical")
  
  # 保存临床数据
  write.csv(clinical, "results/tables/clinical_data_complete.csv", row.names = FALSE)
  
  log_message(paste("临床数据已下载:", nrow(clinical), "个病例"))
  
  # 生成临床特征统计
  if (nrow(clinical) > 0) {
    # 创建统计报告
    clinical_summary <- list(
      total_patients = nrow(clinical),
      gender = table(clinical$gender, useNA = "always"),
      vital_status = table(clinical$vital_status, useNA = "always"),
      race = table(clinical$race, useNA = "always")
    )
    
    # 如果有年龄数据
    if ("age_at_diagnosis" %in% colnames(clinical)) {
      clinical_summary$age_summary <- summary(clinical$age_at_diagnosis)
    }
    
    # 保存统计摘要
    capture.output(clinical_summary, file = "results/tables/clinical_summary.txt")
    
    log_message("临床数据统计已生成")
  }
  
  TRUE
}, "step12_clinical", progress)

# ========================== 生成分析报告 ==========================
cat("\n")
cat("="*80, "\n")
cat("分析完成！\n")
cat("="*80, "\n\n")

log_message("生成分析报告...")

# 创建报告
report <- paste0(
  "TCGA分析完成报告\n",
  "="*80, "\n",
  "项目: ", PARAMS$project, "\n",
  "分析时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "工作目录: ", WORK_DIR, "\n\n",
  "数据统计:\n",
  "- 总基因数: ", nrow(expr_filtered), "\n",
  "- 总样本数: ", ncol(expr_valid), "\n",
  "- 肿瘤样本: ", sum(sample_info_valid$group == "Tumor"), "\n",
  "- 正常样本: ", sum(sample_info_valid$group == "Normal"), "\n\n",
  "差异分析结果:\n",
  "- 显著差异基因: ", nrow(sig_genes), "\n",
  "- 上调基因: ", sum(sig_genes$regulation == "Up"), "\n",
  "- 下调基因: ", sum(sig_genes$regulation == "Down"), "\n\n",
  "输出文件:\n",
  "results/tables/\n",
  "  - differential_expression_all.csv\n",
  "  - differential_expression_significant.csv\n",
  "  - DEG_gene_list.txt\n",
  "  - clinical_data_complete.csv\n",
  "  - sample_info.csv\n\n",
  "results/figures/\n",
  "  - heatmap_DEGs.pdf\n",
  "  - volcano_plot.pdf\n",
  "  - PCA_plot.pdf\n",
  "  - MA_plot.pdf\n",
  "  - survival_plot.pdf\n\n",
  "="*80, "\n"
)

cat(report)
writeLines(report, "results/analysis_report.txt")

log_message("分析完成！所有结果已保存到 results/ 目录")
log_message(paste("详细日志:", LOG_FILE))
