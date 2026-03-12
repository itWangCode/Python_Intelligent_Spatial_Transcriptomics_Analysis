#!/usr/bin/env Rscript
# ==============================================================================
# 模块1：数据预处理和基础分析
# 功能：数据读取、质控、聚类、Marker分析、细胞注释
# ==============================================================================

# ==============================================================================
# Module 1.1: 数据读取和预处理
# ==============================================================================

read_and_process_10x_data <- function(workDir, output_dirs) {
  message("=== 步骤1.1: 读取10X数据 ===")
  
  # 获取所有样本目录
  dirs <- list.dirs(workDir)
  dirs_sample <- dirs[grepl("^[^.]", basename(dirs)) & dirs != workDir]
  
  if (length(dirs_sample) == 0) {
    # 如果没有子目录，尝试直接读取当前目录
    message("未找到样本子目录，尝试读取当前目录...")
    dirs_sample <- workDir
  }
  
  names(dirs_sample) <- gsub(".+\\/(.+)", "\\1", dirs_sample)
  
  # 创建进度条
  if (requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      format = " 读取进度 [:bar] :percent | 时间: :elapsed | 剩余: :eta",
      total = length(dirs_sample),
      clear = FALSE,
      width = 60
    )
    use_pb <- TRUE
  } else {
    use_pb <- FALSE
  }
  
  # 逐一读取数据
  counts <- list()
  for (i in seq_along(dirs_sample)) {
    message(sprintf("正在读取: %s", names(dirs_sample)[i]))
    tryCatch({
      counts[[i]] <- Read10X(data.dir = dirs_sample[i])
      if (use_pb) pb$tick()
    }, error = function(e) {
      message(sprintf("读取失败: %s - %s", names(dirs_sample)[i], e$message))
      if (use_pb) pb$tick()
    })
  }
  
  # 合并数据
  if (length(counts) == 0) {
    stop("没有成功读取任何数据")
  }
  
  counts_combined <- do.call(cbind, counts)
  
  message(sprintf("✓ 成功读取 %d 个样本，共 %d 个基因 × %d 个细胞",
                  length(counts), nrow(counts_combined), ncol(counts_combined)))
  
  return(list(
    counts = counts_combined,
    sample_names = names(dirs_sample)
  ))
}

# ==============================================================================
# Module 1.2: 创建Seurat对象
# ==============================================================================

create_seurat_with_qc <- function(counts_data, params, output_dirs) {
  message("\n=== 步骤1.2: 创建Seurat对象并质控 ===")
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = counts_data$counts,
    min.cells = params$min_cells,
    min.features = params$min_features,
    project = "SingleCell_Complete_Analysis"
  )
  
  message(sprintf("初始细胞数: %d", ncol(seurat_obj)))
  
  # 添加分组信息
  Type <- gsub("(.*?)\\..*", "\\1", colnames(seurat_obj))
  names(Type) <- colnames(seurat_obj)
  seurat_obj <- AddMetaData(
    object = seurat_obj,
    metadata = Type,
    col.name = "Type"
  )
  
  # 计算线粒体基因比例
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # 计算核糖体基因比例
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  # 生成QC报告
  qc_stats <- data.frame(
    CellID = colnames(seurat_obj),
    nCount_RNA = seurat_obj$nCount_RNA,
    nFeature_RNA = seurat_obj$nFeature_RNA,
    percent.mt = seurat_obj$percent.mt,
    percent.rb = seurat_obj$percent.rb,
    Type = seurat_obj$Type
  )
  
  write.csv(qc_stats,
           file = file.path(output_dirs$module1_qc, "QC_Statistics.csv"),
           row.names = FALSE)
  
  # QC可视化
  pdf(file = file.path(output_dirs$module1_qc, "QC_Violin_Plots.pdf"),
      width = 12, height = 8)
  print(
    VlnPlot(seurat_obj, 
           features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
           ncol = 2, pt.size = 0.1)
  )
  dev.off()
  
  # 散点图QC
  pdf(file = file.path(output_dirs$module1_qc, "QC_Scatter_Plots.pdf"),
      width = 12, height = 5)
  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
  
  message("✓ Seurat对象创建完成，质控图表已保存")
  
  return(seurat_obj)
}

# ==============================================================================
# Module 1.3: 数据归一化和降维
# ==============================================================================

normalize_and_reduce <- function(seurat_obj, params, output_dirs) {
  message("\n=== 步骤1.3: 数据归一化和降维 ===")
  
  # 创建进度指示
  steps <- c("归一化", "高变基因", "标准化", "PCA", "邻居", "聚类", "UMAP")
  
  # 1. 归一化
  message("1/7: 数据归一化...")
  seurat_obj <- NormalizeData(
    object = seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  # 2. 识别高变基因
  message("2/7: 识别高变基因...")
  seurat_obj <- FindVariableFeatures(
    object = seurat_obj,
    selection.method = "vst",
    nfeatures = params$n_variable_features
  )
  
  # 可视化高变基因
  top10 <- head(VariableFeatures(seurat_obj), 10)
  
  pdf(file = file.path(output_dirs$module1_qc, "Variable_Features.pdf"),
      width = 10, height = 6)
  plot1 <- VariableFeaturePlot(seurat_obj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot2)
  dev.off()
  
  # 3. 标准化
  message("3/7: 数据标准化...")
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  
  # 4. PCA分析
  message("4/7: PCA降维...")
  seurat_obj <- RunPCA(
    object = seurat_obj,
    features = VariableFeatures(object = seurat_obj),
    npcs = 50,
    verbose = FALSE
  )
  
  # PCA可视化
  pdf(file = file.path(output_dirs$module1_clustering, "PCA_Analysis.pdf"),
      width = 12, height = 10)
  print(DimPlot(seurat_obj, reduction = "pca"))
  print(ElbowPlot(seurat_obj, ndims = 50))
  print(DimHeatmap(seurat_obj, dims = 1:12, cells = 500, balanced = TRUE))
  dev.off()
  
  # 5. 构建KNN图
  message("5/7: 构建邻居图...")
  seurat_obj <- FindNeighbors(
    object = seurat_obj,
    dims = 1:params$pcSelect
  )
  
  # 6. 聚类
  message("6/7: 聚类分析...")
  seurat_obj <- FindClusters(
    object = seurat_obj,
    resolution = params$cluster_resolution,
    verbose = FALSE
  )
  
  # 7. UMAP降维
  message("7/7: UMAP降维...")
  seurat_obj <- RunUMAP(
    object = seurat_obj,
    dims = 1:params$pcSelect,
    verbose = FALSE
  )
  
  # UMAP可视化
  pdf(file = file.path(output_dirs$module1_clustering, "UMAP_Basic.pdf"),
      width = 10, height = 8)
  p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1.5)
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "Type", pt.size = 1.5)
  print(p1)
  print(p2)
  dev.off()
  
  # 分组UMAP
  pdf(file = file.path(output_dirs$module1_clustering, "UMAP_Split_by_Group.pdf"),
      width = 14, height = 6)
  print(
    DimPlot(seurat_obj, reduction = "umap", split.by = "Type",
           label = TRUE, pt.size = 1.2, ncol = 2)
  )
  dev.off()
  
  # 保存聚类信息
  write.table(
    seurat_obj$seurat_clusters,
    file = file.path(output_dirs$module1_clustering, "Cluster_Assignment.txt"),
    quote = FALSE,
    sep = "\t",
    col.names = FALSE
  )
  
  message("✓ 归一化和降维完成")
  
  return(seurat_obj)
}

# ==============================================================================
# Module 1.4: Marker基因分析
# ==============================================================================

find_all_markers <- function(seurat_obj, params, output_dirs) {
  message("\n=== 步骤1.4: Marker基因分析 ===")
  
  # 寻找所有marker基因
  message("识别marker基因（这可能需要一些时间）...")
  all_markers <- FindAllMarkers(
    object = seurat_obj,
    only.pos = FALSE,
    min.pct = 0.25,
    logfc.threshold = params$logFCfilter,
    verbose = FALSE
  )
  
  # 筛选显著marker
  sig_markers <- all_markers %>%
    filter(abs(avg_log2FC) > params$logFCfilter & p_val_adj < params$adjPvalFilter)
  
  message(sprintf("找到 %d 个显著marker基因", nrow(sig_markers)))
  
  # 保存结果
  write.csv(all_markers,
           file = file.path(output_dirs$module1_markers, "All_Markers.csv"),
           row.names = FALSE)
  
  write.csv(sig_markers,
           file = file.path(output_dirs$module1_markers, "Significant_Markers.csv"),
           row.names = FALSE)
  
  # Top marker基因
  top_markers <- sig_markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
  
  write.csv(top_markers,
           file = file.path(output_dirs$module1_markers, "Top10_Markers_Per_Cluster.csv"),
           row.names = FALSE)
  
  # 热图可视化
  if (nrow(top_markers) > 0) {
    pdf(file = file.path(output_dirs$module1_heatmaps, "Marker_Heatmap.pdf"),
        width = 14, height = 12)
    print(
      DoHeatmap(seurat_obj, features = unique(top_markers$gene),
               size = 3, angle = 45) +
        scale_fill_gradientn(colors = c("blue", "white", "red"))
    )
    dev.off()
  }
  
  # 小提琴图（Top基因）
  top_genes <- head(top_markers$gene, 12)
  if (length(top_genes) > 0) {
    pdf(file = file.path(output_dirs$module1_markers, "Top_Markers_Violin.pdf"),
        width = 15, height = 12)
    print(
      VlnPlot(seurat_obj, features = top_genes, ncol = 3, pt.size = 0.1)
    )
    dev.off()
  }
  
  # 特征图
  if (length(top_genes) > 0) {
    pdf(file = file.path(output_dirs$module1_markers, "Top_Markers_Feature.pdf"),
        width = 15, height = 12)
    print(
      FeaturePlot(seurat_obj, features = top_genes[1:min(9, length(top_genes))],
                 ncol = 3, pt.size = 0.5)
    )
    dev.off()
  }
  
  message("✓ Marker基因分析完成")
  
  return(list(
    all_markers = all_markers,
    sig_markers = sig_markers,
    top_markers = top_markers
  ))
}

# ==============================================================================
# Module 1.5: 自动细胞类型注释（SingleR）
# ==============================================================================

automatic_cell_annotation <- function(seurat_obj, output_dirs) {
  message("\n=== 步骤1.5: 自动细胞类型注释 ===")
  
  tryCatch({
    # 加载参考数据库
    message("加载参考数据库...")
    ref <- celldex::HumanPrimaryCellAtlasData()
    
    # 提取表达矩阵
    test_data <- GetAssayData(seurat_obj, slot = "data")
    
    # 运行SingleR
    message("运行SingleR注释（这可能需要几分钟）...")
    singler_results <- SingleR(
      test = test_data,
      ref = ref,
      labels = ref$label.main
    )
    
    # 添加到Seurat对象
    seurat_obj$singler_labels <- singler_results$labels
    
    # 保存结果
    annotation_df <- data.frame(
      CellID = colnames(seurat_obj),
      Cluster = seurat_obj$seurat_clusters,
      SingleR_Annotation = singler_results$labels,
      SingleR_Score = singler_results$scores
    )
    
    write.csv(annotation_df,
             file = file.path(output_dirs$module1_annotation, "SingleR_Annotation.csv"),
             row.names = FALSE)
    
    # 聚类与注释的对应关系
    cluster_celltype <- table(seurat_obj$seurat_clusters, singler_results$labels)
    write.csv(cluster_celltype,
             file = file.path(output_dirs$module1_annotation, "Cluster_CellType_Table.csv"))
    
    # 可视化
    pdf(file = file.path(output_dirs$module1_annotation, "SingleR_UMAP.pdf"),
        width = 12, height = 10)
    p1 <- DimPlot(seurat_obj, group.by = "singler_labels", label = TRUE,
                 repel = TRUE, pt.size = 1) +
      ggtitle("SingleR Automatic Annotation")
    p2 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE,
                 pt.size = 1) +
      ggtitle("Seurat Clusters")
    print(p1)
    print(p2)
    dev.off()
    
    # 热图：聚类 vs 细胞类型
    pdf(file = file.path(output_dirs$module1_annotation, "Cluster_CellType_Heatmap.pdf"),
        width = 10, height = 8)
    pheatmap::pheatmap(
      as.matrix(cluster_celltype),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      color = colorRampPalette(c("white", "blue"))(50),
      main = "Cluster vs Cell Type"
    )
    dev.off()
    
    message("✓ 自动注释完成")
    
    return(seurat_obj)
    
  }, error = function(e) {
    warning(sprintf("自动注释失败: %s", e$message))
    return(seurat_obj)
  })
}

# ==============================================================================
# Module 1.6: 手动细胞类型注释
# ==============================================================================

manual_cell_annotation <- function(seurat_obj, cluster_annotation, output_dirs) {
  message("\n=== 步骤1.6: 应用手动细胞类型注释 ===")
  
  if (is.null(cluster_annotation)) {
    message("未提供手动注释，使用SingleR注释或保持原聚类标签")
    
    if ("singler_labels" %in% colnames(seurat_obj@meta.data)) {
      Idents(seurat_obj) <- seurat_obj$singler_labels
      seurat_obj$cell_type <- Idents(seurat_obj)
    } else {
      Idents(seurat_obj) <- seurat_obj$seurat_clusters
      seurat_obj$cell_type <- paste0("Cluster_", Idents(seurat_obj))
    }
  } else {
    # 应用手动注释
    current_clusters <- levels(Idents(seurat_obj))
    new_ids <- sapply(current_clusters, function(x) {
      if (x %in% names(cluster_annotation)) {
        return(cluster_annotation[x])
      } else {
        return(paste0("Unknown_", x))
      }
    })
    
    names(new_ids) <- current_clusters
    seurat_obj <- RenameIdents(seurat_obj, new_ids)
    seurat_obj$cell_type <- Idents(seurat_obj)
    
    message(sprintf("✓ 已应用 %d 个聚类的手动注释", length(cluster_annotation)))
  }
  
  # 保存注释映射
  annotation_map <- data.frame(
    Cluster = seurat_obj$seurat_clusters,
    CellType = seurat_obj$cell_type
  )
  
  write.csv(unique(annotation_map),
           file = file.path(output_dirs$module1_annotation, "Final_Cell_Type_Mapping.csv"),
           row.names = FALSE)
  
  # 注释后的UMAP
  pdf(file = file.path(output_dirs$module1_annotation, "Annotated_UMAP.pdf"),
      width = 12, height = 10)
  p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1.5,
               label.size = 5) +
    ggtitle("Cell Type Annotation") +
    theme_minimal()
  print(p1)
  dev.off()
  
  # 分组注释UMAP
  pdf(file = file.path(output_dirs$module1_annotation, "Annotated_UMAP_Split.pdf"),
      width = 16, height = 7)
  p2 <- DimPlot(seurat_obj, reduction = "umap", split.by = "Type",
               label = TRUE, pt.size = 1.2, ncol = 2) +
    ggtitle("Cell Type Distribution by Group")
  print(p2)
  dev.off()
  
  message("✓ 细胞类型注释完成")
  
  return(seurat_obj)
}

# ==============================================================================
# Module 1.7: 差异表达分析
# ==============================================================================

perform_differential_expression <- function(seurat_obj, params, output_dirs) {
  message("\n=== 步骤1.7: 差异表达分析 ===")
  
  tryCatch({
    # 检查Type分组
    if (!("Type" %in% colnames(seurat_obj@meta.data))) {
      warning("没有Type分组，跳过差异分析")
      return(seurat_obj)
    }
    
    unique_types <- unique(seurat_obj$Type)
    if (length(unique_types) < 2) {
      warning("Type分组少于2组，跳过差异分析")
      return(seurat_obj)
    }
    
    message(sprintf("比较组别: %s vs %s", unique_types[1], unique_types[2]))
    
    # 设置分组
    Idents(seurat_obj) <- seurat_obj$Type
    
    # 差异分析
    diff_genes <- FindMarkers(
      seurat_obj,
      ident.1 = unique_types[1],
      ident.2 = unique_types[2],
      logfc.threshold = params$logFCfilter,
      min.pct = 0.25,
      verbose = FALSE
    )
    
    diff_genes$gene <- rownames(diff_genes)
    
    # 筛选显著差异基因
    sig_diff <- diff_genes %>%
      filter(abs(avg_log2FC) > params$logFCfilter & p_val_adj < params$adjPvalFilter)
    
    message(sprintf("找到 %d 个显著差异基因", nrow(sig_diff)))
    
    # 保存结果
    write.csv(diff_genes,
             file = file.path(output_dirs$differential, "All_Differential_Genes.csv"),
             row.names = FALSE)
    
    write.csv(sig_diff,
             file = file.path(output_dirs$differential, "Significant_Differential_Genes.csv"),
             row.names = FALSE)
    
    # 火山图
    volcano_data <- diff_genes
    volcano_data$significance <- "NS"
    volcano_data$significance[volcano_data$avg_log2FC > 1 & volcano_data$p_val_adj < 0.05] <- "Up"
    volcano_data$significance[volcano_data$avg_log2FC < -1 & volcano_data$p_val_adj < 0.05] <- "Down"
    
    pdf(file = file.path(output_dirs$differential, "Volcano_Plot.pdf"),
        width = 10, height = 8)
    
    p <- ggplot(volcano_data, aes(x = avg_log2FC, y = -log10(p_val_adj), 
                                  color = significance)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
      theme_minimal() +
      labs(title = "Volcano Plot", x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed")
    
    # 标注top基因
    top_genes <- volcano_data %>%
      filter(significance != "NS") %>%
      arrange(p_val_adj) %>%
      head(20)
    
    if (nrow(top_genes) > 0 & requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = 3,
        max.overlaps = 20
      )
    }
    
    print(p)
    dev.off()
    
    # 恢复原始ident
    Idents(seurat_obj) <- seurat_obj$cell_type
    
    # 保存到对象
    seurat_obj@misc$diff_genes <- diff_genes
    seurat_obj@misc$sig_diff_genes <- sig_diff
    
    message("✓ 差异表达分析完成")
    
    return(seurat_obj)
    
  }, error = function(e) {
    warning(sprintf("差异表达分析失败: %s", e$message))
    Idents(seurat_obj) <- seurat_obj$cell_type
    return(seurat_obj)
  })
}

# ==============================================================================
# Module 1.8: 细胞比例分析
# ==============================================================================

analyze_proportions <- function(seurat_obj, output_dirs) {
  message("\n=== 步骤1.8: 细胞比例分析 ===")
  
  tryCatch({
    # 计算比例
    cell_counts <- table(Idents(seurat_obj), seurat_obj$Type)
    proportions <- prop.table(cell_counts, margin = 2) * 100
    
    # 数据框
    prop_df <- as.data.frame(proportions)
    colnames(prop_df) <- c("CellType", "Group", "Percentage")
    
    # 柱状图
    p <- ggplot(prop_df, aes(x = CellType, y = Percentage, fill = Group)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Cell Type Proportions", x = "Cell Type",
           y = "Percentage (%)", fill = "Group") +
      scale_fill_brewer(palette = "Set2")
    
    ggsave(p,
          filename = file.path(output_dirs$proportions, "Cell_Proportions.pdf"),
          width = 10, height = 6)
    
    # 保存数据
    write.csv(prop_df,
             file = file.path(output_dirs$proportions, "Cell_Proportions.csv"),
             row.names = FALSE)
    
    # 统计检验
    if (requireNamespace("rstatix", quietly = TRUE)) {
      chi_test <- chisq.test(cell_counts)
      
      test_result <- data.frame(
        Test = "Chi-square test",
        Statistic = chi_test$statistic,
        P_value = chi_test$p.value
      )
      
      write.csv(test_result,
               file = file.path(output_dirs$proportions, "Chi_Square_Test.csv"),
               row.names = FALSE)
    }
    
    message("✓ 细胞比例分析完成")
    
    return(prop_df)
    
  }, error = function(e) {
    warning(sprintf("细胞比例分析失败: %s", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 1.9: 保存结果
# ==============================================================================

save_module1_results <- function(seurat_obj, output_dirs, workDir) {
  message("\n=== 步骤1.9: 保存模块1结果 ===")
  
  # 保存Seurat对象和注释
  cellAnn <- as.character(Idents(seurat_obj))
  pbmc <- seurat_obj  # 保持命名一致性
  
  save(pbmc, cellAnn,
       file = file.path(output_dirs$rdata_storage, "Seurat.RData"))
  
  message("✓ Seurat.RData已保存")
  
  # 保存RDS文件
  saveRDS(seurat_obj,
         file = file.path(output_dirs$rdata_storage, "Seurat_Object.rds"))
  
  message("✓ Seurat_Object.rds已保存")
  
  # 保存表达矩阵
  expr_matrix <- GetAssayData(seurat_obj, slot = "counts")
  saveRDS(expr_matrix,
         file = file.path(output_dirs$rdata_storage, "Expression_Matrix.rds"))
  
  message("✓ Expression_Matrix.rds已保存")
  
  # 保存元数据
  metadata <- seurat_obj@meta.data
  write.csv(metadata,
           file = file.path(output_dirs$rdata_storage, "Metadata.csv"),
           row.names = TRUE)
  
  message("✓ Metadata.csv已保存")
  
  # 生成总结报告
  summary_stats <- data.frame(
    Metric = c("Total Cells", "Total Genes", "Number of Clusters",
               "Number of Cell Types", "Number of Groups"),
    Value = c(ncol(seurat_obj), nrow(seurat_obj),
             length(unique(seurat_obj$seurat_clusters)),
             length(unique(seurat_obj$cell_type)),
             length(unique(seurat_obj$Type)))
  )
  
  write.csv(summary_stats,
           file = file.path(output_dirs$statistics, "Module1_Summary.csv"),
           row.names = FALSE)
  
  message("\n=== 模块1完成总结 ===")
  message(sprintf("总细胞数: %d", ncol(seurat_obj)))
  message(sprintf("总基因数: %d", nrow(seurat_obj)))
  message(sprintf("聚类数: %d", length(unique(seurat_obj$seurat_clusters))))
  message(sprintf("细胞类型数: %d", length(unique(seurat_obj$cell_type))))
  
  return(seurat_obj)
}

# ==============================================================================
# Module 1: 主执行函数
# ==============================================================================

run_module1_analysis <- function(workDir, data_file, marker_file, params,
                                 output_dirs, tracker, cluster_annotation = NULL) {
  
  tracker$set_module("Module 1: Core Analysis")
  
  # 1. 读取数据
  if (file.exists(data_file)) {
    message("从RDS文件加载数据...")
    single_cell_data <- safe_execute(readRDS, "load_rds_data", tracker,
                                     file = data_file)
    counts_data <- list(counts = single_cell_data)
  } else {
    counts_data <- safe_execute(read_and_process_10x_data, "read_10x_data",
                               tracker, workDir = workDir,
                               output_dirs = output_dirs)
  }
  
  if (is.null(counts_data)) {
    stop("数据读取失败")
  }
  
  # 2. 创建Seurat对象
  seurat_obj <- safe_execute(create_seurat_with_qc, "create_seurat", tracker,
                             counts_data = counts_data, params = params,
                             output_dirs = output_dirs)
  
  # 3. 归一化和降维
  seurat_obj <- safe_execute(normalize_and_reduce, "normalize_reduce", tracker,
                             seurat_obj = seurat_obj, params = params,
                             output_dirs = output_dirs)
  
  # 4. Marker基因分析
  marker_results <- safe_execute(find_all_markers, "find_markers", tracker,
                                seurat_obj = seurat_obj, params = params,
                                output_dirs = output_dirs)
  
  # 5. 自动注释
  seurat_obj <- safe_execute(automatic_cell_annotation, "auto_annotation", tracker,
                             seurat_obj = seurat_obj, output_dirs = output_dirs)
  
  # 6. 手动注释
  seurat_obj <- safe_execute(manual_cell_annotation, "manual_annotation", tracker,
                             seurat_obj = seurat_obj,
                             cluster_annotation = cluster_annotation,
                             output_dirs = output_dirs)
  
  # 7. 差异分析
  seurat_obj <- safe_execute(perform_differential_expression, "differential_analysis",
                             tracker, seurat_obj = seurat_obj, params = params,
                             output_dirs = output_dirs)
  
  # 8. 比例分析
  safe_execute(analyze_proportions, "proportion_analysis", tracker,
              seurat_obj = seurat_obj, output_dirs = output_dirs)
  
  # 9. 保存结果
  seurat_obj <- safe_execute(save_module1_results, "save_results", tracker,
                             seurat_obj = seurat_obj, output_dirs = output_dirs,
                             workDir = workDir)
  
  message("\n✓✓✓ 模块1：核心分析完成 ✓✓✓\n")
  
  return(list(
    seurat_obj = seurat_obj,
    marker_results = marker_results
  ))
}
