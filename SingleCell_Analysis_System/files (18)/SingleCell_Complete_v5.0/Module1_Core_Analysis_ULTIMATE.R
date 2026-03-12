# ============================================================================
# Module 1: 核心分析 - 终极完整版
# ✓ 修复所有GetAssayData错误
# ✓ 内存优化
# ✓ 完整断点续传
# ✓ 所有缺失图表
# ✓ SCI级别可视化
# ============================================================================

# ============================================================================
# 函数1: 读取和处理10X数据
# ============================================================================

read_and_process_10x_data <- function(data_dir, output_dirs) {
  checkpoint_file <- "checkpoint_read_data.rds"
  
  if (file.exists(checkpoint_file)) {
    cat("⚡ 检测到断点，加载已保存的数据...\n")
    return(readRDS(checkpoint_file))
  }
  
  tryCatch({
    cat("=== 步骤1.1: 读取10X数据 ===\n")
    
    all_dirs <- list.dirs(data_dir, recursive = FALSE)
    
    valid_dirs <- all_dirs[sapply(all_dirs, function(d) {
      files <- list.files(d)
      any(grepl("barcodes|genes|matrix", files, ignore.case = TRUE))
    })]
    
    if (length(valid_dirs) == 0) {
      stop("未找到有效的10X数据目录")
    }
    
    cat(sprintf("找到 %d 个数据目录\n", length(valid_dirs)))
    
    data_list <- list()
    pb <- txtProgressBar(min = 0, max = length(valid_dirs), style = 3)
    
    for (i in seq_along(valid_dirs)) {
      dir_path <- valid_dirs[i]
      sample_name <- basename(dir_path)
      
      tryCatch({
        data <- Read10X(data.dir = dir_path)
        data_list[[sample_name]] <- data
        cat(sprintf("\n✓ %s: %d genes × %d cells\n", 
                   sample_name, nrow(data), ncol(data)))
      }, error = function(e) {
        cat(sprintf("\n✗ %s 读取失败: %s\n", sample_name, e$message))
      })
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    if (length(data_list) == 0) {
      stop("没有成功读取任何数据")
    }
    
    # 合并数据
    if (length(data_list) > 1) {
      all_genes <- Reduce(intersect, lapply(data_list, rownames))
      combined_data <- do.call(cbind, lapply(data_list, function(x) x[all_genes, ]))
    } else {
      combined_data <- data_list[[1]]
    }
    
    result <- list(counts = combined_data)
    
    # 保存断点
    saveRDS(result, checkpoint_file)
    cat(sprintf("\n✓ 合并完成: %d genes × %d cells\n", 
               nrow(combined_data), ncol(combined_data)))
    
    return(result)
    
  }, error = function(e) {
    stop(sprintf("读取数据失败: %s", e$message))
  })
}

# ============================================================================
# 函数2: 创建Seurat对象并质控 - 完整图表版
# ============================================================================

create_seurat_with_qc <- function(counts_data, params, output_dirs) {
  checkpoint_file <- "checkpoint_seurat_qc.rds"
  
  if (file.exists(checkpoint_file)) {
    cat("⚡ 加载已保存的Seurat对象...\n")
    return(readRDS(checkpoint_file))
  }
  
  tryCatch({
    cat("\n=== 步骤1.2: 创建Seurat对象并质控 ===\n")
    
    pbmc <- CreateSeuratObject(
      counts = counts_data$counts,
      project = "SingleCell",
      min.cells = params$min_cells,
      min.features = params$min_features
    )
    
    cat(sprintf("初始细胞数: %d\n", ncol(pbmc)))
    
    # 计算QC指标
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    pbmc[["percent.rb"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
    
    qc_dir <- output_dirs$module1_qc
    
    # ========== SCI级别QC可视化 ==========
    
    # 1. 小提琴图 - 柔和配色
    cat("生成QC小提琴图...\n")
    p1 <- VlnPlot(pbmc, 
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, 
                  pt.size = 0.1,
                  cols = c("#B8E6F0")) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none"
      )
    
    ggsave(file.path(qc_dir, "01_QC_ViolinPlot.pdf"), p1,
           width = 12, height = 5, dpi = 600)
    ggsave(file.path(qc_dir, "01_QC_ViolinPlot.png"), p1,
           width = 12, height = 5, dpi = 300)
    
    # 2. 散点图 - 特征关联
    cat("生成QC散点图...\n")
    p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",
                         pt.size = 0.5) +
      geom_smooth(method = "lm", color = "#D67280", se = TRUE, alpha = 0.2) +
      theme_minimal() +
      labs(title = "UMI Count vs Mitochondrial Content") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    p3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                         pt.size = 0.5) +
      geom_smooth(method = "lm", color = "#7BA7BC", se = TRUE, alpha = 0.2) +
      theme_minimal() +
      labs(title = "UMI Count vs Feature Count") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    p_scatter <- p2 + p3
    ggsave(file.path(qc_dir, "02_QC_FeatureScatter.pdf"), p_scatter,
           width = 12, height = 5, dpi = 600)
    
    # 3. 密度图
    cat("生成QC密度图...\n")
    qc_data <- pbmc@meta.data
    
    p4 <- ggplot(qc_data, aes(x = nFeature_RNA)) +
      geom_density(fill = "#A8D8E8", alpha = 0.7, color = "#4A90A4", size = 1) +
      geom_vline(xintercept = params$min_features, linetype = "dashed", 
                 color = "#D64545", size = 1) +
      geom_vline(xintercept = params$max_features, linetype = "dashed", 
                 color = "#D64545", size = 1) +
      theme_minimal() +
      labs(title = "Feature Distribution", 
           x = "Number of Features", 
           y = "Density") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)
      )
    
    ggsave(file.path(qc_dir, "03_QC_Feature_Density.pdf"), p4,
           width = 8, height = 6, dpi = 600)
    
    # 4. 质控前后统计对比
    cat("生成统计对比...\n")
    qc_stats_before <- data.frame(
      Metric = c("Total Cells", "Mean Features", "Mean Counts", "Mean MT%"),
      Before = c(
        ncol(pbmc),
        round(mean(pbmc$nFeature_RNA), 2),
        round(mean(pbmc$nCount_RNA), 2),
        round(mean(pbmc$percent.mt), 2)
      )
    )
    
    # 质控过滤
    pbmc <- subset(pbmc, subset = 
                     nFeature_RNA > params$min_features & 
                     nFeature_RNA < params$max_features & 
                     percent.mt < params$mt_cutoff)
    
    qc_stats_after <- data.frame(
      After = c(
        ncol(pbmc),
        round(mean(pbmc$nFeature_RNA), 2),
        round(mean(pbmc$nCount_RNA), 2),
        round(mean(pbmc$percent.mt), 2)
      )
    )
    
    qc_comparison <- cbind(qc_stats_before, qc_stats_after)
    qc_comparison$Change <- sprintf("%.1f%%", 
                                    (qc_comparison$After / qc_comparison$Before - 1) * 100)
    
    write.csv(qc_comparison, 
              file.path(qc_dir, "QC_Statistics_Comparison.csv"),
              row.names = FALSE)
    
    # 5. 可视化统计对比
    qc_long <- tidyr::pivot_longer(qc_comparison[,1:3], 
                                    cols = c("Before", "After"),
                                    names_to = "Stage",
                                    values_to = "Value")
    
    p5 <- ggplot(qc_long, aes(x = Metric, y = Value, fill = Stage)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.7) +
      scale_fill_manual(values = c("Before" = "#E8B4B8", "After" = "#88CAE0")) +
      theme_minimal() +
      labs(title = "QC Statistics: Before vs After Filtering",
           x = "", y = "Value") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        legend.title = element_text(face = "bold")
      )
    
    ggsave(file.path(qc_dir, "04_QC_Statistics_Comparison.pdf"), p5,
           width = 10, height = 7, dip = 600)
    
    cat(sprintf("✓ 质控完成，保留 %d 细胞 (%.1f%%)\n", 
                ncol(pbmc),
                ncol(pbmc) / nrow(qc_stats_before) * 100))
    cat("✓ 所有QC图表已保存\n")
    
    # 保存断点
    saveRDS(pbmc, checkpoint_file)
    
    return(pbmc)
    
  }, error = function(e) {
    stop(sprintf("创建Seurat对象失败: %s", e$message))
  })
}

# ============================================================================
# 函数3: 归一化和降维 - 完整图表版
# ============================================================================

normalize_and_reduce <- function(pbmc, params, output_dirs) {
  checkpoint_file <- "checkpoint_normalized.rds"
  
  if (file.exists(checkpoint_file)) {
    cat("⚡ 加载已归一化的数据...\n")
    return(readRDS(checkpoint_file))
  }
  
  tryCatch({
    cat("\n=== 步骤1.3: 数据归一化和降维 ===\n")
    
    cat("1/7: 数据归一化...\n")
    pbmc <- NormalizeData(pbmc)
    
    cat("2/7: 识别高变基因...\n")
    pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", 
                                 nfeatures = params$n_features)
    
    # 高变基因可视化
    cat("生成高变基因图...\n")
    top10 <- head(VariableFeatures(pbmc), 10)
    
    p1 <- VariableFeaturePlot(pbmc) +
      theme_minimal() +
      labs(title = "Highly Variable Features") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "top"
      )
    
    p2 <- LabelPoints(plot = p1, points = top10, repel = TRUE,
                     xnudge = 0, ynudge = 0) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    ggsave(file.path(output_dirs$module1_basic, "01_VariableFeatures.pdf"),
           p2, width = 10, height = 7, dpi = 600)
    ggsave(file.path(output_dirs$module1_basic, "01_VariableFeatures.png"),
           p2, width = 10, height = 7, dpi = 300)
    
    cat("3/7: 数据标准化...\n")
    all.genes <- rownames(pbmc)
    pbmc <- ScaleData(pbmc, features = all.genes)
    
    cat("4/7: PCA降维...\n")
    pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
    
    # PCA可视化
    cat("生成PCA图...\n")
    p3 <- DimPlot(pbmc, reduction = "pca") +
      theme_minimal() +
      labs(title = "PCA Plot") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
    
    ggsave(file.path(output_dirs$module1_basic, "02_PCA_Plot.pdf"),
           p3, width = 8, height = 7, dpi = 600)
    
    # PCA热图
    cat("生成PCA热图...\n")
    pdf(file.path(output_dirs$module1_basic, "03_PCA_Heatmap.pdf"), 
        width = 12, height = 10)
    DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE, ncol = 3)
    dev.off()
    
    # Elbow plot
    cat("生成Elbow图...\n")
    p4 <- ElbowPlot(pbmc, ndims = 50) +
      theme_minimal() +
      geom_vline(xintercept = params$n_pcs, linetype = "dashed", 
                 color = "red", size = 1) +
      labs(title = "PC Selection: Elbow Plot") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
    
    ggsave(file.path(output_dirs$module1_basic, "04_ElbowPlot.pdf"),
           p4, width = 8, height = 6, dpi = 600)
    
    cat("5/7: 构建邻居图...\n")
    pbmc <- FindNeighbors(pbmc, dims = 1:params$n_pcs)
    
    cat("6/7: 聚类分析...\n")
    pbmc <- FindClusters(pbmc, resolution = params$resolution)
    
    cat("7/7: UMAP降维...\n")
    pbmc <- RunUMAP(pbmc, dims = 1:params$n_pcs)
    
    # UMAP可视化 - 多种风格
    cat("生成UMAP图...\n")
    
    # 基础UMAP
    colors_sci <- colorRampPalette(c("#E8F4F8", "#A8D8E8", "#7BA7BC", 
                                     "#5B8FA4", "#D67280"))(
                                       length(unique(pbmc$seurat_clusters)))
    
    p5 <- DimPlot(pbmc, reduction = "umap", label = TRUE, 
                  label.size = 6, pt.size = 0.5,
                  cols = colors_sci) +
      theme_minimal() +
      labs(title = "UMAP Clustering") +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.text = element_text(size = 10)
      )
    
    ggsave(file.path(output_dirs$module1_clustering, "01_UMAP_Clusters.pdf"),
           p5, width = 10, height = 8, dpi = 600)
    ggsave(file.path(output_dirs$module1_clustering, "01_UMAP_Clusters.png"),
           p5, width = 10, height = 8, dpi = 300)
    
    # UMAP with split by orig.ident (if exists)
    if ("orig.ident" %in% colnames(pbmc@meta.data)) {
      p6 <- DimPlot(pbmc, reduction = "umap", split.by = "orig.ident",
                    label = TRUE, label.size = 4, pt.size = 0.3,
                    cols = colors_sci) +
        theme_minimal() +
        labs(title = "UMAP by Sample") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      ggsave(file.path(output_dirs$module1_clustering, "02_UMAP_BySample.pdf"),
             p6, width = 14, height = 7, dpi = 600)
    }
    
    cat("✓ 归一化和降维完成\n")
    cat("✓ 所有降维图表已保存\n")
    
    # 保存断点
    saveRDS(pbmc, checkpoint_file)
    
    return(pbmc)
    
  }, error = function(e) {
    stop(sprintf("归一化和降维失败: %s", e$message))
  })
}

# ============================================================================
# 函数4: Marker基因分析 - 内存优化+完整图表
# ============================================================================

find_all_markers <- function(pbmc, params, output_dirs) {
  checkpoint_file <- "checkpoint_markers.rds"
  
  if (file.exists(checkpoint_file)) {
    cat("⚡ 加载已计算的marker基因...\n")
    return(readRDS(checkpoint_file))
  }
  
  tryCatch({
    cat("\n=== 步骤1.4: Marker基因分析 ===\n")
    cat("识别marker基因（内存优化模式）...\n")
    
    # 内存优化：限制每个聚类的细胞数
    all_markers <- FindAllMarkers(
      pbmc,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      max.cells.per.ident = 500  # 限制细胞数以节省内存
    )
    
    cat(sprintf("找到 %d 个显著marker基因\n", nrow(all_markers)))
    
    # 保存完整结果
    write.csv(all_markers,
              file.path(output_dirs$module1_markers, "All_Markers.csv"),
              row.names = FALSE)
    
    # Top markers
    top_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) %>%
      arrange(cluster, desc(avg_log2FC))
    
    write.csv(top_markers,
              file.path(output_dirs$module1_markers, "Top10_Markers_Per_Cluster.csv"),
              row.names = FALSE)
    
    # ========== 可视化：Marker基因热图 ==========
    
    cat("生成Marker基因热图...\n")
    top5_markers <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = 5, wt = avg_log2FC)
    
    tryCatch({
      # 修复GetAssayData错误 - 使用layer参数
      expr_data <- GetAssayData(pbmc, assay = "RNA", layer = "scale.data")
      
      genes_to_plot <- unique(top5_markers$gene)
      genes_to_plot <- intersect(genes_to_plot, rownames(expr_data))
      
      if (length(genes_to_plot) > 50) {
        genes_to_plot <- head(genes_to_plot, 50)
      }
      
      if (length(genes_to_plot) > 0) {
        p_heat <- DoHeatmap(
          pbmc,
          features = genes_to_plot,
          size = 3,
          angle = 0,
          hjust = 0.5
        ) +
          scale_fill_gradientn(
            colors = c("#7BA7BC", "white", "#D67280"),
            na.value = "white"
          ) +
          labs(title = "Top Marker Genes Heatmap") +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            axis.text.y = element_text(size = 7)
          )
        
        ggsave(file.path(output_dirs$module1_heatmaps, "01_Markers_Heatmap.pdf"),
               p_heat, width = 12, height = 14, dpi = 600)
        
        cat("✓ Marker基因热图已保存\n")
      }
    }, error = function(e) {
      cat(sprintf("⚠ 热图生成失败: %s\n", e$message))
    })
    
    # ========== 可视化：Dot plot ==========
    
    cat("生成Dot plot...\n")
    tryCatch({
      top3_per_cluster <- all_markers %>%
        group_by(cluster) %>%
        top_n(n = 3, wt = avg_log2FC) %>%
        pull(gene) %>%
        unique()
      
      if (length(top3_per_cluster) > 30) {
        top3_per_cluster <- head(top3_per_cluster, 30)
      }
      
      p_dot <- DotPlot(pbmc, features = top3_per_cluster) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
        scale_color_gradient2(low = "#E8F4F8", mid = "white", high = "#D67280") +
        labs(title = "Marker Genes Expression") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
      
      ggsave(file.path(output_dirs$module1_markers, "02_Markers_DotPlot.pdf"),
             p_dot, width = 14, height = 8, dpi = 600)
      
      cat("✓ Marker基因Dot plot已保存\n")
    }, error = function(e) {
      cat(sprintf("⚠ Dot plot生成失败: %s\n", e$message))
    })
    
    # ========== 可视化：Feature plots for top markers ==========
    
    cat("生成Top marker特征图...\n")
    tryCatch({
      top_overall <- all_markers %>%
        group_by(cluster) %>%
        top_n(n = 1, wt = avg_log2FC) %>%
        pull(gene) %>%
        unique() %>%
        head(9)  # 3x3 grid
      
      p_features <- FeaturePlot(pbmc, features = top_overall, ncol = 3,
                                pt.size = 0.5) &
        scale_color_gradientn(colors = c("lightgrey", "#A8D8E8", "#D67280")) &
        theme_minimal() &
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
      
      ggsave(file.path(output_dirs$module1_markers, "03_Top_Markers_FeaturePlot.pdf"),
             p_features, width = 12, height = 12, dpi = 600)
      
      cat("✓ Top marker特征图已保存\n")
    }, error = function(e) {
      cat(sprintf("⚠ 特征图生成失败: %s\n", e$message))
    })
    
    cat("✓ Marker基因分析完成\n")
    
    # 保存断点
    saveRDS(all_markers, checkpoint_file)
    
    return(all_markers)
    
  }, error = function(e) {
    stop(sprintf("Marker基因分析失败: %s", e$message))
  })
}

# ============================================================================
# 函数5: 自动细胞类型注释
# ============================================================================

automatic_cell_annotation <- function(pbmc, output_dirs) {
  checkpoint_file <- "checkpoint_annotation.rds"
  
  if (file.exists(checkpoint_file)) {
    cat("⚡ 加载已有的注释...\n")
    return(readRDS(checkpoint_file))
  }
  
  tryCatch({
    cat("\n=== 步骤1.5: 自动细胞类型注释 ===\n")
    cat("加载参考数据库...\n")
    
    if (!requireNamespace("SingleR", quietly = TRUE)) {
      BiocManager::install("SingleR")
    }
    if (!requireNamespace("celldex", quietly = TRUE)) {
      BiocManager::install("celldex")
    }
    
    library(SingleR)
    library(celldex)
    
    ref <- celldex::HumanPrimaryCellAtlasData()
    
    sce <- as.SingleCellExperiment(pbmc)
    
    pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
    
    pbmc$singler_labels <- pred$labels
    
    # 可视化
    p1 <- DimPlot(pbmc, reduction = "umap", group.by = "singler_labels",
                  label = TRUE, repel = TRUE, label.size = 3, pt.size = 0.5) +
      theme_minimal() +
      labs(title = "Automatic Cell Type Annotation (SingleR)") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "bottom"
      )
    
    ggsave(file.path(output_dirs$module1_annotation, "01_Auto_Annotation_UMAP.pdf"),
           p1, width = 12, height = 10, dpi = 600)
    ggsave(file.path(output_dirs$module1_annotation, "01_Auto_Annotation_UMAP.png"),
           p1, width = 12, height = 10, dpi = 300)
    
    # 保存注释映射
    annotation_table <- table(pbmc$seurat_clusters, pbmc$singler_labels)
    write.csv(annotation_table,
              file.path(output_dirs$module1_annotation, "Cluster_CellType_Mapping.csv"))
    
    # 堆叠柱状图显示每个聚类的细胞类型分布
    mapping_df <- as.data.frame(annotation_table)
    colnames(mapping_df) <- c("Cluster", "CellType", "Count")
    
    p2 <- ggplot(mapping_df, aes(x = Cluster, y = Count, fill = CellType)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_fill_manual(values = colorRampPalette(
        c("#E8F4F8", "#A8D8E8", "#7BA7BC", "#D67280"))(
          length(unique(mapping_df$CellType)))) +
      theme_minimal() +
      labs(title = "Cell Type Composition per Cluster",
           y = "Proportion") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    
    ggsave(file.path(output_dirs$module1_annotation, "02_CellType_Composition.pdf"),
           p2, width = 12, height = 7, dpi = 600)
    
    cat("✓ 自动注释完成\n")
    
    # 保存断点
    saveRDS(pbmc, checkpoint_file)
    
    return(pbmc)
    
  }, error = function(e) {
    cat(sprintf("⚠ 自动注释失败: %s\n", e$message))
    return(pbmc)
  })
}

# ============================================================================
# 函数6: 手动细胞类型注释
# ============================================================================

manual_cell_annotation <- function(pbmc, annotation_map, output_dirs) {
  tryCatch({
    cat("\n=== 步骤1.6: 应用手动细胞类型注释 ===\n")
    
    current_clusters <- levels(pbmc$seurat_clusters)
    
    # 填充缺失的注释
    for (cluster in current_clusters) {
      if (!cluster %in% names(annotation_map)) {
        annotation_map[cluster] <- paste0("Cluster_", cluster)
      }
    }
    
    new.cluster.ids <- annotation_map[current_clusters]
    names(new.cluster.ids) <- current_clusters
    pbmc <- RenameIdents(pbmc, new.cluster.ids)
    pbmc$cell_type <- Idents(pbmc)
    
    cat(sprintf("✓ 已应用 %d 个聚类的手动注释\n", length(current_clusters)))
    
    # 可视化
    colors_celltype <- colorRampPalette(c("#E8F4F8", "#A8D8E8", "#7BA7BC", 
                                          "#5B8FA4", "#D67280"))(
                                            length(unique(pbmc$cell_type)))
    
    p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, repel = TRUE,
                  label.size = 4, pt.size = 0.5, cols = colors_celltype) +
      theme_minimal() +
      labs(title = "Manual Cell Type Annotation") +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.text = element_text(size = 10)
      )
    
    ggsave(file.path(output_dirs$module1_annotation, "03_Manual_Annotation_UMAP.pdf"),
           p1, width = 12, height = 9, dpi = 600)
    ggsave(file.path(output_dirs$module1_annotation, "03_Manual_Annotation_UMAP.png"),
           p1, width = 12, height = 9, dpi = 300)
    
    # 保存注释映射
    mapping_df <- data.frame(
      Cluster = names(new.cluster.ids),
      CellType = as.character(new.cluster.ids),
      stringsAsFactors = FALSE
    )
    
    write.csv(mapping_df,
              file.path(output_dirs$module1_annotation, "Final_Cell_Type_Mapping.csv"),
              row.names = FALSE)
    
    cat("✓ 细胞类型注释完成\n")
    
    return(pbmc)
    
  }, error = function(e) {
    stop(sprintf("手动注释失败: %s", e$message))
  })
}

cat("✓ Module 1 终极完整版已加载\n")
cat("  - ✓ 修复所有GetAssayData错误\n")
cat("  - ✓ 内存优化\n")
cat("  - ✓ 完整断点续传\n")
cat("  - ✓ 所有可视化图表\n")
cat("  - ✓ SCI级别配色\n")
