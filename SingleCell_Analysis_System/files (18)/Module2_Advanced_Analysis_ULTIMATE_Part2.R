# ==============================================================================
# Module 2: 高级分析 - Part 2/2
# 包含: 基因可视化, 高级图表, 综合分析
# ==============================================================================

# ==============================================================================
# Module 2.4: 基因可视化
# ==============================================================================

visualize_gene_of_interest <- function(seurat_obj, params, output_dirs) {
  cat("\n=== 模块2.4: 基因可视化 ===\n")
  
  tryCatch({
    gene_name <- params$gene_of_interest
    cat(sprintf("可视化基因: %s\n", gene_name))
    
    # 检查基因是否存在
    if (!(gene_name %in% rownames(seurat_obj))) {
      cat(sprintf("✗ 基因 %s 不在数据集中\n", gene_name))
      return(NULL)
    }
    
    # === 散点图（FeaturePlot）===
    cat("生成散点图...\n")
    p1 <- FeaturePlot(object = seurat_obj, features = gene_name,
                     cols = c("lightgrey", "#D67280"), pt.size = 1.5) +
      ggtitle(sprintf("%s Expression", gene_name)) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "right"
      )
    
    ggsave(file.path(output_dirs$module2_gene_view, "Gene_FeaturePlot.pdf"),
           p1, width = 8, height = 7, dpi = 600)
    ggsave(file.path(output_dirs$module2_gene_view, "Gene_FeaturePlot.png"),
           p1, width = 8, height = 7, dpi = 300)
    
    # === 小提琴图 ===
    cat("生成小提琴图...\n")
    if ("Type" %in% colnames(seurat_obj@meta.data)) {
      p2 <- VlnPlot(object = seurat_obj, features = gene_name, split.by = "Type",
                   pt.size = 0.1, cols = c("#A8D8E8", "#D67280")) +
        ggtitle(sprintf("%s Expression by Group", gene_name)) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title.x = element_blank()
        )
    } else {
      p2 <- VlnPlot(object = seurat_obj, features = gene_name,
                   pt.size = 0.1) +
        ggtitle(sprintf("%s Expression by Cell Type", gene_name)) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    }
    
    ggsave(file.path(output_dirs$module2_gene_view, "Gene_ViolinPlot.pdf"),
           p2, width = 10, height = 6, dpi = 600)
    
    # === 点图 ===
    cat("生成点图...\n")
    p3 <- DotPlot(object = seurat_obj, features = gene_name) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      ggtitle(sprintf("%s Expression Across Cell Types", gene_name)) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.title.x = element_blank()
      )
    
    ggsave(file.path(output_dirs$module2_gene_view, "Gene_DotPlot.pdf"),
           p3, width = 8, height = 6, dpi = 600)
    
    # === 如果有Type分组，生成分组比较图 ===
    if ("Type" %in% colnames(seurat_obj@meta.data)) {
      cat("生成分组比较图...\n")
      p4 <- FeaturePlot(object = seurat_obj, features = gene_name,
                       split.by = "Type", cols = c("lightgrey", "#D67280"),
                       pt.size = 1) +
        ggtitle(sprintf("%s Expression by Group", gene_name)) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
      
      ggsave(file.path(output_dirs$module2_gene_view, "Gene_FeaturePlot_Split.pdf"),
             p4, width = 12, height = 5, dpi = 600)
    }
    
    # === 岭线图（如果有ggridges包）===
    if (requireNamespace("ggridges", quietly = TRUE)) {
      cat("生成岭线图...\n")
      library(ggridges)
      
      # 修复GetAssayData
      plot_data <- data.frame(
        Expression = GetAssayData(seurat_obj, layer = "data")[gene_name, ],  # ✓ 修复
        CellType = Idents(seurat_obj)
      )
      
      p5 <- ggplot(plot_data, aes(x = Expression, y = CellType, fill = CellType)) +
        geom_density_ridges(alpha = 0.7) +
        scale_fill_manual(values = colorRampPalette(c("#E8F4F8", "#A8D8E8",
                                                      "#7BA7BC", "#D67280"))(
                                                        length(unique(plot_data$CellType)))) +
        theme_minimal() +
        labs(title = sprintf("%s Expression Distribution", gene_name),
             x = "Expression Level", y = "Cell Type") +
        theme(
          legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title = element_text(face = "bold", size = 12)
        )
      
      ggsave(file.path(output_dirs$module2_gene_view, "Gene_RidgePlot.pdf"),
             p5, width = 10, height = 8, dpi = 600)
    }
    
    # === 热图（基因表达在不同细胞类型）===
    cat("生成热图...\n")
    tryCatch({
      p6 <- DoHeatmap(seurat_obj, features = gene_name, size = 4) +
        scale_fill_gradientn(colors = c("#7BA7BC", "white", "#D67280")) +
        labs(title = sprintf("%s Expression Heatmap", gene_name)) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
      
      ggsave(file.path(output_dirs$module2_gene_view, "Gene_Heatmap.pdf"),
             p6, width = 10, height = 6, dpi = 600)
    }, error = function(e) {
      cat("⚠ 热图生成失败\n")
    })
    
    cat("✓ 基因可视化完成\n")
    cat(sprintf("  - 散点图: Gene_FeaturePlot.pdf/png\n"))
    cat(sprintf("  - 小提琴图: Gene_ViolinPlot.pdf\n"))
    cat(sprintf("  - 点图: Gene_DotPlot.pdf\n"))
    
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("✗ 基因可视化失败: %s\n", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.5: 高级可视化（综合图表）
# ==============================================================================

create_advanced_visualizations <- function(seurat_obj, output_dirs) {
  cat("\n=== 模块2.5: 高级综合可视化 ===\n")
  
  tryCatch({
    # === 1. 细胞类型标记基因气泡图 ===
    cat("生成细胞类型标记基因气泡图...\n")
    
    cell_type_markers <- list(
      "T cells" = c("CD3D", "CD3E", "CD4", "CD8A", "IL7R", "TRAC"),
      "NK cells" = c("GNLY", "NKG7", "FGFBP2", "GZMB", "GZMA", "PRF1"),
      "B cells" = c("CD19", "MS4A1", "CD79A", "CD79B", "IGHM", "IGHD"),
      "Monocytes" = c("CD14", "LYZ", "S100A8", "S100A9", "CST3", "FCGR1A"),
      "Dendritic Cells" = c("FCER1A", "CD1C", "CLEC9A", "IDO1", "THBD"),
      "Plasma cells" = c("MZB1", "IGHG1", "IGHA1", "IGHM", "POU2AF1", "XBP1")
    )
    
    # 获取存在的标记基因
    valid_markers <- list()
    for (celltype in names(cell_type_markers)) {
      valid_genes <- intersect(cell_type_markers[[celltype]], rownames(seurat_obj))
      if (length(valid_genes) > 0) {
        valid_markers[[celltype]] <- valid_genes
      }
    }
    
    if (length(valid_markers) > 0) {
      all_markers <- unique(unlist(valid_markers))
      
      # DotPlot
      p1 <- DotPlot(seurat_obj, features = all_markers, dot.scale = 8, scale = TRUE,
                   cols = c("#E8F4F8", "#7BA7BC", "#D67280")) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.title.x = element_blank()
        ) +
        labs(title = "Cell Type Marker Expression Profile")
      
      ggsave(file.path(output_dirs$module2_advanced_plots, "CellType_Markers_DotPlot.pdf"),
             p1, width = 14, height = 9, dpi = 600)
      
      cat("✓ 细胞类型标记基因气泡图已生成\n")
    }
    
    # === 2. 细胞比例堆叠图 ===
    if ("Type" %in% colnames(seurat_obj@meta.data)) {
      cat("生成细胞比例堆叠图...\n")
      
      prop_data <- seurat_obj@meta.data %>%
        group_by(Type, seurat_clusters) %>%
        summarise(count = n(), .groups = 'drop') %>%
        group_by(Type) %>%
        mutate(proportion = count / sum(count) * 100)
      
      p2 <- ggplot(prop_data, aes(x = Type, y = proportion, fill = seurat_clusters)) +
        geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
        scale_fill_manual(values = colorRampPalette(c("#E8F4F8", "#A8D8E8",
                                                      "#7BA7BC", "#5B8FA4", "#D67280"))(
                                                        length(unique(prop_data$seurat_clusters)))) +
        theme_minimal() +
        labs(
          title = "Cell Type Proportions by Group",
          x = "Group",
          y = "Proportion (%)",
          fill = "Cluster"
        ) +
        theme(
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          axis.title = element_text(face = "bold", size = 12),
          legend.title = element_text(face = "bold")
        )
      
      ggsave(file.path(output_dirs$module2_advanced_plots, "Cell_Proportions_Stacked.pdf"),
             p2, width = 10, height = 7, dpi = 600)
      
      cat("✓ 细胞比例堆叠图已生成\n")
    }
    
    # === 3. UMAP with Density ===
    cat("生成UMAP密度图...\n")
    if (requireNamespace("ggpointdensity", quietly = TRUE)) {
      library(ggpointdensity)
      
      umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
      
      p3 <- ggplot(umap_coords, aes(x = UMAP_1, y = UMAP_2)) +
        geom_pointdensity(size = 1) +
        scale_color_viridis_c() +
        theme_minimal() +
        labs(title = "UMAP with Cell Density",
             color = "Density") +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.position = "right"
        )
      
      ggsave(file.path(output_dirs$module2_advanced_plots, "UMAP_Density.pdf"),
             p3, width = 9, height = 7, dpi = 600)
      
      cat("✓ UMAP密度图已生成\n")
    }
    
    # === 4. 细胞周期评分（如果有细胞周期基因）===
    cat("生成细胞周期图...\n")
    tryCatch({
      # 加载细胞周期基因
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
      
      # 计算评分
      seurat_obj <- CellCycleScoring(seurat_obj,
                                     s.features = s.genes,
                                     g2m.features = g2m.genes,
                                     set.ident = FALSE)
      
      # 可视化
      p4 <- DimPlot(seurat_obj, reduction = "umap", group.by = "Phase",
                   cols = c("G1" = "#A8D8E8", "S" = "#D67280", "G2M" = "#7BA7BC")) +
        labs(title = "Cell Cycle Phase Distribution") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14))
      
      ggsave(file.path(output_dirs$module2_advanced_plots, "Cell_Cycle_UMAP.pdf"),
             p4, width = 9, height = 7, dpi = 600)
      
      # 统计
      phase_counts <- table(seurat_obj$Phase)
      phase_df <- data.frame(
        Phase = names(phase_counts),
        Count = as.numeric(phase_counts),
        Percentage = round(as.numeric(phase_counts) / sum(phase_counts) * 100, 1)
      )
      
      write.csv(phase_df,
               file.path(output_dirs$module2_advanced_plots, "Cell_Cycle_Statistics.csv"),
               row.names = FALSE)
      
      cat("✓ 细胞周期图已生成\n")
      
    }, error = function(e) {
      cat("⚠ 细胞周期分析失败（可能缺少细胞周期基因）\n")
    })
    
    # === 5. QC指标综合图 ===
    cat("生成QC指标综合图...\n")
    qc_data <- seurat_obj@meta.data
    
    p5_1 <- ggplot(qc_data, aes(x = nCount_RNA, y = nFeature_RNA)) +
      geom_point(alpha = 0.5, color = "#7BA7BC") +
      geom_smooth(method = "lm", color = "#D67280") +
      theme_minimal() +
      labs(title = "UMI vs Feature Count",
           x = "UMI Count", y = "Feature Count") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    p5_2 <- ggplot(qc_data, aes(x = nCount_RNA, y = percent.mt)) +
      geom_point(alpha = 0.5, color = "#7BA7BC") +
      geom_smooth(method = "lm", color = "#D67280") +
      theme_minimal() +
      labs(title = "UMI vs Mitochondrial %",
           x = "UMI Count", y = "Mitochondrial %") +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    p5_combined <- p5_1 + p5_2
    ggsave(file.path(output_dirs$module2_advanced_plots, "QC_Metrics_Combined.pdf"),
           p5_combined, width = 14, height = 6, dpi = 600)
    
    cat("✓ QC指标综合图已生成\n")
    
    # === 6. 细胞类型比例饼图 ===
    cat("生成细胞类型饼图...\n")
    celltype_counts <- table(Idents(seurat_obj))
    celltype_df <- data.frame(
      CellType = names(celltype_counts),
      Count = as.numeric(celltype_counts),
      Percentage = round(as.numeric(celltype_counts) / sum(celltype_counts) * 100, 1)
    )
    celltype_df$Label <- paste0(celltype_df$CellType, "\n",
                                celltype_df$Percentage, "%")
    
    p6 <- ggplot(celltype_df, aes(x = "", y = Count, fill = CellType)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = colorRampPalette(c("#E8F4F8", "#A8D8E8",
                                                    "#7BA7BC", "#5B8FA4", "#D67280"))(
                                                      nrow(celltype_df))) +
      geom_text(aes(label = Percentage),
               position = position_stack(vjust = 0.5),
               size = 4, fontface = "bold") +
      theme_void() +
      labs(title = "Cell Type Distribution",
           fill = "Cell Type") +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "right"
      )
    
    ggsave(file.path(output_dirs$module2_advanced_plots, "CellType_PieChart.pdf"),
           p6, width = 10, height = 7, dpi = 600)
    
    cat("✓ 细胞类型饼图已生成\n")
    
    cat("✓ 高级综合可视化完成\n")
    cat("  - 标记基因气泡图\n")
    cat("  - 细胞比例图\n")
    cat("  - UMAP密度图\n")
    cat("  - 细胞周期图\n")
    cat("  - QC综合图\n")
    cat("  - 细胞类型饼图\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("✗ 高级可视化失败: %s\n", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.6: 运行所有Module2分析
# ==============================================================================

run_all_module2_analyses <- function(seurat_obj, params, output_dirs) {
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("开始Module 2: 高级分析\n")
  cat(strrep("=", 80), "\n\n")
  
  results <- list()
  
  # 2.1 CellChat
  if (params$run_cellchat) {
    results$cellchat <- run_cellchat_analysis(seurat_obj, params, output_dirs)
  }
  
  # 2.2 虚拟敲除
  if (params$run_knockout) {
    results$knockout <- run_virtual_knockout(seurat_obj, params, output_dirs)
  }
  
  # 2.3 GO富集
  if (params$run_go_enrichment) {
    results$go <- run_go_enrichment(seurat_obj, params, output_dirs)
  }
  
  # 2.4 基因可视化
  if (params$run_gene_visualization) {
    results$gene_viz <- visualize_gene_of_interest(seurat_obj, params, output_dirs)
  }
  
  # 2.5 高级可视化
  if (params$run_advanced_plots) {
    results$advanced_plots <- create_advanced_visualizations(seurat_obj, output_dirs)
  }
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("✓ Module 2 高级分析全部完成\n")
  cat(strrep("=", 80), "\n")
  
  return(results)
}

cat("✓ Module 2 (Part 2/2) 已加载：基因可视化, 高级图表\n")
cat("✓ Module 2 完整版已准备就绪\n")
