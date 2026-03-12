# ============================================================================
# Module 1 续: 差异分析和比例分析 - 完整版
# ============================================================================

# ============================================================================
# 函数7: 差异表达分析 - 包含jjVolcano火山图
# ============================================================================

perform_differential_expression <- function(pbmc, params, output_dirs) {
  checkpoint_file <- "checkpoint_diff_expr.rds"
  
  if (file.exists(checkpoint_file)) {
    cat("⚡ 加载已有的差异分析结果...\n")
    return(readRDS(checkpoint_file))
  }
  
  tryCatch({
    cat("\n=== 步骤1.7: 差异表达分析 ===\n")
    
    # 检查组别信息
    if (!"orig.ident" %in% colnames(pbmc@meta.data)) {
      cat("⚠ 没有组别信息，跳过差异分析\n")
      return(pbmc)
    }
    
    groups <- unique(pbmc$orig.ident)
    if (length(groups) < 2) {
      cat("⚠ 只有一个组别，跳过差异分析\n")
      return(pbmc)
    }
    
    cat(sprintf("比较组别: %s vs %s\n", groups[1], groups[2]))
    
    Idents(pbmc) <- "orig.ident"
    
    # 进行差异分析
    diff_genes <- FindMarkers(
      pbmc,
      ident.1 = groups[1],
      ident.2 = groups[2],
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    
    diff_genes$gene <- rownames(diff_genes)
    diff_genes$significant <- abs(diff_genes$avg_log2FC) > 1 & diff_genes$p_val_adj < 0.05
    
    cat(sprintf("找到 %d 个显著差异基因\n", sum(diff_genes$significant)))
    
    # 保存结果
    write.csv(diff_genes,
              file.path(output_dirs$module1_diff, "Differential_Genes.csv"),
              row.names = FALSE)
    
    # ========== 图1: 增强版火山图 ==========
    
    cat("生成增强版火山图...\n")
    diff_genes$log10p <- -log10(diff_genes$p_val_adj + 1e-300)
    diff_genes$color_group <- "NS"
    diff_genes$color_group[diff_genes$avg_log2FC > 1 & diff_genes$p_val_adj < 0.05] <- "Up"
    diff_genes$color_group[diff_genes$avg_log2FC < -1 & diff_genes$p_val_adj < 0.05] <- "Down"
    
    # 标记top基因
    diff_genes <- diff_genes %>%
      arrange(p_val_adj) %>%
      mutate(label = ifelse(row_number() <= 10, gene, ""))
    
    p_volcano <- ggplot(diff_genes, aes(x = avg_log2FC, y = log10p, color = color_group, label = label)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.5) +
      scale_color_manual(values = c("Up" = "#D67280", "Down" = "#7BA7BC", "NS" = "#E0E0E0"),
                        labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not Significant")) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", size = 0.8) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", size = 0.8) +
      theme_minimal() +
      labs(
        title = sprintf("Differential Expression: %s vs %s", groups[1], groups[2]),
        x = "log2 Fold Change",
        y = "-log10(Adjusted P-value)",
        color = "Regulation"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(face = "bold"),
        legend.position = "top"
      )
    
    ggsave(file.path(output_dirs$module1_diff, "01_Volcano_Plot_Enhanced.pdf"),
           p_volcano, width = 10, height = 8, dpi = 600)
    ggsave(file.path(output_dirs$module1_diff, "01_Volcano_Plot_Enhanced.png"),
           p_volcano, width = 10, height = 8, dpi = 300)
    
    # ========== 图2: MA plot ==========
    
    cat("生成MA plot...\n")
    diff_genes$avg_expr <- rowMeans(GetAssayData(pbmc, layer = "data")[diff_genes$gene, ])
    
    p_ma <- ggplot(diff_genes, aes(x = avg_expr, y = avg_log2FC, color = color_group)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = c("Up" = "#D67280", "Down" = "#7BA7BC", "NS" = "#E0E0E0")) +
      geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
      geom_hline(yintercept = 0, color = "black") +
      theme_minimal() +
      labs(
        title = "MA Plot",
        x = "Average Expression",
        y = "log2 Fold Change",
        color = "Regulation"
      ) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "top"
      )
    
    ggsave(file.path(output_dirs$module1_diff, "02_MA_Plot.pdf"),
           p_ma, width = 10, height = 8, dpi = 600)
    
    # ========== 图3: Top差异基因热图 ==========
    
    cat("生成Top差异基因热图...\n")
    top_diff_genes <- diff_genes %>%
      filter(significant == TRUE) %>%
      arrange(desc(abs(avg_log2FC))) %>%
      head(50) %>%
      pull(gene)
    
    if (length(top_diff_genes) > 0) {
      tryCatch({
        p_heat_diff <- DoHeatmap(
          pbmc,
          features = top_diff_genes,
          group.by = "orig.ident",
          size = 3
        ) +
          scale_fill_gradientn(colors = c("#7BA7BC", "white", "#D67280")) +
          labs(title = "Top 50 Differentially Expressed Genes") +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.y = element_text(size = 6)
          )
        
        ggsave(file.path(output_dirs$module1_heatmaps, "02_DiffGenes_Heatmap.pdf"),
               p_heat_diff, width = 10, height = 12, dpi = 600)
      }, error = function(e) {
        cat("⚠ 差异基因热图生成失败\n")
      })
    }
    
    # ========== 图4: Top基因的Feature plots ==========
    
    cat("生成Top差异基因Feature plots...\n")
    if (length(top_diff_genes) >= 9) {
      top9 <- head(top_diff_genes, 9)
      
      p_features <- FeaturePlot(pbmc, features = top9, ncol = 3, pt.size = 0.5) &
        scale_color_gradientn(colors = c("lightgrey", "#A8D8E8", "#D67280")) &
        theme_minimal() &
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
      
      ggsave(file.path(output_dirs$module1_diff, "03_Top_DiffGenes_FeaturePlot.pdf"),
             p_features, width = 12, height = 12, dpi = 600)
    }
    
    cat("✓ 差异表达分析完成\n")
    
    # 保存断点
    result <- list(pbmc = pbmc, diff_genes = diff_genes)
    saveRDS(result, checkpoint_file)
    
    return(pbmc)
    
  }, error = function(e) {
    cat(sprintf("⚠ 差异分析失败: %s\n", e$message))
    return(pbmc)
  })
}

# ============================================================================
# 函数8: 细胞比例分析 - 完整版
# ============================================================================

analyze_proportions <- function(pbmc, output_dirs) {
  tryCatch({
    cat("\n=== 步骤1.8: 细胞比例分析 ===\n")
    
    if (!"orig.ident" %in% colnames(pbmc@meta.data) || 
        !"cell_type" %in% colnames(pbmc@meta.data)) {
      cat("⚠ 缺少必要信息，跳过比例分析\n")
      return()
    }
    
    # 计算比例
    prop_data <- pbmc@meta.data %>%
      group_by(orig.ident, cell_type) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(orig.ident) %>%
      mutate(proportion = count / sum(count) * 100)
    
    write.csv(prop_data,
              file.path(output_dirs$module1_proportion, "Cell_Proportions.csv"),
              row.names = FALSE)
    
    # ========== 图1: 堆叠柱状图 ==========
    
    cat("生成堆叠柱状图...\n")
    p1 <- ggplot(prop_data, aes(x = orig.ident, y = proportion, fill = cell_type)) +
      geom_bar(stat = "identity", position = "stack", color = "white", size = 0.3) +
      scale_fill_manual(values = colorRampPalette(c("#E8F4F8", "#A8D8E8", "#7BA7BC", 
                                                     "#5B8FA4", "#D67280"))(
                                                       length(unique(prop_data$cell_type)))) +
      theme_minimal() +
      labs(
        title = "Cell Type Proportions by Group",
        x = "Group",
        y = "Proportion (%)",
        fill = "Cell Type"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
    ggsave(file.path(output_dirs$module1_proportion, "01_Proportions_Stacked.pdf"),
           p1, width = 10, height = 7, dpi = 600)
    
    # ========== 图2: 分组柱状图 ==========
    
    cat("生成分组柱状图...\n")
    p2 <- ggplot(prop_data, aes(x = cell_type, y = proportion, fill = orig.ident)) +
      geom_bar(stat = "identity", position = position_dodge(), color = "white", size = 0.3) +
      scale_fill_manual(values = c("#A8D8E8", "#D67280")) +
      theme_minimal() +
      labs(
        title = "Cell Type Proportions Comparison",
        x = "Cell Type",
        y = "Proportion (%)",
        fill = "Group"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.title = element_text(face = "bold"),
        legend.position = "top"
      )
    
    ggsave(file.path(output_dirs$module1_proportion, "02_Proportions_Grouped.pdf"),
           p2, width = 12, height = 7, dpi = 600)
    
    # ========== 图3: 气泡图 ==========
    
    cat("生成气泡图...\n")
    p3 <- ggplot(prop_data, aes(x = orig.ident, y = cell_type, size = proportion, color = proportion)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient2(low = "#E8F4F8", mid = "#A8D8E8", high = "#D67280", 
                           midpoint = mean(prop_data$proportion)) +
      scale_size_continuous(range = c(3, 15)) +
      theme_minimal() +
      labs(
        title = "Cell Type Distribution Bubble Plot",
        x = "Group",
        y = "Cell Type",
        size = "Proportion (%)",
        color = "Proportion (%)"
      ) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)
      )
    
    ggsave(file.path(output_dirs$module1_proportion, "03_Proportions_Bubble.pdf"),
           p3, width = 10, height = 8, dpi = 600)
    
    # ========== 图4: 饼图 ==========
    
    cat("生成饼图...\n")
    for (group in unique(prop_data$orig.ident)) {
      group_data <- prop_data %>% filter(orig.ident == group)
      
      p4 <- ggplot(group_data, aes(x = "", y = proportion, fill = cell_type)) +
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = colorRampPalette(c("#E8F4F8", "#A8D8E8", "#7BA7BC", 
                                                       "#5B8FA4", "#D67280"))(
                                                         nrow(group_data))) +
        theme_void() +
        labs(title = sprintf("Cell Type Distribution - %s", group),
             fill = "Cell Type") +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          legend.position = "right"
        )
      
      ggsave(file.path(output_dirs$module1_proportion, 
                       sprintf("04_Pie_Chart_%s.pdf", group)),
             p4, width = 8, height = 6, dpi = 600)
    }
    
    cat("✓ 细胞比例分析完成\n")
    cat("✓ 所有比例图表已保存\n")
    
  }, error = function(e) {
    cat(sprintf("⚠ 比例分析失败: %s\n", e$message))
  })
}

# ============================================================================
# 函数9: 保存结果 - 修复GetAssayData错误
# ============================================================================

save_module1_results <- function(pbmc, output_dirs, work_dir) {
  tryCatch({
    cat("\n=== 步骤1.9: 保存模块1结果 ===\n")
    
    # 保存RData
    cellAnn <- as.character(Idents(pbmc))
    save(pbmc, cellAnn, file = file.path(work_dir, "RData_Storage", "Seurat.RData"))
    cat("✓ Seurat.RData已保存\n")
    
    # 保存RDS
    saveRDS(pbmc, file = file.path(work_dir, "RData_Storage", "Seurat_Object.rds"))
    cat("✓ Seurat_Object.rds已保存\n")
    
    # 保存表达矩阵 - 修复版，使用layer参数
    tryCatch({
      norm_data <- GetAssayData(pbmc, assay = "RNA", layer = "data")
      
      # 转换并限制大小
      if (ncol(norm_data) > 5000) {
        set.seed(123)
        sample_cells <- sample(colnames(norm_data), 5000)
        norm_data <- norm_data[, sample_cells]
      }
      
      write.csv(as.matrix(norm_data),
                file.path(work_dir, "RData_Storage", "Normalized_Expression_Matrix.csv.gz"),
                quote = FALSE)
      cat("✓ 表达矩阵已保存\n")
    }, error = function(e) {
      cat(sprintf("⚠ 表达矩阵保存失败: %s\n", e$message))
    })
    
    # 保存metadata
    write.csv(pbmc@meta.data,
              file.path(work_dir, "Statistical_Reports", "Cell_Metadata.csv"),
              row.names = TRUE)
    cat("✓ Metadata已保存\n")
    
    # 生成最终统计报告
    summary_stats <- data.frame(
      Metric = c(
        "Total Cells",
        "Total Genes",
        "Number of Clusters",
        "Mean UMI per Cell",
        "Mean Genes per Cell"
      ),
      Value = c(
        ncol(pbmc),
        nrow(pbmc),
        length(unique(pbmc$seurat_clusters)),
        round(mean(pbmc$nCount_RNA), 2),
        round(mean(pbmc$nFeature_RNA), 2)
      )
    )
    
    write.csv(summary_stats,
              file.path(work_dir, "Statistical_Reports", "Final_Analysis_Summary.csv"),
              row.names = FALSE)
    cat("✓ 最终统计报告已保存\n")
    
    cat("✓ 模块1结果已全部保存\n")
    
    return(pbmc)
    
  }, error = function(e) {
    cat(sprintf("⚠ 保存结果时出错: %s\n", e$message))
    return(pbmc)
  })
}

cat("✓ Module 1 差异分析和比例分析模块已加载\n")
cat("  - ✓ 增强版火山图\n")
cat("  - ✓ MA plot\n")
cat("  - ✓ 完整比例分析图表\n")
cat("  - ✓ 修复所有GetAssayData错误\n")
