# ============================================================================
# 高级可视化模块：jjVolcano + Gene Cluster Enrichment Heatmap
# 整合您提供的完整代码功能
# ============================================================================

library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

# ============================================================================
# 函数1: jjVolcano火山图生成器
# ============================================================================

generate_jjVolcano_plot <- function(diff_data_dir, output_dir) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("生成jjVolcano火山图\n")
  cat(strrep("=", 80), "\n\n")
  
  # 安装scRNAtoolVis包
  if (!require("scRNAtoolVis", quietly = TRUE)) {
    cat("正在安装scRNAtoolVis包...\n")
    tryCatch({
      install.packages("scRNAtoolVis", repos = "http://cran.us.r-project.org")
      library(scRNAtoolVis)
      cat("✓ scRNAtoolVis包安装成功\n")
    }, error = function(e) {
      cat(sprintf("✗ scRNAtoolVis安装失败: %s\n", e$message))
      return(NULL)
    })
  }
  
  # 查找差异基因文件
  diff_files <- c(
    list.files(path = diff_data_dir, pattern = "_diffGenes\\.txt$", full.names = TRUE, recursive = TRUE),
    list.files(path = diff_data_dir, pattern = "\\.diffGene\\.txt$", full.names = TRUE, recursive = TRUE)
  )
  diff_files <- unique(diff_files)
  
  if (length(diff_files) == 0) {
    cat("⚠ 未找到差异基因文件，跳过jjVolcano图\n")
    return(NULL)
  }
  
  cat(sprintf("✓ 找到 %d 个差异基因文件\n", length(diff_files)))
  
  # 读取并合并数据
  volcano_data_list <- list()
  
  for (file_path in diff_files) {
    file_name <- basename(file_path)
    
    cell_type <- file_name
    cell_type <- str_replace(cell_type, "_diffGenes\\.txt$", "")
    cell_type <- str_replace(cell_type, "\\.diffGene\\.txt$", "")
    cell_type <- str_replace(cell_type, "^\\d+\\.", "")
    cell_type <- str_replace_all(cell_type, "_", " ")
    
    cat(sprintf("  读取: %s\n", file_name))
    
    df <- tryCatch({
      read.table(file_path, header = TRUE, sep = "\t", 
                 stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(df)) next
    
    required_cols <- c("Gene", "avg_log2FC", "p_val_adj")
    if (!all(required_cols %in% colnames(df))) {
      next
    }
    
    df$gene <- df$Gene
    df$cluster <- cell_type
    
    volcano_data_list[[cell_type]] <- df
  }
  
  if (length(volcano_data_list) == 0) {
    cat("⚠ 没有成功读取任何数据\n")
    return(NULL)
  }
  
  volcano_data <- bind_rows(volcano_data_list)
  
  # 过滤数据
  volcano_data <- volcano_data %>%
    filter(
      !str_detect(gene, "\\.[0-9]$"),
      !str_detect(gene, "^LINC")
    )
  
  cat(sprintf("✓ 数据准备完成: %d 个基因\n", nrow(volcano_data)))
  
  # 准备配色
  unique_clusters <- unique(volcano_data$cluster)
  nColors <- length(unique_clusters)
  myColors <- colorRampPalette(brewer.pal(min(12, nColors), "Set3"))(nColors)
  names(myColors) <- unique_clusters
  
  volcano_data$cluster <- factor(volcano_data$cluster, levels = unique_clusters)
  
  # 绘制jjVolcano火山图
  cat("绘制jjVolcano火山图...\n")
  
  tryCatch({
    p_volcano <- jjVolcano(
      diffData = volcano_data,
      log2FC.cutoff = 1,
      pvalue.cutoff = 0.05,
      tile.col = myColors,
      size = 4.5,
      topGeneN = 5,
      legend.position = "top"
    ) +
      labs(
        title = "Volcano Plot of Differentially Expressed Genes",
        subtitle = "Treatment vs Control",
        x = "Cell Type",
        y = "-log10(p-value)"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", 
                                  margin = ggplot2::margin(b = 10)),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray50", 
                                     margin = ggplot2::margin(b = 15)),
        legend.position = "top",
        legend.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.x = element_text(size = 14, face = "bold", angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 13, face = "bold")
      )
    
    # 保存
    volcano_pdf <- file.path(output_dir, "jjVolcano_Plot.pdf")
    volcano_png <- file.path(output_dir, "jjVolcano_Plot.png")
    
    ggsave(volcano_pdf, p_volcano, width = 14, height = 9, dpi = 600)
    ggsave(volcano_png, p_volcano, width = 14, height = 9, dpi = 300)
    
    cat(sprintf("✓ jjVolcano火山图已保存: %s\n", basename(volcano_pdf)))
    
    # 保存数据
    volcano_file <- file.path(output_dir, "jjVolcano_data.csv")
    write.csv(volcano_data, volcano_file, row.names = FALSE)
    
    # 统计摘要
    summary_stats <- volcano_data %>%
      group_by(cluster) %>%
      summarise(
        Total_Genes = n(),
        Significant = sum(abs(avg_log2FC) > 1 & p_val_adj < 0.05),
        Upregulated = sum(avg_log2FC > 1 & p_val_adj < 0.05),
        Downregulated = sum(avg_log2FC < -1 & p_val_adj < 0.05),
        .groups = 'drop'
      )
    
    summary_file <- file.path(output_dir, "jjVolcano_summary.txt")
    write.table(summary_stats, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("\n火山图统计摘要:\n")
    print(summary_stats)
    
    return(volcano_data)
    
  }, error = function(e) {
    cat(sprintf("✗ 绘制火山图失败: %s\n", e$message))
    return(NULL)
  })
}

# ============================================================================
# 函数2: Gene Cluster Enrichment Heatmap生成器
# ============================================================================

generate_gene_cluster_heatmap <- function(diff_data_dir, output_dir) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("生成Gene Cluster Enrichment Heatmap\n")
  cat(strrep("=", 80), "\n\n")
  
  # 安装必要的包
  required_packages <- list(
    list(name = "ClusterGVis", method = "github", repo = "junjunlab/ClusterGVis"),
    list(name = "ggsci", method = "cran"),
    list(name = "ComplexHeatmap", method = "bioc"),
    list(name = "org.Hs.eg.db", method = "bioc"),
    list(name = "clusterProfiler", method = "bioc")
  )
  
  for (pkg_info in required_packages) {
    pkg_name <- pkg_info$name
    
    if (!require(pkg_name, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("正在安装%s...\n", pkg_name))
      
      tryCatch({
        if (pkg_info$method == "cran") {
          install.packages(pkg_name)
        } else if (pkg_info$method == "bioc") {
          if (!require("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(pkg_name)
        } else if (pkg_info$method == "github") {
          if (!require("remotes", quietly = TRUE)) {
            install.packages("remotes")
          }
          remotes::install_github(pkg_info$repo)
        }
        library(pkg_name, character.only = TRUE)
        cat(sprintf("✓ %s安装成功\n", pkg_name))
      }, error = function(e) {
        cat(sprintf("⚠ %s安装失败\n", pkg_name))
      })
    }
  }
  
  # 检查ClusterGVis是否可用
  if (!require("ClusterGVis", quietly = TRUE)) {
    cat("⚠ ClusterGVis不可用，跳过基因聚类热图\n")
    return(NULL)
  }
  
  library(ClusterGVis)
  library(ggsci)
  library(ComplexHeatmap)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  
  # 查找差异基因文件
  diff_files <- c(
    list.files(path = diff_data_dir, pattern = "_diffGenes\\.txt$", full.names = TRUE, recursive = TRUE),
    list.files(path = diff_data_dir, pattern = "\\.diffGene\\.txt$", full.names = TRUE, recursive = TRUE)
  )
  diff_files <- unique(diff_files)
  
  if (length(diff_files) == 0) {
    cat("⚠ 未找到差异基因文件\n")
    return(NULL)
  }
  
  cat("准备基因聚类数据...\n")
  
  # 读取所有数据
  all_diff_data <- list()
  
  for (file_path in diff_files) {
    file_name <- basename(file_path)
    
    cell_type <- file_name
    cell_type <- str_replace(cell_type, "_diffGenes\\.txt$", "")
    cell_type <- str_replace(cell_type, "\\.diffGene\\.txt$", "")
    cell_type <- str_replace(cell_type, "^\\d+\\.", "")
    cell_type <- str_replace_all(cell_type, "_", " ")
    
    df <- tryCatch({
      read.table(file_path, header = TRUE, sep = "\t", 
                 stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(df) && all(c("Gene", "avg_log2FC", "p_val_adj") %in% colnames(df))) {
      df$CellType <- cell_type
      all_diff_data[[cell_type]] <- df
    }
  }
  
  combined_diff_data <- bind_rows(all_diff_data)
  
  # 过滤
  combined_diff_data <- combined_diff_data %>%
    filter(
      !str_detect(Gene, "\\.[0-9]$"),
      !str_detect(Gene, "^LINC")
    )
  
  cat(sprintf("✓ 总共 %d 个基因用于聚类分析\n", nrow(combined_diff_data)))
  
  # 选择显著差异基因
  cat("筛选显著差异基因...\n")
  
  sig_genes <- combined_diff_data %>%
    filter(abs(avg_log2FC) > 1, p_val_adj < 0.05) %>%
    group_by(CellType) %>%
    arrange(desc(abs(avg_log2FC)), p_val_adj) %>%
    slice_head(n = 50) %>%
    ungroup()
  
  all_sig_genes <- unique(sig_genes$Gene)
  
  # 限制基因数量
  if (length(all_sig_genes) > 500) {
    top_genes_per_celltype <- sig_genes %>%
      group_by(CellType) %>%
      arrange(p_val_adj) %>%
      slice_head(n = 30) %>%
      pull(Gene) %>%
      unique()
    
    all_sig_genes <- head(top_genes_per_celltype, 500)
    cat(sprintf("基因数量限制为 %d 个\n", length(all_sig_genes)))
  }
  
  cat(sprintf("✓ 选择 %d 个显著差异基因进行聚类\n", length(all_sig_genes)))
  
  # 创建表达矩阵
  cat("创建基因表达矩阵...\n")
  
  unique_celltypes <- unique(combined_diff_data$CellType)
  
  expr_matrix <- matrix(0, nrow = length(all_sig_genes), ncol = length(unique_celltypes))
  rownames(expr_matrix) <- all_sig_genes
  colnames(expr_matrix) <- unique_celltypes
  
  for (i in 1:nrow(sig_genes)) {
    gene <- sig_genes$Gene[i]
    celltype <- sig_genes$CellType[i]
    
    if (gene %in% all_sig_genes && celltype %in% unique_celltypes) {
      expr_matrix[gene, celltype] <- sig_genes$avg_log2FC[i]
    }
  }
  
  # 添加少量噪音以避免聚类问题
  set.seed(123)
  expr_matrix <- expr_matrix + matrix(rnorm(length(expr_matrix), 0, 0.1), 
                                      nrow = nrow(expr_matrix))
  
  cat(sprintf("✓ 表达矩阵: %d 基因 × %d 细胞类型\n", 
              nrow(expr_matrix), ncol(expr_matrix)))
  
  # 基因聚类
  cat("进行基因聚类分析...\n")
  
  expr_df <- as.data.frame(expr_matrix)
  
  cluster_num <- min(8, max(3, round(nrow(expr_df) / 50)))
  cat(sprintf("聚类数量: %d\n", cluster_num))
  
  gene_clusters <- tryCatch({
    clusterData(
      obj = expr_df,
      cluster.method = "kmeans",
      cluster.num = cluster_num
    )
  }, error = function(e) {
    cat("kmeans失败，尝试hclust...\n")
    clusterData(
      obj = expr_df,
      cluster.method = "hclust",
      cluster.num = cluster_num
    )
  })
  
  cat("✓ 基因聚类完成\n")
  
  # GO富集分析
  cat("进行GO富集分析...\n")
  
  gene_enrich <- NULL
  
  tryCatch({
    gene_enrich <- enrichCluster(
      object = gene_clusters,
      OrgDb = org.Hs.eg.db,
      type = "BP",
      organism = "hsa",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      topn = 5,
      seed = 5201314
    )
    
    if (!is.null(gene_enrich) && nrow(gene_enrich) > 0) {
      cat(sprintf("✓ GO富集分析完成: %d 个富集项\n", nrow(gene_enrich)))
    } else {
      gene_enrich <- NULL
    }
  }, error = function(e) {
    cat(sprintf("⚠ GO富集分析失败: %s\n", e$message))
    gene_enrich <- NULL
  })
  
  # 选择标记基因
  cat("选择标记基因...\n")
  
  known_markers <- c(
    "CD3D", "CD3E", "CD4", "CD8A", "CD19", "CD14",
    "LYZ", "MS4A1", "GNLY", "NKG7", "FCGR3A", "CST3",
    "IL7R", "CCR7", "GZMB", "PRF1", "FCER1A", "CD1C"
  )
  
  mark_genes <- intersect(known_markers, rownames(expr_matrix))
  
  if (length(mark_genes) < 10) {
    top_var_genes <- sig_genes %>%
      arrange(p_val_adj) %>%
      filter(!Gene %in% mark_genes) %>%
      pull(Gene) %>%
      unique() %>%
      head(20)
    
    mark_genes <- c(mark_genes, top_var_genes)
  }
  
  mark_genes <- intersect(mark_genes, rownames(expr_matrix))
  mark_genes <- head(mark_genes, 20)
  
  cat(sprintf("✓ 选择 %d 个标记基因\n", length(mark_genes)))
  
  # 绘制完整热图
  cat("绘制基因聚类富集热图...\n")
  
  heatmap_complete_pdf <- file.path(output_dir, "Gene_Cluster_Enrichment_Heatmap.pdf")
  
  tryCatch({
    pdf(heatmap_complete_pdf, width = 16, height = 12, onefile = FALSE)
    
    if (!is.null(gene_enrich) && nrow(gene_enrich) > 0) {
      num_terms <- nrow(gene_enrich)
      go_colors <- rep(pal_d3()(8), length.out = num_terms)
      
      visCluster(
        object = gene_clusters,
        plot.type = "both",
        column_names_rot = 45,
        show_row_dend = FALSE,
        markGenes = mark_genes,
        markGenes.side = "left",
        annoTerm.data = gene_enrich,
        line.side = "left",
        go.col = go_colors,
        go.size = "pval",
        add.bar = TRUE
      )
    } else {
      visCluster(
        object = gene_clusters,
        plot.type = "heatmap",
        column_names_rot = 45,
        show_row_dend = FALSE,
        markGenes = mark_genes,
        markGenes.side = "left"
      )
    }
    
    dev.off()
    cat(sprintf("✓ 完整热图已保存: %s\n", basename(heatmap_complete_pdf)))
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    cat(sprintf("✗ 绘制完整热图失败: %s\n", e$message))
  })
  
  # 绘制简版热图
  cat("绘制简版热图...\n")
  
  heatmap_simple_pdf <- file.path(output_dir, "Gene_Cluster_Heatmap_Simple.pdf")
  
  tryCatch({
    pdf(heatmap_simple_pdf, width = 12, height = 10, onefile = FALSE)
    
    visCluster(
      object = gene_clusters,
      plot.type = "heatmap",
      column_names_rot = 45,
      show_row_dend = FALSE,
      markGenes = mark_genes,
      markGenes.side = "left"
    )
    
    dev.off()
    cat(sprintf("✓ 简版热图已保存: %s\n", basename(heatmap_simple_pdf)))
  }, error = function(e) {
    if (dev.cur() != 1) dev.off()
    cat(sprintf("✗ 绘制简版热图失败: %s\n", e$message))
  })
  
  # 保存结果
  cat("保存分析结果...\n")
  
  cluster_results_file <- file.path(output_dir, "Gene_Cluster_Results.csv")
  write.csv(gene_clusters$wide.res, cluster_results_file, row.names = TRUE)
  
  expr_file <- file.path(output_dir, "Expression_Matrix.csv")
  write.csv(expr_matrix, expr_file, row.names = TRUE)
  
  if (!is.null(gene_enrich) && nrow(gene_enrich) > 0) {
    enrich_file <- file.path(output_dir, "GO_Enrichment_Results.csv")
    write.csv(gene_enrich, enrich_file, row.names = FALSE)
  }
  
  cat("✓ 基因聚类热图分析完成\n")
  
  return(list(
    gene_clusters = gene_clusters,
    expr_matrix = expr_matrix,
    gene_enrich = gene_enrich
  ))
}

# ============================================================================
# 主函数：运行完整的高级可视化
# ============================================================================

run_advanced_visualizations <- function(work_dir, data_dir = NULL) {
  
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("高级可视化模块\n")
  cat("jjVolcano + Gene Cluster Enrichment Heatmap\n")
  cat(strrep("=", 80), "\n\n")
  
  # 设置目录
  if (is.null(data_dir)) {
    data_dir <- work_dir
  }
  
  output_dir <- file.path(work_dir, "Advanced_Visualizations")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("工作目录: %s\n", work_dir))
  cat(sprintf("数据目录: %s\n", data_dir))
  cat(sprintf("输出目录: %s\n\n", output_dir))
  
  # 运行jjVolcano
  volcano_result <- generate_jjVolcano_plot(data_dir, output_dir)
  
  # 运行Gene Cluster Heatmap
  heatmap_result <- generate_gene_cluster_heatmap(data_dir, output_dir)
  
  # 最终报告
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("高级可视化完成！\n")
  cat(strrep("=", 80), "\n")
  cat(sprintf("输出目录: %s\n", output_dir))
  cat("\n主要输出文件:\n")
  cat("  1. jjVolcano_Plot.pdf - jjVolcano火山图\n")
  cat("  2. jjVolcano_Plot.png - jjVolcano火山图 (PNG)\n")
  cat("  3. Gene_Cluster_Enrichment_Heatmap.pdf - 完整聚类富集热图\n")
  cat("  4. Gene_Cluster_Heatmap_Simple.pdf - 简版聚类热图\n")
  cat("  5. GO_Enrichment_Results.csv - GO富集结果\n")
  cat("  6. Gene_Cluster_Results.csv - 聚类结果\n")
  cat(strrep("=", 80), "\n")
  
  return(list(
    volcano = volcano_result,
    heatmap = heatmap_result
  ))
}

cat("✓ 高级可视化模块已加载\n")
cat("  - ✓ jjVolcano火山图\n")
cat("  - ✓ Gene Cluster Enrichment Heatmap\n")
cat("  - ✓ 完整功能整合\n")
