# ==============================================================================
# Module 2: 高级分析 - 完整终极版 v5.0
# 包含: CellChat, 虚拟敲除, GO富集, 基因可视化, 高级图表
# 修复: 所有GetAssayData错误 (slot → layer)
# ==============================================================================

# ==============================================================================
# Module 2.1: CellChat 细胞通讯分析
# ==============================================================================

run_cellchat_analysis <- function(seurat_obj, params, output_dirs) {
  cat("\n=== 模块2.1: CellChat 细胞通讯分析 ===\n")
  
  tryCatch({
    # 检查CellChat包
    if (!requireNamespace("CellChat", quietly = TRUE)) {
      cat("⚠ CellChat包未安装，尝试安装...\n")
      tryCatch({
        remotes::install_github("sqjin/CellChat")
      }, error = function(e) {
        cat("✗ CellChat安装失败，跳过此分析\n")
        return(NULL)
      })
    }
    
    library(CellChat)
    
    # 准备数据 - 修复GetAssayData
    cat("准备CellChat数据...\n")
    expMatrix <- as.matrix(GetAssayData(seurat_obj, layer = "data"))  # ✓ 修复
    
    meta <- data.frame(
      labels = as.character(Idents(seurat_obj)),
      row.names = colnames(seurat_obj)
    )
    
    # 创建CellChat对象
    cat("创建CellChat对象...\n")
    cellchat <- createCellChat(object = expMatrix, meta = meta, group.by = "labels")
    cellchat <- setIdent(cellchat, ident.use = "labels")
    groupSize <- as.numeric(table(cellchat@idents))
    
    # 加载配体-受体数据库
    cat(sprintf("加载数据库: %s\n", params$cellchat_db))
    CellChatDB <- CellChatDB.human
    
    if (params$cellchat_db == "Secreted Signaling") {
      CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
    } else {
      CellChatDB.use <- CellChatDB
    }
    
    cellchat@DB <- CellChatDB.use
    
    # 查看数据库分类
    pdf(file = file.path(output_dirs$module2_cellchat, "Database_Category.pdf"),
        width = 7, height = 5)
    showDatabaseCategory(CellChatDB)
    dev.off()
    
    # 数据预处理
    cat("数据预处理...\n")
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    
    # 计算通讯概率
    cat("计算细胞通讯概率...\n")
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = params$cellchat_min_cells)
    
    # 提取通讯网络
    df.net <- subsetCommunication(cellchat)
    write.table(df.net,
               file = file.path(output_dirs$module2_cellchat, "Communication_Network.txt"),
               sep = "\t", row.names = FALSE, quote = FALSE)
    
    # 计算通路水平的通讯
    cat("计算信号通路水平通讯...\n")
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # === 网络可视化 - 交互数量 ===
    cat("生成网络图（交互数量）...\n")
    pdf(file = file.path(output_dirs$module2_cellchat, "Network_Count.pdf"),
        width = 7, height = 6)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                    weight.scale = TRUE, label.edge = FALSE,
                    title.name = "Number of interactions")
    dev.off()
    
    # === 网络可视化 - 交互强度 ===
    cat("生成网络图（交互强度）...\n")
    pdf(file = file.path(output_dirs$module2_cellchat, "Network_Weight.pdf"),
        width = 7, height = 6)
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                    weight.scale = TRUE, label.edge = FALSE,
                    title.name = "Interaction strength")
    dev.off()
    
    # === 单个细胞类型的网络 ===
    cat("生成单细胞类型网络图...\n")
    weight_mat <- cellchat@net$weight
    cell_types <- unique(cellchat@idents)
    
    pdf(file = file.path(output_dirs$module2_cellchat, "Individual_Cell_Networks.pdf"),
        width = 12, height = 10)
    par(mfrow = c(2, 3), xpd = TRUE)
    for (i in 1:min(6, length(cell_types))) {
      mat2 <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat),
                    dimnames = dimnames(weight_mat))
      mat2[cell_types[i], ] <- weight_mat[cell_types[i], ]
      netVisual_circle(mat2, vertex.weight = groupSize,
                      weight.scale = TRUE, edge.weight.max = max(weight_mat),
                      title.name = paste0(cell_types[i], " signaling"))
    }
    dev.off()
    
    # === 信号通路分析 ===
    cat("分析信号通路...\n")
    pathways.show <- cellchat@netP$pathways[1:min(5, length(cellchat@netP$pathways))]
    
    # 层级图
    for (pathway in pathways.show) {
      pdf(file = file.path(output_dirs$module2_cellchat,
                          sprintf("Pathway_%s_Hierarchy.pdf", pathway)),
          width = 10, height = 8)
      tryCatch({
        netVisual_aggregate(cellchat, signaling = pathway, layout = "hierarchy")
      }, error = function(e) {
        cat(sprintf("⚠ 通路 %s 层级图生成失败\n", pathway))
      })
      dev.off()
    }
    
    # 和弦图
    for (pathway in pathways.show) {
      pdf(file = file.path(output_dirs$module2_cellchat,
                          sprintf("Pathway_%s_Chord.pdf", pathway)),
          width = 8, height = 8)
      tryCatch({
        netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
      }, error = function(e) {
        cat(sprintf("⚠ 通路 %s 和弦图生成失败\n", pathway))
      })
      dev.off()
    }
    
    # === 气泡图 - 配体-受体对 ===
    cat("生成配体-受体对气泡图...\n")
    pdf(file = file.path(output_dirs$module2_cellchat, "LR_Pairs_BubblePlot.pdf"),
        width = 12, height = 8)
    netVisual_bubble(cellchat, sources.use = 1:length(groupSize),
                    targets.use = 1:length(groupSize),
                    remove.isolate = FALSE)
    dev.off()
    
    # 保存CellChat对象
    saveRDS(cellchat,
           file = file.path(output_dirs$module2_cellchat, "CellChat_Object.rds"))
    
    cat("✓ CellChat分析完成\n")
    cat(sprintf("  - 通讯网络文件: Communication_Network.txt\n"))
    cat(sprintf("  - 网络图: 2个PDF\n"))
    cat(sprintf("  - 通路图: %d 个通路\n", length(pathways.show)))
    
    return(cellchat)
    
  }, error = function(e) {
    cat(sprintf("✗ CellChat分析失败: %s\n", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.2: 虚拟基因敲除 (Virtual Gene Knockout)
# ==============================================================================

run_virtual_knockout <- function(seurat_obj, params, output_dirs) {
  cat("\n=== 模块2.2: 虚拟基因敲除分析 ===\n")
  
  tryCatch({
    # 检查scTenifoldKnk包
    if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
      cat("⚠ scTenifoldKnk包未安装，尝试安装...\n")
      tryCatch({
        remotes::install_github("cailab-tamu/scTenifoldKnk")
      }, error = function(e) {
        cat("✗ scTenifoldKnk安装失败，跳过此分析\n")
        return(NULL)
      })
    }
    
    library(scTenifoldKnk)
    
    target_gene <- params$knockout_target_gene
    cat(sprintf("敲除目标基因: %s\n", target_gene))
    
    # 检查基因是否存在
    if (!(target_gene %in% rownames(seurat_obj))) {
      cat(sprintf("✗ 基因 %s 不在数据集中\n", target_gene))
      return(NULL)
    }
    
    # 子集数据（如果细胞太多）
    if (ncol(seurat_obj) > params$knockout_max_cells) {
      cat(sprintf("细胞数 %d 超过限制 %d，进行采样...\n",
                  ncol(seurat_obj), params$knockout_max_cells))
      set.seed(123)
      pbmc_subset <- subset(seurat_obj,
                           cells = sample(colnames(seurat_obj),
                                        params$knockout_max_cells))
    } else {
      pbmc_subset <- seurat_obj
    }
    
    # 提取表达矩阵 - 修复GetAssayData
    cat("提取表达矩阵...\n")
    countMat <- GetAssayData(pbmc_subset, layer = "counts")  # ✓ 修复
    
    # 识别高变基因
    cat("识别高变基因...\n")
    pbmc_subset <- FindVariableFeatures(object = pbmc_subset,
                                       selection.method = "vst",
                                       nfeatures = 10000)
    hvgs <- VariableFeatures(pbmc_subset)
    
    # 准备数据
    data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])
    
    cat(sprintf("数据维度: %d 基因 × %d 细胞\n", nrow(data), ncol(data)))
    
    # 执行虚拟敲除
    cat("执行虚拟敲除分析（这可能需要较长时间）...\n")
    result <- scTenifoldKnk(
      countMatrix = data,
      gKO = target_gene,
      qc = TRUE,
      qc_mtThreshold = params$knockout_qc_mtThreshold,
      qc_minLSize = params$knockout_qc_minLSize,
      nc_nNet = params$knockout_nc_nNet,
      nc_nCells = params$knockout_nc_nCells
    )
    
    # 提取差异结果
    df <- result$diffRegulation
    df <- df[df$gene != target_gene, ]
    
    # 保存显著差异基因
    sig_diff <- df[df$p.adj < 0.05, ]
    write.table(sig_diff,
               file = file.path(output_dirs$module2_knockout, "Significant_Diff_Genes.txt"),
               sep = "\t", quote = FALSE, row.names = FALSE)
    
    write.table(df,
               file = file.path(output_dirs$module2_knockout, "All_Diff_Genes.txt"),
               sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat(sprintf("✓ 发现 %d 个显著差异基因\n", nrow(sig_diff)))
    
    # === 柱状图 ===
    cat("生成柱状图...\n")
    top_genes <- head(df[order(-df$FC), ], 20)
    
    p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
      geom_bar(stat = 'identity', fill = '#7BA7BC') +
      coord_flip() +
      labs(title = sprintf("Top 20 Genes After %s Knockout", target_gene),
           x = "Gene", y = "Fold Change") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10)
      )
    
    ggsave(p1,
          filename = file.path(output_dirs$module2_knockout, "Barplot.pdf"),
          width = 8, height = 6, dpi = 600)
    
    # === 火山图 ===
    cat("生成火山图...\n")
    df$log_p.adj <- -log10(df$p.adj)
    df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
    label_genes <- subset(df, p.adj < 0.05)
    
    y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)
    
    p2 <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_manual(values = c("Significant" = "#D67280",
                                    "Not significant" = "#E0E0E0")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                color = "#D67280", size = 1) +
      labs(title = sprintf("Volcano Plot: %s Knockout", target_gene),
           x = "Z-score", y = "-log10(p.adj)") +
      theme_minimal() +
      coord_cartesian(ylim = c(0, y_upper)) +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12)
      )
    
    # 添加基因标签
    if (nrow(label_genes) > 0 & requireNamespace("ggrepel", quietly = TRUE)) {
      library(ggrepel)
      p2 <- p2 + geom_text_repel(
        data = head(label_genes, 20),
        aes(label = gene),
        size = 3,
        max.overlaps = 50,
        box.padding = 0.5
      )
    }
    
    ggsave(p2,
          filename = file.path(output_dirs$module2_knockout, "Volcano_Plot.pdf"),
          width = 10, height = 8, dpi = 600)
    
    # === 饼图 ===
    cat("生成饼图...\n")
    sig_count <- table(df$significant)
    sig_df <- as.data.frame(sig_count)
    colnames(sig_df) <- c("category", "count")
    sig_df$percentage <- paste0(round(sig_df$count / sum(sig_df$count) * 100, 1), "%")
    
    p3 <- ggplot(sig_df, aes(x = "", y = count, fill = category)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(label = percentage),
               position = position_stack(vjust = 0.5), size = 5, fontface = "bold") +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = c("Significant" = "#D67280",
                                   "Not significant" = "#E8E8E8")) +
      labs(title = sprintf("Proportion of Significant Genes (%s KO)", target_gene),
           fill = "") +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )
    
    ggsave(p3,
          filename = file.path(output_dirs$module2_knockout, "Pie_Chart.pdf"),
          width = 8, height = 6, dpi = 600)
    
    # 保存结果对象
    saveRDS(result,
           file = file.path(output_dirs$module2_knockout, "scTenifoldKnk_Result.rds"))
    
    cat("✓ 虚拟敲除分析完成\n")
    cat(sprintf("  - 显著差异基因: %d\n", nrow(sig_diff)))
    cat(sprintf("  - 柱状图: Barplot.pdf\n"))
    cat(sprintf("  - 火山图: Volcano_Plot.pdf\n"))
    cat(sprintf("  - 饼图: Pie_Chart.pdf\n"))
    
    return(result)
    
  }, error = function(e) {
    cat(sprintf("✗ 虚拟敲除分析失败: %s\n", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.3: GO富集分析
# ==============================================================================

run_go_enrichment <- function(seurat_obj, params, output_dirs, gene_file = NULL) {
  cat("\n=== 模块2.3: GO富集分析 ===\n")
  
  tryCatch({
    # 检查必要的包
    required_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2")
    for (pkg in required_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        cat(sprintf("⚠ 缺少包 %s，跳过GO富集分析\n", pkg))
        return(NULL)
      }
    }
    
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    
    # 读取基因列表
    if (!is.null(gene_file) && file.exists(gene_file)) {
      cat(sprintf("从文件读取基因: %s\n", gene_file))
      rt <- read.table(gene_file, header = TRUE, sep = "\t", check.names = FALSE)
      genes <- unique(as.vector(rt[, 1]))
    } else if (!is.null(seurat_obj@misc$sig_diff_genes)) {
      cat("使用差异基因进行GO分析...\n")
      genes <- seurat_obj@misc$sig_diff_genes$gene
    } else {
      cat("⚠ 没有基因列表，跳过GO分析\n")
      return(NULL)
    }
    
    cat(sprintf("基因数量: %d\n", length(genes)))
    
    # 基因ID转换
    cat("转换基因ID...\n")
    entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
    entrezIDs <- as.character(entrezIDs)
    
    # 移除NA
    valid_ids <- entrezIDs[entrezIDs != "NA"]
    cat(sprintf("有效Entrez ID: %d\n", length(valid_ids)))
    
    if (length(valid_ids) == 0) {
      cat("⚠ 没有有效的Entrez ID\n")
      return(NULL)
    }
    
    # GO富集分析
    cat("执行GO富集分析...\n")
    kk <- enrichGO(
      gene = valid_ids,
      OrgDb = org.Hs.eg.db,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      ont = "all",
      readable = TRUE
    )
    
    # 保存结果
    write.table(kk,
               file = file.path(output_dirs$module2_go, "GO_Enrichment_Results.txt"),
               sep = "\t", quote = FALSE, row.names = FALSE)
    
    cat(sprintf("✓ GO富集分析完成，发现 %d 个富集项\n", nrow(kk)))
    
    if (nrow(kk) == 0) {
      cat("⚠ 没有显著富集的GO terms\n")
      return(NULL)
    }
    
    # === 柱状图 ===
    cat("生成柱状图...\n")
    showNum <- min(30, nrow(kk))
    
    pdf(file = file.path(output_dirs$module2_go, "GO_Barplot.pdf"),
        width = 10, height = 8)
    p1 <- barplot(kk, showCategory = showNum,
                 title = "GO Enrichment Analysis")
    print(p1)
    dev.off()
    
    # === 气泡图 ===
    cat("生成气泡图...\n")
    pdf(file = file.path(output_dirs$module2_go, "GO_DotPlot.pdf"),
        width = 10, height = 8)
    p2 <- dotplot(kk, showCategory = showNum,
                 title = "GO Enrichment Dotplot")
    print(p2)
    dev.off()
    
    # === GO关系图 ===
    if (nrow(kk) >= 10) {
      cat("生成GO关系图...\n")
      kk2 <- pairwise_termsim(kk)
      
      pdf(file = file.path(output_dirs$module2_go, "GO_Network.pdf"),
          width = 12, height = 10)
      p3 <- emapplot(kk2, showCategory = min(30, nrow(kk2)))
      print(p3)
      dev.off()
    }
    
    # === 分类气泡图 ===
    cat("生成分类气泡图...\n")
    pdf(file = file.path(output_dirs$module2_go, "GO_Category_DotPlot.pdf"),
        width = 10, height = 8)
    p4 <- dotplot(kk, split = "ONTOLOGY", showCategory = 10) +
      facet_grid(ONTOLOGY~., scales = "free") +
      theme(strip.text = element_text(face = "bold", size = 12))
    print(p4)
    dev.off()
    
    # === Circos图（如果有circlize包）===
    if (requireNamespace("circlize", quietly = TRUE) && nrow(kk) >= 10) {
      cat("生成Circos图...\n")
      library(circlize)
      
      tryCatch({
        # 准备数据
        df <- as.data.frame(kk)
        df <- df[order(-df$Count), ]
        df <- head(df, 20)
        
        # 创建床位数据
        BgGene <- max(df$Count) * 1.1
        
        bed1 <- data.frame(
          chr = paste0("GO:", 1:nrow(df)),
          start = 0,
          end = BgGene,
          stringsAsFactors = FALSE
        )
        
        bed2 <- data.frame(
          chr = paste0("GO:", 1:nrow(df)),
          start = BgGene - df$Count,
          end = BgGene,
          Description = substring(df$Description, 1, 30),
          stringsAsFactors = FALSE
        )
        
        bed3 <- data.frame(
          chr = paste0("GO:", 1:nrow(df)),
          start = BgGene - df$Count - 20,
          end = BgGene - df$Count - 10,
          GO_ID = df$ID,
          stringsAsFactors = FALSE
        )
        
        bed4 <- data.frame(
          chr = paste0("GO:", 1:nrow(df)),
          start = 0,
          end = df$Count,
          value = -log10(df$p.adjust),
          color = ifelse(df$ONTOLOGY == "BP", "#E8B4B8",
                        ifelse(df$ONTOLOGY == "CC", "#88CAE0", "#B8E6B8")),
          stringsAsFactors = FALSE
        )
        
        # 绘制Circos图
        pdf(file = file.path(output_dirs$module2_go, "GO_Circos_Plot.pdf"),
            width = 10, height = 10)
        
        circos.clear()
        circos.par(gap.degree = c(rep(2, nrow(df) - 1), 10),
                  start.degree = 90)
        
        main.col <- rainbow(nrow(bed1), alpha = 0.5)
        
        circos.genomicInitialize(bed1, plotType = NULL)
        
        circos.track(ylim = c(0, 1), track.height = 0.08,
                    bg.col = main.col, bg.border = NA,
                    panel.fun = function(x, y) {
                      xlim <- get.cell.meta.data("xlim")
                      ylim <- get.cell.meta.data("ylim")
                      circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8,
                                 facing = "bending.inside", niceFacing = TRUE)
                    })
        
        for (si in get.all.sector.index()) {
          circos.axis(h = "top", labels.cex = 0.6, sector.index = si,
                     track.index = 1, major.at = seq(0, max(BgGene), by = 100),
                     labels.facing = "clockwise")
        }
        
        circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1,
                           bg.border = "white",
                           panel.fun = function(region, value, ...) {
                             circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                              col = value[, 2], border = NA, ...)
                             circos.genomicText(region, value, y = 0.4,
                                              labels = value[, 1], adj = 0, cex = 0.8, ...)
                           })
        
        circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1,
                           bg.border = "white",
                           panel.fun = function(region, value, ...) {
                             circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                              col = '#BA55D3', border = NA, ...)
                             circos.genomicText(region, value, y = 0.4,
                                              labels = value[, 1], cex = 0.9, adj = 0, ...)
                           })
        
        circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35,
                           bg.border = "white", bg.col = "grey90",
                           panel.fun = function(region, value, ...) {
                             cell.xlim <- get.cell.meta.data("cell.xlim")
                             cell.ylim <- get.cell.meta.data("cell.ylim")
                             for (j in 1:9) {
                               y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                               circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                             }
                             circos.genomicRect(region, value, ytop = 0,
                                              ybottom = value[, 1], col = value[, 2],
                                              border = NA, ...)
                           })
        
        circos.clear()
        dev.off()
        
        cat("✓ Circos图已生成\n")
        
      }, error = function(e) {
        cat(sprintf("⚠ Circos图生成失败: %s\n", e$message))
      })
    }
    
    cat("✓ GO富集分析完成\n")
    cat(sprintf("  - 富集项数: %d\n", nrow(kk)))
    cat(sprintf("  - 柱状图: GO_Barplot.pdf\n"))
    cat(sprintf("  - 气泡图: GO_DotPlot.pdf\n"))
    cat(sprintf("  - 网络图: GO_Network.pdf\n"))
    
    return(kk)
    
  }, error = function(e) {
    cat(sprintf("✗ GO富集分析失败: %s\n", e$message))
    return(NULL)
  })
}

cat("✓ Module 2 (Part 1/2) 已加载：CellChat, 虚拟敲除, GO富集\n")
