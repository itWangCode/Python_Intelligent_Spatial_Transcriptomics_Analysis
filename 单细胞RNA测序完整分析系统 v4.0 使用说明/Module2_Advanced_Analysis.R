#!/usr/bin/env Rscript
# ==============================================================================
# 模块2：高级分析和可视化
# 功能：CellChat细胞通讯、虚拟敲除、GO富集分析、基因可视化
# ==============================================================================

# ==============================================================================
# Module 2.1: CellChat 细胞通讯分析
# ==============================================================================

run_cellchat_analysis <- function(seurat_obj, params, output_dirs) {
  message("\n=== 模块2.1: CellChat 细胞通讯分析 ===")
  
  tryCatch({
    # 检查CellChat包
    if (!requireNamespace("CellChat", quietly = TRUE)) {
      warning("CellChat包未安装，跳过此分析")
      return(NULL)
    }
    
    library(CellChat)
    
    # 准备数据
    message("准备CellChat数据...")
    expMatrix <- as.matrix(GetAssayData(seurat_obj, slot = "data"))
    meta <- data.frame(
      labels = as.character(Idents(seurat_obj)),
      row.names = colnames(seurat_obj)
    )
    
    # 创建CellChat对象
    message("创建CellChat对象...")
    cellchat <- createCellChat(object = expMatrix, meta = meta, group.by = "labels")
    cellchat <- setIdent(cellchat, ident.use = "labels")
    groupSize <- as.numeric(table(cellchat@idents))
    
    # 加载配体-受体数据库
    message(sprintf("加载数据库: %s", params$cellchat_db))
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
    message("数据预处理...")
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, PPI.human)
    
    # 计算通讯概率
    message("计算细胞通讯概率...")
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = params$cellchat_min_cells)
    
    # 提取通讯网络
    df.net <- subsetCommunication(cellchat)
    write.table(df.net,
               file = file.path(output_dirs$module2_cellchat, "Communication_Network.txt"),
               sep = "\t", row.names = FALSE, quote = FALSE)
    
    # 计算通路水平的通讯
    message("计算信号通路水平通讯...")
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    
    # 网络可视化 - 交互数量
    message("生成网络图...")
    pdf(file = file.path(output_dirs$module2_cellchat, "Network_Count.pdf"),
        width = 7, height = 6)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                    weight.scale = TRUE, label.edge = FALSE,
                    title.name = "Number of interactions")
    dev.off()
    
    # 网络可视化 - 交互强度
    pdf(file = file.path(output_dirs$module2_cellchat, "Network_Weight.pdf"),
        width = 7, height = 6)
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                    weight.scale = TRUE, label.edge = FALSE,
                    title.name = "Interaction strength")
    dev.off()
    
    # 单个细胞类型的网络
    message("生成单细胞类型网络图...")
    weight_mat <- cellchat@net$weight
    cell_types <- unique(cellchat@idents)
    
    pdf(file = file.path(output_dirs$module2_cellchat, "Individual_Cell_Networks.pdf"),
        width = 12, height = 10)
    par(mfrow = c(ceiling(length(cell_types)/3), 3), xpd = TRUE)
    for (cel in cell_types) {
      cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat),
                       dimnames = dimnames(weight_mat))
      cir_mat[cel, ] <- weight_mat[cel, ]
      netVisual_circle(cir_mat, vertex.weight = groupSize, weight.scale = TRUE,
                      edge.weight.max = max(weight_mat), vertex.label.cex = 0.8,
                      title.name = cel)
    }
    dev.off()
    
    # 为每个细胞类型单独保存
    for (cel in cell_types) {
      cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat),
                       dimnames = dimnames(weight_mat))
      cir_mat[cel, ] <- weight_mat[cel, ]
      
      pdf(file = file.path(output_dirs$module2_cellchat,
                          paste0("Network_", gsub(" ", "_", cel), ".pdf")),
          width = 6.5, height = 5.5)
      netVisual_circle(cir_mat, vertex.weight = groupSize, weight.scale = TRUE,
                      edge.weight.max = max(weight_mat), vertex.label.cex = 0.8,
                      title.name = cel)
      dev.off()
    }
    
    # 气泡图
    message("生成气泡图...")
    pdf(file = file.path(output_dirs$module2_cellchat, "Bubble_Plot.pdf"),
        width = 12, height = 10)
    netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45, font.size = 6)
    dev.off()
    
    # 和弦图（通路级别）
    if (nrow(df.net) > 0) {
      pathways <- unique(df.net$pathway_name)[1:min(10, length(unique(df.net$pathway_name)))]
      
      for (pathway in pathways) {
        tryCatch({
          pdf(file = file.path(output_dirs$module2_cellchat,
                              paste0("Pathway_", gsub(" ", "_", pathway), ".pdf")),
              width = 8, height = 8)
          netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
          dev.off()
        }, error = function(e) {
          message(sprintf("跳过通路 %s: %s", pathway, e$message))
        })
      }
    }
    
    # 保存CellChat对象
    saveRDS(cellchat,
           file = file.path(output_dirs$module2_cellchat, "CellChat_Object.rds"))
    
    message("✓ CellChat分析完成")
    
    return(cellchat)
    
  }, error = function(e) {
    warning(sprintf("CellChat分析失败: %s", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.2: scTenifoldKnk 虚拟基因敲除分析
# ==============================================================================

run_virtual_knockout <- function(seurat_obj, params, output_dirs) {
  message("\n=== 模块2.2: 虚拟基因敲除分析 ===")
  
  tryCatch({
    # 检查包
    if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
      warning("scTenifoldKnk包未安装，跳过此分析")
      return(NULL)
    }
    
    library(scTenifoldKnk)
    set.seed(123)
    
    target_gene <- params$knockout_target
    message(sprintf("目标基因: %s", target_gene))
    
    # 检查基因是否存在
    if (!(target_gene %in% rownames(seurat_obj))) {
      warning(sprintf("基因 %s 不在数据集中", target_gene))
      return(NULL)
    }
    
    # 只保留肿瘤组（如果有Type分组）
    if ("Type" %in% colnames(seurat_obj@meta.data)) {
      tumor_types <- grep("tumor|disease|treat", unique(seurat_obj$Type),
                         ignore.case = TRUE, value = TRUE)
      if (length(tumor_types) > 0) {
        message(sprintf("筛选肿瘤样本: %s", tumor_types[1]))
        pbmc_subset <- subset(seurat_obj, subset = Type == tumor_types[1])
      } else {
        pbmc_subset <- seurat_obj
      }
    } else {
      pbmc_subset <- seurat_obj
    }
    
    # 提取表达矩阵
    message("提取表达矩阵...")
    countMat <- GetAssayData(pbmc_subset, slot = "counts")
    
    # 识别高变基因
    message("识别高变基因...")
    pbmc_subset <- FindVariableFeatures(object = pbmc_subset,
                                       selection.method = "vst",
                                       nfeatures = 10000)
    hvgs <- VariableFeatures(pbmc_subset)
    
    # 准备数据
    data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])
    
    message(sprintf("数据维度: %d 基因 × %d 细胞", nrow(data), ncol(data)))
    
    # 执行虚拟敲除
    message("执行虚拟敲除分析（这可能需要较长时间）...")
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
    
    message(sprintf("发现 %d 个显著差异基因", nrow(sig_diff)))
    
    # === 柱状图 ===
    message("生成柱状图...")
    top_genes <- head(df[order(-df$FC), ], 20)
    
    p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
      geom_bar(stat = 'identity', fill = '#5A9BD4') +
      coord_flip() +
      labs(title = sprintf("Top 20 Genes After %s Knockout", target_gene),
           x = "Gene", y = "Fold Change") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(p1,
          filename = file.path(output_dirs$module2_knockout, "Barplot.pdf"),
          width = 8, height = 6)
    
    # === 火山图 ===
    message("生成火山图...")
    df$log_p.adj <- -log10(df$p.adj)
    df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
    label_genes <- subset(df, p.adj < 0.05)
    
    y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)
    
    p2 <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
      geom_point(alpha = 0.7, size = 1.5) +
      scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      labs(title = sprintf("Volcano Plot: %s Knockout", target_gene),
           x = "Z-score", y = "-log10(p.adj)") +
      theme_classic() +
      coord_cartesian(ylim = c(0, y_upper)) +
      theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # 添加基因标签
    if (nrow(label_genes) > 0 & requireNamespace("ggrepel", quietly = TRUE)) {
      p2 <- p2 + ggrepel::geom_text_repel(
        data = head(label_genes, 20),
        aes(label = gene),
        size = 3,
        max.overlaps = 50
      )
    }
    
    ggsave(p2,
          filename = file.path(output_dirs$module2_knockout, "Volcano_Plot.pdf"),
          width = 10, height = 8)
    
    # === 饼图 ===
    message("生成饼图...")
    sig_count <- table(df$significant)
    sig_df <- as.data.frame(sig_count)
    colnames(sig_df) <- c("category", "count")
    sig_df$percentage <- paste0(round(sig_df$count / sum(sig_df$count) * 100, 1), "%")
    
    p3 <- ggplot(sig_df, aes(x = "", y = count, fill = category)) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), size = 5) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = c("Significant" = "red", "Not significant" = "lightgray")) +
      labs(title = sprintf("Proportion of Significant Genes (%s KO)", target_gene),
           fill = "") +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    ggsave(p3,
          filename = file.path(output_dirs$module2_knockout, "Pie_Chart.pdf"),
          width = 8, height = 6)
    
    # 保存结果对象
    saveRDS(result,
           file = file.path(output_dirs$module2_knockout, "scTenifoldKnk_Result.rds"))
    
    message("✓ 虚拟敲除分析完成")
    
    return(result)
    
  }, error = function(e) {
    warning(sprintf("虚拟敲除分析失败: %s", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.3: GO富集分析
# ==============================================================================

run_go_enrichment <- function(seurat_obj, params, output_dirs, gene_file = NULL) {
  message("\n=== 模块2.3: GO富集分析 ===")
  
  tryCatch({
    # 检查必要的包
    required_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
    for (pkg in required_pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        warning(sprintf("缺少包 %s，跳过GO富集分析", pkg))
        return(NULL)
      }
    }
    
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    
    # 读取基因列表
    if (!is.null(gene_file) && file.exists(gene_file)) {
      message(sprintf("从文件读取基因: %s", gene_file))
      rt <- read.table(gene_file, header = TRUE, sep = "\t", check.names = FALSE)
      genes <- unique(as.vector(rt[, 1]))
    } else if (!is.null(seurat_obj@misc$sig_diff_genes)) {
      message("使用差异基因进行GO分析...")
      genes <- seurat_obj@misc$sig_diff_genes$gene
    } else {
      warning("没有基因列表，跳过GO分析")
      return(NULL)
    }
    
    message(sprintf("基因数量: %d", length(genes)))
    
    # 基因ID转换
    message("转换基因ID...")
    entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
    entrezIDs <- as.character(entrezIDs)
    
    # 移除NA
    valid_ids <- entrezIDs[entrezIDs != "NA"]
    message(sprintf("有效Entrez ID: %d", length(valid_ids)))
    
    if (length(valid_ids) == 0) {
      warning("没有有效的Entrez ID")
      return(NULL)
    }
    
    # GO富集分析
    message("执行GO富集分析...")
    kk <- enrichGO(
      gene = valid_ids,
      OrgDb = org.Hs.eg.db,
      pvalueCutoff = 1,
      qvalueCutoff = 1,
      ont = "all",
      readable = TRUE
    )
    
    GO <- as.data.frame(kk)
    
    # 筛选显著结果
    GO_sig <- GO[GO$pvalue < params$go_pvalue_filter &
                GO$p.adjust < params$go_adj_pval_filter, ]
    
    message(sprintf("显著GO term: %d", nrow(GO_sig)))
    
    # 保存结果
    write.table(GO,
               file = file.path(output_dirs$module2_enrichment, "GO_All_Results.txt"),
               sep = "\t", quote = FALSE, row.names = FALSE)
    
    write.table(GO_sig,
               file = file.path(output_dirs$module2_enrichment, "GO_Significant_Results.txt"),
               sep = "\t", quote = FALSE, row.names = FALSE)
    
    # 判断用于颜色的列
    colorSel <- "p.adjust"
    if (params$go_adj_pval_filter > 0.05) {
      colorSel <- "pvalue"
    }
    
    # 柱状图
    message("生成柱状图...")
    pdf(file = file.path(output_dirs$module2_enrichment, "GO_Barplot.pdf"),
        width = 10, height = 8)
    bar <- barplot(kk, drop = TRUE, showCategory = 10, label_format = 100,
                  split = "ONTOLOGY", color = colorSel) +
      facet_grid(ONTOLOGY ~ ., scale = 'free')
    print(bar)
    dev.off()
    
    # 气泡图
    message("生成气泡图...")
    pdf(file = file.path(output_dirs$module2_enrichment, "GO_Bubble.pdf"),
        width = 10, height = 8)
    bub <- dotplot(kk, showCategory = 10, orderBy = "GeneRatio",
                  label_format = 100, split = "ONTOLOGY", color = colorSel) +
      facet_grid(ONTOLOGY ~ ., scale = 'free')
    print(bub)
    dev.off()
    
    # === GO圈图 ===
    if (nrow(GO_sig) > 0 && requireNamespace("circlize", quietly = TRUE)) {
      message("生成GO圈图...")
      
      library(circlize)
      library(ComplexHeatmap)
      
      ontology.col <- c("#6A89C2FF", "#D6616BFF", "#67B88BFF")
      
      # 准备数据
      data <- GO[order(GO$pvalue), ]
      datasig <- data[data$pvalue < 0.05, , drop = FALSE]
      
      BP <- datasig[datasig$ONTOLOGY == "BP", , drop = FALSE]
      CC <- datasig[datasig$ONTOLOGY == "CC", , drop = FALSE]
      MF <- datasig[datasig$ONTOLOGY == "MF", , drop = FALSE]
      
      BP <- head(BP, 6)
      CC <- head(CC, 6)
      MF <- head(MF, 6)
      
      data <- rbind(BP, CC, MF)
      
      if (nrow(data) > 0) {
        main.col <- ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
        
        # 准备圈图数据
        BgGene <- as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
        Gene <- as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
        ratio <- Gene / BgGene
        logpvalue <- -log(data$pvalue, 10)
        
        logpvalue.col <- brewer.pal(n = 6, name = "Reds")
        f <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10), colors = logpvalue.col)
        BgGene.col <- f(logpvalue)
        
        df <- data.frame(GO = data$ID, start = 1, end = max(BgGene))
        rownames(df) <- df$GO
        
        bed2 <- data.frame(GO = data$ID, start = 1, end = BgGene,
                          BgGene = BgGene, BgGene.col = BgGene.col)
        bed3 <- data.frame(GO = data$ID, start = 1, end = Gene, BgGene = Gene)
        bed4 <- data.frame(GO = data$ID, start = 1, end = max(BgGene),
                          ratio = ratio, col = main.col)
        bed4$ratio <- bed4$ratio / max(bed4$ratio) * 9.5
        
        # 绘制圈图
        pdf(file = file.path(output_dirs$module2_enrichment, "GO_Circlize.pdf"),
            width = 10, height = 10)
        
        par(omi = c(0.1, 0.1, 0.1, 1.5))
        circos.par(track.margin = c(0.01, 0.01))
        circos.genomicInitialize(df, plotType = "none")
        
        circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
          sector.index <- get.cell.meta.data("sector.index")
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")
          circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8,
                     facing = "bending.inside", niceFacing = TRUE)
        }, track.height = 0.08, bg.border = NA, bg.col = main.col)
        
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
      }
    }
    
    message("✓ GO富集分析完成")
    
    return(kk)
    
  }, error = function(e) {
    warning(sprintf("GO富集分析失败: %s", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.4: 基因可视化
# ==============================================================================

visualize_gene_of_interest <- function(seurat_obj, params, output_dirs) {
  message("\n=== 模块2.4: 基因可视化 ===")
  
  tryCatch({
    gene_name <- params$gene_of_interest
    message(sprintf("可视化基因: %s", gene_name))
    
    # 检查基因是否存在
    if (!(gene_name %in% rownames(seurat_obj))) {
      warning(sprintf("基因 %s 不在数据集中", gene_name))
      return(NULL)
    }
    
    # 散点图（FeaturePlot）
    message("生成散点图...")
    pdf(file = file.path(output_dirs$module2_gene_view, "Gene_FeaturePlot.pdf"),
        width = 8, height = 7)
    p1 <- FeaturePlot(object = seurat_obj, features = gene_name,
                     cols = c("lightgrey", "red"), pt.size = 1.5) +
      ggtitle(sprintf("%s Expression", gene_name)) +
      theme_minimal()
    print(p1)
    dev.off()
    
    # 小提琴图
    message("生成小提琴图...")
    pdf(file = file.path(output_dirs$module2_gene_view, "Gene_ViolinPlot.pdf"),
        width = 10, height = 6)
    p2 <- VlnPlot(object = seurat_obj, features = gene_name, split.by = "Type",
                 pt.size = 0.1) +
      ggtitle(sprintf("%s Expression by Group", gene_name))
    print(p2)
    dev.off()
    
    # 点图
    message("生成点图...")
    pdf(file = file.path(output_dirs$module2_gene_view, "Gene_DotPlot.pdf"),
        width = 8, height = 6)
    p3 <- DotPlot(object = seurat_obj, features = gene_name) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      ggtitle(sprintf("%s Expression Across Cell Types", gene_name))
    print(p3)
    dev.off()
    
    # 如果有Type分组，生成分组比较图
    if ("Type" %in% colnames(seurat_obj@meta.data)) {
      pdf(file = file.path(output_dirs$module2_gene_view, "Gene_FeaturePlot_Split.pdf"),
          width = 12, height = 5)
      p4 <- FeaturePlot(object = seurat_obj, features = gene_name,
                       split.by = "Type", cols = c("lightgrey", "red"),
                       pt.size = 1) +
        ggtitle(sprintf("%s Expression by Group", gene_name))
      print(p4)
      dev.off()
    }
    
    # 岭线图（如果有ggridges包）
    if (requireNamespace("ggridges", quietly = TRUE)) {
      library(ggridges)
      
      plot_data <- data.frame(
        Expression = GetAssayData(seurat_obj, slot = "data")[gene_name, ],
        CellType = Idents(seurat_obj)
      )
      
      pdf(file = file.path(output_dirs$module2_gene_view, "Gene_RidgePlot.pdf"),
          width = 10, height = 8)
      p5 <- ggplot(plot_data, aes(x = Expression, y = CellType, fill = CellType)) +
        geom_density_ridges() +
        theme_minimal() +
        labs(title = sprintf("%s Expression Distribution", gene_name),
             x = "Expression Level", y = "Cell Type") +
        theme(legend.position = "none")
      print(p5)
      dev.off()
    }
    
    message("✓ 基因可视化完成")
    
    return(TRUE)
    
  }, error = function(e) {
    warning(sprintf("基因可视化失败: %s", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2.5: 高级可视化（综合图表）
# ==============================================================================

create_advanced_visualizations <- function(seurat_obj, output_dirs) {
  message("\n=== 模块2.5: 高级可视化 ===")
  
  tryCatch({
    # 细胞类型标记基因气泡图
    message("生成细胞类型标记基因气泡图...")
    
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
      pdf(file = file.path(output_dirs$module2_advanced_plots, "CellType_Markers_DotPlot.pdf"),
          width = 14, height = 9)
      p1 <- DotPlot(seurat_obj, features = all_markers, dot.scale = 8, scale = TRUE,
                   cols = c("lightgrey", "#002868", "#BF0A30")) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
          axis.text.y = element_text(size = 10, face = "bold")
        ) +
        labs(title = "Cell Type Marker Expression Profile")
      print(p1)
      dev.off()
    }
    
    # 细胞比例堆叠图
    if ("Type" %in% colnames(seurat_obj@meta.data)) {
      message("生成细胞比例堆叠图...")
      
      prop_data <- as.data.frame(prop.table(table(Idents(seurat_obj), seurat_obj$Type), margin = 2) * 100)
      colnames(prop_data) <- c("CellType", "Group", "Percentage")
      
      pdf(file = file.path(output_dirs$module2_advanced_plots, "Cell_Proportion_Stacked.pdf"),
          width = 10, height = 7)
      p2 <- ggplot(prop_data, aes(x = Group, y = Percentage, fill = CellType)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        labs(title = "Cell Type Proportions by Group",
             x = "Group", y = "Percentage (%)", fill = "Cell Type") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      print(p2)
      dev.off()
    }
    
    # UMAP美化版
    message("生成美化版UMAP...")
    pdf(file = file.path(output_dirs$module2_advanced_plots, "UMAP_Enhanced.pdf"),
        width = 12, height = 10)
    p3 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 1.5,
                 label.size = 5, repel = TRUE) +
      ggtitle("Enhanced UMAP Visualization") +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        legend.position = "right"
      )
    print(p3)
    dev.off()
    
    message("✓ 高级可视化完成")
    
    return(TRUE)
    
  }, error = function(e) {
    warning(sprintf("高级可视化失败: %s", e$message))
    return(NULL)
  })
}

# ==============================================================================
# Module 2: 主执行函数
# ==============================================================================

run_module2_analysis <- function(workDir, params, output_dirs, tracker) {
  
  tracker$set_module("Module 2: Advanced Analysis")
  
  # 加载Seurat对象
  seurat_file <- file.path(output_dirs$rdata_storage, "Seurat.RData")
  
  if (!file.exists(seurat_file)) {
    stop("未找到Seurat.RData文件，请先运行模块1")
  }
  
  message("加载Seurat对象...")
  load(seurat_file)
  seurat_obj <- pbmc
  
  message(sprintf("✓ 加载完成: %d 细胞 × %d 基因", ncol(seurat_obj), nrow(seurat_obj)))
  
  # 1. CellChat分析
  cellchat_result <- safe_execute(run_cellchat_analysis, "cellchat_analysis",
                                  tracker, seurat_obj = seurat_obj,
                                  params = params, output_dirs = output_dirs)
  
  # 2. 虚拟敲除分析
  knockout_result <- safe_execute(run_virtual_knockout, "virtual_knockout",
                                 tracker, seurat_obj = seurat_obj,
                                 params = params, output_dirs = output_dirs)
  
  # 3. GO富集分析
  # 首先尝试使用虚拟敲除结果
  if (!is.null(knockout_result)) {
    sig_gene_file <- file.path(output_dirs$module2_knockout, "Significant_Diff_Genes.txt")
    if (file.exists(sig_gene_file)) {
      go_result <- safe_execute(run_go_enrichment, "go_enrichment", tracker,
                               seurat_obj = seurat_obj, params = params,
                               output_dirs = output_dirs, gene_file = sig_gene_file)
    }
  } else {
    # 使用模块1的差异基因
    go_result <- safe_execute(run_go_enrichment, "go_enrichment", tracker,
                             seurat_obj = seurat_obj, params = params,
                             output_dirs = output_dirs, gene_file = NULL)
  }
  
  # 4. 基因可视化
  gene_viz_result <- safe_execute(visualize_gene_of_interest, "gene_visualization",
                                  tracker, seurat_obj = seurat_obj,
                                  params = params, output_dirs = output_dirs)
  
  # 5. 高级可视化
  advanced_viz_result <- safe_execute(create_advanced_visualizations,
                                      "advanced_visualization", tracker,
                                      seurat_obj = seurat_obj,
                                      output_dirs = output_dirs)
  
  message("\n✓✓✓ 模块2：高级分析完成 ✓✓✓\n")
  
  return(list(
    cellchat = cellchat_result,
    knockout = knockout_result,
    go_enrichment = go_result,
    gene_viz = gene_viz_result,
    advanced_viz = advanced_viz_result
  ))
}
