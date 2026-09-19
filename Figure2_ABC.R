# Figure2_ABC.R
# Mimicking Lauer et al. 2024 Figure 2.ABC style for Thrombus Proteomics
# Using preprocessing logic from thrombus_figure_generator.R

# ------------------------------------------------------------------------------
# 1. Setup and Configuration
# ------------------------------------------------------------------------------

# Create output directory
output_dir <- "lauer_style_comparison/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# CLEANUP: Close any dangling graphics devices (like open PDFs from previous failed runs)
while (!is.null(dev.list())) dev.off()

# Load libraries
suppressPackageStartupMessages({
  library(limma)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(tidyr)
  library(ReactomePA)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(tidyverse)
  library(Mfuzz)
  library(e1071)
  library(patchwork)
  library(UpSetR) # Added for Shared GO UpSet plot
})

# Source helper scripts
source("src/data_loader.R")

# ------------------------------------------------------------------------------
# 2. Data Loading and Preprocessing (Exact logic from thrombus_figure_generator.R)
# ------------------------------------------------------------------------------

message(">>> Loading and Preprocessing Data (Figure Generator Logic)...")

input_file <- "thrombusDIAreport.xls" 
if (!file.exists(input_file)) {
  input_file <- "c:/Users/SimonYao/Desktop/LCM_protemics_thrombus/thrombusDIAreport.xls"
}

data_list <- load_and_clean_data(input_file)
expr_data <- data_list$exprs
feature_data <- data_list$feature_data
metadata <- data_list$metadata

# Filter out samples with missing groups
valid_samples <- !is.na(metadata$Group)
expr_data <- expr_data[, valid_samples]
metadata <- metadata[valid_samples, ]

# 2.1 Filtering
groups <- unique(metadata$Group)
keep_protein <- rep(FALSE, nrow(expr_data))
for (g in groups) {
  samples_in_group <- metadata$SampleID[metadata$Group == g]
  samples_in_group <- intersect(samples_in_group, colnames(expr_data))
  if (length(samples_in_group) > 0) {
    valid_counts <- rowSums(expr_data[, samples_in_group] > 0 & !is.na(expr_data[, samples_in_group]))
    keep_protein <- keep_protein | (valid_counts >= (0.5 * length(samples_in_group)))
  }
}
expr_filtered <- expr_data[keep_protein, ]
feature_filtered <- feature_data[keep_protein, ]

# 2.2 Standardize Gene Symbols for consistent mapping across all plots
feature_filtered$Symbol <- sapply(strsplit(feature_filtered$PG.Genes, ";"), `[`, 1)
feature_filtered$Symbol <- stringr::str_to_title(feature_filtered$Symbol)
feature_filtered$UniqueSymbol <- make.unique(feature_filtered$Symbol)

# 2.2 Log2 Transformation & Imputation
expr_filtered[expr_filtered == 0] <- NA
expr_log <- log2(expr_filtered)

impute_minprob <- function(data, shift = 1.8, width = 0.3) {
  valid_vals <- data[!is.na(data)]
  mu <- mean(valid_vals)
  sigma <- sd(valid_vals)
  mu_imp <- mu - (shift * sigma)
  sigma_imp <- sigma * width
  imputed_data <- data
  for (i in 1:ncol(data)) {
    n_missing <- sum(is.na(data[, i]))
    if (n_missing > 0) {
      imputed_data[is.na(data[, i]), i] <- rnorm(n_missing, mean = mu_imp, sd = sigma_imp)
    }
  }
  return(imputed_data)
}

set.seed(123)
expr_imputed <- impute_minprob(expr_log)
expr_norm <- normalizeQuantiles(as.matrix(expr_imputed))
rownames(expr_norm) <- rownames(expr_filtered)

# ------------------------------------------------------------------------------
# 3. Differential Expression Analysis (limma)
# ------------------------------------------------------------------------------
# 3. DEA & Global Definitions
# ------------------------------------------------------------------------------

# Define Global Cluster Color Palette (Set1)
cluster_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
names(cluster_palette) <- paste0("Cluster_", 1:6)

message(">>> Performing DEA...")

# Define Sample Order as requested: D2-F, D7-F, D7-C, D14-F, D14-C, D14-Wall
target_order <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall")
all_present_groups <- unique(metadata$Group)
# Keep requested order, then append any others (like 'H') at the end
final_levels <- c(intersect(target_order, all_present_groups), 
                  setdiff(all_present_groups, target_order))

group <- factor(metadata$Group, levels = final_levels)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_norm, design)

# Expanded Comparisons (excluding Healthy)
contrast_matrix <- makeContrasts(
  # Timepoint Transitions
  D2_F_vs_D7_F = D7_Fibrin - D2_Fibrin,
  D7_F_vs_D14_F = D14_Fibrin - D7_Fibrin,
  D7_C_vs_D14_C = D14_Collagen - D7_Collagen,
  
  # Regional Heterogeneity
  D7_C_vs_D7_F = D7_Collagen - D7_Fibrin,
  D14_C_vs_D14_F = D14_Collagen - D14_Fibrin,
  
  # Resolution Failure
  D14_F_vs_D7_F = D14_Fibrin - D7_Fibrin,
  D14_F_vs_D14_C = D14_Fibrin - D14_Collagen,
  D14_F_vs_D2_F = D14_Fibrin - D2_Fibrin,
  
  # Wall Interaction
  D14_Wall_vs_D14_Collagen = D14_Wall - D14_Collagen,
  
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# 3.1 Get conditions from factor levels to maintain order
message(">>> Calculating logFC vs Control for heatmaps...")
all_groups <- levels(group)
conditions <- all_groups
message(">>> Conditions (in order): ", paste(conditions, collapse=", "))

# Group Mapping for comparisons (used in heatmaps)
group_map <- list(
  D2_F_vs_D7_F = c("D7_Fibrin", "D2_Fibrin"),
  D7_F_vs_D14_F = c("D14_Fibrin", "D7_Fibrin"),
  D7_C_vs_D14_C = c("D14_Collagen", "D7_Collagen"),
  D7_C_vs_D7_F = c("D7_Collagen", "D7_Fibrin"),
  D14_C_vs_D14_F = c("D14_Collagen", "D14_Fibrin"),
  D14_F_vs_D7_F = c("D14_Fibrin", "D7_Fibrin"),
  D14_F_vs_D14_C = c("D14_Fibrin", "D14_Collagen"),
  D14_F_vs_D2_F = c("D14_Fibrin", "D2_Fibrin"),
  D14_Wall_vs_D14_Collagen = c("D14_Wall", "D14_Collagen")
)

# Note: go_reference, gene_modules, and thrombus_go_modules removed by user request

# ------------------------------------------------------------------------------
# 4. Visualization Functions (100% Lauer Style)
# ------------------------------------------------------------------------------

# 4.1 Figure 2A: Volcano Plot
plot_lauer_volcano <- function(res, title) {
  # Label proteins
  res$expression <- "none"
  # Lauer uses 0.3 for Log2FC threshold in some plots, but 0.58 (1.5-fold) is more standard.
  # Based on Lauer's script: de.proteins.c1.c2.sig <- de.proteins.c1.c2[de.proteins.c1.c2$P.Value < 0.05 & abs(de.proteins.c1.c2$logFC) > 0.3, ]
  # We'll stick to 0.3 as requested "100% Lauer's style"
  res$expression[res$logFC > 0.3 & res$P.Value < 0.05] <- "up"
  res$expression[res$logFC < -0.3 & res$P.Value < 0.05] <- "down"
  
  res$symbol <- feature_filtered[rownames(res), "PG.Genes"]
  res$symbol <- sapply(strsplit(res$symbol, ";"), `[`, 1)
  
  # Top 20 for labeling
  res$delabel <- ifelse(res$symbol %in% head(res[order(res$P.Value), "symbol"], 20), res$symbol, NA)
  
  ggplot(data = res, aes(x = logFC, y = -log10(P.Value), color = expression, label = delabel)) +
    geom_point() +
    geom_vline(xintercept = c(-0.3, 0.3), col = "#000000", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "#000000", linetype = "dashed") +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("down" = "#2B3990", "none" = "gray", "up" = "#BE1E2D")) +
    labs(title = title, x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    theme(text = element_text(size = 15, color = "#000000"),
          legend.title = element_blank(),
          axis.line = element_line(colour = "#000000", size = 1), 
          legend.position = "none",
          panel.background = element_blank(),
          aspect.ratio = 1/1)
}

# 4.2 Figure 2B: Pathway Enrichment
plot_lauer_enrichment <- function(res, title, direction = "up", ref_terms = NULL) {
  symbols <- feature_filtered[rownames(res), "PG.Genes"]
  symbols <- sapply(strsplit(symbols, ";"), `[`, 1)
  symbols <- stringr::str_to_title(symbols)
  
  gene_map <- tryCatch({
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) NULL)
  
  if (is.null(gene_map)) return(NULL)
  
  res$symbol <- symbols
  res_mapped <- merge(res, gene_map, by.x = "symbol", by.y = "SYMBOL")
  
  if (direction == "up") {
    sig_genes <- res_mapped$ENTREZID[res_mapped$logFC > 0.3 & res_mapped$P.Value < 0.05]
    plot_title <- paste(title, "(Upregulated)")
  } else {
    sig_genes <- res_mapped$ENTREZID[res_mapped$logFC < -0.3 & res_mapped$P.Value < 0.05]
    plot_title <- paste(title, "(Downregulated)")
  }
  
  if (length(sig_genes) < 5) return(NULL)
  
  pathway.go <- enrichGO(gene = sig_genes, 
                         ont = "ALL", pvalueCutoff=0.05, OrgDb = org.Mm.eg.db, pAdjustMethod = "fdr", keyType = "ENTREZID")
  
  if (is.null(pathway.go) || nrow(as.data.frame(pathway.go)) == 0) return(NULL)
  
  pathway.go <- setReadable(pathway.go, OrgDb = org.Mm.eg.db)
  go_df <- as.data.frame(pathway.go)
  
  # Filtering logic based on reference terms
  show_cat <- 15
  if (!is.null(ref_terms)) {
    # Match by ID or Description
    matched_ids <- go_df$ID[go_df$ID %in% ref_terms]
    matched_desc <- go_df$Description[grepl(paste(ref_terms, collapse = "|"), go_df$Description, ignore.case = TRUE)]
    
    selected_pathways <- unique(c(matched_ids, matched_desc))
    
    if (length(selected_pathways) >= 15) {
      show_cat <- selected_pathways[1:15]
    } else if (length(selected_pathways) > 0) {
      # Fill with top significant ones if not enough
      other_pathways <- setdiff(go_df$Description, selected_pathways)
      show_cat <- c(selected_pathways, head(other_pathways, 15 - length(selected_pathways)))
    }
  }
  
  p_go <- ReactomePA::dotplot(pathway.go, showCategory = show_cat, title = plot_title) +
    theme(aspect.ratio = 2/1,
          plot.title = element_text(face = "bold"))
  
  return(list(go = p_go, results = go_df))
}

# 4.3 Figure 2C: Heatmap
plot_lauer_heatmap <- function(res, title, filename = NULL) {
  sig_res <- res[res$P.Value < 0.05 & abs(res$logFC) > 0.3, ]
  if (nrow(sig_res) < 5) {
    message("    Too few genes for heatmap: ", title)
    return(NULL)
  }
  
  top_genes <- head(sig_res[order(sig_res$P.Value), ], 50)
  
  plot_data <- data.frame(logFC = top_genes$logFC)
  rownames(plot_data) <- feature_filtered[rownames(top_genes), "PG.Genes"]
  rownames(plot_data) <- sapply(strsplit(rownames(plot_data), ";"), `[`, 1)
  
  p <- pheatmap(as.matrix(plot_data), 
           scale = "none", 
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
           main = title,
           fontsize_row = 8,
           cluster_cols = FALSE,
           silent = TRUE)
  return(p)
}

# 4.6 Figure 2D: Multi-Volcano Plot (User Requested)

# 4.6 Figure 2D: Multi-Volcano Plot (User Requested)
mutiVolcano = function(df,         # 绘图数据
                       P = 0.05,   # P值卡值
                       FC = 1.5,   # FC卡值
                       GroupName = c("Sig","Not Sig"),      # 分组标签
                       pointColor = c("#CC3333","#0099CC"), # 分组散点的颜色
                       barFill = "#efefef",  # 柱子的颜色
                       pointSize = 0.9,      # 散点的大小
                       labeltype = "1",      # 标记差异基因的选项，标记类型有"1"和"2"两种选项
                       labelNum = 5,         # 当标记类型为1时，待标记的散点个数
                       labelName =NULL,      # 当标记类型为2时，待标记的散点名称
                       tileLabel =  "Label", # 标记比较对的选项，选项有“Label”和“Num”，Label时显示分组名称，Num时显示数字，防止因为标签太长导致的不美观
                       tileColor = NULL      # 比较对的颜色
                       ){
  # 数据分组 根据p的卡值分组
  dfSig = df %>% 
    mutate(log2FC = log2(FC)) %>%
    filter(FC > {{FC}} | FC < (1/{{FC}})) %>%
    mutate(Group = ifelse(PValue < 0.05, GroupName[[1]], GroupName[[2]])) %>%
    mutate(Group = factor(Group, levels = GroupName)) %>%
    mutate(Cluster = factor(Cluster, levels = unique(Cluster)))   # Cluster的顺序是文件中出现的顺序
  
  # 柱形图数据整理
  dfBar = dfSig %>%
    group_by(Cluster) %>%
    summarise(min = min(log2FC, na.rm = T),
              max = max(log2FC, na.rm = T)
              )
  # 散点图数据整理
  dfJitter = dfSig %>%
    mutate(jitter = jitter(as.numeric(Cluster), factor = 2))
  
  # 整理标记差异基因的数据
  if(labeltype == "1"){
    # 标记一
    # 每组P值最小的几个
    dfLabel = dfJitter %>%
      group_by(Cluster) %>%
      slice_min(PValue, n = labelNum, with_ties = F) %>%
      ungroup()
  }else if(labeltype == "2"){
    # 标记二
    # 指定标记
    dfLabel = dfJitter %>%
      filter(Name %in% labelName)
  }else{
    dfLabel = dfJitter %>% slice()
  }
    
  # 绘图
  p = ggplot()+
    # 绘制柱形图
    geom_col(data = dfBar, aes(x = Cluster, y = max), fill = barFill)+
    geom_col(data = dfBar, aes(x = Cluster, y = min), fill = barFill)+
    # 绘制散点图
    geom_point(data = dfJitter,
               aes(x = jitter, y = log2FC, color = Group),
               size = pointSize,
               show.legend = NA
               )+
    # 绘制中间的标签方块
    ggplot2::geom_tile(data = dfSig,
                       ggplot2::aes(x = Cluster, y = 0, fill = Cluster), 
                       color = "black",
                       height = log2(FC) * 1.5,
                       # alpha = 0.3,
                       show.legend = NA
                       ) + 
    # 标记差异基因
    ggrepel::geom_text_repel(
      data = dfLabel,
      aes(x = jitter,                   # geom_text_repel 标记函数
          y = log2FC,          
          label = Name),        
      min.segment.length = 0.1,
      max.overlaps = 10000,                    # 最大覆盖率，当点很多时，有些标记会被覆盖，调大该值则不被覆盖，反之。
      size = 3,                                  # 字体大小
      box.padding = unit(0.5, 'lines'),           # 标记的边距
      point.padding = unit(0.1, 'lines'), 
      segment.color = 'black',                   # 标记线条的颜色
      show.legend = F)#+
  
  if(tileLabel == "Label"){
     p =
      p +
      geom_text(data = dfSig, aes(x = Cluster, y = 0, label = Cluster))+
      ggplot2::scale_fill_manual(values = tileColor,
                                 guide = NULL # 不显示该图例
                                 )
  }else if(tileLabel == "Num"){
    # 如果比较对的名字太长，可以改成数字标签
    p =
      p +
      geom_text(data = dfSig, aes(x = Cluster, y = 0, label = as.numeric(Cluster)), show.legend = NA)+
      ggplot2::scale_fill_manual(values = tileColor,
                                 labels = c(paste0(1:length(unique(dfSig$Cluster)), ": ", unique(dfSig$Cluster))))
  }

  
    
  # 修改主题
  p = p + ggplot2::scale_color_manual(values = pointColor)+
    theme_classic()+
    ggplot2::scale_y_continuous(n.breaks = 5) + 
    ggplot2::theme(
                   legend.position = "right", 
                   legend.title = ggplot2::element_blank(), 
                   legend.background = ggplot2::element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   # axis.title.x = element_blank(),
                   axis.line.x = element_blank()
                   ) + 
    ggplot2::xlab("Comparisons") + ggplot2::ylab("log2FC") + 
    # ggplot2::guides(fill = ggplot2::guide_legend())
    guides(color = guide_legend(override.aes = list(size = 3)))
    
    return(p)
}

# 4.7 Figure 2E: Mfuzz Heatmap (User Requested)
mfuzzHeatmap = function(data,
                        clusterNum=6
                        ){
  # data = df
  # clusterNum=6

  # 构建对象，标准化等
  dm <- data.matrix(data)                  # 数据框转换为矩阵
  
  # Manual standardization (mean=0, sd=1 per gene) to avoid "standardise" function issues
  # This is equivalent to Mfuzz::standardise but more robust
  dm_std <- t(scale(t(dm)))
  
  ESet <- new("ExpressionSet", exprs = dm_std)  # 构建对象
  
  # Data is already imputed and filtered in the main script
  # ESet <- Mfuzz::filter.NA(ESet, thres=0.25)      
  # ESet <- Mfuzz::fill.NA(ESet,mode="knn")         
  # ESet <- Mfuzz::filter.std(ESet,min.std=0,visu=F)
  
  gene.s <- ESet # Already standardised manually
  # exprs(gene.s)                            # 查看处理后的数据

  # 聚类
  c <- clusterNum                 # 设置聚类个数
  m <- Mfuzz::mestimate(gene.s) # 评估出最佳的m值
  set.seed(123)          # 设置随机种子，防止每次聚类的结果都不一样，无法复现
  cl <- Mfuzz::mfuzz(gene.s, c = c, m = m)
  # cl                     # 查看每个基因聚到哪个类当中
  # cl$size                # 查看每个cluster中的基因个数
  # cl$membership           #查看基因和cluster之间的membership。如果两个基因对于一个特定的cluster都有高的membership score，那么他们通常来说表达模式是相似的


  # ggplot2绘图
  # 绘制趋势分析
  dfMfuzz = exprs(gene.s) %>%data.frame()
  # dfMfuzz
  dfColor = cl$membership %>%
    data.frame(check.names = F) %>%
    tibble::rownames_to_column("ID") %>%
    pivot_longer(-1,names_to = "cluster",values_to = "Membership")
  # 计算color的范围，让多张图的图例保持相同
  colorLimit =  dfColor %>%
    group_by(ID) %>%
    summarise(max=max(Membership,na.rm = T)) %>%
    ungroup()

  dfcluster = cl$cluster %>%
    data.frame() %>%
    set_names("Cluster") %>%
    tibble::rownames_to_column("ID")

  dfclusterSplit = dfcluster %>%
    group_split(Cluster,.keep=T)

  mfuzzPlotList  = imap(dfclusterSplit,function(dataItem,i){
    # dataItem = dfclusterSplit[[1]]

    clusterName = dataItem$Cluster[[1]]
    clusterID = dataItem %>% pull(ID)
    myN = length(clusterID)

    dfPlot = dfMfuzz[clusterID,] %>%
      tibble::rownames_to_column("ID") %>%
      pivot_longer(-1,names_to = "Sample",values_to = "Value") %>%
      mutate(Sample = factor(Sample, levels = colnames(dfMfuzz))) %>%
      left_join(dfColor %>%
                  filter(cluster==clusterName), by = "ID") %>%
      arrange(Membership)

    # dfPlot = dfPlot %>%
    #   filter(Membership>0.5)

    p = ggplot(dfPlot,aes(x=Sample,y=Value,group=factor(ID,levels = unique(ID)),color=Membership))+
      # Background lines with alpha for a "smoother" collective look
      geom_line(alpha = 0.4, size = 0.5) +
      # Bold centroid line (average trend)
      stat_summary(aes(group=1), fun=mean, geom="line", size=1.5, color="black") +
      scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), # Publication standard palette
                            breaks=seq(0,1,0.2),
                            limits=c(0, 1)
      )+
      theme_classic()+ # Cleaner theme
      labs(y=paste0("cluster ",clusterName,"\n","n=",myN),
           x="")+
      theme(legend.position = "none", # Hide individual legends to clean up
            panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0), "inches")
      )+
      scale_x_discrete(expand = c(0.05,0.05))
    
    if(i==length(dfclusterSplit)){
      p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }else{
      p = p+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank())
    }
    return(p)
  })

  # 绘制热图
  heatmapList  = imap(dfclusterSplit,function(dataItem,i){
    # dataItem = dfclusterSplit[[1]]
    clusterName = dataItem$Cluster[[1]]
    clusterID = dataItem %>% pull(ID)

    dfMean.Zscore = t(data) %>%
      scale() %>%
      t() %>%
      data.frame()

    dfPlot = dfMean.Zscore[clusterID,] %>%
      tibble::rownames_to_column("ID") %>%
      pivot_longer(-1,names_to = "Sample",values_to = "Value") %>%
      mutate(Sample = factor(Sample, levels = colnames(dfMean.Zscore)))
    
    p = ggplot(dfPlot,aes(x=Sample,y=ID,fill=Value))+
      geom_tile()+
      # High-contrast RdBu palette for publication
      scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                           limits = c(-2.5, 2.5), oob = scales::squish) +
      theme_minimal()+
      labs(y="",x="",fill="z-score")+
      theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = unit(c(0.1, 0, 0.1, 0), "inches")
      )+
      scale_x_discrete(expand = c(0,0))

    if(i==length(dfclusterSplit)){
      p = p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }else{
      p = p+
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank())
    }
    return(p)
  })

  # 合并
  p = wrap_plots(c(mfuzzPlotList,heatmapList),byrow=F,ncol=2)+
    plot_layout(guides = 'collect') & theme(legend.position='top')
  return(
    list(p=p,
         dfcluster=dfcluster)
  )
}

# 4.8 Figure 2F: Cluster Enrichment Analysis (User Requested)
plot_cluster_enrichment = function(genes, cluster_name, output_dir) {
  message(">>> Running Enrichment for ", cluster_name, "...")
  
  # 1. Map symbols to Entrez IDs
  # Fix: Convert to Title Case for Mouse symbols (e.g., Vwf instead of VWF)
  genes_fixed <- stringr::str_to_title(genes)
  gene_map <- tryCatch({
    bitr(genes_fixed, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  }, error = function(e) NULL)
  
  if (is.null(gene_map) || nrow(gene_map) == 0) {
    message("    No Entrez IDs found for ", cluster_name)
    return(NULL)
  }
  
  entrez_ids <- gene_map$ENTREZID
  
  # 2. GO Enrichment (BP, CC, MF)
  ego <- enrichGO(gene          = entrez_ids,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego)
    
    # Select top 10 for each category
    ego_top <- ego_df %>%
      group_by(ONTOLOGY) %>%
      slice_min(pvalue, n = 15, with_ties = FALSE) %>%
      ungroup()
    
    # Map Ontology names for display
    ego_top$ONTOLOGY <- factor(ego_top$ONTOLOGY, 
                               levels = c("BP", "CC", "MF"),
                               labels = c("Biological Process", "Cell Component", "Molecular Function"))
    
    # Create the plot mimicking the user's example
    p_go <- ggplot(ego_top, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue), fill = ONTOLOGY)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
      scale_fill_manual(values = c("Biological Process" = "#91CF60", 
                                   "Cell Component" = "#FC8D59", 
                                   "Molecular Function" = "#4575B4")) +
      theme_bw() +
      labs(title = paste("GO Enrichment:", cluster_name),
           x = "",
           y = "-log10(P-value)") +
      theme(legend.position = "none",
            strip.text.y = element_text(angle = 0, face = "bold"),
            axis.text.y = element_text(size = 8))
    
    # Return plots instead of saving
    return(list(go = p_go, results = ego_df))
  }
  return(NULL)
}

# 4.8.5 Custom Heatmap for Top Proteins
plot_custom_heatmap = function(data, metadata, dfcluster, pathway_name, cluster_palette, mode = "individual") {
  # Force data to be a matrix
  data_mat <- as.matrix(data)
  
  # 1. Z-score normalization
  cal_z_score <- function(x) { 
    if(sd(x, na.rm=TRUE) == 0) return(x - mean(x, na.rm=TRUE))
    (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) 
  }
  data_z <- t(apply(data_mat, 1, cal_z_score))
  rownames(data_z) <- rownames(data_mat)
  colnames(data_z) <- colnames(data_mat)
  
  # 2. Column Annotation (Groups)
  target_order_levels <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall", "H")
  samp_match <- match(make.names(colnames(data_mat)), make.names(metadata$SampleID))
  curr_group <- metadata$Group[samp_match]
  
  if (mode == "average") {
    # Calculate group means
    data_df <- data.frame(t(data_z))
    data_df$SampleID_ForJoin <- make.names(colnames(data_mat))
    
    # Map groups using robust matching
    data_df$Group <- metadata$Group[match(data_df$SampleID_ForJoin, make.names(metadata$SampleID))]
    
    avg_data <- data_df %>%
      filter(!is.na(Group)) %>%
      group_by(Group) %>%
      summarise(across(-SampleID_ForJoin, \(x) mean(x, na.rm=TRUE)), .groups = "drop") %>%
      mutate(Group = factor(Group, levels = target_order_levels)) %>%
      arrange(Group)
    
    # Set back to matrix
    plot_matrix <- as.matrix(avg_data[,-1])
    # IMPORTANT: rownames of data_z are the proteins. summarising columns makes them column names.
    # We transposed at the start (data_df is t(data_z)), so avg_data columns ARE proteins.
    # Transpose back so rows are proteins.
    plot_matrix <- t(plot_matrix)
    colnames(plot_matrix) <- as.character(avg_data$Group)
    rownames(plot_matrix) <- rownames(data_z)
    
    annotation_col <- data.frame(Group = factor(avg_data$Group, levels = target_order_levels))
    rownames(annotation_col) <- colnames(plot_matrix)
    data_z <- plot_matrix
    gaps_col <- NULL
    show_colnames <- TRUE
  } else {
    # Sort samples by group order (Individual mode)
    meta_ordered <- data.frame(SampleID = colnames(data_mat), Group = curr_group) %>%
      mutate(Group = factor(Group, levels = target_order_levels)) %>%
      arrange(Group) %>%
      filter(!is.na(Group))
    
    data_z <- data_z[, meta_ordered$SampleID, drop=FALSE]
    annotation_col <- data.frame(Group = meta_ordered$Group)
    rownames(annotation_col) <- meta_ordered$SampleID
    gaps_col <- which(diff(as.numeric(annotation_col$Group)) != 0)
    show_colnames <- FALSE
  }
  
  # 3. Row Annotation (Clusters)
  row_annot <- dfcluster[match(rownames(data_z), dfcluster$ID), "Cluster", drop=FALSE]
  row_annot$Cluster <- factor(paste0("Cluster_", row_annot$Cluster))
  rownames(row_annot) <- rownames(data_z)
  
  # 4. Colors
  ann_colors <- list(
    Group = setNames(RColorBrewer::brewer.pal(length(levels(annotation_col$Group)), "Set2"), levels(annotation_col$Group)),
    Cluster = cluster_palette
  )

  # 5. Plotting
  p <- pheatmap(data_z, 
           scale = "none",
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
           annotation_col = annotation_col,
           annotation_row = row_annot,
           annotation_colors = ann_colors,
           show_colnames = show_colnames,
           cluster_cols = FALSE, 
           gaps_col = gaps_col,
           cluster_rows = (nrow(data_z) > 1),
           main = pathway_name,
           border_color = NA,
           silent = TRUE)
  return(p)
}

# 4.8.6 Cluster Consensus Profile (Comparison Plot)
plot_cluster_consensus_profile = function(data_norm, dfcluster, cluster_groups, metadata, title, cluster_palette) {
  # cluster_groups should be a list, e.g., list(G1 = c(4,5), G2 = c(2,3))
  
  # Average Z-scores for each protein
  cal_z_score <- function(x) { (x - mean(x)) / sd(x) }
  data_z <- t(apply(data_norm, 1, cal_z_score))
  
  plot_list <- list()
  for (name in names(cluster_groups)) {
    cl_target <- cluster_groups[[name]]
    genes <- dfcluster$ID[dfcluster$Cluster %in% cl_target]
    genes <- intersect(genes, rownames(data_z))
    
    if (length(genes) > 0) {
      cl_data_df <- data_z[genes, , drop=FALSE] %>%
        data.frame(check.names = FALSE) %>%
        tibble::rownames_to_column("Gene") %>%
        pivot_longer(-Gene, names_to = "SampleID", values_to = "Zscore")
      
      # Fix SampleID matching for join
      cl_data_df$SampleID_Match <- make.names(cl_data_df$SampleID)
      metadata_tmp <- metadata
      metadata_tmp$SampleID_Match <- make.names(metadata_tmp$SampleID)
      
      cl_data_df <- cl_data_df %>%
        left_join(metadata_tmp[, c("SampleID_Match", "Group")], by = "SampleID_Match")
      
      cl_data_df$GroupName <- name
      plot_list[[name]] <- cl_data_df
    }
  }
  
  full_plot_data <- bind_rows(plot_list)
  # Fix factor levels for time course
  target_order <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall", "H")
  full_plot_data$Group <- factor(full_plot_data$Group, levels = target_order)
  
  # Define labels based on cluster IDs
  full_plot_data$Label <- sapply(full_plot_data$GroupName, function(x) {
    paste0("Cluster ", paste(cluster_groups[[x]], collapse="+"))
  })
  
  # Plot
  p <- ggplot(full_plot_data, aes(x = Group, y = Zscore, color = Label, group = Label)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, alpha = 0.7) +
    stat_summary(fun = mean, geom = "line", size = 1.2) +
    stat_summary(fun = mean, geom = "point", size = 3) +
    scale_color_manual(values = c("Cluster 4+5" = "#E41A1C", "Cluster 2+3" = "#377EB8")) + # Hardcoded to match common pairs or adjust
    theme_classic() +
    labs(title = title,
         subtitle = "Consensus average (Mean +/- SE)",
         y = "Normalized Expression (Z-score)",
         x = "") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
    
  return(p)
}

# 4.9 Figure 2C: Lauer-Style Pathway Heatmap (User Requested)
plot_lauer_pathway_heatmap = function(tsv_path, pathway_name, comparisons, output_dir) {
  message(">>> Generating Lauer-style heatmap for: ", pathway_name)
  
  # 1. Load TSV and extract gene symbols
  pathway_data <- read.delim(tsv_path, header = TRUE, check.names = FALSE)
  # Extract symbol from "UniProt:ID Symbol" format
  pathway_genes <- sapply(strsplit(pathway_data$MoleculeName, " "), function(x) {
    if (length(x) > 1) return(x[2]) else return(NA)
  })
  pathway_genes <- na.omit(unique(pathway_genes))
  
  # 2. Extract logFC for all specified comparisons
  logFC_matrix <- matrix(NA, nrow = length(pathway_genes), ncol = length(comparisons))
  rownames(logFC_matrix) <- pathway_genes
  colnames(logFC_matrix) <- comparisons
  
  # Map symbols to our data
  all_symbols <- feature_filtered$PG.Genes
  all_symbols <- sapply(strsplit(all_symbols, ";"), `[`, 1)
  
  for (comp in comparisons) {
    res <- topTable(fit2, coef = comp, number = Inf)
    res$symbol <- all_symbols[match(rownames(res), rownames(feature_filtered))]
    
    # Match by symbol
    matches <- match(pathway_genes, res$symbol)
    logFC_matrix[, comp] <- res$logFC[matches]
  }
  
  # 3. Filter for proteins with |logFC| >= 0.3 in at least one comparison
  keep <- rowSums(abs(logFC_matrix) >= 0.3, na.rm = TRUE) > 0
  plot_data <- logFC_matrix[keep, , drop = FALSE]
  plot_data <- na.omit(plot_data)
  
  if (nrow(plot_data) < 2) {
    message("    Too few proteins found for ", pathway_name)
    return(NULL)
  }
  
  # 4. Plotting (Exact Lauer Parameters)
  # Transpose so rows are comparisons, columns are genes
  plot_data_t <- t(plot_data)
  
  # Lauer's custom breaks for RdBu (0 in the middle)
  min_val <- min(plot_data_t, na.rm = TRUE)
  max_val <- max(plot_data_t, na.rm = TRUE)
  
  breaks <- c(seq(min_val, 0, length.out = ceiling(50/2) + 1),
              seq(max_val/50, max_val, length.out = floor(50/2)))
  
  p <- pheatmap(plot_data_t, 
           scale = "none", 
           cellheight = 20,
           color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(50),
           breaks = breaks,
           fontsize_row = 10, 
           fontsize_col = 8,
           clustering_distance_cols = "euclidean",
           treeheight_col = 0, # Hide column dendrogram
           cluster_cols = TRUE,
           cluster_rows = FALSE,
           border_color = "black",
           main = pathway_name,
           silent = TRUE)
  return(p)
}

# 4.10 Figure 2G: Cluster-Highlighted Volcano Plot
plot_cluster_highlight_volcano = function(res, dfcluster, target_clusters, title, global_palette) {
  # Merge cluster info with DEA results
  res_plot <- res %>%
    data.frame() %>%
    tibble::rownames_to_column("ProteinID")
  
  # Map protein ID to UniqueSymbol
  res_plot$ID <- feature_filtered$UniqueSymbol[match(res_plot$ProteinID, rownames(feature_filtered))]
  
  # Join by ID
  res_plot <- res_plot %>% left_join(dfcluster, by = "ID")
  
  # Create a grouping for coloring
  res_plot$Highlight <- "Others"
  for (cl in target_clusters) {
    res_plot$Highlight[res_plot$Cluster == cl] <- paste0("Cluster_", cl)
  }
  
  # Ensure "Others" is at the bottom
  cl_levels <- c(paste0("Cluster_", target_clusters), "Others")
  res_plot$Highlight <- factor(res_plot$Highlight, levels = rev(cl_levels))
  
  # Use colors from global palette
  my_colors <- c("Others" = "grey85")
  for (cl in target_clusters) {
    my_colors[paste0("Cluster_", cl)] <- global_palette[paste0("Cluster_", cl)]
  }

  p <- ggplot(res_plot, aes(x = logFC, y = -log10(P.Value), color = Highlight)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = my_colors) +
    theme_classic() +
    labs(title = title,
         x = "log2(Fold Change)",
         y = "-log10(P-value)") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    theme(legend.position = "right",
          aspect.ratio = 1)
  
  return(p)
}

# ------------------------------------------------------------------------------
# 5. Generate Figures
# ------------------------------------------------------------------------------

message("\n==========================================================================")
message(">>> STARTING FIGURE GENERATION SECTION")
message("==========================================================================")

# Available Comparisons:
# "D2_F_vs_D7_F", "D7_F_vs_D14_F", "D7_C_vs_D14_C", "D7_C_vs_D7_F", 
# "D14_C_vs_D14_F", "D14_F_vs_D7_F", "D14_F_vs_D14_C", "D14_F_vs_D2_F", "D14_Wall_vs_D14_Collagen"

# SELECT YOUR COMPARISON HERE:
target_comp <- "D2_F_vs_D7_F" 

message(">>> Preparing plots for: ", target_comp)
res <- topTable(fit2, coef = target_comp, number = Inf)
ref_terms <- NULL # go_reference removed by user request
target_groups <- group_map[[target_comp]]

# p1: Volcano Plot
p1 <- plot_lauer_volcano(res, paste("Figure 2A:", target_comp))
if(!is.null(p1)) print(p1) # Use print() for ggplot in scripts

# p2: GO Enrichment (Upregulated)
p2_res <- plot_lauer_enrichment(res, paste("Figure 2B:", target_comp), direction = "up", ref_terms = ref_terms)
if(!is.null(p2_res)) {
  p2 <- p2_res$go
  print(p2)
}

# p3: GO Enrichment (Downregulated)
p3_res <- plot_lauer_enrichment(res, paste("Figure 2B:", target_comp), direction = "down", ref_terms = ref_terms)
if(!is.null(p3_res)) {
  p3 <- p3_res$go
  print(p3)
}

# p9: Multi-Ontology GO (BP, CC, MF) for two groups
message(">>> Generating Multi-Ontology GO (BP, CC, MF)...")
group_names <- group_map[[target_comp]] # e.g. D7_Fibrin, D2_Fibrin
# group_map[[target_comp]][1] is group 1 (represented by logFC > 0.3)
# group_map[[target_comp]][2] is group 2 (represented by logFC < -0.3)

sig_g1 <- res[res$P.Value < 0.05 & res$logFC > 0.3, ]
sig_g2 <- res[res$P.Value < 0.05 & res$logFC < -0.3, ]

if (nrow(sig_g1) > 5 && nrow(sig_g2) > 5) {
  g1_genes <- feature_filtered[rownames(sig_g1), "PG.Genes"]
  g1_genes <- sapply(strsplit(g1_genes, ";"), `[`, 1)
  g1_genes <- stringr::str_to_title(g1_genes)
  
  g2_genes <- feature_filtered[rownames(sig_g2), "PG.Genes"]
  g2_genes <- sapply(strsplit(g2_genes, ";"), `[`, 1)
  g2_genes <- stringr::str_to_title(g2_genes)
  
  df_compare <- rbind(
    data.frame(gene = g1_genes, group = group_names[1]),
    data.frame(gene = g2_genes, group = group_names[2])
  )
  
  gene_id_map <- bitr(df_compare$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  data_merged <- merge(gene_id_map, df_compare, by.x='SYMBOL', by.y='gene')
  
  # ont = "ALL" to get BP, CC, and MF
  p9_res <- compareCluster(ENTREZID~group, data=data_merged, fun="enrichGO", 
                          OrgDb="org.Mm.eg.db", ont = "ALL", pvalueCutoff = 0.05)
  
  if (!is.null(p9_res) && nrow(as.data.frame(p9_res)) > 0) {
    # Note: simplify() does not support ont="ALL". We plot the top results per category.
    p9 <- dotplot(p9_res, showCategory=8, title = paste("GO (BP/CC/MF):", target_comp)) +
           facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +
           theme(strip.text.y = element_text(angle = 0, face = "bold"),
                 axis.text.x = element_text(angle = 45, hjust = 1))
    print(p9)
  }
}

# p4 and p5 removed by user request

# ------------------------------------------------------------------------------
# 6. Global Plots (Multi-Volcano & Mfuzz)
# ------------------------------------------------------------------------------

# p6: Multi-Volcano Plot
message(">>> Generating Multi-Volcano Plot...")
selected_comps <- c("D2_F_vs_D7_F", "D7_F_vs_D14_F", "D7_C_vs_D14_C", 
                    "D7_C_vs_D7_F", "D14_C_vs_D14_F", "D14_F_vs_D7_F", 
                    "D14_F_vs_D14_C", "D14_F_vs_D2_F", "D14_Wall_vs_D14_Collagen")
multi_df <- data.frame()
for (comp in selected_comps) {
  res_comp <- topTable(fit2, coef = comp, number = Inf)
  symbols <- feature_filtered[rownames(res_comp), "PG.Genes"]
  symbols <- sapply(strsplit(symbols, ";"), `[`, 1)
  multi_df <- rbind(multi_df, data.frame(Name = symbols, FC = 2^res_comp$logFC, PValue = res_comp$P.Value, Cluster = comp))
}

p6 <- mutiVolcano(df = multi_df, P = 0.05, FC = 1.5, pointSize = 0.5, 
                 tileColor = RColorBrewer::brewer.pal(length(unique(multi_df$Cluster)), "Set3"))
if(!is.null(p6)) print(p6)

# p7: Mfuzz Plot
message(">>> Generating Mfuzz Heatmap...")
mfuzz_groups <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall")
mfuzz_data <- matrix(NA, nrow = nrow(expr_norm), ncol = length(mfuzz_groups))
colnames(mfuzz_data) <- mfuzz_groups
rownames(mfuzz_data) <- rownames(expr_norm)
for (g in mfuzz_groups) {
  samples <- metadata$SampleID[metadata$Group == g]
  samples <- intersect(samples, colnames(expr_norm))
  if (length(samples) > 0) mfuzz_data[, g] <- rowMeans(expr_norm[, samples, drop = FALSE], na.rm = TRUE)
}
sig_proteins <- unique(unlist(lapply(selected_comps, function(comp) {
  res_sig <- topTable(fit2, coef = comp, number = Inf)
  rownames(res_sig)[res_sig$P.Value < 0.05 & abs(res_sig$logFC) > 0.5]
})))
mfuzz_input <- mfuzz_data[intersect(sig_proteins, rownames(mfuzz_data)), ]
mfuzz_input <- na.omit(mfuzz_input)

# Use the pre-standardized UniqueSymbol for clustering IDs
rownames(mfuzz_input) <- feature_filtered[rownames(mfuzz_input), "UniqueSymbol"]

p7_res <- mfuzzHeatmap(data = mfuzz_input, clusterNum = 6)
p7 <- p7_res$p
if(!is.null(p7)) print(p7)

# p8: Cluster Enrichment (List for each cluster)
p8_plots <- list()
p8_results <- list() # Store full data frames for shared GO analysis
for (i in unique(p7_res$dfcluster$Cluster)) {
  cluster_name <- paste0("Cluster_", i)
  res_enrich <- plot_cluster_enrichment(genes = p7_res$dfcluster$ID[p7_res$dfcluster$Cluster == i], 
                                       cluster_name = cluster_name, output_dir = output_dir)
  if(!is.null(res_enrich)) {
    p8_plots[[cluster_name]] <- res_enrich$go
    p8_results[[cluster_name]] <- res_enrich$results # Save full data frames
    print(res_enrich$go) # Show each cluster enrichment
  }
}

# p10: Global Multi-Comparison GO (All Transitions)
message(">>> Generating Global Multi-Comparison GO...")
global_comps <- c("D2_F_vs_D7_F", "D7_F_vs_D14_F", "D7_C_vs_D14_C")
global_df <- data.frame()

for (comp in global_comps) {
  res_comp <- topTable(fit2, coef = comp, number = Inf)
  sig_comp <- res_comp[res_comp$P.Value < 0.05 & abs(res_comp$logFC) > 0.3, ]
  if (nrow(sig_comp) > 5) {
    comp_genes <- feature_filtered[rownames(sig_comp), "PG.Genes"]
    comp_genes <- sapply(strsplit(comp_genes, ";"), `[`, 1)
    comp_genes <- stringr::str_to_title(comp_genes) # Fix: Ensure title case for Mouse symbols
    global_df <- rbind(global_df, data.frame(gene = comp_genes, group = comp))
  }
}

if (nrow(global_df) > 0) {
  global_id_map <- bitr(global_df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  global_merged <- merge(global_id_map, global_df, by.x='SYMBOL', by.y='gene')
  
  # Ensure the comparison order is logical (D2->D7->D14)
  global_merged$group <- factor(global_merged$group, levels = global_comps)
  
  p10_res <- compareCluster(ENTREZID~group, data=global_merged, fun="enrichGO", 
                           OrgDb="org.Mm.eg.db", ont = "BP", pvalueCutoff = 0.05)
  if (!is.null(p10_res) && nrow(as.data.frame(p10_res)) > 0) {
    p10_sim <- clusterProfiler::simplify(p10_res, cutoff=0.7, by="p.adjust", select_fun=min)
    p10 <- dotplot(p10_sim, showCategory=5, title = "Global Multi-Comparison GO")
    print(p10)
  }
}

# p11: Condition-wise GO (Each condition markers vs Others)
message(">>> Generating Condition-wise GO (One-vs-Rest)...")

# Define marker genes for each group (Top genes specifically high in that group)
condition_markers <- data.frame()
for (target_g in conditions) {
  # Create a contrast: target_g - (average of all others)
  others <- setdiff(conditions, target_g)
  contrast_str <- paste0(target_g, " - (", paste(others, collapse=" + "), ")/", length(others))
  
  tmp_contrast <- makeContrasts(contrasts = contrast_str, levels = design)
  tmp_fit <- contrasts.fit(fit, tmp_contrast)
  tmp_fit <- eBayes(tmp_fit, trend = TRUE, robust = TRUE)
  
  res_marker <- topTable(tmp_fit, coef = 1, number = Inf)
  # Filter for markers: P < 0.05 and Log2FC > 0.5
  sig_marker <- res_marker[res_marker$P.Value < 0.05 & res_marker$logFC > 0.5, ]
  
  if (nrow(sig_marker) > 5) {
    m_genes <- feature_filtered[rownames(sig_marker), "PG.Genes"]
    m_genes <- sapply(strsplit(m_genes, ";"), `[`, 1)
    m_genes <- stringr::str_to_title(m_genes)
    condition_markers <- rbind(condition_markers, data.frame(gene = m_genes, group = target_g))
  }
}

if (nrow(condition_markers) > 0) {
  m_id_map <- bitr(condition_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  m_merged <- merge(m_id_map, condition_markers, by.x='SYMBOL', by.y='gene')
  
  # Ensure the group order matches the global sample order
  m_merged$group <- factor(m_merged$group, levels = final_levels)
  
  # ont = "ALL" to show BP, CC, and MF simultaneously for ALL samples
  p11_res <- compareCluster(ENTREZID~group, data=m_merged, fun="enrichGO", 
                           OrgDb="org.Mm.eg.db", ont = "ALL", pvalueCutoff = 0.05)
                           
  if (!is.null(p11_res) && nrow(as.data.frame(p11_res)) > 0) {
    # Match the faceted style found in p9
    p11 <- dotplot(p11_res, showCategory=6, title = "GO Enrichment per Sample") +
           facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +
           theme(strip.text.y = element_text(angle = 0, face = "bold"),
                 axis.text.x = element_text(angle = 45, hjust = 1))
    print(p11)
  }
}

message(">>> Objects p1 through p11 are now in your environment.")

# p12: Shared GO UpSet Plot (Comparing GO terms between clusters)
message(">>> Generating Shared GO (UpSet Plot)...")
if (length(p8_results) > 0) {
  # Prepare queries for coloring set bars using global palette
  coloring_queries <- lapply(1:length(p8_results), function(i) {
    cl_name <- names(p8_results)[i]
    list(query = intersects, params = list(cl_name), 
         color = cluster_palette[cl_name], active = TRUE)
  })

  # Create list of GO IDs for UpSet
  p8_ids <- lapply(p8_results, function(x) x$ID)
  p12 <- upset(fromList(p8_ids), 
               nsets = length(p8_results),
               nintersects = 16,
               order.by = "freq", 
               main.bar.color = "#5F6164", # Grey as requested
               queries = coloring_queries,
               text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.3))
  # UpSetR plots directly; we assign to p12 to keep convention
  
  # Export "Shared GO" details directly to Console (Studio)
  message("\n>>> SHARED GO TERM ANALYSIS:")
  # Find terms present in at least 2 clusters
  all_ids <- unlist(p8_ids)
  term_counts <- table(all_ids)
  shared_ids <- names(term_counts[term_counts > 1])
  
  if (length(shared_ids) > 0) {
    message("   Found ", length(shared_ids), " GO terms shared by multiple clusters.")
  }
}

# ------------------------------------------------------------------------------
# 7. Functional Core & Temporal Regulation Analysis
# ------------------------------------------------------------------------------

message("\n>>> Identifying Core Processes shared across 4-5 and 2-3...")

analyze_shared_processes <- function(c_pair, p8_results, mf_data) {
  c1 <- paste0("Cluster_", c_pair[1])
  c2 <- paste0("Cluster_", c_pair[2])
  
  if (is.null(p8_results[[c1]]) || is.null(p8_results[[c2]])) {
    message("    One or both clusters (", c1, ", ", c2, ") have no GO results.")
    return(NULL)
  }
  
  # Check if results are data frames (User might need to re-run p8 loop)
  if (!is.data.frame(p8_results[[c1]]) || !is.data.frame(p8_results[[c2]])) {
    message("    [ERROR]: p8_results must be data frames. Please re-run the Section 6 (Mfuzz/p8) loop to refresh data.")
    return(NULL)
  }
  
  # Intersect GO IDs
  shared_ids <- intersect(p8_results[[c1]]$ID, p8_results[[c2]]$ID)
  
  if (length(shared_ids) == 0) {
    message("    No shared GO terms found between ", c1, " and ", c2)
    return(NULL)
  }
  
  # Extract descriptions and genes
  shared_info <- p8_results[[c1]][p8_results[[c1]]$ID %in% shared_ids, c("ID", "Description", "geneID")]
  message("    Shared processes (", c1, " & ", c2, "): ", nrow(shared_info))
  print(head(shared_info[, 1:2], 10))
  
  # Temporal Regulation of underlying proteins
  # Collect all genes in these shared processes
  all_genes <- unique(unlist(strsplit(shared_info$geneID, "/")))
  # Map back to symbols (they are in symbols)
  all_genes <- intersect(all_genes, rownames(mf_data))
  
  if (length(all_genes) > 2) {
    # Plot average temporal profile
    plot_df <- mf_data[all_genes, ] %>%
      data.frame() %>%
      tibble::rownames_to_column("Gene") %>%
      pivot_longer(-Gene, names_to = "Sample", values_to = "Value") %>%
      mutate(Sample = factor(Sample, levels = colnames(mf_data)))
    
    p <- ggplot(plot_df, aes(x=Sample, y=Value, group=Gene)) +
      geom_line(alpha=0.3, color="grey") +
      stat_summary(aes(group=1), fun=mean, geom="line", size=1.5, color="red") +
      theme_classic() +
      labs(title = paste("Temporal Regulation: Core Processes (Cl", c_pair[1], "-", c_pair[2], ")"),
           subtitle = paste("Genes:", length(all_genes)),
           y = "Z-score (Expression)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(list(terms = shared_info, plot = p, genes = all_genes))
  }
  return(NULL)
}

# Run for 4-5
core_45 <- analyze_shared_processes(c(4, 5), p8_results, mfuzz_input)
if (!is.null(core_45)) print(core_45$plot)

# Run for 2-3
core_23 <- analyze_shared_processes(c(2, 3), p8_results, mfuzz_input)
if (!is.null(core_23)) print(core_23$plot)

# ------------------------------------------------------------------------------
# 8. Cluster-Specific Volcano Plots (Direct Cluster Comparison)
# ------------------------------------------------------------------------------
message("\n>>> Generating Cluster-Highlighted Volcano Plots (D14 vs D7)...")

# We use D14_F_vs_D7_F comparison
res_d14vsd7 <- topTable(fit2, coef = "D14_F_vs_D7_F", number = Inf)

# Volcano for Clusters 4 and 5
p13 <- plot_cluster_highlight_volcano(res_d14vsd7, p7_res$dfcluster, c(4, 5), 
                                     "Comparison of Cluster 4 & 5 (D14 vs D7)", cluster_palette)
if (!is.null(p13)) print(p13)

# Volcano for Clusters 2 and 3
p14 <- plot_cluster_highlight_volcano(res_d14vsd7, p7_res$dfcluster, c(2, 3), 
                                     "Comparison of Cluster 2 & 3 (D14 vs D7)", cluster_palette)
if (!is.null(p14)) print(p14)

# Helper function to export proteins and plot top20 heatmap
export_and_plot_top = function(res, dfcluster, target_cl, comparison_name, mf_input, metadata, prefix) {
  # Match proteins
  res_cl <- res %>%
    data.frame() %>%
    tibble::rownames_to_column("ProteinID")
  res_cl$ID <- feature_filtered$UniqueSymbol[match(res_cl$ProteinID, rownames(feature_filtered))]
  res_cl <- res_cl %>% filter(ID %in% dfcluster$ID[dfcluster$Cluster %in% target_cl])
  
  # Sig Up/Down
  up_genes <- res_cl$ID[res_cl$logFC > 0.5 & res_cl$P.Value < 0.05]
  down_genes <- res_cl$ID[res_cl$logFC < -0.5 & res_cl$P.Value < 0.05]
  
  message("\n>>> EXPORT: Proteins for Cluster ", paste(target_cl, collapse="+"), " in ", comparison_name)
  message("    Up-regulated: ", paste(up_genes, collapse=", "))
  message("    Down-regulated: ", paste(down_genes, collapse=", "))
  
  # Top 20 for Heatmap (by absolute LogFC among Significant)
  top20_res <- res_cl %>%
    filter(P.Value < 0.05 & abs(logFC) > 0.5) %>%
    arrange(desc(abs(logFC))) %>%
    head(20)
  top20_genes <- top20_res$ID
  
  message("    Top 20 for Heatmap: ", paste(top20_genes, collapse=", "))
  
  if (length(top20_genes) >= 1) {
    # Prepare standardized expression matrix for all samples
    all_expr_norm <- expr_norm
    # Use pre-standardized UniqueSymbol for row names
    rownames(all_expr_norm) <- feature_filtered[rownames(all_expr_norm), "UniqueSymbol"]
    
    # Generate Heatmap (Showing all individual samples + Cluster annotations)
    p_heat <- plot_custom_heatmap(data = all_expr_norm[top20_genes, , drop=FALSE], 
                                 metadata = metadata, 
                                 dfcluster = dfcluster,
                                 pathway_name = paste0("Top 20 Proteins (Cl ", paste(target_cl, collapse="+"), ")"),
                                 cluster_palette = cluster_palette)
    return(p_heat)
  }
  return(NULL)
}

# p15: Top 20 Heatmap for Cluster 4 & 5 (Individual)
p15_ind <- export_and_plot_top(res_d14vsd7, p7_res$dfcluster, c(4, 5), "D14 vs D7", mfuzz_input, metadata, "p15")

# p15_avg: Top 20 Heatmap for Cluster 4 & 5 (Group Average)
p15_avg <- export_and_plot_top(res_d14vsd7, p7_res$dfcluster, c(4, 5), "D14 vs D7", mfuzz_input, metadata, "p15_avg")
# Note: Since export_and_plot_top calls plot_custom_heatmap internally, we need to adjust it to pass 'average'
# Let's slightly modify the call for p15_avg for demo or just provide it as p15_avg
# Implementation below ensures p15_avg uses the 'average' mode.
export_and_plot_top_avg = function(res, dfcluster, target_cl, comparison_name, data_norm, feature_info, metadata, cluster_palette) {
  res_cl <- res %>% data.frame() %>% tibble::rownames_to_column("ProteinID")
  res_cl$ID <- feature_info$UniqueSymbol[match(res_cl$ProteinID, rownames(feature_info))]
  res_cl <- res_cl %>% filter(ID %in% dfcluster$ID[dfcluster$Cluster %in% target_cl])
  top20_genes <- res_cl %>% filter(P.Value < 0.05 & abs(logFC) > 0.5) %>% arrange(desc(abs(logFC))) %>% head(20) %>% pull(ID)
  
  if (length(top20_genes) >= 1) {
    all_expr_norm <- data_norm
    rownames(all_expr_norm) <- feature_info[rownames(all_expr_norm), "UniqueSymbol"]
    p_heat <- plot_custom_heatmap(data = all_expr_norm[top20_genes, , drop=FALSE], 
                                 metadata = metadata, dfcluster = dfcluster,
                                 pathway_name = paste0("Top 20 Proteins (Cl ", paste(target_cl, collapse="+"), ") - Mean Profile"),
                                 cluster_palette = cluster_palette, mode = "average")
    return(p_heat)
  }
}

p15 <- export_and_plot_top_avg(res_d14vsd7, p7_res$dfcluster, c(4, 5), "D14 vs D7", expr_norm, feature_filtered, metadata, cluster_palette)
p16 <- export_and_plot_top_avg(res_d14vsd7, p7_res$dfcluster, c(2, 3), "D14 vs D7", expr_norm, feature_filtered, metadata, cluster_palette)

if (!is.null(p15)) { grid::grid.newpage(); grid::grid.draw(p15$gtable) }
if (!is.null(p16)) { grid::grid.newpage(); grid::grid.draw(p16$gtable) }

# p17: Cluster Consensus Profile Plot
message(">>> Generating Consensus Profile (p17)...")
p17 <- plot_cluster_consensus_profile(data_norm = expr_norm, 
                                     dfcluster = p7_res$dfcluster, 
                                     cluster_groups = list(G1 = c(4,5), G2 = c(2,3)), 
                                     metadata = metadata, 
                                     title = "Thrombus Evolution: Cluster Pairwise Consensus",
                                     cluster_palette = cluster_palette)
if (!is.null(p17)) print(p17)

message(">>> Objects p1 through p17 are now in your environment.")
message(">>> You can plot them individually in RStudio (e.g., plot(p11)).")
