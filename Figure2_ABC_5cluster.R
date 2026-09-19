# Figure2_ABC.R
# Mimicking Lauer et al. 2024 Figure 2.ABC style for Thrombus Proteomics
# Using preprocessing logic from thrombus_figure_generator.R

# ------------------------------------------------------------------------------
# 1. Setup and Configuration
# ------------------------------------------------------------------------------

# Create output directory
output_dir <- "lauer_style_comparison/output20260205"
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
  library(UpSetR)
  library(VennDiagram) # For Matrisome Intersections
  library(pathlinkR)   # For Network and signaling analysis
  library(tidygraph)   # For network manipulation
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

message("    --- Data Filtering Summary ---")
message("    Total Proteins in Header: ", nrow(expr_data))
message("    Proteins Passing 50% Filter: ", nrow(expr_filtered))

# 2.2 Standardize Gene Symbols for consistent mapping across all plots
feature_filtered$Symbol <- as.character(feature_filtered$PG.Genes)
missing_gene <- is.na(feature_filtered$Symbol) | feature_filtered$Symbol == ""

if (any(missing_gene)) {
  descriptions <- feature_filtered$PG.ProteinDescriptions[missing_gene]
  extracted_genes <- stringr::str_extract(descriptions, "GN=[^ ]+")
  extracted_genes <- gsub("GN=", "", extracted_genes)
  feature_filtered$Symbol[missing_gene] <- ifelse(!is.na(extracted_genes), extracted_genes, 
                                                sapply(strsplit(feature_filtered$PG.ProteinGroups[missing_gene], ";"), `[`, 1))
}
feature_filtered$Symbol <- stringr::str_to_title(feature_filtered$Symbol)
feature_filtered$UniqueSymbol <- make.unique(feature_filtered$Symbol)

# 2.3 Log2 Transformation & Imputation
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

# Define Global Color Palettes & Themes (Imported from thrombus_proteomics_analysis.R)
theme_elegant <- function() {
  theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      text = element_text(family = "sans"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

time_colors <- c("D2" = "#D5F5E3", "D7" = "#52BE80", "D14" = "#1E8449")
region_colors <- c("Fibrin" = "#F1948A", "Collagen" = "#5DADE2", "Wall" = "#8E44AD")
group_colors <- c(
  "D2_Fibrin"     = "#FADBD8", 
  "D7_Fibrin"     = "#EC7063", 
  "D14_Fibrin"    = "#B03A2E", 
  "D7_Collagen"   = "#AED6F1", 
  "D14_Collagen"  = "#2E86C1", 
  "D14_Wall"      = "#8E44AD"
)

cluster_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
names(cluster_palette) <- paste0("Cluster_", 1:6)

message(">>> Performing DEA...")

# Define Sample Order as requested: D2-F, D7-F, D7-C, D14-F, D14-C, D14-Wall,"D14_Fibrin"
target_order <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Collagen", "D14_Wall")
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
  
  # Regional Heterogeneity
  D7_C_vs_D7_F = D7_Collagen - D7_Fibrin,
  D14_C_vs_D14_F = D14_Collagen - D14_Fibrin,
  
  # Combined Maturation & Transition (User Requested)
  D14_C_vs_D7_F = D14_Collagen - D7_Fibrin,
  D14_C_vs_D7_C = D14_Collagen - D7_Collagen,
  
  # Resolution Failure Context
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
    theme_elegant() +
    theme(text = element_text(size = 15, color = "#000000"),
          legend.title = element_blank(),
          axis.line = element_line(colour = "#000000", size = 1), 
          legend.position = "none",
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
    theme_elegant() +
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
           color = colorRampPalette(c("#2E86B9", "white", "#B03A2E"))(50),
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
    theme_elegant()+
    ggplot2::scale_y_continuous(n.breaks = 5) + 
    ggplot2::theme(
                   legend.position = "right", 
                   legend.title = ggplot2::element_blank(), 
                   legend.background = ggplot2::element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.line.x = element_blank()
                   ) + 
    ggplot2::xlab("Comparisons") + ggplot2::ylab("log2FC") + 
    guides(color = guide_legend(override.aes = list(size = 3)))
    
    return(p)
}

# 4.7 Figure 2E: Mfuzz Heatmap (User Requested)
mfuzzHeatmap = function(data,
                        clusterNum=4
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
      stat_summary(aes(group=1), fun=mean, geom="line", size=1.0, color="red") +
      scale_color_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "Spectral")), # Publication standard palette
                            breaks=seq(0,1,0.2),
                            limits=c(0, 1)
      )+
      theme_elegant()+ # Applied elegant theme
      labs(y=paste0("cluster ",clusterName,"\n","n=",myN),
           x="")+
      theme(legend.position = "none", # Hide individual legends to clean up
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
      scale_fill_gradientn(colors = colorRampPalette(c("#2E86B9", "white", "#B03A2E"))(100),
                           limits = c(-2.5, 2.5), oob = scales::squish) +
      theme_elegant()+
      labs(y="",x="",fill="z-score")+
      theme(legend.position = "none",
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
                               labels = c("BP", "CC", "MF"))
    
    # Create the plot mimicking the user's example
    p_go <- ggplot(ego_top, aes(x = reorder(Description, -log10(pvalue)), y = -log10(pvalue), fill = ONTOLOGY)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free") +
      scale_fill_manual(values = c("BP" = "#91CF60", 
                                   "CC" = "#FC8D59", 
                                   "MF" = "#4575B4")) +
      theme_elegant() +
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
  # Strictly match rownames to avoid NAs in annotation
  matching_indices <- match(rownames(data_z), dfcluster$ID)
  # Filter out proteins not in cluster df (though they should be)
  valid_row_indices <- which(!is.na(matching_indices))
  
  if (length(valid_row_indices) < nrow(data_z)) {
    data_z <- data_z[valid_row_indices, , drop = FALSE]
    matching_indices <- matching_indices[valid_row_indices]
  }
  
  row_annot <- dfcluster[matching_indices, "Cluster", drop=FALSE]
  row_annot$Cluster <- factor(paste0("Cluster_", row_annot$Cluster))
  rownames(row_annot) <- rownames(data_z)
  
  # 4. Colors
  ann_colors <- list(
    Group = group_colors,
    Cluster = cluster_palette
  )

  # 5. Plotting
  p <- pheatmap(data_z, 
           scale = "none",
           color = colorRampPalette(c("#2E86B9", "white", "#B03A2E"))(100),
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
    theme_elegant() +
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
  
# 4.9.5 Lauer-Style Expression Heatmap (Group Average)
plot_lauer_style_expression_heatmap = function(data, metadata, title, order_levels) {
  # 1. Group Averaging
  data_df <- as.data.frame(t(data))
  data_df$Group <- metadata$Group[match(make.names(rownames(data_df)), make.names(metadata$SampleID))]
  
  avg_data <- data_df %>%
    filter(!is.na(Group)) %>%
    group_by(Group) %>%
    summarise(across(everything(), \(x) mean(x, na.rm=TRUE)), .groups = "drop") %>%
    mutate(Group = factor(Group, levels = order_levels)) %>%
    filter(!is.na(Group)) %>%  # Safety: Remove any groups not mentioned in order_levels
    arrange(Group) %>%
    tibble::column_to_rownames("Group")
  
  # 2. Z-score Scaling (on average data)
  plot_mat <- t(as.matrix(avg_data))
  plot_mat_z <- t(apply(plot_mat, 1, function(x) {
    if(sd(x, na.rm=T) == 0) return(x - mean(x, na.rm=T))
    (x - mean(x, na.rm=T)) / sd(x, na.rm=T)
  }))
  
  # 3. Final NA and Name Cleaning
  # Ensure NO NAs in names before pheatmap
  valid_rows <- which(apply(plot_mat_z, 1, function(x) !any(is.na(x))) & !is.na(rownames(plot_mat_z)))
  if (length(valid_rows) < 2) return(NULL)
  plot_mat_z <- plot_mat_z[valid_rows, , drop=FALSE]
  
  # 4. Plotting using Lauer parameters
  p <- pheatmap(plot_mat_z, 
           scale = "none", 
           color = colorRampPalette(c("#2E86B9", "white", "#B03A2E"))(100),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           border_color = "black",
           cellheight = 12,
           fontsize_row = 8,
           fontsize_col = 10,
           main = title,
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
    theme_elegant() +
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
selected_comps <- colnames(contrast_matrix)
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
mfuzz_groups <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Collagen", "D14_Wall")
mfuzz_data <- matrix(NA, nrow = nrow(expr_norm), ncol = length(mfuzz_groups))
colnames(mfuzz_data) <- mfuzz_groups
rownames(mfuzz_data) <- rownames(expr_norm)
for (g in mfuzz_groups) {
  samples <- metadata$SampleID[metadata$Group == g]
  samples <- intersect(samples, colnames(expr_norm))
  if (length(samples) > 0) mfuzz_data[, g] <- rowMeans(expr_norm[, samples, drop = FALSE], na.rm = TRUE)
}
# Identify significant proteins for clustering
# Robust selection: Check each coefficient in fit2
sig_proteins <- character(0)
for (comp in selected_comps) {
  if (comp %in% colnames(fit2)) {
    res_sig <- topTable(fit2, coef = comp, number = Inf)
    if (!is.null(res_sig) && nrow(res_sig) > 0) {
      sig_ids <- rownames(res_sig)[res_sig$P.Value < 0.05 & abs(res_sig$logFC) > 0.5]
      sig_proteins <- c(sig_proteins, sig_ids)
    }
  }
}
sig_proteins <- unique(sig_proteins)
message("    Significant proteins found for Mfuzz: ", length(sig_proteins))

mfuzz_input <- mfuzz_data[intersect(sig_proteins, rownames(mfuzz_data)), , drop = FALSE]
mfuzz_input <- na.omit(mfuzz_input)

# Use the pre-standardized UniqueSymbol for clustering IDs
rownames(mfuzz_input) <- feature_filtered[rownames(mfuzz_input), "UniqueSymbol"]

p7_res <- mfuzzHeatmap(data = mfuzz_input, clusterNum = 5)
p7 <- p7_res$p
if(!is.null(p7)) print(p7)

# p8: Cluster Enrichment (List for each cluster)
p8_plots <- list()
p8_results <- list() # Store full data frames for shared GO analysis

if (!is.null(p7_res$dfcluster) && nrow(p7_res$dfcluster) > 0) {
    for (i in sort(unique(p7_res$dfcluster$Cluster))) {
      cluster_name <- paste0("Cluster_", i)
      res_enrich <- plot_cluster_enrichment(genes = p7_res$dfcluster$ID[p7_res$dfcluster$Cluster == i], 
                                           cluster_name = cluster_name, output_dir = output_dir)
      if(!is.null(res_enrich)) {
        p8_plots[[cluster_name]] <- res_enrich$go
        p8_results[[cluster_name]] <- res_enrich$results # Save full data frames
        print(res_enrich$go) # Show each cluster enrichment
      }
    }
} else {
    message("    Note: No clusters found/generated by Mfuzz.")
}
p8
# p8b: Cluster-Condition Bubble Plot (User Requested)
# This shows GO enrichment per cluster, broken down by condition (D2-F, D7-C, etc.)
message(">>> Generating Cluster-Condition Bubble Plots (p8b)...")
p8b_plots <- list()

# Define the specific requested order and original-to-label mapping
p8b_order_map <- c(
  "D2_Fibrin"   = "D2-F",
  "D7_Fibrin"   = "D7-F",
  "D7_Collagen" = "D7-C",
  "D14_Fibrin"  = "D14-F",
  "D14_Collagen" = "D14-C",
  "D14_Wall"     = "D14-Wall"
)

# Iterate through each cluster to create an individual plot
for (i in sort(unique(p7_res$dfcluster$Cluster))) {
  cl_name <- paste0("Cluster_", i)
  cl_genes <- p7_res$dfcluster$ID[p7_res$dfcluster$Cluster == i]
  
  # For this cluster, find which genes are markers for each condition
  cl_condition_markers <- data.frame()
  
  for (target_g in intersect(names(p8b_order_map), conditions)) {
    # Marker logic: Genes in this cluster that are UP in this condition vs average of all others
    others <- setdiff(conditions, target_g)
    contrast_str <- paste0("`", target_g, "` - (", paste(paste0("`", others, "`"), collapse=" + "), ")/", length(others))
    tmp_contrast <- makeContrasts(contrasts = contrast_str, levels = design)
    tmp_fit <- contrasts.fit(fit, tmp_contrast)
    tmp_fit <- eBayes(tmp_fit, trend = TRUE, robust = TRUE)
    res_marker <- topTable(tmp_fit, coef = 1, number = Inf)
    
    # Sig in this condition (P < 0.05, LogFC > 0.3) AND belongs to this cluster
    sig_protein_ids <- rownames(res_marker)[res_marker$P.Value < 0.05 & res_marker$logFC > 0.3]
    m_genes <- feature_filtered$UniqueSymbol[match(sig_protein_ids, rownames(feature_filtered))]
    m_genes <- intersect(m_genes, cl_genes)
    
    if (length(m_genes) >= 3) {
      # Map to ENTREZID for enrichment
      m_entrez_df <- bitr(m_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
      if (nrow(m_entrez_df) > 0) {
        cl_condition_markers <- rbind(cl_condition_markers, 
                                     data.frame(ENTREZID = m_entrez_df$ENTREZID, 
                                                group = p8b_order_map[target_g]))
      }
    }
  }
  
  if (nrow(cl_condition_markers) > 0) {
    # Ensure group is a factor with requested order
    cl_condition_markers$group <- factor(cl_condition_markers$group, levels = unname(p8b_order_map))
    
    # Run compareCluster for this specific cluster's conditions
    ck_cl <- compareCluster(ENTREZID ~ group, data = cl_condition_markers, fun = "enrichGO", 
                            OrgDb = "org.Mm.eg.db", ont = "ALL", pvalueCutoff = 0.05)
    
    if (!is.null(ck_cl) && nrow(as.data.frame(ck_cl)) > 1) {
      p_cl <- dotplot(ck_cl, showCategory = 5) + 
              facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +
              theme_elegant() +
              theme(strip.text.y = element_text(angle = 0, face = "bold"),
                    axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(title = paste("GO Enrichment by Condition:", cl_name),
                   subtitle = "Restricted to cluster-specific marker proteins") +
              scale_color_gradient(low = "blue", high = "red")
      
      p8b_plots[[cl_name]] <- p_cl
      print(p_cl) # Display in RStudio
    } else {
      message("    No significant GO terms found for ", cl_name, " by condition subset.")
    }
  } else {
    message("    Too few marker genes found in ", cl_name, " for condition-wise breakdown.")
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

# p11: Condition-wise GO (Each condition markers vs Others)
message(">>> Generating Condition-wise GO (One-vs-Rest)...")

# 1. Define the specific requested order and original-to-label mapping
# This explicitly excludes H_Inflammation (High inflammation)
p11_order_map <- c(
  "D2_Fibrin"   = "D2-F",
  "D7_Fibrin"   = "D7-F",
  "D7_Collagen" = "D7-C",
  "D14_Fibrin"  = "D14-F",
  "D14_Collagen" = "D14-C",
  "D14_Wall"     = "D14-Wall"
)

# Only process groups that exist in our mapping and the data
p11_targets <- intersect(names(p11_order_map), conditions)

# 2. Define marker genes for each group (Top genes specifically high in that group)
condition_markers <- data.frame()
for (target_g in p11_targets) {
  # Create a contrast: target_g - (average of all others)
  # 'others' includes all groups in the original 'conditions' (including H_Inflammation)
  others <- setdiff(conditions, target_g)
  
  # Use backticks in contrast string to handle special characters correctly
  contrast_str <- paste0("`", target_g, "` - (", paste(paste0("`", others, "`"), collapse=" + "), ")/", length(others))
  
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

# 3. Perform Enrichment and Generate Plot
if (nrow(condition_markers) > 0) {
  # Map symbols to Entrez IDs
  m_id_map <- bitr(condition_markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  m_merged <- merge(m_id_map, condition_markers, by.x='SYMBOL', by.y='gene')
  
  # Rename the groups to the requested labels (D2-F, etc.) and set the order
  m_merged$group <- factor(p11_order_map[m_merged$group], levels = unname(p11_order_map))
  
  # ont = "ALL" to show BP, CC, and MF simultaneously
  p11_res <- compareCluster(ENTREZID~group, data=m_merged, fun="enrichGO", 
                            OrgDb="org.Mm.eg.db", ont = "ALL", pvalueCutoff = 0.05)
  
  if (!is.null(p11_res) && nrow(as.data.frame(p11_res)) > 0) {
    p11 <- dotplot(p11_res, showCategory=6, title = "GO Enrichment per Sample") +
      facet_grid(ONTOLOGY ~ ., scales = "free", space = "free") +
      theme(strip.text.y = element_text(angle = 0, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    print(p11)
  }
}


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
p12
# ------------------------------------------------------------------------------
# 7. Pathway-Filtered Functional Core Analysis (Cluster 2 & 3)
# ------------------------------------------------------------------------------

# ==============================================================================
# 7.5 Matrisome Intersection & Venn Diagrams (User Requested)
# ==============================================================================
message("\n>>> Generating Matrisome vs Cluster Intersections (Venn Diagrams)...")

# Provide/Load Matrisome Database Lists
matrisome_dir <- "C:/Users/SimonYao/Desktop/LCM_protemics_thrombus/Matrisome"
core_paths <- list.files(file.path(matrisome_dir, "Core Matrisome"), pattern="\\.tsv$|\\.txt$|\\.csv$", full.names=TRUE)
assoc_paths <- list.files(file.path(matrisome_dir, "Matrisome associated"), pattern="\\.tsv$|\\.txt$|\\.csv$", full.names=TRUE)

# Helper to normalize identifiers (UniProt/Gene) - Exact logic from thrombus_proteomics_analysis.R
normalize_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("^UNIPROT\\:[[:space:]]*|^SP\\|", "", x)
  x <- gsub("(_MOUSE|_HUMAN)$", "", x)
  x <- gsub("-[0-9]+$", "", x) # Strip isoforms
  return(x)
}

load_local_db_proteins_as_list <- function(file_paths) {
  db_list <- list()
  for (f in file_paths) {
    if (file.exists(f)) {
        name <- gsub("\\.tsv$|\\.txt$|\\.csv$", "", basename(f))
        
        # Robust reading for manually organized files (Gene \t UniProt)
        df <- tryCatch({
            read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                       fill = TRUE, quote = "", comment.char = "")
        }, error = function(e) return(NULL))
        
        if (!is.null(df) && nrow(df) > 0) {
            ids <- character(0)
            # Use columns 1 and 2 (Gene and UniProt)
            cols_to_use <- intersect(1:ncol(df), 1:2)
            for (i in cols_to_use) {
                raw_vals <- as.character(df[[i]])
                # Clean up: skip common headers if present
                clean_vals <- raw_vals[!grepl("^Gene$|^UniProt$|^Symbol$|^Protein$", raw_vals, ignore.case = TRUE)]
                normalized <- normalize_id(clean_vals)
                ids <- c(ids, normalized)
            }
            db_list[[name]] <- unique(ids)
            message("      - Loaded Category [", name, "]: ", length(unique(ids)), " identifiers")
        }
    }
  }
  return(db_list)
}

# --- Perform Matrisome Intersection with Clusters ---
message("\n>>> Generating Matrisome vs Cluster Intersections (Venn Diagrams)...")

# User's manually organized Matrisome directory
matrisome_dir <- "c:/Users/SimonYao/Desktop/LCM_protemics_thrombus/lauer_style_comparison/Matrisome"
core_paths <- list.files(file.path(matrisome_dir, "Core Matrisome"), pattern="\\.tsv$|\\.txt$|\\.csv$", full.names=TRUE)
assoc_paths <- list.files(file.path(matrisome_dir, "Matrisome associated"), pattern="\\.tsv$|\\.txt$|\\.csv$", full.names=TRUE)

matrisome_full_db <- load_local_db_proteins_as_list(c(core_paths, assoc_paths))

# Helper to map identifiers to our features robustly - RETURNS UniqueSymbol
get_genes_from_ids <- function(id_list) {
  valid_db_ids <- unique(id_list) 
  if (length(valid_db_ids) == 0) return(character(0))
  
  # Prepare searchable tokens for each protein in our dataset
  # We match against Gene Symbols and Protein Accessions
  match_mask <- rep(FALSE, nrow(feature_filtered))
  
  for (i in 1:nrow(feature_filtered)) {
    # Extract identifiers from multiple possible columns
    raw_ids <- c(
      feature_filtered$PG.Genes[i],
      feature_filtered$PG.ProteinGroups[i],
      feature_filtered$Symbol[i],
      rownames(feature_filtered)[i]
    )
    
    # Split by common delimiters: ; , | / and space
    tokens <- unlist(strsplit(as.character(raw_ids), "[; ,/|]"))
    tokens <- tokens[!is.na(tokens) & tokens != ""]
    
    # Normalize (Upper, strip sp|, etc.)
    clean_tokens <- normalize_id(tokens)
    
    # Check for match
    if (any(clean_tokens %in% valid_db_ids)) {
      match_mask[i] <- TRUE
    }
  }
  
  # Return UniqueSymbol to match Mfuzz cluster IDs
  return(feature_filtered$UniqueSymbol[match_mask])
}

if (length(matrisome_full_db) > 0) {
  total_matrisome_ids <- unique(unlist(matrisome_full_db))
  
  # Map recognized IDs to our dataset genes
  matrisome_detected_genes <- get_genes_from_ids(total_matrisome_ids)
  message("    Detected Matrisome Proteins in dataset: ", length(matrisome_detected_genes))
}

# Print detected Matrisome genes clearly for future reference
if (length(matrisome_detected_genes) > 0) {
    cat("\n--- DETECTED MATRISOME GENES (N=", length(matrisome_detected_genes), ") ---\n")
    cat(paste(sort(matrisome_detected_genes), collapse = ", "), "\n")
    cat("--------------------------------------------------\n\n")

  # We will do this per cluster using the p7_res$dfcluster object from Mfuzz
  cluster_genes <- split(p7_res$dfcluster$ID, p7_res$dfcluster$Cluster)
  
  # Function to draw advanced Venn (blue/grey minimalist style)
  draw_advanced_mat_venn_cluster <- function(venn_list, title, filename) {
    if (length(venn_list) < 2) return(NULL)
    # Refined palette: Grey and Deep Blue as requested
    venn_colors <- c("#EBEBEB", "#005596") 
    
    venn_obj <- VennDiagram::venn.diagram(
      x = venn_list,
      category.names = names(venn_list),
      filename = NULL,
      fill = venn_colors[1:length(venn_list)],
      alpha = 0.8,
      col = "black", # Black borders
      lwd = 1.5,
      cex = 1.5,
      cat.cex = 1.2,
      cat.fontface = "bold",
      main = title,
      main.fontface = "bold",
      main.cex = 1.5,
      margin = 0.05
    )
    
    pdf(filename, width = 8, height = 8)
    grid::grid.draw(venn_obj)
    
    # Text summary of intersection
    shared_prots <- intersect(venn_list[[1]], venn_list[[2]])
    if (length(shared_prots) > 0) {
      grid::grid.text(paste("Intersect:", length(shared_prots)), x = 0.5, y = 0.1, gp = grid::gpar(fontface="bold", cex=1.2))
      # Top 15 display
      display_genes <- head(shared_prots, 15)
      grid::grid.text(paste(display_genes, collapse=", "), x = 0.5, y = 0.05, gp = grid::gpar(fontsize=8))
    }
    dev.off()
    
    # Print to console for RStudio
    grid::grid.newpage()
    grid::grid.draw(venn_obj)
  }
  
  # Generate Venn, CSV and Heatmap for each cluster
  dir.create(file.path(output_dir, "Matrisome_Venns"), showWarnings = FALSE)
  for (cl_name in names(cluster_genes)) {
    c_genes <- cluster_genes[[cl_name]]
    int_genes <- intersect(c_genes, matrisome_detected_genes)
    
    if (length(int_genes) > 0) {
      # 1. Draw Venn
      v_list <- list(Cluster = c_genes, Matrisome = matrisome_detected_genes)
      names(v_list)[1] <- paste0("Cluster ", cl_name)
      draw_advanced_mat_venn_cluster(
        v_list, 
        paste("Intersection: Cluster", cl_name, "& Matrisome"), 
        file.path(output_dir, "Matrisome_Venns", paste0("Venn_Matrisome_Cluster", cl_name, ".pdf"))
      )
      
      # 2. Export CSV with Sub-category Annotation
      # Match UniqueSymbols in intersection back to their features
      export_df <- feature_filtered[feature_filtered$UniqueSymbol %in% int_genes, ]
      
      # For each gene, find its Matrisome sub-sets
      export_df$Matrisome_Categories <- sapply(export_df$UniqueSymbol, function(gs) {
        row_idx <- which(feature_filtered$UniqueSymbol == gs)[1]
        # Match using identifiers only (avoiding descriptions to maintain precision)
        ids <- normalize_id(unlist(strsplit(paste(feature_filtered$PG.Genes[row_idx], 
                                                 feature_filtered$PG.ProteinGroups[row_idx], 
                                                 feature_filtered$Symbol[row_idx], sep=";"), "[; ,/|]")))
        
        # Check against each Matrisome sub-database
        cats <- names(matrisome_full_db)[sapply(matrisome_full_db, function(db_ids) any(ids %in% db_ids))]
        return(paste(cats, collapse="; "))
      })
      
      # Reorder for better CSV readability
      export_df <- export_df %>% select(UniqueSymbol, Matrisome_Categories, everything())
      csv_path <- file.path(output_dir, "Matrisome_Venns", paste0("Cluster", cl_name, "_Matrisome_Intersection.csv"))
      write.csv(export_df, csv_path, row.names=FALSE)
      message("      [CSV Export] Saved: ", normalizePath(csv_path))
      
      # 3. Group-Level Heatmaps (User Style)
      # Prepare Annotation for Group-Level Heatmaps
      target_order_levels <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall")
      heat_anno_col <- data.frame(Group = factor(target_order_levels, levels = target_order_levels))
      rownames(heat_anno_col) <- heat_anno_col$Group
      heat_anno_colors <- list(Group = group_colors)
      
      # 3a. Top 20 Individual Proteins (Group Averages)
      top_20_int <- export_df %>%
        mutate(MeanExpr = rowMeans(expr_norm[rownames(export_df), ], na.rm=TRUE)) %>%
        arrange(desc(MeanExpr)) %>%
        head(20) %>%
        pull(UniqueSymbol)
      
      if (length(top_20_int) > 0) {
        top_ids <- rownames(feature_filtered)[feature_filtered$UniqueSymbol %in% top_20_int]
        # Calculate Group Averages
        heat_mat <- matrix(0, nrow = length(top_ids), ncol = length(target_order_levels))
        colnames(heat_mat) <- target_order_levels
        rownames(heat_mat) <- feature_filtered$UniqueSymbol[match(top_ids, rownames(feature_filtered))]
        
        for (g in colnames(heat_mat)) {
          s_ids <- metadata$SampleID[metadata$Group == g]
          if (length(s_ids) > 0) {
            heat_mat[, g] <- rowMeans(expr_norm[top_ids, s_ids, drop=FALSE], na.rm=TRUE)
          }
        }
        
        # Plot to RStudio
        ph1 <- pheatmap(heat_mat, 
                 scale = "row", 
                 cluster_cols = FALSE,
                 show_colnames = FALSE, 
                 annotation_col = heat_anno_col, 
                 annotation_colors = heat_anno_colors,
                 main = paste0("Cluster ", cl_name, " Matrisome: Top 20 Proteins"), 
                 color = colorRampPalette(c("#2E86B9", "white", "#B03A2E"))(100), 
                 border_color = NA,
                 fontsize_row = 10)
        print(ph1)
        
        # Save to PDF
        pdf(file.path(output_dir, "Matrisome_Venns", paste0("Heatmap_Cluster", cl_name, "_Matrisome_Top20_GroupAvg.pdf")), width=8, height=8)
        grid::grid.draw(ph1$gtable)
        dev.off()
      }
      
      # 3b. Categorical Group Heatmaps (Category Averages)
      # For rows, we use the 6 Matrisome categories
      cat_names <- names(matrisome_full_db)
      cat_heat_mat <- matrix(0, nrow = length(cat_names), ncol = length(target_order_levels))
      colnames(cat_heat_mat) <- target_order_levels
      rownames(cat_heat_mat) <- cat_names
      
      for (cat in cat_names) {
        # Find genes in this cluster AND in this matrisome category
        cat_db_ids <- matrisome_full_db[[cat]]
        
        # Check matching for each gene in the intersection
        cat_mask <- sapply(int_genes, function(gs) {
          row_idx <- which(feature_filtered$UniqueSymbol == gs)[1]
          ids <- normalize_id(unlist(strsplit(paste(feature_filtered$PG.Genes[row_idx], 
                                                   feature_filtered$PG.ProteinGroups[row_idx], 
                                                   feature_filtered$Symbol[row_idx], sep=";"), "[; ,/|]")))
          any(ids %in% cat_db_ids)
        })
        cat_genes <- int_genes[cat_mask]
        
        if (length(cat_genes) > 0) {
          cat_ids <- rownames(feature_filtered)[feature_filtered$UniqueSymbol %in% cat_genes]
          for (g in colnames(cat_heat_mat)) {
            s_ids <- metadata$SampleID[metadata$Group == g]
            if (length(s_ids) > 0) {
              # Average of all proteins in this category across all samples in this group
              cat_heat_mat[cat, g] <- mean(expr_norm[cat_ids, s_ids], na.rm=TRUE)
            }
          }
        }
      }
      
      # Remove categories with no data
      cat_heat_mat <- cat_heat_mat[rowSums(cat_heat_mat != 0, na.rm=TRUE) > 0, , drop=FALSE]
      
      if (nrow(cat_heat_mat) > 0) {
        # Plot to RStudio
        ph2 <- pheatmap(cat_heat_mat, 
                 scale = "row", 
                 cluster_cols = FALSE,
                 show_colnames = FALSE, 
                 annotation_col = heat_anno_col, 
                 annotation_colors = heat_anno_colors,
                 main = paste0("Cluster ", cl_name, " Matrisome: Category Trends"), 
                 color = colorRampPalette(c("#2E86B9", "white", "#B03A2E"))(100), 
                 border_color = NA,
                 fontsize_row = 10)
        print(ph2)
        
        # Save to PDF
        pdf(file.path(output_dir, "Matrisome_Venns", paste0("Heatmap_Cluster", cl_name, "_Matrisome_Categories_GroupAvg.pdf")), width=8, height=6)
        grid::grid.draw(ph2$gtable)
        dev.off()
      }
      message("    - Cluster ", cl_name, ": Exported Venn, CSV and Group-Level Heatmaps")
    }
  }
  # --- 6.1 Cluster 1 x Matrisome Highlight Volcano Plots ---
  # Intersection specifically for Cluster 1 (the main Matrisome cluster)
  c1_genes <- cluster_genes[["1"]]
  c1_mat_intersection <- intersect(c1_genes, matrisome_detected_genes)
  
  if (length(c1_mat_intersection) > 0) {
    message("\n>>> Generating Highlight Volcano Plots for Cluster 1 Matrisome intersection...")
    
    # Comparisons of interest
    volcano_contrasts <- c("D14_C_vs_D7_F", "D14_Wall_vs_D14_Collagen")
    
    for (con in volcano_contrasts) {
      if (con %in% colnames(fit2)) {
        res_comp <- topTable(fit2, coef = con, number = Inf)
        res_comp$UniqueSymbol <- feature_filtered$UniqueSymbol[match(rownames(res_comp), rownames(feature_filtered))]
        
        # Color coding
        res_comp$Color <- "Neutral"
        res_comp$Color[res_comp$P.Value < 0.05 & res_comp$logFC > 0.5] <- "Up"
        res_comp$Color[res_comp$P.Value < 0.05 & res_comp$logFC < -0.5] <- "Down"
        res_comp$Color[res_comp$UniqueSymbol %in% c1_mat_intersection] <- "Highlight"
        
        # Plotting
        v_title <- gsub("_", " ", con)
        p_vol <- ggplot(res_comp, aes(x = logFC, y = -log10(P.Value))) +
          geom_point(data = filter(res_comp, Color == "Neutral"), color = "grey80", alpha = 0.4, size = 1) +
          geom_point(data = filter(res_comp, Color == "Up"), color = "#F1948A", alpha = 0.4, size = 1) +
          geom_point(data = filter(res_comp, Color == "Down"), color = "#5DADE2", alpha = 0.4, size = 1) +
          geom_point(data = filter(res_comp, Color == "Highlight"), color = "#D4AC0D", alpha = 0.9, size = 2.5) +
          geom_text_repel(data = filter(res_comp, Color == "Highlight" & (abs(logFC) > 1 | -log10(P.Value) > 2)),
                          aes(label = UniqueSymbol), size = 3, force = 2, max.overlaps = 100) +
          theme_elegant() +
          labs(title = paste0("Highlight: Cluster 1 Matrisome (", v_title, ")"),
               x = "log2 Fold Change", y = "-log10 P-Value") +
          geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50")
        
        print(p_vol)
        ggsave(file.path(output_dir, paste0("Volcano_Highlight_C1_Matrisome_", con, ".pdf")), p_vol, width = 7, height = 6)
      }
    }
  }
} else {
    message("    Warning: No Matrisome proteins detected.")
}

# --- 6.2 pathlinkR Analysis for Maturation Intersection ---
message("\n>>> Reached pathlinkR Analysis block...")
# We use the Cluster 1 Matrisome intersection for network signaling analysis
if (length(c1_mat_intersection) > 0) {
  message("\n>>> Performing pathlinkR analysis on Cluster 1 Matrisome intersection...")
  
  # Prepare pathlinkR input: Map Mouse Symbols to Human ENSEMBL IDs (ENSG...)
  # pathlinkR 1.6.0 strictly validates for ENSG... IDs even in Sigora/PPI modules
  c1_mat_features <- feature_filtered[feature_filtered$UniqueSymbol %in% c1_mat_intersection, ]
  
  # 1. Load the internal mapping file as our ortholog bridge
  data("mappingFile", package = "pathlinkR", envir = environment())
  m_file <- as.data.frame(get("mappingFile"))
  
  # 2. Convert Mouse Symbols to Human Uppercase Symbols for mapping
  c1_mat_features$HumanSymbol <- toupper(c1_mat_features$Symbol)
  
  # 3. Join with internal mappingFile to get ENSG IDs
  mapping_bridge <- merge(c1_mat_features, m_file, by.x = "HumanSymbol", by.y = "hgncSymbol", all.x = TRUE)
  mapping_bridge <- mapping_bridge[!is.na(mapping_bridge$ensemblGeneId), ]
  mapping_bridge <- mapping_bridge[!duplicated(mapping_bridge$ensemblGeneId), ]
  
  message("    - Human ENSEMBL (ENSG) bridge success: ", nrow(mapping_bridge), " / ", nrow(c1_mat_features), " genes.")
  
    if (nrow(mapping_bridge) > 0) {
    # Create the mock result with ENSG row names (Required for pathlinkR validity)
    res_mock <- data.frame(
      log2FoldChange = rep(2.0, nrow(mapping_bridge)), 
      padj = rep(0.01, nrow(mapping_bridge)),
      hgncSymbol = as.character(mapping_bridge$Symbol), # Explicitly character for case_when
      stringsAsFactors = FALSE
    )
    rownames(res_mock) <- mapping_bridge$ensemblGeneId
    
    # --- pathlinkR Signaling Flow (Manual 6.1/6.2 Pattern) ---
    tryCatch({
      # 2. Build PPI Network (InnateDB)
      message("    - Building PPI Network (InnateDB)...")
      exNetwork <- ppiBuildNetwork(
        rnaseqResult = res_mock,
        filterInput = FALSE,
        columnFC = "log2FoldChange",
        columnP = "padj",
        order = "zero"
      )
      
      if (!is.null(exNetwork)) {
        # 3. Enrich Network for Pathways (Sigora)
        message("    - Performing Network Enrichment...")
        pLink_results <- ppiEnrichNetwork(exNetwork, analysis = "sigora")
        
        # Display Results in Console
        if (!is.null(pLink_results) && nrow(pLink_results) > 0) {
          message("    - Pathway enrichment: Found ", nrow(pLink_results), " enriched pathways.")
          print(pLink_results) 
          write.csv(pLink_results, file.path(output_dir, "pathlinkR_Enriched_Pathways_C1_Matrisome.csv"), row.names = FALSE)
        }
        
        # 4. Plot and Display in RStudio Plots Window
        # Fixed: Manual Attribute Injection for robust labeling
        # pathlinkR 1.6.0 often fails to carry over attributes during ppiBuildNetwork
        if (requireNamespace("igraph", quietly = TRUE)) {
          # 1. Get current Ensembl IDs from the network nodes
          node_ids <- igraph::V(exNetwork)$name
          # 2. Match these IDs back to our res_mock mapping table
          symbol_map <- res_mock$hgncSymbol[match(node_ids, rownames(res_mock))]
          # 3. Explicitly set hgncSymbol as a character attribute
          igraph::vertex_attr(exNetwork, "hgncSymbol") <- as.character(symbol_map)
          message("    - Diagnostic: Network labels localized for ", length(node_ids), " nodes.")
        }

        pLink_net <- ppiPlotNetwork(
          network = exNetwork,
          title = "Cluster 1 Matrisome PPI Network (Signal-Mapped)",
          fillColumn = log2FoldChange,
          fillType = "foldChange",
          label = TRUE,
          labelColumn = hgncSymbol,
          legend = TRUE
        )
        
        if (!is.null(pLink_net)) {
          print(pLink_net) # Show in RStudio Plot window
          ggsave(file.path(output_dir, "pathlinkR_Network_C1_Matrisome.pdf"), pLink_net, width = 10, height = 8)
          message("    - pathlinkR analysis completed. Network and table displayed.")
        }
      }
    }, error = function(e) {
      message("    Warning: pathlinkR signaling flow encountered an error: ", e$message)
    })
  } else {
    message("    Warning: No genes could be mapped to pathlinkR standard IDs (ENSG). Check mapping bridge.")
  }
}

# --- 6.3 pathlinkR Analysis for Full Cluster 1 ---
message("\n>>> Performing pathlinkR analysis on Full Cluster 1...")
# Full Cluster 1 (all proteins)
if (exists("cluster_genes") && "1" %in% names(cluster_genes)) {
  c1_genes <- cluster_genes[["1"]]
  c1_full_features <- feature_filtered[feature_filtered$UniqueSymbol %in% c1_genes, ]
  message("    - Full Cluster 1 size: ", nrow(c1_full_features), " genes.")

  # 1. Use the internal mapping file bridge
  if (!exists("m_file")) {
    data("mappingFile", package = "pathlinkR", envir = environment())
    m_file <- as.data.frame(get("mappingFile"))
  }
  
  c1_full_features$HumanSymbol <- toupper(c1_full_features$Symbol)
  mapping_bridge_full <- merge(c1_full_features, m_file, by.x = "HumanSymbol", by.y = "hgncSymbol", all.x = TRUE)
  mapping_bridge_full <- mapping_bridge_full[!is.na(mapping_bridge_full$ensemblGeneId), ]
  mapping_bridge_full <- mapping_bridge_full[!duplicated(mapping_bridge_full$ensemblGeneId), ]

  message("    - Full Cluster 1 Human ENSEMBL (ENSG) bridge success: ", nrow(mapping_bridge_full), " / ", nrow(c1_full_features), " genes.")

  if (nrow(mapping_bridge_full) > 0) {
    res_mock_full <- data.frame(
      log2FoldChange = rep(2.0, nrow(mapping_bridge_full)), 
      padj = rep(0.01, nrow(mapping_bridge_full)),
      hgncSymbol = as.character(mapping_bridge_full$Symbol),
      stringsAsFactors = FALSE
    )
    rownames(res_mock_full) <- mapping_bridge_full$ensemblGeneId
    
    tryCatch({
      # 2. Build PPI Network
      message("    - Building Full Cluster 1 PPI Network...")
      exNetworkFull <- ppiBuildNetwork(
        rnaseqResult = res_mock_full,
        filterInput = FALSE,
        columnFC = "log2FoldChange",
        columnP = "padj",
        order = "zero"
      )
      
      if (!is.null(exNetworkFull)) {
        # 3. Enrich Network
        message("    - Performing Full Cluster 1 Network Enrichment...")
        pLink_results_full <- ppiEnrichNetwork(exNetworkFull, analysis = "sigora")
        
        if (!is.null(pLink_results_full) && nrow(pLink_results_full) > 0) {
          message("    - Full Cluster 1 Pathway enrichment: Found ", nrow(pLink_results_full), " enriched pathways.")
          print(pLink_results_full)
          write.csv(pLink_results_full, file.path(output_dir, "pathlinkR_Enriched_Pathways_C1_Full.csv"), row.names = FALSE)
        } else {
          message("    - Notice: No enriched pathways found for Full Cluster 1.")
        }
        
        # Fixed: Manual Attribute Injection for robust labeling
        if (requireNamespace("igraph", quietly = TRUE)) {
          node_ids_full <- igraph::V(exNetworkFull)$name
          symbol_map_full <- res_mock_full$hgncSymbol[match(node_ids_full, rownames(res_mock_full))]
          igraph::vertex_attr(exNetworkFull, "hgncSymbol") <- as.character(symbol_map_full)
          message("    - Diagnostic: Full network labels localized for ", length(node_ids_full), " nodes.")
        }

        pLink_net_full <- ppiPlotNetwork(
          network = exNetworkFull,
          title = "Full Cluster 1 PPI Network (Signal-Mapped)",
          fillColumn = log2FoldChange,
          fillType = "foldChange",
          label = TRUE,
          labelColumn = hgncSymbol,
          legend = TRUE
        )
        
        if (!is.null(pLink_net_full)) {
          print(pLink_net_full)
          ggsave(file.path(output_dir, "pathlinkR_Network_C1_Full.pdf"), pLink_net_full, width = 10, height = 8)
          message("    - Full Cluster 1 pathlinkR basic analysis completed.")
        }
        
        # --- Advanced Diagnostics: Collagen Highlight & Subnetwork ---
        # 1. Highlight Collagen Pathway (R-HSA-1650814) in Full C1 Graph
        col_row <- pLink_results_full[as.character(pLink_results_full$pathwayId) == "R-HSA-1650814", ]
        if (nrow(col_row) > 0) {
          message("    - Highlighting Collagen Pathway in Full C1 Network...")
          col_genes <- unlist(strsplit(as.character(col_row$genes), ";"))
          
          # Inject 'Category' for highlighting
          if (requireNamespace("tidygraph", quietly = TRUE)) {
            exNetworkFull <- exNetworkFull %>%
              tidygraph::mutate(Highlight = ifelse(hgncSymbol %in% col_genes, "Collagen Pathway", "Other"))
          }
          
          p_highlight <- ppiPlotNetwork(
            network = exNetworkFull,
            title = "Cluster 1: Collagen Pathway Highlighted in Full PPI",
            fillColumn = Highlight,
            fillType = "categorical", # Changed according to snippet
            catFillColours = c("Collagen Pathway" = "#D4AC0D", "Other" = "grey92"),
            label = TRUE,
            labelColumn = hgncSymbol,
            legend = TRUE
          )
          print(p_highlight)
          ggsave(file.path(output_dir, "pathlinkR_Network_C1_Full_Collagen_Highlighted.pdf"), p_highlight, width = 10, height = 8)

          # 2. Extract and Plot Collagen Subnetwork
          message("    - Extracting Collagen Subnetwork (R-HSA-1650814)...")
          # Using the more robust extraction method from user's snippet if possible
          exSubnet <- tryCatch({
            ppiExtractSubnetwork(
              network = exNetworkFull,
              pathwayEnrichmentResult = pLink_results_full,
              pathwayToExtract = as.character(col_row$pathwayName[1])
            )
          }, error = function(e) { NULL })
          
          if (!is.null(exSubnet)) {
            # Ensure labels are character
            node_ids_sub <- igraph::V(exSubnet)$name
            igraph::vertex_attr(exSubnet, "hgncSymbol") <- as.character(res_mock_full$hgncSymbol[match(node_ids_sub, rownames(res_mock_full))])
            
            p_sub <- ppiPlotNetwork(
              network = exSubnet,
              title = "Cluster 1: Collagen Biosynthesis Module",
              fillColumn = log2FoldChange,
              fillType = "foldChange",
              label = TRUE,
              labelColumn = hgncSymbol,
              legend = TRUE
            )
            print(p_sub)
            ggsave(file.path(output_dir, "pathlinkR_Network_C1_Collagen_Subnetwork.pdf"), p_sub, width = 8, height = 6)
          }
        }
        
        # --- Advanced Diagnostics: Hub Protein Analysis ---
        message("    - Identifying Transition Hubs via Centrality Analysis...")
        if (requireNamespace("igraph", quietly = TRUE)) {
          # Calculate Centrality on the underlying igraph for reporting
          deg <- igraph::degree(exNetworkFull)
          bet <- igraph::betweenness(exNetworkFull)
          # Normalized score (Combined Degree and Betweenness)
          hub_score <- (deg/max(deg)) + (bet/max(bet))
          hub_df <- data.frame(name = names(deg), HubScore = hub_score) %>% arrange(desc(HubScore))
          top_hubs <- head(hub_df, 10)
          
          # Map Hub IDs to Symbols
          top_hubs$Symbol <- as.character(res_mock_full$hgncSymbol[match(top_hubs$name, rownames(res_mock_full))])
          message("    - Top 10 Transition Hubs identified: ", paste(top_hubs$Symbol, collapse=", "))
          write.csv(top_hubs, file.path(output_dir, "pathlinkR_Hub_Proteins_C1.csv"), row.names = FALSE)
          
          # Note: pathlinkR ppiPlotNetwork handles blue label highlighting internally for Hubs
          # when hubMeasure is used correctly, but we can also manually label them
        }
      }
    }, error = function(e) {
      message("    Warning: Full Cluster 1 advanced analysis failed: ", e$message)
    })
  } else {
    message("    Warning: No Full Cluster 1 genes could be mapped to pathlinkR standard IDs (ENSG).")
  }
}

# --- 6.4 pathlinkR Comparative Analysis for Cluster 4 ---
message("\n>>> Performing pathlinkR Comparative Analysis for Cluster 4...")
if (exists("cluster_genes") && "4" %in% names(cluster_genes)) {
  c4_genes <- cluster_genes[["4"]]
  c4_features <- feature_filtered[feature_filtered$UniqueSymbol %in% c4_genes, ]
  
  # C4 Bridge
  c4_features$HumanSymbol <- toupper(c4_features$Symbol)
  mapping_bridge_c4 <- merge(c4_features, m_file, by.x = "HumanSymbol", by.y = "hgncSymbol", all.x = TRUE)
  mapping_bridge_c4 <- mapping_bridge_c4[!is.na(mapping_bridge_c4$ensemblGeneId), ]
  mapping_bridge_c4 <- mapping_bridge_c4[!duplicated(mapping_bridge_c4$ensemblGeneId), ]
  
  if (nrow(mapping_bridge_c4) > 0) {
    res_mock_c4 <- data.frame(
      log2FoldChange = rep(1.5, nrow(mapping_bridge_c4)), 
      padj = rep(0.01, nrow(mapping_bridge_c4)),
      hgncSymbol = as.character(mapping_bridge_c4$Symbol),
      stringsAsFactors = FALSE
    )
    rownames(res_mock_c4) <- mapping_bridge_c4$ensemblGeneId
    
    tryCatch({
      message("    - Building Cluster 4 PPI Network...")
      # Fixed: order='first' (valid options: zero, first, minSimple)
      exNetworkC4 <- ppiBuildNetwork(
        res_mock_c4, 
        filterInput = FALSE, 
        order = "first",
        hubMeasure = "betweenness" # Added for hub identification
      ) 
      
      if (!is.null(exNetworkC4)) {
        # Manual Label Inject for C4
        node_ids_c4 <- igraph::V(exNetworkC4)$name
        igraph::vertex_attr(exNetworkC4, "hgncSymbol") <- as.character(res_mock_c4$hgncSymbol[match(node_ids_c4, rownames(res_mock_c4))])
        
        # Enrichment for C4
        pLink_results_c4 <- ppiEnrichNetwork(exNetworkC4, analysis = "sigora")
        if (!is.null(pLink_results_c4) && nrow(pLink_results_c4) > 0) {
          message("    - Cluster 4 Pathway enrichment: Found ", nrow(pLink_results_c4), " pathways.")
          print(head(pLink_results_c4))
          write.csv(pLink_results_c4, file.path(output_dir, "pathlinkR_Enriched_Pathways_C4.csv"), row.names = FALSE)
        }
        
        # Plot C4
        p_c4 <- ppiPlotNetwork(
          network = exNetworkC4,
          title = "Full Cluster 4 PPI Network (Comparative - Order 1)",
          fillColumn = log2FoldChange,
          fillType = "foldChange",
          label = TRUE,
          labelColumn = hgncSymbol,
          legend = TRUE
        )
        print(p_c4)
        ggsave(file.path(output_dir, "pathlinkR_Network_C4_Full.pdf"), p_c4, width = 10, height = 8)

        # --- Cluster 4 Refinements: NA Clarification & TGFB1 Axis ---
        message("    - Refining Cluster 4 Visualization & TGFB1 logic...")
        
        # 1. Distinguish Proteomics vs Database (NA) nodes
        if (requireNamespace("tidygraph", quietly = TRUE)) {
          exNetworkC4 <- exNetworkC4 %>%
            tidygraph::mutate(DataSource = ifelse(is.na(log2FoldChange), "Database Bridge", "Proteomics Dataset"))
        }

        p_c4_source <- ppiPlotNetwork(
          network = exNetworkC4,
          title = "Cluster 4: PPI Network with Signal Bridges (Grey = DB)",
          fillColumn = DataSource,
          fillType = "categorical",
          catFillColours = c("Proteomics Dataset" = "#C0392B", "Database Bridge" = "grey85"),
          label = TRUE,
          labelColumn = hgncSymbol,
          legend = TRUE
        )
        print(p_c4_source)
        ggsave(file.path(output_dir, "pathlinkR_Network_C4_Refined_Sources.pdf"), p_c4_source, width = 10, height = 8)

        # 2. TGFB1 Axis Extraction
        tgfb1_row <- pLink_results_c4[grep("TGFB1", as.character(pLink_results_c4$genes)), ]
        if (nrow(tgfb1_row) > 0) {
          message("    - Extracting TGFB1 Signaling Subnetwork...")
          exSubnet_TGFB1 <- tryCatch({
            ppiExtractSubnetwork(
              network = exNetworkC4,
              pathwayEnrichmentResult = pLink_results_c4,
              pathwayToExtract = as.character(tgfb1_row$pathwayName[1])
            )
          }, error = function(e) { NULL })
          
          if (!is.null(exSubnet_TGFB1)) {
            # Ensure labels for subnetwork
            node_ids_tgf <- igraph::V(exSubnet_TGFB1)$name
            igraph::vertex_attr(exSubnet_TGFB1, "hgncSymbol") <- as.character(res_mock_c4$hgncSymbol[match(node_ids_tgf, rownames(res_mock_c4))])
            
            p_tgf <- ppiPlotNetwork(
              network = exSubnet_TGFB1,
              title = paste0("Cluster 4: TGFB1 Axis (", tgfb1_row$pathwayName[1], ")"),
              fillType = "oneSided",
              fillColumn = degree,
              label = TRUE,
              labelColumn = hgncSymbol,
              legendTitle = "Degree"
            )
            print(p_tgf)
            ggsave(file.path(output_dir, "pathlinkR_Network_C4_TGFB1_Subnetwork.pdf"), p_tgf, width = 8, height = 6)
          }
        }
      }
    }, error = function(e) {
      message("    Warning: Cluster 4 analysis failed: ", e$message)
    })
  }
}


message("\n>>> Starting Pathway-Filtered Functional Core Analysis...")
pathway_files <- c(
  "lauer_style_comparison/Participating_Molecules_[R-MMU-109582].tsv",
  "lauer_style_comparison/Participating_Molecules_[R-MMU-1280215].tsv",
  "lauer_style_comparison/Participating_Molecules_[R-MMU-1474244].tsv"
)

# Robust loader for symbols in "UniProt:ID Symbol" format
load_pathway_symbols <- function(files) {
  all_symbols <- c()
  for (f in files) {
    if (file.exists(f)) {
      df <- read.delim(f, check.names = FALSE)
      # Extract symbols: MoleculeName usually like "UniProt:Q80T91 Symbol"
      raw_names <- df$MoleculeName
      symbols <- sapply(strsplit(as.character(raw_names), " "), function(x) {
        if (length(x) > 1) return(x[2]) else return(NA)
      })
      all_symbols <- c(all_symbols, na.omit(symbols))
    }
  }
  return(unique(all_symbols))
}

pathway_ref_list <- load_pathway_symbols(pathway_files)
message("    Total Pathway Reference Symbols loaded: ", length(pathway_ref_list))

# 2. Extract Shared Functional Core (Cl 2 & 3)
analyze_shared_processes <- function(c_pair, p8_results, mf_data) {
  c1 <- paste0("Cluster_", c_pair[1]); c2 <- paste0("Cluster_", c_pair[2])
  if (is.null(p8_results[[c1]]) || is.null(p8_results[[c2]])) return(NULL)
  shared_ids <- intersect(p8_results[[c1]]$ID, p8_results[[c2]]$ID)
  if (length(shared_ids) == 0) return(NULL)
  
  shared_info <- p8_results[[c1]][p8_results[[c1]]$ID %in% shared_ids, c("ID", "Description", "geneID")]
  all_genes <- unique(unlist(strsplit(shared_info$geneID, "/")))
  all_genes <- intersect(all_genes, rownames(mf_data))
  return(all_genes)
}

core_genes_23 <- analyze_shared_processes(c(2, 3), p8_results, mfuzz_input)

# 3. Intersection: Functional Core INTERSECT Reactome Molecules
final_transition_drivers <- intersect(core_genes_23, pathway_ref_list)
message("    Final Pathway-Filtered Transition Drivers: ", length(final_transition_drivers))
message("    Drivers: ", paste(final_transition_drivers, collapse=", "))

# 4. Dual Highlighted Volcano Plot Function
plot_pathway_highlight_volcano <- function(res, highlighted_genes, title) {
  res_plot <- res %>% data.frame() %>% tibble::rownames_to_column("ProteinID")
  res_plot$ID <- feature_filtered$UniqueSymbol[match(res_plot$ProteinID, rownames(feature_filtered))]
  
  # Map Mfuzz clusters for background context
  res_plot <- res_plot %>% left_join(p7_res$dfcluster, by = "ID")
  
  # Define highlight group
  res_plot$HighlightGroup <- "Background"
  res_plot$HighlightGroup[res_plot$ID %in% highlighted_genes] <- "Pathway_Core"
  
  # Volcano Plot Core Logic
  ggplot(res_plot, aes(x = logFC, y = -log10(P.Value))) +
    # Background: Show Mfuzz clusters with low alpha
    geom_point(aes(color = factor(Cluster)), alpha = 0.15, size = 0.8) +
    scale_color_manual(values = cluster_palette, na.value = "grey85") +
    # Highlighted: Pathway Core proteins
    geom_point(data = filter(res_plot, HighlightGroup == "Pathway_Core"), 
               color = "black", fill = "#E41A1C", shape = 21, size = 2, alpha = 0.9) +
    theme_elegant() +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(P-value)") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    ggrepel::geom_text_repel(data = filter(res_plot, HighlightGroup == "Pathway_Core"),
                            aes(label = ID), size = 3.5, fontface = "bold", 
                            box.padding = 0.5, max.overlaps = 20) +
    theme(aspect.ratio = 1, legend.position = "none")
}


# Generate Plots
if (length(final_transition_drivers) > 0) {
  # Transition: D14_C vs D7_F
  # Ensure the contrast exists in fit2 (in case user ran partial code)
  if (!"D14_C_vs_D7_F" %in% colnames(fit2)) {
    stop("Contrast 'D14_C_vs_D7_F' not found in fit2. Please re-run the Section 3 setup.")
  }
  
  p19a <- plot_pathway_highlight_volcano(topTable(fit2, coef="D14_C_vs_D7_F", number=Inf), 
                                        final_transition_drivers, 
                                        "Transition Core: D14C vs D7F (Pathway Filtered)")
  
  # Maturation: D14_C vs D7_C
  p19b <- plot_pathway_highlight_volcano(topTable(fit2, coef="D14_C_vs_D7_C", number=Inf), 
                                        final_transition_drivers, 
                                        "Maturation Core: D14C vs D7C (Pathway Filtered)")
  
  print(p19a)
  print(p19b)
  
  # 5. Pathway-Specific Hub Heatmaps (Lauer Style - Group Average)
  message(">>> Generating Pathway-Specific Hub Heatmaps (GO Core Intersection)...")
  
  pathway_config <- list(
    hemostasis = list(file = "lauer_style_comparison/Participating_Molecules_[R-MMU-109582].tsv", name = "Hemostasis Hubs (GO Core)"),
    cytokine = list(file = "lauer_style_comparison/Participating_Molecules_[R-MMU-1280215].tsv", name = "Cytokine Signaling Hubs (GO Core)"),
    ecm = list(file = "lauer_style_comparison/Participating_Molecules_[R-MMU-1474244].tsv", name = "ECM Organization Hubs (GO Core)")
  )
  
  # Target levels for group order
  target_order_levels <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall")
  
  for (p_key in names(pathway_config)) {
    p_info <- pathway_config[[p_key]]
    message("    - Processing: ", p_info$name)
    
    # 1. Load symbols from this specific TSV
    spec_path_genes <- load_pathway_symbols(list(p_info$file))
    
    # 2. Intersect with the Shared GO Core (core_genes_23)
    path_hubs <- intersect(core_genes_23, spec_path_genes)
    path_hubs <- na.omit(as.character(path_hubs))
    path_hubs <- path_hubs[path_hubs != "NA"]
    
    if (length(path_hubs) >= 2) {
      # 3. Find protein indices in our data
      # Filter UniqueSymbol for matches and ensure no NA
      match_mask <- (feature_filtered$UniqueSymbol %in% path_hubs) & (!is.na(feature_filtered$UniqueSymbol))
      match_indices <- which(match_mask)
      
      if (length(match_indices) >= 2) {
        # 4. Extract Expression Data
        p_expr_mat <- expr_norm[match_indices, , drop = FALSE]
        rownames(p_expr_mat) <- as.character(feature_filtered$UniqueSymbol[match_indices])
        
        # 5. Call the Lauer-style expression heatmap function
        p_heat <- plot_lauer_style_expression_heatmap(
          data = p_expr_mat, 
          metadata = metadata, 
          title = p_info$name,
          order_levels = target_order_levels
        )
        
        if (!is.null(p_heat)) {
           grid::grid.newpage()
           grid::grid.draw(p_heat$gtable)
           assign(paste0("p21_", p_key), p_heat, envir = .GlobalEnv)
        }
      } else {
        message("      ! Notice: Too few matching proteins in data matrix.")
      }
    } else {
      message("      ! Notice: No hub proteins found for this pathway.")
    }
  }
}

# ------------------------------------------------------------------------------
# 8. ECM Sub-pathway Abundance Trend Analysis
# ------------------------------------------------------------------------------
message("\n>>> Starting ECM Sub-pathway Abundance Trend Analysis...")

# 1. Configuration: Mapping Reactome IDs to Pathway Names
ecm_sub_pathways <- list(
  "Collagen formation" = "R-MMU-1474290",
  "Fibronectin matrix formation" = "R-MMU-1566977",
  "Elastic fibre formation" = "R-MMU-1566948",
  "Laminin interactions" = "R-MMU-3000157",
  "Non-integrin membrane-ECM interactions" = "R-MMU-3000171",
  "ECM proteoglycans" = "R-MMU-3000178",
  "Degradation of the extracellular matrix" = "R-MMU-1474228",
  "Integrin cell surface interactions" = "R-MMU-216083",
  "Invadopodia formation" = "R-MMU-8941237"
)

# 2. Extract Protein Lists for each Sub-pathway
# Try multiple potential paths to find the ECM reactome folder
possible_ecm_paths <- c(
  "ECM reactome", 
  "lauer_style_comparison/ECM reactome",
  "C:/Users/SimonYao/Desktop/LCM_protemics_thrombus/lauer_style_comparison/ECM reactome"
)
ecm_folder <- possible_ecm_paths[dir.exists(possible_ecm_paths)][1]

if (is.na(ecm_folder)) {
  stop(">>> Error: Could not find 'ECM reactome' folder. Please check your working directory.")
} else {
  message(">>> Using ECM folder: ", ecm_folder)
}

ecm_pathway_proteins <- list()

for (p_name in names(ecm_sub_pathways)) {
  p_id <- ecm_sub_pathways[[p_name]]
  f_name <- paste0("Participating_Molecules_[", p_id, "].tsv")
  f_path <- file.path(ecm_folder, f_name)
  
  if (file.exists(f_path)) {
    # Extract SYMBOLS using the existing robust loader
    ecm_pathway_proteins[[p_name]] <- load_pathway_symbols(list(f_path))
  } else {
    message("    ! Warning: File not found for ", p_name, " (", f_path, ")")
  }
}

# 3. Calculate Average Expression (Z-score) per Pathway
# Use expr_norm (log2, normalized) for relative abundance trends
abundance_list <- list()
for (p_name in names(ecm_pathway_proteins)) {
  # Skip specific category as requested
  if (p_name == "Non-integrin membrane-ECM interactions") next
  
  target_symbols <- ecm_pathway_proteins[[p_name]]
  matching_indices <- which(feature_filtered$UniqueSymbol %in% target_symbols)
  
  if (length(matching_indices) > 0) {
    # Calculate average log2 expression for this pathway across samples
    p_expr_avg <- colMeans(expr_norm[matching_indices, , drop = FALSE], na.rm = TRUE)
    abundance_list[[p_name]] <- p_expr_avg
    message("    - ", p_name, ": matched ", length(matching_indices), " proteins.")
  }
}

# 4. Prepare Data for Plotting
message(">>> Assembling trend data...")
if (length(abundance_list) == 0) {
  stop("!!! Error: No abundance data calculated.")
}

abundance_df <- bind_rows(lapply(names(abundance_list), function(n) {
  vec <- abundance_list[[n]]
  data.frame(SampleID = names(vec), 
             Value = as.numeric(vec), 
             Pathway = n,
             stringsAsFactors = FALSE)
}))

# 5. Join with metadata and Filter Groups
# Explicitly filter out High Inflammation (which often appears as NA in target_order_levels)
plot_data_ecm <- abundance_df %>%
  left_join(metadata %>% select(SampleID, Group), by = "SampleID") %>%
  filter(!is.na(Group)) %>%
  filter(Group %in% target_order_levels) %>% # Only keep the 6 core groups (D2-F to D14-Wall)
  mutate(Group = factor(Group, levels = target_order_levels)) %>%
  group_by(Group, Pathway) %>%
  summarise(Mean_Val = mean(Value, na.rm = TRUE),
            SE = sd(Value, na.rm = TRUE) / sqrt(n()), .groups = "drop")

# Calculate Z-score per Pathway to show relative trends (since absolute abundance is low)
plot_data_ecm <- plot_data_ecm %>%
  group_by(Pathway) %>%
  mutate(Mean_Z = (Mean_Val - mean(Mean_Val, na.rm=TRUE)) / sd(Mean_Val, na.rm=TRUE),
         SE_Z = SE / sd(Mean_Val, na.rm=TRUE)) %>%
  ungroup()

# 6. Visualization: Two Targeted ECM Trend Plots
message(">>> Generating Split ECM Trend Plots (p22a, p22b)...")

# Biologically relevant grouping:
# Group A: Structural Scaffold (Physical components)
scaffold_pathways <- c("Collagen formation", "Fibronectin matrix formation", 
                      "Elastic fibre formation", "ECM proteoglycans")

# Group B: Functional Remodeling & Signaling (Active processes)
functional_pathways <- c("Degradation of the extracellular matrix", 
                        "Integrin cell surface interactions", 
                        "Laminin interactions", "Invadopodia formation")

# Use original palette for consistency
ecm_colors <- c(
  "Collagen formation" = "#E41A1C",
  "Fibronectin matrix formation" = "#377EB8",
  "Elastic fibre formation" = "#4DAF4A",
  "Laminin interactions" = "#984EA3",
  "ECM proteoglycans" = "#FFFF33",
  "Degradation of the extracellular matrix" = "#A65628",
  "Integrin cell surface interactions" = "#F781BF",
  "Invadopodia formation" = "#999999"
)

# Plot A: Structural Scaffold
p22a <- ggplot(filter(plot_data_ecm, Pathway %in% scaffold_pathways), 
               aes(x = Group, y = Mean_Z, color = Pathway, group = Pathway)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = ecm_colors) +
  theme_elegant() +
  labs(title = "ECM Structural Scaffold Trends",
       y = "Relative Abundance (Z-score)", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Plot B: Functional Remodeling & Signaling
p22b <- ggplot(filter(plot_data_ecm, Pathway %in% functional_pathways), 
               aes(x = Group, y = Mean_Z, color = Pathway, group = Pathway)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = ecm_colors) +
  theme_elegant() +
  labs(title = "ECM Remodeling & Interaction Trends",
       y = "Relative Abundance (Z-score)", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"))

print(p22a)
print(p22b)

# Save separately
ggsave(file.path(output_dir, "Figure2_ECM_Structural_Trends.pdf"), p22a, width = 8, height = 6)
ggsave(file.path(output_dir, "Figure2_ECM_Functional_Trends.pdf"), p22b, width = 8, height = 6)

# ------------------------------------------------------------------------------
# 9. One-Click Results Export (CSV)
# ------------------------------------------------------------------------------
message("\n>>> Exporting Analysis Results to CSV...")

export_dir <- file.path(output_dir, "CSV_Results")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Export Pathway-Validated Transition Drivers
if (exists("final_transition_drivers")) {
  drivers_df <- feature_filtered %>%
    filter(UniqueSymbol %in% final_transition_drivers) %>%
    select(UniqueSymbol, ProteinID = PG.ProteinGroups, Description = PG.ProteinDescriptions) %>%
    mutate(Analysis_Group = "Pathway_Validated_Transition_Driver")
  
  write.csv(drivers_df, file.path(export_dir, "Table1_Transition_Drivers.csv"), row.names = FALSE)
  message("    - Exported: Table1_Transition_Drivers.csv")
}

# 2. Export Pathway-Specific Hubs (Heatmap Proteins)
if (exists("pathway_config")) {
  all_hubs <- list()
  for (p_key in names(pathway_config)) {
    spec_path_genes <- load_pathway_symbols(list(pathway_config[[p_key]]$file))
    hubs <- intersect(core_genes_23, spec_path_genes)
    if (length(hubs) > 0) {
      all_hubs[[p_key]] <- data.frame(Pathway = pathway_config[[p_key]]$name, 
                                     Symbol = hubs)
    }
  }
  hubs_df <- bind_rows(all_hubs)
  write.csv(hubs_df, file.path(export_dir, "Table2_Pathway_Hubs.csv"), row.names = FALSE)
  message("    - Exported: Table2_Pathway_Hubs.csv")
}

# 3. Export ECM Sub-pathway Trend Data (Z-scores)
if (exists("plot_data_ecm")) {
  write.csv(plot_data_ecm, file.path(export_dir, "Table3_ECM_Trend_Data.csv"), row.names = FALSE)
  message("    - Exported: Table3_ECM_Trend_Data.csv")
}

# 4. Export DEA Results for Key Contrasts (Top 500 each)
# Note: Using fit2 from Section 7
if (exists("fit2")) {
  key_contrasts <- c("D14_C_vs_D7_F", "D14_C_vs_D7_C")
  for (con in key_contrasts) {
    if (con %in% colnames(fit2)) {
      res <- topTable(fit2, coef = con, number = 500)
      write.csv(res, file.path(export_dir, paste0("DEA_", con, ".csv")), row.names = TRUE)
    }
  }
  message("    - Exported: Key DEA Results (Top 500)")
}

message("\n>>> ALL EXPORTS COMPLETED. Files are located in: ", export_dir)
