# D14 Spatial Heterogeneity Analysis: Fibrin vs Collagen vs Wall

# ------------------------------------------------------------------------------
# 1. Setup and Configuration
# ------------------------------------------------------------------------------
output_dir <- "lauer_style_comparison/output_D14_Spatial"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load libraries
suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(patchwork)
  library(ReactomePA)
  library(gridExtra)
})

# Source helper scripts
source("src/data_loader.R")

# ------------------------------------------------------------------------------
# 2. Data Loading & Global Preprocessing (Bias Avoidance)
# ------------------------------------------------------------------------------
message(">>> Loading Data...")
input_file <- "thrombusDIAreport.xls"
if (!file.exists(input_file)) input_file <- "c:/Users/SimonYao/Desktop/LCM_protemics_thrombus/thrombusDIAreport.xls"

data_list <- load_and_clean_data(input_file)
expr_data <- data_list$exprs
feature_data <- data_list$feature_data
metadata <- data_list$metadata

# 2.1 GLOBAL FILTERING (Consistent with Main Analysis)
# This ensures we define "valid proteins" based on the whole experiment, 
# preventing bias from subset-specific missingness patterns.
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

# 2.2 Standardize Gene Symbols (Robust extraction)
# 1. Identify the best column for Gene Symbols
gene_col_idx <- grep("PG.Genes|Genes|Gene|Symbol", names(feature_filtered), ignore.case = TRUE)[1]
if (!is.na(gene_col_idx)) {
  feature_filtered$Symbol <- as.character(feature_filtered[[gene_col_idx]])
} else {
  feature_filtered$Symbol <- NA
}

# 2. Fill missing symbols from Descriptions (GN= pattern)
missing_gene <- is.na(feature_filtered$Symbol) | feature_filtered$Symbol == "" | feature_filtered$Symbol == "NA"
if (any(missing_gene)) {
  descriptions <- feature_filtered$PG.ProteinDescriptions[missing_gene]
  extracted_genes <- stringr::str_extract(descriptions, "GN=[^ ]+")
  extracted_genes <- gsub("GN=", "", extracted_genes)
  # Fallback to Protein ID if still missing
  feature_filtered$Symbol[missing_gene] <- ifelse(!is.na(extracted_genes), extracted_genes, 
                                                 sapply(strsplit(as.character(rownames(feature_filtered)[missing_gene]), ";"), `[`, 1))
}

feature_filtered$Symbol <- stringr::str_to_title(as.character(feature_filtered$Symbol))
feature_filtered$UniqueSymbol <- make.unique(feature_filtered$Symbol)

# 2.3 GLOBAL LOG2 & IMPUTATION
# Imputing using global distribution parameters prevents "shifting" the D14 data 
# relative to D7/D2 if they were processed in isolation.
expr_filtered[expr_filtered == 0] <- NA
expr_log <- log2(expr_filtered)

set.seed(123)
impute_minprob <- function(data) {
  # Calculate stats on the FULL dataset
  valid_vals <- data[!is.na(data)]
  mu <- mean(valid_vals) - 1.8 * sd(valid_vals)
  sig <- 0.3 * sd(valid_vals)
  
  data_imp <- data
  is_na <- is.na(data)
  data_imp[is_na] <- rnorm(sum(is_na), mean = mu, sd = sig)
  return(data_imp)
}

expr_imp <- impute_minprob(expr_log)

# 2.3 GLOBAL NORMALIZATION
# Normalizing D14 against the global set ensures comparability.
expr_norm_global <- normalizeQuantiles(expr_imp)
rownames(expr_norm_global) <- rownames(expr_filtered)

# 2.4 SUBSETTING FOR SPATIAL-TEMPORAL ANALYSIS
message(">>> Subsetting for D2, D7, and D14 (Fibrin, Collagen, Wall)...")
target_groups <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall")
target_mask <- metadata$Group %in% target_groups

expr_subset <- expr_norm_global[, target_mask]
metadata_subset <- metadata[target_mask, ]
feature_subset <- feature_filtered # Same rows

# Re-factor groups
metadata_subset$Group <- factor(metadata_subset$Group, levels = target_groups)
# Rename for convenience in downstream code
expr_norm <- expr_subset 

message("    - Spatial-Temporal Dataset: ", ncol(expr_norm), " samples, ", nrow(expr_norm), " proteins")

# ------------------------------------------------------------------------------
# 3. Visualization Setup (Colors & Themes)
# ------------------------------------------------------------------------------
# Color Scheme: Fibrin, Collagen, Wall (Standardized)
spatial_colors <- c(
  "D2_Fibrin"     = "#FADBD8",
  "D7_Fibrin"     = "#E6B0AA",
  "D14_Fibrin"    = "#B03A2E", 
  "D7_Collagen"   = "#AED6F1",
  "D14_Collagen"  = "#2E86C1", 
  "D14_Wall"      = "#8E44AD"
)

theme_elegant <- function() {
  theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

theme_paper <- theme_elegant()

# --- Perform Matrisome Intersection with Clusters ---
message("\n>>> Generating Matrisome vs Cluster Intersections (Venn Diagrams)...")

# CRITICAL: Using the ROOT Matrisome directory which contains the FULL database (60KB files)
matrisome_dir <- "C:/Users/SimonYao/Desktop/LCM_protemics_thrombus/Matrisome"
core_paths <- list.files(file.path(matrisome_dir, "Core Matrisome"), pattern="\\.tsv$|\\.txt$|\\.csv$", full.names=TRUE)
assoc_paths <- list.files(file.path(matrisome_dir, "Matrisome associated"), pattern="\\.tsv$|\\.txt$|\\.csv$", full.names=TRUE)

# Helper to normalize identifiers (UniProt/Gene)
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

matrisome_full_db <- load_local_db_proteins_as_list(c(core_paths, assoc_paths))

# Helper to map identifiers to our features robustly - RETURNS UniqueSymbol
get_genes_from_ids <- function(id_list) {
  valid_db_ids <- unique(id_list) 
  if (length(valid_db_ids) == 0) return(character(0))
  
  # Prepare searchable tokens for each protein in our dataset
  match_mask <- rep(FALSE, nrow(feature_filtered))
  
  for (i in 1:nrow(feature_filtered)) {
    # Extract identifiers from multiple possible columns
    raw_ids <- c(
      as.character(feature_filtered$PG.Genes[i]),
      as.character(feature_filtered$PG.ProteinGroups[i]),
      as.character(feature_filtered$Symbol[i]),
      as.character(rownames(feature_filtered)[i])
    )
    
    # Split by common delimiters: ; , | / and space
    tokens <- unlist(strsplit(raw_ids, "[; ,/|]"))
    tokens <- tokens[!is.na(tokens) & tokens != ""]
    
    # Normalize and Check match
    clean_tokens <- normalize_id(tokens)
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

# ------------------------------------------------------------------------------
# 4. PCA Analysis (D14 Only)
# ------------------------------------------------------------------------------
message(">>> Running PCA...")
pca_res <- prcomp(t(expr_norm), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- metadata_subset$Group
pca_df$SampleID <- metadata_subset$SampleID

# Calculate variance exp
var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = SampleID)) +
  geom_point(size = 4, alpha = 0.8) +
  # geom_text_repel(size = 3, show.legend = FALSE) + # Optional: Hide labels for cleaner look
  stat_ellipse(level = 0.95, linetype = "dashed", size = 0.5) +
  scale_color_manual(values = spatial_colors) +
  labs(title = "PCA: D14 Spatial Heterogeneity",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_paper

print(p_pca)
ggsave(file.path(output_dir, "1_D14_Spatial_PCA.pdf"), p_pca, width = 6, height = 5)

# ------------------------------------------------------------------------------
# 5. Differential Expression (Limma: Pairwise & Global)
# ------------------------------------------------------------------------------
message(">>> Running Limma...")
design <- model.matrix(~ 0 + Group, data = metadata_subset)
colnames(design) <- levels(metadata_subset$Group)

fit <- lmFit(expr_norm, design)
contrast_matrix <- makeContrasts(
  C_vs_F = D14_Collagen - D14_Fibrin,
  W_vs_C = D14_Wall - D14_Collagen,
  W_vs_F = D14_Wall - D14_Fibrin,
  # Requested Evolution Trajectories (Target_vs_Reference)
  D7F_vs_D2F   = D7_Fibrin - D2_Fibrin,
  D7C_vs_D7F   = D7_Collagen - D7_Fibrin,
  D14C_vs_D7C  = D14_Collagen - D7_Collagen,
  D14W_vs_D14C = D14_Wall - D14_Collagen,
  D14F_vs_D7F  = D14_Fibrin - D7_Fibrin,
  D14C_vs_D14F = D14_Collagen - D14_Fibrin,
  D7C_vs_D2F   = D7_Collagen - D2_Fibrin,
  D14C_vs_D7F  = D14_Collagen - D7_Fibrin,
  D14W_vs_D14F = D14_Wall - D14_Fibrin,
  D14F_vs_D2F  = D14_Fibrin - D2_Fibrin,
  D14C_vs_D2F  = D14_Collagen - D2_Fibrin,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Global Test (F-statistic for any difference across 3 groups)
# Equivalent to testing if all 3 contrasts are zero
fit_global <- eBayes(fit) # Re-run standard? No, construct F-test from fit2
# Actually, topTableF from fit2 works for all contrasts
res_global <- topTable(fit2, number = Inf, p.value = 1) 
# Note: topTable on fit2 with multiple coefficients does an F-test
message("    - Global Analysis: ", sum(res_global$adj.P.Val < 0.05), " significant proteins (F-test FDR < 0.05)")

# ------------------------------------------------------------------------------
# 6. Pairwise Volcano Plots (Standardized)
# ------------------------------------------------------------------------------
message(">>> Generating Pairwise Volcano Plots...")
plot_volcano_custom <- function(res, title, color_up, color_down, contrast_label) {
  res$Significance <- "NS"
  res$Significance[res$logFC > 0.58 & res$adj.P.Val < 0.05] <- "Up"
  res$Significance[res$logFC < -0.58 & res$adj.P.Val < 0.05] <- "Down"
  
  # Map Symbols mapping
  res$Symbol <- feature_subset$Symbol[match(rownames(res), rownames(feature_subset))]
  
  # Top labels selection (Top 10 up, Top 10 down)
  top_up <- res %>% filter(Significance == "Up") %>% arrange(adj.P.Val) %>% head(10)
  top_down <- res %>% filter(Significance == "Down") %>% arrange(adj.P.Val) %>% head(10)
  label_df <- rbind(top_up, top_down)
  
  ggplot(res, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(data = subset(res, Significance == "NS"), color = "grey85", alpha = 0.4, size = 1.2) +
    geom_point(data = subset(res, Significance == "Up"), color = color_up, alpha = 0.7, size = 2) +
    geom_point(data = subset(res, Significance == "Down"), color = color_down, alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey40", size = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.3) +
    geom_text_repel(data = label_df, aes(label = Symbol), size = 3, fontface = "bold", 
                    box.padding = 0.5, max.overlaps = 50) +
    labs(title = title, x = expression(bold(Log[2]~Fold~Change)), y = expression(bold(-Log[10]~Adj.P.Value))) +
    theme_elegant() +
    theme(legend.position = "none")
}

# 1. C vs F (D14 Collagen vs Fibrin)
res_CF <- topTable(fit2, coef = "C_vs_F", number = Inf)
p_v1_new <- plot_volcano_custom(res_CF, "D14: Collagen vs Fibrin", spatial_colors["D14_Collagen"], spatial_colors["D14_Fibrin"], "C_vs_F")

# 2. W vs C (D14 Wall vs Collagen)
res_WC <- topTable(fit2, coef = "W_vs_C", number = Inf)
p_v2_new <- plot_volcano_custom(res_WC, "D14: Wall vs Collagen", spatial_colors["D14_Wall"], spatial_colors["D14_Collagen"], "W_vs_C")

# 3. W vs F (D14 Wall vs Fibrin)
res_WF <- topTable(fit2, coef = "W_vs_F", number = Inf)
p_v3_new <- plot_volcano_custom(res_WF, "D14: Wall vs Fibrin", spatial_colors["D14_Wall"], spatial_colors["D14_Fibrin"], "W_vs_F")

print(p_v1_new); print(p_v2_new); print(p_v3_new)
ggsave(file.path(output_dir, "2_Volcano_D14_C_vs_F.pdf"), p_v1_new, width = 6, height = 5.5)
ggsave(file.path(output_dir, "3_Volcano_D14_W_vs_C.pdf"), p_v2_new, width = 6, height = 5.5)
ggsave(file.path(output_dir, "4_Volcano_D14_W_vs_F.pdf"), p_v3_new, width = 6, height = 5.5)

# ------------------------------------------------------------------------------
# 7. Optimized Ternary Plot (Niche & Boundary)
# ------------------------------------------------------------------------------
message(">>> Generating Optimized Ternary Plot...")

# 1. Coordinate & Abundance Calculation
# Calculate Mean Expression per Group (Linear scale for barycentric coords)
expr_lin <- 2^expr_norm
group_means <- data.frame(row.names = rownames(expr_lin))
for (g in target_groups) {
  group_means[[g]] <- rowMeans(expr_lin[, metadata_d14$Group == g], na.rm = TRUE)
}

# Ternary Transformation (Normalize to sum to 1)
row_sums <- rowSums(group_means)
ternary_df <- group_means / row_sums

# Coordinate Transform (Barycentric to Cartesian)
# Top=Wall, Left=Fibrin, Right=Collagen
ternary_df$x <- ternary_df$D14_Collagen + 0.5 * ternary_df$D14_Wall
ternary_df$y <- (sqrt(3) / 2) * ternary_df$D14_Wall

# 2. Assign Metadata for Plotting
ternary_df$Gene <- feature_d14$Symbol[match(rownames(ternary_df), rownames(feature_d14))]
ternary_df$Dominant <- apply(group_means, 1, function(x) target_groups[which.max(x)])
ternary_df$adj.P.Val <- res_global$adj.P.Val[match(rownames(ternary_df), rownames(res_global))]
ternary_df$MeanLog2 <- rowMeans(expr_norm)

# Scale MeanLog2 for Alpha (0.3 to 1.0)
ternary_df$PlotAlpha <- 0.3 + 0.7 * (ternary_df$MeanLog2 - min(ternary_df$MeanLog2)) / (max(ternary_df$MeanLog2) - min(ternary_df$MeanLog2))

# 3. Identify Boundary Proteins (Shared by two niches, low in third)
ratios <- apply(group_means, 1, function(x) {
  xs <- sort(x, decreasing = TRUE)
  return(c(Max_to_Mid = xs[1]/xs[2], Max_to_Min = xs[1]/xs[3]))
})
ternary_df$Max_to_Mid <- ratios[1, ]
ternary_df$Max_to_Min <- ratios[2, ]

# Classification
ternary_df$Niche_Status <- "Shared/NS"
ternary_df$Niche_Status[ternary_df$adj.P.Val < 0.05 & ternary_df$Max_to_Mid > 1.5] <- "Niche Specific"
ternary_df$Niche_Status[ternary_df$adj.P.Val < 0.05 & ternary_df$Max_to_Mid <= 1.5 & ternary_df$Max_to_Min > 2] <- "Boundary"

# 4. Refine Plot Colors
ext_colors_opt <- c(spatial_colors, "Background" = "#D3D3D3")
ternary_df$FinalColor <- "Background"
ternary_df$FinalColor[ternary_df$Niche_Status != "Shared/NS"] <- as.character(ternary_df$Dominant[ternary_df$Niche_Status != "Shared/NS"])
ternary_df$FinalColor <- factor(ternary_df$FinalColor, levels = c(target_groups, "Background"))

# 5. Label Selection
# Label Top 25 niche-specific markers AND ALL Boundary markers
top_specific <- ternary_df %>% filter(Niche_Status == "Niche Specific") %>% arrange(adj.P.Val) %>% head(25)
all_boundary <- ternary_df %>% filter(Niche_Status == "Boundary")
label_set <- unique(c(rownames(top_specific), rownames(all_boundary)))
ternary_df$Label <- ifelse(rownames(ternary_df) %in% label_set, as.character(ternary_df$Gene), NA)

p_ternary_opt <- ggplot(ternary_df, aes(x = x, y = y)) +
  # Triangle Borders & Separators (Using annotate to avoid length errors)
  annotate("segment", x = 0, y = 0, xend = 1, yend = 0, color = "black", size = 0.6) +
  annotate("segment", x = 1, y = 0, xend = 0.5, yend = sqrt(3)/2, color = "black", size = 0.6) +
  annotate("segment", x = 0.5, y = sqrt(3)/2, xend = 0, yend = 0, color = "black", size = 0.6) +
  annotate("segment", x = 0.5, y = sqrt(3)/6, xend = 0.5, yend = 0, linetype = "dotted", color = "grey70") +
  annotate("segment", x = 0.5, y = sqrt(3)/6, xend = 0.25, yend = sqrt(3)/4, linetype = "dotted", color = "grey70") +
  annotate("segment", x = 0.5, y = sqrt(3)/6, xend = 0.75, yend = sqrt(3)/4, linetype = "dotted", color = "grey70") +
  # Points: Uniform shape, Color by Dominance, Alpha by Abundance
  geom_point(aes(color = FinalColor, alpha = PlotAlpha), size = 2.5, shape = 16) +
  # Highlight Boundary proteins with a gold ring
  geom_point(data = filter(ternary_df, Niche_Status == "Boundary"), color = "gold", shape = 1, size = 3, stroke = 0.8) +
  # Labels (Labeling ALL boundary markers and TOP specific ones)
  geom_text_repel(aes(label = Label), size = 3, fontface = "bold", box.padding = 0.6, 
                  max.overlaps = 100, segment.color = "grey50") +
  # Scales
  scale_color_manual(values = ext_colors_opt) +
  scale_alpha_identity() +
  # Vertices Labels
  annotate("text", x = -0.05, y = 0, label = "Fibrin", fontface = "bold", color = spatial_colors[1], hjust=1) +
  annotate("text", x = 1.05, y = 0, label = "Collagen", fontface = "bold", color = spatial_colors[2], hjust=0) +
  annotate("text", x = 0.5, y = sqrt(3)/2 + 0.05, label = "Wall", fontface = "bold", color = spatial_colors[3]) +
  theme_void() +
  labs(title = "D14 Spatial Heterogeneity: Functional Landscapes",
       subtitle = "Dots: Individual Proteins | Opacity: Abundance Level",
       caption = "Legend: Colored regions = Niche Specificity | Gold Ring = Boundary Protein (Shared Signature)") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
        plot.caption = element_text(hjust = 0.5, color = "grey40", face = "italic", size = 9),
        legend.position = "none",
        plot.margin = margin(20, 20, 20, 20))

print(p_ternary_opt)
ggsave(file.path(output_dir, "5_Ternary_Plot_Optimized.pdf"), p_ternary_opt, width = 7, height = 7)

# ------------------------------------------------------------------------------
# 8. Boxplots of Significant Proteins (Top 6 Global)
# ------------------------------------------------------------------------------
message(">>> Generating Boxplots...")
top_6_ids <- head(rownames(res_global[order(res_global$F, decreasing = TRUE), ]), 20)

plot_list <- list()
for (id in top_6_ids) {
  gene_name <- feature_subset[id, "PG.Genes"]
  gene_name <- sapply(strsplit(gene_name, ";"), `[`, 1)
  
  df_prot <- data.frame(
    Expression = as.numeric(expr_norm[id, ]),
    Group = metadata_subset$Group
  )
  
  p_box <- ggplot(df_prot, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.8) +
    scale_fill_manual(values = spatial_colors) +
    labs(title = gene_name, x = "", y = "Log2 Intensity") +
    theme_paper +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  plot_list[[gene_name]] <- p_box
}

p_combined_box <- wrap_plots(plot_list, ncol = 5)
print(p_combined_box)
ggsave(file.path(output_dir, "6_Top_Global_Boxplots.pdf"), p_combined_box, width = 10, height = 7)


# ------------------------------------------------------------------------------
# 9. Functional Enrichment (GSEA)
# ------------------------------------------------------------------------------

# Helper 1: Comparative Side-by-Side Dotplot (Grouped by Category)
plot_comparative_dotplot <- function(gsea_res_list, title, group_up, group_down) {
  # Merge BP, CC, MF if it's a list, otherwise handle single result
  if (is.list(gsea_res_list) && !inherits(gsea_res_list, "gseaResult")) {
    df_list <- lapply(names(gsea_res_list), function(ont) {
      d <- as.data.frame(gsea_res_list[[ont]])
      if (nrow(d) > 0) d$Ontology <- ont
      return(d)
    })
    df <- do.call(rbind, df_list)
  } else {
    df <- as.data.frame(gsea_res_list)
    df$Ontology <- "Other"
  }
  
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  # Assign group names based on NES sign
  df$Enriched_In <- ifelse(df$NES > 0, group_up, group_down)
  
  # Pick top terms for each side and each ontology
  df_top <- df %>%
    group_by(Ontology, Enriched_In) %>%
    slice_max(order_by = abs(NES), n = 5) %>%
    ungroup() %>%
    unique()
  
  # Factorize Ontology to match common scientific order
  df_top$Ontology <- factor(df_top$Ontology, levels = c("MF", "CC", "BP", "Reactome", "WikiPathways", "Other"))
  # Human-readable labels for facets
  levels(df_top$Ontology) <- c("Molecular function", "Cellular component", "Biological process", "Reactome", "WikiPathways", "Other")
  
  # Factorize Group for X-axis
  df_top$Enriched_In <- factor(df_top$Enriched_In, levels = c(group_down, group_up))
  
  # Order Descriptions for vertical alignment
  df_top$Description <- factor(df_top$Description, levels = rev(unique(df_top$Description[order(df_top$Ontology, df_top$NES)])))
  
  p <- ggplot(df_top, aes(x = Enriched_In, y = Description, size = setSize, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#BE1E2D", high = "#2B3990", name = "Adj.P") +
    facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL, size = "Set Size") +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 10),
          axis.text.x = element_text(face = "bold", color = "black"),
          axis.text.y = element_text(size = 9),
          panel.grid.major.x = element_line(color = "grey90", linetype = "dotted"),
          legend.position = "right")
  
  return(p)
}

# Helper 2: Enrichment Volcano (Consolidated)
plot_enrichment_volcano <- function(results_list, title, use_pvalue = FALSE) {
  plot_df <- list()
  if (!is.null(results_list$Reactome)) {
    df <- as.data.frame(results_list$Reactome); df$Ontology <- "Reactome"; plot_df$Reactome <- df
  }
  if (!is.null(results_list$GO_BP)) {
    df <- as.data.frame(results_list$GO_BP); df$Ontology <- "GO:BP"; plot_df$GO_BP <- df
  }
  if (!is.null(results_list$Wiki)) {
    df <- as.data.frame(results_list$Wiki); df$Ontology <- "WikiPathways"; plot_df$Wiki <- df
  }
  
  if (length(plot_df) == 0) return(NULL)
  combined_df <- bind_rows(plot_df)
  
  # Adaptive Y-axis: Use pvalue for sparse results to show distribution
  if (use_pvalue) {
    combined_df$logP <- -log10(combined_df$pvalue)
    y_label <- "-Log10 p-value"
  } else {
    combined_df$logP <- -log10(combined_df$p.adjust)
    y_label <- "-Log10 q-value"
  }
  
  combined_df$logP[is.infinite(combined_df$logP)] <- max(combined_df$logP[is.finite(combined_df$logP)], na.rm=TRUE) + 1
  
  # Labeling Selection: Force top results if using p-value (to ensure informative plot)
  if (use_pvalue) {
    label_df <- combined_df %>%
      group_by(Ontology, sign(NES)) %>%
      slice_max(logP, n = 3) %>%
      ungroup() %>%
      unique()
  } else {
    label_df <- combined_df %>%
      filter(p.adjust < 0.05) %>%
      group_by(Ontology, sign(NES)) %>%
      slice_max(logP, n = 5) %>%
      ungroup() %>%
      unique()
  }
  
  label_df$Description <- ifelse(nchar(label_df$Description) > 50, paste0(substr(label_df$Description, 1, 47), "..."), label_df$Description)

  p <- ggplot(combined_df, aes(x = NES, y = logP, color = Ontology, size = setSize)) +
    theme_bw() +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("GO:BP" = "#74ADD1", "Reactome" = "#D73027", "WikiPathways" = "#91CF60")) +
    # Reference Lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey90") +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted", color = "grey80") +
    # Labels
    geom_text_repel(data = label_df, aes(label = Description), size = 2.5, max.overlaps = 30, box.padding = 0.8,
                    point.padding = 0.3, segment.color = "grey30", segment.alpha = 0.5, force = 2, show.legend = FALSE) +
    labs(title = title, x = "NES", y = y_label, size = "Set Size", color = "Ontology") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Add footnote if p-value is used
  if (use_pvalue) {
    p <- p + labs(caption = "Note: Y-axis uses p-value for enhanced visual dispersion; dashed line indicates q=0.05 threshold.")
  }
  
  return(p)
}

# Master Function: Run GSEA for a contrast
run_pairwise_gsea <- function(contrast_name, title) {
  # Safety Check: Ensure contrast exists in the model
  if (!contrast_name %in% colnames(fit2)) {
    message("    - Skip GSEA for ", contrast_name, ": Contrast not found in Limma model.")
    return(NULL)
  }
  
  # Parse group names
  groups <- strsplit(contrast_name, "_vs_")[[1]]
  group_up <- groups[1]
  group_down <- groups[2]
  
  # Prepare Ranked Gene List (Using Index to prevent subscript errors)
  target_idx <- which(colnames(fit2) == contrast_name)
  res <- topTable(fit2, coef = target_idx, number = Inf)
  
  # Ensure feature_subset is used consistently
  res$Gene <- feature_subset[rownames(res), "PG.Genes"]
  res$Gene <- sapply(strsplit(as.character(res$Gene), ";"), `[`, 1)
  res$Gene <- stringr::str_to_title(res$Gene)
  
  res <- res %>% filter(!is.na(Gene)) %>% group_by(Gene) %>% filter(abs(t) == max(abs(t))) %>% ungroup() %>% arrange(desc(t))
  gene_list <- res$t
  names(gene_list) <- res$Gene
  
  # Entrez Mapping
  gene_map <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  matched_idx <- match(gene_map$SYMBOL, names(gene_list))
  entrez_list <- gene_list[matched_idx]
  names(entrez_list) <- gene_map$ENTREZID
  entrez_list <- sort(entrez_list, decreasing = TRUE)
  
  results_list <- list()
  
  # A. Reactome
  results_list$Reactome <- gsePathway(entrez_list, organism = "mouse", pvalueCutoff = 0.2, verbose = FALSE)
  
  # B. GO (BP, CC, MF)
  message("    - Running GO BP...")
  gsea_bp <- gseGO(entrez_list, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.2, verbose = FALSE)
  if (!is.null(gsea_bp)) results_list$GO_BP <- simplify(gsea_bp)
  
  message("    - Running GO CC...")
  gsea_cc <- gseGO(entrez_list, OrgDb = org.Mm.eg.db, ont = "CC", pvalueCutoff = 0.2, verbose = FALSE)
  if (!is.null(gsea_cc)) results_list$GO_CC <- simplify(gsea_cc)
  
  message("    - Running GO MF...")
  gsea_mf <- gseGO(entrez_list, OrgDb = org.Mm.eg.db, ont = "MF", pvalueCutoff = 0.2, verbose = FALSE)
  if (!is.null(gsea_mf)) results_list$GO_MF <- simplify(gsea_mf)
  
  # C. WikiPathways
  results_list$Wiki <- tryCatch({
    gseWP(entrez_list, organism = "Mus musculus", pvalueCutoff = 0.2, verbose = FALSE)
  }, error = function(e) NULL)
  
  # Create Plots
  plots <- list()
  
  # 1. Combined GO Dotplot (Side-by-Side as per reference)
  go_list <- list(BP = results_list$GO_BP, CC = results_list$GO_CC, MF = results_list$GO_MF)
  plots$GO_Comparative <- plot_comparative_dotplot(go_list, paste(title, "- GO Categories"), group_up, group_down)
  
  # 2. Reactome Comparative
  plots$Reactome <- plot_comparative_dotplot(results_list$Reactome, paste(title, "- Reactome"), group_up, group_down)
  
  # 3. Wiki Comparative
  plots$Wiki <- plot_comparative_dotplot(results_list$Wiki, paste(title, "- WikiPathways"), group_up, group_down)
  
  # 4. Enrichment Volcano
  # Boost aesthetics for Wall vs Collagen (W_vs_C) as requested
  use_boost <- (contrast_name == "W_vs_C")
  plots$Volcano <- plot_enrichment_volcano(results_list, paste(title, "- Global Enrichment"), use_pvalue = use_boost)
  
  return(plots)
}

# Execute Analysis
message(">>> Starting Pairwise GSEA Comparisons...")
p_gsea_CF <- run_pairwise_gsea("C_vs_F", "D14 Collagen vs Fibrin")
p_gsea_WC <- run_pairwise_gsea("W_vs_C", "D14 Wall vs Collagen")
p_gsea_WF   <- run_pairwise_gsea("W_vs_F", "D14 Wall vs Fibrin")

# Execution for Requested Trajectories
message(">>> Running Requested Trajectory GSEA...")
conts_to_run <- c("D7F_vs_D2F", "D7C_vs_D7F", "D14C_vs_D7C", "D14W_vs_D14C", "D14F_vs_D7F", 
                  "D14C_vs_D14F", "D7C_vs_D2F", "D14C_vs_D7F", "D14W_vs_D14F", "D14F_vs_D2F", "D14C_vs_D2F")
for (cc in conts_to_run) {
  p_temp <- run_pairwise_gsea(cc, cc)
  if (!is.null(p_temp)) {
    save_gsea_results(p_temp, paste0("GSEA_", cc, ".pdf"))
  }
}

# Helper to save GSEA plots
save_gsea_results <- function(p_list, filename) {
  if (length(p_list) > 0) {
    v <- p_list$Volcano
    d <- wrap_plots(p_list[names(p_list) != "Volcano"], ncol = 1)
    g <- (v / d) + plot_layout(heights = c(1, 1))
    ggsave(file.path(output_dir, filename), g, width = 10, height = 15)
  }
}

# Save initial results
save_gsea_results(p_gsea_CF, "7_GSEA_Collagen_vs_Fibrin.pdf")
save_gsea_results(p_gsea_WC, "8_GSEA_Wall_vs_Collagen.pdf")
save_gsea_results(p_gsea_WF, "9_GSEA_Wall_vs_Fibrin.pdf")

# ------------------------------------------------------------------------------
# 8. Unified IPA Data Export (Multi-sheet Excel)
# ------------------------------------------------------------------------------
message("\n>>> Exporting Unified IPA Excel Data (all 9 requested contrasts)...")

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

ipa_wb <- createWorkbook()

# Function to extract and format for IPA (All Proteins)
extract_ipa_sheet <- function(contrast_name) {
  if (!contrast_name %in% colnames(fit2)) return(NULL)
  
  res <- topTable(fit2, coef = contrast_name, number = Inf)
  
  # Format as per User Request: ID, GeneName, logFC, P
  ipa_df <- data.frame(
    ProteinID = rownames(res),
    Symbol = feature_subset$Symbol[match(rownames(res), rownames(feature_subset))],
    logFC = res$logFC,
    P_value = res$P.Value,
    P_value_adjust = res$adj.P.Val,
    check.names = FALSE
  )
  return(ipa_df)
}

# Add sheets for all trajectories (Displaying as Hyphenated labels for IPA, e.g., D2F-D7F)
# Also export individual Excel files for each trajectory
for (cc in conts_to_run) {
    sheet_data <- extract_ipa_sheet(cc)
    if (!is.null(sheet_data)) {
        # 1. Add to Master Workbook
        sheet_label <- gsub("([A-Z0-9]+)_vs_([A-Z0-9]+)", "\\2-\\1", cc)
        addWorksheet(ipa_wb, sheet_label)
        writeData(ipa_wb, sheet_label, sheet_data, startCol = 1, startRow = 1)
        setColWidths(ipa_wb, sheet_label, cols = 1:5, widths = "auto")
        
        # 2. Export Standalone Excel File
        standalone_wb <- createWorkbook()
        addWorksheet(standalone_wb, sheet_label)
        writeData(standalone_wb, sheet_label, sheet_data)
        setColWidths(standalone_wb, sheet_label, cols = 1:5, widths = "auto")
        saveWorkbook(standalone_wb, file.path(output_dir, paste0("IPA_Trajectory_", sheet_label, ".xlsx")), overwrite = TRUE)
    }
}

saveWorkbook(ipa_wb, file.path(output_dir, "IPA_Master_Export_All_Proteins.xlsx"), overwrite = TRUE)
message("    - Completed: 11 Individual Excel exports and Master Export (All Proteins)")

# ------------------------------------------------------------------------------
# 8. IPA Data Export (All Contrasts)
# ------------------------------------------------------------------------------
message("\n>>> Exporting IPA Data for all contrasts...")

export_ipa_data <- function(contrast_name, filename) {
  if (!contrast_name %in% colnames(fit2)) return(NULL)
  
  # Extract DE results
  res <- topTable(fit2, coef = contrast_name, number = Inf)
  
  # Join with Feature Data metadata
  res$Symbol <- feature_subset[rownames(res), "Symbol"]
  res$ProteinID <- rownames(res)
  res$Description <- feature_subset[rownames(res), "PG.ProteinDescriptions"]
  
  # Compile dataframe for IPA consumption
  ipa_df <- res[, c("Symbol", "logFC", "P.Value", "adj.P.Val", "ProteinID", "Description")]
  
  write.csv(ipa_df, file.path(output_dir, filename), row.names = FALSE)
  message("    - Exported IPA data: ", filename)
}

# 1. Spatial Contrasts (D14)
export_ipa_data("C_vs_F", "IPA_D14_Collagen_vs_Fibrin.csv")
export_ipa_data("W_vs_C", "IPA_D14_Wall_vs_Collagen.csv")
export_ipa_data("W_vs_F", "IPA_D14_Wall_vs_Fibrin.csv")

# 2. Longitudinal/Temporal Contrasts
export_ipa_data("D7C_vs_D2F", "IPA_D7C_vs_D2F.csv")   # Was D2F-D7C
export_ipa_data("D7C_vs_D7F", "IPA_D7C_vs_D7F.csv")   # Was D7F-D7C
export_ipa_data("D7F_vs_D2F", "IPA_D7F_vs_D2F.csv")   # Was D2F-D7F
export_ipa_data("D14C_vs_D7C", "IPA_D14C_vs_D7C.csv") # Was D7C-D14C

# ------------------------------------------------------------------------------
# 8. PPI Signaling Networks (pathlinkR)
# ------------------------------------------------------------------------------
message("\n>>> Starting pathlinkR Signaling Analysis...")

# Helper to process pathlinkR for a given comparison
run_pathlinkr_niche <- function(contrast_name, title, output_fname) {
  if (!contrast_name %in% colnames(fit2)) return(NULL)
  
  res <- topTable(fit2, coef = contrast_name, number = Inf)
  # Ensure Symbol is character and not NULL
  if (!"Symbol" %in% colnames(res)) {
     res$Symbol <- feature_subset$Symbol[match(rownames(res), rownames(feature_subset))]
  }
  
  # Filter for significantly regulated genes (p < 0.05, |logFC| > 0.5)
  sig_res <- res[res$P.Value < 0.05 & abs(res$logFC) > 0.5, ]
  if (nrow(sig_res) < 5) {
    message("    - Notice: Too few significant genes for ", title, " network.")
    return(NULL)
  }
  
  # Use top 200 genes for a readable network if many are significant
  sig_res <- head(sig_res[order(sig_res$P.Value), ], 200)
  
  # 1. Map Mouse Symbols to Human ENSEMBL IDs (Required by pathlinkR 1.6.0)
  if (!exists("mappingFile")) {
    data("mappingFile", package = "pathlinkR", envir = environment())
    m_file <<- as.data.frame(get("mappingFile"))
  }
  
  # Ensure sig_res has Symbol
  if (!"Symbol" %in% colnames(sig_res) || is.null(sig_res$Symbol)) {
      message("    - Warning: Symbol column missing in sig_res for ", title)
      return(NULL)
  }
  
  sig_res$HumanSymbol <- toupper(as.character(sig_res$Symbol))
  
  # Filter out any NAs from Symbol
  sig_res <- sig_res[!is.na(sig_res$HumanSymbol), ]
  
  mapping_bridge <- merge(sig_res, m_file, by.x = "HumanSymbol", by.y = "hgncSymbol", all.x = TRUE)
  mapping_bridge <- mapping_bridge[!is.na(mapping_bridge$ensemblGeneId), ]
  mapping_bridge <- mapping_bridge[!duplicated(mapping_bridge$ensemblGeneId), ]
  
  if (nrow(mapping_bridge) == 0) {
    message("    - Warning: No genes from ", title, " mapped to Human ENSEMBL IDs.")
    return(NULL)
  }
  
  # Prepare pathlinkR input
  res_mock <- data.frame(
    log2FoldChange = mapping_bridge$logFC,
    padj = mapping_bridge$adj.P.Val,
    hgncSymbol = as.character(mapping_bridge$Symbol),
    stringsAsFactors = FALSE
  )
  rownames(res_mock) <- as.character(mapping_bridge$ensemblGeneId)
  
  message("    - Building PPI Network for ", title, " (", nrow(res_mock), " genes)...")
  
  tryCatch({
    # 2. Build PPI Network (InnateDB)
    exNetwork <- ppiBuildNetwork(
      rnaseqResult = res_mock,
      filterInput = FALSE,
      columnFC = "log2FoldChange",
      columnP = "padj",
      order = "zero"
    )
    
    if (!is.null(exNetwork)) {
        # 3. Label Injection (Fix for pathlinkR 1.6.0 label mapping)
        if (requireNamespace("igraph", quietly = TRUE)) {
          node_ids <- igraph::V(exNetwork)$name
          symbol_map <- res_mock$hgncSymbol[match(node_ids, rownames(res_mock))]
          igraph::vertex_attr(exNetwork, "hgncSymbol") <- as.character(symbol_map)
        }

        # 4. Plot Network
        p_net <- ppiPlotNetwork(
          network = exNetwork,
          title = paste0("Niche Signaling: ", title),
          fillColumn = log2FoldChange,
          fillType = "foldChange",
          label = TRUE,
          labelColumn = hgncSymbol,
          legend = TRUE
        )
        
        if (!is.null(p_net)) {
          print(p_net)
          ggsave(file.path(output_dir, output_fname), p_net, width = 10, height = 8)
          message("    - Exported: ", output_fname)
        }
    }
  }, error = function(e) {
    message("    - Warning: pathlinkR failed for ", title, ": ", e$message)
  })
}

# Run Signaling Analysis for all 3 spatial contrasts
run_pathlinkr_niche("C_vs_F", "Collagen vs Fibrin Niche", "10_pathlinkR_Collagen_vs_Fibrin.pdf")
run_pathlinkr_niche("W_vs_C", "Wall vs Collagen Niche", "11_pathlinkR_Wall_vs_Collagen.pdf")
run_pathlinkr_niche("W_vs_F", "Wall vs Fibrin Niche", "12_pathlinkR_Wall_vs_Fibrin.pdf")

# ------------------------------------------------------------------------------
# 11. Deep Matrisome Characterization
# ------------------------------------------------------------------------------
message("\n>>> Starting Deep Matrisome Characterization...")

if (exists("matrisome_detected_genes") && length(matrisome_detected_genes) > 0) {
  # 1. Map Proteins to Categories (Using row indices from feature_filtered)
  matrisome_rows <- which(feature_filtered$UniqueSymbol %in% matrisome_detected_genes)
  mat_df <- data.frame(
    ProteinID = rownames(feature_filtered)[matrisome_rows],
    UniqueSymbol = feature_filtered$UniqueSymbol[matrisome_rows],
    stringsAsFactors = FALSE
  )
  mat_df$Category <- "Unknown"
  
  for (cat in names(matrisome_full_db)) {
    cat_genes <- get_genes_from_ids(matrisome_full_db[[cat]])
    mat_df$Category[mat_df$UniqueSymbol %in% cat_genes] <- cat
  }
  
  # 2. Composition Heatmap (Mean Expression by Niche)
  # Subset expression using ProteinID
  expr_mat <- expr_norm[mat_df$ProteinID, ]
  rownames(expr_mat) <- mat_df$UniqueSymbol # For better labeling
  
  # Calculate per-niche means
  n_groups <- levels(metadata_subset$Group)
  niche_means <- matrix(NA, nrow = nrow(expr_mat), ncol = length(n_groups))
  rownames(niche_means) <- rownames(expr_mat)
  colnames(niche_means) <- n_groups
  
  for (i in 1:length(n_groups)) {
    g_mask <- metadata_subset$Group == n_groups[i]
    if (sum(g_mask) > 0) {
      niche_means[, i] <- rowMeans(expr_mat[, g_mask, drop = FALSE], na.rm = TRUE)
    }
  }
  
  # Remove any rows with all NAs (safety)
  valid_rows <- rowSums(!is.na(niche_means)) > 0
  niche_means <- niche_means[valid_rows, , drop = FALSE]
  mat_df_sub <- mat_df[valid_rows, ]
  
  if (nrow(niche_means) > 1) {
    # Z-score scaling for heatmap
    niche_z <- t(apply(niche_means, 1, scale))
    colnames(niche_z) <- colnames(niche_means)
    
    # Plot Heatmap
    library(pheatmap)
    mat_annotation <- data.frame(Category = mat_df_sub$Category, row.names = mat_df_sub$UniqueSymbol)
    
    # Sort by Category for better visualization
    sort_idx <- order(mat_df_sub$Category)
    
    p_mat_heat <- pheatmap(niche_z[sort_idx, ], 
                          annotation_row = mat_annotation,
                          cluster_rows = FALSE,
                          cluster_cols = FALSE,
                          show_colnames = TRUE,
                          main = "Matrisome Composition by Spatial Niche",
                          color = colorRampPalette(c("#2B3990", "white", "#BE1E2D"))(100),
                          silent = TRUE)
    
    save_pheatmap_pdf <- function(x, filename, width=7, height=10) {
      pdf(filename, width=width, height=height)
      grid::grid.newpage()
      grid::grid.draw(x$gtable)
      dev.off()
    }
    save_pheatmap_pdf(p_mat_heat, file.path(output_dir, "13_Matrisome_Composition_Heatmap.pdf"))
  }
  
  # 3. Categorical Distribution Pie Chart (Optimized based on Reference)
  cat_counts <- as.data.frame(table(mat_df$Category))
  colnames(cat_counts) <- c("Category", "Count")
  cat_counts <- cat_counts %>%
    mutate(Percentage = Count / sum(Count) * 100,
           Label = paste0(Category, " (", round(Percentage, 0), "%)")) %>%
    arrange(desc(Percentage)) # Sort by size for cleaner pie
  
  # Professional palette similar to reference (muted, high-quality)
  mat_palette <- c(
    "#54AD8F", # Green (Large)
    "#9DA6C9", # Soft Blue
    "#CDCDCD", # Grey
    "#D8B38D", # Sand
    "#E6C15C", # Gold
    "#A8C97F", # Soft Green
    "#D589AD", # Pink
    "#8B99C1", # Steel Blue
    "#E89275"  # Muted Orange
  )
  
  p_pie <- ggplot(cat_counts, aes(x = "", y = Count, fill = reorder(Category, -Percentage))) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0.3) +
    coord_polar("y", start = 0) +
    # Use geom_text_repel for cleaner labels if many categories, or simple geom_text
    geom_text(aes(label = Label, x = 1.3), 
              position = position_stack(vjust = 0.5), 
              size = 3.5, fontface = "bold") +
    scale_fill_manual(values = mat_palette) +
    theme_void() +
    labs(title = "Matrisome Categorical Distribution (N=80)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "none")
  
  if (nrow(cat_counts) > 0) {
    ggsave(file.path(output_dir, "14_Matrisome_Categories_Pie.pdf"), p_pie, width = 7, height = 7)
    print(p_pie)
  }
}
# ------------------------------------------------------------------------------
# 12. Spatial Abundance Profiling (6 Groups)
# ------------------------------------------------------------------------------
message("\n>>> Generating 6-Group Spatial Abundance Profiling...")

# 1. Define Groups & Order (consistent with dataset architecture)
all_target_groups <- c("D2_Fibrin", "D7_Fibrin", "D7_Collagen", "D14_Fibrin", "D14_Collagen", "D14_Wall")
all_mask <- metadata$Group %in% all_target_groups
expr_6g <- expr_norm_global[, all_mask]
meta_6g <- metadata[all_mask, ]
meta_6g$Group <- factor(meta_6g$Group, levels = all_target_groups)

# 2. Target Markers (Robust Matching)
markers <- c("S100a4", "Tgfb1", "Ogn", "Col1a1", "Col4a2", "Col3a1", "Lox", "Serpinh1", "P4ha2", "Bgn", "Cspg4", "Lgals1", "Lgals3", "Pf4;Pf4v1", "Fbln7")

# Map markers to UniqueSymbols in dataset (Case-insensitive & Partial matching)
marker_map <- list()
for (m in markers) {
  # Try case-insensitive match against both Symbol and PG.Genes/PG.ProteinGroups
  idx <- which(toupper(feature_filtered$Symbol) == toupper(m))[1]
  if (is.na(idx)) {
    # Try regex matching on both PG.Genes and PG.ProteinGroups
    # This catches things like "Pf4" in "Pf4;Pf4v1" or vice versa
    search_term <- gsub(";.*", "", m) # Use base gene name before semicolon
    idx <- grep(paste0("(^|;)", search_term, "($|;)"), feature_filtered$PG.Genes, ignore.case = TRUE)[1]
    if (is.na(idx)) {
        idx <- grep(paste0("(^|;)", search_term, "($|;)"), feature_filtered$Symbol, ignore.case = TRUE)[1]
    }
  }
  
  if (!is.na(idx)) {
    # Store the actual Rowname (ProteinGroup ID) for indexing
    marker_map[[m]] <- rownames(feature_filtered)[idx]
  } else {
    message("    - Note: Marker ", m, " not found in dataset.")
  }
}

# 3. Extract Data (Using IMPUTED Global Data: expr_norm_global)
plot_data_list <- list()
for (m in names(marker_map)) {
  id_prot <- marker_map[[m]]
  # id_prot is the ProteinGroup ID (vignette: rownames of the matrix)
  df_m <- data.frame(
    Expression = as.numeric(expr_norm_global[id_prot, all_mask]),
    Group = meta_6g$Group,
    Gene = m
  )
  plot_data_list[[m]] <- df_m
}
plot_df <- do.call(rbind, plot_data_list)

# Summary Stats
summary_df <- plot_df %>%
  group_by(Gene, Group) %>%
  summarize(
    Mean = mean(Expression, na.rm = TRUE),
    SEM = sd(Expression, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

# 4. Plot Line Profiles
# Use niche colors for grouping
group_colors <- c(
  "D2_Fibrin" = "#B03A2E", "D7_Fibrin" = "#B03A2E", "D14_Fibrin" = "#B03A2E",
  "D7_Collagen" = "#2E86C1", "D14_Collagen" = "#2E86C1",
  "D14_Wall" = "#8E44AD"
)

p_lines <- ggplot(summary_df, aes(x = Group, y = Mean, group = Gene, color = Group)) +
  geom_line(color = "grey60", size = 0.5) + # Connector line
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
  geom_point(size = 2.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 5) +
  scale_color_manual(values = group_colors) +
  theme_elegant() +
  labs(title = "Spatial and Temporal Evolution of Key Markers",
       x = NULL, y = "Log2 Normalized Intensity",
       caption = "Groups: D2-D14 Fibrin, D7-D14 Collagen, D14 Wall") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"),
        legend.position = "none")

print(p_lines)
ggsave(file.path(output_dir, "15_Key_Markers_6Group_Lines.pdf"), p_lines, width = 12, height = 10)

message("    - 6-Group profiles exported: 15_Key_Markers_6Group_Lines.pdf")

# ------------------------------------------------------------------------------
# 10. Statistical Methods Summary
# ------------------------------------------------------------------------------
sink(file.path(output_dir, "Statistical_Methods.txt"))
cat("Statistical Analytical Methods:\n\n")
cat("1. Data Preprocessing:\n")
cat("   - Log2 transformation\n")
cat("   - Missing value imputation: MinProb method (Mean shifted by -1.8sd, width 0.3sd)\n")
cat("   - Normalization: Quantile normalization\n\n")
cat("2. Differential Expression Analysis (Limma):\n")
cat("   - Linear model fitted for spatial-temporal subgroups (D2-D14 Fibrin, D7-D14 Collagen, D14 Wall)\n")
cat("   - Pairwise contrasts: C vs F, W vs C, W vs F evaluated using empirical Bayes t-tests.\n")
cat("   - Global significance evaluated using F-test across all groups.\n")
cat("   - Significance threshold: adj.P.Value < 0.05 and |log2FC| > 0.58 (1.5-fold).\n\n")
cat("3. Functional Enrichment (GSEA):\n")
cat("   - Ranked list metric: t-statistic from Limma.\n")
cat("   - Methods: gsePathway (Reactome), gseGO (GO Biological Process).\n")
cat("   - P-value cutoff: 0.2 (exploratory) with BH correction.\n")
sink()

message("\n>>> D14 Spatial Analysis Completed. Check output folder: ", output_dir)
