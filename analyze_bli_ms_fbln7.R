# ==============================================================================
# BLI-MS (Biolayer Interferometry - Mass Spectrometry) Data Analysis Script
# Project: FBLN7 (Fibulin-7) Molecular Fishing & Background Protein Filtering
# Reference Paper: BLI-MS Technology (DOI: 10.1002/pmic.202100031)
# Date: 2026-08-31 (Updated with Medium-Confidence Labels & Robust Barplot)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Environment Setup & Package Auto-Installation
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "dplyr", "tidyr", "readr", "pheatmap", "ggrepel", "scales")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing required package:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# 2. File Path Configuration & Data Loading
# ------------------------------------------------------------------------------
input_file <- "蛋白鉴定列表.csv"
if (!file.exists(input_file)) {
  input_file <- "c:/Users/SimonYao/Desktop/分子垂钓/蛋白鉴定列表.csv"
}

cat("=== Loading Dataset:", input_file, "===\n")

df_raw <- tryCatch({
  read.csv(input_file, fileEncoding = "UTF-8", stringsAsFactors = FALSE, check.names = FALSE)
}, error = function(e) {
  read.csv(input_file, fileEncoding = "GBK", stringsAsFactors = FALSE, check.names = FALSE)
})

cat("Loaded:", nrow(df_raw), "proteins x", ncol(df_raw), "columns.\n")

# ------------------------------------------------------------------------------
# 3. Step 1 — Remove Contaminants, Reverse Hits, and IDs-only Entries
# ------------------------------------------------------------------------------
n_before <- nrow(df_raw)

# Remove reverse decoy sequences
if ("Reverse" %in% names(df_raw)) {
  df_raw <- df_raw[df_raw[["Reverse"]] != "+", ]
}
# Remove common contaminants (keratins, trypsin, BSA, etc.)
if ("Potential contaminant" %in% names(df_raw)) {
  df_raw <- df_raw[df_raw[["Potential contaminant"]] != "+", ]
}
# Remove entries without gene/protein names
df_raw <- df_raw[!(is.na(df_raw[["Protein IDs"]]) | df_raw[["Protein IDs"]] == ""), ]

cat(sprintf("After contaminant/reverse filtering: %d -> %d proteins.\n", n_before, nrow(df_raw)))

# ------------------------------------------------------------------------------
# 4. Define Replicates & Metadata Columns
# ------------------------------------------------------------------------------
fishing_cols     <- c("A5", "B5", "C5", "D5", "E5")
blank_cols       <- c("F5", "G5", "H5")
uniq_pep_fishing <- c("Unique peptides A5", "Unique peptides B5", "Unique peptides C5",
                      "Unique peptides D5", "Unique peptides E5")

# Convert intensity columns to numeric
for (col in c(fishing_cols, blank_cols)) {
  df_raw[[col]] <- suppressWarnings(as.numeric(as.character(df_raw[[col]])))
  df_raw[[col]][is.na(df_raw[[col]])] <- 0
}

# Convert unique peptide count columns to numeric
for (col in uniq_pep_fishing) {
  if (col %in% names(df_raw)) {
    df_raw[[col]] <- suppressWarnings(as.numeric(as.character(df_raw[[col]])))
    df_raw[[col]][is.na(df_raw[[col]])] <- 0
  }
}

# ------------------------------------------------------------------------------
# 5. Background Filtering & Statistical Calculations
# ------------------------------------------------------------------------------
cat("=== Performing Background Protein Subtraction & Statistical Testing ===\n")

df_analysis <- df_raw %>%
  mutate(
    n_fishing = rowSums(select(., all_of(fishing_cols)) > 0),
    n_blank   = rowSums(select(., all_of(blank_cols)) > 0),
    mean_fishing_raw = rowMeans(select(., all_of(fishing_cols))),
    mean_blank_raw   = rowMeans(select(., all_of(blank_cols)))
  )

# Calculate max unique peptides across fishing replicates
upep_cols_present <- uniq_pep_fishing[uniq_pep_fishing %in% names(df_analysis)]
if (length(upep_cols_present) > 0) {
  df_analysis$max_unique_pep_fishing <- apply(df_analysis[, upep_cols_present], 1, max, na.rm = TRUE)
} else {
  df_analysis$max_unique_pep_fishing <- 1
  warning("Unique peptide columns not found. Unique-peptide filtering disabled.")
}

# LOD imputation for zeros (0.5 × global minimum detected intensity)
all_intensities <- unlist(df_analysis[, c(fishing_cols, blank_cols)])
min_detected    <- min(all_intensities[all_intensities > 0], na.rm = TRUE)
lod_value       <- min_detected * 0.5
cat(sprintf("Minimum Detected Intensity: %.2e | LOD Imputation Value: %.2e\n", min_detected, lod_value))

imputed_fishing <- df_analysis[, fishing_cols]
imputed_fishing[imputed_fishing == 0] <- lod_value

imputed_blank <- df_analysis[, blank_cols]
imputed_blank[imputed_blank == 0] <- lod_value

# Log2 transformation
log2_fishing <- log2(imputed_fishing)
log2_blank   <- log2(imputed_blank)

# Welch's t-test (unequal variance)
p_values <- apply(cbind(log2_fishing, log2_blank), 1, function(row) {
  f_vals <- row[1:length(fishing_cols)]
  b_vals <- row[(length(fishing_cols) + 1):length(row)]
  if (sd(f_vals) == 0 && sd(b_vals) == 0) return(1.0)
  test <- tryCatch(t.test(f_vals, b_vals, var.equal = FALSE), error = function(e) NULL)
  if (is.null(test)) return(1.0) else return(test$p.value)
})

df_analysis$log2_fishing_mean <- rowMeans(log2_fishing)
df_analysis$log2_blank_mean   <- rowMeans(log2_blank)
df_analysis$log2_FC           <- df_analysis$log2_fishing_mean - df_analysis$log2_blank_mean
df_analysis$Fold_Change       <- 2 ^ df_analysis$log2_FC

df_analysis$pvalue <- p_values
df_analysis$fdr    <- p.adjust(p_values, method = "BH")

# Net Specificity Score S = (Fishing - Blank) / Fishing
df_analysis$Specificity_Score <- ifelse(
  df_analysis$mean_fishing_raw > 0,
  (df_analysis$mean_fishing_raw - df_analysis$mean_blank_raw) / df_analysis$mean_fishing_raw,
  0
)

# Detection label for volcano plot
df_analysis$label_detected <- paste0(
  ifelse(!is.na(df_analysis[["Gene names"]]) & df_analysis[["Gene names"]] != "", 
         df_analysis[["Gene names"]], df_analysis[["Protein IDs"]]),
  " (", df_analysis$n_fishing, "/5)"
)

# ------------------------------------------------------------------------------
# 6. Multi-Tier Classification
# ------------------------------------------------------------------------------
df_analysis <- df_analysis %>%
  mutate(
    Is_FBLN7 = grepl("FBLN7|Fibulin-7", `Gene names`, ignore.case = TRUE) |
               grepl("FBLN7|Fibulin-7", `Protein names`, ignore.case = TRUE),

    Category = case_when(
      Is_FBLN7 ~ "Bait Protein (FBLN7)",
      log2_FC >= 2.0 & pvalue < 0.05 & n_fishing >= 3 & max_unique_pep_fishing >= 2 ~ "High-Confidence Interactor",
      log2_FC >= 1.0 & pvalue < 0.05 & n_fishing >= 2 ~ "Medium-Confidence Interactor",
      mean_blank_raw > mean_fishing_raw ~ "High Background Protein",
      TRUE ~ "Non-Specific / Noise"
    )
  )

# Print Summary
cat("\n=== Protein Classification Summary ===\n")
print(table(df_analysis$Category))

# ------------------------------------------------------------------------------
# 7. Build Output Table (Include GO Annotations)
# ------------------------------------------------------------------------------
go_cols_present <- intersect(
  c("Gene Ontology (biological process)",
    "Gene Ontology (cellular component)",
    "Gene Ontology (molecular function)"),
  names(df_analysis)
)

output_cols <- c("Protein IDs", "Protein names", "Gene names",
                 "log2_FC", "Fold_Change", "pvalue", "fdr",
                 "n_fishing", "n_blank", "max_unique_pep_fishing",
                 "mean_fishing_raw", "mean_blank_raw", "Specificity_Score",
                 "Category", go_cols_present)

high_conf_interactors <- df_analysis %>%
  filter(Category %in% c("Bait Protein (FBLN7)", "High-Confidence Interactor")) %>%
  arrange(desc(log2_FC))

all_interactors <- df_analysis %>%
  filter(Category %in% c("Bait Protein (FBLN7)", "High-Confidence Interactor", "Medium-Confidence Interactor")) %>%
  arrange(desc(log2_FC))

cat("\nTop 15 Specific FBLN7 Interactors:\n")
print(all_interactors %>%
        select(any_of(c("Protein IDs", "Gene names", "Protein names",
                        "log2_FC", "Fold_Change", "pvalue", "Category"))) %>%
        head(15))

# Export CSVs
write.csv(df_analysis[, output_cols], "BLI_MS_FBLN7_all_proteins_annotated.csv", row.names = FALSE)
write.csv(all_interactors[, output_cols], "BLI_MS_FBLN7_interactors_filtered.csv", row.names = FALSE)
cat("\nExported annotated CSV files.\n")

# ------------------------------------------------------------------------------
# 8. Data Visualization
# ------------------------------------------------------------------------------
cat("=== Generating Visualizations ===\n")

# --- Plot 1: Volcano Plot (Labels include Bait, High-Confidence, AND Medium-Confidence Interactors) ---
volcano_plot <- ggplot(df_analysis, aes(x = log2_FC, y = -log10(pvalue), color = Category)) +
  geom_point(alpha = 0.75, size = 2.2) +
  geom_vline(xintercept = c(1, 2), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c(
    "Bait Protein (FBLN7)"          = "#D95F02",
    "High-Confidence Interactor"    = "#E41A1C",
    "Medium-Confidence Interactor"  = "#377EB8",
    "High Background Protein"       = "#999999",
    "Non-Specific / Noise"          = "#D9D9D9"
  )) +
  # Now labeling Bait, High-Confidence, AND Medium-Confidence interactors!
  geom_text_repel(
    data = subset(df_analysis, Category %in% c("Bait Protein (FBLN7)", "High-Confidence Interactor", "Medium-Confidence Interactor")),
    aes(label = label_detected),
    size = 2.8, max.overlaps = 40, box.padding = 0.35, point.padding = 0.2, segment.color = "gray60"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "BLI-MS Volcano Plot: FBLN7 Fishing vs Blank Control",
    subtitle = "Labels show Gene Name + Detection Frequency (n_fishing / 5 replicates) for High & Medium confidence interactors",
    x = expression(log[2] ~ "(Fold Change: Fishing / Blank)"),
    y = expression(-log[10] ~ "(p-value)"),
    color = "Protein Group"
  ) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("Volcano_Plot_FBLN7_BLI_MS.png", volcano_plot, width = 11, height = 7.5, dpi = 300)
ggsave("Volcano_Plot_FBLN7_BLI_MS.pdf", volcano_plot, width = 11, height = 7.5)

# --- Plot 2: Top 20 Interactors Bar Plot (Pulls from all High & Medium confidence interactors) ---
bar_candidates <- df_analysis %>%
  filter(Category %in% c("High-Confidence Interactor", "Medium-Confidence Interactor") | (log2_FC >= 1.0 & pvalue < 0.05 & !Is_FBLN7)) %>%
  filter(!Is_FBLN7) %>%
  arrange(desc(log2_FC))

if (nrow(bar_candidates) == 0) {
  bar_candidates <- df_analysis %>% filter(!Is_FBLN7) %>% arrange(desc(log2_FC))
}

top20_df <- bar_candidates %>% head(20) %>%
  mutate(Display_Name = ifelse(!is.na(`Gene names`) & `Gene names` != "", `Gene names`, `Protein IDs`))

barplot_top20 <- ggplot(top20_df, aes(x = reorder(Display_Name, log2_FC), y = log2_FC, fill = Category)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c(
    "High-Confidence Interactor"   = "#E41A1C",
    "Medium-Confidence Interactor" = "#377EB8",
    "Non-Specific / Noise"         = "#4DAF4A"
  )) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Top 20 FBLN7 Interactors Enriched Over Blank",
    subtitle = "Ranked by Log2 Fold Change (Fishing / Blank)",
    x = "Gene / Protein Name",
    y = expression(log[2] ~ "(Fold Enrichment over Blank)"),
    fill = "Confidence Class"
  ) +
  theme(panel.grid.minor = element_blank(), legend.position = "right")

ggsave("Barplot_Top20_Interactors.png", barplot_top20, width = 8.5, height = 6.5, dpi = 300)
ggsave("Barplot_Top20_Interactors.pdf", barplot_top20, width = 8.5, height = 6.5)

# --- Plot 3: Heatmap of Top Interactors ---
if (nrow(all_interactors) >= 2) {
  top_heatmap_proteins <- all_interactors %>% head(30)
} else {
  top_heatmap_proteins <- df_analysis %>% arrange(desc(log2_FC)) %>% head(30)
}

heatmap_matrix <- as.matrix(top_heatmap_proteins[, c(fishing_cols, blank_cols)])
rownames(heatmap_matrix) <- ifelse(
  !is.na(top_heatmap_proteins[["Gene names"]]) & top_heatmap_proteins[["Gene names"]] != "",
  top_heatmap_proteins[["Gene names"]],
  top_heatmap_proteins[["Protein IDs"]]
)

heatmap_matrix_z <- t(scale(t(log2(heatmap_matrix + lod_value))))
heatmap_matrix_z[is.na(heatmap_matrix_z)] <- 0

annotation_col <- data.frame(Group = factor(c(rep("FBLN7 Fishing", 5), rep("Blank Control", 3))))
rownames(annotation_col) <- c(fishing_cols, blank_cols)

pheatmap(
  heatmap_matrix_z,
  annotation_col = annotation_col,
  main = "Heatmap: Top FBLN7 Interactors (Z-score Log2 Intensity)",
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
  fontsize_row = 8,
  cluster_rows = (nrow(heatmap_matrix_z) > 1),
  cluster_cols = TRUE,
  filename = "Heatmap_Top_Interactors.png",
  width = 7.5, height = 9.5
)

# --- Plot 4: PCA for Sample QC ---
pca_matrix <- cbind(as.matrix(log2_fishing), as.matrix(log2_blank))
pca_data   <- t(pca_matrix)

var_per_col <- apply(pca_data, 2, var)
pca_data_filtered <- pca_data[, var_per_col > 0, drop = FALSE]

pca_res  <- prcomp(pca_data_filtered, scale. = TRUE)
pca_df   <- data.frame(
  Sample = c(fishing_cols, blank_cols),
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  Group  = factor(c(rep("FBLN7 Fishing (A5-E5)", 5), rep("Blank Control (F5-H5)", 3)))
)
percent_var <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 4) +
  scale_color_manual(values = c("FBLN7 Fishing (A5-E5)" = "#E41A1C", "Blank Control (F5-H5)" = "#377EB8")) +
  theme_minimal(base_size = 13) +
  labs(
    title = "PCA of BLI-MS Samples (Quality Control)",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  )

ggsave("PCA_Sample_Clustering.png", pca_plot, width = 7, height = 5, dpi = 300)

cat("\n=== BLI-MS Analysis Completed Successfully! ===\n")
cat("Output Files:\n")
cat("  Tables:  BLI_MS_FBLN7_interactors_filtered.csv\n")
cat("           BLI_MS_FBLN7_all_proteins_annotated.csv\n")
cat("  Figures: Volcano_Plot_FBLN7_BLI_MS.png/.pdf\n")
cat("           Barplot_Top20_Interactors.png/.pdf\n")
cat("           Heatmap_Top_Interactors.png\n")
cat("           PCA_Sample_Clustering.png\n")
