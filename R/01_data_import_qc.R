#!/usr/bin/env Rscript
source("R/pipeline_utils.R")

load_required(c("ggplot2"))
cfg <- load_config()
ensure_dirs(cfg)
set.seed(cfg$analysis$random_seed %||% 123)

x <- load_proteomics(cfg)
a <- align_expr_metadata(x$expr, x$metadata)
expr <- a$expr
metadata <- a$metadata
feature <- x$feature[rownames(expr), , drop = FALSE]

expr_num <- expr
expr_num[expr_num == 0] <- NA
expr_log2 <- log2(expr_num)

sample_detected <- colSums(!is.na(expr_num))
protein_detect_rate <- rowMeans(!is.na(expr_num))
zero_count <- sum(expr == 0, na.rm = TRUE)
na_count <- sum(is.na(expr))

dup_pg <- if ("PG.ProteinGroups" %in% names(feature)) sum(duplicated(feature$PG.ProteinGroups)) else NA_integer_
dup_gene <- if ("PG.Genes" %in% names(feature)) sum(duplicated(feature$PG.Genes)) else NA_integer_

group_counts <- as.data.frame(table(metadata$Group), stringsAsFactors = FALSE)
colnames(group_counts) <- c("Group", "N")

sample_total_intensity <- colSums(expr, na.rm = TRUE)
cor_mat <- suppressWarnings(stats::cor(expr_log2, use = "pairwise.complete.obs", method = "pearson"))

pc <- stats::prcomp(t(expr_log2), center = TRUE, scale. = TRUE)
pc_df <- data.frame(SampleID = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], stringsAsFactors = FALSE)
pc_df <- merge(pc_df, metadata, by = "SampleID", all.x = TRUE)

hc <- stats::hclust(dist(t(expr_log2)))

all_qc <- data.frame(
  metric = c("n_proteins", "n_samples", "zero_count", "missing_count", "duplicated_proteingroup", "duplicated_gene", "mean_sample_detected", "mean_protein_detect_rate"),
  value = c(nrow(expr), ncol(expr), zero_count, na_count, dup_pg, dup_gene, mean(sample_detected), mean(protein_detect_rate)),
  stringsAsFactors = FALSE
)

utils::write.csv(all_qc, "/home/runner/work/DVT/DVT/results/tables/All_Protein_QC_Summary.csv", row.names = FALSE)
utils::write.csv(metadata, "/home/runner/work/DVT/DVT/results/tables/Sample_Metadata_Clean.csv", row.names = FALSE)
utils::write.csv(feature, "/home/runner/work/DVT/DVT/results/tables/Protein_Feature_Clean.csv", row.names = TRUE)
utils::write.table(expr, "/home/runner/work/DVT/DVT/results/tables/Expression_Matrix_Clean.tsv", sep = "\t", quote = FALSE, col.names = NA)

safe_pdf("/home/runner/work/DVT/DVT/results/figures/Sample_Correlation.pdf", {
  graphics::par(mar = c(6, 6, 2, 1))
  graphics::image(1:ncol(cor_mat), 1:ncol(cor_mat), cor_mat[nrow(cor_mat):1, ], xaxt = "n", yaxt = "n", col = grDevices::colorRampPalette(c("navy", "white", "firebrick"))(100))
  graphics::axis(1, at = 1:ncol(cor_mat), labels = colnames(cor_mat), las = 2, cex.axis = 0.6)
  graphics::axis(2, at = 1:ncol(cor_mat), labels = rev(colnames(cor_mat)), las = 2, cex.axis = 0.6)
})

safe_pdf("/home/runner/work/DVT/DVT/results/figures/PCA_Global.pdf", {
  print(ggplot2::ggplot(pc_df, ggplot2::aes(PC1, PC2, color = Group, label = RawSampleName)) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::geom_text(size = 2, vjust = -0.8, show.legend = FALSE) +
    ggplot2::theme_bw())
})

safe_pdf("/home/runner/work/DVT/DVT/results/figures/Missingness_Overview.pdf", {
  hist_df <- data.frame(detect_rate = protein_detect_rate)
  print(ggplot2::ggplot(hist_df, ggplot2::aes(detect_rate)) + ggplot2::geom_histogram(bins = 30, fill = "steelblue") + ggplot2::theme_bw())
})

outlier_cut <- mean(sample_detected) - 2 * stats::sd(sample_detected)
outliers <- names(sample_detected)[sample_detected < outlier_cut]

qc_md <- c(
  "# QC Report",
  paste0("- Input expression dimensions: ", nrow(expr), " proteins x ", ncol(expr), " samples"),
  paste0("- Groups: ", paste(sprintf("%s=%s", group_counts$Group, group_counts$N), collapse = ", ")),
  paste0("- Zero count: ", zero_count),
  paste0("- Missing count (NA): ", na_count),
  paste0("- Duplicated ProteinGroup: ", dup_pg),
  paste0("- Duplicated Gene symbol: ", dup_gene),
  paste0("- Potential outlier samples: ", ifelse(length(outliers) == 0, "None", paste(outliers, collapse = ", "))),
  "- Batch effect: placeholder (requires explicit batch metadata)",
  "- Abnormal samples flagged by low detected-protein count only in this module."
)
write_md("/home/runner/work/DVT/DVT/results/logs/QC_Report.md", qc_md)
save_session(cfg, "01_data_import_qc")
