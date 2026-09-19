#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
load_required(c("ggplot2"))
cfg <- load_config(); ensure_dirs(cfg)

expr <- utils::read.delim("/home/runner/work/DVT/DVT/results/tables/Normalized_Expression_Matrix.tsv", check.names = FALSE, row.names = 1)
expr <- as.matrix(expr); mode(expr) <- "numeric"
meta <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Sample_Metadata_Clean.csv")
mat <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/All_Matrisome_Annotation.csv")
mat_ids <- intersect(mat$ProteinGroup_ID[mat$Matrisome_status == "matched"], rownames(expr))
expr_m <- expr[mat_ids, , drop = FALSE]
meta <- meta[match(colnames(expr_m), meta$SampleID), , drop = FALSE]

expand <- expand.grid(Time = unique(meta$Time), Spatial_region = unique(meta$Spatial_region), stringsAsFactors = FALSE)
expand$Matrisome_protein_n <- apply(expand, 1, function(r) {
  s <- meta$SampleID[meta$Time == r[["Time"]] & meta$Spatial_region == r[["Spatial_region"]]]
  if (length(s) == 0) return(0)
  sum(rowMeans(!is.na(expr_m[, s, drop = FALSE])) > 0)
})
expand$mean_abundance <- apply(expand, 1, function(r) {
  s <- meta$SampleID[meta$Time == r[["Time"]] & meta$Spatial_region == r[["Spatial_region"]]]
  if (length(s) == 0) return(NA_real_)
  mean(expr_m[, s, drop = FALSE], na.rm = TRUE)
})

cat_ab <- merge(mat[, c("ProteinGroup_ID", "Matrisome_category")], data.frame(ProteinGroup_ID = rownames(expr_m), mean_expr = rowMeans(expr_m, na.rm = TRUE), stringsAsFactors = FALSE), by = "ProteinGroup_ID", all.y = TRUE)
cat_sum <- aggregate(mean_expr ~ Matrisome_category, cat_ab, mean)

utils::write.csv(expand, "/home/runner/work/DVT/DVT/results/tables/Matrisome_Time_Space_Abundance.csv", row.names = FALSE)
utils::write.csv(cat_sum, "/home/runner/work/DVT/DVT/results/tables/Matrisome_Category_Abundance.csv", row.names = FALSE)

safe_pdf("/home/runner/work/DVT/DVT/results/figures/Matrisome_Time_Space_Heatmap.pdf", {
  z <- t(scale(t(expr_m))); z[is.na(z)] <- 0
  graphics::image(t(z[nrow(z):1, , drop = FALSE]), axes = FALSE, col = grDevices::colorRampPalette(c("navy", "white", "firebrick"))(100))
})
safe_pdf("/home/runner/work/DVT/DVT/results/figures/Matrisome_Category_Composition.pdf", {
  print(ggplot2::ggplot(cat_sum, ggplot2::aes(reorder(Matrisome_category, mean_expr), mean_expr)) + ggplot2::geom_col(fill = "steelblue") + ggplot2::coord_flip() + ggplot2::theme_bw())
})
safe_pdf("/home/runner/work/DVT/DVT/results/figures/Matrisome_Trajectory.pdf", {
  p <- ggplot2::ggplot(expand, ggplot2::aes(Time, mean_abundance, group = Spatial_region, color = Spatial_region)) + ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::theme_bw()
  print(p)
})
save_session(cfg, "07_matrisome_time_space")
