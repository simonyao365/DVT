#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
cfg <- load_config(); ensure_dirs(cfg)

expr <- utils::read.delim("/home/runner/work/DVT/DVT/results/tables/Normalized_Expression_Matrix.tsv", check.names = FALSE, row.names = 1)
expr <- as.matrix(expr); mode(expr) <- "numeric"
meta <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Sample_Metadata_Clean.csv")
mat <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/All_Matrisome_Annotation.csv")
mat_ids <- intersect(mat$ProteinGroup_ID[mat$Matrisome_status == "matched"], rownames(expr))

time_order <- c("D2", "D7", "D14")
meta$Time <- factor(meta$Time, levels = time_order)

traj <- lapply(mat_ids, function(pid) {
  tmp <- data.frame(SampleID = colnames(expr), abundance = as.numeric(expr[pid, ]), stringsAsFactors = FALSE)
  tmp <- merge(tmp, meta[, c("SampleID", "Time", "Spatial_region", "Group")], by = "SampleID", all.x = TRUE)
  agg <- aggregate(abundance ~ Time + Spatial_region, tmp, mean)
  out <- data.frame(ProteinGroup_ID = pid, stringsAsFactors = FALSE)
  out$trend_positive <- all(diff(aggregate(abundance ~ Time, agg, mean)$abundance) > 0, na.rm = TRUE)
  out$det_rate <- mean(!is.na(tmp$abundance))
  out$effect_size <- max(tmp$abundance, na.rm = TRUE) - min(tmp$abundance, na.rm = TRUE)
  out
})
traj_df <- do.call(rbind, traj)
traj_df$trajectory_cluster <- ifelse(traj_df$trend_positive, "persistent_up", "non_monotonic")

utils::write.csv(traj_df[traj_df$trend_positive & traj_df$det_rate >= 0.5, ], "/home/runner/work/DVT/DVT/results/tables/Persistent_Matrisome_Proteins.csv", row.names = FALSE)
utils::write.csv(traj_df[!traj_df$trend_positive, ], "/home/runner/work/DVT/DVT/results/tables/NonMonotonic_Matrisome_Proteins.csv", row.names = FALSE)
utils::write.csv(traj_df[, c("ProteinGroup_ID", "trajectory_cluster", "effect_size", "det_rate")], "/home/runner/work/DVT/DVT/results/tables/Matrisome_Trajectory_Clusters.csv", row.names = FALSE)

safe_pdf("/home/runner/work/DVT/DVT/results/figures/Persistent_Matrisome_Heatmap.pdf", {
  keep <- traj_df$ProteinGroup_ID[traj_df$trend_positive]
  if (length(keep) == 0) keep <- mat_ids
  m <- expr[intersect(keep, rownames(expr)), , drop = FALSE]
  z <- t(scale(t(m))); z[is.na(z)] <- 0
  graphics::image(t(z[nrow(z):1, , drop = FALSE]), axes = FALSE, col = grDevices::colorRampPalette(c("navy", "white", "firebrick"))(100))
})
save_session(cfg, "08_persistent_matrisome")
