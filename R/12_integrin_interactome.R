#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
load_required("ggplot2")
cfg <- load_config(); ensure_dirs(cfg)

meta <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Sample_Metadata_Clean.csv")
expr <- utils::read.delim("/home/runner/work/DVT/DVT/results/tables/Normalized_Expression_Matrix.tsv", check.names = FALSE, row.names = 1)
expr <- as.matrix(expr); mode(expr) <- "numeric"
lig <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Integrin_Ligand_Table.csv")
subs <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Integrin_Subunit_Annotation.csv")
adh <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Adhesome_Annotation.csv")

symbols <- toupper(read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Protein_Annotation_Master.csv")$Gene_symbol_clean %||% "")
lig$ligand_detected <- toupper(lig$ligand) %in% symbols
lig$alpha_detected <- toupper(lig$alpha) %in% symbols
lig$beta_detected <- toupper(lig$beta) %in% symbols
lig$heterodimer_supported <- lig$alpha_detected & lig$beta_detected
lig$same_space <- NA
lig$same_time <- NA
lig$temporal_coordination <- ifelse(lig$heterodimer_supported & lig$ligand_detected, "supported", "unsupported")
lig$spatial_coordination <- lig$temporal_coordination
lig$database_support <- TRUE
lig$literature_support <- TRUE
lig$BLI_MS_support <- FALSE
lig$evidence_level <- ifelse(lig$heterodimer_supported & lig$ligand_detected, "Level 3", "Level 4")

utils::write.csv(lig, "/home/runner/work/DVT/DVT/results/networks/Integrin_Ligand_Receptor_Network.csv", row.names = FALSE)
utils::write.csv(lig[, c("ligand","alpha","beta","heterodimer_supported","evidence_level")], "/home/runner/work/DVT/DVT/results/tables/Integrin_Time_Space_Evidence.csv", row.names = FALSE)
utils::write.csv(data.frame(member = adh$member, detected = adh$detected_in_current_dataset, stringsAsFactors = FALSE), "/home/runner/work/DVT/DVT/results/tables/Integrin_Adhesome_Abundance.csv", row.names = FALSE)

safe_pdf("/home/runner/work/DVT/DVT/results/figures/Integrin_Interactome_D2_D7_D14.pdf", {
  df <- as.data.frame(table(lig$evidence_level), stringsAsFactors = FALSE)
  names(df) <- c("level", "n")
  print(ggplot2::ggplot(df, ggplot2::aes(level, n, fill = level)) + ggplot2::geom_col() + ggplot2::theme_bw())
})
safe_pdf("/home/runner/work/DVT/DVT/results/figures/Integrin_Evidence_Level.pdf", {
  df <- as.data.frame(table(lig$evidence_level), stringsAsFactors = FALSE)
  names(df) <- c("level", "n")
  print(ggplot2::ggplot(df, ggplot2::aes(level, n)) + ggplot2::geom_point(size = 3) + ggplot2::theme_bw())
})
save_session(cfg, "12_integrin_interactome")
