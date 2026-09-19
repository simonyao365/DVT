#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
cfg <- load_config(); ensure_dirs(cfg)

prot <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Protein_Annotation_Master.csv")
lig <- if (file.exists(cfg$._resolved$integrin_ligand_receptor_table)) read_delim_auto(cfg$._resolved$integrin_ligand_receptor_table) else data.frame(ligand=character(),alpha=character(),beta=character())

integrin_alpha <- paste0("ITGA", 1:11)
integrin_beta <- paste0("ITGB", c(1:8))
hetero <- expand.grid(alpha = integrin_alpha, beta = integrin_beta, stringsAsFactors = FALSE)
hetero$heterodimer <- paste0(hetero$alpha, "/", hetero$beta)

prot$symbol_upper <- toupper(prot$Gene_symbol_clean %||% "")
sub_anno <- data.frame(
  gene = c(integrin_alpha, integrin_beta),
  class = c(rep("alpha", length(integrin_alpha)), rep("beta", length(integrin_beta))),
  detected_in_current_dataset = c(integrin_alpha, integrin_beta) %in% prot$symbol_upper,
  stringsAsFactors = FALSE
)

if (nrow(lig) == 0) {
  lig <- data.frame(ligand = c("FN1", "COL1A1"), alpha = c("ITGA5", "ITGA2"), beta = c("ITGB1", "ITGB1"), stringsAsFactors = FALSE)
}
lig$heterodimer <- paste0(lig$alpha, "/", lig$beta)
lig$evidence_type <- "curated"
lig$species <- cfg$project$species_primary %||% "mouse"
lig$source <- "config_or_default"
lig$direct_ligand_status <- TRUE
lig$detected_in_current_dataset <- toupper(lig$ligand) %in% prot$symbol_upper

adhesome <- data.frame(member = c("TLN1","TLN2","FERMT1","FERMT2","FERMT3","FLNA","FLNB","PTK2","SRC","PXN","VCL","ILK","ACTN1","RHOA","ROCK1","PIK3CA","AKT1","MAPK1","YAP1","WWTR1","TGFB1"), stringsAsFactors = FALSE)
adhesome$detected_in_current_dataset <- adhesome$member %in% prot$symbol_upper

utils::write.csv(sub_anno, "/home/runner/work/DVT/DVT/results/tables/Integrin_Subunit_Annotation.csv", row.names = FALSE)
utils::write.csv(hetero, "/home/runner/work/DVT/DVT/results/tables/Integrin_Heterodimer_Table.csv", row.names = FALSE)
utils::write.csv(lig, "/home/runner/work/DVT/DVT/results/tables/Integrin_Ligand_Table.csv", row.names = FALSE)
utils::write.csv(adhesome, "/home/runner/work/DVT/DVT/results/tables/Adhesome_Annotation.csv", row.names = FALSE)
utils::write.csv(data.frame(resource = c("integrin_ligand_receptor","adhesome_core"), version = "TO_BE_FILLED", download_date = "TO_BE_FILLED"), "/home/runner/work/DVT/DVT/results/logs/Annotation_Source_Version.csv", row.names = FALSE)
save_session(cfg, "11_integrin_annotation")
