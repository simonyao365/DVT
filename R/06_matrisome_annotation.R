#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
cfg <- load_config(); ensure_dirs(cfg)

prot <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Protein_Annotation_Master.csv")
mat <- if (file.exists(cfg$._resolved$matrisome_database)) read_delim_auto(cfg$._resolved$matrisome_database) else data.frame()
if (nrow(mat) > 0) {
  id_col <- names(mat)[1]
  cat_col <- if ("Matrisome_category" %in% names(mat)) "Matrisome_category" else names(mat)[min(2, ncol(mat))]
  key <- toupper(trimws(as.character(mat[[id_col]])))
  m <- match(toupper(trimws(prot$Gene_symbol_clean %||% "")), key)
  prot$Matrisome_status <- ifelse(!is.na(m), "matched", "unmapped")
  prot$Matrisome_category <- ifelse(!is.na(m), as.character(mat[[cat_col]][m]), NA_character_)
} else {
  prot$Matrisome_status <- "unmapped"
  prot$Matrisome_category <- NA_character_
}
prot$Matrisome_subcategory <- prot$Matrisome_category
prot$matched_ID <- ifelse(prot$Matrisome_status == "matched", prot$Gene_symbol_clean, NA)
prot$mapping_method <- "symbol_exact_uppercase"
prot$mapping_confidence <- ifelse(prot$Matrisome_status == "matched", "high", "none")
prot$multiple_category_status <- FALSE

utils::write.csv(prot, "/home/runner/work/DVT/DVT/results/tables/All_Matrisome_Annotation.csv", row.names = FALSE)
utils::write.csv(prot[, c("ProteinGroup_ID", "Matrisome_status", "matched_ID", "mapping_method", "mapping_confidence")], "/home/runner/work/DVT/DVT/results/logs/Matrisome_Mapping_Log.csv", row.names = FALSE)
utils::write.csv(prot[prot$Matrisome_status != "matched", ], "/home/runner/work/DVT/DVT/results/logs/Unmapped_Matrisome_Proteins.csv", row.names = FALSE)
utils::write.csv(as.data.frame(table(prot$Matrisome_category, useNA = "ifany")), "/home/runner/work/DVT/DVT/results/tables/Matrisome_Category_Summary.csv", row.names = FALSE)
save_session(cfg, "06_matrisome_annotation")
