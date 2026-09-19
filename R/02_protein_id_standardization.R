#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
cfg <- load_config(); ensure_dirs(cfg)

x <- load_proteomics(cfg)
feature <- x$feature

split_ids <- function(v) unique(unlist(strsplit(as.character(v %||% ""), "[;,| ]+")))
extract_gn <- function(desc) {
  y <- stringr::str_extract(as.character(desc), "GN=[^ ]+")
  gsub("GN=", "", y)
}

if (!requireNamespace("stringr", quietly = TRUE)) stop("Missing package: stringr")
raw_symbol <- if ("PG.Genes" %in% names(feature)) as.character(feature$PG.Genes) else rep(NA_character_, nrow(feature))
gn_from_desc <- if ("PG.ProteinDescriptions" %in% names(feature)) extract_gn(feature$PG.ProteinDescriptions) else rep(NA_character_, nrow(feature))
gene_clean <- ifelse(is.na(raw_symbol) | raw_symbol == "", gn_from_desc, raw_symbol)

df <- data.frame(
  ProteinGroup_ID = rownames(feature),
  Gene_symbol_raw = raw_symbol,
  Gene_symbol_clean = gene_clean,
  UniProt_ID = if ("PG.ProteinGroups" %in% names(feature)) sapply(feature$PG.ProteinGroups, function(x) paste(split_ids(x), collapse = ";")) else NA_character_,
  Ensembl_ID = NA_character_,
  Entrez_ID = NA_character_,
  Protein_name = if ("PG.ProteinDescriptions" %in% names(feature)) as.character(feature$PG.ProteinDescriptions) else NA_character_,
  Species = cfg$project$species_primary %||% "mouse",
  Isoform = grepl("-[0-9]+", if ("PG.ProteinGroups" %in% names(feature)) feature$PG.ProteinGroups else ""),
  Mapping_status = ifelse(is.na(gene_clean) | gene_clean == "", "unmapped", "mapped"),
  Mapping_source = ifelse(is.na(raw_symbol) | raw_symbol == "", "description_GN", "PG.Genes"),
  Ambiguous_mapping = grepl(";|,|\\|", gene_clean %||% ""),
  Unique_protein_key = make.unique(rownames(feature)),
  stringsAsFactors = FALSE
)

utils::write.csv(df, "/home/runner/work/DVT/DVT/results/tables/Protein_Annotation_Master.csv", row.names = FALSE)
utils::write.csv(df[, c("ProteinGroup_ID", "Mapping_status", "Mapping_source", "Ambiguous_mapping")], "/home/runner/work/DVT/DVT/results/logs/Protein_ID_Mapping_Log.csv", row.names = FALSE)
utils::write.csv(df[df$Mapping_status == "unmapped", ], "/home/runner/work/DVT/DVT/results/logs/Unmapped_Proteins.csv", row.names = FALSE)
utils::write.csv(df[duplicated(df$Gene_symbol_clean) & !is.na(df$Gene_symbol_clean) & df$Gene_symbol_clean != "", ], "/home/runner/work/DVT/DVT/results/logs/Duplicated_Gene_Symbols.csv", row.names = FALSE)
save_session(cfg, "02_protein_id_standardization")
