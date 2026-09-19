#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
cfg <- load_config(); ensure_dirs(cfg)

expr <- utils::read.delim("/home/runner/work/DVT/DVT/results/tables/Normalized_Expression_Matrix.tsv", check.names = FALSE, row.names = 1)
expr <- as.matrix(expr); mode(expr) <- "numeric"
mat <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/All_Matrisome_Annotation.csv")
ids <- intersect(mat$ProteinGroup_ID[mat$Matrisome_status == "matched"], rownames(expr))

build_edges <- function(ids) {
  if (length(ids) < 2) return(data.frame())
  cmb <- utils::combn(ids, 2)
  data.frame(Protein_A = cmb[1, ], Protein_B = cmb[2, ], stringsAsFactors = FALSE)
}

if (file.exists(cfg$._resolved$ppi_database)) {
  ppi <- read_delim_auto(cfg$._resolved$ppi_database)
  names(ppi)[1:2] <- c("Protein_A", "Protein_B")
  edges <- ppi
} else {
  edges <- build_edges(ids)
}

edges$Species <- cfg$project$species_primary %||% "mouse"
edges$Database <- ifelse(file.exists(cfg$._resolved$ppi_database), "Configured_PPI", "CoDetected_Exploratory")
edges$Evidence_type <- ifelse(edges$Database == "Configured_PPI", "database", "co-detected")
edges$Experimental_system <- "NA"
edges$Confidence <- ifelse(edges$Database == "Configured_PPI", 0.8, 0.2)
edges$Direct_or_indirect <- ifelse(edges$Database == "Configured_PPI", "mixed", "indirect")
edges$Current_dataset_support <- edges$Protein_A %in% ids & edges$Protein_B %in% ids
edges$BLI_MS_support <- FALSE
edges$Literature_support <- NA

utils::write.csv(edges, "/home/runner/work/DVT/DVT/results/networks/ECM_PPI_All_Network.csv", row.names = FALSE)
utils::write.csv(edges[edges$Confidence >= 0.7, ], "/home/runner/work/DVT/DVT/results/networks/ECM_PPI_HighConfidence.csv", row.names = FALSE)
utils::write.csv(edges[edges$Confidence < 0.7, ], "/home/runner/work/DVT/DVT/results/networks/ECM_PPI_Exploratory.csv", row.names = FALSE)
utils::write.csv(edges[, c("Protein_A", "Protein_B", "Database", "Evidence_type", "Confidence")], "/home/runner/work/DVT/DVT/results/networks/ECM_PPI_Database_Evidence.csv", row.names = FALSE)
utils::write.csv(data.frame(node = unique(c(edges$Protein_A, edges$Protein_B))), "/home/runner/work/DVT/DVT/results/networks/Cytoscape_Nodes.csv", row.names = FALSE)
utils::write.csv(edges[, c("Protein_A", "Protein_B", "Confidence", "Evidence_type")], "/home/runner/work/DVT/DVT/results/networks/Cytoscape_Edges.csv", row.names = FALSE)
utils::write.csv(data.frame(metric = c("n_nodes", "n_edges"), value = c(length(unique(c(edges$Protein_A, edges$Protein_B))), nrow(edges))), "/home/runner/work/DVT/DVT/results/networks/ECM_Network_Stats.csv", row.names = FALSE)
save_session(cfg, "09_ecm_ppi_network")
