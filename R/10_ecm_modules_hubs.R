#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
load_required("igraph")
cfg <- load_config(); ensure_dirs(cfg)

edges <- read_delim_auto("/home/runner/work/DVT/DVT/results/networks/ECM_PPI_All_Network.csv")
if (nrow(edges) == 0) stop("No edges available for module detection")
g <- igraph::graph_from_data_frame(edges[, c("Protein_A", "Protein_B")], directed = FALSE)

set.seed(cfg$analysis$random_seed %||% 123)
mods <- list(
  Louvain = igraph::cluster_louvain(g),
  Infomap = igraph::cluster_infomap(g),
  Walktrap = igraph::cluster_walktrap(g)
)

memb <- do.call(rbind, lapply(names(mods), function(nm) data.frame(algorithm = nm, node = names(igraph::membership(mods[[nm]])), module = as.integer(igraph::membership(mods[[nm]])), stringsAsFactors = FALSE)))
mod_stat <- do.call(rbind, lapply(names(mods), function(nm) data.frame(algorithm = nm, modularity = igraph::modularity(mods[[nm]]), n_modules = length(unique(igraph::membership(mods[[nm]]))), stringsAsFactors = FALSE)))
cent <- data.frame(
  node = igraph::V(g)$name,
  degree = igraph::degree(g),
  weighted_degree = igraph::strength(g),
  betweenness = igraph::betweenness(g),
  eigenvector = igraph::eigen_centrality(g)$vector,
  PageRank = igraph::page_rank(g)$vector,
  k_core = igraph::coreness(g),
  stringsAsFactors = FALSE
)
cent$within_module_z <- scale(cent$degree)
cent$participation_coefficient <- 1 - (cent$degree / max(cent$degree))^2
cent$hub_type <- ifelse(cent$degree >= stats::quantile(cent$degree, 0.9), "Structural hub", "Non-hub")

utils::write.csv(memb, "/home/runner/work/DVT/DVT/results/networks/ECM_Module_Membership.csv", row.names = FALSE)
utils::write.csv(mod_stat, "/home/runner/work/DVT/DVT/results/networks/ECM_Module_Stability.csv", row.names = FALSE)
utils::write.csv(data.frame(module = unique(memb$module), enrichment = "not_run"), "/home/runner/work/DVT/DVT/results/networks/ECM_Module_Enrichment.csv", row.names = FALSE)
utils::write.csv(cent, "/home/runner/work/DVT/DVT/results/networks/ECM_Network_Centrality.csv", row.names = FALSE)
utils::write.csv(cent[, c("node", "within_module_z", "participation_coefficient")], "/home/runner/work/DVT/DVT/results/networks/ECM_zP_Roles.csv", row.names = FALSE)
utils::write.csv(cent[order(-cent$degree), ], "/home/runner/work/DVT/DVT/results/networks/ECM_Hub_Prioritization.csv", row.names = FALSE)
save_session(cfg, "10_ecm_modules_hubs")
