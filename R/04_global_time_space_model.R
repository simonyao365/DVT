#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
cfg <- load_config(); ensure_dirs(cfg)

load_required("limma")
expr <- as.matrix(read_delim_auto(cfg$._resolved$tables_dir %||% "/home/runner/work/DVT/DVT/results/tables/Imputed_Expression_Matrix.tsv"))
if (colnames(expr)[1] == "") { rownames(expr) <- expr[,1]; expr <- expr[,-1,drop=FALSE] }
mode(expr) <- "numeric"
meta <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Sample_Metadata_Clean.csv")
if (!"SampleID" %in% names(meta)) stop("Sample metadata missing SampleID")
meta <- meta[meta$SampleID %in% colnames(expr), , drop = FALSE]
expr <- expr[, meta$SampleID, drop = FALSE]
meta$Time <- factor(meta$Time)
meta$Spatial_region <- factor(meta$Spatial_region)
meta$Batch <- factor(meta$Batch %||% "Batch1")

has_interaction <- nlevels(meta$Time) > 1 && nlevels(meta$Spatial_region) > 1 && nrow(meta) >= (nlevels(meta$Time) * nlevels(meta$Spatial_region))
form <- if (has_interaction) ~ 0 + Time + Spatial_region + Time:Spatial_region + Batch else ~ 0 + Time + Spatial_region + Batch
design <- model.matrix(form, data = meta)
fit <- limma::eBayes(limma::lmFit(expr, design))
res <- limma::topTable(fit, number = Inf, sort.by = "P")
res$FDR <- res$adj.P.Val
res$CI_low <- res$logFC - 1.96 * res$SE
res$CI_high <- res$logFC + 1.96 * res$SE

utils::write.csv(res, "/home/runner/work/DVT/DVT/results/tables/Global_Time_Space_Model.csv", row.names = TRUE)
var_df <- data.frame(term = colnames(design), variance = apply(design, 2, var), stringsAsFactors = FALSE)
utils::write.csv(var_df, "/home/runner/work/DVT/DVT/results/tables/Global_Variance_Decomposition.csv", row.names = FALSE)
utils::write.csv(data.frame(model = deparse(form), has_interaction = has_interaction, stringsAsFactors = FALSE), "/home/runner/work/DVT/DVT/results/logs/Global_Model_Design.csv", row.names = FALSE)
write_md("/home/runner/work/DVT/DVT/results/logs/Global_Model_Interpretation.md", c("# Global model interpretation", "- Results include combined coefficient-level effects.", "- Time, spatial, and interaction interpretations must follow coefficient naming.", "- Mixed contrasts remain exploratory only."))
save_session(cfg, "04_global_time_space_model")
