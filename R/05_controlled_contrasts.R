#!/usr/bin/env Rscript
source("R/pipeline_utils.R")
load_required(c("limma", "ggplot2"))
cfg <- load_config(); ensure_dirs(cfg)

expr <- utils::read.delim("/home/runner/work/DVT/DVT/results/tables/Imputed_Expression_Matrix.tsv", check.names = FALSE, row.names = 1)
expr <- as.matrix(expr); mode(expr) <- "numeric"
meta <- read_delim_auto("/home/runner/work/DVT/DVT/results/tables/Sample_Metadata_Clean.csv")
meta <- meta[match(colnames(expr), meta$SampleID), , drop = FALSE]
if (any(is.na(meta$SampleID))) stop("Metadata does not cover all expression samples")

meta$Group <- factor(meta$Group)
design <- model.matrix(~ 0 + Group, data = meta)
colnames(design) <- levels(meta$Group)
fit <- limma::lmFit(expr, design)
cons <- read_contrasts(cfg)
if (nrow(cons) == 0) stop("No contrasts found")

mk <- function(n, d) paste0("`", n, "`-`", d, "`")
contrast_strings <- setNames(vapply(seq_len(nrow(cons)), function(i) mk(cons$numerator[i], cons$denominator[i]), character(1)), cons$contrast_id)
cmat <- limma::makeContrasts(contrasts = contrast_strings, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cmat))

for (cid in colnames(cmat)) {
  res <- limma::topTable(fit2, coef = cid, number = Inf, sort.by = "P")
  res$FDR <- res$adj.P.Val
  res$effect_size <- res$logFC
  res$B_statistic <- res$B
  res$contrast_id <- cid
  out_csv <- file.path("/home/runner/work/DVT/DVT/results/tables", paste0("Contrast_", cid, "_limma.csv"))
  utils::write.csv(res, out_csv, row.names = TRUE)

  p1 <- ggplot2::ggplot(res, ggplot2::aes(logFC, -log10(P.Value), color = FDR < (cfg$analysis$thresholds$fdr %||% 0.05))) + ggplot2::geom_point(alpha = 0.5) + ggplot2::theme_bw() + ggplot2::labs(title = paste("Volcano", cid))
  p2 <- ggplot2::ggplot(res, ggplot2::aes(AveExpr, logFC, color = FDR < (cfg$analysis$thresholds$fdr %||% 0.05))) + ggplot2::geom_point(alpha = 0.5) + ggplot2::theme_bw() + ggplot2::labs(title = paste("MA", cid))
  safe_pdf(file.path("/home/runner/work/DVT/DVT/results/figures", paste0("Volcano_", cid, ".pdf")), { print(p1) })
  safe_pdf(file.path("/home/runner/work/DVT/DVT/results/figures", paste0("MA_", cid, ".pdf")), { print(p2) })
}

utils::write.csv(cons, "/home/runner/work/DVT/DVT/results/logs/contrast_manifest.csv", row.names = FALSE)
write_md("/home/runner/work/DVT/DVT/results/logs/Contrast_Interpretation_Limits.md", c("# Contrast interpretation limits", "- Primary conclusions should use controlled time/space contrasts.", "- Mixed-variable contrasts are exploratory composite-state differences only."))
save_session(cfg, "05_controlled_contrasts")
