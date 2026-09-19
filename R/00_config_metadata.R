#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required. Install dependencies listed in README before running.")
  }
})

required_config_fields <- c(
  "paths.raw_proteomics_file",
  "paths.metadata_file",
  "paths.matrisome_database",
  "paths.integrin_ligand_receptor_table",
  "paths.ppi_database",
  "paths.bli_ms_file",
  "paths.contrasts",
  "analysis.random_seed",
  "analysis.thresholds.fdr",
  "analysis.thresholds.abs_log2fc"
)

`%||%` <- function(x, y) if (is.null(x)) y else x

get_nested <- function(x, key) {
  parts <- strsplit(key, "\\.")[[1]]
  out <- x
  for (p in parts) {
    if (is.null(out[[p]])) return(NULL)
    out <- out[[p]]
  }
  out
}

resolve_path <- function(project_root, p) {
  if (is.null(p) || is.na(p) || p == "") return(NA_character_)
  if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(project_root, p), winslash = "/", mustWork = FALSE)
}

safe_read_delim <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt", "xls")) {
    return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

check_required_fields <- function(cfg) {
  missing <- required_config_fields[vapply(required_config_fields, function(f) is.null(get_nested(cfg, f)), logical(1))]
  if (length(missing) > 0) {
    stop("Missing required config fields: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

check_contrast_schema <- function(contrast_path) {
  required_cols <- c(
    "contrast_id", "numerator", "denominator", "variable_type", "controlled_variable",
    "biological_question", "interpretation_limit", "primary_or_exploratory"
  )
  x <- utils::read.csv(contrast_path, stringsAsFactors = FALSE, check.names = FALSE)
  miss <- setdiff(required_cols, names(x))
  if (length(miss) > 0) {
    stop("contrasts.csv missing columns: ", paste(miss, collapse = ", "))
  }
  x
}

check_metadata_expression_match <- function(metadata_path, expr_path, logs_dir) {
  if (!file.exists(metadata_path) || !file.exists(expr_path)) {
    return(data.frame(
      check = c("metadata_exists", "expression_exists", "sample_overlap"),
      status = c(file.exists(metadata_path), file.exists(expr_path), NA),
      detail = c(metadata_path, expr_path, "Skipped: missing metadata or expression file"),
      stringsAsFactors = FALSE
    ))
  }

  meta <- safe_read_delim(metadata_path)
  expr <- safe_read_delim(expr_path)

  sample_col <- if ("SampleID" %in% names(meta)) "SampleID" else names(meta)[1]
  meta_ids <- unique(as.character(meta[[sample_col]]))
  expr_ids <- colnames(expr)

  overlap <- intersect(meta_ids, expr_ids)
  only_meta <- setdiff(meta_ids, expr_ids)
  only_expr <- setdiff(expr_ids, meta_ids)

  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    data.frame(source = "metadata_only", sample_id = only_meta, stringsAsFactors = FALSE),
    file.path(logs_dir, "SampleID_Metadata_Only.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(source = "expression_only", sample_id = only_expr, stringsAsFactors = FALSE),
    file.path(logs_dir, "SampleID_Expression_Only.csv"),
    row.names = FALSE
  )

  data.frame(
    check = c("sample_overlap", "metadata_only_n", "expression_only_n"),
    status = c(length(overlap) > 0, length(only_meta) == 0, length(only_expr) == 0),
    detail = c(length(overlap), length(only_meta), length(only_expr)),
    stringsAsFactors = FALSE
  )
}

main <- function(config_path = "config/config.yml") {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  cfg <- yaml::read_yaml(resolve_path(project_root, config_path))
  check_required_fields(cfg)

  paths <- cfg$paths
  logs_dir <- resolve_path(project_root, paths$logs_dir %||% "results/logs")
  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  resolved <- list(
    raw_proteomics_file = resolve_path(project_root, paths$raw_proteomics_file),
    metadata_file = resolve_path(project_root, paths$metadata_file),
    matrisome_database = resolve_path(project_root, paths$matrisome_database),
    integrin_ligand_receptor_table = resolve_path(project_root, paths$integrin_ligand_receptor_table),
    ppi_database = resolve_path(project_root, paths$ppi_database),
    bli_ms_file = resolve_path(project_root, paths$bli_ms_file),
    contrasts = resolve_path(project_root, paths$contrasts)
  )

  file_check <- data.frame(
    key = names(resolved),
    path = unlist(resolved, use.names = FALSE),
    exists = vapply(resolved, file.exists, logical(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(file_check, file.path(logs_dir, "Config_Path_Check.csv"), row.names = FALSE)

  contrasts <- check_contrast_schema(resolved$contrasts)
  utils::write.csv(contrasts, file.path(logs_dir, "Contrast_Manifest_Checked.csv"), row.names = FALSE)

  if (file.exists(resolved$metadata_file)) {
    meta <- safe_read_delim(resolved$metadata_file)
    if ("Group" %in% names(meta)) {
      groups <- unique(as.character(meta$Group))
      contrasts$numerator_exists <- contrasts$numerator %in% groups
      contrasts$denominator_exists <- contrasts$denominator %in% groups
      utils::write.csv(contrasts, file.path(logs_dir, "Contrast_Metadata_Validation.csv"), row.names = FALSE)
    }
  }

  meta_expr_check <- check_metadata_expression_match(resolved$metadata_file, resolved$raw_proteomics_file, logs_dir)
  utils::write.csv(meta_expr_check, file.path(logs_dir, "Metadata_Expression_Consistency.csv"), row.names = FALSE)

  session_info <- utils::capture.output(sessionInfo())
  writeLines(session_info, con = file.path(logs_dir, "sessionInfo_00_config_metadata.txt"))

  message("00_config_metadata completed. Logs written to: ", logs_dir)
}

args <- commandArgs(trailingOnly = TRUE)
config_arg <- if (length(args) > 0) args[[1]] else "config/config.yml"
main(config_arg)
