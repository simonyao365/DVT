#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

load_required <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) {
    stop("Missing required packages: ", paste(miss, collapse = ", "), ". Install dependencies listed in README.")
  }
  invisible(TRUE)
}

load_config <- function(config_path = "config/config.yml") {
  load_required("yaml")
  cfg <- yaml::read_yaml(config_path)
  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  resolve <- function(p) {
    if (is.null(p) || is.na(p) || p == "") return(NA_character_)
    if (grepl("^(/|[A-Za-z]:[/\\\\])", p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    normalizePath(file.path(root, p), winslash = "/", mustWork = FALSE)
  }
  cfg$._resolved <- list(
    raw_proteomics_file = resolve(cfg$paths$raw_proteomics_file),
    metadata_file = resolve(cfg$paths$metadata_file),
    matrisome_database = resolve(cfg$paths$matrisome_database),
    integrin_ligand_receptor_table = resolve(cfg$paths$integrin_ligand_receptor_table),
    ppi_database = resolve(cfg$paths$ppi_database),
    bli_ms_file = resolve(cfg$paths$bli_ms_file),
    output_root = resolve(cfg$paths$output_root),
    logs_dir = resolve(cfg$paths$logs_dir),
    tables_dir = resolve(cfg$paths$tables_dir),
    figures_dir = resolve(cfg$paths$figures_dir),
    networks_dir = resolve(cfg$paths$networks_dir),
    contrasts = resolve(cfg$paths$contrasts)
  )
  cfg
}

ensure_dirs <- function(cfg) {
  dirs <- unlist(cfg$._resolved[c("output_root", "logs_dir", "tables_dir", "figures_dir", "networks_dir")], use.names = FALSE)
  for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

save_session <- function(cfg, module_id) {
  ensure_dirs(cfg)
  out <- file.path(cfg$._resolved$logs_dir, paste0("sessionInfo_", module_id, ".txt"))
  writeLines(capture.output(sessionInfo()), out)
}

read_delim_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt", "xls")) {
    return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

infer_metadata_from_expr <- function(expr_df) {
  qty_cols <- grep("PG\\.Quantity$", names(expr_df), value = TRUE)
  if (length(qty_cols) == 0) stop("No quantity columns detected (pattern: PG.Quantity).")
  raw_name <- gsub("^\\[[0-9]+\\]\\s*", "", qty_cols)
  raw_name <- gsub("\\.raw\\.PG\\.Quantity$", "", raw_name)

  guess_time <- function(x) {
    if (grepl("(^|-)14($|[A-Z])", x)) return("D14")
    if (grepl("(^|-)7($|[A-Z])", x)) return("D7")
    if (grepl("(^|-)2($|[A-Z])", x)) return("D2")
    NA_character_
  }
  guess_region <- function(x) {
    if (grepl("[A-Z]C$|C$", x)) return("Collagen")
    if (grepl("[A-Z]F$|F$", x)) return("Fibrin")
    if (grepl("[A-Z]V$|V$", x)) return("Wall")
    if (grepl("[A-Z]H$|H$", x)) return("Healthy")
    NA_character_
  }
  time <- vapply(raw_name, guess_time, character(1))
  region <- vapply(raw_name, guess_region, character(1))
  group <- ifelse(is.na(time) | is.na(region), NA_character_, paste(time, region, sep = "_"))

  data.frame(
    SampleID = qty_cols,
    RawSampleName = raw_name,
    Time = time,
    Spatial_region = region,
    Group = group,
    Batch = "Batch1",
    stringsAsFactors = FALSE
  )
}

load_proteomics <- function(cfg) {
  expr_raw <- read_delim_auto(cfg$._resolved$raw_proteomics_file)
  qty_cols <- grep("PG\\.Quantity$", names(expr_raw), value = TRUE)
  if (length(qty_cols) == 0) stop("No quantity columns found in proteomics input.")

  expr <- as.matrix(expr_raw[, qty_cols, drop = FALSE])
  mode(expr) <- "numeric"

  if ("PG.ProteinGroups" %in% names(expr_raw)) {
    row_ids <- as.character(expr_raw$PG.ProteinGroups)
  } else {
    row_ids <- paste0("ProteinGroup_", seq_len(nrow(expr_raw)))
  }
  row_ids[is.na(row_ids) | row_ids == ""] <- paste0("ProteinGroup_", which(is.na(row_ids) | row_ids == ""))
  rownames(expr) <- make.unique(row_ids)

  feature <- expr_raw
  rownames(feature) <- rownames(expr)

  metadata <- if (file.exists(cfg$._resolved$metadata_file)) {
    md <- read_delim_auto(cfg$._resolved$metadata_file)
    if (!"SampleID" %in% names(md)) stop("Metadata must contain 'SampleID'.")
    md
  } else {
    infer_metadata_from_expr(expr_raw)
  }

  list(expr = expr, feature = feature, metadata = metadata)
}

align_expr_metadata <- function(expr, metadata) {
  common <- intersect(colnames(expr), metadata$SampleID)
  if (length(common) == 0) stop("No overlap between expression columns and metadata SampleID.")
  meta2 <- metadata[match(common, metadata$SampleID), , drop = FALSE]
  expr2 <- expr[, common, drop = FALSE]
  list(expr = expr2, metadata = meta2)
}

safe_pdf <- function(path, code) {
  grDevices::pdf(path)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
}

write_md <- function(path, lines) {
  writeLines(lines, con = path)
}

read_contrasts <- function(cfg) {
  if (!file.exists(cfg$._resolved$contrasts)) return(data.frame())
  utils::read.csv(cfg$._resolved$contrasts, stringsAsFactors = FALSE, check.names = FALSE)
}

classify_evidence_level <- function(score) {
  cut(score, breaks = c(-Inf, 1, 2, 3, 4, Inf), labels = c("Level 5", "Level 4", "Level 3", "Level 2", "Level 1"))
}
