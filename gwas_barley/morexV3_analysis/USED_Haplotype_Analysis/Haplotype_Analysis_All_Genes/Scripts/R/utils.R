suppressPackageStartupMessages({
  library(tools)
  library(data.table)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

eps_tag <- function(x) {
  gsub("\\.", "p", as.character(x))
}

now_ts <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

file_ts <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

as_numeric_epsilon <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

make_gene_name <- function(gene_file) {
  sub("\\.vcf\\.gz$", "", gene_file)
}

make_base_name <- function(trait) {
  paste0(trait, "_corrected_V3")
}

make_raw_vcf_dir <- function(vcf_root, trait) {
  file.path(vcf_root, paste0(make_base_name(trait), "_targets_gene_220000bp_HARMONIZED"))
}

make_imputed_vcf_dir <- function(vcf_root, trait) {
  file.path(vcf_root, paste0(make_base_name(trait), "_imputed_VCFs_for_LD_targets_gene_220000bp_HARMONIZED"))
}

make_pheno_file <- function(pheno_root, trait) {
  file.path(pheno_root, paste0(make_base_name(trait), ".pheno"))
}

list_trait_targets <- function(vcf_root, trait) {
  raw_dir <- make_raw_vcf_dir(vcf_root, trait)
  files <- list.files(raw_dir, pattern = "\\.vcf\\.gz$", full.names = FALSE)
  files <- sort(unique(files))
  data.frame(
    trait = trait,
    gene_file = files,
    title = files,
    stringsAsFactors = FALSE
  )
}

write_targets_yaml <- function(targets_df, out_path) {
  lines <- c("targets:")
  for (i in seq_len(nrow(targets_df))) {
    lines <- c(
      lines,
      paste0("  - trait: \"", targets_df$trait[i], "\""),
      paste0("    gene_file: \"", targets_df$gene_file[i], "\""),
      paste0("    title: \"", targets_df$title[i], "\""),
      ""
    )
  }
  writeLines(lines, out_path)
  invisible(out_path)
}

init_event_log <- function(path) {
  fwrite(
    data.table(
      timestamp = character(),
      level = character(),
      trait = character(),
      gene_file = character(),
      MGmin = integer(),
      epsilon = character(),
      stage = character(),
      event = character(),
      message = character()
    ),
    path
  )
  invisible(path)
}

append_event_log <- function(path, level = "INFO", trait = NA_character_, gene_file = NA_character_,
                             MGmin = NA_integer_, epsilon = NA_character_, stage = NA_character_,
                             event = NA_character_, message = NA_character_) {
  row <- data.table(
    timestamp = now_ts(),
    level = as.character(level),
    trait = as.character(trait),
    gene_file = as.character(gene_file),
    MGmin = as.integer(MGmin),
    epsilon = as.character(epsilon),
    stage = as.character(stage),
    event = as.character(event),
    message = as.character(message)
  )
  fwrite(row, path, append = TRUE)
  invisible(row)
}

log_msg <- function(..., level = "INFO", trait = NULL, gene_file = NULL, MGmin = NULL,
                    epsilon = NULL, stage = NULL, event = NULL, event_log_path = NULL) {
  parts <- c(
    paste0("[", now_ts(), "]"),
    paste0("[", level, "]"),
    if (!is.null(stage) && nzchar(stage)) paste0("[stage=", stage, "]") else NULL,
    if (!is.null(trait) && nzchar(trait)) paste0("[trait=", trait, "]") else NULL,
    if (!is.null(gene_file) && nzchar(gene_file)) paste0("[gene=", gene_file, "]") else NULL,
    if (!is.null(MGmin) && length(MGmin) == 1 && !is.na(MGmin)) paste0("[MGmin=", MGmin, "]") else NULL,
    if (!is.null(epsilon) && nzchar(as.character(epsilon))) paste0("[eps=", epsilon, "]") else NULL,
    paste(...)
  )
  line <- paste(parts, collapse = " ")
  message(line)

  if (!is.null(event_log_path)) {
    append_event_log(
      path = event_log_path,
      level = level,
      trait = trait %||% NA_character_,
      gene_file = gene_file %||% NA_character_,
      MGmin = MGmin %||% NA_integer_,
      epsilon = as.character(epsilon %||% NA_character_),
      stage = stage %||% NA_character_,
      event = event %||% NA_character_,
      message = paste(..., collapse = " ")
    )
  }

  invisible(line)
}

init_failures_df <- function() {
  data.table(
    trait = character(),
    gene_file = character(),
    gene_name = character(),
    MGmin = integer(),
    epsilon = character(),
    stage = character(),
    status = character(),
    message = character()
  )
}

add_failure <- function(failures_dt, trait, gene_file, MGmin = NA_integer_, epsilon = NA_character_,
                        stage, status, message) {
  gene_name <- make_gene_name(gene_file)
  rbind(
    failures_dt,
    data.table(
      trait = trait,
      gene_file = gene_file,
      gene_name = gene_name,
      MGmin = as.integer(MGmin),
      epsilon = as.character(epsilon),
      stage = stage,
      status = status,
      message = message
    ),
    fill = TRUE
  )
}

safe_p_adjust <- function(p, method) {
  if (length(p) == 0) return(numeric(0))
  stats::p.adjust(p, method = method)
}

group_size_signature <- function(indfile) {
  if (is.null(indfile) || nrow(indfile) == 0) return(NA_character_)
  x <- indfile
  x$hap <- as.character(x$hap)
  x$Pheno <- suppressWarnings(as.numeric(x$Pheno))
  x <- x[!is.na(x$Pheno) & x$hap != "0", , drop = FALSE]
  if (nrow(x) == 0) return(NA_character_)
  sizes <- sort(as.integer(table(x$hap)))
  paste(sizes, collapse = "|")
}

compute_valid_kw_test <- function(HapObject, label) {
  ind <- HapObject[[label]]$Indfile
  if (is.null(ind) || nrow(ind) == 0) {
    return(list(valid = FALSE, reason = "Indfile missing or empty"))
  }

  ind2 <- ind
  ind2$hap <- as.character(ind2$hap)
  ind2$Pheno <- suppressWarnings(as.numeric(ind2$Pheno))
  ind2 <- ind2[ind2$hap != "0" & !is.na(ind2$Pheno), , drop = FALSE]

  if (nrow(ind2) == 0) {
    return(list(valid = FALSE, reason = "No rows after filtering hap!=0 and non-missing phenotype"))
  }

  group_sizes <- as.integer(sort(table(ind2$hap)))
  n_groups <- length(group_sizes)
  n_ind <- nrow(ind2)

  if (n_groups <= 1) {
    return(list(valid = FALSE, reason = "Not enough groups after filtering"))
  }

  p <- tryCatch(
    kruskal.test(Pheno ~ hap, data = ind2)$p.value,
    error = function(e) NA_real_
  )

  if (is.na(p)) {
    return(list(valid = FALSE, reason = "Kruskal-Wallis p-value could not be computed"))
  }

  list(
    valid = TRUE,
    reason = NA_character_,
    kw_p_raw = as.numeric(p),
    n_ind = as.integer(n_ind),
    n_groups = as.integer(n_groups),
    group_sizes = paste(group_sizes, collapse = "|"),
    sorted_group_sizes = paste(group_sizes, collapse = "|")
  )
}