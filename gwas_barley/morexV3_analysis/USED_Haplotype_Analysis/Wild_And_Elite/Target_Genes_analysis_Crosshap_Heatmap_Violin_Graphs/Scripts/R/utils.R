suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

now_ts <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

file_ts <- function() {
  format(Sys.time(), "%Y%m%d_%H%M%S")
}

make_gene_name <- function(gene_file) {
  sub("\\.vcf\\.gz$", "", gene_file)
}

make_base_name <- function(trait) {
  paste0(trait, "_corrected_V3")
}

make_imputed_vcf_dir <- function(ld_vcf_root, trait) {
  file.path(ld_vcf_root, paste0(make_base_name(trait), "_imputed_VCFs_for_LD_targets_gene_220000bp_HARMONIZED"))
}

make_pheno_file <- function(pheno_root, trait) {
  file.path(pheno_root, paste0(make_base_name(trait), "_with_elite.pheno"))
}

make_merged_vcf_path <- function(merged_vcf_dir, merged_vcf_file) {
  file.path(merged_vcf_dir, merged_vcf_file)
}

make_imputed_vcf_path <- function(ld_vcf_root, trait, gene_file) {
  file.path(make_imputed_vcf_dir(ld_vcf_root, trait), gene_file)
}

normalize_target <- function(target) {
  target$trait <- as.character(target$trait)
  target$gene_file <- as.character(target$gene_file)
  target$merged_vcf_file <- as.character(target$merged_vcf_file)
  target$title <- as.character(target$title)
  target$epsilon <- as.numeric(target$epsilon)
  target$MGmin <- as.integer(target$MGmin)
  if (is.null(target$force_ref) || is.na(target$force_ref) || identical(target$force_ref, "null")) {
    target$force_ref <- NULL
  } else {
    target$force_ref <- as.character(target$force_ref)
  }
  target
}

init_event_log <- function(path) {
  fwrite(
    data.table(
      timestamp = character(),
      level = character(),
      trait = character(),
      gene_file = character(),
      stage = character(),
      event = character(),
      message = character()
    ),
    path
  )
  invisible(path)
}

append_event_log <- function(path, level = "INFO", trait = NA_character_, gene_file = NA_character_,
                             stage = NA_character_, event = NA_character_, message = NA_character_) {
  fwrite(
    data.table(
      timestamp = now_ts(),
      level = as.character(level),
      trait = as.character(trait),
      gene_file = as.character(gene_file),
      stage = as.character(stage),
      event = as.character(event),
      message = as.character(message)
    ),
    path,
    append = TRUE
  )
}

log_msg <- function(..., level = "INFO", trait = NULL, gene_file = NULL, stage = NULL,
                    event = NULL, event_log_path = NULL) {
  parts <- c(
    paste0("[", now_ts(), "]"),
    paste0("[", level, "]"),
    if (!is.null(stage) && nzchar(stage)) paste0("[stage=", stage, "]") else NULL,
    if (!is.null(trait) && nzchar(trait)) paste0("[trait=", trait, "]") else NULL,
    if (!is.null(gene_file) && nzchar(gene_file)) paste0("[gene=", gene_file, "]") else NULL,
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
    stage = character(),
    status = character(),
    message = character()
  )
}

add_failure <- function(failures_dt, trait, gene_file, stage, status, message) {
  rbind(
    failures_dt,
    data.table(
      trait = trait,
      gene_file = gene_file,
      gene_name = make_gene_name(gene_file),
      stage = stage,
      status = status,
      message = message
    ),
    fill = TRUE
  )
}

load_pheno_wild_elite <- function(pheno_file) {
  raw_pheno <- data.table::fread(pheno_file, header = FALSE, fill = TRUE)
  colnames(raw_pheno) <- c("FID", "IID", "trait")

  raw_pheno$Ind <- ifelse(grepl("^HS", raw_pheno$FID), paste(raw_pheno$FID, raw_pheno$IID, sep = "_"), raw_pheno$FID)
  raw_pheno$Pheno <- suppressWarnings(as.numeric(raw_pheno$trait))

  out <- unique(raw_pheno[, .(Ind, Pheno)])
  out
}

format_p_plain <- function(p, digits = 4, min_print = 1e-4) {
  if (is.null(p) || length(p) == 0 || is.na(p)) return("NA")
  if (p < min_print) return(paste0("<", formatC(min_print, format = "f", digits = digits)))
  formatC(p, format = "f", digits = digits)
}

p_to_signif_symbol <- function(p) {
  if (is.na(p)) return("ns")
  if (p <= 0.0001) return("****")
  if (p <= 0.001) return("***")
  if (p <= 0.01) return("**")
  if (p <= 0.05) return("*")
  "ns"
}
