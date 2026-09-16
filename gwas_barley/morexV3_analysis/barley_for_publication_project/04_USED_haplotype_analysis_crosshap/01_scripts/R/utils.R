suppressPackageStartupMessages({
  library(tools)
  library(data.table)
})

# ---------------------------------------------------------------------------
# utils.R -- shared helpers for the step-04 crosshap pipeline.
# Ported from the proven all-genes tool; path/sample logic updated for the new
# bcftools-cut, single-IID, gene+/-1000bp inputs. Statistics helpers unchanged.
# ---------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

eps_tag         <- function(x) gsub("\\.", "p", as.character(x))
now_ts          <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
file_ts         <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
as_numeric_epsilon <- function(x) suppressWarnings(as.numeric(as.character(x)))
make_gene_name  <- function(gene_file) sub("\\.vcf\\.gz$", "", gene_file)

# ---- Input path helpers (config-driven; new layout) ----
raw_vcf_path     <- function(cfg, trait, gene_file) file.path(cfg$raw_vcf_dir, trait, gene_file)
imputed_vcf_path <- function(cfg, trait, gene_file) file.path(cfg$imputed_vcf_dir, trait, gene_file)
pheno_path       <- function(cfg, trait) file.path(cfg$pheno_root, paste0(trait, cfg$pheno_suffix))

# ---- Targets from gene_windows.tsv ----
# Returns one row per gene with gene_file (= gene_id.vcf.gz), a human title that
# embeds the origin lead SNP + class, and all origin/annotation columns.
read_gene_windows <- function(gw_path) {
  gw <- read.delim(gw_path, stringsAsFactors = FALSE)
  gw$gene_file <- paste0(gw$gene_id, ".vcf.gz")
  gw$title <- sprintf("%s | %s | lead %s (%s)",
                      gw$gene_id, gw$trait, gw$lead_SNP, gw$class)
  gw
}

# ---- Event log ----
init_event_log <- function(path) {
  fwrite(data.table(timestamp = character(), level = character(), trait = character(),
                    gene_file = character(), MGmin = integer(), epsilon = character(),
                    stage = character(), event = character(), message = character()), path)
  invisible(path)
}

append_event_log <- function(path, level = "INFO", trait = NA_character_, gene_file = NA_character_,
                             MGmin = NA_integer_, epsilon = NA_character_, stage = NA_character_,
                             event = NA_character_, message = NA_character_) {
  fwrite(data.table(timestamp = now_ts(), level = as.character(level), trait = as.character(trait),
                    gene_file = as.character(gene_file), MGmin = as.integer(MGmin),
                    epsilon = as.character(epsilon), stage = as.character(stage),
                    event = as.character(event), message = as.character(message)),
         path, append = TRUE)
  invisible(NULL)
}

log_msg <- function(..., level = "INFO", trait = NULL, gene_file = NULL, MGmin = NULL,
                    epsilon = NULL, stage = NULL, event = NULL, event_log_path = NULL) {
  parts <- c(paste0("[", now_ts(), "]"), paste0("[", level, "]"),
             if (!is.null(stage) && nzchar(stage)) paste0("[stage=", stage, "]"),
             if (!is.null(trait) && nzchar(trait)) paste0("[trait=", trait, "]"),
             if (!is.null(gene_file) && nzchar(gene_file)) paste0("[gene=", gene_file, "]"),
             if (!is.null(MGmin) && length(MGmin) == 1 && !is.na(MGmin)) paste0("[MGmin=", MGmin, "]"),
             if (!is.null(epsilon) && nzchar(as.character(epsilon))) paste0("[eps=", epsilon, "]"),
             paste(...))
  line <- paste(parts, collapse = " ")
  message(line)
  if (!is.null(event_log_path)) {
    append_event_log(event_log_path, level = level, trait = trait %||% NA_character_,
                     gene_file = gene_file %||% NA_character_, MGmin = MGmin %||% NA_integer_,
                     epsilon = as.character(epsilon %||% NA_character_), stage = stage %||% NA_character_,
                     event = event %||% NA_character_, message = paste(..., collapse = " "))
  }
  invisible(line)
}

# ---- Failures table ----
init_failures_df <- function() {
  data.table(trait = character(), gene_file = character(), gene_name = character(),
             MGmin = integer(), epsilon = character(), stage = character(),
             status = character(), message = character())
}

add_failure <- function(failures_dt, trait, gene_file, MGmin = NA_integer_, epsilon = NA_character_,
                        stage, status, message) {
  rbind(failures_dt,
        data.table(trait = trait, gene_file = gene_file, gene_name = make_gene_name(gene_file),
                   MGmin = as.integer(MGmin), epsilon = as.character(epsilon),
                   stage = stage, status = status, message = message), fill = TRUE)
}

safe_p_adjust <- function(p, method) {
  if (length(p) == 0) return(numeric(0))
  stats::p.adjust(p, method = method)
}

# ---- Omnibus Kruskal-Wallis validity (UNCHANGED science) ----
# Valid test = >1 haplotype group after dropping hap "0" and missing phenotype,
# and a computable KW p-value.
compute_valid_kw_test <- function(HapObject, label) {
  ind <- HapObject[[label]]$Indfile
  if (is.null(ind) || nrow(ind) == 0)
    return(list(valid = FALSE, reason = "Indfile missing or empty"))

  ind2 <- ind
  ind2$hap   <- as.character(ind2$hap)
  ind2$Pheno <- suppressWarnings(as.numeric(ind2$Pheno))
  ind2 <- ind2[ind2$hap != "0" & !is.na(ind2$Pheno), , drop = FALSE]
  if (nrow(ind2) == 0)
    return(list(valid = FALSE, reason = "No rows after filtering hap!=0 and non-missing phenotype"))

  group_sizes <- as.integer(sort(table(ind2$hap)))
  n_groups <- length(group_sizes); n_ind <- nrow(ind2)
  if (n_groups <= 1)
    return(list(valid = FALSE, reason = "Not enough groups after filtering"))

  p <- tryCatch(kruskal.test(Pheno ~ hap, data = ind2)$p.value, error = function(e) NA_real_)
  if (is.na(p))
    return(list(valid = FALSE, reason = "Kruskal-Wallis p-value could not be computed"))

  list(valid = TRUE, reason = NA_character_, kw_p_raw = as.numeric(p),
       n_ind = as.integer(n_ind), n_groups = as.integer(n_groups),
       group_sizes = paste(group_sizes, collapse = "|"),
       sorted_group_sizes = paste(group_sizes, collapse = "|"))
}
