#!/usr/bin/env Rscript
# 03_run_crosshap_pipeline.R
# Orchestrator + statistics engine for step 04.
#
# For every gene (from gene_windows.tsv) x MGmin {2,3}: run crosshap (all epsilon),
# render the screening combined PDF + per-epsilon heatmaps, and compute the omnibus
# Kruskal-Wallis statistics. The statistics workflow is preserved exactly from the
# proven tool:
#   per gene x MGmin x epsilon omnibus KW (drop hap 0 / missing pheno / need >1 group)
#   -> duplicate collapse within gene (signature = raw p + n_groups + sorted sizes;
#      keep smallest epsilon, then larger MGmin)
#   -> within-gene Holm -> representative = smallest Holm p (deterministic tie-break)
#   -> per-trait across-gene BH (FDR) + Bonferroni at alpha.
#
# Outputs under 04_runs/<run_id>/: Cache, CombinedPDF, Heatmaps, Logs, Stats.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table); library(dplyr) })

argv <- commandArgs(trailingOnly = FALSE)
file_arg <- argv[grepl("^--file=", argv)]
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()

source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "run_crosshap.R"))
source(file.path(script_dir, "R", "plot_combined_pdf.R"))
source(file.path(script_dir, "R", "plot_heatmaps.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)
gw_path <- if (length(args) >= 2) args[2] else file.path(cfg$output_root, "00_config", "gene_windows.tsv")

# Defaults (mirror config)
cfg$mgmin_values   <- cfg$mgmin_values %||% c(2, 3)
cfg$epsilon_vector <- cfg$epsilon_vector %||% c(0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 0.85)
cfg$minHap         <- cfg$minHap %||% 9
cfg$hetmiss_as     <- cfg$hetmiss_as %||% "allele"
cfg$keep_outliers  <- cfg$keep_outliers %||% FALSE
cfg$overwrite_cache<- cfg$overwrite_cache %||% FALSE
cfg$continue_on_error <- cfg$continue_on_error %||% TRUE
cfg$fail_fast      <- cfg$fail_fast %||% FALSE
cfg$stats          <- cfg$stats %||% list()
cfg$stats$alpha    <- cfg$stats$alpha %||% 0.05
cfg$stats$within_gene_method     <- cfg$stats$within_gene_method %||% "holm"
cfg$stats$across_gene_fdr_method <- cfg$stats$across_gene_fdr_method %||% "BH"

run_root      <- file.path(cfg$output_root, "04_runs", cfg$run_id)
cache_root    <- file.path(run_root, "Cache")
combined_root <- file.path(run_root, "CombinedPDF")
heatmaps_root <- file.path(run_root, "Heatmaps")
logs_root     <- file.path(run_root, "Logs")
stats_root    <- file.path(run_root, "Stats")
tmp_root      <- file.path(run_root, "tmp")
for (d in list(run_root, cache_root, combined_root, heatmaps_root, logs_root, stats_root, tmp_root)) ensure_dir(d)

event_log_path <- file.path(logs_root, paste0("event_log_", cfg$run_id, "_", file_ts(), ".csv"))
init_event_log(event_log_path)

status_from_invalid_reason <- function(reason) {
  if (is.na(reason) || !nzchar(reason)) return("no_valid_kw_pvalue")
  if (grepl("Not enough groups|No rows after filtering|empty", reason, ignore.case = TRUE))
    return("no_valid_haplotype_grouping")
  "no_valid_kw_pvalue"
}

empty_per_test_stats <- function() data.table(
  trait=character(), gene_file=character(), gene_name=character(), title=character(),
  MGmin=integer(), epsilon=character(), test_signature=character(), n_ind=integer(),
  n_groups=integer(), group_sizes=character(), kw_p_raw=numeric(),
  kw_p_holm_within_gene=numeric(), is_best_within_gene=logical(),
  combined_pdf_path=character(), heatmap_pdf_path=character(), notes=character())

empty_gene_summary <- function() data.table(
  trait=character(), gene_file=character(), gene_name=character(), title=character(),
  lead_SNP=character(), class=character(), lead_pos=integer(), lead_neg_log10p=numeric(),
  dist_to_lead_bp=integer(), locus_id=character(), description=character(), n_snps_window=integer(),
  n_total_valid_tests_before_collapse=integer(), n_unique_valid_tests=integer(),
  best_MGmin=integer(), best_epsilon=character(), best_n_ind=integer(), best_n_groups=integer(),
  best_group_sizes=character(), best_kw_p_raw=numeric(), best_kw_p_holm_within_gene=numeric(),
  trait_fdr_p=numeric(), trait_bonferroni_p=numeric(),
  significant_fdr=logical(), significant_bonferroni=logical())

empty_successful_runs <- function() data.table(
  trait=character(), gene_file=character(), gene_name=character(), title=character(),
  MGmin=integer(), cache_rds=character(), combined_pdf_path=character())

log_msg("Run started.", stage="run", event="run_start", event_log_path=event_log_path)
log_msg(paste0("Run ID: ", cfg$run_id, " | window_bp=", cfg$window_bp,
               " | MGmin=", paste(cfg$mgmin_values, collapse=","),
               " | eps=", paste(cfg$epsilon_vector, collapse=",")),
        stage="run", event="run_config", event_log_path=event_log_path)

targets_dt <- as.data.table(read_gene_windows(gw_path))
setorder(targets_dt, trait, gene_file)
if (nrow(targets_dt) == 0) stop("No targets in ", gw_path)
log_msg(paste0("Targets: ", nrow(targets_dt), " genes from ", gw_path),
        stage="discovery", event="targets_loaded", event_log_path=event_log_path)

failures_dt        <- init_failures_df()
per_test_dt        <- empty_per_test_stats()
gene_summary_dt    <- empty_gene_summary()
successful_runs_dt <- empty_successful_runs()

for (i in seq_len(nrow(targets_dt))) {
  tr        <- targets_dt$trait[i]
  gene_file <- targets_dt$gene_file[i]
  title     <- targets_dt$title[i]
  gene_name <- make_gene_name(gene_file)
  gene_n_common <- NA_integer_

  log_msg(paste0("Processing gene ", i, "/", nrow(targets_dt)),
          trait=tr, gene_file=gene_file, stage="gene", event="gene_start", event_log_path=event_log_path)

  gene_valid_dt <- empty_per_test_stats()[0]

  for (MGmin in cfg$mgmin_values) {
    log_msg("Starting MGmin analysis.", trait=tr, gene_file=gene_file, MGmin=MGmin,
            stage="mgmin", event="mgmin_start", event_log_path=event_log_path)

    cache_dir <- file.path(cache_root, tr, gene_name, paste0("MGmin_", MGmin)); ensure_dir(cache_dir)
    cache_rds <- file.path(cache_dir, "HapObject.rds")
    res <- NULL

    if (file.exists(cache_rds) && !isTRUE(cfg$overwrite_cache)) {
      res <- tryCatch(readRDS(cache_rds), error = function(e) {
        failures_dt <<- add_failure(failures_dt, tr, gene_file, MGmin, NA_character_, "cache",
                                    "cache_read_failed", conditionMessage(e))
        log_msg(paste0("Cache read failed: ", conditionMessage(e)), level="ERROR", trait=tr,
                gene_file=gene_file, MGmin=MGmin, stage="cache", event="failure_encountered",
                event_log_path=event_log_path); NULL })
      if (!is.null(res)) log_msg(paste0("Using cached HapObject: ", cache_rds), trait=tr,
                                 gene_file=gene_file, MGmin=MGmin, stage="cache", event="cache_used",
                                 event_log_path=event_log_path)
    } else {
      log_msg("Running CrossHap.", trait=tr, gene_file=gene_file, MGmin=MGmin, stage="crosshap",
              event="crosshap_start", event_log_path=event_log_path)
      res <- tryCatch(
        run_crosshap(cfg = cfg, trait = tr, gene_file = gene_file, MGmin = MGmin, tmp_dir = tmp_root),
        error = function(e) {
          failures_dt <<- add_failure(failures_dt, tr, gene_file, MGmin, NA_character_, "crosshap",
                                      "crosshap_failed", conditionMessage(e))
          log_msg(paste0("CrossHap failed: ", conditionMessage(e)), level="ERROR", trait=tr,
                  gene_file=gene_file, MGmin=MGmin, stage="crosshap", event="failure_encountered",
                  event_log_path=event_log_path); NULL })
      if (!is.null(res)) {
        saveRDS(res, cache_rds)
        log_msg(paste0("Saved cache: ", cache_rds), trait=tr, gene_file=gene_file, MGmin=MGmin,
                stage="cache", event="cache_written", event_log_path=event_log_path)
      }
    }

    if (is.null(res)) {
      if (isTRUE(cfg$fail_fast) && !isTRUE(cfg$continue_on_error)) stop("fail_fast and not continue_on_error.")
      next
    }
    if (is.na(gene_n_common)) gene_n_common <- as.integer(res$n_common %||% NA_integer_)

    comb_dir <- file.path(combined_root, tr, gene_name, paste0("MGmin_", MGmin)); ensure_dir(comb_dir)
    combined_pdf_path <- file.path(comb_dir, paste0("CrosshapTree+Violin__", gene_name,
                                   "__MGmin", MGmin, "__EpsVec", length(cfg$epsilon_vector), ".pdf"))

    successful_runs_dt <- rbind(successful_runs_dt, data.table(trait=tr, gene_file=gene_file,
      gene_name=gene_name, title=title, MGmin=as.integer(MGmin), cache_rds=cache_rds,
      combined_pdf_path=combined_pdf_path), fill = TRUE)

    tryCatch(write_combined_pdf(HapObject=res$HapObject, out_pdf=combined_pdf_path, title=title,
      trait=tr, gene_file=gene_file, gene_name=gene_name, MGmin=MGmin,
      epsilon_vector=cfg$epsilon_vector, mgmin_test_stats=NULL, gene_summary_row=NULL),
      error = function(e) {
        failures_dt <<- add_failure(failures_dt, tr, gene_file, MGmin, NA_character_, "combined_pdf",
                                    "combined_pdf_failed", conditionMessage(e))
        log_msg(paste0("Combined PDF failed: ", conditionMessage(e)), level="ERROR", trait=tr,
                gene_file=gene_file, MGmin=MGmin, stage="combined_pdf", event="failure_encountered",
                event_log_path=event_log_path) })

    heatmap_dir <- file.path(heatmaps_root, tr, gene_name, paste0("MGmin_", MGmin)); ensure_dir(heatmap_dir)

    for (eps in cfg$epsilon_vector) {
      label <- paste0("Haplotypes_MGmin", MGmin, "_E", eps)
      heatmap_pdf_path <- file.path(heatmap_dir, paste0("Heatmaps__", gene_name, "__MGmin", MGmin,
                                    "__Eps", eps_tag(eps), ".pdf"))
      if (is.null(res$HapObject[[label]])) {
        failures_dt <- add_failure(failures_dt, tr, gene_file, MGmin, as.character(eps),
                                   "epsilon_validation", "no_valid_haplotype_grouping",
                                   "No HapObject entry for this epsilon.")
        next
      }
      tryCatch(write_heatmaps_for_eps(HapObject=res$HapObject, label=label, gene_file=gene_file,
        title=title, trait=tr, eps=eps, MGmin=MGmin, raw_path=res$raw_path,
        common_ids=res$common_ids, vcf_raw_ids=res$vcf_raw_ids, out_pdf=heatmap_pdf_path),
        error = function(e) {
          failures_dt <<- add_failure(failures_dt, tr, gene_file, MGmin, as.character(eps), "heatmap",
                                      "heatmap_failed", conditionMessage(e))
          log_msg(paste0("Heatmap failed: ", conditionMessage(e)), level="ERROR", trait=tr,
                  gene_file=gene_file, MGmin=MGmin, epsilon=as.character(eps), stage="heatmap",
                  event="failure_encountered", event_log_path=event_log_path) })

      kw_info <- compute_valid_kw_test(res$HapObject, label)
      if (!isTRUE(kw_info$valid)) {
        failures_dt <- add_failure(failures_dt, tr, gene_file, MGmin, as.character(eps),
                                   "kw_validation", status_from_invalid_reason(kw_info$reason), kw_info$reason)
        next
      }
      test_signature <- paste(formatC(kw_info$kw_p_raw, digits=15, format="fg", flag="#"),
                              kw_info$n_groups, kw_info$sorted_group_sizes, sep="__")
      gene_valid_dt <- rbind(gene_valid_dt, data.table(trait=tr, gene_file=gene_file,
        gene_name=gene_name, title=title, MGmin=as.integer(MGmin), epsilon=as.character(eps),
        test_signature=test_signature, n_ind=as.integer(kw_info$n_ind),
        n_groups=as.integer(kw_info$n_groups), group_sizes=kw_info$group_sizes,
        kw_p_raw=as.numeric(kw_info$kw_p_raw), kw_p_holm_within_gene=NA_real_,
        is_best_within_gene=FALSE, combined_pdf_path=combined_pdf_path,
        heatmap_pdf_path=heatmap_pdf_path, notes=""), fill = TRUE)
      log_msg(paste0("Valid KW test. p=", format_p_plain(kw_info$kw_p_raw), " groups=",
              kw_info$n_groups, " sizes=", kw_info$group_sizes), trait=tr, gene_file=gene_file,
              MGmin=MGmin, epsilon=as.character(eps), stage="kw_validation",
              event="valid_kw_test_found", event_log_path=event_log_path)
    }
  }

  n_before <- nrow(gene_valid_dt)
  if (n_before == 0) {
    failures_dt <- add_failure(failures_dt, tr, gene_file, NA_integer_, NA_character_, "gene_stats",
                               "no_valid_tests_for_gene", "Zero valid omnibus tests across MGmin/epsilon.")
    log_msg("Gene produced zero valid omnibus tests.", level="WARN", trait=tr, gene_file=gene_file,
            stage="gene_stats", event="failure_encountered", event_log_path=event_log_path)
    next
  }

  # ---- Duplicate collapse within gene (smallest epsilon, then larger MGmin) ----
  gene_valid_dt[, epsilon_num := as_numeric_epsilon(epsilon)]
  setorder(gene_valid_dt, epsilon_num, -MGmin, test_signature)
  gene_unique_dt <- copy(gene_valid_dt[!duplicated(gene_valid_dt$test_signature)])
  gene_valid_dt[, epsilon_num := NULL]; gene_unique_dt[, epsilon_num := NULL]

  dup_counts <- gene_valid_dt[, .N, by = test_signature]
  gene_unique_dt <- merge(gene_unique_dt, dup_counts, by="test_signature", all.x=TRUE, sort=FALSE)
  setnames(gene_unique_dt, "N", "duplicate_count")
  for (k in seq_len(nrow(gene_unique_dt))) {
    sig <- gene_unique_dt$test_signature[k]; dup_rows <- gene_valid_dt[test_signature == sig]
    if (nrow(dup_rows) > 1) {
      gene_unique_dt$notes[k] <- paste0("Duplicate signature collapsed from ",
        paste0("MGmin", dup_rows$MGmin, ":eps", dup_rows$epsilon, collapse=";"))
      log_msg(paste0("Collapsed duplicates. Signature=", sig), trait=tr, gene_file=gene_file,
              MGmin=gene_unique_dt$MGmin[k], epsilon=gene_unique_dt$epsilon[k], stage="collapse",
              event="duplicate_collapsed", event_log_path=event_log_path)
    }
  }

  # ---- Within-gene Holm ----
  gene_unique_dt$kw_p_holm_within_gene <- safe_p_adjust(gene_unique_dt$kw_p_raw, method=cfg$stats$within_gene_method)
  log_msg(paste0("Within-gene Holm done. before=", n_before, " unique=", nrow(gene_unique_dt)),
          trait=tr, gene_file=gene_file, stage="within_gene_correction", event="holm_correction_complete",
          event_log_path=event_log_path)

  # ---- Representative test (smallest Holm p; deterministic tie-break) ----
  best_holm_p <- min(gene_unique_dt$kw_p_holm_within_gene, na.rm = TRUE)
  best_idx <- which(gene_unique_dt$kw_p_holm_within_gene == best_holm_p)
  if (length(best_idx) != 1) {
    tied <- copy(gene_unique_dt[best_idx]); tied[, en := as_numeric_epsilon(epsilon)]
    setorder(tied, kw_p_raw, en, -MGmin); chosen <- tied[1]
    log_msg(paste0("AMBIGUOUS_BEST_TEST_RESOLVED: holm_p=", best_holm_p, " chosen MGmin=", chosen$MGmin,
            " eps=", chosen$epsilon), level="WARN", trait=tr, gene_file=gene_file, MGmin=chosen$MGmin,
            epsilon=chosen$epsilon, stage="within_gene_summary", event="ambiguous_best_test_resolved",
            event_log_path=event_log_path)
    best_idx <- which(gene_unique_dt$MGmin==chosen$MGmin & gene_unique_dt$epsilon==chosen$epsilon &
                      gene_unique_dt$kw_p_raw==chosen$kw_p_raw &
                      gene_unique_dt$kw_p_holm_within_gene==chosen$kw_p_holm_within_gene)[1]
  }
  gene_unique_dt$is_best_within_gene[best_idx] <- TRUE
  gene_unique_dt$duplicate_count <- NULL
  per_test_dt <- rbind(per_test_dt, gene_unique_dt, fill = TRUE)

  best <- gene_unique_dt[best_idx]
  gene_summary_dt <- rbind(gene_summary_dt, data.table(
    trait=tr, gene_file=gene_file, gene_name=gene_name, title=title,
    lead_SNP=targets_dt$lead_SNP[i], class=targets_dt$class[i],
    lead_pos=as.integer(targets_dt$lead_pos[i]), lead_neg_log10p=as.numeric(targets_dt$lead_neg_log10p[i]),
    dist_to_lead_bp=as.integer(targets_dt$dist_to_lead_bp[i]), locus_id=targets_dt$locus_id[i],
    description=as.character(targets_dt$description[i]), n_snps_window=gene_n_common,
    n_total_valid_tests_before_collapse=as.integer(n_before),
    n_unique_valid_tests=as.integer(nrow(gene_unique_dt)),
    best_MGmin=as.integer(best$MGmin), best_epsilon=as.character(best$epsilon),
    best_n_ind=as.integer(best$n_ind), best_n_groups=as.integer(best$n_groups),
    best_group_sizes=as.character(best$group_sizes), best_kw_p_raw=as.numeric(best$kw_p_raw),
    best_kw_p_holm_within_gene=as.numeric(best$kw_p_holm_within_gene),
    trait_fdr_p=NA_real_, trait_bonferroni_p=NA_real_, significant_fdr=NA, significant_bonferroni=NA), fill=TRUE)

  log_msg(paste0("Representative: MGmin=", best$MGmin, " eps=", best$epsilon, " raw_p=",
          format_p_plain(best$kw_p_raw), " holm_p=", format_p_plain(best$kw_p_holm_within_gene)),
          trait=tr, gene_file=gene_file, MGmin=best$MGmin, epsilon=best$epsilon,
          stage="within_gene_summary", event="best_test_selected", event_log_path=event_log_path)
}

# ---- Per-trait across-gene correction ----
if (nrow(gene_summary_dt) > 0) {
  gene_summary_dt$trait_fdr_p <- NA_real_; gene_summary_dt$trait_bonferroni_p <- NA_real_
  gene_summary_dt$significant_fdr <- FALSE; gene_summary_dt$significant_bonferroni <- FALSE
  for (tn in unique(gene_summary_dt$trait)) {
    idx <- which(gene_summary_dt$trait == tn)
    fdr  <- safe_p_adjust(gene_summary_dt$best_kw_p_holm_within_gene[idx], method=cfg$stats$across_gene_fdr_method)
    bonf <- safe_p_adjust(gene_summary_dt$best_kw_p_holm_within_gene[idx], method="bonferroni")
    gene_summary_dt$trait_fdr_p[idx] <- fdr
    gene_summary_dt$trait_bonferroni_p[idx] <- bonf
    gene_summary_dt$significant_fdr[idx] <- fdr <= cfg$stats$alpha
    gene_summary_dt$significant_bonferroni[idx] <- bonf <= cfg$stats$alpha
    log_msg(paste0("Trait correction done. genes=", length(idx), " fdr=", cfg$stats$across_gene_fdr_method,
            " alpha=", cfg$stats$alpha), trait=tn, stage="trait_correction",
            event="trait_level_correction_complete", event_log_path=event_log_path)
  }
}

# ---- Regenerate combined PDFs with corrected values ----
if (nrow(successful_runs_dt) > 0) {
  successful_runs_dt <- unique(successful_runs_dt); setorder(successful_runs_dt, trait, gene_file, MGmin)
  for (i in seq_len(nrow(successful_runs_dt))) {
    tr <- successful_runs_dt$trait[i]; gene_file <- successful_runs_dt$gene_file[i]
    gene_name <- successful_runs_dt$gene_name[i]; title <- successful_runs_dt$title[i]
    MGmin <- successful_runs_dt$MGmin[i]; cache_rds <- successful_runs_dt$cache_rds[i]
    combined_pdf_path <- successful_runs_dt$combined_pdf_path[i]
    res <- tryCatch(readRDS(cache_rds), error=function(e) NULL)
    if (is.null(res)) next
    mgmin_test_stats <- as.data.frame(per_test_dt[trait==tr & gene_file==gene_file & MGmin==MGmin], stringsAsFactors=FALSE)
    gene_summary_row <- as.data.frame(gene_summary_dt[trait==tr & gene_file==gene_file], stringsAsFactors=FALSE)
    tryCatch(write_combined_pdf(HapObject=res$HapObject, out_pdf=combined_pdf_path, title=title,
      trait=tr, gene_file=gene_file, gene_name=gene_name, MGmin=MGmin, epsilon_vector=cfg$epsilon_vector,
      mgmin_test_stats=mgmin_test_stats, gene_summary_row=gene_summary_row),
      error=function(e) {
        failures_dt <<- add_failure(failures_dt, tr, gene_file, MGmin, NA_character_,
                                    "combined_pdf_postprocess", "combined_pdf_postprocess_failed", conditionMessage(e))
        log_msg(paste0("Corrected PDF regen failed: ", conditionMessage(e)), level="ERROR", trait=tr,
                gene_file=gene_file, MGmin=MGmin, stage="combined_pdf_postprocess",
                event="failure_encountered", event_log_path=event_log_path) })
  }
}

# ---- Write tables ----
if (nrow(per_test_dt) > 0)     setorder(per_test_dt, trait, gene_file, MGmin, epsilon)
if (nrow(gene_summary_dt) > 0) setorder(gene_summary_dt, trait, best_kw_p_holm_within_gene, best_kw_p_raw)
if (nrow(failures_dt) > 0)     setorder(failures_dt, trait, gene_file, MGmin, epsilon, stage)

fwrite(per_test_dt,     file.path(stats_root, "per_test_stats.csv"))
fwrite(gene_summary_dt, file.path(stats_root, "gene_summary.csv"))
fwrite(failures_dt,     file.path(stats_root, "failures_notes.csv"))

log_msg(paste0("Run finished. targets=", nrow(targets_dt), " per_test=", nrow(per_test_dt),
        " gene_summary=", nrow(gene_summary_dt), " failures=", nrow(failures_dt)),
        stage="run", event="run_end", event_log_path=event_log_path)
cat(sprintf("\nDONE. gene_summary rows=%d, per_test rows=%d, failures=%d\nStats: %s\n",
            nrow(gene_summary_dt), nrow(per_test_dt), nrow(failures_dt), stats_root))
