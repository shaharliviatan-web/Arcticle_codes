#!/usr/bin/env Rscript
# =============================================================================
# 03_run_crosshap_pipeline.R  --  orchestrator + statistics engine for step 04.
#
# REWRITTEN 2026-09-08. The v1 version is in
#   _ARCHIVE_v1_epsilon_sweep_2026-09-08/03_run_crosshap_pipeline.R
#
# WHAT IT DOES
#   For every gene in 00_config/gene_windows.tsv:
#     1. run crosshap ONCE  (MGmin = 2, epsilon = 0.6, both fixed in config.yaml)
#     2. render the combined tree+violin PDF and the heatmap PDF
#     3. run ONE Kruskal-Wallis test:  phenotype ~ haplotype group
#     4. compute effect sizes (per-group means, eta-squared, top-vs-bottom in SD)
#   then, per trait, apply Benjamini-Hochberg FDR to the RAW p-values.
#
# WHAT CHANGED FROM v1, AND WHY
#   v1 swept 7 epsilon x 2 MGmin per gene, collapsed duplicate results, applied a
#   within-gene Holm correction, took the smallest Holm p as the gene's result, and
#   then applied BOTH BH and Bonferroni across genes. That was removed because:
#
#   (a) DOUBLE CORRECTION. v1 fed Holm-ADJUSTED p-values into BH. BH applied to
#       FWER-adjusted values does not control FDR at any interpretable rate, so the
#       reported "FDR" was not an FDR. Here BH receives the RAW p.
#   (b) ARBITRARY PENALTY. The Holm penalty a gene paid equalled its number of unique
#       epsilon results, which ranged 1 to 13 -- a property of the gene's SNP structure,
#       not of the hypothesis. 12 genes paid no penalty at all; 8 of those were called
#       significant.
#   (c) SELECTION ON THE OUTCOME. Taking the epsilon with the smallest p is a forking
#       path. Fixing epsilon on genotype-only grounds (crosshap never uses the phenotype
#       to build groups) removes it. For 3 of the 4 v1 publication genes the epsilon
#       choice changed nothing at all -- they were flat across the whole sweep.
#
#   BH is retained despite genes within a locus being in LD: BH controls FDR under
#   positive regression dependence (Benjamini & Yekutieli 2001), which LD-induced
#   correlation satisfies. The caveat is one of REPORTING, not validity -- significant
#   genes cluster in loci, so locus_summary.tsv reports loci alongside genes.
#
# UNCHANGED, DELIBERATELY (the proven core -- see R/run_crosshap.R)
#   LD comes from the IMPUTED per-gene VCF (complete r2 matrix); haplotyping runs on the
#   RAW per-gene VCF (real genotypes); the two are intersected on CHROM:POS. Caching,
#   event logging and failure capture are v1's and are kept as-is.
#
# OUTPUTS  04_runs/<run_id>/Stats/
#   gene_results.tsv       one row per TESTED gene: grouping, KW p, BH q, effect sizes
#   haplotype_groups.tsv   one row per gene x haplotype group: n, mean, sd, median
#   genes_not_tested.tsv   genes with no valid grouping, with the reason
#   locus_summary.tsv      one row per locus: genes tested / significant
#   per_trait_summary.tsv  one row per trait
# =============================================================================

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
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)
gw_path <- if (length(args) >= 2) args[2] else file.path(cfg$output_root, "00_config", "gene_windows.tsv")

## ---- Fixed parameters: exactly one MGmin and one epsilon ----------------------
cfg$mgmin_values   <- cfg$mgmin_values   %||% 2
cfg$epsilon_vector <- cfg$epsilon_vector %||% 0.6
if (length(cfg$mgmin_values) != 1L || length(cfg$epsilon_vector) != 1L) {
  stop("This pipeline is single-configuration by design: config.yaml must give exactly ",
       "one mgmin_values and one epsilon_vector entry (got ",
       length(cfg$mgmin_values), " and ", length(cfg$epsilon_vector), "). ",
       "To explore parameters, use the v1 script in _ARCHIVE_v1_epsilon_sweep_2026-09-08/.")
}
MGMIN <- as.integer(cfg$mgmin_values[[1]])
EPS   <- as.numeric(cfg$epsilon_vector[[1]])
cfg$minHap        <- cfg$minHap        %||% 9
cfg$hetmiss_as    <- cfg$hetmiss_as    %||% "allele"
cfg$keep_outliers <- cfg$keep_outliers %||% FALSE
cfg$overwrite_cache   <- cfg$overwrite_cache   %||% FALSE
cfg$continue_on_error <- cfg$continue_on_error %||% TRUE
cfg$stats         <- cfg$stats         %||% list()
cfg$stats$alpha   <- cfg$stats$alpha   %||% 0.05
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

log_msg("Run started.", stage="run", event="run_start", event_log_path=event_log_path)
log_msg(sprintf("run_id=%s window_bp=%s MGmin=%d epsilon=%s alpha=%s",
                cfg$run_id, cfg$window_bp, MGMIN, EPS, cfg$stats$alpha),
        stage="run", event="run_config", event_log_path=event_log_path)

## =============================================================================
## One Kruskal-Wallis test + effect sizes for one haplotype grouping.
##
## Filtering (unchanged from v1): drop hap "0" (crosshap's unassigned label) and any
## accession with a missing phenotype; require >= 2 remaining groups.
##
## eta_squared = (H - k + 1) / (n - k)   -- the standard KW effect size. It is the
##   proportion of phenotype rank variance explained by haplotype group. Reported
##   because a KW p-value at a locus that was SELECTED for association with this very
##   phenotype is close to guaranteed; the effect size is what a reader can judge.
## delta_top_bottom_sd = (mean of the highest group - mean of the lowest) / SD of all
##   phenotype values used, i.e. the spread between haplotypes in phenotype SD units.
## =============================================================================
kw_test_with_effects <- function(HapObject, label) {
  ind <- HapObject[[label]]$Indfile
  if (is.null(ind) || nrow(ind) == 0)
    return(list(valid = FALSE, reason = "Indfile missing or empty"))

  d <- data.table(hap = as.character(ind$hap),
                  Pheno = suppressWarnings(as.numeric(ind$Pheno)),
                  Ind = as.character(ind$Ind))
  n_unassigned <- sum(d$hap == "0", na.rm = TRUE)
  d <- d[hap != "0" & !is.na(Pheno)]
  if (nrow(d) == 0)
    return(list(valid = FALSE, reason = "No accessions left after dropping hap 0 and missing phenotypes"))
  k <- uniqueN(d$hap)
  if (k < 2)
    return(list(valid = FALSE, reason = sprintf("Only %d haplotype group after filtering (need >= 2)", k)))

  kt <- tryCatch(kruskal.test(Pheno ~ hap, data = d), error = function(e) NULL)
  if (is.null(kt) || is.na(kt$p.value))
    return(list(valid = FALSE, reason = "Kruskal-Wallis p-value could not be computed"))

  n <- nrow(d); H <- unname(kt$statistic)
  eta2 <- if (n > k) (H - k + 1) / (n - k) else NA_real_
  eta2 <- if (is.na(eta2)) NA_real_ else max(0, min(1, eta2))

  grp <- d[, .(n = .N, mean = mean(Pheno), median = median(Pheno),
               sd = if (.N > 1) sd(Pheno) else NA_real_), by = hap][order(-mean)]
  sd_all <- sd(d$Pheno)
  delta_sd <- if (!is.na(sd_all) && sd_all > 0) (grp$mean[1] - grp$mean[nrow(grp)]) / sd_all else NA_real_

  list(valid = TRUE, reason = NA_character_,
       kw_p = as.numeric(kt$p.value), kw_H = H, df = unname(kt$parameter),
       n_ind = n, n_groups = k, n_unassigned = n_unassigned,
       group_sizes = paste(sort(grp$n, decreasing = TRUE), collapse = "|"),
       eta_squared = eta2, delta_top_bottom_sd = delta_sd,
       top_group = grp$hap[1], top_n = grp$n[1], top_mean = grp$mean[1],
       bottom_group = grp$hap[nrow(grp)], bottom_n = grp$n[nrow(grp)],
       bottom_mean = grp$mean[nrow(grp)],
       groups = grp)
}

## =============================================================================
## Main loop -- one crosshap run per gene.
## =============================================================================
## Raw per-gene SNP counts, so a not-tested gene can still report how much data it had.
raw_manifest <- file.path(cfg$output_root, "03_per_gene_vcfs", "raw_1000bp_manifest.tsv")
raw_snp_count <- if (file.exists(raw_manifest)) {
  m <- fread(raw_manifest); setNames(as.list(m$n_snps), m$gene_id)
} else list()

targets <- as.data.table(read_gene_windows(gw_path))
setorder(targets, trait, gene_file)
if (nrow(targets) == 0) stop("No targets in ", gw_path)
log_msg(sprintf("Targets: %d genes from %s", nrow(targets), gw_path),
        stage="discovery", event="targets_loaded", event_log_path=event_log_path)

failures    <- init_failures_df()
results     <- list()
groups_out  <- list()
not_tested  <- list()
label       <- paste0("Haplotypes_MGmin", MGMIN, "_E", EPS)

for (i in seq_len(nrow(targets))) {
  tr <- targets$trait[i]; gene_file <- targets$gene_file[i]
  title <- targets$title[i]; gene_name <- make_gene_name(gene_file)
  cat(sprintf("[%2d/%d] %-11s %s\n", i, nrow(targets), tr, gene_name))
  log_msg(sprintf("Processing gene %d/%d", i, nrow(targets)), trait=tr, gene_file=gene_file,
          stage="gene", event="gene_start", event_log_path=event_log_path)

  cache_dir <- file.path(cache_root, tr, gene_name, paste0("MGmin_", MGMIN)); ensure_dir(cache_dir)
  cache_rds <- file.path(cache_dir, "HapObject.rds")
  res <- NULL

  if (file.exists(cache_rds) && !isTRUE(cfg$overwrite_cache)) {
    res <- tryCatch(readRDS(cache_rds), error = function(e) NULL)
    if (!is.null(res)) log_msg("Using cached HapObject.", trait=tr, gene_file=gene_file,
                               MGmin=MGMIN, stage="cache", event="cache_used",
                               event_log_path=event_log_path)
  }
  crosshap_err <- NA_character_
  if (is.null(res)) {
    res <- tryCatch(run_crosshap(cfg=cfg, trait=tr, gene_file=gene_file, MGmin=MGMIN, tmp_dir=tmp_root),
      error = function(e) { msg <- conditionMessage(e); crosshap_err <<- msg
        failures <<- add_failure(failures, tr, gene_file, MGMIN, as.character(EPS), "crosshap",
                                 "crosshap_failed", msg)
        log_msg(paste0("CrossHap failed: ", msg), level="ERROR", trait=tr, gene_file=gene_file,
                MGmin=MGMIN, stage="crosshap", event="failure_encountered",
                event_log_path=event_log_path); NULL })
    if (!is.null(res)) saveRDS(res, cache_rds)
  }

  if (is.null(res)) {
    # Classify the failure rather than lumping everything together. crosshap raises
    # "object 'Varfile' not found" when DBSCAN finds NO marker groups at this epsilon
    # (all SNPs are noise) -- that is a parameter-coverage outcome, not a data problem,
    # and it is the dominant reason on a sparse gene set. Distinguish it from genuinely
    # having too few variants, so the two are countable separately in the paper.
    n_raw <- suppressWarnings(as.integer(raw_snp_count[[gene_name]]))
    reason <- if (!is.na(crosshap_err) && grepl("Varfile", crosshap_err, fixed = TRUE))
                sprintf("no_marker_groups_at_epsilon_%s", EPS)
              else if (!is.na(crosshap_err) && grepl("Too few common variants", crosshap_err))
                "fewer_than_MGmin_variants_after_raw_imputed_intersection"
              else if (!is.na(crosshap_err) && grepl("0 variants", crosshap_err))
                "no_variants_in_window"
              else paste0("crosshap_error: ", substr(crosshap_err %||% "unknown", 1, 120))
    not_tested[[length(not_tested)+1]] <- data.table(
      trait=tr, gene_id=gene_name, locus_id=targets$locus_id[i], lead_SNP=targets$lead_SNP[i],
      n_snps_window=if (length(n_raw)) n_raw else NA_integer_, reason=reason)
    next
  }

  n_common <- as.integer(res$n_common %||% NA_integer_)

  ## ---- Figures (unchanged renderers) ----
  comb_dir <- file.path(combined_root, tr, gene_name); ensure_dir(comb_dir)
  combined_pdf <- file.path(comb_dir, sprintf("CrosshapTree+Violin__%s__MGmin%d__Eps%s.pdf",
                                              gene_name, MGMIN, eps_tag(EPS)))
  tryCatch(write_combined_pdf(HapObject=res$HapObject, out_pdf=combined_pdf, title=title,
      trait=tr, gene_file=gene_file, gene_name=gene_name, MGmin=MGMIN,
      epsilon_vector=EPS, mgmin_test_stats=NULL, gene_summary_row=NULL),
    error=function(e) { failures <<- add_failure(failures, tr, gene_file, MGMIN, as.character(EPS),
      "combined_pdf", "combined_pdf_failed", conditionMessage(e)) })

  hm_dir <- file.path(heatmaps_root, tr, gene_name); ensure_dir(hm_dir)
  heatmap_pdf <- file.path(hm_dir, sprintf("Heatmaps__%s__MGmin%d__Eps%s.pdf",
                                           gene_name, MGMIN, eps_tag(EPS)))
  if (!is.null(res$HapObject[[label]])) {
    tryCatch(write_heatmaps_for_eps(HapObject=res$HapObject, label=label, gene_file=gene_file,
        title=title, trait=tr, eps=EPS, MGmin=MGMIN, raw_path=res$raw_path,
        common_ids=res$common_ids, vcf_raw_ids=res$vcf_raw_ids, out_pdf=heatmap_pdf),
      error=function(e) { failures <<- add_failure(failures, tr, gene_file, MGMIN, as.character(EPS),
        "heatmap", "heatmap_failed", conditionMessage(e)) })
  }

  ## ---- The single test ----
  if (is.null(res$HapObject[[label]])) {
    failures <- add_failure(failures, tr, gene_file, MGMIN, as.character(EPS), "epsilon_validation",
                            "no_valid_haplotype_grouping", "No HapObject entry at this epsilon.")
    not_tested[[length(not_tested)+1]] <- data.table(trait=tr, gene_id=gene_name,
      locus_id=targets$locus_id[i], lead_SNP=targets$lead_SNP[i], n_snps_window=n_common,
      reason="no_haplotype_grouping_at_this_epsilon")
    next
  }
  kw <- kw_test_with_effects(res$HapObject, label)
  if (!isTRUE(kw$valid)) {
    failures <- add_failure(failures, tr, gene_file, MGMIN, as.character(EPS), "kw_validation",
                            "no_valid_kw_pvalue", kw$reason)
    not_tested[[length(not_tested)+1]] <- data.table(trait=tr, gene_id=gene_name,
      locus_id=targets$locus_id[i], lead_SNP=targets$lead_SNP[i], n_snps_window=n_common,
      reason=kw$reason)
    log_msg(paste0("Not testable: ", kw$reason), level="WARN", trait=tr, gene_file=gene_file,
            stage="kw_validation", event="failure_encountered", event_log_path=event_log_path)
    next
  }

  results[[length(results)+1]] <- data.table(
    trait=tr, gene_id=gene_name, locus_id=targets$locus_id[i], class=targets$class[i],
    chr=targets$chr[i], gene_start=as.integer(targets$gene_start[i]),
    gene_end=as.integer(targets$gene_end[i]), strand=targets$strand[i],
    lead_SNP=targets$lead_SNP[i], lead_pos=as.integer(targets$lead_pos[i]),
    lead_neg_log10p=as.numeric(targets$lead_neg_log10p[i]),
    dist_to_lead_bp=as.integer(targets$dist_to_lead_bp[i]),
    n_snps_window=n_common, MGmin=MGMIN, epsilon=EPS,
    n_ind=kw$n_ind, n_unassigned=kw$n_unassigned, n_groups=kw$n_groups,
    group_sizes=kw$group_sizes, kw_H=kw$kw_H, kw_df=kw$df, kw_p_raw=kw$kw_p,
    eta_squared=kw$eta_squared, delta_top_bottom_sd=kw$delta_top_bottom_sd,
    top_group=kw$top_group, top_n=kw$top_n, top_mean=kw$top_mean,
    bottom_group=kw$bottom_group, bottom_n=kw$bottom_n, bottom_mean=kw$bottom_mean,
    description=as.character(targets$description[i]),
    combined_pdf_path=combined_pdf, heatmap_pdf_path=heatmap_pdf)

  g <- copy(kw$groups); g[, `:=`(trait=tr, gene_id=gene_name)]
  setnames(g, "hap", "haplotype_group")
  groups_out[[length(groups_out)+1]] <- g[, .(trait, gene_id, haplotype_group,
                                              n, mean_pheno=mean, median_pheno=median, sd_pheno=sd)]

  log_msg(sprintf("KW p=%s groups=%d eta2=%.3f", format_p_plain(kw$kw_p), kw$n_groups, kw$eta_squared),
          trait=tr, gene_file=gene_file, MGmin=MGMIN, epsilon=as.character(EPS),
          stage="kw", event="kw_test_done", event_log_path=event_log_path)
}

## =============================================================================
## Per-trait BH on the RAW p-values. Genes that produced no valid grouping are NOT
## in this denominator -- no test was performed on them -- but they are reported in
## genes_not_tested.tsv so the real search space stays visible.
## =============================================================================
res_dt <- if (length(results)) rbindlist(results) else data.table()
if (nrow(res_dt)) {
  res_dt[, fdr_q := p.adjust(kw_p_raw, method = cfg$stats$across_gene_fdr_method), by = trait]
  res_dt[, significant_fdr := fdr_q <= cfg$stats$alpha]
  res_dt[, n_tests_in_trait_family := .N, by = trait]
  setorder(res_dt, trait, kw_p_raw)
  for (tn in unique(res_dt$trait))
    log_msg(sprintf("BH done: %d tests, %d significant at alpha=%s",
                    sum(res_dt$trait==tn), sum(res_dt$significant_fdr[res_dt$trait==tn]),
                    cfg$stats$alpha), trait=tn, stage="trait_correction",
            event="trait_level_correction_complete", event_log_path=event_log_path)
}

nt_dt <- if (length(not_tested)) rbindlist(not_tested) else data.table()
gr_dt <- if (length(groups_out)) rbindlist(groups_out) else data.table()

## ---- locus_summary.tsv: report loci, not only genes ------------------------
loc_dt <- if (nrow(res_dt)) res_dt[, .(
    chr = chr[1], lead_SNP = lead_SNP[1], class = class[1],
    lead_neg_log10p = lead_neg_log10p[1],
    n_genes_tested = .N, n_genes_significant = sum(significant_fdr),
    best_gene = gene_id[which.min(kw_p_raw)], best_kw_p = min(kw_p_raw),
    best_fdr_q = fdr_q[which.min(kw_p_raw)],
    max_eta_squared = max(eta_squared, na.rm = TRUE)
  ), by = .(trait, locus_id)][order(trait, locus_id)] else data.table()
if (nrow(nt_dt) && nrow(loc_dt)) {
  nt_by_loc <- nt_dt[, .(n_genes_not_tested = .N), by = .(trait, locus_id)]
  loc_dt <- merge(loc_dt, nt_by_loc, by = c("trait","locus_id"), all.x = TRUE)
  loc_dt[is.na(n_genes_not_tested), n_genes_not_tested := 0L]
}

## ---- per_trait_summary.tsv -------------------------------------------------
tr_dt <- data.table()
if (nrow(res_dt)) {
  tr_dt <- res_dt[, .(n_genes_tested = .N, n_significant_fdr = sum(significant_fdr),
                      n_loci_tested = uniqueN(locus_id),
                      n_loci_with_significant = uniqueN(locus_id[significant_fdr]),
                      median_eta_squared = round(median(eta_squared, na.rm=TRUE), 4),
                      min_fdr_q = min(fdr_q)), by = trait]
  if (nrow(nt_dt)) {
    tr_dt <- merge(tr_dt, nt_dt[, .(n_genes_not_tested = .N), by = trait], by = "trait", all.x = TRUE)
    tr_dt[is.na(n_genes_not_tested), n_genes_not_tested := 0L]
  }
  setorder(tr_dt, trait)
}

## ---- Write ------------------------------------------------------------------
w <- function(d, f) if (nrow(d)) fwrite(d, file.path(stats_root, f), sep = "\t", na = "NA")
w(res_dt, "gene_results.tsv")
w(gr_dt,  "haplotype_groups.tsv")
w(nt_dt,  "genes_not_tested.tsv")
w(loc_dt, "locus_summary.tsv")
w(tr_dt,  "per_trait_summary.tsv")
if (nrow(failures)) fwrite(failures, file.path(stats_root, "failures_notes.tsv"), sep = "\t", na = "NA")

cat("\n================ SUMMARY ================\n")
cat(sprintf("Genes in         : %d\n", nrow(targets)))
cat(sprintf("Genes tested     : %d\n", nrow(res_dt)))
cat(sprintf("Genes not tested : %d\n", nrow(nt_dt)))
if (nrow(res_dt)) {
  cat(sprintf("Significant (BH q <= %s): %d genes in %d loci\n", cfg$stats$alpha,
              sum(res_dt$significant_fdr), uniqueN(res_dt$locus_id[res_dt$significant_fdr])))
  cat("\nPer trait:\n"); print(tr_dt, row.names = FALSE)
}
cat(sprintf("\nStats written to %s\n", stats_root))
log_msg("Run finished.", stage="run", event="run_end", event_log_path=event_log_path)
