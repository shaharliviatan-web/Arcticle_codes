#!/usr/bin/env Rscript
# 04_make_shortlist.R
# Build the single decision sheet for picking significant genes + epsilon + MGmin.
#
# One row per candidate gene (ALL genes in gene_windows.tsv, tested or not), with:
#   - origin: lead_SNP, class (significant|marginal), lead_neg_log10p, dist_to_lead_bp, locus_id
#   - n_snps_window
#   - status: SIGNIFICANT_FDR / tested_ns / no_valid_test
#   - best test: MGmin, epsilon, n_groups, group_sizes, raw p, within-gene Holm p
#   - trait-level FDR + Bonferroni p and significance flags
#   - description, and paths to the best combined PDF + heatmap dir
# Sorted by trait, then significant first, then trait_fdr_p ascending.
#
# Output: 04_runs/<run_id>/Stats/gene_shortlist.csv

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)

run_root   <- file.path(cfg$output_root, "04_runs", cfg$run_id)
stats_root <- file.path(run_root, "Stats")
gw_path    <- file.path(cfg$output_root, "00_config", "gene_windows.tsv")
raw_manifest <- file.path(cfg$output_root, "03_per_gene_vcfs", "raw_1000bp_manifest.tsv")

gw <- fread(gw_path)
gw[, gene_file := paste0(gene_id, ".vcf.gz")]

gs_path <- file.path(stats_root, "gene_summary.csv")
pt_path <- file.path(stats_root, "per_test_stats.csv")
fa_path <- file.path(stats_root, "failures_notes.csv")
gs <- if (file.exists(gs_path)) fread(gs_path) else data.table()
pt <- if (file.exists(pt_path)) fread(pt_path) else data.table()
fa <- if (file.exists(fa_path)) fread(fa_path) else data.table()
man <- if (file.exists(raw_manifest)) fread(raw_manifest) else data.table(gene_id=character(), n_snps=integer())

# n_snps in window from the raw manifest
nsnp <- setNames(man$n_snps, man$gene_id)

# best combined PDF + best-epsilon heatmap path (from the is_best_within_gene row)
best_pt <- if (nrow(pt) > 0) pt[is_best_within_gene == TRUE,
                               .(combined_pdf_path = combined_pdf_path[1],
                                 heatmap_pdf_path = heatmap_pdf_path[1]),
                               by = .(trait, gene_file)] else
  data.table(trait=character(), gene_file=character(),
             combined_pdf_path=character(), heatmap_pdf_path=character())

# per-gene failure reason (for untested genes): take the gene-level status if present
fa_reason <- if (nrow(fa) > 0) fa[, .(fail_status = paste(sort(unique(status)), collapse=";")),
                                  by = .(trait, gene_file)] else
  data.table(trait=character(), gene_file=character(), fail_status=character())

# ---- Assemble: start from all genes (gw), left-join gene_summary ----
gs_cols <- c("trait","gene_file","best_MGmin","best_epsilon","best_n_ind","best_n_groups",
             "best_group_sizes","best_kw_p_raw","best_kw_p_holm_within_gene",
             "trait_fdr_p","trait_bonferroni_p","significant_fdr","significant_bonferroni",
             "n_unique_valid_tests","n_total_valid_tests_before_collapse")
if (nrow(gs) > 0) gs_use <- gs[, ..gs_cols] else
  gs_use <- setNames(data.table(matrix(nrow=0, ncol=length(gs_cols))), gs_cols)

sl <- merge(gw[, .(trait, gene_id, gene_file, chr, gene_start, gene_end, strand,
                   lead_SNP, class, lead_neg_log10p, dist_to_lead_bp, locus_id, description)],
            gs_use, by = c("trait","gene_file"), all.x = TRUE, sort = FALSE)
sl <- merge(sl, best_pt, by = c("trait","gene_file"), all.x = TRUE, sort = FALSE)
sl <- merge(sl, fa_reason, by = c("trait","gene_file"), all.x = TRUE, sort = FALSE)

sl[, n_snps_window := as.integer(nsnp[gene_id])]
sl[is.na(n_snps_window), n_snps_window := 0L]

tested <- !is.na(sl$best_kw_p_raw)
sl[, status := fifelse(!tested, "no_valid_test",
                fifelse(significant_fdr %in% TRUE & significant_bonferroni %in% TRUE, "SIGNIFICANT_FDR_and_BONF",
                fifelse(significant_fdr %in% TRUE, "SIGNIFICANT_FDR_only", "tested_ns")))]
sl[!tested & (is.na(fail_status) | fail_status==""), fail_status := "no_valid_haplotype_grouping"]
sl[tested, fail_status := ""]

# heatmap dir for the gene (browse all epsilons)
sl[, heatmap_dir := file.path(run_root, "Heatmaps", trait, gene_id)]

# ---- Order columns + sort ----
setcolorder(sl, c("trait","gene_id","chr","gene_start","gene_end","strand",
                  "lead_SNP","class","lead_neg_log10p","dist_to_lead_bp","locus_id",
                  "n_snps_window","status","fail_status",
                  "best_MGmin","best_epsilon","best_n_ind","best_n_groups","best_group_sizes",
                  "best_kw_p_raw","best_kw_p_holm_within_gene","trait_fdr_p","trait_bonferroni_p",
                  "significant_fdr","significant_bonferroni",
                  "n_unique_valid_tests","n_total_valid_tests_before_collapse",
                  "description","combined_pdf_path","heatmap_pdf_path","heatmap_dir"))

# sort: trait, tested-before-untested, FDR p asc (NA last)
sl[, .sortp := fifelse(is.na(trait_fdr_p), Inf, trait_fdr_p)]
setorder(sl, trait, .sortp, -lead_neg_log10p)
sl[, .sortp := NULL]

out_path <- file.path(stats_root, "gene_shortlist.csv")
fwrite(sl, out_path)

# ---- Console summary ----
cat(sprintf("gene_shortlist.csv: %d genes\n", nrow(sl)))
cat("Per-trait status counts:\n")
print(sl[, .N, by = .(trait, status)][order(trait, status)])
cat("\nSignificant (FDR) genes per trait:\n")
print(sl[significant_fdr %in% TRUE, .N, by = trait])
cat(sprintf("\nWrote %s\n", out_path))
