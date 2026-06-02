#!/usr/bin/env Rscript
# 06b_make_readable_summary.R
# Compact, human-readable per-gene view derived from gene_annotation_review.csv.
# One row per gene (per trait): trait, gene, passed FDR (and Bonferroni), the p-value,
# and the annotation ("NO ANNOTATION" when none). Sorted by trait, then p (most
# significant first; untested genes last).
#
# Output: 05_results/tables/genes_readable_summary.csv

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)
tab <- file.path(cfg$output_root, "05_results", "tables")

ar <- fread(file.path(tab, "gene_annotation_review.csv"))

fmtp <- function(p) ifelse(is.na(p), "", formatC(as.numeric(p), format = "g", digits = 3))

out <- data.table(
  trait            = ar$trait,
  gene             = ar$gene_id,
  passed_FDR       = fifelse(ar$significant_fdr %in% TRUE, "YES", "NO"),
  passed_Bonferroni= fifelse(ar$significant_bonferroni %in% TRUE, "YES", "NO"),
  KW_p_raw         = fmtp(ar$best_kw_p_raw),
  FDR_p            = fmtp(ar$trait_fdr_p),
  Bonferroni_p     = fmtp(ar$trait_bonferroni_p),
  tested           = fifelse(is.na(ar$best_kw_p_raw), "no_valid_test", "tested"),
  annotation       = fifelse(ar$annotation_status == "NEED_ANNOTATION" | is.na(ar$annotation) |
                             ar$annotation == "NEED TO ADD ANNOTATION", "NO ANNOTATION", ar$annotation)
)

# sort: trait, then FDR p ascending (most significant first), NA/untested last
out[, .sortp := fifelse(is.na(ar$trait_fdr_p),
                        fifelse(is.na(ar$best_kw_p_raw), Inf, as.numeric(ar$best_kw_p_raw)),
                        as.numeric(ar$trait_fdr_p))]
setorder(out, trait, .sortp)
out[, .sortp := NULL]

fwrite(out, file.path(tab, "genes_readable_summary.csv"))
cat(sprintf("Wrote genes_readable_summary.csv (%d rows)\n", nrow(out)))
cat("Per trait: passed_FDR / total:\n")
print(out[, .(genes = .N, passed_FDR = sum(passed_FDR == "YES"),
              passed_Bonferroni = sum(passed_Bonferroni == "YES"),
              no_annotation = sum(annotation == "NO ANNOTATION")), by = trait])
