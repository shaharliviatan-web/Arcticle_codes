#!/usr/bin/env Rscript
# 06_make_annotation_review.R
# Build the per-trait gene annotation review (run AFTER the pipeline + 04_make_shortlist.R).
#
# For each trait: genes that PASSED the statistical filter (significant_fdr) first, sorted most
# significant -> least; then a clearly separated block of genes that did NOT pass (tested but not
# significant), then genes with no valid haplotype test. Each gene shows its name/annotation.
# Annotation source priority: (1) curated lookup from the previous tool (prev_annotations_lookup.tsv,
# "Relevance to ..." already stripped), (2) Morex V3 GFF description, (3) none -> NEED ANNOTATION.
#
# Output (05_results/tables/):
#   gene_annotation_review.csv   machine-readable; genes lacking an annotation show the literal
#                                "NEED TO ADD ANNOTATION" in the annotation column.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)

stats_root <- file.path(cfg$output_root, "04_runs", cfg$run_id, "Stats")
tab_root   <- file.path(cfg$output_root, "05_results", "tables")
dir.create(tab_root, showWarnings = FALSE, recursive = TRUE)

sl_path  <- file.path(stats_root, "gene_shortlist.csv")
lk_path  <- file.path(tab_root, "prev_annotations_lookup.tsv")
stopifnot(file.exists(sl_path))
sl <- fread(sl_path)
lk <- if (file.exists(lk_path)) fread(lk_path, header = FALSE, col.names = c("gene_id","annotation")) else
  data.table(gene_id = character(), annotation = character())
lk_map <- setNames(lk$annotation, lk$gene_id)

# genes whose curated note differed across traits in the source (best one kept) -> verify
ver_path <- file.path(tab_root, "prev_annotations_needs_verification.txt")
ver_ids <- if (file.exists(ver_path)) readLines(ver_path) else character(0)
ver_ids <- ver_ids[nzchar(ver_ids)]

# ---- Resolve annotation + source ----
gff_desc <- sl$description
sl[, curated := lk_map[gene_id]]
sl[, annotation := fifelse(!is.na(curated) & nzchar(curated), curated,
                    fifelse(!is.na(gff_desc) & gff_desc != "NA" & nzchar(gff_desc), gff_desc, NA_character_))]
sl[, annotation_source := fifelse(!is.na(curated) & nzchar(curated), "curated",
                           fifelse(!is.na(gff_desc) & gff_desc != "NA" & nzchar(gff_desc), "gff", "none"))]
sl[, annotation_status := fifelse(annotation_source == "none", "NEED_ANNOTATION", "HAVE")]
# Make the gap explicit in the CSV itself (no colour available in CSV)
sl[annotation_source == "none", annotation := "NEED TO ADD ANNOTATION"]
# Flag genes whose source note differed across traits (best kept; please verify)
sl[, needs_verification := gene_id %in% ver_ids]

# ---- Pass status + ordering ----
sl[, pass_group := fifelse(significant_fdr %in% TRUE, "1_PASSED_FDR",
                    fifelse(!is.na(best_kw_p_raw), "2_did_not_pass", "3_no_valid_test"))]
# Combined filter status: FDR is the primary gene filter; Bonferroni (0.05) shown alongside.
sl[, filter_status := fifelse(pass_group == "3_no_valid_test", "no_valid_test",
                       fifelse(significant_fdr %in% TRUE & significant_bonferroni %in% TRUE, "PASS_FDR_and_BONF",
                       fifelse(significant_fdr %in% TRUE, "PASS_FDR_only", "not_significant")))]
sl[, sortp := fifelse(is.na(trait_fdr_p), fifelse(is.na(best_kw_p_raw), Inf, best_kw_p_raw), trait_fdr_p)]
setorder(sl, trait, pass_group, sortp, -lead_neg_log10p)
sl[, rank_in_block := rowid(trait, pass_group)]

# ---- CSV ----
csv_cols <- c("trait","pass_group","filter_status","rank_in_block","gene_id","lead_SNP","class",
              "dist_to_lead_bp","n_snps_window","best_MGmin","best_epsilon","best_n_groups",
              "best_group_sizes","best_kw_p_raw","best_kw_p_holm_within_gene","trait_fdr_p",
              "trait_bonferroni_p","significant_fdr","significant_bonferroni",
              "annotation_status","annotation_source","needs_verification","annotation",
              "combined_pdf_path","heatmap_dir")
fwrite(sl[, ..csv_cols], file.path(tab_root, "gene_annotation_review.csv"))


# ---- Console summary ----
cat(sprintf("Wrote gene_annotation_review.csv (%d genes)\n", nrow(sl)))
cat("Per trait: passed / need-annotation:\n")
print(sl[, .(genes = .N,
             passed_fdr = sum(significant_fdr %in% TRUE),
             passed_bonferroni = sum(significant_bonferroni %in% TRUE),
             need_annotation = sum(annotation_status == "NEED_ANNOTATION")), by = trait])
