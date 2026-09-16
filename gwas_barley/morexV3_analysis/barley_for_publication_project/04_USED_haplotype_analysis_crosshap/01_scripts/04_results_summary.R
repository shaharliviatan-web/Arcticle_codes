#!/usr/bin/env Rscript
# =============================================================================
# 04_results_summary.R
#
# WHAT   Turn the Stats tables into paste-ready prose + a parameters table for the
#        manuscript. Reads only; computes nothing new.
#
# READS  04_runs/<run_id>/Stats/{gene_results,genes_not_tested,locus_summary,
#                                per_trait_summary}.tsv
# WRITES 04_runs/<run_id>/Stats/analysis_parameters.tsv
#        04_runs/<run_id>/Stats/results_chapter_numbers.txt
#
# Replaces v1's 04_make_shortlist.R -- gene_results.tsv already IS the decision
# sheet, so the separate shortlist was pure duplication.
# =============================================================================
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)
S <- file.path(cfg$output_root, "04_runs", cfg$run_id, "Stats")

gr <- fread(file.path(S, "gene_results.tsv"))
nt <- if (file.exists(file.path(S,"genes_not_tested.tsv"))) fread(file.path(S,"genes_not_tested.tsv")) else data.table()
lo <- fread(file.path(S, "locus_summary.tsv"))
pt <- fread(file.path(S, "per_trait_summary.tsv"))
alpha <- cfg$stats$alpha
EPS <- cfg$epsilon_vector[[1]]; MG <- cfg$mgmin_values[[1]]

## ---- parameters table -------------------------------------------------------
p <- data.table(parameter = c(
  "run_id","candidate_genes_source","n_genes_in","per_gene_window_bp",
  "haplotyping_tool","LD_input","haplotyping_input","MGmin_dbscan_minPts",
  "epsilon_dbscan_eps","epsilon_selection_basis","minHap","hetmiss_as","keep_outliers",
  "test","test_filtering","effect_size","multiple_testing","family","alpha",
  "not_tested_handling","n_genes_tested","n_genes_not_tested",
  "n_significant_genes","n_loci_with_significant_gene","date_run"),
  value = c(
  cfg$run_id, basename(cfg$candidate_genes_tsv), as.character(nrow(gr)+nrow(nt)),
  as.character(cfg$window_bp),
  "crosshap (DBSCAN on a PLINK --r2 square matrix)",
  "imputed per-gene VCF (complete r2 matrix)",
  "raw per-gene VCF (observed genotypes)",
  as.character(MG), as.character(EPS),
  "genotype-only (gene coverage + accession assignment rate); the phenotype is never used to build haplotype groups",
  as.character(cfg$minHap), as.character(cfg$hetmiss_as), as.character(cfg$keep_outliers),
  "Kruskal-Wallis, phenotype (BLUP) ~ haplotype group, one test per gene",
  "drop hap 0 (unassigned) and missing phenotypes; require >= 2 groups",
  "eta-squared = (H - k + 1)/(n - k); plus top-vs-bottom group difference in phenotype SD",
  "Benjamini-Hochberg on RAW p-values", "per trait", as.character(alpha),
  "genes with no valid grouping are excluded from the BH denominator and listed in genes_not_tested.tsv",
  as.character(nrow(gr)), as.character(nrow(nt)),
  as.character(sum(gr$significant_fdr)),
  as.character(uniqueN(gr$locus_id[gr$significant_fdr])),
  format(Sys.Date(), "%Y-%m-%d")))
fwrite(p, file.path(S, "analysis_parameters.tsv"), sep = "\t")

## ---- prose ------------------------------------------------------------------
con <- file(file.path(S, "results_chapter_numbers.txt"), "w"); wl <- function(...) writeLines(paste0(...), con)
wl("================================================================")
wl("  Haplotype analysis at the candidate genes -- numbers for the manuscript")
wl("  run_id: ", cfg$run_id)
wl("  generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " by 04_results_summary.R")
wl("================================================================")
wl("")
wl("METHOD")
wl(sprintf("  crosshap haplotyping per gene (gene +/- %s bp), MGmin = %s, epsilon = %s, both FIXED.",
           cfg$window_bp, MG, EPS))
wl("  LD comes from the imputed per-gene VCF; haplotypes are called on the raw (observed)")
wl("  genotypes; the two variant sets are intersected on CHROM:POS.")
wl("  One Kruskal-Wallis test per gene: BLUP phenotype ~ haplotype group, dropping")
wl("  unassigned accessions (hap 0) and missing phenotypes, requiring >= 2 groups.")
wl(sprintf("  Benjamini-Hochberg FDR on the RAW p-values, within each trait, alpha = %s.", alpha))
wl("")
wl("HEADLINE NUMBERS")
wl(sprintf("  Candidate genes in       : %d", nrow(gr) + nrow(nt)))
wl(sprintf("  Genes tested             : %d", nrow(gr)))
wl(sprintf("  Genes not testable       : %d", nrow(nt)))
wl(sprintf("  Significant (BH q <= %s) : %d genes, in %d loci",
           alpha, sum(gr$significant_fdr), uniqueN(gr$locus_id[gr$significant_fdr])))
wl(sprintf("  Median haplotype groups  : %g", median(gr$n_groups)))
wl(sprintf("  Median eta-squared       : %.3f (significant genes: %.3f)",
           median(gr$eta_squared, na.rm=TRUE),
           median(gr$eta_squared[gr$significant_fdr], na.rm=TRUE)))
wl("")
wl("  NOTE ON COUNTING: significant genes cluster inside loci -- genes in one locus are")
wl("  in LD and are not independent discoveries. Report loci alongside genes.")
wl("")
if (nrow(nt)) {
  wl("WHY GENES WERE NOT TESTABLE")
  for (r in nt[, .N, by = reason][order(-N)]$reason)
    wl(sprintf("  %-58s %d", substr(r,1,58), nt[reason == r, .N]))
  wl("")
}
wl("PER TRAIT")
for (i in seq_len(nrow(pt))) with(pt[i,], wl(sprintf(
  "  %-11s %2d tested / %2d not testable | %2d significant in %d of %d loci | median eta2 %.3f",
  trait, n_genes_tested, n_genes_not_tested, n_significant_fdr,
  n_loci_with_significant, n_loci_tested, median_eta_squared)))
wl("")
wl("SIGNIFICANT GENES (BH q <= ", alpha, "), strongest first")
s <- gr[significant_fdr == TRUE][order(kw_p_raw)]
for (i in seq_len(nrow(s))) with(s[i,], wl(sprintf(
  "  %-11s %-28s %-14s q=%-9.2g eta2=%.3f  %d groups  top-vs-bottom %+.2f SD  (%s n=%d vs %s n=%d)",
  trait, gene_id, locus_id, fdr_q, eta_squared, n_groups, delta_top_bottom_sd,
  top_group, top_n, bottom_group, bottom_n)))
wl("")
wl("PER LOCUS")
for (tr in sort(unique(lo$trait))) {
  wl(sprintf("  --- %s ---", tr))
  for (i in which(lo$trait == tr)) with(lo[i,], wl(sprintf(
    "  %-14s %s lead %-14s  %d/%d genes significant  best %s (q=%.2g, eta2=%.3f)",
    locus_id, chr, lead_SNP, n_genes_significant, n_genes_tested,
    sub("HORVU.MOREX.r3.","",best_gene), best_fdr_q, max_eta_squared)))
}
close(con)

cat(sprintf("Wrote analysis_parameters.tsv (%d rows) and results_chapter_numbers.txt\n", nrow(p)))
cat(sprintf("%d genes tested, %d significant in %d loci\n",
            nrow(gr), sum(gr$significant_fdr), uniqueN(gr$locus_id[gr$significant_fdr])))
