#!/usr/bin/env Rscript
# 11_chosen_config_tables_v3.R
# Numeric / TSV package for the CHOSEN configuration.
#
#   Configuration (locked 2026-08-17, after reviewing the 24-cell sensitivity grid):
#     phenotype value   BLUP
#     PCs as covariates 3
#     alpha             0.10
#     Bonferroni        -log10p = 6.0454   (alpha / 111,017 LD-pruned SNPs)
#
# SCOPE (decision 2026-08-17):
#   * SIGNIFICANT SNPs ONLY -- no marginal / suggestive class.
#   * SNP LEVEL ONLY -- no LD clumping, no loci, no gene-search windows.
#     LD decay is to be recomputed per sub-ecosystem later; until that value is
#     settled, nothing here depends on a window, so nothing here goes stale.
#     Locus grouping and search windows will be added as a separate step.
#
# Outputs -> results/00_FINAL_BLUP_3PC/
#   tables/significant_snps.tsv       every SNP above the threshold, full stats
#   tables/per_trait_summary.tsv      per-trait counts, lambda, top signal
#   tables/analysis_parameters.tsv    every constant used, one row per parameter
#   tables/top15_per_chr__<trait>.tsv top 15 SNPs per chromosome
#   pc_selection/Table_S2_pc_variance.tsv
#   README.md                         all numbers in prose, for the manuscript
#
# Created 2026-08-17.
#
# BETA SIGN CONVENTION (fixed 2026-09-22): `beta` in significant_snps.tsv and `top_beta` in
# per_trait_summary.tsv are the effect of A1, the minor allele (PLINK A1). EMMAX's .ps beta is
# the effect of the allele coded "2" in the --recode12 tped (PLINK A2, the major allele), so it
# is sign-flipped below. Verification: see 36_paper_tables_loci.R. The top15_per_chr__* files
# are copied unchanged from results/tables/top_snps_per_chr/ and keep the raw EMMAX (A2) sign.
# CAUTION: this script also rewrites results/00_FINAL_BLUP_3PC/README.md with an old
# auto-generated text; that README has since been curated by hand. On 2026-09-22 the script was
# re-run with that one write suppressed (writeLines masked for that path).

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))

PIPE   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
INTER  <- file.path(PIPE, "intermediates")
PS_DIR <- file.path(PIPE, "results", "emmax_ps")
TAB_SRC<- file.path(PIPE, "results", "tables")
SET    <- file.path(PIPE, "results", "00_FINAL_BLUP_3PC")
OUT_T  <- file.path(SET, "00_snp_level"); OUT_P <- file.path(SET, "03_pc_selection")
for (d in c(OUT_T, OUT_P)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

PHENO <- "BLUP"; N_PCS <- 3L; PC_TAG <- "pc3"; ALPHA <- 0.10
PRUNE_IN <- file.path(INTER, "morexV3_pruned_for_covs.prune.in")
N_PRUNED <- length(readLines(PRUNE_IN))
BONF     <- -log10(ALPHA / N_PRUNED)
N_SAMPLES<- 290L
TRAITS   <- c("betaglucan", "fiber", "protein", "starch")
TRAIT_LAB<- c(betaglucan = "beta-glucan", fiber = "Fiber", protein = "Protein", starch = "Starch")
cat(sprintf("[11] %s x %d PCs | alpha=%.2f | N_pruned=%s | Bonf=%.4f | SNP level only\n",
            PHENO, N_PCS, ALPHA, format(N_PRUNED, big.mark=","), BONF))

snp_map <- fread(file.path(INTER, "morexV3_290.bim"), header = FALSE,
                 col.names = c("CHR_raw","SNP","cm","BP","A1","A2"))
snp_map[, CHR := as.integer(sub("H$","",CHR_raw))]
N_ALL <- nrow(snp_map)
maf <- fread(file.path(INTER, "morexV3_290_freq.frq"))[, .(SNP, MAF)]
lam <- fread(file.path(TAB_SRC, "lambda_table.tsv"))[pheno_type == PHENO & n_PCs == N_PCS,
              .(trait, lambda_GC)]
stopifnot(nrow(lam) == 4)

sig_l <- list(); summ_l <- list()
for (tr in TRAITS) {
  gw <- fread(file.path(PS_DIR, sprintf("morexV3__%s__%s__%s.ps", tr, PHENO, PC_TAG)),
              header = FALSE, col.names = c("SNP","beta","SE","P"), showProgress = FALSE)
  stopifnot(nrow(gw) == N_ALL)
  d <- data.table(CHR = snp_map$CHR, SNP = snp_map$SNP, BP = snp_map$BP,
                  A1 = snp_map$A1, A2 = snp_map$A2,
                  beta = -as.numeric(gw$beta),   # -beta: EMMAX A2 effect -> A1 effect
                  SE = as.numeric(gw$SE), P = as.numeric(gw$P))
  d <- d[!is.na(P) & P > 0 & P <= 1][, nlp := -log10(P)]
  d <- merge(d, maf, by = "SNP", all.x = TRUE)

  sig <- d[nlp >= BONF][order(-nlp)]
  top <- d[which.max(nlp)]
  sig_l[[tr]] <- sig[, .(trait = tr, SNP_id = SNP, chr = paste0(CHR,"H"), position_bp = BP,
                         A1, A2, MAF = round(MAF,4), beta = round(beta,6), SE = round(SE,6),
                         p_value = signif(P,4), neg_log10_p = round(nlp,4),
                         Bonf_threshold = round(BONF,4))]
  summ_l[[tr]] <- data.table(
    trait = tr, trait_label = TRAIT_LAB[[tr]],
    n_SNPs_tested = nrow(d), lambda_GC = lam[trait == tr, lambda_GC],
    Bonf_threshold_neglog10p = round(BONF,4), Bonf_threshold_pvalue = signif(ALPHA/N_PRUNED,4),
    n_significant_SNPs = nrow(sig),
    n_chromosomes_with_signal = uniqueN(sig$CHR),
    top_SNP = top$SNP, top_chr = paste0(top$CHR,"H"), top_pos = top$BP,
    top_neg_log10_p = round(top$nlp,4), top_p_value = signif(top$P,4),
    top_MAF = round(top$MAF,4), top_beta = round(top$beta,6))
  cat(sprintf("[11] %-11s significant SNPs=%3d  top=%.4f  lambda=%.4f\n",
              tr, nrow(sig), top$nlp, lam[trait == tr, lambda_GC]))
  rm(d, gw); gc(verbose = FALSE)
}

sig <- rbindlist(sig_l); summ <- rbindlist(summ_l)
fwrite(sig,  file.path(OUT_T, "significant_snps.tsv"), sep = "\t")
fwrite(summ, file.path(OUT_T, "per_trait_summary.tsv"), sep = "\t")

scree <- fread(file.path(INTER, "morexV3_pca_scree_data.tsv"))
params <- data.table(
  parameter = c("phenotype_value","n_PCs_covariates","alpha","n_SNPs_tested",
                "n_SNPs_LD_pruned","Bonferroni_neglog10p","Bonferroni_pvalue",
                "LD_pruning_call","n_samples","kinship","GWAS_software",
                "PC1_pct_var","PC1_3_cum_pct_var","lambda_GC_min","lambda_GC_max",
                "n_significant_SNPs_total"),
  value = c(PHENO, N_PCS, ALPHA, format(N_ALL, big.mark=","),
            format(N_PRUNED, big.mark=","), sprintf("%.4f", BONF),
            format(signif(ALPHA/N_PRUNED,4), scientific = TRUE),
            "--indep-pairwise 1000kb 1 0.2", N_SAMPLES,
            "EMMAX aIBS (290x290)", "EMMAX", sprintf("%.2f%%", scree$pct_variance[1]),
            sprintf("%.2f%%", scree$cum_pct[3]),
            sprintf("%.4f", min(summ$lambda_GC)), sprintf("%.4f", max(summ$lambda_GC)),
            sum(summ$n_significant_SNPs)))
fwrite(params, file.path(OUT_T, "analysis_parameters.tsv"), sep = "\t")

for (tr in TRAITS) {
  src <- file.path(TAB_SRC, "top_snps_per_chr",
                   sprintf("top15perChr__%s__%s__%s.tsv", tr, PHENO, PC_TAG))
  if (file.exists(src)) file.copy(src, file.path(OUT_T, sprintf("top15_per_chr__%s.tsv", tr)), overwrite = TRUE)
}
fwrite(scree, file.path(OUT_P, "Table_S2_pc_variance.tsv"), sep = "\t")

md <- c("# GWAS results - chosen configuration", "",
  sprintf("Generated %s. **SNP-level results only** - no LD clumping, no loci, no gene-search windows.", Sys.Date()), "",
  "## Configuration", "", "| parameter | value |", "|---|---|",
  sprintf("| Phenotype value | %s |", PHENO),
  sprintf("| PCs as EMMAX covariates | %d |", N_PCS),
  sprintf("| Kinship | EMMAX aIBS, %d x %d |", N_SAMPLES, N_SAMPLES),
  sprintf("| Samples | %d wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*), Southern Levant |", N_SAMPLES),
  sprintf("| SNPs tested | %s |", format(N_ALL, big.mark=",")),
  sprintf("| LD pruning (for PCA/kinship) | `--indep-pairwise 1000kb 1 0.2` -> %s SNPs |", format(N_PRUNED, big.mark=",")),
  sprintf("| Significance | Bonferroni, alpha = %.2f over %s independent tests |", ALPHA, format(N_PRUNED, big.mark=",")),
  sprintf("| Threshold | **-log10p = %.4f** (p < %s) |", BONF, format(signif(ALPHA/N_PRUNED,4), scientific=TRUE)),
  sprintf("| PC1 variance | %.2f%% |", scree$pct_variance[1]),
  sprintf("| PC1-PC3 cumulative | %.2f%% |", scree$cum_pct[3]),
  sprintf("| lambda_GC range | %.4f - %.4f |", min(summ$lambda_GC), max(summ$lambda_GC)), "",
  "**Significant SNPs only.** No marginal / suggestive class is used.", "",
  "## Results per trait", "",
  "| trait | lambda_GC | significant SNPs | chromosomes | top SNP | top -log10p | top p |",
  "|---|---|---|---|---|---|---|")
for (i in seq_len(nrow(summ))) md <- c(md, sprintf("| %s | %.4f | %d | %d | `%s` | %.4f | %s |",
  summ$trait_label[i], summ$lambda_GC[i], summ$n_significant_SNPs[i],
  summ$n_chromosomes_with_signal[i], summ$top_SNP[i], summ$top_neg_log10_p[i],
  format(summ$top_p_value[i], scientific = TRUE)))
md <- c(md, "", sprintf("**Total: %d significant SNPs.**", sum(summ$n_significant_SNPs)), "",
  "## Significant SNPs", "")
for (tr in TRAITS) {
  d <- sig[trait == tr]
  md <- c(md, sprintf("### %s (%d SNPs)", TRAIT_LAB[[tr]], nrow(d)), "",
    "| SNP | chr | position (bp) | A1/A2 | MAF | beta | SE | p | -log10p |",
    "|---|---|---|---|---|---|---|---|---|",
    sprintf("| `%s` | %s | %s | %s/%s | %.4f | %+.4f | %.4f | %s | %.4f |", d$SNP_id, d$chr,
            format(d$position_bp, big.mark=","), d$A1, d$A2, d$MAF, d$beta, d$SE,
            format(d$p_value, scientific=TRUE), d$neg_log10_p), "")
}
md <- c(md, "## Files", "", "| file | contents |", "|---|---|",
  "| `tables/significant_snps.tsv` | every SNP above the threshold, full stats |",
  "| `tables/per_trait_summary.tsv` | per-trait counts, lambda, top signal |",
  "| `tables/analysis_parameters.tsv` | every constant used in this analysis |",
  "| `tables/top15_per_chr__<trait>.tsv` | top 15 SNPs per chromosome |",
  "| `pc_selection/Table_S2_pc_variance.tsv` | PC variance table (PC1-PC20) |", "",
  "## Not included yet", "",
  "LD clumping into loci, gene-search windows around peaks, and Manhattan figures.",
  "These depend on the LD decay distance, which is being recomputed per sub-ecosystem.", "")
writeLines(md, file.path(SET, "README.md"))

cat("\n---------------------------------\n")
cat(sprintf("[11] CHECKPOINT: significant SNPs = %d\n", nrow(sig)))
cat(sprintf("[11] CHECKPOINT: files in tables/ = %d\n", length(list.files(OUT_T))))
print(summ[, .(trait, lambda_GC, n_significant_SNPs, top_SNP, top_neg_log10_p)])
cat(sprintf("[11] OK: %s\n", SET))
