args <- commandArgs(trailingOnly = TRUE)

results_dir <- if (length(args) >= 1) args[[1]] else "gwas_barley/morexV3_analysis/USED_data/results_maf0.000001_geno1_thin0_pc5"
out_file <- if (length(args) >= 2) args[[2]] else "gwas_barley/morexV3_analysis/FDR_check_for_manahatan_plots/FDR_diagnostics/combined_fdr_summary.csv"
bim_file <- file.path(results_dir, "snp_map.bim")

ps_files <- sort(list.files(results_dir, pattern = "^morexV3_.*_gwas\\.ps$", full.names = TRUE))
if (length(ps_files) == 0) {
  stop(sprintf("No GWAS .ps files found in %s", results_dir))
}

bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)

extract_trait_name <- function(path) {
  sub("_gwas\\.ps$", "", sub("^morexV3_", "", basename(path)))
}

build_summary_row <- function(ps_file) {
  gwas <- read.table(
    ps_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    col.names = c("SNP_ID_DUMMY", "BETA", "SE", "P")
  )

  if (nrow(gwas) != nrow(bim)) {
    stop(sprintf(
      "Row mismatch for %s: GWAS=%d BIM=%d",
      basename(ps_file),
      nrow(gwas),
      nrow(bim)
    ))
  }

  pvals <- suppressWarnings(as.numeric(gwas$P))
  valid <- is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1
  pvals <- pvals[valid]
  qvals <- p.adjust(pvals, method = "BH")
  m <- length(pvals)

  bh05 <- if (any(qvals <= 0.05)) max(pvals[qvals <= 0.05]) else NA_real_
  bh10 <- if (any(qvals <= 0.10)) max(pvals[qvals <= 0.10]) else NA_real_

  data.frame(
    trait_name = extract_trait_name(ps_file),
    total_snps = nrow(gwas),
    valid_snps = m,
    min_raw_p_value = min(pvals),
    max_minus_log10_p = max(-log10(pvals)),
    min_bh_q_value = min(qvals),
    n_q_le_0_05 = sum(qvals <= 0.05),
    n_q_le_0_10 = sum(qvals <= 0.10),
    n_bonferroni = sum(pvals <= (0.05 / m)),
    n_bonferroni_0_10 = sum(pvals <= (0.10 / m)),
    bonferroni_0_05_cutoff_raw_p = 0.05 / m,
    bonferroni_0_10_cutoff_raw_p = 0.10 / m,
    bh_q_0_05_cutoff_raw_p = bh05,
    bh_q_0_05_cutoff_minus_log10_p = if (is.na(bh05)) NA_real_ else -log10(bh05),
    bh_q_0_10_cutoff_raw_p = bh10,
    bh_q_0_10_cutoff_minus_log10_p = if (is.na(bh10)) NA_real_ else -log10(bh10),
    lambda_gc = median(qchisq(1 - pvals, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, lapply(ps_files, build_summary_row))
write.csv(summary_df, out_file, row.names = FALSE, quote = TRUE)
cat(sprintf("Wrote %s with %d rows\n", out_file, nrow(summary_df)))
print(summary_df[, c("trait_name", "n_q_le_0_05", "n_q_le_0_10", "n_bonferroni")])
