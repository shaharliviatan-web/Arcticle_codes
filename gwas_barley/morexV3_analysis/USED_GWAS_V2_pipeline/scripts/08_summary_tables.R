#!/usr/bin/env Rscript
# 08_summary_tables.R
# Build the master 24-row summary table from the 24 EMMAX runs:
#   trait | pheno_type | n_PCs | lambda_GC | n_valid_pvals
#   n_SNPs_ge_Bonf_all | n_SNPs_ge_Bonf_pruned | n_SNPs_ge_FDR05 | n_SNPs_ge_FDR10
#   FDR05_threshold_neglog10p | FDR10_threshold_neglog10p
#   top_SNP_1H .. top_SNP_7H   (each: "SNPID:neglog10p")
# Plus a compact per-config rollup printed to stdout/log for the manual pick.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# --- Paths ---
PIPE       <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
INTER_DIR  <- file.path(PIPE, "intermediates")
PS_DIR     <- file.path(PIPE, "results", "emmax_ps")
TABLES_DIR <- file.path(PIPE, "results", "tables")
BIM_FILE   <- file.path(INTER_DIR, "snp_map.bim")
PRUNE_IN   <- file.path(INTER_DIR, "morexV3_pruned_for_covs.prune.in")
LAMBDA_TSV <- file.path(TABLES_DIR, "lambda_table.tsv")
OUT_TSV    <- file.path(TABLES_DIR, "summary_24runs.tsv")

# --- Load SNP map ONCE (CHR + SNP + BP; positional join with .ps files) ---
cat("[08] Loading snp_map.bim ...\n")
snp_map <- fread(BIM_FILE, header = FALSE,
                 select = c(1, 2, 4),
                 col.names = c("CHR_raw", "SNP", "BP"))
snp_map[, CHR := as.integer(sub("H$", "", sub("^chr", "", as.character(CHR_raw))))]
snp_map[is.na(CHR), CHR := as.integer(CHR_raw) - 30L]
stopifnot(all(snp_map$CHR %in% 1:7))

N_ALL    <- nrow(snp_map)
N_PRUNED <- length(readLines(PRUNE_IN))
BONF_ALL    <- -log10(0.05 / N_ALL)
BONF_PRUNED <- -log10(0.05 / N_PRUNED)
cat(sprintf("[08] N_all=%s  N_pruned=%s  Bonf_all=%.3f  Bonf_pruned=%.3f\n",
            format(N_ALL, big.mark = ","), format(N_PRUNED, big.mark = ","),
            BONF_ALL, BONF_PRUNED))

# --- lambda_table.tsv anchor (24 rows from step 06) ---
lambda_dt <- fread(LAMBDA_TSV)
stopifnot(nrow(lambda_dt) == 24)

# --- Discover + parse .ps files ---
ps_files <- list.files(PS_DIR, pattern = "^morexV3__.+__.+__pc\\d+\\.ps$",
                       full.names = TRUE)
stopifnot(length(ps_files) == 24)
ps_dt <- data.table(path = ps_files, base = sub("\\.ps$", "", basename(ps_files)))
ps_dt[, c("prefix", "trait", "pheno_type", "pcs_str") := tstrsplit(base, "__", fixed = TRUE)]
ps_dt[, n_PCs := as.integer(sub("^pc", "", pcs_str))]

# --- Per-cell summary ---
summarize_one <- function(i) {
  row <- ps_dt[i, ]
  gw <- fread(row$path, select = 4, col.names = "P", showProgress = FALSE)
  stopifnot(nrow(gw) == N_ALL)

  d <- data.table(CHR = snp_map$CHR, SNP = snp_map$SNP, P = gw$P)
  d <- d[!is.na(P) & P > 0 & P <= 1]
  d[, nlp := -log10(P)]

  n_bonf_all    <- d[nlp >= BONF_ALL,    .N]
  n_bonf_pruned <- d[nlp >= BONF_PRUNED, .N]

  q <- p.adjust(d$P, method = "BH")
  i05 <- which(q <= 0.05); i10 <- which(q <= 0.10)
  fdr05_thr <- if (length(i05) > 0) -log10(max(d$P[i05], na.rm = TRUE)) else NA_real_
  fdr10_thr <- if (length(i10) > 0) -log10(max(d$P[i10], na.rm = TRUE)) else NA_real_

  # Top SNP per chromosome (min P)
  top <- d[, .SD[which.min(P)], by = CHR][order(CHR)]
  top[, label := sprintf("%s:%.2f", SNP, nlp)]
  top_v <- setNames(rep(NA_character_, 7L), paste0("top_SNP_", 1:7, "H"))
  for (k in seq_len(nrow(top))) top_v[paste0("top_SNP_", top$CHR[k], "H")] <- top$label[k]

  c(list(trait = row$trait, pheno_type = row$pheno_type, n_PCs = row$n_PCs,
         n_SNPs_ge_Bonf_all = n_bonf_all, n_SNPs_ge_Bonf_pruned = n_bonf_pruned,
         n_SNPs_ge_FDR05 = length(i05), n_SNPs_ge_FDR10 = length(i10),
         FDR05_threshold_neglog10p = round(fdr05_thr, 4),
         FDR10_threshold_neglog10p = round(fdr10_thr, 4)),
    as.list(top_v))
}

cat(sprintf("[08] Summarizing %d cells (mclapply mc.cores=8)...\n", nrow(ps_dt)))
res <- mclapply(seq_len(nrow(ps_dt)), summarize_one, mc.cores = 8)

err <- which(vapply(res, inherits, logical(1), what = "try-error"))
if (length(err) > 0) {
  for (i in err) { cat(sprintf("[08] FAIL cell %d: %s\n", i, basename(ps_dt$path[i]))); print(res[[i]]) }
  stop(sprintf("[08] FAIL: %d cells errored", length(err)))
}

summary_dt <- rbindlist(res, fill = TRUE)
summary_dt <- merge(summary_dt,
                    lambda_dt[, .(trait, pheno_type, n_PCs, lambda_GC, n_valid_pvals)],
                    by = c("trait", "pheno_type", "n_PCs"), all.x = TRUE)

col_order <- c("trait", "pheno_type", "n_PCs", "lambda_GC", "n_valid_pvals",
               "n_SNPs_ge_Bonf_all", "n_SNPs_ge_Bonf_pruned",
               "n_SNPs_ge_FDR05", "n_SNPs_ge_FDR10",
               "FDR05_threshold_neglog10p", "FDR10_threshold_neglog10p",
               paste0("top_SNP_", 1:7, "H"))
setcolorder(summary_dt, col_order)
setorder(summary_dt, pheno_type, n_PCs, trait)

fwrite(summary_dt, OUT_TSV, sep = "\t", na = "NA")
cat(sprintf("\n[08] Wrote %s (%d rows x %d cols)\n", OUT_TSV, nrow(summary_dt), ncol(summary_dt)))

# --- Compact per-config rollup for manual pick ---
cat("\n[08] Per-config rollup (sorted by mean |lambda - 1|; for inspection only):\n")
cfg <- summary_dt[, .(
  mean_abs_lambda_minus_1 = round(mean(abs(lambda_GC - 1)), 4),
  min_lambda = round(min(lambda_GC), 4),
  max_lambda = round(max(lambda_GC), 4),
  hits_Bonf_all    = sum(n_SNPs_ge_Bonf_all),
  hits_Bonf_pruned = sum(n_SNPs_ge_Bonf_pruned),
  hits_FDR05       = sum(n_SNPs_ge_FDR05),
  hits_FDR10       = sum(n_SNPs_ge_FDR10)
), by = .(pheno_type, n_PCs)]
setorder(cfg, mean_abs_lambda_minus_1)
print(cfg)

# --- Checkpoints ---
cat("\n---------------------------------\n")
stopifnot(nrow(summary_dt) == 24, ncol(summary_dt) == length(col_order))
cat(sprintf("[08] CHECKPOINT: summary_24runs.tsv = %d rows x %d cols (expected 24 x %d)\n",
            nrow(summary_dt), ncol(summary_dt), length(col_order)))
cat("[08] OK\n")
