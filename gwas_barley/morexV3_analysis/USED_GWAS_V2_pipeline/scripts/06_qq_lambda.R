#!/usr/bin/env Rscript
# 06_qq_lambda.R
# For each of the 24 EMMAX .ps files:
#   - compute genomic inflation lambda (median-chi-square method)
#   - generate a QQ plot in BOTH vector PDF and 300dpi PNG (TAG spec)
# Aggregates a single results/tables/lambda_table.tsv with 24 rows.

# --- Environment FIRST (before data.table loads, per project rule) ---
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# --- Paths ---
PIPE       <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
PS_DIR     <- file.path(PIPE, "results", "emmax_ps")
QQ_DIR     <- file.path(PIPE, "results", "qq")
TABLES_DIR <- file.path(PIPE, "results", "tables")
SCRIPT_DIR <- file.path(PIPE, "scripts")

dir.create(QQ_DIR,     showWarnings = FALSE, recursive = TRUE)
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(SCRIPT_DIR, "helpers", "tag_theme.R"))

LAMBDA_TSV <- file.path(TABLES_DIR, "lambda_table.tsv")

# --- Discover and parse .ps files ---
ps_files <- list.files(PS_DIR,
                       pattern = "^morexV3__.+__.+__pc\\d+\\.ps$",
                       full.names = TRUE)
if (length(ps_files) != 24) {
  stop(sprintf("[06] FAIL: expected 24 .ps files, found %d", length(ps_files)))
}
parts <- as.data.table(do.call(rbind,
            strsplit(basename(ps_files), "__|\\.", perl = FALSE)))
setnames(parts, c("prefix", "trait", "pheno_type", "pcs_str", "ext"))
parts[, pcs := as.integer(sub("^pc", "", pcs_str))]
parts[, path := ps_files]

# --- Lambda via median-chi-square (same formula as original QQ script) ---
calc_lambda <- function(pvec) {
  chi <- qchisq(1 - pvec, df = 1)
  median(chi, na.rm = TRUE) / qchisq(0.5, df = 1)
}

# --- Process one cell ---
process_one <- function(i) {
  row <- parts[i, ]
  d <- fread(row$path, select = 4, col.names = "p", showProgress = FALSE)
  d <- d[!is.na(p) & p > 0 & p <= 1]
  n <- nrow(d)
  lambda <- calc_lambda(d$p)

  # QQ data: expected from uniform quantiles
  setorder(d, p)
  expected_all <- -log10(ppoints(n))
  observed_all <- -log10(d$p)

  # Subsample for plot (preserve all small-p tail, subsample bulk)
  TAIL_KEEP <- 50000
  BULK_KEEP <- 150000
  if (n > TAIL_KEEP + BULK_KEEP) {
    idx_tail <- seq_len(TAIL_KEEP)
    idx_bulk <- sort(sample(seq(TAIL_KEEP + 1, n), BULK_KEEP))
    idx <- c(idx_tail, idx_bulk)
  } else {
    idx <- seq_len(n)
  }
  expected <- expected_all[idx]
  observed <- observed_all[idx]
  xmax <- max(expected, na.rm = TRUE)
  ymax <- max(observed, na.rm = TRUE)

  title <- sprintf("%s | %s | %d PCs   (lambda=%.4f, n=%s)",
                   row$trait, row$pheno_type, row$pcs,
                   lambda, format(n, big.mark = ","))

  draw_qq <- function() {
    par(mar = c(4, 4.2, 2.2, 0.8), las = 1,
        cex.axis = 0.85, cex.lab = 0.95, cex.main = 0.95)
    plot(expected, observed,
         xlim = c(0, xmax),
         ylim = c(0, max(ymax, xmax)),
         pch = 16, cex = 0.35, col = "grey30",
         xlab = expression(Expected ~~ -log[10](p)),
         ylab = expression(Observed ~~ -log[10](p)),
         main = title)
    abline(0, 1, col = "red", lwd = TAG_LINE_LWD)
    legend("topleft",
           legend = c(sprintf("lambda[GC] = %.4f", lambda),
                      sprintf("n = %s SNPs", format(n, big.mark = ","))),
           bty = "n", text.col = c("red", "grey20"),
           cex = 0.85, inset = c(0.02, 0.02))
  }

  base <- sprintf("qq__%s__%s__pc%d", row$trait, row$pheno_type, row$pcs)
  tag_plot_dual(
    file.path(QQ_DIR, paste0(base, ".pdf")),
    file.path(QQ_DIR, paste0(base, ".png")),
    draw_qq,
    width_in = 4, height_in = 4
  )

  list(trait = row$trait, pheno_type = row$pheno_type,
       n_PCs = row$pcs, lambda_GC = lambda, n_valid_pvals = n)
}

# --- Run in parallel (mclapply); 8 cores * 4 BLAS threads is fine here ---
N_CORES <- min(8, nrow(parts))
cat(sprintf("[06] Processing %d cells with %d cores in parallel...\n",
            nrow(parts), N_CORES))
res <- mclapply(seq_len(nrow(parts)), process_one, mc.cores = N_CORES)

# Surface any errors from mclapply
err_idx <- which(vapply(res, inherits, logical(1), what = "try-error"))
if (length(err_idx) > 0) {
  for (i in err_idx) {
    cat(sprintf("[06] FAIL at cell %d: %s\n", i, basename(parts$path[i])))
    print(res[[i]])
  }
  stop(sprintf("[06] FAIL: %d cells errored out", length(err_idx)))
}

# --- Aggregate and write lambda_table.tsv ---
lambda_dt <- rbindlist(res)
setorder(lambda_dt, pheno_type, n_PCs, trait)
fwrite(lambda_dt, LAMBDA_TSV, sep = "\t")

cat("\n[06] lambda_table.tsv:\n")
print(lambda_dt)

# --- Checkpoints ---
cat("\n---------------------------------\n")
qq_files <- list.files(QQ_DIR, pattern = "^qq__.+\\.(pdf|png)$")
cat(sprintf("[06] CHECKPOINT: %d rows in lambda_table.tsv (expected 24)\n",
            nrow(lambda_dt)))
cat(sprintf("[06] CHECKPOINT: %d QQ files in results/qq/ (expected 48)\n",
            length(qq_files)))
if (nrow(lambda_dt) != 24 || length(qq_files) != 48) {
  stop("[06] FAIL: counts off; halting.")
}

# Range check on lambda (sanity)
lo <- min(lambda_dt$lambda_GC); hi <- max(lambda_dt$lambda_GC)
cat(sprintf("[06] CHECKPOINT: lambda range = [%.4f, %.4f]\n", lo, hi))
if (lo < 0.85 || hi > 1.20) {
  cat("[06] WARNING: lambda outside typical [0.85, 1.20] range -- inspect manually.\n")
}

cat("[06] OK: 24 QQ figures (PDF + PNG) + lambda_table.tsv ready\n")
