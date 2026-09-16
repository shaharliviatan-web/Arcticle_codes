#!/usr/bin/env Rscript
# 06b_qq_compare_v2_v3.R
# Side-by-side QQ comparison of the two population-structure corrections:
#   v2 = PCA + aIBS kinship built on --indep-pairwise 50 5 0.2   (590,462 SNPs)
#   v3 = PCA + aIBS kinship built on --indep-pairwise 1000kb 1 0.2 (111,017 SNPs)
#
# Both corrections were applied to the SAME 7.1M-SNP tped and the same 24-cell
# grid (4 traits x {BLUP,BLUE} x {3,5,10} PCs), so any difference in the QQ
# curves is attributable to the correction alone.
#
# Outputs (results/comparison_v2_vs_v3/):
#   qq_overlay_BLUP.{pdf,png}   12-panel grid, v2 vs v3 overlaid
#   qq_overlay_BLUE.{pdf,png}   12-panel grid, v2 vs v3 overlaid
#   lambda_shift.{pdf,png}      dumbbell plot: lambda movement toward 1.0
#   lambda_compare.tsv          per-cell lambda, both versions, delta
#
# Created 2026-08-16 for the v3 re-run review gate.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

PIPE    <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
PS_V3   <- file.path(PIPE, "results", "emmax_ps")
ARCHIVE <- file.path(PIPE, "results", "_archive", "_archive_v2_win50snp_2026-08-16")
PS_V2   <- file.path(ARCHIVE, "emmax_ps")
OUT_DIR <- file.path(PIPE, "results", "comparison_v2_vs_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(PIPE, "scripts", "helpers", "tag_theme.R"))

N_PRUNED_V2 <- 590462L
N_PRUNED_V3 <- 111017L
ALPHA       <- 0.10
BONF_V2     <- -log10(ALPHA / N_PRUNED_V2)   # 6.7712
BONF_V3     <- -log10(ALPHA / N_PRUNED_V3)   # 6.0454

COL_V2 <- "#9E9E9E"   # grey  - old correction
COL_V3 <- "#1F6FB4"   # blue  - new correction
COL_ID <- "#D62728"   # red   - y = x

TRAITS <- c("betaglucan", "fiber", "protein", "starch")
PCS    <- c(3L, 5L, 10L)

calc_lambda <- function(p) median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

# Read one .ps, return sorted -log10(p) thinned for plotting + lambda on the full set.
load_cell <- function(dir, trait, pheno, pc) {
  f <- file.path(dir, sprintf("morexV3__%s__%s__pc%d.ps", trait, pheno, pc))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, select = 4, col.names = "p", showProgress = FALSE)
  p <- d$p[!is.na(d$p) & d$p > 0 & d$p <= 1]
  n <- length(p)
  lambda <- calc_lambda(p)
  p <- sort(p)
  exp_all <- -log10(ppoints(n))
  obs_all <- -log10(p)
  # keep the whole tail, subsample the bulk (same policy as 06_qq_lambda.R)
  TAIL_KEEP <- 50000L; BULK_KEEP <- 100000L
  idx <- if (n > TAIL_KEEP + BULK_KEEP) {
    c(seq_len(TAIL_KEEP), sort(sample(seq(TAIL_KEEP + 1L, n), BULK_KEEP)))
  } else seq_len(n)
  list(expected = exp_all[idx], observed = obs_all[idx], lambda = lambda, n = n)
}

# ---- Load all 48 cells (24 per version) in parallel ----
grid <- CJ(trait = TRAITS, pheno = c("BLUP", "BLUE"), pc = PCS, sorted = FALSE)
cat(sprintf("[06b] Loading %d cells x 2 versions ...\n", nrow(grid)))

load_pair <- function(i) {
  g <- grid[i, ]
  list(key = sprintf("%s|%s|%d", g$trait, g$pheno, g$pc),
       trait = g$trait, pheno = g$pheno, pc = g$pc,
       v3 = load_cell(PS_V3, g$trait, g$pheno, g$pc),
       v2 = load_cell(PS_V2, g$trait, g$pheno, g$pc))
}
cells <- mclapply(seq_len(nrow(grid)), load_pair, mc.cores = 8)

err <- which(vapply(cells, inherits, logical(1), what = "try-error"))
if (length(err)) stop(sprintf("[06b] FAIL: %d cells errored", length(err)))

names(cells) <- vapply(cells, `[[`, character(1), "key")

# ---- Panel drawing ----
draw_panel <- function(cell) {
  v2 <- cell$v2; v3 <- cell$v3
  if (is.null(v3)) { plot.new(); title(main = "missing"); return(invisible()) }

  xmax <- max(v3$expected, if (!is.null(v2)) v2$expected else 0, na.rm = TRUE)
  ymax <- max(v3$observed, if (!is.null(v2)) v2$observed else 0, na.rm = TRUE)
  lim  <- max(xmax, ymax)

  par(mar = c(3.4, 3.6, 2.4, 0.6), las = 1, mgp = c(2.1, 0.6, 0),
      cex.axis = 0.72, cex.lab = 0.8, cex.main = 0.82)
  plot(NA, xlim = c(0, xmax), ylim = c(0, lim),
       xlab = expression(Expected ~~ -log[10](p)),
       ylab = expression(Observed ~~ -log[10](p)),
       main = sprintf("%s | %s | %d PCs", cell$trait, cell$pheno, cell$pc))

  # Bonferroni reference lines (dotted, version-coloured)
  abline(h = BONF_V2, col = COL_V2, lty = 3, lwd = TAG_LINE_LWD)
  abline(h = BONF_V3, col = COL_V3, lty = 3, lwd = TAG_LINE_LWD)

  if (!is.null(v2)) points(v2$expected, v2$observed, pch = 16, cex = 0.28, col = COL_V2)
  points(v3$expected, v3$observed, pch = 16, cex = 0.28, col = COL_V3)
  abline(0, 1, col = COL_ID, lwd = TAG_LINE_LWD)

  legend("topleft", bty = "n", cex = 0.68, inset = c(-0.02, -0.01),
         legend = c(sprintf("v3  lambda = %.4f", v3$lambda),
                    if (!is.null(v2)) sprintf("v2  lambda = %.4f", v2$lambda) else NULL),
         text.col = c(COL_V3, COL_V2))
}

draw_grid <- function(pheno) {
  # 4 trait rows x 3 PC columns, plus a figure-level legend row
  layout(rbind(matrix(1:12, nrow = 4, byrow = TRUE), rep(13, 3)),
         heights = c(1, 1, 1, 1, 0.22))
  for (tr in TRAITS) for (pc in PCS) {
    draw_panel(cells[[sprintf("%s|%s|%d", tr, pheno, pc)]])
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", horiz = TRUE, bty = "n", cex = 0.9,
         legend = c(sprintf("v3: 1000kb window, %s SNPs (Bonf %.3f)",
                            format(N_PRUNED_V3, big.mark = ","), BONF_V3),
                    sprintf("v2: 50-SNP window, %s SNPs (Bonf %.3f)",
                            format(N_PRUNED_V2, big.mark = ","), BONF_V2),
                    "y = x (null)"),
         col = c(COL_V3, COL_V2, COL_ID), pch = c(16, 16, NA),
         lty = c(NA, NA, 1), lwd = c(NA, NA, 1.2))
}

for (ph in c("BLUP", "BLUE")) {
  cat(sprintf("[06b] Drawing %s overlay grid ...\n", ph))
  tag_plot_dual(file.path(OUT_DIR, sprintf("qq_overlay_%s.pdf", ph)),
                file.path(OUT_DIR, sprintf("qq_overlay_%s.png", ph)),
                function() draw_grid(ph),
                width_in = 10.5, height_in = 13)
}

# ---- Lambda comparison table ----
lam <- rbindlist(lapply(cells, function(c0) data.table(
  trait = c0$trait, pheno_type = c0$pheno, n_PCs = c0$pc,
  lambda_v2 = if (is.null(c0$v2)) NA_real_ else c0$v2$lambda,
  lambda_v3 = c0$v3$lambda)))
lam[, delta := lambda_v3 - lambda_v2]
lam[, dist_to_1_v2 := abs(lambda_v2 - 1)]
lam[, dist_to_1_v3 := abs(lambda_v3 - 1)]
lam[, improved := dist_to_1_v3 < dist_to_1_v2]
setorder(lam, pheno_type, n_PCs, trait)
fwrite(lam, file.path(OUT_DIR, "lambda_compare.tsv"), sep = "\t")

# ---- Lambda dumbbell plot ----
draw_dumbbell <- function() {
  d <- copy(lam)
  setorder(d, -pheno_type, -n_PCs, -trait)
  d[, lab := sprintf("%s  %s  %dPC", trait, pheno_type, n_PCs)]
  n <- nrow(d)
  xr <- range(c(d$lambda_v2, d$lambda_v3, 1), na.rm = TRUE)
  xr <- xr + c(-1, 1) * diff(xr) * 0.12

  par(mar = c(4, 9.5, 2.6, 1), las = 1, mgp = c(2.3, 0.6, 0),
      cex.axis = 0.7, cex.lab = 0.85)
  plot(NA, xlim = xr, ylim = c(0.5, n + 0.5), yaxt = "n",
       xlab = expression(lambda[GC]), ylab = "",
       main = "Genomic inflation: v2 -> v3 correction")
  abline(v = 1, col = COL_ID, lwd = 1.2)
  axis(2, at = seq_len(n), labels = d$lab, tick = FALSE, cex.axis = 0.62)
  segments(d$lambda_v2, seq_len(n), d$lambda_v3, seq_len(n),
           col = "grey65", lwd = 1.6)
  points(d$lambda_v2, seq_len(n), pch = 16, cex = 0.85, col = COL_V2)
  points(d$lambda_v3, seq_len(n), pch = 16, cex = 0.85, col = COL_V3)
  legend("bottomleft", bty = "n", cex = 0.72, pch = 16,
         col = c(COL_V2, COL_V3, COL_ID),
         legend = c("v2 (50-SNP window)", "v3 (1000kb window)", "lambda = 1 (no inflation)"))
}
tag_plot_dual(file.path(OUT_DIR, "lambda_shift.pdf"),
              file.path(OUT_DIR, "lambda_shift.png"),
              draw_dumbbell, width_in = 7.5, height_in = 8)

# ---- Checkpoints ----
cat("\n---------------------------------\n")
cat(sprintf("[06b] CHECKPOINT: cells compared = %d (expected 24)\n", nrow(lam)))
cat(sprintf("[06b] CHECKPOINT: cells closer to lambda=1 under v3 = %d / %d\n",
            sum(lam$improved, na.rm = TRUE), nrow(lam)))
cat(sprintf("[06b] v2 lambda range = [%.4f, %.4f];  mean |lambda-1| = %.4f\n",
            min(lam$lambda_v2), max(lam$lambda_v2), mean(lam$dist_to_1_v2)))
cat(sprintf("[06b] v3 lambda range = [%.4f, %.4f];  mean |lambda-1| = %.4f\n",
            min(lam$lambda_v3), max(lam$lambda_v3), mean(lam$dist_to_1_v3)))
cat(sprintf("[06b] OK: figures + lambda_compare.tsv in %s\n", OUT_DIR))
