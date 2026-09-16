#!/usr/bin/env Rscript
# 06c_qq_lambda_v3_grid.R
# v3-ONLY versions of the QQ / lambda summary figures (no v2 overlay), so the new
# run can be read on its own terms. The v2-overlaid counterparts are produced by
# 06b_qq_compare_v2_v3.R.
#
# Writes into results/comparison_v2_vs_v3/02_qq_lambda/v3_only/:
#   qq_grid_BLUP.{pdf,png}    12-panel QQ grid, v3 only
#   qq_grid_BLUE.{pdf,png}
#   lambda_v3.{pdf,png}       lambda per cell vs the lambda = 1 line
#   lambda_v3.tsv
#
# Created 2026-08-16 for the v3 re-run review gate.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(parallel) })

PIPE    <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
PS_DIR  <- file.path(PIPE, "results", "emmax_ps")
OUT_DIR <- file.path(PIPE, "results", "comparison_v2_vs_v3", "02_qq_lambda", "v3_only")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(PIPE, "scripts", "helpers", "tag_theme.R"))

N_PRUNED <- 111017L
ALPHA    <- 0.10
BONF     <- -log10(ALPHA / N_PRUNED)          # 6.0454
BONF_ALL <- -log10(ALPHA / 7110996)           # 7.8519

COL_V3 <- "#1F6FB4"; COL_ID <- "#D62728"; COL_ALL <- "#7F7F7F"
TRAITS <- c("betaglucan", "fiber", "protein", "starch")
PCS    <- c(3L, 5L, 10L)

calc_lambda <- function(p) median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)

load_cell <- function(trait, pheno, pc) {
  f <- file.path(PS_DIR, sprintf("morexV3__%s__%s__pc%d.ps", trait, pheno, pc))
  d <- fread(f, select = 4, col.names = "p", showProgress = FALSE)
  p <- sort(d$p[!is.na(d$p) & d$p > 0 & d$p <= 1])
  n <- length(p)
  idx <- if (n > 150000L) c(seq_len(50000L), sort(sample(seq(50001L, n), 100000L))) else seq_len(n)
  list(expected = -log10(ppoints(n))[idx], observed = -log10(p)[idx],
       lambda = calc_lambda(p), n = n, trait = trait, pheno = pheno, pc = pc)
}

grid <- CJ(trait = TRAITS, pheno = c("BLUP", "BLUE"), pc = PCS, sorted = FALSE)
cat(sprintf("[06c] Loading %d cells ...\n", nrow(grid)))
cells <- mclapply(seq_len(nrow(grid)),
                  function(i) load_cell(grid$trait[i], grid$pheno[i], grid$pc[i]),
                  mc.cores = 8)
names(cells) <- sprintf("%s|%s|%d", grid$trait, grid$pheno, grid$pc)

draw_panel <- function(cl) {
  xmax <- max(cl$expected); ymax <- max(cl$observed)
  par(mar = c(3.4, 3.6, 2.4, 0.6), las = 1, mgp = c(2.1, 0.6, 0),
      cex.axis = 0.72, cex.lab = 0.8, cex.main = 0.82)
  plot(cl$expected, cl$observed, pch = 16, cex = 0.28, col = COL_V3,
       xlim = c(0, xmax), ylim = c(0, max(xmax, ymax)),
       xlab = expression(Expected ~~ -log[10](p)),
       ylab = expression(Observed ~~ -log[10](p)),
       main = sprintf("%s | %s | %d PCs", cl$trait, cl$pheno, cl$pc))
  abline(0, 1, col = COL_ID, lwd = TAG_LINE_LWD)
  abline(h = BONF,     col = COL_V3,  lty = 3, lwd = TAG_LINE_LWD)
  abline(h = BONF_ALL, col = COL_ALL, lty = 2, lwd = TAG_LINE_LWD)
  legend("topleft", bty = "n", cex = 0.7, text.col = COL_V3,
         legend = sprintf("lambda = %.4f", cl$lambda), inset = c(-0.02, -0.01))
}

for (ph in c("BLUP", "BLUE")) {
  cat(sprintf("[06c] Drawing %s grid ...\n", ph))
  tag_plot_dual(
    file.path(OUT_DIR, sprintf("qq_grid_%s.pdf", ph)),
    file.path(OUT_DIR, sprintf("qq_grid_%s.png", ph)),
    function() {
      layout(rbind(matrix(1:12, nrow = 4, byrow = TRUE), rep(13, 3)),
             heights = c(1, 1, 1, 1, 0.22))
      for (tr in TRAITS) for (pc in PCS) draw_panel(cells[[sprintf("%s|%s|%d", tr, ph, pc)]])
      par(mar = c(0, 0, 0, 0)); plot.new()
      legend("center", horiz = TRUE, bty = "n", cex = 0.9,
             legend = c(sprintf("Bonferroni, LD-pruned n=%s  (%.3f)",
                                format(N_PRUNED, big.mark = ","), BONF),
                        sprintf("Bonferroni, all SNPs  (%.3f)", BONF_ALL),
                        "y = x (null)"),
             col = c(COL_V3, COL_ALL, COL_ID), lty = c(3, 2, 1), lwd = c(1.2, 1.2, 1.2))
    }, width_in = 10.5, height_in = 13)
}

lam <- rbindlist(lapply(cells, function(c0) data.table(
  trait = c0$trait, pheno_type = c0$pheno, n_PCs = c0$pc, lambda_GC = c0$lambda)))
setorder(lam, pheno_type, n_PCs, trait)
fwrite(lam, file.path(OUT_DIR, "lambda_v3.tsv"), sep = "\t")

tag_plot_dual(file.path(OUT_DIR, "lambda_v3.pdf"), file.path(OUT_DIR, "lambda_v3.png"),
  function() {
    d <- copy(lam); setorder(d, -pheno_type, -n_PCs, -trait)
    d[, lab := sprintf("%s  %s  %dPC", trait, pheno_type, n_PCs)]
    n <- nrow(d); xr <- range(c(d$lambda_GC, 1)); xr <- xr + c(-1, 1) * diff(xr) * 0.15
    par(mar = c(4, 9.5, 2.6, 1), las = 1, mgp = c(2.3, 0.6, 0), cex.axis = 0.7)
    plot(NA, xlim = xr, ylim = c(0.5, n + 0.5), yaxt = "n",
         xlab = expression(lambda[GC]), ylab = "",
         main = "Genomic inflation, v3 correction")
    abline(v = 1, col = COL_ID, lwd = 1.2)
    axis(2, at = seq_len(n), labels = d$lab, tick = FALSE, cex.axis = 0.62)
    segments(1, seq_len(n), d$lambda_GC, seq_len(n), col = "grey75", lwd = 1.4)
    points(d$lambda_GC, seq_len(n), pch = 16, cex = 0.9, col = COL_V3)
    legend("bottomleft", bty = "n", cex = 0.72, pch = c(16, NA), lty = c(NA, 1),
           col = c(COL_V3, COL_ID),
           legend = c("v3 (1000kb window)", "lambda = 1 (no inflation)"))
  }, width_in = 7.5, height_in = 8)

cat(sprintf("\n[06c] CHECKPOINT: cells = %d (expected 24)\n", nrow(lam)))
cat(sprintf("[06c] lambda range = [%.4f, %.4f]\n", min(lam$lambda_GC), max(lam$lambda_GC)))
cat(sprintf("[06c] OK: outputs in %s\n", OUT_DIR))
