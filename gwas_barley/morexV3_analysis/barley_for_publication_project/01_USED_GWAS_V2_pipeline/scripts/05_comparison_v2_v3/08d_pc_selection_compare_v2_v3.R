#!/usr/bin/env Rscript
# 08d_pc_selection_compare_v2_v3.R
# PC-selection diagnostics comparing the two population-structure corrections:
#   v2 = PCA + aIBS kinship on --indep-pairwise 50 5 0.2     (590,462 SNPs)
#   v3 = PCA + aIBS kinship on --indep-pairwise 1000kb 1 0.2 (111,017 SNPs)
#
# Writes into results/comparison_v2_vs_v3/01_pc_selection/v2_vs_v3/:
#   scree_compare.{pdf,png}          % variance per PC, both versions
#   scree_cumulative_compare.{pdf,png}
#   scree_elbow_compare.{pdf,png}    line/elbow style, PC1..PC20
#   pca_scatter_compare.{pdf,png}    PC1-PC2 / PC1-PC3 / PC2-PC3, v2 vs v3 side by side
#   pc_score_correlation.{pdf,png}   |r| between v2 and v3 sample scores, PC1..PC10
#   lambda_vs_npcs_compare.{pdf,png} lambda vs n_PCs, both versions
#   pc_variance_compare.tsv
#   pc_score_correlation.tsv
#
# Inputs (all already on disk):
#   v3: intermediates/morexV3_pca_scree_data.tsv, morexV3_pca.eigenvec
#       results/tables/lambda_table.tsv
#   v2: intermediates/_archive_v2_win50snp_2026-08-16/... (same names)
#       results/_archive/_archive_v2_win50snp_2026-08-16/tables/lambda_table.tsv
#
# Created 2026-08-16 for the v3 re-run review gate.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))

PIPE    <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
INTER   <- file.path(PIPE, "intermediates")
ARC_I   <- file.path(INTER, "_archive_v2_win50snp_2026-08-16")
ARC_R   <- file.path(PIPE, "results", "_archive", "_archive_v2_win50snp_2026-08-16")
OUT_DIR <- file.path(PIPE, "results", "comparison_v2_vs_v3", "01_pc_selection", "v2_vs_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_V2 <- "#9E9E9E"; COL_V3 <- "#1F6FB4"; COL_ACC <- "#D62728"
LAB_V2 <- "v2: 50-SNP window (590,462 SNPs)"
LAB_V3 <- "v3: 1000kb window (111,017 SNPs)"

dual <- function(base, draw_fn, w = 8, h = 5) {
  cairo_pdf(file.path(OUT_DIR, paste0(base, ".pdf")), width = w, height = h,
            family = "sans", pointsize = 10); draw_fn(); dev.off()
  png(file.path(OUT_DIR, paste0(base, ".png")), width = w, height = h,
      units = "in", res = 300, family = "sans", pointsize = 10); draw_fn(); dev.off()
  cat(sprintf("[08d] OK  %s.{pdf,png}\n", base))
}

# ---------------- Load ----------------
s3 <- fread(file.path(INTER, "morexV3_pca_scree_data.tsv"))
s2 <- fread(file.path(ARC_I, "morexV3_pca_scree_data.tsv"))
e3 <- fread(file.path(INTER, "morexV3_pca.eigenvec"))
e2 <- fread(file.path(ARC_I, "morexV3_pca.eigenvec"))
l3 <- fread(file.path(PIPE,  "results", "tables", "lambda_table.tsv"))
l2 <- fread(file.path(ARC_R, "tables", "lambda_table.tsv"))

NPC <- min(nrow(s2), nrow(s3), 20L)
s2 <- s2[1:NPC]; s3 <- s3[1:NPC]

# ---------------- Variance table ----------------
vt <- data.table(PC = s3$PC,
                 eigenval_v2 = s2$eigenval, eigenval_v3 = s3$eigenval,
                 pct_v2 = s2$pct_variance,  pct_v3 = s3$pct_variance,
                 cum_v2 = s2$cum_pct,       cum_v3 = s3$cum_pct)
vt[, pct_delta := pct_v3 - pct_v2]
fwrite(vt, file.path(OUT_DIR, "pc_variance_compare.tsv"), sep = "\t")

# ---------------- Scree: % variance per PC (grouped bars) ----------------
dual("scree_compare", function() {
  par(mar = c(4.2, 4.4, 3, 1), las = 1, mgp = c(2.6, 0.6, 0))
  m <- rbind(vt$pct_v2, vt$pct_v3)
  bp <- barplot(m, beside = TRUE, names.arg = vt$PC,
                col = c(COL_V2, COL_V3), border = NA, space = c(0, 0.6),
                xlab = "Principal component", ylab = "% variance explained",
                main = "Scree: variance per PC, v2 vs v3 correction",
                ylim = c(0, max(m) * 1.18))
  legend("topright", bty = "n", fill = c(COL_V2, COL_V3), border = NA,
         legend = c(LAB_V2, LAB_V3), cex = 0.85)
})

# ---------------- Scree: cumulative ----------------
dual("scree_cumulative_compare", function() {
  par(mar = c(4.2, 4.4, 3, 1), las = 1, mgp = c(2.6, 0.6, 0))
  plot(vt$PC, vt$cum_v2, type = "o", pch = 16, cex = 0.8, lwd = 1.8, col = COL_V2,
       ylim = c(0, max(vt$cum_v2, vt$cum_v3) * 1.1),
       xlab = "Principal component", ylab = "Cumulative % variance explained",
       main = "Cumulative variance, v2 vs v3 correction")
  lines(vt$PC, vt$cum_v3, type = "o", pch = 16, cex = 0.8, lwd = 1.8, col = COL_V3)
  for (k in c(3, 5, 10)) abline(v = k, col = COL_ACC, lty = 3, lwd = 0.8)
  text(c(3, 5, 10), rep(max(vt$cum_v2) * 1.06, 3), labels = paste0(c(3, 5, 10), " PCs"),
       col = COL_ACC, cex = 0.7)
  legend("bottomright", bty = "n", lwd = 1.8, pch = 16, col = c(COL_V2, COL_V3),
         legend = c(LAB_V2, LAB_V3), cex = 0.85)
})

# ---------------- Scree: elbow (line only, 10b style) ----------------
dual("scree_elbow_compare", function() {
  par(mar = c(4.2, 4.4, 3, 1), las = 1, mgp = c(2.6, 0.6, 0))
  plot(vt$PC, vt$pct_v2, type = "o", pch = 16, cex = 0.9, lwd = 1.8, col = COL_V2,
       ylim = c(0, max(vt$pct_v2, vt$pct_v3) * 1.12),
       xlab = "Principal component", ylab = "% variance explained",
       main = "Scree elbow, v2 vs v3 correction")
  lines(vt$PC, vt$pct_v3, type = "o", pch = 16, cex = 0.9, lwd = 1.8, col = COL_V3)
  legend("topright", bty = "n", lwd = 1.8, pch = 16, col = c(COL_V2, COL_V3),
         legend = c(LAB_V2, LAB_V3), cex = 0.85)
})

# ---------------- PCA scatters: v2 vs v3 side by side ----------------
# eigenvec cols: FID IID PC1 PC2 ... -> PCk is column k+2
pcv <- function(dt, k) as.numeric(dt[[k + 2L]])
dual("pca_scatter_compare", function() {
  pairs_list <- list(c(1, 2), c(1, 3), c(2, 3))
  par(mfrow = c(2, 3), mar = c(4, 4.2, 2.6, 0.8), las = 1, mgp = c(2.4, 0.6, 0),
      cex.axis = 0.8, cex.lab = 0.9, cex.main = 0.95)
  for (ver in c("v2", "v3")) {
    dt  <- if (ver == "v2") e2 else e3
    sc  <- if (ver == "v2") s2 else s3
    col <- if (ver == "v2") COL_V2 else COL_V3
    for (pr in pairs_list) {
      plot(pcv(dt, pr[1]), pcv(dt, pr[2]), pch = 16, cex = 0.6, col = col,
           xlab = sprintf("PC%d (%.2f%%)", pr[1], sc$pct_variance[pr[1]]),
           ylab = sprintf("PC%d (%.2f%%)", pr[2], sc$pct_variance[pr[2]]),
           main = sprintf("%s  PC%d vs PC%d", ver, pr[1], pr[2]))
      abline(h = 0, v = 0, col = "grey85", lwd = 0.6)
    }
  }
}, w = 11, h = 7.5)

# ---------------- PC score correlation v2 vs v3 ----------------
K <- min(10L, ncol(e3) - 2L, ncol(e2) - 2L)
cors <- data.table(PC = 1:K,
                   r = vapply(1:K, function(k) cor(pcv(e2, k), pcv(e3, k)), numeric(1)))
cors[, abs_r := abs(r)]
fwrite(cors, file.path(OUT_DIR, "pc_score_correlation.tsv"), sep = "\t")

dual("pc_score_correlation", function() {
  par(mar = c(4.2, 4.4, 3, 1), las = 1, mgp = c(2.6, 0.6, 0))
  bp <- barplot(cors$abs_r, names.arg = cors$PC, col = COL_V3, border = NA,
                ylim = c(0, 1.08), xlab = "Principal component",
                ylab = expression("|r| between v2 and v3 sample scores"),
                main = "Do v2 and v3 recover the same population structure?")
  abline(h = c(0.9, 1), col = c(COL_ACC, "grey60"), lty = c(2, 1), lwd = c(1, 0.8))
  text(bp, cors$abs_r + 0.035, sprintf("%.3f", cors$abs_r), cex = 0.7)
  legend("bottomleft", bty = "n", lty = 2, col = COL_ACC, cex = 0.78,
         legend = "|r| = 0.9")
})

# ---------------- Cross-correlation heatmap: every v2 PC vs every v3 PC ----------------
# The element-wise (diagonal) correlation above understates the agreement, because
# PCs with near-equal eigenvalues can swap rank between the two marker sets.
# This panel shows the full |r| matrix so a swap is visible as an off-diagonal hit.
M2 <- as.matrix(e2[, 3:(2 + K)]); M3 <- as.matrix(e3[, 3:(2 + K)])
CC <- abs(cor(M2, M3))
best <- data.table(v3_PC = 1:K,
                   best_v2_PC = apply(CC, 2, which.max),
                   best_abs_r = apply(CC, 2, max),
                   diag_abs_r = diag(CC))
fwrite(best, file.path(OUT_DIR, "pc_best_match.tsv"), sep = "\t")

dual("pc_crosscorrelation", function() {
  par(mar = c(4.2, 4.4, 3.2, 3.5), las = 1, mgp = c(2.6, 0.6, 0))
  pal <- colorRampPalette(c("white", "#C6DBEF", COL_V3, "#08306B"))(100)
  image(1:K, 1:K, CC, col = pal, zlim = c(0, 1), axes = FALSE,
        xlab = "v2 PC (50-SNP window)", ylab = "v3 PC (1000kb window)",
        main = "PC score agreement |r|: swaps appear off-diagonal")
  axis(1, at = 1:K, cex.axis = 0.8); axis(2, at = 1:K, cex.axis = 0.8)
  box()
  for (i in 1:K) for (j in 1:K) {
    if (CC[i, j] >= 0.45)
      text(i, j, sprintf("%.2f", CC[i, j]), cex = 0.62,
           col = if (CC[i, j] > 0.75) "white" else "grey15")
  }
  # outline the rank-matched cell of each v3 PC
  for (j in 1:K) {
    i <- best$best_v2_PC[j]
    rect(i - .5, j - .5, i + .5, j + .5, border = COL_ACC, lwd = 1.6)
  }
  # legend sits inside the (empty) top-left corner of the matrix, clear of the title
  legend("topleft", bty = "n", cex = 0.72, border = COL_ACC, fill = NA,
         legend = "best match for that v3 PC", inset = c(0.01, 0.01))
}, w = 7.5, h = 6.5)

# ---------------- lambda vs n_PCs, both versions ----------------
l3[, ver := "v3"]; l2[, ver := "v2"]
lam <- rbind(l3, l2, fill = TRUE)
dual("lambda_vs_npcs_compare", function() {
  traits <- sort(unique(lam$trait))
  par(mfrow = c(2, 4), mar = c(3.8, 4, 2.6, 0.7), las = 1, mgp = c(2.3, 0.6, 0),
      cex.axis = 0.78, cex.lab = 0.85, cex.main = 0.9)
  for (ph in c("BLUP", "BLUE")) for (tr in traits) {
    d <- lam[trait == tr & pheno_type == ph]
    yr <- range(c(d$lambda_GC, 1)); yr <- yr + c(-1, 1) * diff(yr) * 0.25
    plot(NA, xlim = c(2, 11), ylim = yr, xaxt = "n",
         xlab = "number of PCs", ylab = expression(lambda[GC]),
         main = sprintf("%s | %s", tr, ph))
    axis(1, at = c(3, 5, 10))
    abline(h = 1, col = COL_ACC, lwd = 1)
    for (v in c("v2", "v3")) {
      dv <- d[ver == v]; setorder(dv, n_PCs)
      lines(dv$n_PCs, dv$lambda_GC, type = "o", pch = 16, cex = 0.9, lwd = 1.7,
            col = if (v == "v2") COL_V2 else COL_V3)
    }
  }
}, w = 13, h = 7)

# ---------------- Checkpoints ----------------
cat("\n---------------------------------\n")
cat(sprintf("[08d] CHECKPOINT: PCs compared = %d\n", NPC))
cat(sprintf("[08d] v2 cum%% at PC3/PC5/PC10 = %.2f / %.2f / %.2f\n",
            vt$cum_v2[3], vt$cum_v2[5], vt$cum_v2[10]))
cat(sprintf("[08d] v3 cum%% at PC3/PC5/PC10 = %.2f / %.2f / %.2f\n",
            vt$cum_v3[3], vt$cum_v3[5], vt$cum_v3[10]))
cat(sprintf("[08d] PC score |r| v2~v3, PC1..PC%d: min=%.3f  median=%.3f\n",
            K, min(cors$abs_r), median(cors$abs_r)))
cat(sprintf("[08d] OK: outputs in %s\n", OUT_DIR))
