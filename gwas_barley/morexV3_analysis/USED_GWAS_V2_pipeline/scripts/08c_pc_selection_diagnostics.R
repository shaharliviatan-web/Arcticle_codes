#!/usr/bin/env Rscript
# 08c_pc_selection_diagnostics.R
# Build the PC-count justification kit for the paper methods:
#   - scree_plot.{pdf,png}              : % variance per PC + cumulative %
#   - lambda_vs_npcs.{pdf,png}          : lambda_GC vs n_PCs per (trait, pheno_type) -- HEADLINE
#   - pca_scatter_PC1_PC2.{pdf,png}     : PC1 vs PC2 sample scatter (popstruct visualization)
#   - pca_scatter_PC1_PC3.{pdf,png}     : PC1 vs PC3
#   - pca_scatter_PC2_PC3.{pdf,png}     : PC2 vs PC3
#   - pc_variance_table.tsv             : PC, eigenval, pct_variance, cum_pct (10 rows)
#   - lambda_vs_npcs_long.tsv           : long-format reshape of lambda_table (24 rows)
#   - lambda_vs_npcs_wide.tsv           : wide pivot: (trait, pheno_type) x (pc3, pc5, pc10) -- 8 rows
#
# Outputs to results/pc_selection/. Inputs (already on disk, unchanged):
#   intermediates/morexV3_pca.eigenvec, morexV3_pca_scree_data.tsv
#   results/tables/lambda_table.tsv
#
# Literature anchors (cite in methods):
#   Price et al. 2006 (PCA correction in GWAS)
#   Patterson, Price, Reich 2006 (PCA in popgen; Tracy-Widom)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
})

PIPE      <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
INTER_DIR <- file.path(PIPE, "intermediates")
TAB_DIR   <- file.path(PIPE, "results", "tables")
OUT_DIR   <- file.path(PIPE, "results", "pc_selection")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EIGENVEC   <- file.path(INTER_DIR, "morexV3_pca.eigenvec")
SCREE_TSV  <- file.path(INTER_DIR, "morexV3_pca_scree_data.tsv")
LAMBDA_TSV <- file.path(TAB_DIR,   "lambda_table.tsv")

stopifnot(file.exists(EIGENVEC), file.exists(SCREE_TSV), file.exists(LAMBDA_TSV))

# Plot device wrappers (TAG-aware sizing; per-cell working size for inspection)
W_IN  <- 8
H_IN  <- 5
PT    <- 10
RES   <- 300
dual <- function(file_base, draw_fn, w = W_IN, h = H_IN) {
  cairo_pdf(file.path(OUT_DIR, paste0(file_base, ".pdf")),
            width = w, height = h, pointsize = PT); draw_fn(); dev.off()
  png(file.path(OUT_DIR, paste0(file_base, ".png")),
      width = w, height = h, units = "in", res = RES, pointsize = PT,
      type = "cairo-png"); draw_fn(); dev.off()
}

# ============================================================================
# 1. Scree + cumulative variance
# ============================================================================
cat("[08c] Building scree plot (PC1..PC20) ...\n")
scree <- fread(SCREE_TSV)
stopifnot(nrow(scree) == 20)
fwrite(scree, file.path(OUT_DIR, "pc_variance_table.tsv"), sep = "\t")

draw_scree <- function() {
  par(mar = c(4.2, 4.4, 2.5, 4.4), las = 1, cex.axis = 0.85, cex.lab = 1.0)
  # bars: % variance per PC
  bp <- barplot(scree$pct_variance, names.arg = paste0("PC", scree$PC),
                ylab = "% variance explained per PC",
                col = "steelblue3", border = "white",
                ylim = c(0, max(scree$pct_variance) * 1.15),
                main = "Scree: variance captured per PC (LD-pruned ~590K SNPs, top 20 of 289 PCs)")
  # mark the PC counts under consideration as covariates (k = 3, 5, 10)
  for (k in c(3, 5, 10)) abline(v = bp[k], col = "grey60", lty = 3)
  text(x = bp[c(3, 5, 10)], y = max(scree$pct_variance) * 1.10,
       labels = paste0("k=", c(3, 5, 10)), col = "grey30", cex = 0.85)
  # secondary axis: cumulative % variance (line, capped at the max we display)
  par(new = TRUE)
  plot(bp, scree$cum_pct, type = "o", pch = 16, col = "firebrick3", lwd = 1.6,
       axes = FALSE, xlab = "", ylab = "",
       ylim = c(0, max(50, ceiling(max(scree$cum_pct) / 10) * 10)))
  axis(side = 4, col.axis = "firebrick3", col = "firebrick3")
  mtext("Cumulative % variance", side = 4, line = 2.6, col = "firebrick3", las = 0)
  legend("topright",
         legend = c("% variance per PC (bars, left axis)",
                    "Cumulative % variance (line, right axis)"),
         col = c("steelblue3", "firebrick3"),
         pch = c(15, 16), lty = c(NA, 1), pt.cex = c(1.6, 1),
         bty = "n", cex = 0.8, inset = c(0.01, 0.01))
}
# wider device for 20 bars
dual("scree_plot", draw_scree, w = 12, h = 5.5)

# ============================================================================
# 2. lambda_GC vs n_PCs   (the headline PC-count justification plot)
# ============================================================================
cat("[08c] Building lambda_vs_npcs plot ...\n")
lam <- fread(LAMBDA_TSV)
stopifnot(nrow(lam) == 24)

# Long-format and wide-format tables
fwrite(lam[, .(trait, pheno_type, n_PCs, lambda_GC)],
       file.path(OUT_DIR, "lambda_vs_npcs_long.tsv"), sep = "\t")
wide <- dcast(lam, trait + pheno_type ~ n_PCs,
              value.var = "lambda_GC")
setnames(wide, c("3", "5", "10"), c("lambda_pc3", "lambda_pc5", "lambda_pc10"))
fwrite(wide, file.path(OUT_DIR, "lambda_vs_npcs_wide.tsv"), sep = "\t")

draw_lambda <- function() {
  par(mar = c(4.2, 4.6, 2.8, 1.2), las = 1, cex.axis = 0.9, cex.lab = 1.0)
  pc_levels <- c(3, 5, 10)
  yrng <- range(lam$lambda_GC) + c(-0.01, 0.01)
  yrng <- range(c(yrng, 0.94, 1.02))
  plot(NA, xlim = c(2.5, 10.5), ylim = yrng,
       xaxt = "n",
       xlab = "Number of PCs included as covariates",
       ylab = expression(lambda[GC]),
       main = "Genomic inflation vs number of PCs (per trait x phenotype)")
  axis(1, at = pc_levels, labels = pc_levels)
  abline(h = 1.00, col = "grey30", lty = 1, lwd = 1.0)
  abline(h = 0.95, col = "firebrick3", lty = 2, lwd = 0.9)
  text(x = 10.4, y = 0.95, labels = "deflation floor",
       col = "firebrick3", pos = 3, cex = 0.7, offset = 0.2)

  # 4 traits x 2 pheno types = 8 lines
  traits <- c("betaglucan", "fiber", "protein", "starch")
  cols <- setNames(c("forestgreen", "dodgerblue3", "darkorange2", "purple3"), traits)
  for (tr in traits) for (pt in c("BLUP", "BLUE")) {
    sub <- lam[trait == tr & pheno_type == pt][order(n_PCs)]
    lines(sub$n_PCs, sub$lambda_GC, col = cols[tr],
          lty = ifelse(pt == "BLUP", 1, 2), lwd = 1.6,
          type = "o", pch = ifelse(pt == "BLUP", 16, 1))
  }
  legend("topleft", legend = traits, col = cols[traits],
         pch = 16, lwd = 1.6, lty = 1, bty = "n", cex = 0.8,
         title = "Trait")
  legend("topright", legend = c("BLUP", "BLUE"), col = "grey25",
         pch = c(16, 1), lty = c(1, 2), lwd = 1.6, bty = "n", cex = 0.8,
         title = "Pheno correction")
}
dual("lambda_vs_npcs", draw_lambda, w = 8, h = 5)

# ============================================================================
# 3. PC scatter plots (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3) -- visualize structure
# ============================================================================
cat("[08c] Building PC scatter plots ...\n")
eig <- fread(EIGENVEC, header = FALSE)
# eigenvec format: FID IID PC1..PC10 (no header)
setnames(eig, c("FID", "IID", paste0("PC", 1:(ncol(eig) - 2L))))
stopifnot(nrow(eig) == 290)

draw_scatter <- function(xpc, ypc) {
  function() {
    xv <- scree$pct_variance[xpc]; yv <- scree$pct_variance[ypc]
    par(mar = c(4.2, 4.4, 2.5, 1.0), las = 1, cex.axis = 0.9, cex.lab = 1.0)
    plot(eig[[paste0("PC", xpc)]], eig[[paste0("PC", ypc)]],
         pch = 16, cex = 0.85, col = "steelblue3",
         xlab = sprintf("PC%d (%.1f%% variance)", xpc, xv),
         ylab = sprintf("PC%d (%.1f%% variance)", ypc, yv),
         main = sprintf("PCA: PC%d vs PC%d (n=290 wild barley)", xpc, ypc))
    abline(h = 0, v = 0, col = "grey80", lty = 3)
  }
}
dual("pca_scatter_PC1_PC2", draw_scatter(1, 2), w = 6, h = 6)
dual("pca_scatter_PC1_PC3", draw_scatter(1, 3), w = 6, h = 6)
dual("pca_scatter_PC2_PC3", draw_scatter(2, 3), w = 6, h = 6)

# ============================================================================
# Checkpoints
# ============================================================================
cat("\n---------------------------------\n")
expected_pdf <- c("scree_plot.pdf", "lambda_vs_npcs.pdf",
                  "pca_scatter_PC1_PC2.pdf", "pca_scatter_PC1_PC3.pdf",
                  "pca_scatter_PC2_PC3.pdf")
expected_png <- sub("\\.pdf$", ".png", expected_pdf)
expected_tsv <- c("pc_variance_table.tsv",
                  "lambda_vs_npcs_long.tsv", "lambda_vs_npcs_wide.tsv")
all_ok <- TRUE
for (f in c(expected_pdf, expected_png, expected_tsv)) {
  p <- file.path(OUT_DIR, f)
  if (file.exists(p) && file.info(p)$size > 0) {
    cat(sprintf("[08c] OK  %-32s (%s)\n", f,
                ifelse(file.info(p)$size > 1024,
                       sprintf("%.1f KB", file.info(p)$size / 1024),
                       sprintf("%d B", file.info(p)$size))))
  } else {
    cat(sprintf("[08c] FAIL: %s missing/empty\n", f)); all_ok <- FALSE
  }
}
if (!all_ok) stop("[08c] FAIL")
cat("\n[08c] OK: PC-selection diagnostics ready in results/pc_selection/\n")
cat("[08c] Cite in methods: Price et al. 2006 (PCA + GWAS); Patterson, Price, Reich 2006 (PCA + Tracy-Widom).\n")
