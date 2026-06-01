#!/usr/bin/env Rscript
# 10b_extra_scree_plots.R
# Adds classic "elbow" scree plots (% variance explained on y-axis, PC on x-axis)
# in the line-and-points style of the reference photo, to the publication folder.
#
# Writes into the chosen-config publication dir
# (results/publication_BonfOnly_BLUP_3PC/pc_selection/):
#   Figure_S1b_PC_scree_withoutPC1.{pdf,png}        -- PC2..PC20 bar + cumulative line
#   Figure_S1c_PC_scree_line.{pdf,png}              -- PC1..PC20 elbow
#   Figure_S1d_PC_scree_line_withoutPC1.{pdf,png}   -- PC2..PC20 elbow (PC1 dropped)
#
# Style:
#   - elbow (S1c/S1d): standard white background, black axes, single teal line + filled circles.
#   - bar (S1b): per-PC % variance bars + red cumulative-% overlay on a right axis.
#   - no plot titles and no k-PC reference lines (clean for publication).
#
# Reads: intermediates/morexV3_pca_scree_data.tsv (PC1..PC20, % computed against
# the full 289-eigenvalue trace).

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
})

PIPE      <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
SCREE_TSV <- file.path(PIPE, "intermediates", "morexV3_pca_scree_data.tsv")
DEST_DIRS <- c(
  file.path(PIPE, "results", "publication_BonfOnly_BLUP_3PC", "pc_selection")
)
stopifnot(file.exists(SCREE_TSV))
for (d in DEST_DIRS) dir.create(d, showWarnings = FALSE, recursive = TRUE)

scree <- fread(SCREE_TSV)
stopifnot(all(c("PC", "eigenval", "pct_variance", "cum_pct") %in% names(scree)))

# ---- Style constants (standard publication: white background, black axes) ----
LINE_COL    <- "#1F7A8C"   # teal accent for the line/points
POINT_COL   <- "#1F7A8C"
W_IN <- 8; H_IN <- 5
PNG_RES <- 600

draw_elbow <- function(include_pc1 = TRUE) {
  d <- if (include_pc1) scree[PC <= 20] else scree[PC >= 2 & PC <= 20]

  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar = c(4.4, 5.0, 1.5, 1.5), las = 1,
      cex.axis = 1.0, cex.lab = 1.15,
      bg = "white", col.axis = "black", col.lab = "black",
      family = "sans")

  ymax <- max(d$pct_variance) * 1.10
  plot(d$PC, d$pct_variance, type = "n",
       xlim = c(min(d$PC) - 0.4, max(d$PC) + 0.4),
       ylim = c(0, ymax),
       xlab = "Principal component", ylab = "% variance explained",
       main = "", xaxt = "n", bty = "l")
  axis(1, at = d$PC, labels = d$PC, cex.axis = 0.9)
  lines(d$PC, d$pct_variance, col = LINE_COL, lwd = 2.6)
  points(d$PC, d$pct_variance, pch = 21, bg = POINT_COL, col = LINE_COL,
         cex = 1.5, lwd = 1.2)
}

# Bar style (PC2..PC20): % variance per PC as bars + cumulative-% line on right axis.
draw_bar_noPC1 <- function() {
  d <- scree[PC >= 2 & PC <= 20]
  op <- par(no.readonly = TRUE); on.exit(par(op))
  par(mar = c(4.2, 4.6, 2.8, 4.6), las = 1, cex.axis = 0.95, cex.lab = 1.05, bg = "white")
  bp <- barplot(d$pct_variance, names.arg = d$PC,
                col = "#4C9AC4", border = "#2E6A8E",
                ylim = c(0, max(d$pct_variance) * 1.15),
                xlab = "Principal component (PC1 excluded)",
                ylab = "% variance explained",
                main = "")
  par(new = TRUE)
  plot(bp, d$cum_pct, type = "o", pch = 19, col = "#B22222",
       axes = FALSE, xlab = "", ylab = "", ylim = c(0, 100), xlim = range(bp))
  axis(4, col.axis = "#B22222", col = "#B22222")
  mtext("Cumulative % variance", side = 4, line = 3, col = "#B22222", cex = 1.05, las = 0)
  legend("topright", bty = "n", cex = 0.85,
         legend = c("% per PC", "cumulative %"),
         col = c("#4C9AC4", "#B22222"),
         pch = c(15, 19), lty = c(NA, 1), lwd = c(NA, 1.5))
}

emit <- function(file_base, drawer) {
  for (dest in DEST_DIRS) {
    pdf_path <- file.path(dest, paste0(file_base, ".pdf"))
    png_path <- file.path(dest, paste0(file_base, ".png"))
    cairo_pdf(pdf_path, width = W_IN, height = H_IN); drawer(); dev.off()
    png(png_path, width = W_IN, height = H_IN, units = "in",
        res = PNG_RES, type = "cairo-png"); drawer(); dev.off()
    cat(sprintf("[10b] wrote %s.{pdf,png}\n", file.path(dest, file_base)))
  }
}

emit("Figure_S1b_PC_scree_withoutPC1",      function() draw_bar_noPC1())
emit("Figure_S1c_PC_scree_line",            function() draw_elbow(include_pc1 = TRUE))
emit("Figure_S1d_PC_scree_line_withoutPC1", function() draw_elbow(include_pc1 = FALSE))

cat("[10b] OK\n")
