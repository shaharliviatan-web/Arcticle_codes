#!/usr/bin/env Rscript
# 07c_manhattan_v3_and_mirror_tagstyle.R
# Manhattan figures in the SAME visual style as 07_manhattan.R (the v2 house style):
# rasterized point cloud composited under vector axes, alternating blue4/orange3
# chromosome colours, single solid black Bonferroni line, legend top-right.
#
# Produces two things 07b (the plain-mirror script) did not:
#   1. v3-ONLY stacked panels  -> 03_manhattan/v3_only/
#        manhattan_v3_<pheno>_pc<N>.{pdf,png}   4 traits stacked, one config per file
#   2. v2-vs-v3 MIRROR in house style -> 03_manhattan/v2_vs_v3/
#        mirror_tagstyle_<pheno>_pc<N>.{pdf,png}
#        v3 plotted upward, v2 downward on a shared genomic axis; each half carries
#        its own Bonferroni line (v3 6.0454, v2 6.7712).
#
# The earlier flat-colour mirrors from 07b are kept alongside as mirror_<pheno>_pc<N>.*
#
# Created 2026-08-16 for the v3 re-run review gate.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({
  library(data.table); library(png); library(parallel)
})

PIPE     <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
PS_V3    <- file.path(PIPE, "results", "emmax_ps")
PS_V2    <- file.path(PIPE, "results", "_archive", "_archive_v2_win50snp_2026-08-16", "emmax_ps")
BIM      <- file.path(PIPE, "intermediates", "snp_map.bim")
CMP      <- file.path(PIPE, "results", "comparison_v2_vs_v3", "03_manhattan")
OUT_V3   <- file.path(CMP, "v3_only")
OUT_MIR  <- file.path(CMP, "v2_vs_v3")
for (d in c(OUT_V3, OUT_MIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- House style constants (copied from 07_manhattan.R so the look matches) ----
PANEL_W_IN  <- 12
PANEL_H_IN  <- 3.6      # per trait panel; 4 stacked -> 14.4in tall figure
POINT_CEX   <- 0.85
RASTER_DPI  <- 900
PNG_OUT_DPI <- 600
CHR_COLS    <- c("blue4", "orange3")

BONF_V3 <- -log10(0.10 / 111017L)   # 6.0454
BONF_V2 <- -log10(0.10 / 590462L)   # 6.7712
TRAITS  <- c("betaglucan", "fiber", "protein", "starch")

# ---- Genomic map ----
cat("[07c] Loading snp_map.bim ...\n")
snp_map <- fread(BIM, select = c(1, 2, 4), col.names = c("CHR", "SNP", "BP"))
snp_map[, CHR := as.integer(sub("H$", "", CHR))]
setorder(snp_map, CHR, BP)
chr_max <- snp_map[, .(mx = max(BP)), by = CHR][order(CHR)]
chr_max[, off := cumsum(as.numeric(shift(mx, fill = 0)))]
snp_map <- merge(snp_map, chr_max[, .(CHR, off)], by = "CHR"); setorder(snp_map, CHR, BP)
snp_map[, x := as.numeric(BP) + off]
snp_map[, col := ifelse(CHR %% 2L == 0L, CHR_COLS[2], CHR_COLS[1])]
chr_cent <- snp_map[, .(c = (min(x) + max(x)) / 2), by = CHR][order(CHR)]
XLIM <- c(min(snp_map$x), max(snp_map$x))

read_cell <- function(dir, trait, pheno, pc) {
  f <- file.path(dir, sprintf("morexV3__%s__%s__pc%d.ps", trait, pheno, pc))
  if (!file.exists(f)) return(NULL)
  gw <- fread(f, select = 4, col.names = "P", showProgress = FALSE)
  if (nrow(gw) != nrow(snp_map))
    stop(sprintf("[07c] FAIL: %s has %d rows, snp_map has %d", basename(f), nrow(gw), nrow(snp_map)))
  d <- data.table(x = snp_map$x, col = snp_map$col, P = gw$P)
  d <- d[!is.na(P) & P > 0 & P <= 1]
  d[, nlp := -log10(P)]
  d[]
}

# Rasterize a point cloud to a transparent PNG and read it back (house technique).
raster_of <- function(x, y, ylim, w_in, h_in) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = w_in, height = h_in, units = "in", res = RASTER_DPI,
      bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(x, y, xlim = XLIM, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "", ylab = "", pch = 16, cex = POINT_CEX, col = attr(y, "col"))
  dev.off()
  img <- readPNG(tmp); unlink(tmp); img
}

# ---------------- v3-only stacked panels ----------------
draw_v3_panel <- function(d, trait, pheno, pc, img, ymax) {
  n_above <- sum(d$nlp >= BONF_V3)
  par(mar = c(4.2, 4.8, 2.8, 1.2), las = 1, cex.axis = 0.95, cex.lab = 1.0)
  plot(NA, xlim = XLIM, ylim = c(0, ymax), xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "Chromosome", ylab = expression(-log[10](italic(p))),
       main = sprintf("Manhattan: %s | %s | %d PCs   (v3, Bonf alpha=0.10 LD-pruned)",
                      tools::toTitleCase(trait), pheno, pc))
  rasterImage(img, XLIM[1], 0, XLIM[2], ymax, interpolate = FALSE)
  axis(2); axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE); box()
  abline(h = BONF_V3, lty = 1, lwd = 1.4, col = "black")
  legend("topright", bty = "n", cex = 0.8, inset = c(0.01, 0.02),
         legend = c(sprintf("Bonf (LD-pruned, n=111,017) = %.2f", BONF_V3),
                    sprintf("n above = %d", n_above)),
         col = c("black", NA), lty = c(1, NA), lwd = c(1.4, NA))
}

# ---------------- v2-vs-v3 mirror, house style ----------------
draw_mirror_panel <- function(trait, pheno, pc, img3, img2, ymax, n3, n2) {
  par(mar = c(4.2, 4.8, 2.8, 1.2), las = 1, cex.axis = 0.95, cex.lab = 1.0)
  plot(NA, xlim = XLIM, ylim = c(-ymax, ymax), xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "Chromosome", ylab = expression(-log[10](italic(p))),
       main = sprintf("%s | %s | %d PCs      v3 (up)  vs  v2 (down)",
                      tools::toTitleCase(trait), pheno, pc))
  rasterImage(img3, XLIM[1],     0, XLIM[2],  ymax, interpolate = FALSE)
  rasterImage(img2, XLIM[1], -ymax, XLIM[2],     0, interpolate = FALSE)
  yt <- pretty(c(0, ymax), 4); yt <- yt[yt <= ymax & yt > 0]
  axis(2, at = c(-rev(yt), 0, yt), labels = c(rev(yt), 0, yt))
  axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE); box()
  abline(h = 0, col = "grey35", lwd = 0.8)
  abline(h =  BONF_V3, lty = 1, lwd = 1.4, col = "black")
  abline(h = -BONF_V2, lty = 1, lwd = 1.4, col = "black")
  legend("topright", bty = "n", cex = 0.78, inset = c(0.01, 0.02),
         legend = c(sprintf("v3 Bonf = %.2f   (n above = %d)", BONF_V3, n3)),
         col = "black", lty = 1, lwd = 1.4)
  legend("bottomright", bty = "n", cex = 0.78, inset = c(0.01, 0.02),
         legend = c(sprintf("v2 Bonf = %.2f   (n above = %d)", BONF_V2, n2)),
         col = "black", lty = 1, lwd = 1.4)
}

for (ph in c("BLUP", "BLUE")) for (pc in c(3L, 5L, 10L)) {
  cat(sprintf("[07c] %s pc%d ...\n", ph, pc))

  dat <- mclapply(TRAITS, function(tr) {
    d3 <- read_cell(PS_V3, tr, ph, pc)
    d2 <- read_cell(PS_V2, tr, ph, pc)
    list(trait = tr, d3 = d3, d2 = d2)
  }, mc.cores = 4)

  # ---- v3-only figure ----
  ymax3 <- ceiling(max(vapply(dat, function(x) max(x$d3$nlp), numeric(1)),
                       BONF_V3, 6)) + 1
  imgs3 <- lapply(dat, function(x) {
    y <- x$d3$nlp; attr(y, "col") <- x$d3$col
    raster_of(x$d3$x, y, c(0, ymax3), PANEL_W_IN, PANEL_H_IN)
  })
  base3 <- sprintf("manhattan_v3_%s_pc%d", ph, pc)
  for (fmt in c("pdf", "png")) {
    f <- file.path(OUT_V3, paste0(base3, ".", fmt))
    if (fmt == "pdf") cairo_pdf(f, width = PANEL_W_IN, height = PANEL_H_IN * 4)
    else png(f, width = PANEL_W_IN, height = PANEL_H_IN * 4, units = "in",
             res = PNG_OUT_DPI, type = "cairo-png")
    par(mfrow = c(4, 1))
    for (i in seq_along(dat))
      draw_v3_panel(dat[[i]]$d3, dat[[i]]$trait, ph, pc, imgs3[[i]], ymax3)
    dev.off()
  }
  cat(sprintf("[07c] OK  v3_only/%s.{pdf,png}\n", base3))

  # ---- mirror figure ----
  ymaxm <- ceiling(max(
    vapply(dat, function(x) max(x$d3$nlp), numeric(1)),
    vapply(dat, function(x) max(x$d2$nlp), numeric(1)), BONF_V2, 6)) + 1
  imgs_m3 <- lapply(dat, function(x) {
    y <- x$d3$nlp; attr(y, "col") <- x$d3$col
    raster_of(x$d3$x, y, c(0, ymaxm), PANEL_W_IN, PANEL_H_IN)
  })
  imgs_m2 <- lapply(dat, function(x) {
    # draw v2 upward into its own raster, then place it flipped below the axis
    y <- x$d2$nlp; attr(y, "col") <- x$d2$col
    img <- raster_of(x$d2$x, y, c(0, ymaxm), PANEL_W_IN, PANEL_H_IN)
    img[rev(seq_len(dim(img)[1])), , , drop = FALSE]   # vertical flip
  })
  basem <- sprintf("mirror_tagstyle_%s_pc%d", ph, pc)
  for (fmt in c("pdf", "png")) {
    f <- file.path(OUT_MIR, paste0(basem, ".", fmt))
    if (fmt == "pdf") cairo_pdf(f, width = PANEL_W_IN, height = PANEL_H_IN * 4)
    else png(f, width = PANEL_W_IN, height = PANEL_H_IN * 4, units = "in",
             res = PNG_OUT_DPI, type = "cairo-png")
    par(mfrow = c(4, 1))
    for (i in seq_along(dat))
      draw_mirror_panel(dat[[i]]$trait, ph, pc, imgs_m3[[i]], imgs_m2[[i]], ymaxm,
                        sum(dat[[i]]$d3$nlp >= BONF_V3),
                        sum(dat[[i]]$d2$nlp >= BONF_V2))
    dev.off()
  }
  cat(sprintf("[07c] OK  v2_vs_v3/%s.{pdf,png}\n", basem))
  rm(dat, imgs3, imgs_m3, imgs_m2); gc(verbose = FALSE)
}

cat("\n---------------------------------\n")
cat(sprintf("[07c] CHECKPOINT: v3_only files  = %d (expected 12)\n",
            length(list.files(OUT_V3, pattern = "^manhattan_v3_.*\\.(pdf|png)$"))))
cat(sprintf("[07c] CHECKPOINT: tagstyle mirrors = %d (expected 12)\n",
            length(list.files(OUT_MIR, pattern = "^mirror_tagstyle_.*\\.(pdf|png)$"))))
cat("[07c] OK\n")
