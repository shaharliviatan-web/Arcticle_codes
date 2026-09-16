#!/usr/bin/env Rscript
# 09b_manhattan_grids_v3_and_mirror.R
# 4 x 6 Manhattan contact sheets in the same layout as
# 09_comparison_views.R / manhattan_grid_6configs_BonfPruned010.png:
#   rows = betaglucan / fiber / protein / starch
#   cols = BLUP pc3, pc5, pc10  then  BLUE pc3, pc5, pc10
# so every trait x config cell is visible in one image.
#
# Produces two sheets:
#   1. v3 ONLY   -> 03_manhattan/v3_only/manhattan_grid_6configs_v3_BonfPruned010.{pdf,png}
#      Tiles the per-cell PNGs already written by 07_manhattan.R (no recompute).
#   2. v3 vs v2 MIRRORED -> 03_manhattan/v2_vs_v3/manhattan_grid_6configs_mirror.{pdf,png}
#      Needs per-cell mirror panels, which do not exist yet, so they are rendered
#      first into 03_manhattan/v2_vs_v3/_tiles_mirror/ (kept: they are useful on
#      their own and make re-tiling free).
#
# Mirror panels use the house style of 07_manhattan.R: rasterized point cloud,
# alternating blue4/orange3 chromosomes, solid black Bonferroni line.
# v3 is drawn upward against its own line (6.0454); v2 downward against its own
# line (6.7712).
#
# Tile rendering is deliberately at TILE_DPI (not the 600/900 dpi used for
# standalone figures) because each tile occupies ~4 x 2.4 in of a 150 dpi sheet;
# rendering higher only costs time and disk.
#
# Created 2026-08-16 for the v3 re-run review gate.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({
  library(data.table); library(png); library(grid); library(gridExtra); library(parallel)
})

PIPE     <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
PS_V3    <- file.path(PIPE, "results", "emmax_ps")
PS_V2    <- file.path(PIPE, "results", "_archive", "_archive_v2_win50snp_2026-08-16", "emmax_ps")
MAN_DIR  <- file.path(PIPE, "results", "manhattan")
BIM      <- file.path(PIPE, "intermediates", "snp_map.bim")
CMP      <- file.path(PIPE, "results", "comparison_v2_vs_v3", "03_manhattan")
OUT_V3   <- file.path(CMP, "v3_only")
OUT_MIR  <- file.path(CMP, "v2_vs_v3")
TILE_MIR <- file.path(OUT_MIR, "_tiles_mirror")
for (d in c(OUT_V3, OUT_MIR, TILE_MIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

TRAITS  <- c("betaglucan", "fiber", "protein", "starch")
CONFIGS <- list(c("BLUP","pc3"), c("BLUP","pc5"), c("BLUP","pc10"),
                c("BLUE","pc3"), c("BLUE","pc5"), c("BLUE","pc10"))

BONF_V3 <- -log10(0.10 / 111017L)
BONF_V2 <- -log10(0.10 / 590462L)
CHR_COLS <- c("blue4", "orange3")
TILE_W_IN <- 12; TILE_H_IN <- 7; TILE_DPI <- 200; POINT_CEX <- 0.85

# MIRROR_FLOOR: drop SNPs below this -log10(p) from the mirror panels.
#   0 (default) = house style, every point drawn from the axis up.
#   2           = "sheet" variant. In the 4 x 6 contact sheet each mirror panel is
#                 squeezed to ~4 x 2.4 in AND spans a doubled y-range (-ymax..+ymax),
#                 so the dense 0-2 bulk of both halves collapses into one solid slab
#                 and the peaks lose contrast. Clipping the bulk restores it.
# Set via env var so both variants come from one script:  MIRROR_FLOOR=2 Rscript ...
MIRROR_FLOOR <- as.numeric(Sys.getenv("MIRROR_FLOOR", "0"))
SUFFIX       <- if (MIRROR_FLOOR > 0) sprintf("_floor%g", MIRROR_FLOOR) else ""
TILE_MIR     <- paste0(TILE_MIR, SUFFIX)
dir.create(TILE_MIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[09b] MIRROR_FLOOR = %g   tiles -> %s\n", MIRROR_FLOOR, basename(TILE_MIR)))

# ---------------- Shared grid assembler (same as 09_comparison_views.R) ----------------
make_grid <- function(files, out_dir, out_base, title, tile_w_in = 4.0,
                      tile_h_in = 2.4, res = 150) {
  stopifnot(all(file.exists(files)))
  grobs <- lapply(files, function(f) rasterGrob(readPNG(f), interpolate = TRUE))
  arranged <- arrangeGrob(grobs = grobs, nrow = 4, ncol = 6,
                          top = textGrob(title, gp = gpar(fontsize = 16, fontface = "bold")))
  W <- tile_w_in * 6; H <- tile_h_in * 4 + 0.4
  png_f <- file.path(out_dir, paste0(out_base, ".png"))
  pdf_f <- file.path(out_dir, paste0(out_base, ".pdf"))
  png(png_f, width = W, height = H, units = "in", res = res, type = "cairo-png")
  grid.draw(arranged); dev.off()
  cairo_pdf(pdf_f, width = W, height = H)
  grid.draw(arranged); dev.off()
  cat(sprintf("[09b] %s -> %.0f x %.0f in  (PNG %.1f MB, PDF %.1f MB)\n", out_base, W, H,
              file.info(png_f)$size/1024/1024, file.info(pdf_f)$size/1024/1024))
}

# ---------------- 1. v3-only sheet, from existing per-cell PNGs ----------------
cat("[09b] Building v3-only grid from existing results/manhattan/ PNGs ...\n")
for (tag in c("BonfPruned010", "BonfAll010")) {
  f <- unlist(lapply(TRAITS, function(tr) vapply(CONFIGS, function(cf)
        file.path(MAN_DIR, sprintf("manhattan__%s__%s__%s__%s.png", tr, cf[1], cf[2], tag)),
        character(1))))
  miss <- f[!file.exists(f)]
  if (length(miss)) { cat(sprintf("[09b] SKIP %s: %d panels missing\n", tag, length(miss))); next }
  ttl <- if (tag == "BonfPruned010")
    "v3 Manhattan (Bonf alpha=0.10, LD-pruned 111,017)  |  rows: betaglucan / fiber / protein / starch  |  cols: BLUP pc3,5,10  then BLUE pc3,5,10"
  else
    "v3 Manhattan (Bonf alpha=0.10, all 7.1M SNPs)  |  rows: betaglucan / fiber / protein / starch  |  cols: BLUP pc3,5,10  then BLUE pc3,5,10"
  make_grid(f, OUT_V3, sprintf("manhattan_grid_6configs_v3_%s", tag), ttl)
}

# ---------------- 2. Per-cell mirror tiles ----------------
cat("[09b] Loading snp_map.bim ...\n")
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

read_nlp <- function(dir, trait, pheno, pc) {
  f <- file.path(dir, sprintf("morexV3__%s__%s__%s.ps", trait, pheno, pc))
  if (!file.exists(f)) return(NULL)
  gw <- fread(f, select = 4, col.names = "P", showProgress = FALSE)
  if (nrow(gw) != nrow(snp_map)) stop(sprintf("[09b] FAIL: row mismatch in %s", basename(f)))
  keep <- !is.na(gw$P) & gw$P > 0 & gw$P <= 1
  list(x = snp_map$x[keep], col = snp_map$col[keep], nlp = -log10(gw$P[keep]))
}

raster_layer <- function(x, y, col, ylim) {
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = TILE_W_IN, height = TILE_H_IN, units = "in", res = TILE_DPI,
      bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(x, y, xlim = XLIM, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "", ylab = "", pch = 16, cex = POINT_CEX, col = col)
  dev.off()
  img <- readPNG(tmp); unlink(tmp); img
}

make_tile <- function(i) {
  tr <- TILE_JOBS$trait[i]; ph <- TILE_JOBS$pheno[i]; pc <- TILE_JOBS$pc[i]
  out <- file.path(TILE_MIR, sprintf("mirror__%s__%s__%s.png", tr, ph, pc))
  if (file.exists(out) && file.info(out)$size > 0) return(sprintf("cached %s", basename(out)))

  d3 <- read_nlp(PS_V3, tr, ph, pc); d2 <- read_nlp(PS_V2, tr, ph, pc)
  # Hit counts are taken BEFORE any floor is applied, so the legend is always exact.
  n3 <- sum(d3$nlp >= BONF_V3); n2 <- sum(d2$nlp >= BONF_V2)
  ymax <- ceiling(max(d3$nlp, d2$nlp, BONF_V2, 6)) + 1
  ymin <- MIRROR_FLOOR
  if (MIRROR_FLOOR > 0) {
    k3 <- d3$nlp >= MIRROR_FLOOR; d3 <- list(x = d3$x[k3], col = d3$col[k3], nlp = d3$nlp[k3])
    k2 <- d2$nlp >= MIRROR_FLOOR; d2 <- list(x = d2$x[k2], col = d2$col[k2], nlp = d2$nlp[k2])
  }
  img3 <- raster_layer(d3$x, d3$nlp, d3$col, c(ymin, ymax))
  img2 <- raster_layer(d2$x, d2$nlp, d2$col, c(ymin, ymax))
  img2 <- img2[rev(seq_len(dim(img2)[1])), , , drop = FALSE]   # flip for the lower half

  png(out, width = TILE_W_IN, height = TILE_H_IN, units = "in", res = TILE_DPI,
      type = "cairo-png")
  par(mar = c(4.5, 4.8, 3, 1.2), las = 1, cex.axis = 0.95, cex.lab = 1.0)
  plot(NA, xlim = XLIM, ylim = c(-ymax, ymax), xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "Chromosome", ylab = expression(-log[10](italic(p))),
       main = sprintf("%s | %s | %s PCs    v3 up / v2 down",
                      tools::toTitleCase(tr), ph, sub("^pc", "", pc)))
  rasterImage(img3, XLIM[1],  ymin, XLIM[2],  ymax, interpolate = FALSE)
  rasterImage(img2, XLIM[1], -ymax, XLIM[2], -ymin, interpolate = FALSE)
  yt <- pretty(c(0, ymax), 4); yt <- yt[yt > 0 & yt <= ymax]
  axis(2, at = c(-rev(yt), 0, yt), labels = c(rev(yt), 0, yt))
  axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE); box()
  abline(h = 0, col = "grey35", lwd = 0.8)
  abline(h =  BONF_V3, lty = 1, lwd = 1.4, col = "black")
  abline(h = -BONF_V2, lty = 1, lwd = 1.4, col = "black")
  legend("topright", bty = "n", cex = 0.8, inset = c(0.01, 0.02), lty = 1, lwd = 1.4,
         col = "black", legend = sprintf("v3 Bonf %.2f  (n=%d)", BONF_V3, n3))
  legend("bottomright", bty = "n", cex = 0.8, inset = c(0.01, 0.02), lty = 1, lwd = 1.4,
         col = "black", legend = sprintf("v2 Bonf %.2f  (n=%d)", BONF_V2, n2))
  dev.off()
  rm(d3, d2, img3, img2); gc(verbose = FALSE)
  sprintf("wrote %s", basename(out))
}

TILE_JOBS <- rbindlist(lapply(TRAITS, function(tr)
  rbindlist(lapply(CONFIGS, function(cf) data.table(trait = tr, pheno = cf[1], pc = cf[2])))))
cat(sprintf("[09b] Rendering %d per-cell mirror tiles (4 at a time) ...\n", nrow(TILE_JOBS)))
res <- mclapply(seq_len(nrow(TILE_JOBS)), function(i)
         tryCatch(make_tile(i), error = function(e) paste("ERROR:", conditionMessage(e))),
       mc.cores = 4)
for (r in res) cat(sprintf("[09b]   %s\n", r))
bad <- grep("^ERROR", unlist(res), value = TRUE)
if (length(bad)) stop(sprintf("[09b] FAIL: %d tiles errored:\n%s", length(bad), paste(bad, collapse = "\n")))

# ---------------- 3. Mirror sheet ----------------
mir_files <- unlist(lapply(TRAITS, function(tr) vapply(CONFIGS, function(cf)
  file.path(TILE_MIR, sprintf("mirror__%s__%s__%s.png", tr, cf[1], cf[2])), character(1))))
make_grid(mir_files, OUT_MIR, paste0("manhattan_grid_6configs_mirror", SUFFIX),
  sprintf("v3 (up) vs v2 (down)%s  |  rows: betaglucan / fiber / protein / starch  |  cols: BLUP pc3,5,10  then BLUE pc3,5,10",
          if (MIRROR_FLOOR > 0) sprintf("  [-log10p >= %g shown]", MIRROR_FLOOR) else ""))

cat("\n---------------------------------\n")
cat(sprintf("[09b] CHECKPOINT: mirror tiles = %d (expected 24)\n",
            length(list.files(TILE_MIR, pattern = "^mirror__.*\\.png$"))))
cat("[09b] OK\n")
