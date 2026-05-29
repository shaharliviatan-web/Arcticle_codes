#!/usr/bin/env Rscript
# 09b_bonferroni_only_grid.R
# Companion to 09_comparison_views.R: build a 4x6 Manhattan comparison grid
# showing ONLY the 2 Bonferroni reference lines (Bonf_all + Bonf_pruned), with
# NO FDR line. Does NOT touch any existing per-cell or grid files.
#
# Per-cell tiles are re-rendered at the SAME quality as step 07 (12x7 in,
# raster 900 dpi, cex 0.85, vector axes/labels) -- just with the FDR line
# omitted -- so the grid PDF/PNG matches the quality and size of the existing
# FDR 0.05 / 0.10 grids. Per-cell PNGs go to a temp dir and are cleaned at end.
#
# Output:
#   results/comparison/manhattan_grid_6configs_BonfOnly.png
#   results/comparison/manhattan_grid_6configs_BonfOnly.pdf

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(png)
  library(grid)
  library(gridExtra)
  library(parallel)
})

PIPE      <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
INTER_DIR <- file.path(PIPE, "intermediates")
PS_DIR    <- file.path(PIPE, "results", "emmax_ps")
OUT_DIR   <- file.path(PIPE, "results", "comparison")
BIM_FILE  <- file.path(INTER_DIR, "snp_map.bim")
PRUNE_IN  <- file.path(INTER_DIR, "morexV3_pruned_for_covs.prune.in")

# --- Same per-cell rendering style as step 07 ---
PLOT_W_IN   <- 12      # device width (inches)
PLOT_H_IN   <- 7       # device height (inches)
POINT_CEX   <- 0.85
RASTER_DPI  <- 900     # point-cloud raster resolution
PNG_OUT_DPI <- 600     # output PNG dpi
CHR_COLS    <- c("blue4", "orange3")

# --- Grid layout (same tile geometry as step 09) ---
TILE_W_IN <- 4.0
TILE_H_IN <- 2.4

# --- Bonferroni thresholds ---
n_all       <- 7110996
n_pruned    <- length(readLines(PRUNE_IN))
bonf_all    <- -log10(0.05 / n_all)
bonf_pruned <- -log10(0.05 / n_pruned)
cat(sprintf("[09b] Bonf_all=%.3f  Bonf_pruned=%.3f\n", bonf_all, bonf_pruned))

# --- SNP map ---
snp_map <- fread(BIM_FILE, header = FALSE, select = c(1, 2, 4),
                 col.names = c("CHR_raw", "SNP", "BP"))
snp_map[, CHR := as.integer(sub("H$", "", sub("^chr", "", as.character(CHR_raw))))]
stopifnot(all(snp_map$CHR %in% 1:7))

# --- .ps files + parsing ---
ps_files <- list.files(PS_DIR, pattern = "^morexV3__.+__.+__pc\\d+\\.ps$", full.names = TRUE)
stopifnot(length(ps_files) == 24)
parts <- as.data.table(do.call(rbind, strsplit(basename(ps_files), "__|\\.", perl = FALSE)))
setnames(parts, c("prefix", "trait", "pheno_type", "pcs_str", "ext"))
parts[, pcs := as.integer(sub("^pc", "", pcs_str))]
parts[, path := ps_files]

# Temp tile dir under our designated TMPDIR
TMP_DIR <- tempfile(pattern = "bonf_tiles_", tmpdir = Sys.getenv("TMPDIR"))
dir.create(TMP_DIR)

# --- Per-cell tile render: matches step 07 exactly except NO FDR line ---
render_tile <- function(i) {
  row <- parts[i, ]
  gw <- fread(row$path, select = 4, col.names = "P", showProgress = FALSE)
  if (nrow(gw) != nrow(snp_map)) stop("row count mismatch")
  d <- data.table(CHR = snp_map$CHR, SNP = snp_map$SNP, BP = snp_map$BP, P = gw$P)
  d <- d[!is.na(P) & P > 0 & P <= 1]
  d[, nlp := -log10(P)]

  # Manhattan layout (same as step 07)
  setorder(d, CHR, BP)
  chr_max <- d[, .(mx = max(BP)), by = CHR][order(CHR)]
  chr_max[, off := cumsum(as.numeric(shift(mx, fill = 0)))]
  d <- merge(d, chr_max[, .(CHR, off)], by = "CHR"); setorder(d, CHR, BP)
  d[, x := as.numeric(BP) + off]
  d[, col := ifelse(CHR %% 2L == 0L, CHR_COLS[2], CHR_COLS[1])]
  chr_cent <- d[, .(c = (min(x) + max(x)) / 2), by = CHR][order(CHR)]
  xlim <- c(min(d$x), max(d$x))
  max_obs <- max(d$nlp, na.rm = TRUE)
  ymax <- ceiling(max(max_obs, bonf_all, bonf_pruned, 6, na.rm = TRUE)) + 1
  ylim <- c(0, ymax)
  pheno_disp <- sprintf("%s | %s | %d PCs",
                        tools::toTitleCase(row$trait), row$pheno_type, row$pcs)

  # Render point cloud ONCE to transparent raster (same as step 07)
  pts_tmp <- tempfile(fileext = ".png", tmpdir = TMP_DIR)
  png(pts_tmp, width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
      res = RASTER_DPI, bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(d$x, d$nlp, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "", ylab = "", pch = 16, cex = POINT_CEX, col = d$col)
  dev.off()
  img <- readPNG(pts_tmp)
  unlink(pts_tmp)

  # Composite: vector axes + raster points + 2 Bonferroni lines (NO FDR LINE)
  out <- file.path(TMP_DIR,
                   sprintf("tile__%s__%s__pc%d.png", row$trait, row$pheno_type, row$pcs))
  draw <- function() {
    title <- sprintf("Manhattan: %s   (Bonferroni only)", pheno_disp)
    par(mar = c(4.5, 4.8, 3, 1.2), las = 1, cex.axis = 0.95, cex.lab = 1.0)
    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "Chromosome", ylab = expression(-log[10](italic(p))), main = title)
    rasterImage(img, xlim[1], ylim[1], xlim[2], ylim[2], interpolate = FALSE)
    axis(2)
    axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE)
    box()
    abline(h = bonf_all,    lty = 1, lwd = 1.4, col = "black")
    abline(h = bonf_pruned, lty = 2, lwd = 1.2, col = "black")
    legend("topright",
           legend = c(sprintf("Bonf_all = %.2f",    bonf_all),
                      sprintf("Bonf_pruned = %.2f", bonf_pruned)),
           col = c("black", "black"), lty = c(1, 2), lwd = 1.2,
           cex = 0.8, bty = "n", inset = c(0.01, 0.02))
  }
  png(out, width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
      res = PNG_OUT_DPI, type = "cairo-png")
  draw(); dev.off()
  out
}

cat("[09b] Rendering 24 per-cell tiles at step-07 quality (12x7 in, 900 dpi raster)...\n")
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(parts)), render_tile, mc.cores = 8)
cat(sprintf("[09b] tile render elapsed: %.1f sec\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

err <- which(vapply(res, inherits, logical(1), what = "try-error"))
if (length(err) > 0) {
  for (i in err) { cat(sprintf("[09b] FAIL at cell %d: %s\n", i, basename(parts$path[i]))); print(res[[i]]) }
  stop(sprintf("[09b] FAIL: %d tiles errored", length(err)))
}

# --- Reorder into trait-major, config-minor (rows: 4 traits, cols: BLUP pc3/5/10 then BLUE pc3/5/10) ---
TRAITS  <- c("betaglucan", "fiber", "protein", "starch")
CONFIGS <- list(c("BLUP","pc3"), c("BLUP","pc5"), c("BLUP","pc10"),
                c("BLUE","pc3"), c("BLUE","pc5"), c("BLUE","pc10"))
ordered <- character(0)
for (tr in TRAITS) for (cf in CONFIGS) {
  ordered <- c(ordered, file.path(TMP_DIR, sprintf("tile__%s__%s__%s.png", tr, cf[1], cf[2])))
}
stopifnot(all(file.exists(ordered)))

# --- Compose grid (same approach as step 09) ---
cat("[09b] Composing grid (this embeds 24 full-res rasters; ~3-5 min)...\n")
t1 <- Sys.time()
grobs <- lapply(ordered, function(f) rasterGrob(readPNG(f), interpolate = TRUE))
arranged <- arrangeGrob(grobs = grobs, nrow = 4, ncol = 6,
                        top = textGrob("Manhattan (Bonferroni only)  |  rows: betaglucan / fiber / protein / starch  |  cols: BLUP pc3,5,10  then BLUE pc3,5,10",
                                       gp = gpar(fontsize = 16, fontface = "bold")))
W <- TILE_W_IN * 6
H <- TILE_H_IN * 4 + 0.4

png_f <- file.path(OUT_DIR, "manhattan_grid_6configs_BonfOnly.png")
pdf_f <- file.path(OUT_DIR, "manhattan_grid_6configs_BonfOnly.pdf")
png(png_f, width = W, height = H, units = "in", res = 150, type = "cairo-png")
grid.draw(arranged); dev.off()
cairo_pdf(pdf_f, width = W, height = H)
grid.draw(arranged); dev.off()

cat(sprintf("[09b] grid compose elapsed: %.1f sec\n",
            as.numeric(difftime(Sys.time(), t1, units = "secs"))))
cat(sprintf("\n[09b] Wrote %s  (%.1f MB)\n[09b] Wrote %s  (%.1f MB)\n",
            basename(png_f), file.info(png_f)$size / 1024 / 1024,
            basename(pdf_f), file.info(pdf_f)$size / 1024 / 1024))

unlink(TMP_DIR, recursive = TRUE)
cat("[09b] OK: BonfOnly grid ready (matched step-09 quality); per-cell temp tiles cleaned\n")
