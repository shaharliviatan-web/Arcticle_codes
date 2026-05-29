#!/usr/bin/env Rscript
# 09_comparison_views.R  (v2: Bonferroni alpha=0.10, no FDR)
# Compose the per-cell QQ and Manhattan PNGs into at-a-glance comparison contact sheets.
# Pure visual aid for the manual configuration pick -- NO logic, NO recommendation.
# Layout: 4 traits (rows) x 6 configs (cols: BLUP pc3/5/10, then BLUE pc3/5/10).
# Each source panel is already self-titled.
# Uses R (png + grid + gridExtra) to avoid ImageMagick's policy.xml resource ceilings.
#
# v2 outputs in results/comparison/:
#   qq_grid_6configs.{pdf,png}
#   manhattan_grid_6configs_BonfAll010.{pdf,png}      (one Bonf line each, alpha=0.10, all SNPs)
#   manhattan_grid_6configs_BonfPruned010.{pdf,png}   (one Bonf line each, alpha=0.10, LD-pruned)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(png)
  library(grid)
  library(gridExtra)
})

PIPE    <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
QQ_DIR  <- file.path(PIPE, "results", "qq")
MAN_DIR <- file.path(PIPE, "results", "manhattan")
OUT_DIR <- file.path(PIPE, "results", "comparison")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TRAITS  <- c("betaglucan", "fiber", "protein", "starch")
# Column order: BLUP pc3,5,10 then BLUE pc3,5,10
CONFIGS <- list(c("BLUP","pc3"), c("BLUP","pc5"), c("BLUP","pc10"),
                c("BLUE","pc3"), c("BLUE","pc5"), c("BLUE","pc10"))

# Build ordered file vectors (trait-major, config-minor = row-major)
qq_files          <- character(0)
man_files_all010  <- character(0)
man_files_prn010  <- character(0)
for (tr in TRAITS) for (cf in CONFIGS) {
  qq_files         <- c(qq_files,         file.path(QQ_DIR,  sprintf("qq__%s__%s__%s.png", tr, cf[1], cf[2])))
  man_files_all010 <- c(man_files_all010, file.path(MAN_DIR, sprintf("manhattan__%s__%s__%s__BonfAll010.png",    tr, cf[1], cf[2])))
  man_files_prn010 <- c(man_files_prn010, file.path(MAN_DIR, sprintf("manhattan__%s__%s__%s__BonfPruned010.png", tr, cf[1], cf[2])))
}
stopifnot(all(file.exists(qq_files)),
          all(file.exists(man_files_all010)),
          all(file.exists(man_files_prn010)))
cat(sprintf("[09] QQ panels: %d   Manhattan BonfAll010: %d   Manhattan BonfPruned010: %d (expect 24 each)\n",
            length(qq_files), length(man_files_all010), length(man_files_prn010)))

make_grid <- function(files, out_base, title, tile_w_in, tile_h_in, res = 150) {
  grobs <- lapply(files, function(f) rasterGrob(readPNG(f), interpolate = TRUE))
  arranged <- arrangeGrob(grobs = grobs, nrow = 4, ncol = 6,
                          top = textGrob(title, gp = gpar(fontsize = 16, fontface = "bold")))
  W <- tile_w_in * 6
  H <- tile_h_in * 4 + 0.4   # +0.4 for the title strip
  png_f <- file.path(OUT_DIR, paste0(out_base, ".png"))
  pdf_f <- file.path(OUT_DIR, paste0(out_base, ".pdf"))
  png(png_f, width = W, height = H, units = "in", res = res, type = "cairo-png")
  grid.draw(arranged); dev.off()
  cairo_pdf(pdf_f, width = W, height = H)
  grid.draw(arranged); dev.off()
  cat(sprintf("[09] %s  -> %.1f x %.1f in  (PNG %.1f MB, PDF %.1f MB)\n",
              out_base, W, H,
              file.info(png_f)$size/1024/1024, file.info(pdf_f)$size/1024/1024))
}

cat("[09] Building QQ grid...\n")
make_grid(qq_files, "qq_grid_6configs",
          "QQ plots  |  rows: betaglucan / fiber / protein / starch  |  cols: BLUP pc3,5,10  then BLUE pc3,5,10",
          tile_w_in = 2.6, tile_h_in = 2.6, res = 150)

cat("[09] Building Manhattan grid (Bonferroni alpha=0.10, all SNPs)...\n")
make_grid(man_files_all010, "manhattan_grid_6configs_BonfAll010",
          "Manhattan (Bonf alpha=0.10, all 7.1M SNPs)  |  rows: betaglucan / fiber / protein / starch  |  cols: BLUP pc3,5,10  then BLUE pc3,5,10",
          tile_w_in = 4.0, tile_h_in = 2.4, res = 150)

cat("[09] Building Manhattan grid (Bonferroni alpha=0.10, LD-pruned 590K)...\n")
make_grid(man_files_prn010, "manhattan_grid_6configs_BonfPruned010",
          "Manhattan (Bonf alpha=0.10, LD-pruned 590K)  |  rows: betaglucan / fiber / protein / starch  |  cols: BLUP pc3,5,10  then BLUE pc3,5,10",
          tile_w_in = 4.0, tile_h_in = 2.4, res = 150)

# Checkpoints
cat("\n---------------------------------\n")
expected <- c("qq_grid_6configs.png", "qq_grid_6configs.pdf",
              "manhattan_grid_6configs_BonfAll010.png",    "manhattan_grid_6configs_BonfAll010.pdf",
              "manhattan_grid_6configs_BonfPruned010.png", "manhattan_grid_6configs_BonfPruned010.pdf")
ok <- TRUE
for (f in expected) {
  p <- file.path(OUT_DIR, f)
  if (file.exists(p) && file.info(p)$size > 0) {
    cat(sprintf("[09] OK  %s  (%.1f MB)\n", f, file.info(p)$size/1024/1024))
  } else { cat(sprintf("[09] FAIL: %s missing/empty\n", f)); ok <- FALSE }
}
if (!ok) stop("[09] FAIL: some contact sheets missing")
cat("[09] OK: comparison contact sheets ready in results/comparison/\n")
