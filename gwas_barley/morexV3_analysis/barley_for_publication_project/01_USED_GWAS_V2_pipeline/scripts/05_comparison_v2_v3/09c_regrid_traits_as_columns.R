#!/usr/bin/env Rscript
# 09c_regrid_traits_as_columns.R
# Re-tile the existing Manhattan panels into a TRANSPOSED contact sheet:
#
#   columns = trait          (betaglucan, fiber, protein, starch)   -> 4 cols
#   rows    = config, in order BLUP 3, BLUP 5, BLUP 10,
#                              BLUE 3, BLUE 5, BLUE 10              -> 6 rows
#
# i.e. one trait runs straight down a single column, so the 3/5/10 PC and
# BLUP/BLUE variants of that trait sit stacked together.
#
# This is the transpose of 09b's sheets (4 traits as rows x 6 configs as cols);
# both layouts are kept. Nothing is deleted or overwritten -- outputs carry the
# "_traitcols" suffix.
#
# No re-rendering: every panel PNG already exists, this only re-arranges them.
#   v3 panels      : results/manhattan/manhattan__<trait>__<pheno>__<pc>__<tag>.png
#   mirror panels  : 03_manhattan/v2_vs_v3/_tiles_mirror{,_floor2}/mirror__*.png
#
# Outputs:
#   03_manhattan/v3_only/manhattan_grid_traitcols_v3_<tag>.{pdf,png}
#   03_manhattan/v2_vs_v3/manhattan_grid_traitcols_mirror.{pdf,png}
#   03_manhattan/v2_vs_v3/manhattan_grid_traitcols_mirror_floor2.{pdf,png}
#
# Created 2026-08-16.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(png); library(grid); library(gridExtra) })

PIPE    <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
MAN_DIR <- file.path(PIPE, "results", "manhattan")
CMP     <- file.path(PIPE, "results", "comparison_v2_vs_v3", "03_manhattan")
OUT_V3  <- file.path(CMP, "v3_only")
OUT_MIR <- file.path(CMP, "v2_vs_v3")

TRAITS  <- c("betaglucan", "fiber", "protein", "starch")          # -> columns
CONFIGS <- list(c("BLUP","pc3"), c("BLUP","pc5"), c("BLUP","pc10"),
                c("BLUE","pc3"), c("BLUE","pc5"), c("BLUE","pc10"))  # -> rows

# arrangeGrob fills ROW-MAJOR, so emit config-major: for each config row, all 4 traits.
order_files <- function(fn) unlist(lapply(CONFIGS, function(cf)
  vapply(TRAITS, function(tr) fn(tr, cf[1], cf[2]), character(1))))

make_grid <- function(files, out_dir, out_base, title,
                      tile_w_in = 4.5, tile_h_in = 3.0, res = 150) {
  miss <- files[!file.exists(files)]
  if (length(miss)) { cat(sprintf("[09c] SKIP %s: %d panels missing (first: %s)\n",
                                  out_base, length(miss), basename(miss[1]))); return(invisible()) }
  grobs <- lapply(files, function(f) rasterGrob(readPNG(f), interpolate = TRUE))
  arranged <- arrangeGrob(grobs = grobs, nrow = 6, ncol = 4,
                          top = textGrob(title, gp = gpar(fontsize = 15, fontface = "bold")))
  W <- tile_w_in * 4; H <- tile_h_in * 6 + 0.4
  png_f <- file.path(out_dir, paste0(out_base, ".png"))
  pdf_f <- file.path(out_dir, paste0(out_base, ".pdf"))
  png(png_f, width = W, height = H, units = "in", res = res, type = "cairo-png")
  grid.draw(arranged); dev.off()
  cairo_pdf(pdf_f, width = W, height = H)
  grid.draw(arranged); dev.off()
  cat(sprintf("[09c] OK  %s -> %.0f x %.0f in  (PNG %.1f MB, PDF %.1f MB)\n", out_base, W, H,
              file.info(png_f)$size/1024/1024, file.info(pdf_f)$size/1024/1024))
}

SUB <- "cols: betaglucan / fiber / protein / starch   |   rows: BLUP pc3,5,10 then BLUE pc3,5,10"

# ---- v3-only, both Bonferroni variants ----
for (tag in c("BonfPruned010", "BonfAll010")) {
  f <- order_files(function(tr, ph, pc)
        file.path(MAN_DIR, sprintf("manhattan__%s__%s__%s__%s.png", tr, ph, pc, tag)))
  ttl <- sprintf("v3 Manhattan (%s)  |  %s",
                 if (tag == "BonfPruned010") "Bonf alpha=0.10, LD-pruned 111,017"
                 else "Bonf alpha=0.10, all 7.1M SNPs", SUB)
  make_grid(f, OUT_V3, sprintf("manhattan_grid_traitcols_v3_%s", tag), ttl)
}

# ---- mirrors, both the house-style and the floor2 variant ----
for (v in c("", "_floor2")) {
  tdir <- file.path(OUT_MIR, paste0("_tiles_mirror", v))
  f <- order_files(function(tr, ph, pc)
        file.path(tdir, sprintf("mirror__%s__%s__%s.png", tr, ph, pc)))
  ttl <- sprintf("v3 (up) vs v2 (down)%s  |  %s",
                 if (v == "_floor2") "  [-log10p >= 2 shown]" else "", SUB)
  make_grid(f, OUT_MIR, paste0("manhattan_grid_traitcols_mirror", v), ttl)
}

cat("[09c] OK\n")
