#!/usr/bin/env Rscript
# 10d_manhattan_highlight_2x2.R
# Combine the four highlighted per-trait Manhattan plots (BLUP x 3 PCs) into a
# single 2x2 publication figure with enlarged fonts (axis labels, tick labels,
# trait names). Layout matches the reference: top-left beta-glucan, top-right
# Starch, bottom-left Fiber, bottom-right Protein. No QQ plots.
#
# Reuses the highlighting logic of 10c_manhattan_highlight.R:
#   - significant : Bonferroni-passing SNPs (green filled circles), from lead_snps.tsv
#   - marginal    : curated near-Bonferroni SNPs (red filled circles)
#
# Output: results/publication_BonfOnly_BLUP_3PC/manhattan_highlighted/
#         manhattan_highlight_2x2_BLUP_pc3.{pdf,png}

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(png)
})

PIPE     <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/01_USED_GWAS_V2_pipeline"
INTER    <- file.path(PIPE, "intermediates")
PS_DIR   <- file.path(PIPE, "results", "emmax_ps")
LEAD_TSV <- file.path(PIPE, "results", "publication_BonfOnly_BLUP_3PC", "tables", "lead_snps.tsv")
OUT_DIR  <- file.path(PIPE, "results", "publication_BonfOnly_BLUP_3PC", "manhattan_highlighted")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PHENO  <- "BLUP"; PC_TAG <- "pc3"
BONF   <- -log10(0.10 / 590462L)        # 6.7712
CHR_COLS <- c("blue4", "orange3")
SIG_COL  <- "#1B9E77"                    # significant: green filled circle
MARG_COL <- "#E41A1C"                    # marginal:    red filled circle
HL_PCH   <- 21
HL_CEX   <- 1.9
RASTER_DPI <- 900

# ---- enlarged fonts (relative to 10c) ----
CEX_AXIS  <- 1.65   # tick labels
CEX_LAB   <- 2.0    # axis titles ("Chromosome", -log10(p))
CEX_TRAIT <- 2.3    # trait name in panel

# Display order to match the reference image (row-major in a 2x2 mfrow grid)
TRAITS       <- c("betaglucan", "starch", "fiber", "protein")
TRAIT_LABELS <- c(betaglucan = "β-glucan", fiber = "Fiber",
                  protein = "Protein", starch = "Starch")

MARGINAL <- list(
  betaglucan = c("7H:144534576", "5H:224439424", "6H:545947347"),
  fiber      = c("7H:14817657",  "1H:344520079"),
  protein    = character(0),
  starch     = c("6H:525776080", "7H:151110354", "3H:546433616", "7H:573606460")
)

lead <- fread(LEAD_TSV)
SIGNIF <- split(lead$SNP_id, lead$trait)

snp_map <- fread(file.path(INTER, "morexV3_290.bim"), header = FALSE,
                 col.names = c("CHR_raw", "SNP", "cm", "BP", "A1", "A2"))
snp_map[, CHR := as.integer(sub("H$", "", CHR_raw))]
stopifnot(all(snp_map$CHR %in% 1:7))

# Build the raster point cloud + plotting closure for one trait
build_trait <- function(tr) {
  gw <- fread(file.path(PS_DIR, sprintf("morexV3__%s__%s__%s.ps", tr, PHENO, PC_TAG)),
              header = FALSE, col.names = c("SNP", "beta", "SE", "P"), showProgress = FALSE)
  stopifnot(nrow(gw) == nrow(snp_map))
  d <- data.table(CHR = snp_map$CHR, SNP = snp_map$SNP, BP = snp_map$BP, P = as.numeric(gw$P))
  d <- d[!is.na(P) & P > 0 & P <= 1]
  d[, nlp := -log10(P)]

  setorder(d, CHR, BP)
  chr_max <- d[, .(mx = max(BP)), by = CHR][order(CHR)]
  chr_max[, off := cumsum(as.numeric(shift(mx, fill = 0)))]
  d <- merge(d, chr_max[, .(CHR, off)], by = "CHR"); setorder(d, CHR, BP)
  d[, x := as.numeric(BP) + off]
  d[, col := ifelse(CHR %% 2L == 0L, CHR_COLS[2], CHR_COLS[1])]
  chr_cent <- d[, .(c = (min(x) + max(x)) / 2), by = CHR][order(CHR)]
  xlim <- c(min(d$x), max(d$x))
  ymax <- ceiling(max(max(d$nlp), BONF, 6)) + 1
  ylim <- c(0, ymax)

  sig_ids  <- SIGNIF[[tr]]; if (is.null(sig_ids)) sig_ids <- character(0)
  marg_ids <- MARGINAL[[tr]]
  sig  <- d[SNP %in% sig_ids]
  marg <- d[SNP %in% marg_ids]
  stopifnot(nrow(sig) == length(sig_ids), nrow(marg) == length(marg_ids))

  # raster point cloud (2:1 panel aspect)
  tmp <- tempfile(fileext = ".png", tmpdir = Sys.getenv("TMPDIR"))
  png(tmp, width = 10, height = 5, units = "in", res = RASTER_DPI,
      bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(d$x, d$nlp, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "", ylab = "", pch = 16, cex = 0.85, col = d$col)
  dev.off()
  img <- readPNG(tmp); unlink(tmp)

  function() {
    par(mar = c(5.2, 6.2, 2.0, 1.0), las = 1, cex.axis = CEX_AXIS, cex.lab = CEX_LAB,
        mgp = c(3.8, 0.9, 0))
    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "Chromosome", ylab = expression(-log[10](italic(p))), main = "")
    rasterImage(img, xlim[1], ylim[1], xlim[2], ylim[2], interpolate = FALSE)
    axis(2, lwd = 1.4)
    axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE)
    box(lwd = 1.4)
    abline(h = BONF, col = "black", lty = 1, lwd = 1.8)
    if (nrow(marg) > 0) points(marg$x, marg$nlp, pch = HL_PCH, bg = MARG_COL, col = "black", cex = HL_CEX, lwd = 1.2)
    if (nrow(sig)  > 0) points(sig$x,  sig$nlp,  pch = HL_PCH, bg = SIG_COL,  col = "black", cex = HL_CEX, lwd = 1.2)
    text(xlim[1] + 0.015 * diff(xlim), ylim[2] * 0.96, labels = TRAIT_LABELS[[tr]],
         adj = c(0, 1), font = 2, cex = CEX_TRAIT)
  }
}

cat("[10d] building 4 panels...\n")
panels <- lapply(TRAITS, build_trait)

# Full-figure dimensions: 2x2 of 10x5 panels
W_IN <- 20; H_IN <- 10

draw_grid <- function() {
  par(mfrow = c(2, 2), oma = c(0, 0, 0, 0))
  for (p in panels) p()
}

base <- "manhattan_highlight_2x2_BLUP_pc3"
cairo_pdf(file.path(OUT_DIR, paste0(base, ".pdf")), width = W_IN, height = H_IN, pointsize = 12)
draw_grid(); dev.off()
png(file.path(OUT_DIR, paste0(base, ".png")), width = W_IN, height = H_IN, units = "in",
    res = 400, pointsize = 12, type = "cairo-png")
draw_grid(); dev.off()

cat(sprintf("[10d] wrote %s.{pdf,png} in %s\n", base, OUT_DIR))
cat("[10d] OK\n")
