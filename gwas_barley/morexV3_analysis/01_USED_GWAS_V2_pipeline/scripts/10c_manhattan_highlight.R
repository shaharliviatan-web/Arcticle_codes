#!/usr/bin/env Rscript
# 10c_manhattan_highlight.R
# Per-trait Manhattan plots (chosen config: BLUP x 3 PCs) that HIGHLIGHT, with
# distinct shape+colour, the SNPs classified as:
#   - significant : Bonferroni-passing (-log10p >= 6.7712), read from lead_snps.tsv
#   - marginal    : curated near-Bonferroni SNPs (hard-coded list below)
# so the two classes are distinguishable from the Manhattan plot alone.
#
# Output: results/publication_BonfOnly_BLUP_3PC/manhattan_highlighted/
#         manhattan_highlight_BLUP_pc3__{trait}.{pdf,png}   (4 traits x 2 formats)
#   also: results/publication_BonfOnly_BLUP_3PC/tables/marginal_snps.tsv
#
# Style: same raster-composite Manhattan as step 10 (raster point cloud + vector
# axes/line), no title. Significant SNPs = green filled circles; marginal SNPs =
# red filled circles (both black-outlined, same shape). No per-class legend.

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
HL_PCH   <- 21                           # same shape for both classes
HL_CEX   <- 1.5
W_IN <- 10; H_IN <- 5; RASTER_DPI <- 900; PNG_DPI <- 600

TRAITS       <- c("betaglucan", "fiber", "protein", "starch")
TRAIT_LABELS <- c(betaglucan = "β-glucan", fiber = "Fiber",
                  protein = "Protein", starch = "Starch")

# Curated marginal SNPs (BLUP x 3 PCs); 7H:144534576 (β-glucan) and 7H:14817657
# (fiber) were Bonferroni-significant under BLUE x 3 PCs but fall just below 6.7712 under BLUP.
MARGINAL <- list(
  betaglucan = c("7H:144534576", "5H:224439424", "6H:545947347"),
  fiber      = c("7H:14817657",  "1H:344520079"),
  protein    = character(0),
  starch     = c("6H:525776080", "7H:151110354", "3H:546433616", "7H:573606460")
)

# Significant SNPs per trait, from the published lead-SNP table
lead <- fread(LEAD_TSV)
SIGNIF <- split(lead$SNP_id, lead$trait)

# SNP map
snp_map <- fread(file.path(INTER, "morexV3_290.bim"), header = FALSE,
                 col.names = c("CHR_raw", "SNP", "cm", "BP", "A1", "A2"))
snp_map[, CHR := as.integer(sub("H$", "", CHR_raw))]
stopifnot(all(snp_map$CHR %in% 1:7))

draw_trait <- function(tr) {
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
  # sanity: every highlighted SNP must be found
  stopifnot(nrow(sig) == length(sig_ids), nrow(marg) == length(marg_ids))

  # raster point cloud
  tmp <- tempfile(fileext = ".png", tmpdir = Sys.getenv("TMPDIR"))
  png(tmp, width = W_IN, height = H_IN, units = "in", res = RASTER_DPI,
      bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(d$x, d$nlp, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "", ylab = "", pch = 16, cex = 0.85, col = d$col)
  dev.off()
  img <- readPNG(tmp); unlink(tmp)

  draw <- function() {
    par(mar = c(4.2, 4.6, 1.6, 0.8), las = 1, cex.axis = 0.95, cex.lab = 1.05,
        mgp = c(2.7, 0.6, 0))
    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "Chromosome", ylab = expression(-log[10](italic(p))), main = "")
    rasterImage(img, xlim[1], ylim[1], xlim[2], ylim[2], interpolate = FALSE)
    axis(2); axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE); box()
    abline(h = BONF, col = "black", lty = 1, lwd = 1.6)
    # Significant = red, marginal = purple (same circle shape). No per-class legend.
    if (nrow(marg) > 0) points(marg$x, marg$nlp, pch = HL_PCH, bg = MARG_COL, col = "black", cex = HL_CEX, lwd = 1.1)
    if (nrow(sig)  > 0) points(sig$x,  sig$nlp,  pch = HL_PCH, bg = SIG_COL,  col = "black", cex = HL_CEX, lwd = 1.1)
    # trait label (top-left)
    text(xlim[1] + 0.015 * diff(xlim), ylim[2] * 0.96, labels = TRAIT_LABELS[[tr]],
         adj = c(0, 1), font = 2, cex = 1.15)
  }

  base <- sprintf("manhattan_highlight_%s_%s__%s", PHENO, PC_TAG, tr)
  cairo_pdf(file.path(OUT_DIR, paste0(base, ".pdf")), width = W_IN, height = H_IN, pointsize = 10)
  draw(); dev.off()
  png(file.path(OUT_DIR, paste0(base, ".png")), width = W_IN, height = H_IN, units = "in",
      res = PNG_DPI, pointsize = 10, type = "cairo-png")
  draw(); dev.off()
  cat(sprintf("[10c] %-11s sig=%d marg=%d -> %s.{pdf,png}\n", tr, nrow(sig), nrow(marg), base))

  # marginal rows for the table (NA-safe if a trait has none)
  if (nrow(marg) == 0L) return(NULL)
  data.table(trait = tr, SNP_id = marg$SNP, chr = paste0(marg$CHR, "H"),
             position_bp = marg$BP, neg_log10_p = round(marg$nlp, 4),
             p_value = signif(marg$P, 4), class = "marginal")
}

cat(sprintf("[10c] Bonferroni = %.4f ; highlighting significant + marginal SNPs (BLUP x 3 PCs)\n", BONF))
marg_tbl <- rbindlist(lapply(TRAITS, draw_trait))
setorder(marg_tbl, trait, -neg_log10_p)
TAB_OUT <- file.path(PIPE, "results", "publication_BonfOnly_BLUP_3PC", "tables", "marginal_snps.tsv")
fwrite(marg_tbl, TAB_OUT, sep = "\t")
cat(sprintf("[10c] wrote %s (%d marginal SNPs)\n", TAB_OUT, nrow(marg_tbl)))
cat("[10c] OK\n")
