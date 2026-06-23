#!/usr/bin/env Rscript
# 10e_manhattan_gene_annotated.R
# Like 10d (highlighted 2x2 Manhattan, BLUP x 3 PCs) but additionally marks, on
# each trait panel, the GWAS peak from which the crosshap candidate gene(s) were
# identified: a dashed black vertical line at the lead/anchor SNP of that peak,
# with the gene symbol(s) labelled (italic) at the top of the line.
#
# Gene -> trait -> peak SNP (from 04_.../06_publication_figures/config/figure_genes.tsv
# and 04_.../00_config/gene_windows.tsv):
#   beta-glucan : AP2/ERF (3HG0299440)            @ 3H:537247620  (significant)
#   starch      : Pho (3HG0301750) + PHT (3HG0301710) @ 3H:546433616 (marginal)
#   fiber       : BAHD (7HG0642350)               @ 7H:14817657   (marginal)
#   protein     : none
#
# Modes:
#   Rscript 10e_manhattan_gene_annotated.R test   -> single beta-glucan panel
#   Rscript 10e_manhattan_gene_annotated.R full    -> 2x2 figure (all four traits)
#
# Output: results/publication_BonfOnly_BLUP_3PC/manhattan_gene_annotated/

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(png)
})

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else "test"
stopifnot(MODE %in% c("test", "full"))

PIPE     <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/01_USED_GWAS_V2_pipeline"
INTER    <- file.path(PIPE, "intermediates")
PS_DIR   <- file.path(PIPE, "results", "emmax_ps")
LEAD_TSV <- file.path(PIPE, "results", "publication_BonfOnly_BLUP_3PC", "tables", "lead_snps.tsv")
OUT_DIR  <- file.path(PIPE, "results", "publication_BonfOnly_BLUP_3PC", "manhattan_gene_annotated")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PHENO  <- "BLUP"; PC_TAG <- "pc3"
BONF   <- -log10(0.10 / 590462L)        # 6.7712
CHR_COLS <- c("blue4", "orange3")
SIG_COL  <- "#1B9E77"                    # significant: green filled circle
MARG_COL <- "#E41A1C"                    # marginal:    red filled circle
HL_PCH   <- 21
HL_CEX   <- 1.9
SHOW_HL  <- FALSE                        # plain Manhattan: do NOT colour-highlight selected SNPs
RASTER_DPI <- 900

# ---- enlarged fonts (same as 10d) ----
CEX_AXIS  <- 1.65
CEX_LAB   <- 2.0
CEX_TRAIT <- 2.3

# Display order for the 2x2 grid (row-major). Fiber sits BELOW Starch (right
# column) so the shared 7H peak can be bracketed by one box across both panels.
TRAITS       <- c("betaglucan", "starch", "protein", "fiber")
TRAIT_LABELS <- c(betaglucan = "β-glucan", fiber = "Fiber",
                  protein = "Protein", starch = "Starch")

MARGINAL <- list(
  betaglucan = c("7H:144534576", "5H:224439424", "6H:545947347"),
  fiber      = c("7H:14817657",  "1H:344520079"),
  protein    = character(0),
  starch     = c("6H:525776080", "7H:151110354", "3H:546433616", "7H:573606460")
)

# ============================================================================
# EDIT HERE  --  gene labels + label text size  (fast, one-line changes)
# ----------------------------------------------------------------------------
#   * To rename a gene  : change the `genes = c(...)` text for that trait.
#                         Two genes on one peak -> c("Pho", "PHT") (joined " / ").
#   * To resize ALL labels at once : change GENE_LABEL_CEX below.
#   * To resize ONE trait's label  : set its `cex = <number>` (NA = use global).
#   * `snp` is the peak coordinate (the dashed line position) -- normally leave as is.
# ============================================================================
GENE_LABEL_CEX <- 2.3    # global gene-label font size (match trait-label size CEX_TRAIT)
GENE_LINE_COL  <- "grey70"  # dashed peak-line colour (light, so it doesn't hide SNPs)
GENE_PEAKS <- list(
  betaglucan = list(snp = "3H:537247620", genes = c("AP2/ERF"),    cex = NA),
  starch     = list(snp = "3H:546433616", genes = c("Pho", "PHT"), cex = NA),
  fiber      = list(snp = "7H:14817657",  genes = c("BAHD"),       cex = NA)
)

# Shared 7H peak (~573.6 Mb) that appears in BOTH starch and fiber, only 154 bp
# apart. In the 2x2 figure (Fiber stacked under Starch) it is highlighted by a
# SINGLE dashed rectangle spanning both right-column panels, same colour as the
# gene-name peak line (GENE_LINE_COL).
SHOW_SHARED7H_BOX <- FALSE  # FALSE -> no shared-peak box (added manually in slides)
SHARED7H        <- list(starch = 573606460, fiber = 573606306)  # 7H positions (bp)
SHARED7H_HW_BP  <- 18e6     # half-width of the box (bp on either side of the peak)
BOX             <- new.env(parent = emptyenv())   # closures stash NDC corners here

lead <- fread(LEAD_TSV)
SIGNIF <- split(lead$SNP_id, lead$trait)

snp_map <- fread(file.path(INTER, "morexV3_290.bim"), header = FALSE,
                 col.names = c("CHR_raw", "SNP", "cm", "BP", "A1", "A2"))
snp_map[, CHR := as.integer(sub("H$", "", CHR_raw))]
stopifnot(all(snp_map$CHR %in% 1:7))

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

  # gene-peak annotation for this trait (if any)
  gp <- GENE_PEAKS[[tr]]
  gene_x <- NA_real_; gene_lab <- NULL; gene_cex <- GENE_LABEL_CEX
  if (!is.null(gp)) {
    pk <- d[SNP == gp$snp]
    stopifnot(nrow(pk) == 1)
    gene_x   <- pk$x
    gene_lab <- gp$genes
    if (!is.null(gp$cex) && !is.na(gp$cex)) gene_cex <- gp$cex
  }

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
    par(mar = c(5.2, 6.2, 6.4, 1.0), las = 1, cex.axis = CEX_AXIS, cex.lab = CEX_LAB,
        mgp = c(3.8, 0.9, 0))
    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "Chromosome", ylab = expression(-log[10](italic(p))), main = "")
    # gene-peak dashed vertical line, extended a little above the top box.
    # Drawn BEFORE the point cloud so the SNP points sit on top of the line.
    line_top <- ylim[2] + 0.05 * diff(ylim)
    if (!is.na(gene_x)) {
      segments(gene_x, ylim[1], gene_x, line_top, col = GENE_LINE_COL, lty = 2, lwd = 2.0, xpd = NA)
    }
    rasterImage(img, xlim[1], ylim[1], xlim[2], ylim[2], interpolate = FALSE)
    axis(2, lwd = 1.4)
    axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE)
    box(lwd = 1.4)
    abline(h = BONF, col = "black", lty = 1, lwd = 1.8)
    # optional colour-highlighting of selected SNPs (off by default -> plain Manhattan)
    if (SHOW_HL) {
      if (nrow(marg) > 0) points(marg$x, marg$nlp, pch = HL_PCH, bg = MARG_COL, col = "black", cex = HL_CEX, lwd = 1.2)
      if (nrow(sig)  > 0) points(sig$x,  sig$nlp,  pch = HL_PCH, bg = SIG_COL,  col = "black", cex = HL_CEX, lwd = 1.2)
    }
    # gene name(s) just above the top of the line (italic); >1 gene = stacked lines
    if (!is.na(gene_x) && length(gene_lab) > 0) {
      lab <- paste(gene_lab, collapse = "\n")
      adj_x <- if (gene_x > mean(xlim)) 1.0 else 0.0   # keep inside the panel
      text(gene_x, line_top + 0.015 * diff(ylim), labels = lab, adj = c(adj_x, 0),
           font = 4, cex = gene_cex, xpd = NA)
    }
    # trait label (top-left)
    text(xlim[1] + 0.015 * diff(xlim), ylim[2] * 0.96, labels = TRAIT_LABELS[[tr]],
         adj = c(0, 1), font = 2, cex = CEX_TRAIT)
    # stash NDC corners for the shared-7H cross-panel box (drawn later, grid only)
    if (tr %in% names(SHARED7H)) {
      pk7 <- d[CHR == 7L][which.min(abs(BP - as.numeric(SHARED7H[[tr]])))]
      bxl <- grconvertX(pk7$x - SHARED7H_HW_BP, "user", "ndc")
      bxr <- grconvertX(pk7$x + SHARED7H_HW_BP, "user", "ndc")
      if (tr == "starch") {            # box TOP edge sits just above the starch peak
        ytop <- min(ylim[2] - 0.08 * diff(ylim), pk7$nlp + 1.3)
        assign("starch", list(xl = bxl, xr = bxr, y = grconvertY(ytop, "user", "ndc")), envir = BOX)
      } else {                          # fiber: box BOTTOM edge sits just below the fiber peak
        ybot <- max(ylim[1] + 0.08 * diff(ylim), pk7$nlp - 1.3)
        assign("fiber", list(xl = bxl, xr = bxr, y = grconvertY(ybot, "user", "ndc")), envir = BOX)
      }
    }
  }
}

if (MODE == "test") {
  cat("[10e] TEST: beta-glucan single panel\n")
  panel <- build_trait("betaglucan")
  base <- "manhattan_gene_annotated_TEST_betaglucan"
  cairo_pdf(file.path(OUT_DIR, paste0(base, ".pdf")), width = 10, height = 5, pointsize = 12)
  panel(); dev.off()
  png(file.path(OUT_DIR, paste0(base, ".png")), width = 10, height = 5, units = "in",
      res = 400, pointsize = 12, type = "cairo-png")
  panel(); dev.off()
  cat(sprintf("[10e] wrote %s.{pdf,png}\n", base))
} else {
  cat("[10e] FULL: building 4 panels...\n")
  panels <- setNames(lapply(TRAITS, build_trait), TRAITS)

  # ---- individual per-trait panels (10x5 each) ----
  for (tr in TRAITS) {
    base <- sprintf("manhattan_gene_annotated_%s", tr)
    cairo_pdf(file.path(OUT_DIR, paste0(base, ".pdf")), width = 10, height = 5, pointsize = 12)
    panels[[tr]](); dev.off()
    png(file.path(OUT_DIR, paste0(base, ".png")), width = 10, height = 5, units = "in",
        res = 400, pointsize = 12, type = "cairo-png")
    panels[[tr]](); dev.off()
    cat(sprintf("[10e] wrote %s.{pdf,png}\n", base))
  }

  # ---- combined 2x2 figure ----
  W_IN <- 20; H_IN <- 10
  draw_grid <- function() {
    if (exists("starch", envir = BOX)) rm("starch", envir = BOX)
    if (exists("fiber",  envir = BOX)) rm("fiber",  envir = BOX)
    par(mfrow = c(2, 2), oma = c(0, 0, 0, 0))
    for (tr in TRAITS) panels[[tr]]()
    # ONE dashed rectangle around the shared 7H peak across the stacked starch/fiber panels
    if (SHOW_SHARED7H_BOX && exists("starch", envir = BOX) && exists("fiber", envir = BOX)) {
      s <- get("starch", envir = BOX); f <- get("fiber", envir = BOX)
      xl <- mean(c(s$xl, f$xl)); xr <- mean(c(s$xr, f$xr))
      par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      rect(xl, f$y, xr, s$y, border = GENE_LINE_COL, lty = 2, lwd = 2.4, xpd = NA)
    }
  }
  base <- "manhattan_gene_annotated_2x2_BLUP_pc3"
  cairo_pdf(file.path(OUT_DIR, paste0(base, ".pdf")), width = W_IN, height = H_IN, pointsize = 12)
  draw_grid(); dev.off()
  png(file.path(OUT_DIR, paste0(base, ".png")), width = W_IN, height = H_IN, units = "in",
      res = 400, pointsize = 12, type = "cairo-png")
  draw_grid(); dev.off()
  cat(sprintf("[10e] wrote %s.{pdf,png}\n", base))
}
cat("[10e] OK\n")
