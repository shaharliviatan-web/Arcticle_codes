#!/usr/bin/env Rscript
# make_figure_3.R — Fig. 3 of the TAG manuscript (Results ch. 2): GWAS Manhattan + QQ plots.
#
# Layout, as in the mini paper: one row per trait (a β-glucan, b fiber, c protein, d starch),
# Manhattan on the left, QQ on the right.
#
# Manhattan
#   * every SNP is one dot, and ALL dots are the same size (no enlarged lead SNPs);
#   * the members of each LD-clumped locus (02_loci_FINAL, 50 kb gap rule) are painted in one
#     colour per locus over an alternating-grey chromosome background, so the physical extent of
#     each association is visible. A locus reduced to its lead SNP by the gap rule is drawn like
#     any other locus (no open circle);
#   * no title, no legend, no locus labels (L01, L02 …): locus identities go to the Online
#     Resource table. Colours only separate neighbouring loci; they carry no identity, so the
#     8-colour palette is reused along the genome (neighbouring loci always differ);
#   * one horizontal line: the Bonferroni threshold (alpha 0.10 over the LD-pruned SNP count).
# QQ
#   * observed vs expected −log10(p) with the y = x line and λGC only (no title, no legend,
#     no threshold lines).
#
# TAG figure spec (10_USED_Paper_writing/TAG_requirements.md): 174 mm wide (full page width),
# height ≤ 234 mm, Arial-metric sans (Liberation Sans) at a TRUE 9 pt at final size (TAG range 8–12),
# lines ≥ 0.3 pt, RGB, 600 dpi (combination art). The point layer holds 4 × 7.1 M SNPs, so the
# figure is raster: drawn once to PNG (docx build), converted to TIFF (submission).
#
# Inputs (read-only, 01_USED_GWAS_V2_pipeline):
#   results/00_FINAL_BLUP_3PC/01_assoc/<trait>.assoc               SNP, P (7,110,996 rows)
#   results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/loci_summary.tsv  one row per locus
#   results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/loci_members.tsv  one row per member SNP
#   intermediates/morexV3_pruned_for_covs.prune.in                  LD-pruned SNPs (threshold)
#   results/tables/lambda_table.tsv                                 λGC cross-check
# Outputs: Figure_3/Fig3.tif (LZW), Figure_3/Fig3.png
#
# Created and approved 2026-09-22 (see README.md). Run: Rscript make_figure_3.R   (~4 min, ~10 GB RAM)
#
# Fixed 2026-09-27 (user decision) — TRUE TEXT SIZE, the same bug found in make_figure_4.R on
# 2026-09-24: layout() with three or more rows silently sets par(cex = 0.66) and nothing reset it,
# so the version approved 2026-09-22 as "12 pt" measured ~7.9 pt on the page, below TAG's 8 pt
# minimum, and the margin arithmetic (LINE_IN, in lines of PT) was wrong for the same reason.
# par(cex = 1) is now set after layout() and PT = 9 — a true 9 pt, matching Fig. 4, ~14% larger
# than the approved look. PT_CEX was rescaled (0.40 -> 0.35) so the SNP dots keep the approved
# size on the page. Nothing else changed: same data, same colours, same layout.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))
set.seed(1)                                   # QQ bulk subsample only

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
PIPE <- file.path(ROOT, "01_USED_GWAS_V2_pipeline")
SET  <- file.path(PIPE, "results", "00_FINAL_BLUP_3PC")
LOC  <- file.path(SET, "02_loci_FINAL", "tables")
OUT  <- file.path(ROOT, "08_USED_creating_figures", "Figure_3")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TRAITS    <- c("betaglucan", "fiber", "protein", "starch")
TRAIT_LAB <- c(betaglucan = "β-glucan", fiber = "Fiber", protein = "Protein", starch = "Starch")
N_SNP     <- 7110996L
N_PRUNED  <- length(readLines(file.path(PIPE, "intermediates", "morexV3_pruned_for_covs.prune.in")))
BONF      <- -log10(0.10 / N_PRUNED)          # 6.0454 for 111,017 pruned SNPs

# ── Style ─────────────────────────────────────────────────────────────────────
FONT     <- "Liberation Sans"                 # metric-compatible with Arial
PT       <- 9                                 # all lettering, TRUE size on the page (TAG 8-12); needs par(cex = 1) after layout()
DPI      <- 600
W_MM     <- 174                               # TAG full width
W_MAN_MM <- 128; W_QQ_MM <- W_MM - W_MAN_MM   # Manhattan 128 mm + QQ 46 mm
PLOT_IN  <- 1.30                              # plot-region height per row (in)
MAR_TOP  <- 1.5                               # lines: panel letter + trait
MAR_BOT  <- 1.5                               # lines: chromosome / tick labels
MAR_BOTX <- 2.9                               # bottom row: + axis title
LINE_IN  <- 1.2 * PT / 72                     # one margin line in inches
PT_CEX   <- 0.35                              # one dot size for every SNP, both plot types (0.40 x 12/9 x 0.66: keeps the approved dot size)
LWD      <- 0.75                              # ≥ 0.3 pt
BG_COLS  <- c("grey74", "grey56")             # odd / even chromosomes
# Validated categorical order (dataviz skill reference palette, light mode):
# adjacent-pair CVD ΔE ≥ 9.1, normal-vision ΔE ≥ 19.6.
PALETTE  <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
              "#e87ba4", "#008300", "#4a3aa7", "#e34948")
COL_QQ   <- "grey20"; COL_NULL <- "#e34948"

# ── Data ──────────────────────────────────────────────────────────────────────
S <- fread(file.path(LOC, "loci_summary.tsv"))
M <- fread(file.path(LOC, "loci_members.tsv"))
LAM_REF <- fread(file.path(PIPE, "results", "tables", "lambda_table.tsv"))[pheno_type == "BLUP" & n_PCs == 3]
calc_lambda <- function(p) median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)

snp <- NULL; G <- list()
for (tr in TRAITS) {
  a <- fread(file.path(SET, "01_assoc", paste0(tr, ".assoc")), showProgress = FALSE)
  stopifnot(nrow(a) == N_SNP, !anyNA(a$P), all(a$P > 0 & a$P <= 1))
  if (is.null(snp)) {
    snp <- a[, .(SNP)]
    snp[, c("chr", "bp") := tstrsplit(SNP, ":", fixed = TRUE)]
    snp[, `:=`(chr = as.integer(sub("H$", "", chr)), bp = as.numeric(bp))]
    stopifnot(!is.unsorted(snp$chr))
    cm <- snp[, .(len = max(bp)), by = chr][order(chr)]
    cm[, off := cumsum(shift(len, fill = 0))]
    snp[cm, x := bp + i.off, on = "chr"]
    BG_IDX <- list(which(snp$chr %% 2L == 1L), which(snp$chr %% 2L == 0L))   # odd / even
    CHR_MID <- snp[, .(mid = (min(x) + max(x)) / 2), by = chr]
    XLIM    <- c(0, sum(cm$len))
  } else stopifnot(identical(a$SNP, snp$SNP))  # same SNP order in every file

  nlp <- -log10(a$P)
  lam <- calc_lambda(a$P)
  stopifnot(abs(lam - LAM_REF[trait == tr, lambda_GC]) < 1e-9)

  # loci: colour in genome order; members painted, looked up by SNP id (same run as S)
  K <- S[trait == tr][order(chr, lead_bp)]
  K[, colour := PALETTE[(seq_len(.N) - 1L) %% length(PALETTE) + 1L]]
  MM <- M[trait == tr][K[, .(locus_id, colour)], on = "locus_id"]
  idx <- match(MM$SNP, snp$SNP); stopifnot(!anyNA(idx))
  MM[, `:=`(x = snp$x[idx], nlp = nlp[idx])]
  sig <- snp$SNP[nlp > BONF]
  stopifnot(all(sig %in% MM$SNP))              # every significant SNP sits in a painted locus

  # QQ: all points in the tail, a random 150k of the bulk (visually identical)
  o <- order(a$P); n <- length(o)
  q_idx <- c(seq_len(50000L), sort(sample(50001L:n, 150000L)))
  G[[tr]] <- list(nlp = nlp, MM = MM, lam = lam, n_loci = nrow(K), n_sig = length(sig),
                  qq_exp = -log10(ppoints(n))[q_idx], qq_obs = nlp[o][q_idx])
  cat(sprintf("[fig3] %-10s lambda %.4f | %2d loci | %4d painted SNPs | %2d significant | max -log10p %.2f\n",
              tr, lam, nrow(K), nrow(MM), length(sig), max(nlp)))
  rm(a, o); gc(verbose = FALSE)
}

# ── Drawing ───────────────────────────────────────────────────────────────────
mm2in <- function(mm) mm / 25.4
ROW_IN  <- PLOT_IN + (MAR_TOP + MAR_BOT) * LINE_IN
LAST_IN <- PLOT_IN + (MAR_TOP + MAR_BOTX) * LINE_IN
W_IN <- mm2in(W_MM); H_IN <- 3 * ROW_IN + LAST_IN
stopifnot(H_IN * 25.4 <= 234)

draw <- function() {
  layout(matrix(1:8, ncol = 2, byrow = TRUE),
         widths = c(W_MAN_MM, W_QQ_MM), heights = c(rep(ROW_IN, 3), LAST_IN))
  par(cex = 1)                                # undo layout()'s automatic 0.66 shrink
  par(family = FONT, las = 1, mgp = c(1.9, 0.45, 0), tcl = -0.25,
      cex.axis = 1, cex.lab = 1, lwd = LWD, xpd = FALSE)
  for (i in seq_along(TRAITS)) {
    tr <- TRAITS[i]; g <- G[[tr]]; last <- i == length(TRAITS)
    bot <- if (last) MAR_BOTX else MAR_BOT

    # Manhattan
    ymax <- ceiling(max(g$nlp, BONF)) + 0.5
    par(mar = c(bot, 3.0, MAR_TOP, 0.4))
    plot(NA, xlim = XLIM, ylim = c(0, ymax), xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "", ylab = expression(-log[10](italic(p))))
    for (k in 1:2) {                          # one call per grey: no per-point colour parsing
      j <- BG_IDX[[k]]; points(snp$x[j], g$nlp[j], pch = 16, cex = PT_CEX, col = BG_COLS[k])
    }
    abline(h = BONF, lwd = LWD, col = "black")
    points(g$MM$x, g$MM$nlp, pch = 16, cex = PT_CEX, col = g$MM$colour)
    axis(2, at = seq(0, ymax, 2), lwd = LWD)
    axis(1, at = CHR_MID$mid, labels = paste0(CHR_MID$chr, "H"), tick = FALSE, line = -0.2)
    box(bty = "l", lwd = LWD)
    if (last) title(xlab = "Chromosome", line = 1.7)
    # panel letter + trait: bottom-left of the label 0.06 in above the plot box, at the row's left edge
    text(grconvertX(0, "nfc", "user"),
         grconvertY(grconvertY(1, "npc", "inches") + 0.06, "inches", "user"),
         bquote(bold(.(letters[i])) ~~ .(TRAIT_LAB[[tr]])), adj = c(0, 0), xpd = NA)

    # QQ
    lim <- max(g$qq_exp, g$qq_obs) * 1.04
    par(mar = c(bot, 3.0, MAR_TOP, 0.5))
    plot(g$qq_exp, g$qq_obs, pch = 16, cex = PT_CEX, col = COL_QQ,
         xlim = c(0, lim), ylim = c(0, lim), xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "", ylab = expression(Observed ~ -log[10](italic(p))))
    abline(0, 1, col = COL_NULL, lwd = LWD)
    axis(1, at = seq(0, lim, 2), lwd = LWD); axis(2, at = seq(0, lim, 2), lwd = LWD)
    box(bty = "l", lwd = LWD)
    if (last) title(xlab = expression(Expected ~ -log[10](italic(p))), line = 1.7)
    text(0.04 * lim, 0.96 * lim, bquote(lambda[GC] == .(sprintf("%.3f", g$lam))), adj = c(0, 1))
  }
}

png(file.path(OUT, "Fig3.png"), width = W_IN, height = H_IN, units = "in", res = DPI,
    pointsize = PT, type = "cairo", family = FONT, bg = "white")
draw(); invisible(dev.off())
# TIFF for submission: same pixels as the PNG, RGB 8 bit/channel, LZW (TAG: RGB, ≥ 600 dpi)
stopifnot(system2("python3", c("-c", shQuote(sprintf(
  "from PIL import Image; Image.open('%s').convert('RGB').save('%s', compression='tiff_lzw', dpi=(%d,%d))",
  file.path(OUT, "Fig3.png"), file.path(OUT, "Fig3.tif"), DPI, DPI)))) == 0)
cat(sprintf("[fig3] OK -> %s  (%.0f x %.0f mm, %d dpi; threshold -log10p %.4f)\n",
            OUT, W_MM, H_IN * 25.4, DPI, BONF))
