#!/usr/bin/env Rscript
# make_figure_4.R — Fig. 4 of the TAG manuscript (Results ch. 4, optional): haplotype structure
# of the GDSL esterase/lipase HORVU.MOREX.r3.7HG0729030 at the shared 7H fiber/starch signal,
# with the elite cultivars drawn as aligned barcodes.
#
# Built on make_figure_3.R (then make_figure_4.R): same drawing code, palette, lettering (true 9 pt), line widths, row
# height and column widths (violins 70 mm | barcodes 104 mm), so Figs. 3 and 4 (then 4 and 5)
# read as one set. Differences from Fig. 3, all decided by the user 2026-09-24:
#   1. LAYOUT — one gene, two traits. Fiber and starch give IDENTICAL haplotype groups (the same
#      212 | 34 accessions; checked below), so their barcodes would be identical. The violins are
#      stacked on the left (a fiber, b starch) and ONE barcode (c) spans both rows on the right.
#   2. STATISTIC — the label above each violin is the raw Kruskal-Wallis P with eta-squared, not a
#      BH q: this gene is a single pre-specified test outside the step-04 BH family.
#   3. RED TRIANGLES under the three genome-wide significant SNP columns (fiber lead
#      7H:573,606,306; starch lead 7H:573,606,460; starch 7H:573,606,491), with a legend entry.
#      The fourth haplotype-defining SNP (7H:573,606,282, -log10P 2.35 / 3.13) is deliberately
#      NOT marked. Kept from the exploratory figure of the replacement analysis, which marked all
#      four marker-group SNPs.
#
# Source analysis: 03_01_7H_branch_Starch_Fiber_shared_signal_explore/ (the branch root since 2026-09-24;
# until then it lived in the subfolder REPLACEMENT_ANALYSIS_MGmin3_eps0.9)
# (crosshap at MGmin = 3, epsilon = 0.9, gene +/- 1 kb; step 04's run_crosshap.R unmodified).
# Only the `shared_sites` version is drawn, as in Fig. 3: a column exists only where the wild and
# the elite call sets both have a record with identical REF and ALT, so no genotype is assumed.
# Elite lines are NOT assigned to a haplotype group — crosshap never saw them.
#
# Wild consensus = step 07's rule (07_.../scripts/03_build_matrices.R): the raw per-gene VCF,
# genotype codes 0/1/2, per-SNP majority over the group's members ignoring missing, exact ties
# -> REF when REF is among the tied states (step 07's rule since 2026-09-24; this gene has no ties).
#
# TAG figure spec (10_USED_Paper_writing/TAG_requirements.md): 174 mm wide, height <= 234 mm,
# Liberation Sans (Arial-metric) at a true 9 pt (par(cex = 1) after layout(), as Fig. 3 since
# 2026-09-24), lines >= 0.3 pt, RGB, 600 dpi. No title in the image.
#
# Inputs (read-only), R = 03_01_7H_branch_Starch_Fiber_shared_signal_explore:
#   R/Cache/<trait>/<gene>/MGmin_3/eps_0.9/HapObject.rds    groups, phenotypes, marker groups
#   R/results/tables/mgmin3_gene_results.tsv                KW P, eta2, group sizes (cross-check)
#   R/results/tables/pairwise_group_tests.tsv               bracket statistics (Wilcoxon, Holm)
#   R/results/tables/Table_site_overlap.tsv                 shared-site count (cross-check)
#   R/results/tables/Table_allele_concordance.tsv           REF/ALT check, re-asserted here
#   R/results/tables/Table_elite_genotypes_wide__<gene>.tsv elite states (cross-check)
#   R/work/raw/<trait>/<gene>.vcf.gz                        wild raw genotypes, gene +/- 1 kb
#   R/work/elite/7HG0729030.elite.vcf.gz                    elite genotypes (DivBrowse, MorexV3)
#   R/archive_v1_MGmin2_eps0.6/results/tables/signal_snps.tsv  genome-wide significance per SNP
#     (an MGmin-independent output of the archived first version, still valid)
#   07_USED_.../config/elite_lines.tsv                      the five cultivars, in display order
# Output: Figure_4/Fig4.{png,tif}  (png for the docx build, tif for submission)
#
# Created 2026-09-24. Run: Rscript make_figure_4.R   (seconds)
#
# RENUMBERED 2026-09-30 (S. Hübner's comments on the figures; user decision): old Fig. 2 was dissolved,
# so this figure is now Fig. 4 (it was Fig. 5) and this script was make_figure_5.R, writing
# Figure_5/Fig5.*. Only the output folder, file names, log tags and this header changed; the
# re-rendered image is pixel-identical to the approved Fig5 (checked 2026-09-30). Hand-over:
# 10_USED_Paper_writing/new_publishing_paper/build/STRUCTURE_CHANGES_2026-10.md

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(vcfR) })

ROOT   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
BRANCH <- file.path(ROOT, "03_01_7H_branch_Starch_Fiber_shared_signal_explore")
REPL   <- BRANCH                               # the MGmin = 3 analysis is the branch root (2026-09-24)
ARCH   <- file.path(BRANCH, "archive_v1_MGmin2_eps0.6")   # archived first version (MGmin 2)
STEP07 <- file.path(ROOT, "07_USED_elite_lines_compariosn_to_wild_lines")
OUT    <- file.path(ROOT, "08_USED_creating_figures", "Figure_4")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

GENE  <- "HORVU.MOREX.r3.7HG0729030"
MGMIN <- 3; EPS <- 0.9
TRAITS <- list(
  list(trait = "fiber",  label = "GDSL", trait_lab = "fiber"),
  list(trait = "starch", label = "GDSL", trait_lab = "starch")
)

# ── Style (identical to make_figure_3.R) ──────────────────────────────────────
FONT   <- "Liberation Sans"
PT     <- 9                                  # TRUE size on the page; needs par(cex = 1) after layout()
DPI    <- 600
W_MM   <- 174
H_MAX  <- 234
LWD    <- 0.75
LINE_IN <- 1.2 * PT / 72

CODE_REF <- 0; CODE_ALT <- 1; CODE_HET <- 2
GT_COL <- c("0" = "#FFFACD", "1" = "#2F4F4F", "2" = "#C46210")
GT_LAB <- c("0" = "Reference", "1" = "Alternate", "2" = "Heterozygous")
COL_MISS <- "grey70"
TILE_H   <- 0.72
HAP_COL  <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")
SIG_COL  <- "red"                            # triangles, as in the exploratory figure
SIG_PCH  <- 17

# ── Data ──────────────────────────────────────────────────────────────────────
tab <- function(f) fread(file.path(REPL, "results", "tables", f))
GR    <- tab("mgmin3_gene_results.tsv")[MGmin == MGMIN & epsilon == EPS]
PW    <- tab("pairwise_group_tests.tsv")[eps == EPS]
OVL   <- tab("Table_site_overlap.tsv")
CONC  <- tab("Table_allele_concordance.tsv")
EWIDE <- tab(paste0("Table_elite_genotypes_wide__", GENE, ".tsv"))
SIGT  <- fread(file.path(ARCH, "results", "tables", "signal_snps.tsv"))
ELITE <- fread(file.path(STEP07, "config", "elite_lines.tsv"), skip = "line_name")

stopifnot(!any(CONC$status %in% c("swapped", "ref_differs")))   # REF/ALT mean the same base

gt_code <- function(gt) {                    # step 07's coding
  out <- rep(NA_real_, length(gt))
  out[gt %in% c("0/0", "0|0")] <- CODE_REF
  out[gt %in% c("1/1", "1|1")] <- CODE_ALT
  out[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- CODE_HET
  out
}
read_gt <- function(path) {
  v   <- suppressMessages(read.vcfR(path, verbose = FALSE))
  fix <- as.data.frame(getFIX(v), stringsAsFactors = FALSE)
  gt  <- extract.gt(v, element = "GT", as.numeric = FALSE)
  key <- paste(sub("^chr", "", fix$CHROM), fix$POS, fix$REF, fix$ALT, sep = ":")
  m   <- apply(gt, 2, gt_code); m <- matrix(as.numeric(m), nrow = nrow(gt), dimnames = list(key, colnames(gt)))
  t(m)                                        # rows = samples, cols = sites
}
consensus_row <- function(m) apply(m, 2, function(x) {   # step 07's rule (ties -> REF, 2026-09-24)
  x <- x[!is.na(x)]; if (!length(x)) return(NA_real_)
  tb <- table(x); top <- names(tb)[tb == max(tb)]
  if (length(top) > 1) return(if (as.character(CODE_REF) %in% top) CODE_REF else NA_real_)
  as.numeric(top)
})

W <- read_gt(file.path(REPL, "work", "raw", "fiber", paste0(GENE, ".vcf.gz")))
stopifnot(identical(W, read_gt(file.path(REPL, "work", "raw", "starch", paste0(GENE, ".vcf.gz")))))
E <- read_gt(file.path(REPL, "work", "elite", "7HG0729030.elite.vcf.gz"))
stopifnot(all(ELITE$sample_id %in% rownames(E)))
E <- E[ELITE$sample_id, , drop = FALSE]; rownames(E) <- ELITE$line_name

shared <- intersect(colnames(W), colnames(E))
spos   <- as.integer(sub("^[^:]+:([0-9]+):.*", "\\1", shared))
shared <- shared[order(spos)]; spos <- sort(spos)
stopifnot(length(shared) == OVL$n_shared)
EL <- E[, shared, drop = FALSE]

# elite states re-derived from the VCF must equal the published table
st <- function(x) unname(ifelse(is.na(x), "Missing", GT_LAB[as.character(x)]))
stopifnot(identical(EWIDE$site_key, shared),
          all(vapply(ELITE$line_name, function(l) identical(st(EL[l, ]), EWIDE[[l]]), logical(1))))

# Haplotype groups per trait, and the check that justifies ONE barcode for both traits
P <- list(); groups <- list()
for (tt in TRAITS) {
  ho  <- readRDS(file.path(REPL, "Cache", tt$trait, GENE, paste0("MGmin_", MGMIN),
                           paste0("eps_", EPS), "HapObject.rds"))
  hh  <- ho$HapObject[[paste0("Haplotypes_MGmin", MGMIN, "_E", EPS)]]
  ind <- as.data.table(hh$Indfile)[, .(Ind = as.character(Ind), hap = as.character(hap), Pheno = as.numeric(Pheno))]
  ind <- ind[hap != "0" & !is.na(Pheno)]
  gs  <- ind[, .(n = .N, mean_pheno = mean(Pheno)), by = hap][order(hap)]
  g   <- GR[trait == tt$trait]
  stopifnot(nrow(g) == 1L, g$haplotype_groups == nrow(gs),
            identical(as.integer(gs$n), as.integer(strsplit(g$group_sizes, "|", fixed = TRUE)[[1]])))
  vf  <- as.data.frame(hh$Varfile)
  br  <- PW[trait == tt$trait][, .(group1 = g1, group2 = g2, y.position = y, p.adj.signif = sym)]
  P[[tt$trait]] <- c(tt, list(ind = ind, gs = gs, br = br, p = g$kw_p_raw, eta2 = g$eta_squared,
                              mg = vf$ID[vf$MGs != "0"]))
  groups[[tt$trait]] <- setNames(ind$hap, ind$Ind)[order(ind$Ind)]
  cat(sprintf("[fig4] %-6s %d groups (%s) | KW P %.2e | eta2 %.3f\n", tt$trait, nrow(gs),
              paste(gs$n, collapse = "|"), g$kw_p_raw, g$eta_squared))
}
stopifnot(identical(groups$fiber, groups$starch))            # same accessions, same groups

gs   <- P$fiber$gs
cons <- t(vapply(gs$hap, function(h) {
  ids <- names(groups$fiber)[groups$fiber == h]
  stopifnot(all(ids %in% rownames(W)))
  consensus_row(W[ids, shared, drop = FALSE])
}, numeric(length(shared))))
rownames(cons) <- gs$hap
stopifnot(!any(cons %in% CODE_HET))                           # wild consensus: REF/ALT/NA only

# Genome-wide significant SNPs (either trait) — the triangles
sig_pos <- SIGT[sig_fiber == TRUE | sig_starch == TRUE, pos]
stopifnot(length(sig_pos) == 3L, all(sig_pos %in% spos),
          all(paste0("7H:", sig_pos) %in% P$fiber$mg))          # all three are in the defining group
sig_col <- match(sig_pos, spos)
cat(sprintf("[fig4] %d shared SNPs | significant SNP columns: %s | %d elite rows\n",
            length(shared), paste(sig_col, collapse = ", "), nrow(EL)))

drawn <- sort(unique(c(as.vector(cons), as.vector(EL)))); drawn <- drawn[!is.na(drawn)]
has_missing <- anyNA(cons) || anyNA(EL)

# ── Drawing helpers (from make_figure_3.R) ────────────────────────────────────
mm2in <- function(mm) mm / 25.4

p_expr <- function(p, eta2) {                # Kruskal-Wallis P, journal style: 2.1 x 10^-3
  e <- floor(log10(p)); mant <- p / 10^e
  bquote(italic(P) == .(sprintf("%.1f", mant)) %*% 10^.(e) * "," ~~ italic(eta)^2 == .(sprintf("%.3f", eta2)))
}

violin_panel <- function(p, ylab, cex_ax = 1) {
  gs  <- p$gs
  ind <- p$ind
  k   <- nrow(gs)

  # Densities are computed FIRST and the y range is taken from them, not from the
  # observations: the first draft truncated each density at the group's extreme values,
  # which cut the violins off flat at top and bottom (2026-09-23). `cut = 2` lets each
  # shape close naturally two bandwidths beyond the data without growing a long tail.
  D <- lapply(gs$hap, function(h) {
    v <- ind$Pheno[ind$hap == h]
    density(v, bw = "nrd0", adjust = 1, cut = 2, n = 512)
  })
  dlo <- min(vapply(D, function(d) min(d$x), 0))
  dhi <- max(vapply(D, function(d) max(d$x), 0))

  # Brackets are stacked one text line apart ABOVE the violins, sized in inches (2026-09-24).
  # Step 07's y.position values were in data units and let the labels sit on the next bracket
  # once the text was drawn at its true size; they now set only the stacking order.
  nb      <- nrow(p$br)
  step_in <- 0.85 * LINE_IN * cex_ax                # one bracket per ~text line
  gap_in  <- 0.06
  plot.new()
  pin_h <- par("pin")[2]
  pad   <- 0.03 * (dhi - dlo)
  core  <- (dhi + pad) - (dlo - pad)
  top_in <- if (nb) gap_in + nb * step_in else 0
  ru    <- core * pin_h / (pin_h - top_in)           # full y range, data units
  ylim  <- c(dlo - pad, dlo - pad + ru)
  plot.window(xlim = c(0.4, k + 0.6), ylim = ylim, xaxs = "i", yaxs = "i")
  title(ylab = ylab)
  for (i in seq_len(k)) {
    v <- ind$Pheno[ind$hap == gs$hap[i]]
    d <- D[[i]]
    w <- 0.40 * d$y / max(d$y)
    polygon(c(i - w, rev(i + w)), c(d$x, rev(d$x)),
            col = HAP_COL[(i - 1) %% length(HAP_COL) + 1], border = "black", lwd = LWD)
    qs <- quantile(v, c(0.25, 0.5, 0.75)); iqr <- qs[3] - qs[1]
    lo <- min(v[v >= qs[1] - 1.5 * iqr]); hi <- max(v[v <= qs[3] + 1.5 * iqr])
    segments(i, lo, i, hi, lwd = LWD)
    rect(i - 0.09, qs[1], i + 0.09, qs[3], col = "white", border = "black", lwd = LWD)
    segments(i - 0.09, qs[2], i + 0.09, qs[2], lwd = LWD * 1.6)
  }
  # brackets: every group against the largest, Holm-corrected within the gene (step 07)
  if (nb) {
    u_in <- ru / pin_h                                 # data units per inch
    tick <- 0.012 * diff(ylim)
    rk   <- rank(p$br$y.position, ties.method = "first")
    yb   <- dhi + pad + tick + (gap_in + (rk - 1) * step_in) * u_in
    xr <- match(p$br$group1, gs$hap); xo <- match(p$br$group2, gs$hap)
    segments(xr, yb, xo, yb, lwd = LWD)
    segments(c(xr, xo), rep(yb, 2) - tick, c(xr, xo), rep(yb, 2), lwd = LWD)
    text((xr + xo) / 2, yb + 0.02 * u_in, p$br$p.adj.signif, adj = c(0.5, 0), cex = cex_ax)
  }
  # y ticks only over the data, not up into the bracket zone
  yt <- pretty(c(dlo, dhi)); yt <- yt[yt >= ylim[1] & yt <= dhi]
  axis(2, at = yt, lwd = LWD, cex.axis = cex_ax)
  # gap.axis = -1: draw every label; at 9 pt R otherwise dropped touching "(n=..)" labels (GH17)
  axis(1, at = seq_len(k), labels = gs$hap, tick = FALSE, line = -0.55, cex.axis = cex_ax, gap.axis = -1)
  axis(1, at = seq_len(k), labels = paste0("(n=", gs$n, ")"), tick = FALSE, line = 0.15,
       cex.axis = cex_ax, gap.axis = -1)
  box(bty = "l", lwd = LWD)
}

# Barcode as in Fig. 3, plus one extra row at the bottom carrying the significance triangles.
TRI_ROWS <- 0.8                              # height of the triangle row, in barcode rows
barcode_panel <- function(cons, el, cex_ax = 1) {
  n_c <- ncol(cons); rows <- nrow(cons) + 1L + nrow(el)
  mat  <- rbind(cons, rep(NA_real_, n_c), el)
  labs <- c(paste("Group", rownames(cons)), "", rownames(el))
  blank <- nrow(cons) + 1L
  plot(NA, xlim = c(0.5, n_c + 0.5), ylim = c(rows + 0.5 + TRI_ROWS, 0.5), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  for (r in seq_len(rows)) {
    if (r == blank) next
    for (cc in seq_len(n_c)) {
      v <- mat[r, cc]
      rect(cc - 0.5, r - TILE_H / 2, cc + 0.5, r + TILE_H / 2,
           col = if (is.na(v)) COL_MISS else GT_COL[[as.character(v)]], border = "white", lwd = 0.25)
    }
  }
  points(sig_col, rep(rows + 0.5 + TRI_ROWS / 2, length(sig_col)), pch = SIG_PCH, col = SIG_COL, cex = 0.9)
  axis(2, at = seq_len(rows)[-blank], labels = labs[-blank], tick = FALSE,
       line = -0.15, las = 1, cex.axis = cex_ax)
}

panel_label <- function(p, letter, name = TRUE, stats = TRUE, dy_in = 0.06, dx_in = 0, cex = 1) {
  ytxt <- grconvertY(grconvertY(1, "npc", "inches") + dy_in, "inches", "user")
  xl   <- grconvertX(grconvertX(0, "nfc", "inches") + dx_in, "inches", "user")
  lab  <- if (name) bquote(bold(.(letter)) ~~ .(p$label) ~ "(" * .(p$trait_lab) * ")") else bquote(bold(.(letter)))
  text(xl, ytxt, lab, adj = c(0, 0), xpd = NA, cex = cex)
  if (stats) text(grconvertX(grconvertX(1, "nfc", "inches") - 0.10, "inches", "user"),
                  ytxt, p_expr(p$p, p$eta2), adj = c(1, 0), xpd = NA, cex = cex)
}

# Legend under the barcode column, drawn with Fig. 3's hlegend() (per-item widths): genotype
# states on the first line, the triangle key on the second — one line does not fit 104 mm at 9 pt.
hlegend <- function(labs, fill = rep(NA, length(labs)), pch = rep(NA, length(labs)),
                    pcol = rep(NA, length(labs)), y = 0.5, cex = 1) {
  ch   <- par("cin")[2] * par("cex") * cex            # character height, inches
  key  <- 0.75 * ch; gap <- 0.35 * ch; sep <- 1.1 * ch
  tw   <- strwidth(labs, units = "inches", cex = cex)
  wid  <- key + gap + tw
  pin  <- par("pin"); x_in <- (pin[1] - (sum(wid) + sep * (length(labs) - 1))) / 2
  ux <- function(i) i / pin[1]; uy <- function(i) i / pin[2]
  for (j in seq_along(labs)) {
    if (!is.na(fill[j])) rect(ux(x_in), y - uy(key / 2), ux(x_in + key), y + uy(key / 2),
                              col = fill[j], border = "black", lwd = LWD, xpd = NA)
    if (!is.na(pch[j]))  points(ux(x_in + key / 2), y, pch = pch[j], col = pcol[j], cex = 0.9, xpd = NA)
    text(ux(x_in + key + gap), y, labs[j], adj = c(0, 0.5), cex = cex, xpd = NA)
    x_in <- x_in + wid[j] + sep
  }
}

# Legend drawn INSIDE the barcode column, directly under the barcode (2026-09-24, user request:
# the separate legend row took height the figure did not need). Same keys and sizes as hlegend(),
# but positioned in device inches: xc_in = centre of the row, y_in = its vertical centre.
hlegend_dev <- function(labs, xc_in, y_in, fill = rep(NA, length(labs)), pch = rep(NA, length(labs)),
                        pcol = rep(NA, length(labs)), cex = 1) {
  ch  <- par("cin")[2] * par("cex") * cex
  key <- 0.75 * ch; gap <- 0.35 * ch; sep <- 1.1 * ch
  wid <- key + gap + strwidth(labs, units = "inches", cex = cex)
  x   <- xc_in - (sum(wid) + sep * (length(labs) - 1)) / 2
  ux <- function(i) grconvertX(i, "inches", "user"); uy <- function(i) grconvertY(i, "inches", "user")
  for (j in seq_along(labs)) {
    if (!is.na(fill[j])) rect(ux(x), uy(y_in - key / 2), ux(x + key), uy(y_in + key / 2),
                              col = fill[j], border = "black", lwd = LWD, xpd = NA)
    if (!is.na(pch[j]))  points(ux(x + key / 2), uy(y_in), pch = pch[j], col = pcol[j], cex = 0.9, xpd = NA)
    text(ux(x + key + gap), uy(y_in), labs[j], adj = c(0, 0.5), cex = cex, xpd = NA)
    x <- x + wid[j] + sep
  }
}

ylab_of <- function(p) paste(if (p$trait == "fiber") "Fiber" else "Starch", "BLUP")

# ── Layout: violins a/b stacked left, one barcode c spanning both rows right ──
ROW_IN <- 2.35                               # as Fig. 3
LEG_BLOCK_IN <- 0.52                         # legend (two lines) under the barcode, inside column c
BAR_ROW_IN <- 0.175                          # barcode row pitch, as Fig. 3

draw_fig4 <- function() {
  M_TOP <- 1.5; M_BOT <- 2.6
  layout(rbind(c(1, 3), c(2, 3)), widths = c(70, 104), heights = c(ROW_IN, ROW_IN))
  par(cex = 1)                                               # undo layout()'s automatic 0.66 shrink
  par(family = FONT, las = 1, mgp = c(2.0, 0.45, 0), tcl = -0.25,
      cex.axis = 1, cex.lab = 1, lwd = LWD, xpd = FALSE)
  for (i in seq_along(P)) {
    p <- P[[i]]
    par(mar = c(M_BOT, 3.8, M_TOP, 0.6))
    violin_panel(p, ylab_of(p))
    panel_label(p, letters[i])                               # a fiber, b starch
  }
  # barcode + legend as one block, centred vertically in column c
  nr    <- nrow(cons) + 1L + nrow(EL) + TRI_ROWS
  bar_h <- nr * BAR_ROW_IN
  free  <- 2 * ROW_IN - bar_h - LEG_BLOCK_IN
  par(mai = c(free / 2 + LEG_BLOCK_IN, 5.2 * LINE_IN, free / 2, 0.6 * LINE_IN))
  barcode_panel(cons, EL)
  panel_label(P$fiber, "c", name = FALSE, stats = FALSE, dx_in = 0.30)
  xc  <- mm2in(70 + 104 / 2)                                 # centre of the barcode column
  bot <- free / 2 + LEG_BLOCK_IN                             # bottom of the barcode plot region
  keys <- as.character(drawn); labs <- unname(GT_LAB[keys]); cols <- unname(GT_COL[keys])
  if (has_missing) { labs <- c(labs, "Missing"); cols <- c(cols, COL_MISS) }
  hlegend_dev(labs, xc, bot - 0.17, fill = cols)
  hlegend_dev("Genome-wide significant SNP", xc, bot - 0.39, pch = SIG_PCH, pcol = SIG_COL)
}

render <- function(name, fun, h_in) {
  stopifnot(h_in * 25.4 <= H_MAX)
  png(file.path(OUT, paste0(name, ".png")), width = mm2in(W_MM), height = h_in,
      units = "in", res = DPI, pointsize = PT, type = "cairo", family = FONT, bg = "white")
  fun(); invisible(dev.off())
  stopifnot(system2("python3", c("-c", shQuote(sprintf(
    "from PIL import Image; Image.open('%s').convert('RGB').save('%s', compression='tiff_lzw', dpi=(%d,%d))",
    file.path(OUT, paste0(name, ".png")), file.path(OUT, paste0(name, ".tif")), DPI, DPI)))) == 0)
  cat(sprintf("[fig4] OK -> %s.{png,tif}  (%.0f x %.1f mm, %d dpi)\n", name, W_MM, h_in * 25.4, DPI))
}

render("Fig4", draw_fig4, 2 * ROW_IN)
