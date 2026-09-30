#!/usr/bin/env Rscript
# make_figure_3.R — Fig. 3 of the TAG manuscript (Results ch. 3): haplotype structure of the
# three candidate genes carried forward, with the elite cultivars drawn as aligned barcodes.
#
# Panels: a GPAT6 (fiber), b GH17 (fiber), c PHT4;3 (starch) — the three tier-1 genes of
# 06_USED_genes_selected_to_present, in the order they are discussed in the text.
#
# Each panel has two parts:
#   * violins of the trait BLUP by wild haplotype group (boxplot inside, group n below),
#     with the Wilcoxon-vs-largest-group brackets of step 07 (Holm-corrected within the gene)
#     and a q + eta-squared label from step 04;
#   * aligned genotype barcodes: one row per wild haplotype group (per-SNP majority consensus),
#     a white gap, then one row per elite cultivar. Columns are SNPs in genomic order.
#
# Only the `shared_sites` version is drawn (decided in build/BLUEPRINT.md): a column appears
# only where the wild and the elite call sets both have a record with identical REF and ALT,
# so no genotype is assumed anywhere in the figure. Elite lines are NOT assigned to a haplotype
# group — crosshap never saw them; the barcodes are aligned and the comparison left to the reader.
#
# LAYOUT — side by side, chosen by the user 2026-09-23: one row per gene, violins left
# (a, b, c) and the aligned barcodes right (d, e, f). A stacked alternative (violins above
# the barcodes, as in the mini paper's Fig. 4) was drawn on 2026-09-22 and rejected: at
# 174 mm it reached 226 mm of the 234 mm height limit, which left the violins cramped.
# `draw_stacked()` is kept so it can be regenerated, but it is no longer rendered.
#
# Changes requested on the first side-by-side draft (2026-09-23):
#   1. violins are drawn untrimmed, so no violin is cut off at the extreme observations;
#   2. lettering raised from 9 pt to 11 pt (TAG allows 8-12) — the 9 pt draft was unreadable
#      at print size;
#   3. the genotype legend moved from the centre of the figure to under the barcode column;
#   4. the barcode panels carry their own panel letters d, e, f.
#
# Fixed 2026-09-24 (user decision) — TRUE TEXT SIZE. layout() with three or more rows silently
# sets par(cex = 0.66), and nothing reset it, so every label was drawn at 0.66 x PT: the version
# approved 2026-09-23 as "11 pt" measured ~7.3 pt on the page, below TAG's 8 pt minimum. It also
# made the margin arithmetic (LINE_IN, in lines of PT) wrong. par(cex = 1) is now set after each
# layout() call, and PT is 9: a true 9 pt, inside TAG's 8-12 pt, ~25% larger than the approved look.
# Same date: the barcode consensus now resolves 50/50 ties to REF (step 07, user decision), which
# turns two GH17 cells in panel e from missing to REF.
#
# TAG figure spec (10_USED_Paper_writing/TAG_requirements.md): 174 mm wide (full page width),
# height <= 234 mm, Arial-metric sans (Liberation Sans) at a true 9 pt at final size (TAG range 8-12),
# lines >= 0.3 pt, RGB, 600 dpi (combination art). No title inside the image: the gene names are
# panel labels, and the windows, SNP counts and colour definitions belong to the caption.
#
# Inputs (read-only):
#   07_USED_elite_lines_compariosn_to_wild_lines/
#     intermediates/matrices/<gene_id>__shared_sites.rds   consensus, elite, group_stats, indfile
#     results/tables/Table_pairwise_group_tests.tsv        bracket statistics
#     results/tables/Table_allele_concordance.tsv          REF/ALT check, re-asserted here
#     config/elite_lines.tsv                               the five cultivars, in display order
#   04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv
#                                                          q, eta2, n_groups cross-check
# Output: Figure_3/Fig3.{png,tif}  (png for the docx build, tif for submission)
#
# Created 2026-09-22. Run: Rscript make_figure_3.R   (seconds)
#
# RENUMBERED 2026-09-30 (S. Hübner's comments on the figures; user decision): old Fig. 2 was dissolved,
# so this figure is now Fig. 3 (it was Fig. 4) and this script was make_figure_4.R, writing
# Figure_4/Fig4.*. Only the output folder, file names, log tags and this header changed; the
# re-rendered image is pixel-identical to the approved Fig4 (checked 2026-09-30). Hand-over:
# 10_USED_Paper_writing/new_publishing_paper/build/STRUCTURE_CHANGES_2026-10.md

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))

ROOT   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
STEP07 <- file.path(ROOT, "07_USED_elite_lines_compariosn_to_wild_lines")
STEP04 <- file.path(ROOT, "04_USED_haplotype_analysis_crosshap", "04_runs", "loci_LDspan_eps06_V4")
OUT    <- file.path(ROOT, "08_USED_creating_figures", "Figure_3")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

VERSION <- "shared_sites"

# Panel order = the order the genes are discussed in Results ch. 3.
GENES <- list(
  list(id = "HORVU.MOREX.r3.3HG0301300", short = "GPAT6",  trait = "fiber",
       label = "GPAT6",  trait_lab = "fiber"),
  list(id = "HORVU.MOREX.r3.5HG0487060", short = "GH17",   trait = "fiber",
       label = "GH17",   trait_lab = "fiber"),
  list(id = "HORVU.MOREX.r3.3HG0301710", short = "PHT4-3", trait = "starch",
       label = "PHT4;3", trait_lab = "starch")
)

# ── Style ─────────────────────────────────────────────────────────────────────
FONT   <- "Liberation Sans"                  # metric-compatible with Arial
PT     <- 9                                  # all lettering, TRUE size on the page (TAG 8-12); needs par(cex = 1) after layout()
DPI    <- 600
W_MM   <- 174                                # TAG full width
H_MAX  <- 234                                # TAG maximum height
LWD    <- 0.75                               # >= 0.3 pt
LINE_IN <- 1.2 * PT / 72                     # one margin line, inches

# Genotype codes, from 07_.../scripts/03_build_matrices.R
CODE_REF <- 0; CODE_ALT <- 1; CODE_HET <- 2; CODE_NORECORD <- 3; CODE_TRI <- 4
# Barcode colours, from 07_.../config/params.sh (themselves inherited from step 04's heatmaps)
GT_COL <- c("0" = "#FFFACD", "1" = "#2F4F4F", "2" = "#C46210",
            "3" = "#FFFFFF", "4" = "#7B3FA0")
GT_LAB <- c("0" = "Reference", "1" = "Alternate", "2" = "Heterozygous",
            "3" = "No elite record", "4" = "Third allele")
COL_MISS <- "grey70"
TILE_H   <- 0.72                             # white gap between barcode rows, as in step 07

# Haplotype-group fills: the CVD-validated categorical palette used for the Manhattan plots (Fig. 2).
HAP_COL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")

# ── Data ──────────────────────────────────────────────────────────────────────
PW    <- fread(file.path(STEP07, "results", "tables", "Table_pairwise_group_tests.tsv"))
CONC  <- fread(file.path(STEP07, "results", "tables", "Table_allele_concordance.tsv"))
G04   <- fread(file.path(STEP04, "Stats", "gene_results.tsv"))
ELITE <- fread(file.path(STEP07, "config", "elite_lines.tsv"), skip = "line_name")

# The whole figure rests on REF and ALT meaning the same base in both call sets; if a shared
# position were swapped, every barcode below would be silently inverted. Step 07 aborts on that,
# and it is re-asserted here so this script cannot draw an inverted figure from a stale matrix.
stopifnot(!any(CONC$status %in% c("swapped", "ref_differs")))

P <- list()
for (g in GENES) {
  m <- readRDS(file.path(STEP07, "intermediates", "matrices",
                         paste0(g$id, "__", VERSION, ".rds")))
  stopifnot(m$version == VERSION, m$gene_id == g$id, m$trait == g$trait)

  s04 <- G04[gene_id == g$id & trait == g$trait]
  stopifnot(nrow(s04) == 1L,
            s04$n_groups == nrow(m$group_stats),          # same grouping as step 04
            identical(as.integer(m$group_stats$n), as.integer(strsplit(s04$group_sizes, "|", fixed = TRUE)[[1]])),
            nrow(m$indfile) == sum(m$group_stats$n))      # every grouped accession is drawn
  # group means recomputed from the individuals, not taken on trust
  gm <- tapply(m$indfile$Pheno, m$indfile$hap, mean)
  stopifnot(max(abs(gm[m$group_stats$hap] - m$group_stats$mean_pheno)) < 1e-6)

  stopifnot(identical(rownames(m$elite), ELITE$line_name),  # five cultivars, display order
            ncol(m$consensus) == m$n_shared, ncol(m$elite) == m$n_shared,
            !any(m$consensus %in% c(CODE_HET, CODE_NORECORD, CODE_TRI)))  # wild: REF/ALT/NA only

  br <- PW[gene_id == g$id][order(match(group2, m$group_stats$hap))]
  P[[g$short]] <- c(g, list(m = m, br = br,
                            q = s04$fdr_q, eta2 = s04$eta_squared,
                            ref_grp = m$group_stats$hap[which.max(m$group_stats$n)]))
  cat(sprintf("[fig3] %-7s %2d groups | %3d shared SNPs | q %.2e | eta2 %.3f | %d elite rows\n",
              g$short, nrow(m$group_stats), m$n_shared, s04$fdr_q, s04$eta_squared, nrow(m$elite)))
}

# States actually drawn anywhere in the figure — the legend lists only these (step 07 convention).
drawn <- sort(unique(unlist(lapply(P, function(p) c(as.vector(p$m$consensus), as.vector(p$m$elite))))))
drawn <- drawn[!is.na(drawn)]
has_missing <- any(vapply(P, function(p) anyNA(p$m$consensus) || anyNA(p$m$elite), logical(1)))

# ── Drawing helpers ───────────────────────────────────────────────────────────
mm2in <- function(mm) mm / 25.4

# Scientific q in the journal's style: 5.3 x 10^-6
q_expr <- function(q, eta2) {
  e <- floor(log10(q)); mant <- q / 10^e
  bquote(italic(q) == .(sprintf("%.1f", mant)) %*% 10^.(e) * "," ~~ italic(eta)^2 == .(sprintf("%.3f", eta2)))
}

violin_panel <- function(p, ylab, cex_ax = 1) {
  gs  <- p$m$group_stats
  ind <- p$m$indfile
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

barcode_panel <- function(p, cex_ax = 1) {
  cons <- p$m$consensus; el <- p$m$elite
  n_c <- ncol(cons); rows <- nrow(cons) + 1L + nrow(el)     # + one blank separator row
  mat  <- rbind(cons, rep(NA_real_, n_c), el)
  labs <- c(paste("Group", rownames(cons)), "", rownames(el))
  blank <- nrow(cons) + 1L

  plot(NA, xlim = c(0.5, n_c + 0.5), ylim = c(rows + 0.5, 0.5), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  for (r in seq_len(rows)) {
    if (r == blank) next
    y <- r
    for (cc in seq_len(n_c)) {
      v <- mat[r, cc]
      col <- if (is.na(v)) COL_MISS else GT_COL[[as.character(v)]]
      rect(cc - 0.5, y - TILE_H / 2, cc + 0.5, y + TILE_H / 2,
           col = col, border = "white", lwd = 0.25)
    }
  }
  axis(2, at = seq_len(rows)[-blank], labels = labs[-blank], tick = FALSE,
       line = -0.15, las = 1, cex.axis = cex_ax)
}

# `name = TRUE` prints the gene and trait after the letter; the barcode panel of the same
# row carries the letter alone, since the row header already names the gene.
panel_label <- function(p, letter, name = TRUE, stats = TRUE,
                        dy_in = 0.06, dx_in = 0, cex = 1) {
  ytxt <- grconvertY(grconvertY(1, "npc", "inches") + dy_in, "inches", "user")
  xl   <- grconvertX(grconvertX(0, "nfc", "inches") + dx_in, "inches", "user")
  lab  <- if (name) bquote(bold(.(letter)) ~~ .(p$label) ~ "(" * .(p$trait_lab) * ")")
          else      bquote(bold(.(letter)))
  text(xl, ytxt, lab, adj = c(0, 0), xpd = NA, cex = cex)
  if (stats) text(grconvertX(grconvertX(1, "nfc", "inches") - 0.10, "inches", "user"),
                  ytxt, q_expr(p$q, p$eta2), adj = c(1, 0), xpd = NA, cex = cex)
}

# The legend sits under the barcode column rather than under the whole figure, so it reads
# as belonging to the barcodes and not to the violins (requested 2026-09-23).
# Drawn by hand since 2026-09-24: R 4.1's legend() gives every item the width of the longest
# label, and at a true 9 pt that no longer fit the column. Each key now takes its own width, and
# the row is centred on the whole barcode column (labels + tiles).
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

draw_legend <- function(cex = 1, pos = "center", mar = c(0, 0, 0, 0)) {
  par(mar = mar); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  keys <- as.character(drawn); labs <- unname(GT_LAB[keys]); cols <- unname(GT_COL[keys])
  if (has_missing) { labs <- c(labs, "Missing"); cols <- c(cols, COL_MISS) }
  hlegend(labs, fill = cols, cex = cex)
}

ylab_of <- function(p) paste(if (p$trait == "fiber") "Fiber" else "Starch", "BLUP")

# ── Layout 1: stacked — one full-width panel per gene ─────────────────────────
draw_stacked <- function() {
  # per gene: violin (plot 1.00 in) over barcode (0.115 in per row)
  bar_in <- vapply(P, function(p) (nrow(p$m$consensus) + 1L + nrow(p$m$elite)) * 0.115, numeric(1))
  VIO_IN <- 1.00
  M_TOP <- 1.5; M_VIO_BOT <- 2.1; M_BAR_TOP <- 0.3; M_BAR_BOT <- 1.0
  M_LEFT <- 5.0; M_RIGHT <- 0.6
  h <- as.vector(rbind(VIO_IN + (M_TOP + M_VIO_BOT) * LINE_IN,
                       bar_in + (M_BAR_TOP + M_BAR_BOT) * LINE_IN))
  layout(matrix(1:7, ncol = 1), heights = c(h, 0.24))
  par(cex = 1)                                         # undo layout()'s automatic 0.66 shrink
  par(family = FONT, las = 1, mgp = c(1.7, 0.4, 0), tcl = -0.22,
      cex.axis = 1, cex.lab = 1, lwd = LWD, xpd = FALSE)
  for (i in seq_along(P)) {
    p <- P[[i]]
    par(mar = c(M_VIO_BOT, M_LEFT, M_TOP, M_RIGHT))
    violin_panel(p, ylab_of(p))
    panel_label(p, letters[i])
    par(mar = c(M_BAR_BOT, M_LEFT, M_BAR_TOP, M_RIGHT))
    barcode_panel(p)
  }
  draw_legend()
}

# ── Layout 2: side by side — violins left, barcodes right ─────────────────────
ROW_IN_SBS <- 2.35     # row height, inches — also used to size the canvas
LEG_IN_SBS <- 0.34

draw_sidebyside <- function() {
  M_TOP <- 1.5; M_BOT <- 2.6
  # panel 7 is left empty so the legend (panel 8) sits under the barcode column only
  layout(rbind(matrix(1:6, ncol = 2, byrow = TRUE), c(7, 8)),
         widths = c(70, 104), heights = c(rep(ROW_IN_SBS, 3), LEG_IN_SBS))
  par(cex = 1)                                         # undo layout()'s automatic 0.66 shrink
  par(family = FONT, las = 1, mgp = c(2.0, 0.45, 0), tcl = -0.25,
      cex.axis = 1, cex.lab = 1, lwd = LWD, xpd = FALSE)
  for (i in seq_along(P)) {
    p <- P[[i]]
    par(mar = c(M_BOT, 3.8, M_TOP, 0.6))
    violin_panel(p, ylab_of(p))
    panel_label(p, letters[i])                       # a, b, c — violins
    # barcode centred vertically in the row, so it does not stretch with the number of groups
    nr  <- nrow(p$m$consensus) + 1L + nrow(p$m$elite)
    pad <- (ROW_IN_SBS - (M_TOP + M_BOT) * LINE_IN - nr * 0.175) / 2 / LINE_IN
    par(mar = c(M_BOT + pad, 5.2, M_TOP + pad, 0.6))
    barcode_panel(p)
    panel_label(p, letters[i + 3L], name = FALSE, stats = FALSE, dx_in = 0.30)  # d, e, f
  }
  par(mar = c(0, 0, 0, 0)); plot.new()               # spacer under the violin column
  draw_legend(mar = c(0, 0, 0, 0))                    # centred on the whole barcode column
}

# ── Render ────────────────────────────────────────────────────────────────────
render <- function(name, fun, h_in) {
  stopifnot(h_in * 25.4 <= H_MAX)
  png(file.path(OUT, paste0(name, ".png")), width = mm2in(W_MM), height = h_in,
      units = "in", res = DPI, pointsize = PT, type = "cairo", family = FONT, bg = "white")
  fun(); invisible(dev.off())
  stopifnot(system2("python3", c("-c", shQuote(sprintf(
    "from PIL import Image; Image.open('%s').convert('RGB').save('%s', compression='tiff_lzw', dpi=(%d,%d))",
    file.path(OUT, paste0(name, ".png")), file.path(OUT, paste0(name, ".tif")), DPI, DPI)))) == 0)
  cat(sprintf("[fig3] OK -> %s.{png,tif}  (%.0f x %.1f mm, %d dpi)\n", name, W_MM, h_in * 25.4, DPI))
}

render("Fig3", draw_sidebyside, 3 * ROW_IN_SBS + LEG_IN_SBS)
