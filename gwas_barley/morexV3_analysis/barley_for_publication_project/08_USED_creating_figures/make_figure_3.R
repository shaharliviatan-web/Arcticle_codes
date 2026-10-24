#!/usr/bin/env Rscript
# make_figure_3.R — Fig. 3 of the TAG manuscript (Results ch. 3): haplotype structure of the
# four candidate genes, with the elite cultivars drawn as aligned barcodes.
#
# COMBINED 2026-10-17 (S. Hübner's recommendation; user decision): Results ch. 3 and ch. 4 are merged
# into one chapter, and old Fig. 3 (three genes) and old Fig. 4 (the 7H GDSL gene,
# make_figure_4.R) into this one figure. Both scripts already shared the drawing code; this one
# keeps Fig. 3's code and data checks and takes over Fig. 4's data loading and checks unchanged.
# Nothing is recomputed: every value is read from the same step-04, step-07 and 7H-branch outputs.
# APPROVED 2026-10-17 by the user (drawn first as Figure_3/Fig3_combined_DRAFT.*, then renamed): it
# replaces the former Fig3.{png,tif}; the former Figure_4/Fig4.* and make_figure_4.R were deleted
# the same day (user request; the script and Fig4.tif remain in git).
#
# Panels, in the order of the text (pipeline order, user decision 2026-10-17):
#   violins  a GPAT6 (fiber)   b GH17 (fiber)   c PHT4;3 (starch)   d GDSL (fiber)   e GDSL (starch)
#   barcodes f GPAT6           g GH17           h PHT4;3            i GDSL
# One row per gene, violins left and the aligned barcodes right. The GDSL row holds two half-width
# violins, fiber and starch side by side, and one barcode: the two traits give identical
# haplotype groups (the same 212 | 34 accessions; asserted below), so their barcodes are identical.
#
# Each violin panel: trait BLUP by wild haplotype group (boxplot inside, group n below), the
# Wilcoxon-vs-largest-group brackets (Holm-corrected within the gene), and above the panel the
# gene-level statistic exactly as in the two former figures: the Benjamini-Hochberg q (step 04)
# for GPAT6, GH17 and PHT4;3; the raw Kruskal-Wallis P for the GDSL gene (a single pre-specified
# test outside the step-04 BH family, as in the former Fig. 4).
# Each barcode panel: one row per wild haplotype group (per-SNP majority consensus, ties -> REF),
# a white gap, then one row per elite cultivar; columns are SNPs in genomic order. The GDSL
# barcode keeps the former Fig. 4's red triangles under the three genome-wide significant SNPs
# (fiber lead 7H:573,606,306; starch lead 7H:573,606,460; starch 7H:573,606,491), with their key
# directly under that barcode; the genotype legend sits under the barcode column.
#
# Only the `shared_sites` version is drawn (decided in build/BLUEPRINT.md): a column appears
# only where the wild and the elite call sets both have a record with identical REF and ALT,
# so no genotype is assumed anywhere in the figure. Elite lines are NOT assigned to a haplotype
# group — crosshap never saw them; the barcodes are aligned and the comparison left to the reader.
#
# LAYOUT (2026-10-17). TAG maximum: 174 mm x 234 mm; the figure uses the full height (rows of
# ~56 mm; the former Fig. 3 had ~60 mm). Text stays at a true 9 pt. Column split 86 | 88 mm (was
# 70 | 104): each GDSL half panel must hold its two-line header (the P label is 39.3 mm at 9 pt),
# and the genotype legend (85.4 mm at 9 pt) must fit under the barcode column. The GDSL panels
# print the gene and trait on one line and the statistic on the line below; the other panels keep
# both on one line.
#
# History of the former Fig. 3 (kept; the first two items still describe this code):
#   2026-09-23 side-by-side layout chosen by the user (a stacked alternative was rejected: at
#     174 mm it reached 226 mm of the 234 mm limit and left the violins cramped; its draw_stacked()
#     was removed on 2026-10-17 — it cannot draw four genes, and it is in git);
#   2026-09-23 violins untrimmed; legend under the barcode column; barcode panels carry letters;
#   2026-09-24 TRUE TEXT SIZE: layout() with three or more rows silently sets par(cex = 0.66);
#     par(cex = 1) is set after each layout() call, and PT is 9 (a true 9 pt, TAG range 8-12);
#     the barcode consensus resolves 50/50 ties to REF (step 07, user decision).
#   2026-09-30 renumbered from Fig. 4 to Fig. 3 (S. Hübner's comments; pixel-identical).
#
# TAG figure spec (10_USED_Paper_writing/TAG_requirements.md): 174 mm wide (full page width),
# height <= 234 mm, Arial-metric sans (Liberation Sans) at a true 9 pt at final size (TAG range 8-12),
# lines >= 0.3 pt, RGB, 600 dpi (combination art). No title inside the image: the gene names are
# panel labels, and the windows, SNP counts and colour definitions belong to the caption.
#
# Inputs (read-only):
#   07_USED_elite_lines_compariosn_to_wild_lines/                  (GPAT6, GH17, PHT4;3)
#     intermediates/matrices/<gene_id>__shared_sites.rds   consensus, elite, group_stats, indfile
#     results/tables/Table_pairwise_group_tests.tsv        bracket statistics
#     results/tables/Table_allele_concordance.tsv          REF/ALT check, re-asserted here
#     config/elite_lines.tsv                               the five cultivars, in display order (all four genes)
#   04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv
#                                                          q, eta2, n_groups cross-check
#   03_01_7H_branch_Starch_Fiber_shared_signal_explore/ = B       (GDSL gene; crosshap MGmin 3, eps 0.9)
#     B/Cache/<trait>/<gene>/MGmin_3/eps_0.9/HapObject.rds    groups, phenotypes, marker groups
#     B/results/tables/mgmin3_gene_results.tsv                KW P, eta2, group sizes (cross-check)
#     B/results/tables/pairwise_group_tests.tsv               bracket statistics (Wilcoxon, Holm)
#     B/results/tables/Table_site_overlap.tsv                 shared-site count (cross-check)
#     B/results/tables/Table_allele_concordance.tsv           REF/ALT check, re-asserted here
#     B/results/tables/Table_elite_genotypes_wide__<gene>.tsv elite states (cross-check)
#     B/work/raw/<trait>/<gene>.vcf.gz                        wild raw genotypes, gene +/- 1 kb
#     B/work/elite/7HG0729030.elite.vcf.gz                    elite genotypes (DivBrowse, MorexV3)
#     B/archive_v1_MGmin2_eps0.6/results/tables/signal_snps.tsv  genome-wide significance per SNP
# Output: Figure_3/Fig3.{png,tif}  (png for the docx build, tif for submission)
#
# Created 2026-09-22 (as make_figure_4.R; renumbered 2026-09-30); combined 2026-10-17.
# Run: Rscript make_figure_3.R   (seconds)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(vcfR) })

ROOT   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
STEP07 <- file.path(ROOT, "07_USED_elite_lines_compariosn_to_wild_lines")
STEP04 <- file.path(ROOT, "04_USED_haplotype_analysis_crosshap", "04_runs", "loci_LDspan_eps06_V4")
BRANCH <- file.path(ROOT, "03_01_7H_branch_Starch_Fiber_shared_signal_explore")   # the GDSL analysis
ARCH   <- file.path(BRANCH, "archive_v1_MGmin2_eps0.6")   # its archived first version (MGmin-independent tables)
OUT    <- file.path(ROOT, "08_USED_creating_figures", "Figure_3")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
OUT_NAME <- "Fig3"                                        # the draft "Fig3_combined_DRAFT" was approved 2026-10-17

VERSION <- "shared_sites"

# Rows 1-3 = the three genes carried forward from the step-04 run, in the order of the text.
GENES <- list(
  list(id = "HORVU.MOREX.r3.3HG0301300", short = "GPAT6",  trait = "fiber",
       label = "GPAT6",  trait_lab = "fiber"),
  list(id = "HORVU.MOREX.r3.5HG0487060", short = "GH17",   trait = "fiber",
       label = "GH17",   trait_lab = "fiber"),
  list(id = "HORVU.MOREX.r3.3HG0301710", short = "PHT4-3", trait = "starch",
       label = "PHT4;3", trait_lab = "starch")
)
# Row 4 = the GDSL gene at the shared 7H signal, both traits.
GDSL_ID <- "HORVU.MOREX.r3.7HG0729030"
MGMIN <- 3; EPS <- 0.9
GDSL_TRAITS <- list(
  list(trait = "fiber",  label = "GDSL", trait_lab = "fiber"),
  list(trait = "starch", label = "GDSL", trait_lab = "starch")
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
SIG_COL  <- "red"                            # significance triangles (former Fig. 4)
SIG_PCH  <- 17
TRI_ROWS <- 0.8                              # height of the triangle row, in barcode rows

# Haplotype-group fills: the CVD-validated categorical palette used for the Manhattan plots (Fig. 2).
HAP_COL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4")

# ── Data: GPAT6, GH17, PHT4;3 (steps 04 and 07) ───────────────────────────────
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
  P[[g$short]] <- c(g, list(m = m, br = br, stat = "q", sval = s04$fdr_q, eta2 = s04$eta_squared,
                            sig_col = NULL,
                            ref_grp = m$group_stats$hap[which.max(m$group_stats$n)]))
  cat(sprintf("[fig3] %-7s %2d groups | %3d shared SNPs | q %.2e | eta2 %.3f | %d elite rows\n",
              g$short, nrow(m$group_stats), m$n_shared, s04$fdr_q, s04$eta_squared, nrow(m$elite)))
}

# ── Data: the GDSL gene (7H branch; loading and checks as in the former make_figure_4.R) ──
tab7 <- function(f) fread(file.path(BRANCH, "results", "tables", f))
GR7    <- tab7("mgmin3_gene_results.tsv")[MGmin == MGMIN & epsilon == EPS]
PW7    <- tab7("pairwise_group_tests.tsv")[eps == EPS]
OVL7   <- tab7("Table_site_overlap.tsv")
CONC7  <- tab7("Table_allele_concordance.tsv")
EWIDE7 <- tab7(paste0("Table_elite_genotypes_wide__", GDSL_ID, ".tsv"))
SIGT7  <- fread(file.path(ARCH, "results", "tables", "signal_snps.tsv"))

stopifnot(!any(CONC7$status %in% c("swapped", "ref_differs")))   # REF/ALT mean the same base

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

W7 <- read_gt(file.path(BRANCH, "work", "raw", "fiber", paste0(GDSL_ID, ".vcf.gz")))
stopifnot(identical(W7, read_gt(file.path(BRANCH, "work", "raw", "starch", paste0(GDSL_ID, ".vcf.gz")))))
E7 <- read_gt(file.path(BRANCH, "work", "elite", "7HG0729030.elite.vcf.gz"))
stopifnot(all(ELITE$sample_id %in% rownames(E7)))
E7 <- E7[ELITE$sample_id, , drop = FALSE]; rownames(E7) <- ELITE$line_name

shared7 <- intersect(colnames(W7), colnames(E7))
spos7   <- as.integer(sub("^[^:]+:([0-9]+):.*", "\\1", shared7))
shared7 <- shared7[order(spos7)]; spos7 <- sort(spos7)
stopifnot(length(shared7) == OVL7$n_shared)
EL7 <- E7[, shared7, drop = FALSE]

# elite states re-derived from the VCF must equal the published table
st <- function(x) unname(ifelse(is.na(x), "Missing", GT_LAB[as.character(x)]))
stopifnot(identical(EWIDE7$site_key, shared7),
          all(vapply(ELITE$line_name, function(l) identical(st(EL7[l, ]), EWIDE7[[l]]), logical(1))))

# Haplotype groups per trait, and the check that justifies ONE barcode for both traits
G7 <- list(); groups7 <- list()
for (tt in GDSL_TRAITS) {
  ho  <- readRDS(file.path(BRANCH, "Cache", tt$trait, GDSL_ID, paste0("MGmin_", MGMIN),
                           paste0("eps_", EPS), "HapObject.rds"))
  hh  <- ho$HapObject[[paste0("Haplotypes_MGmin", MGMIN, "_E", EPS)]]
  ind <- as.data.table(hh$Indfile)[, .(Ind = as.character(Ind), hap = as.character(hap), Pheno = as.numeric(Pheno))]
  ind <- ind[hap != "0" & !is.na(Pheno)]
  gs  <- ind[, .(n = .N, mean_pheno = mean(Pheno)), by = hap][order(hap)]
  g   <- GR7[trait == tt$trait]
  stopifnot(nrow(g) == 1L, g$haplotype_groups == nrow(gs),
            identical(as.integer(gs$n), as.integer(strsplit(g$group_sizes, "|", fixed = TRUE)[[1]])))
  vf  <- as.data.frame(hh$Varfile)
  br  <- PW7[trait == tt$trait][, .(group1 = g1, group2 = g2, y.position = y, p.adj.signif = sym)]
  G7[[tt$trait]] <- c(tt, list(m = list(group_stats = as.data.frame(gs), indfile = as.data.frame(ind)),
                               br = br, stat = "P", sval = g$kw_p_raw, eta2 = g$eta_squared,
                               mg = vf$ID[vf$MGs != "0"]))
  groups7[[tt$trait]] <- setNames(ind$hap, ind$Ind)[order(ind$Ind)]
  cat(sprintf("[fig3] GDSL %-6s %d groups (%s) | KW P %.2e | eta2 %.3f\n", tt$trait, nrow(gs),
              paste(gs$n, collapse = "|"), g$kw_p_raw, g$eta_squared))
}
stopifnot(identical(groups7$fiber, groups7$starch))          # same accessions, same groups

gs7   <- G7$fiber$m$group_stats
cons7 <- t(vapply(gs7$hap, function(h) {
  ids <- names(groups7$fiber)[groups7$fiber == h]
  stopifnot(all(ids %in% rownames(W7)))
  consensus_row(W7[ids, shared7, drop = FALSE])
}, numeric(length(shared7))))
rownames(cons7) <- gs7$hap
stopifnot(!any(cons7 %in% CODE_HET))                          # wild consensus: REF/ALT/NA only

# Genome-wide significant SNPs (either trait) — the triangles
sig_pos <- SIGT7[sig_fiber == TRUE | sig_starch == TRUE, pos]
stopifnot(length(sig_pos) == 3L, all(sig_pos %in% spos7),
          all(paste0("7H:", sig_pos) %in% G7$fiber$mg))       # all three are in the defining group
# one barcode for both traits: the consensus and elite matrices go with the fiber entry
G7$fiber$m$consensus <- cons7; G7$fiber$m$elite <- EL7
G7$fiber$sig_col <- match(sig_pos, spos7)
cat(sprintf("[fig3] GDSL   %d shared SNPs | significant SNP columns: %s | %d elite rows\n",
            length(shared7), paste(G7$fiber$sig_col, collapse = ", "), nrow(EL7)))

# Barcode entries, one per row: the three genes, then the GDSL gene
BARS <- c(P, list(GDSL = G7$fiber))

# States actually drawn anywhere in the figure — the legend lists only these (step 07 convention).
drawn <- sort(unique(unlist(lapply(BARS, function(p) c(as.vector(p$m$consensus), as.vector(p$m$elite))))))
drawn <- drawn[!is.na(drawn)]
has_missing <- any(vapply(BARS, function(p) anyNA(p$m$consensus) || anyNA(p$m$elite), logical(1)))

# ── Drawing helpers ───────────────────────────────────────────────────────────
mm2in <- function(mm) mm / 25.4

# Gene-level statistic in the journal's style: q = 5.3 x 10^-6 (step 04) or P = 2.1 x 10^-3 (GDSL)
stat_expr <- function(p) {
  e <- floor(log10(p$sval)); mant <- p$sval / 10^e
  if (p$stat == "q")
    bquote(italic(q) == .(sprintf("%.1f", mant)) %*% 10^.(e) * "," ~~ italic(eta)^2 == .(sprintf("%.3f", p$eta2)))
  else
    bquote(italic(P) == .(sprintf("%.1f", mant)) %*% 10^.(e) * "," ~~ italic(eta)^2 == .(sprintf("%.3f", p$eta2)))
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

# Barcode; with sig_col (GDSL), one extra row at the bottom carries the significance triangles.
barcode_panel <- function(p, cex_ax = 1) {
  cons <- p$m$consensus; el <- p$m$elite
  n_c <- ncol(cons); rows <- nrow(cons) + 1L + nrow(el)     # + one blank separator row
  mat  <- rbind(cons, rep(NA_real_, n_c), el)
  labs <- c(paste("Group", rownames(cons)), "", rownames(el))
  blank <- nrow(cons) + 1L
  extra <- if (length(p$sig_col)) TRI_ROWS else 0

  plot(NA, xlim = c(0.5, n_c + 0.5), ylim = c(rows + 0.5 + extra, 0.5), xaxs = "i", yaxs = "i",
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
  if (extra) points(p$sig_col, rep(rows + 0.5 + extra / 2, length(p$sig_col)),
                    pch = SIG_PCH, col = SIG_COL, cex = 0.9)
  axis(2, at = seq_len(rows)[-blank], labels = labs[-blank], tick = FALSE,
       line = -0.15, las = 1, cex.axis = cex_ax)
}

# `name = TRUE` prints the gene and trait after the letter; the barcode panel of the same
# row carries the letter alone, since the row header already names the gene.
# `two_line = TRUE` (the half-width GDSL panels): the statistic goes on its own line under the
# name, left-aligned with it, since the two do not fit side by side in 43 mm at 9 pt.
panel_label <- function(p, letter, name = TRUE, stats = TRUE, two_line = FALSE,
                        dy_in = 0.06, dx_in = 0, cex = 1) {
  yin  <- grconvertY(1, "npc", "inches") + dy_in
  ytxt <- grconvertY(yin + if (two_line) 1.1 * LINE_IN else 0, "inches", "user")
  xl   <- grconvertX(grconvertX(0, "nfc", "inches") + dx_in, "inches", "user")
  lab  <- if (name) bquote(bold(.(letter)) ~~ .(p$label) ~ "(" * .(p$trait_lab) * ")")
          else      bquote(bold(.(letter)))
  text(xl, ytxt, lab, adj = c(0, 0), xpd = NA, cex = cex)
  if (stats && two_line)
    text(xl, grconvertY(yin, "inches", "user"), stat_expr(p), adj = c(0, 0), xpd = NA, cex = cex)
  else if (stats)
    text(grconvertX(grconvertX(1, "nfc", "inches") - 0.10, "inches", "user"),
         ytxt, stat_expr(p), adj = c(1, 0), xpd = NA, cex = cex)
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
  pin  <- par("pin"); tot <- sum(wid) + sep * (length(labs) - 1)
  # centred, but never closer to the right page edge than the barcode tiles end (0.6 lines):
  # in the 88-mm column the 85.4-mm row would otherwise end 0.5 mm from the edge (2026-10-17)
  x_in <- min((pin[1] - tot) / 2, pin[1] - tot - 0.6 * LINE_IN)
  ux <- function(i) i / pin[1]; uy <- function(i) i / pin[2]
  for (j in seq_along(labs)) {
    if (!is.na(fill[j])) rect(ux(x_in), y - uy(key / 2), ux(x_in + key), y + uy(key / 2),
                              col = fill[j], border = "black", lwd = LWD, xpd = NA)
    if (!is.na(pch[j]))  points(ux(x_in + key / 2), y, pch = pch[j], col = pcol[j], cex = 0.9, xpd = NA)
    text(ux(x_in + key + gap), y, labs[j], adj = c(0, 0.5), cex = cex, xpd = NA)
    x_in <- x_in + wid[j] + sep
  }
}

# The triangle key, drawn inside the GDSL barcode cell directly under its barcode (as in the
# former Fig. 4): same keys and sizes as hlegend(), positioned in device inches.
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

draw_legend <- function(cex = 1, mar = c(0, 0, 0, 0)) {
  par(mar = mar); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  keys <- as.character(drawn); labs <- unname(GT_LAB[keys]); cols <- unname(GT_COL[keys])
  if (has_missing) { labs <- c(labs, "Missing"); cols <- c(cols, COL_MISS) }
  hlegend(labs, fill = cols, cex = cex)
}

ylab_of <- function(p) paste(if (p$trait == "fiber") "Fiber" else "Starch", "BLUP")

# ── Layout: one row per gene — violins left, barcodes right ───────────────────
H_IN       <- floor(H_MAX / 25.4 * DPI) / DPI   # the full TAG height, whole pixels at 600 dpi (<= 234 mm)
LEG_IN     <- 0.34                               # genotype legend row, as the former Fig. 3
ROW_IN     <- (H_IN - LEG_IN) / 4                # four equal gene rows
COL_MM     <- c(43, 43, 88)                      # violin column 86 mm (two 43-mm halves in the GDSL row) | barcodes 88 mm
BAR_ROW_IN <- 0.175                              # barcode row pitch, as the former Figs. 3 and 4
KEY_IN     <- 0.24                               # triangle key under the GDSL barcode

draw_combined <- function() {
  M_TOP <- 1.5; M_TOP2 <- 2.6; M_BOT <- 2.6      # M_TOP2: the two-line headers of the GDSL panels
  layout(rbind(c(1, 1, 2),
               c(3, 3, 4),
               c(5, 5, 6),
               c(7, 8, 9),
               c(10, 10, 11)),                   # panel 10 empty: the legend sits under the barcode column only
         widths = COL_MM, heights = c(rep(ROW_IN, 4), LEG_IN))
  par(cex = 1)                                   # undo layout()'s automatic 0.66 shrink
  par(family = FONT, las = 1, mgp = c(2.0, 0.45, 0), tcl = -0.25,
      cex.axis = 1, cex.lab = 1, lwd = LWD, xpd = FALSE)

  # rows 1-3: GPAT6, GH17, PHT4;3 — violins a-c, barcodes f-h
  for (i in seq_along(P)) {
    p <- P[[i]]
    par(mar = c(M_BOT, 3.8, M_TOP, 0.6))
    violin_panel(p, ylab_of(p))
    panel_label(p, letters[i])
    # barcode centred vertically in the row, so it does not stretch with the number of groups
    nr  <- nrow(p$m$consensus) + 1L + nrow(p$m$elite)
    pad <- (ROW_IN - (M_TOP + M_BOT) * LINE_IN - nr * BAR_ROW_IN) / 2 / LINE_IN
    par(mar = c(M_BOT + pad, 5.2, M_TOP + pad, 0.6))
    barcode_panel(p)
    panel_label(p, letters[i + 5L], name = FALSE, stats = FALSE, dx_in = 0.30)
  }

  # row 4: GDSL — violins d (fiber) and e (starch) side by side, barcode i
  for (j in seq_along(G7)) {
    p <- G7[[j]]
    par(mar = c(M_BOT, 3.8, M_TOP2, 0.6))
    violin_panel(p, ylab_of(p))
    panel_label(p, letters[3L + j], two_line = TRUE)
  }
  p     <- G7$fiber
  nr    <- nrow(p$m$consensus) + 1L + nrow(p$m$elite) + TRI_ROWS
  bar_h <- nr * BAR_ROW_IN
  free  <- ROW_IN - (M_TOP + M_BOT) * LINE_IN - bar_h - KEY_IN   # barcode + key centred like rows 1-3
  par(mai = c(M_BOT * LINE_IN + free / 2 + KEY_IN, 5.2 * LINE_IN, M_TOP * LINE_IN + free / 2, 0.6 * LINE_IN))
  barcode_panel(p)
  panel_label(p, "i", name = FALSE, stats = FALSE, dx_in = 0.30)
  xc  <- mm2in(COL_MM[1] + COL_MM[2] + COL_MM[3] / 2)            # centre of the barcode column
  bot <- grconvertY(0, "npc", "inches")                          # bottom of the barcode plot region
  hlegend_dev("Genome-wide significant SNP", xc, bot - KEY_IN / 2 - 0.02, pch = SIG_PCH, pcol = SIG_COL)

  par(mar = c(0, 0, 0, 0)); plot.new()           # spacer under the violin column
  draw_legend(mar = c(0, 0, 0, 0))                # centred on the whole barcode column
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

render(OUT_NAME, draw_combined, H_IN)
