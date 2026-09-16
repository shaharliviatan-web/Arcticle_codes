#!/usr/bin/env Rscript
# 07b_manhattan_compare_v2_v3.R
# Manhattan comparison between the two population-structure corrections.
#
# For every trait (headline config BLUP x 3 PCs, plus a full 24-cell sweep) draws
# a MIRRORED Manhattan: v3 upward, v2 downward on a shared genomic axis. A mirror
# is used rather than an overlay because 7.1M overplotted points hide whichever
# series is drawn second; mirroring keeps both readable and makes gained/lost
# peaks obvious at a glance.
#
# Writes into results/comparison_v2_vs_v3/03_manhattan/:
#   v2_vs_v3/mirror_BLUP_pc3.{pdf,png}        headline, 4 traits stacked
#   v2_vs_v3/mirror_<pheno>_pc<N>.{pdf,png}   one per config (6 configs)
#   v2_vs_v3/hit_counts_compare.tsv           hits per trait/config, both versions
#   v3_only/manhattan_grid_<pheno>_pc<N>...   handled by 09_comparison_views.R
#
# Thresholds: each version is drawn against ITS OWN Bonferroni line
#   v2 = -log10(0.10/590462) = 6.7712     v3 = -log10(0.10/111017) = 6.0454
# and the table additionally scores v2 at the v3 line so the effect of the new
# correction can be separated from the effect of the easier threshold.
#
# Created 2026-08-16 for the v3 re-run review gate.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(parallel) })

PIPE    <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
PS_V3   <- file.path(PIPE, "results", "emmax_ps")
PS_V2   <- file.path(PIPE, "results", "_archive", "_archive_v2_win50snp_2026-08-16", "emmax_ps")
BIM     <- file.path(PIPE, "intermediates", "snp_map.bim")
OUT_DIR <- file.path(PIPE, "results", "comparison_v2_vs_v3", "03_manhattan", "v2_vs_v3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(PIPE, "scripts", "helpers", "tag_theme.R"))

BONF_V2 <- -log10(0.10 / 590462L)
BONF_V3 <- -log10(0.10 / 111017L)
COL_V2  <- c("#B0B0B0", "#8A8A8A")   # alternating chromosome shades, grey  = v2
COL_V3  <- c("#5B9BD5", "#1F6FB4")   # alternating chromosome shades, blue  = v3
COL_SIG <- "#D62728"
TRAITS  <- c("betaglucan", "fiber", "protein", "starch")
PLOT_FLOOR <- 2   # only plot SNPs with -log10p >= this (keeps files sane)

# ---- Genomic coordinates ----
cat("[07b] Loading snp_map.bim ...\n")
map <- fread(BIM, select = c(1, 2, 4), col.names = c("chr", "snp", "bp"))
map[, chr := sub("H$", "", chr)][, chr := as.integer(chr)]
setorder(map, chr, bp)
chr_len <- map[, .(len = max(bp)), by = chr][order(chr)]
chr_len[, offset := cumsum(as.numeric(len)) - len]
map <- merge(map, chr_len[, .(chr, offset)], by = "chr", sort = FALSE)
map[, gpos := as.numeric(bp) + offset]
setkey(map, snp)
axis_at  <- chr_len[, offset + len / 2]
axis_lab <- paste0(chr_len$chr, "H")
gmax     <- chr_len[, max(offset + len)]

read_ps <- function(dir, trait, pheno, pc) {
  f <- file.path(dir, sprintf("morexV3__%s__%s__pc%d.ps", trait, pheno, pc))
  if (!file.exists(f)) return(NULL)
  d <- fread(f, select = c(1, 4), col.names = c("snp", "p"), showProgress = FALSE)
  d <- d[!is.na(p) & p > 0 & p <= 1]
  d[, lp := -log10(p)]
  d <- d[lp >= PLOT_FLOOR]
  d <- merge(d, map[, .(snp, chr, gpos)], by = "snp", sort = FALSE)
  setorder(d, chr, gpos)
  d[]
}

draw_mirror <- function(v3, v2, trait, pheno, pc, ymax) {
  par(mar = c(3.2, 4.2, 2.2, 0.8), las = 1, mgp = c(2.6, 0.6, 0),
      cex.axis = 0.75, cex.lab = 0.85, cex.main = 0.9)
  plot(NA, xlim = c(0, gmax), ylim = c(-ymax, ymax), xaxt = "n", yaxt = "n",
       xlab = "", ylab = expression(-log[10](p)),
       main = sprintf("%s | %s | %d PCs      v3 (up)  vs  v2 (down)", trait, pheno, pc))
  axis(1, at = axis_at, labels = axis_lab, tick = FALSE)
  yt <- pretty(c(0, ymax), 4); yt <- yt[yt <= ymax]
  axis(2, at = c(-rev(yt), yt), labels = c(rev(yt), yt))
  abline(h = 0, col = "grey40", lwd = 0.7)

  if (!is.null(v3)) points(v3$gpos,  v3$lp, pch = 16, cex = 0.22, col = COL_V3[v3$chr %% 2 + 1])
  if (!is.null(v2)) points(v2$gpos, -v2$lp, pch = 16, cex = 0.22, col = COL_V2[v2$chr %% 2 + 1])
  abline(h =  BONF_V3, col = COL_SIG, lty = 2, lwd = 0.9)
  abline(h = -BONF_V2, col = COL_SIG, lty = 2, lwd = 0.9)
  text(gmax, BONF_V3,  sprintf(" v3 Bonf %.2f", BONF_V3), adj = c(1, -0.4), cex = 0.62, col = COL_SIG)
  text(gmax, -BONF_V2, sprintf(" v2 Bonf %.2f", BONF_V2), adj = c(1, 1.3),  cex = 0.62, col = COL_SIG)
}

hit_rows <- list()

for (ph in c("BLUP", "BLUE")) for (pc in c(3L, 5L, 10L)) {
  cat(sprintf("[07b] %s pc%d ...\n", ph, pc))
  dat <- mclapply(TRAITS, function(tr)
           list(trait = tr, v3 = read_ps(PS_V3, tr, ph, pc), v2 = read_ps(PS_V2, tr, ph, pc)),
         mc.cores = 4)

  for (x in dat) {
    hit_rows[[length(hit_rows) + 1L]] <- data.table(
      trait = x$trait, pheno_type = ph, n_PCs = pc,
      v3_hits_at_v3thr = if (is.null(x$v3)) NA_integer_ else sum(x$v3$lp >= BONF_V3),
      v2_hits_at_v2thr = if (is.null(x$v2)) NA_integer_ else sum(x$v2$lp >= BONF_V2),
      v2_hits_at_v3thr = if (is.null(x$v2)) NA_integer_ else sum(x$v2$lp >= BONF_V3),
      v3_top = if (is.null(x$v3)) NA_real_ else max(x$v3$lp),
      v2_top = if (is.null(x$v2)) NA_real_ else max(x$v2$lp))
  }

  ymax <- max(unlist(lapply(dat, function(x)
             c(if (is.null(x$v3)) 0 else max(x$v3$lp),
               if (is.null(x$v2)) 0 else max(x$v2$lp)))), BONF_V2) * 1.06

  base <- sprintf("mirror_%s_pc%d", ph, pc)
  tag_plot_dual(file.path(OUT_DIR, paste0(base, ".pdf")),
                file.path(OUT_DIR, paste0(base, ".png")),
    function() {
      par(mfrow = c(4, 1), oma = c(0.5, 0, 1.4, 0))
      for (x in dat) draw_mirror(x$v3, x$v2, x$trait, ph, pc, ymax)
    }, width_in = 11, height_in = 12)
  cat(sprintf("[07b] OK  %s.{pdf,png}\n", base))
}

hits <- rbindlist(hit_rows)
setorder(hits, pheno_type, n_PCs, trait)
fwrite(hits, file.path(OUT_DIR, "hit_counts_compare.tsv"), sep = "\t")

cat("\n---------------------------------\n")
cat(sprintf("[07b] CHECKPOINT: configs drawn = %d (expected 6)\n",
            length(unique(paste(hits$pheno_type, hits$n_PCs)))))
cat(sprintf("[07b] CHECKPOINT: rows in hit_counts_compare.tsv = %d (expected 24)\n", nrow(hits)))
cat("\n[07b] Headline (BLUP x 3 PCs):\n")
print(hits[pheno_type == "BLUP" & n_PCs == 3L])
cat(sprintf("\n[07b] OK: outputs in %s\n", OUT_DIR))
