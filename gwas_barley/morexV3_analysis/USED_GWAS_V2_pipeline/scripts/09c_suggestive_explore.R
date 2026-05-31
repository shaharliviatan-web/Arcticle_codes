#!/usr/bin/env Rscript
# 09c_suggestive_explore.R
# Exploratory Manhattan plots for the chosen config (BLUE x 10 PCs):
#   - Solid line at BonfPruned010 = 6.7712 (main threshold; alpha=0.10, N=590,462)
#   - 4 suggestive dashed lines at -log10p = 6.7, 6.6, 6.5, 6.4
#   - Legend reports per-suggestive cumulative count: SNPs in [X, BonfPruned010)
#
# Outputs to results/suggestive_explore/:
#   manhattan_explore__{trait}__BLUE__pc10.{pdf,png}    (4 plots)
#   suggestive_added_snps_all_traits.tsv                (every SNP in [6.4, 6.7712))

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(png)
  library(parallel)
})

PIPE      <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
INTER_DIR <- file.path(PIPE, "intermediates")
PS_DIR    <- file.path(PIPE, "results", "emmax_ps")
OUT_DIR   <- file.path(PIPE, "results", "suggestive_explore")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
BIM_FILE  <- file.path(INTER_DIR, "snp_map.bim")

# --- Thresholds ---
N_PRUNED  <- 590462L
ALPHA     <- 0.10
BONF      <- -log10(ALPHA / N_PRUNED)            # 6.7712
SUGGS     <- c(6.7, 6.6, 6.5, 6.4)               # suggestive lines (descending)
SUGG_COLS <- c("dodgerblue3", "darkorange2", "purple3", "firebrick3")

# --- Plot constants (same as step 07) ---
PLOT_W_IN <- 12; PLOT_H_IN <- 7
POINT_CEX <- 0.85; RASTER_DPI <- 900; PNG_OUT_DPI <- 600
CHR_COLS  <- c("blue4", "orange3")

# --- SNP map ---
snp_map <- fread(BIM_FILE, header = FALSE, select = c(1, 2, 4),
                 col.names = c("CHR_raw", "SNP", "BP"))
snp_map[, CHR := as.integer(sub("H$", "", as.character(CHR_raw)))]
stopifnot(all(snp_map$CHR %in% 1:7))

TRAITS <- c("betaglucan", "fiber", "protein", "starch")

# --- Per-trait processor ---
process_trait <- function(tr) {
  ps_file <- file.path(PS_DIR, sprintf("morexV3__%s__BLUE__pc10.ps", tr))
  gw <- fread(ps_file, select = 4, col.names = "P", showProgress = FALSE)
  d <- data.table(CHR = snp_map$CHR, SNP = snp_map$SNP, BP = snp_map$BP, P = gw$P)
  d <- d[!is.na(P) & P > 0 & P <= 1]
  d[, nlp := -log10(P)]

  # SNPs above each level (cumulative)
  n_above_bonf <- sum(d$nlp >= BONF)
  n_added <- vapply(SUGGS, function(x) sum(d$nlp >= x & d$nlp < BONF), integer(1))

  # SNPs in [6.4, BonfPruned010) -- the suggestive band
  band <- d[nlp >= 6.4 & nlp < BONF][order(-nlp)]
  band[, trait := tr]
  band[, captured_by := vapply(nlp, function(v) {
    s <- SUGGS[v >= SUGGS]   # lines that capture this SNP
    if (length(s) == 0) NA_character_ else sprintf("%.1f", min(s))  # lowest line wins
  }, character(1))]

  # Manhattan x-coords
  setorder(d, CHR, BP)
  chr_max <- d[, .(mx = max(BP)), by = CHR][order(CHR)]
  chr_max[, off := cumsum(as.numeric(shift(mx, fill = 0)))]
  d <- merge(d, chr_max[, .(CHR, off)], by = "CHR"); setorder(d, CHR, BP)
  d[, x := as.numeric(BP) + off]
  d[, col := ifelse(CHR %% 2L == 0L, CHR_COLS[2], CHR_COLS[1])]
  chr_cent <- d[, .(c = (min(x) + max(x)) / 2), by = CHR][order(CHR)]
  xlim <- c(min(d$x), max(d$x))
  max_obs <- max(d$nlp, na.rm = TRUE)
  ymax <- ceiling(max(max_obs, BONF, 6, na.rm = TRUE)) + 1
  ylim <- c(0, ymax)

  # Render point cloud raster
  tmp_png <- tempfile(fileext = ".png")
  png(tmp_png, width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
      res = RASTER_DPI, bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(d$x, d$nlp, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "", ylab = "", pch = 16, cex = POINT_CEX, col = d$col)
  dev.off()
  img <- readPNG(tmp_png); unlink(tmp_png)

  # Composite
  draw <- function() {
    title <- sprintf("Manhattan: %s | BLUE | 10 PCs   (BonfPruned010 = 6.77; n above = %d)",
                     tools::toTitleCase(tr), n_above_bonf)
    par(mar = c(4.5, 4.8, 3, 1.2), las = 1, cex.axis = 0.95, cex.lab = 1.0)
    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "Chromosome", ylab = expression(-log[10](italic(p))), main = title)
    rasterImage(img, xlim[1], ylim[1], xlim[2], ylim[2], interpolate = FALSE)
    axis(2); axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE); box()
    # Suggestive lines first (so the Bonf line sits on top)
    for (i in seq_along(SUGGS)) {
      abline(h = SUGGS[i], col = SUGG_COLS[i], lty = 2, lwd = 1.4)
    }
    # Main Bonf line: solid black, thick
    abline(h = BONF, col = "black", lty = 1, lwd = 1.8)
    legend("topright",
           legend = c(
             sprintf("BonfPruned010 = %.4f  (n above = %d)", BONF, n_above_bonf),
             sprintf("suggestive %.1f: +%d SNPs above (cumulative)", SUGGS[1], n_added[1]),
             sprintf("suggestive %.1f: +%d SNPs above (cumulative)", SUGGS[2], n_added[2]),
             sprintf("suggestive %.1f: +%d SNPs above (cumulative)", SUGGS[3], n_added[3]),
             sprintf("suggestive %.1f: +%d SNPs above (cumulative)", SUGGS[4], n_added[4])),
           col = c("black", SUGG_COLS),
           lty = c(1, 2, 2, 2, 2),
           lwd = c(1.8, 1.4, 1.4, 1.4, 1.4),
           cex = 0.78, bty = "n", inset = c(0.005, 0.01))
  }

  base <- sprintf("manhattan_explore__%s__BLUE__pc10", tr)
  cairo_pdf(file.path(OUT_DIR, paste0(base, ".pdf")),
            width = PLOT_W_IN, height = PLOT_H_IN); draw(); dev.off()
  png(file.path(OUT_DIR, paste0(base, ".png")),
      width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
      res = PNG_OUT_DPI, type = "cairo-png"); draw(); dev.off()

  list(trait = tr, n_above_bonf = n_above_bonf,
       n_added_67 = n_added[1], n_added_66 = n_added[2],
       n_added_65 = n_added[3], n_added_64 = n_added[4],
       band = band)
}

cat(sprintf("[09c] BonfPruned010 = %.4f (alpha=%.2f / %d)\n", BONF, ALPHA, N_PRUNED))
cat(sprintf("[09c] Suggestive lines: %s\n", paste(SUGGS, collapse = ", ")))
cat("[09c] Rendering 4 Manhattan plots (raster point cloud, 12x7 in)...\n")
t0 <- Sys.time()
res <- mclapply(TRAITS, process_trait, mc.cores = 4)
cat(sprintf("[09c] elapsed: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# --- Per-trait counts table to stdout ---
counts <- rbindlist(lapply(res, function(x) {
  data.table(trait = x$trait, n_above_BonfPruned010 = x$n_above_bonf,
             added_by_6.7 = x$n_added_67, added_by_6.6 = x$n_added_66,
             added_by_6.5 = x$n_added_65, added_by_6.4 = x$n_added_64)
}))
cat("\n[09c] Cumulative added-SNP counts (BLUE x pc10):\n")
print(counts)
fwrite(counts, file.path(OUT_DIR, "suggestive_counts_per_trait.tsv"), sep = "\t")

# --- Full table of SNPs in the [6.4, BonfPruned010) suggestive band ---
all_band <- rbindlist(lapply(res, function(x) x$band))
all_band[, p_value := signif(P, 4)]
all_band[, neg_log10_p := round(nlp, 4)]
out_cols <- c("trait", "CHR", "SNP", "BP", "neg_log10_p", "p_value", "captured_by")
all_band <- all_band[, ..out_cols]
setorder(all_band, trait, -neg_log10_p)
fwrite(all_band, file.path(OUT_DIR, "suggestive_added_snps_all_traits.tsv"), sep = "\t")
cat(sprintf("\n[09c] Wrote suggestive_added_snps_all_traits.tsv (%d rows)\n", nrow(all_band)))
cat("[09c] OK\n")
