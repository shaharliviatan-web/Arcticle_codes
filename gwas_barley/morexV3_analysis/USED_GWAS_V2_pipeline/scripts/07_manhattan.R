#!/usr/bin/env Rscript
# 07_manhattan.R
# For each of the 24 EMMAX runs, generate 2 Manhattan plots (FDR q=0.05 and q=0.10),
# each with BOTH Bonferroni reference lines (all-7.1M and LD-pruned-590K) drawn,
# in BOTH vector PDF and PNG. Total: 24 x 2 x 2 = 96 files.
#
# Rendering: the full 7.1M-point cloud is drawn ONCE per cell to a high-res (900 dpi)
# transparent raster layer (dense + continuous, no thinning artifacts), then composited
# under TRUE-VECTOR axes / labels / threshold lines / legend in the PDF. This keeps files
# small (~0.3 MB PDF) while looking like a published Manhattan. The same raster is reused
# across both FDR variants and both output formats.

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(data.table)
  library(png)
  library(parallel)
})

# --- Paths ---
PIPE       <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_GWAS_V2_pipeline"
INTER_DIR  <- file.path(PIPE, "intermediates")
PS_DIR     <- file.path(PIPE, "results", "emmax_ps")
MAN_DIR    <- file.path(PIPE, "results", "manhattan")
dir.create(MAN_DIR, showWarnings = FALSE, recursive = TRUE)

BIM_FILE   <- file.path(INTER_DIR, "snp_map.bim")
PRUNE_IN   <- file.path(INTER_DIR, "morexV3_pruned_for_covs.prune.in")

# --- Plot style constants (chosen with the user) ---
# Per-cell WORKING Manhattans at 12x7 in (large, easy to inspect). The TAG-spec
# publication figures (174 mm) are produced separately in step 10.
PLOT_W_IN   <- 12      # device width (inches)
PLOT_H_IN   <- 7       # device height (inches)
POINT_CEX   <- 0.85    # dot size
RASTER_DPI  <- 900     # resolution of the rasterized point-cloud layer (crisp)
PNG_OUT_DPI <- 600     # resolution of the standalone PNG output
CHR_COLS    <- c("blue4", "orange3")

# --- Load SNP map ONCE ---
cat("[07] Loading snp_map.bim ...\n")
snp_map <- fread(BIM_FILE, header = FALSE,
                 select = c(1, 2, 4),
                 col.names = c("CHR_raw", "SNP", "BP"))
snp_map[, CHR := as.integer(sub("H$", "", sub("^chr", "", as.character(CHR_raw))))]
snp_map[is.na(CHR), CHR := as.integer(CHR_raw) - 30L]
stopifnot(all(snp_map$CHR %in% 1:7))
cat(sprintf("[07] snp_map: %s SNPs, chr %s\n",
            format(nrow(snp_map), big.mark = ","),
            paste(sort(unique(snp_map$CHR)), collapse = ",")))

# --- Bonferroni thresholds ---
n_all    <- nrow(snp_map)
n_pruned <- length(readLines(PRUNE_IN))
bonf_all    <- -log10(0.05 / n_all)
bonf_pruned <- -log10(0.05 / n_pruned)
cat(sprintf("[07] Bonf_all=%.3f  Bonf_pruned=%.3f\n", bonf_all, bonf_pruned))

# --- Discover + parse .ps files ---
ps_files <- list.files(PS_DIR, pattern = "^morexV3__.+__.+__pc\\d+\\.ps$", full.names = TRUE)
stopifnot(length(ps_files) == 24)
parts <- as.data.table(do.call(rbind, strsplit(basename(ps_files), "__|\\.", perl = FALSE)))
setnames(parts, c("prefix", "trait", "pheno_type", "pcs_str", "ext"))
parts[, pcs := as.integer(sub("^pc", "", pcs_str))]
parts[, path := ps_files]

# --- Per-cell processing ---
process_one <- function(i) {
  row <- parts[i, ]
  base <- sprintf("manhattan__%s__%s__pc%d", row$trait, row$pheno_type, row$pcs)

  gw <- fread(row$path, select = 4, col.names = "P", showProgress = FALSE)
  if (nrow(gw) != nrow(snp_map)) {
    stop(sprintf("[07] FAIL: %s has %d rows, snp_map has %d", basename(row$path), nrow(gw), nrow(snp_map)))
  }
  d <- data.table(CHR = snp_map$CHR, SNP = snp_map$SNP, BP = snp_map$BP, P = gw$P)
  d <- d[!is.na(P) & P > 0 & P <= 1]
  d[, nlp := -log10(P)]

  # FDR via Benjamini-Hochberg (on FULL data)
  q_vals <- p.adjust(d$P, method = "BH")
  fdr_info <- lapply(c(0.05, 0.10), function(qv) {
    sig <- which(q_vals <= qv)
    if (length(sig) > 0) list(q = qv, thr = -log10(max(d$P[sig], na.rm = TRUE)),
                              n_sig = length(sig), valid = TRUE)
    else list(q = qv, thr = NA_real_, n_sig = 0L, valid = FALSE)
  })

  # Manhattan layout: cumulative x positions, chromosome centres, alternating colours
  setorder(d, CHR, BP)
  chr_max <- d[, .(mx = max(BP)), by = CHR][order(CHR)]
  chr_max[, off := cumsum(as.numeric(shift(mx, fill = 0)))]
  d <- merge(d, chr_max[, .(CHR, off)], by = "CHR"); setorder(d, CHR, BP)
  d[, x := as.numeric(BP) + off]
  d[, col := ifelse(CHR %% 2L == 0L, CHR_COLS[2], CHR_COLS[1])]
  chr_cent <- d[, .(c = (min(x) + max(x)) / 2), by = CHR][order(CHR)]
  xlim <- c(min(d$x), max(d$x))
  max_obs <- max(d$nlp, na.rm = TRUE)
  ymax <- ceiling(max(max_obs, bonf_all, bonf_pruned, 6, na.rm = TRUE)) + 1
  ylim <- c(0, ymax)
  pheno_disp <- sprintf("%s | %s | %d PCs", tools::toTitleCase(row$trait), row$pheno_type, row$pcs)

  # --- Render the point cloud ONCE to a transparent high-res raster ---
  tmp_png <- tempfile(fileext = ".png")
  png(tmp_png, width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
      res = RASTER_DPI, bg = "transparent", type = "cairo-png")
  par(mar = c(0, 0, 0, 0))
  plot(d$x, d$nlp, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = "", ylab = "", pch = 16, cex = POINT_CEX, col = d$col)
  dev.off()
  img <- readPNG(tmp_png)

  # --- Composite drawing: vector axes + raster points + vector lines/legend ---
  draw_man <- function(fi) {
    title <- sprintf("Manhattan: %s   (FDR q=%.2f)", pheno_disp, fi$q)
    par(mar = c(4.5, 4.8, 3, 1.2), las = 1, cex.axis = 0.95, cex.lab = 1.0)
    plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i", axes = FALSE,
         xlab = "Chromosome", ylab = expression(-log[10](italic(p))), main = title)
    rasterImage(img, xlim[1], ylim[1], xlim[2], ylim[2], interpolate = FALSE)
    axis(2)
    axis(1, at = chr_cent$c, labels = paste0(chr_cent$CHR, "H"), tick = FALSE)
    box()
    abline(h = bonf_all,    lty = 1, lwd = 1.4, col = "black")
    abline(h = bonf_pruned, lty = 2, lwd = 1.2, col = "black")
    leg  <- c(sprintf("Bonf_all = %.2f", bonf_all),
              sprintf("Bonf_pruned = %.2f", bonf_pruned))
    cols <- c("black", "black"); ltys <- c(1, 2)
    if (fi$valid) {
      abline(h = fi$thr, col = "blue2", lty = 1, lwd = 1.2)
      leg  <- c(leg, sprintf("FDR %.2f = %.2f (n=%d)", fi$q, fi$thr, fi$n_sig))
      cols <- c(cols, "blue2"); ltys <- c(ltys, 1)
    } else {
      leg  <- c(leg, sprintf("FDR %.2f: no SNPs passed", fi$q))
      cols <- c(cols, "blue2"); ltys <- c(ltys, NA)
    }
    legend("topright", legend = leg, col = cols, lty = ltys, lwd = 1.2,
           cex = 0.8, bty = "n", inset = c(0.01, 0.02))
  }

  # Two FDR variants x two formats, reusing the same raster
  for (fi in fdr_info) {
    tag <- sprintf("FDR%s", sub("\\.", "", sprintf("%.2f", fi$q)))
    pdf_f <- file.path(MAN_DIR, sprintf("%s__%s.pdf", base, tag))
    png_f <- file.path(MAN_DIR, sprintf("%s__%s.png", base, tag))
    cairo_pdf(pdf_f, width = PLOT_W_IN, height = PLOT_H_IN); draw_man(fi); dev.off()
    png(png_f, width = PLOT_W_IN, height = PLOT_H_IN, units = "in",
        res = PNG_OUT_DPI, type = "cairo-png"); draw_man(fi); dev.off()
  }

  unlink(tmp_png)
  list(cell = sprintf("%s_%s_pc%d", row$trait, row$pheno_type, row$pcs),
       n_sig_fdr05 = fdr_info[[1]]$n_sig, n_sig_fdr10 = fdr_info[[2]]$n_sig)
}

# --- Run in parallel ---
N_CORES <- min(8, nrow(parts))
cat(sprintf("[07] Generating Manhattan plots for %d cells (%d cores)...\n", nrow(parts), N_CORES))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(parts)), process_one, mc.cores = N_CORES)
cat(sprintf("[07] mclapply elapsed: %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

err_idx <- which(vapply(res, inherits, logical(1), what = "try-error"))
if (length(err_idx) > 0) {
  for (i in err_idx) { cat(sprintf("[07] FAIL at cell %d (%s):\n", i, basename(parts$path[i]))); print(res[[i]]) }
  stop(sprintf("[07] FAIL: %d cells errored out", length(err_idx)))
}

# --- Checkpoints ---
cat("\n---------------------------------\n")
n_pdf <- length(list.files(MAN_DIR, pattern = "^manhattan__.+__FDR\\d+\\.pdf$"))
n_png <- length(list.files(MAN_DIR, pattern = "^manhattan__.+__FDR\\d+\\.png$"))
cat(sprintf("[07] CHECKPOINT: %d Manhattan PDFs (expected 48)\n", n_pdf))
cat(sprintf("[07] CHECKPOINT: %d Manhattan PNGs (expected 48)\n", n_png))
stopifnot(n_pdf == 48, n_png == 48)

cat("\n[07] Per-cell FDR hit counts:\n")
print(rbindlist(res))
cat("\n[07] OK: 96 Manhattan files generated (48 PDF + 48 PNG)\n")
