# =============================================================================
# 02_compute_decay_profile.R
#
# Reads the pre-aggregated per-chromosome TSV files produced by
# 01b_preaggregate_ld_parallel.sh (one row per 1-kb distance bin per chromosome,
# with SUM, SUM-OF-SQUARES, COUNT of pairwise r^2). Computes mean / SD / SE per
# bin (both per-chromosome and genome-wide), applies a pair-count-weighted LOESS
# smoother to the genome-wide profile, and finds the distances at which the
# smoothed curve crosses r^2 = 0.1, 0.2, 0.3.
#
# Outputs (written to ../results/):
#   - ld_decay_profile_genomewide.tsv
#   - ld_decay_profile_per_chromosome.tsv
#   - ld_decay_crossings.tsv
#   - ld_decay_loess_predictions.tsv
#   - ld_decay_loess_fit.rds
#   - ld_decay_summary_for_results_chapter.txt
#
# This script is fast (<1 min) because the heavy single-pass aggregation
# already happened in 01b_preaggregate_ld.sh.
# =============================================================================

# Force temp files off the system disk (small files only here, but safe).
dir.create("/mnt/data/shahar/.tmp", showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp",
           TMP    = "/mnt/data/shahar/.tmp",
           TEMP   = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
    library(data.table)
})

# --- 1. Paths -----------------------------------------------------------------
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE))
if (length(script_path) == 0) script_path <- "02_compute_decay_profile.R"
SCRIPT_DIR  <- normalizePath(dirname(script_path))
PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))

INTER_DIR <- file.path(PROJECT_DIR, "intermediates")
RES_DIR   <- file.path(PROJECT_DIR, "results")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# --- 2. Configuration ---------------------------------------------------------
CHROMS         <- c("1H", "2H", "3H", "4H", "5H", "6H", "7H")
BIN_SIZE_BP    <- 1000
MAX_DIST_BP    <- 2500000
R2_THRESHOLDS  <- c(0.1, 0.2, 0.3)
N_SAMPLES      <- 290

# --- 3. Load pre-aggregated per-chromosome tables -----------------------------
# Defensive reader:
#   * fill=TRUE so fread doesn't halt at malformed rows (a parallel-pipe race
#     in step 1b can occasionally tear ~1 line per chromosome at a byte boundary,
#     producing rows with <4 fields. We want to keep going past them.)
#   * Drop rows with any NA in the four numeric columns.
#   * Drop rows where Bin is not a non-negative multiple of BIN_SIZE_BP or > MAX_DIST_BP.
#   * Dedupe by Bin (sum the partials) in case a valid bin appears in more than
#     one row due to the same tearing artifact.
cat("Loading per-chromosome pre-aggregated tables ...\n")
per_chrom_list <- lapply(CHROMS, function(chr) {
    f <- file.path(INTER_DIR, sprintf("pre_aggregated_chr%s.tsv", chr))
    if (!file.exists(f)) stop(sprintf("Missing pre-aggregated file: %s", f))

    d_raw <- fread(f, fill = TRUE)

    n_raw <- nrow(d_raw)
    d <- d_raw[
        !is.na(Bin) & !is.na(Sum_R2) & !is.na(SumSq_R2) & !is.na(N_Pairs) &
        Bin >= 0 & Bin <= MAX_DIST_BP & (Bin %% BIN_SIZE_BP == 0)
    ]
    n_clean <- nrow(d)
    if (n_raw != n_clean) {
        cat(sprintf("  [%s] dropped %d malformed row(s) of %d\n", chr, n_raw - n_clean, n_raw))
    }

    # Dedupe by Bin (sum across any duplicate rows for the same bin)
    d <- d[, .(
        Sum_R2   = sum(Sum_R2),
        SumSq_R2 = sum(SumSq_R2),
        N_Pairs  = sum(N_Pairs)
    ), by = Bin]

    d[, Chromosome := chr]
    d
})
per_chrom <- rbindlist(per_chrom_list)
cat(sprintf("Loaded %d total bin rows across %d chromosomes.\n",
            nrow(per_chrom), length(CHROMS)))

# Derived per-chromosome statistics
per_chrom[, Avg_R2 := Sum_R2 / N_Pairs]
per_chrom[, Var_R2 := pmax(SumSq_R2 / N_Pairs - Avg_R2^2, 0)]
per_chrom[, SD_R2  := sqrt(Var_R2)]
per_chrom[, SE_R2  := SD_R2 / sqrt(N_Pairs)]
per_chrom[, Distance_kb := Bin / 1000]
setcolorder(per_chrom,
    c("Chromosome", "Bin", "Distance_kb", "Avg_R2", "SD_R2", "SE_R2", "N_Pairs",
      "Sum_R2", "SumSq_R2", "Var_R2"))

# --- 4. Genome-wide aggregation ----------------------------------------------
cat("Aggregating genome-wide profile ...\n")
gw <- per_chrom[, .(
    Sum_R2   = sum(Sum_R2),
    SumSq_R2 = sum(SumSq_R2),
    N_Pairs  = sum(N_Pairs)
), by = Bin]
gw[, Avg_R2 := Sum_R2 / N_Pairs]
gw[, Var_R2 := pmax(SumSq_R2 / N_Pairs - Avg_R2^2, 0)]
gw[, SD_R2  := sqrt(Var_R2)]
gw[, SE_R2  := SD_R2 / sqrt(N_Pairs)]
gw[, Distance_kb := Bin / 1000]
setorder(gw, Bin)
setcolorder(gw,
    c("Bin", "Distance_kb", "Avg_R2", "SD_R2", "SE_R2", "N_Pairs",
      "Sum_R2", "SumSq_R2", "Var_R2"))

# --- 5. LOESS smoothing of the binned mean r^2 -------------------------------
# Non-parametric smoother on the 1-kb-binned genome-wide profile. This is the
# standard approach in plant LD-decay papers (Comadran et al. 2009; Russell et
# al. 2016; etc.). Bins are weighted by their pair count (N_Pairs) so high-count
# bins dominate the local fit.
#
# We deliberately do NOT use the Hill-Weir parametric model here. In selfing
# inbred populations like wild barley the empirical r^2 has a non-zero
# asymptote (~0.10 in this dataset) from variance-component / population-
# structure effects, which the textbook Hill-Weir form (asymptote -> 1/n)
# cannot represent. LOESS makes no such assumption.
LOESS_SPAN <- 0.10
cat(sprintf("Fitting LOESS smoother (span = %.2f) on genome-wide bins ...\n", LOESS_SPAN))

fit_data <- gw[!is.na(Avg_R2) & N_Pairs >= 10 & Bin > 0]
loess_fit <- loess(
    Avg_R2 ~ Distance_kb,
    data    = fit_data,
    weights = N_Pairs,
    span    = LOESS_SPAN,
    degree  = 2,
    control = loess.control(surface = "direct")
)

# Dense 1-kb grid for prediction
pred_grid <- data.table(Distance_bp = seq(BIN_SIZE_BP, MAX_DIST_BP, by = BIN_SIZE_BP))
pred_grid[, Distance_kb := Distance_bp / 1000]
pred_grid[, R2_fit      := predict(loess_fit, newdata = pred_grid[, .(Distance_kb)])]

# --- 6. Crossings of the LOESS curve -----------------------------------------
find_crossing <- function(pred_dt, threshold) {
    below <- pred_dt[!is.na(R2_fit) & R2_fit < threshold]
    if (nrow(below) == 0) return(NA_real_)
    below[which.min(Distance_bp), Distance_kb]
}
crossings <- data.table(
    R2_threshold = R2_THRESHOLDS,
    Distance_kb  = sapply(R2_THRESHOLDS, function(t) find_crossing(pred_grid, t))
)
crossings[, Window_bp        := Distance_kb * 1000]
crossings[, Window_2sided_kb := 2 * Distance_kb]

# Empirical asymptote: mean r^2 over the last 100 kb of the analysed window
r2_baseline   <- mean(gw[Bin >= MAX_DIST_BP - 100000, Avg_R2], na.rm = TRUE)
# r^2 at the first usable bin (LOESS at 1 kb) -- treat as "peak" / r^2_max
r2_max        <- as.numeric(pred_grid[Distance_bp == BIN_SIZE_BP, R2_fit])
# Half-decay: where the LOESS curve falls to r2_max/2
half_decay_kb <- find_crossing(pred_grid, r2_max / 2)

# --- 7. Write outputs --------------------------------------------------------
cat("Writing outputs to ", RES_DIR, " ...\n", sep = "")

fwrite(gw,        file.path(RES_DIR, "ld_decay_profile_genomewide.tsv"),       sep = "\t")
fwrite(per_chrom, file.path(RES_DIR, "ld_decay_profile_per_chromosome.tsv"),   sep = "\t")
fwrite(crossings, file.path(RES_DIR, "ld_decay_crossings.tsv"),                sep = "\t")
fwrite(pred_grid, file.path(RES_DIR, "ld_decay_loess_predictions.tsv"),        sep = "\t")
saveRDS(list(loess_fit = loess_fit, span = LOESS_SPAN, n = N_SAMPLES,
             r2_max = r2_max, r2_baseline = r2_baseline,
             half_decay_kb = half_decay_kb,
             crossings = crossings, pred_grid = pred_grid),
        file.path(RES_DIR, "ld_decay_loess_fit.rds"))

# --- 8. Results-chapter summary block ----------------------------------------
total_pairs <- format(sum(gw$N_Pairs), big.mark = ",")
cross02     <- crossings[R2_threshold == 0.2, Distance_kb]

summary_lines <- c(
    "===========================================================",
    "  LD decay summary  (wild barley, morexV3, all 7 chromosomes)",
    "===========================================================",
    sprintf("  Samples used                          : %d", N_SAMPLES),
    sprintf("  Chromosomes                           : %s", paste(CHROMS, collapse = ", ")),
    sprintf("  Total binned SNP pairs                : %s", total_pairs),
    sprintf("  Physical window analysed              : %d kb", MAX_DIST_BP / 1000),
    sprintf("  Smoothing                             : LOESS, span = %.2f, pair-weighted", LOESS_SPAN),
    sprintf("  Smoothed r^2 at 1 kb (~r^2_max)       : %.3f", r2_max),
    sprintf("  Empirical r^2 asymptote (last 100 kb) : %.3f", r2_baseline),
    sprintf("  Half-decay distance (r^2_max / 2)     : %.1f kb", half_decay_kb),
    "  ----- Crossings of the LOESS curve -----",
    sprintf("  r^2 = 0.30  -> %6.1f kb  (window +/- = %6.1f kb)",
            crossings[R2_threshold == 0.3, Distance_kb],
            crossings[R2_threshold == 0.3, Window_2sided_kb]),
    sprintf("  r^2 = 0.20  -> %6.1f kb  (window +/- = %6.1f kb)  <- recommended GWAS window",
            crossings[R2_threshold == 0.2, Distance_kb],
            crossings[R2_threshold == 0.2, Window_2sided_kb]),
    sprintf("  r^2 = 0.10  -> %6.1f kb  (window +/- = %6.1f kb)",
            crossings[R2_threshold == 0.1, Distance_kb],
            crossings[R2_threshold == 0.1, Window_2sided_kb]),
    "===========================================================",
    "",
    "Paste-ready paragraph for the Results chapter:",
    "-----------------------------------------------",
    sprintf("  Genome-wide linkage disequilibrium (LD) decay was estimated from"),
    sprintf("  pairwise r^2 values computed in PLINK for all SNP pairs within"),
    sprintf("  %.1f Mb across the seven barley chromosomes in the 290 wild barley", MAX_DIST_BP / 1e6),
    sprintf("  accessions. Mean r^2 within 1-kb distance bins was smoothed using"),
    sprintf("  a pair-count-weighted LOESS curve (span = %.2f). The smoothed", LOESS_SPAN),
    sprintf("  genome-wide r^2 decayed from ~%.2f at short distances toward an", r2_max),
    sprintf("  asymptote of ~%.2f at 2.5 Mb, crossing r^2 = 0.2 at ~%.0f kb.", r2_baseline, cross02),
    sprintf("  This distance was adopted as the half-width of the +/- %.0f kb", cross02),
    sprintf("  candidate-gene search window around lead SNPs from the GWAS.")
)
writeLines(summary_lines, file.path(RES_DIR, "ld_decay_summary_for_results_chapter.txt"))

cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
cat("\n--- Done. Outputs in: ", RES_DIR, " ---\n", sep = "")
