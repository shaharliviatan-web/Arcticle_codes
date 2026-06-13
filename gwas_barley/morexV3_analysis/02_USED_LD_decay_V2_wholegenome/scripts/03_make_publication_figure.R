# =============================================================================
# 03_make_publication_figure.R
#
# Publication-ready LD decay figures, built on the outputs of
# 02_compute_decay_profile.R.
#
# Produces:
#   results/figure_ld_decay_genomewide.png       (300 dpi raster)
#   results/figure_ld_decay_genomewide.pdf       (cairo vector, 174 mm)
#   results/figure_ld_decay_perchromosome.pdf    (faceted supplementary)
#   results/figure_ld_decay_perchromosome.png    (raster supplementary)
#
# All annotations (cross-point at r^2 = 0.2, etc.) are read from
# ld_decay_crossings.tsv / ld_decay_loess_fit.rds — nothing is hardcoded.
# =============================================================================

# Force temp files off the system disk. /tmp is small and gets purged.
dir.create("/mnt/data/shahar/.tmp", showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp",
           TMP    = "/mnt/data/shahar/.tmp",
           TEMP   = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

# --- 1. Paths -----------------------------------------------------------------
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE))
if (length(script_path) == 0) script_path <- "03_make_publication_figure.R"
SCRIPT_DIR  <- normalizePath(dirname(script_path))
PROJECT_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
RES_DIR     <- file.path(PROJECT_DIR, "results")

# --- 2. Read inputs -----------------------------------------------------------
gw        <- fread(file.path(RES_DIR, "ld_decay_profile_genomewide.tsv"))
per_chrom <- fread(file.path(RES_DIR, "ld_decay_profile_per_chromosome.tsv"))
crossings <- fread(file.path(RES_DIR, "ld_decay_crossings.tsv"))
pred_grid <- fread(file.path(RES_DIR, "ld_decay_loess_predictions.tsv"))
fit_obj   <- readRDS(file.path(RES_DIR, "ld_decay_loess_fit.rds"))

# Genome-wide bins use Bin (bp) -> Distance_kb
gw[, Distance_kb := Bin / 1000]
per_chrom[, Distance_kb := Bin / 1000]

THRESHOLD     <- 0.2
cross02_kb    <- crossings[R2_threshold == THRESHOLD, Distance_kb]
loess_span    <- fit_obj$span
r2_baseline   <- fit_obj$r2_baseline
half_decay_kb <- fit_obj$half_decay_kb
n_samples     <- fit_obj$n
max_kb        <- max(gw$Distance_kb, na.rm = TRUE)

annot_label <- sprintf("r^2 = %.1f at ~%.0f kb", THRESHOLD, cross02_kb)
cat(sprintf("Annotating cross-point at %.1f kb (r^2 = %.2f)\n", cross02_kb, THRESHOLD))

# --- 3. Common theme ----------------------------------------------------------
publication_theme <- theme_classic(base_size = 15) +
    theme(
        axis.title       = element_text(size = 15),
        axis.text        = element_text(size = 13),
        plot.title       = element_text(size = 16, hjust = 0.5),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text       = element_text(size = 13, face = "bold"),
        plot.margin      = margin(5, 8, 5, 5, "mm")
    )

# --- 4. Main genome-wide figure ----------------------------------------------
main_fig <- ggplot() +
    geom_point(data = gw,
               aes(x = Distance_kb, y = Avg_R2),
               color = "grey45", alpha = 0.35, size = 0.7) +
    geom_line(data = pred_grid,
              aes(x = Distance_kb, y = R2_fit),
              color = "darkred", linewidth = 1.0) +
    geom_hline(yintercept = THRESHOLD,
               linetype = "dashed", color = "steelblue", linewidth = 0.6) +
    geom_vline(xintercept = cross02_kb,
               linetype = "dotted", color = "steelblue", linewidth = 0.5) +
    scale_x_continuous(
        name   = "Physical distance (kb)",
        limits = c(0, max_kb), expand = c(0, 0)
    ) +
    scale_y_continuous(
        name   = expression(Average ~ italic(r)^2),
        limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))
    ) +
    publication_theme

ggsave(file.path(RES_DIR, "figure_ld_decay_genomewide.pdf"),
       plot = main_fig, device = cairo_pdf,
       width = 174, height = 110, units = "mm")
ggsave(file.path(RES_DIR, "figure_ld_decay_genomewide.png"),
       plot = main_fig, width = 174, height = 110, units = "mm", dpi = 300)

# --- 5. Supplementary: faceted per-chromosome --------------------------------
per_chrom_fig <- ggplot(per_chrom, aes(x = Distance_kb, y = Avg_R2)) +
    geom_point(color = "grey45", alpha = 0.35, size = 0.5) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.5,
                se = FALSE, color = "darkred", linewidth = 0.7) +
    geom_hline(yintercept = THRESHOLD,
               linetype = "dashed", color = "steelblue", linewidth = 0.5) +
    facet_wrap(~ Chromosome, ncol = 4, scales = "free_y") +
    scale_x_continuous(
        name   = "Physical distance (kb)",
        limits = c(0, max_kb), expand = c(0, 0),
        breaks = c(0, 1000, 2000)
    ) +
    scale_y_continuous(
        name   = expression(Average ~ italic(r)^2),
        limits = c(0, NA), expand = expansion(mult = c(0.01, 0.05))
    ) +
    labs(title = "LD decay by chromosome (wild barley, n = 290)") +
    publication_theme

ggsave(file.path(RES_DIR, "figure_ld_decay_perchromosome.pdf"),
       plot = per_chrom_fig, device = cairo_pdf,
       width = 260, height = 150, units = "mm")
ggsave(file.path(RES_DIR, "figure_ld_decay_perchromosome.png"),
       plot = per_chrom_fig, width = 260, height = 150, units = "mm", dpi = 300)

cat(sprintf("Saved figures to: %s\n", RES_DIR))

# =============================================================================
# Suggested figure legend (paste into the manuscript):
#
#   Figure X. Genome-wide linkage disequilibrium (LD) decay in 290 wild barley
#   accessions, computed from pairwise r^2 across all seven chromosomes.
#   Grey points: mean r^2 in 1-kb physical distance bins (0 - 2.5 Mb).
#   Red curve: pair-count-weighted LOESS smoother (span = 0.10).
#   Dashed blue line: r^2 = 0.2 threshold; the smoothed curve crosses this
#   threshold at ~[cross02_kb] kb (dotted vertical line), which was adopted
#   as the half-width of the candidate-gene search window around lead SNPs
#   from the GWAS.
# =============================================================================
