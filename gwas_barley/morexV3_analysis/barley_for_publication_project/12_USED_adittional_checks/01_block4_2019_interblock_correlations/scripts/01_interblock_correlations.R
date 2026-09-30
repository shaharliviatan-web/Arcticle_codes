# =============================================================================
# 01_interblock_correlations.R
# -----------------------------------------------------------------------------
# Check 01: inter-block Pearson correlations of the four NIR grain traits
# (protein, starch, beta-glucan, fiber; grain content, %) in the wild-barley
# common garden, per season. Documents why block 4 of the 2019-20 season was
# excluded from the phenotype analysis (00_THIN_.../01_phenotypic_analysis_
# no_GxE_v3_streamlined.R, section 2, "Filter 2").
#
# Input  (read-only): 00_THIN_Generate_Plots_For_Publication/all years barley.csv
# Filters (same as the 00_THIN script, section 2):
#   * drop cultivated checks: short_Tag starting with c / C / M
#   * drop site 04 (Hachola) accessions: short_Tag starting with HS04
#   * block 4 of 2019 is KEPT (it is the object of the check)
#   * the "+8.85 starch in 2021" adjustment is NOT applied and values are NOT
#     centred: a constant within a season does not change within-season
#     correlations.
# Analyses:
#   raw : raw values
#   sd3 : after the project's per-trait +/-3 SD outlier removal
#         (clean_trait_local(): values centred within season, pooled across
#         seasons, NA removed, values outside mean +/- 3 SD dropped).
#         Centring here uses the retained data (block 4 of 2019 included).
#
# Run: /usr/bin/Rscript scripts/01_interblock_correlations.R   (from the check dir)
# Author: Shahar Liviatan
# =============================================================================

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- Paths -----------------------------------------------------------------
PROJ      <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
INPUT_CSV <- file.path(PROJ, "00_THIN_Generate_Plots_For_Publication", "all years barley.csv")
CHECK_DIR <- file.path(PROJ, "12_USED_adittional_checks", "01_block4_2019_interblock_correlations")
OUT_TAB   <- file.path(CHECK_DIR, "results", "tables")
OUT_FIG   <- file.path(CHECK_DIR, "results", "figures")
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)

TRAITS       <- c("ProteinAsis.", "StarchAsis.", "BetaglucansAsis.", "FiberAsis.")
TRAIT_LABEL  <- c(ProteinAsis. = "Protein", StarchAsis. = "Starch",
                  BetaglucansAsis. = "Beta-glucan", FiberAsis. = "Fiber")
SEASON_LABEL <- c("2019" = "2019-20", "2020" = "2020-21", "2021" = "2021-22")

cat("R version:", R.version.string, "\n")
cat("Run date :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ---- 1. Load + filter --------------------------------------------------------
d <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
cat("Rows read:", nrow(d), "\n")
d <- d[!grepl("^[cCM]", d$short_Tag), ]   # cultivated checks (Morex, Clipper)
d <- d[!grepl("^HS04",  d$short_Tag), ]   # site 04 (Hachola)
cat("Rows after filters:", nrow(d), " | unique accessions:",
    length(unique(d$short_Tag)), "\n")
stopifnot(!any(duplicated(d[, c("short_Tag", "season", "Block_no")])))  # 1 plant / accession / block
cat("\nRecords per season x block:\n"); print(table(d$season, d$Block_no))

# ---- 2. +/-3 SD cleaning (clean_trait_local() logic) ---------------------------
# Returns a copy of d with trait values outside mean +/- 3 SD (on season-centred,
# season-pooled values) set to NA.
sd3_clean <- function(d, trait) {
  v  <- d[[trait]]
  vc <- v - ave(v, d$season, FUN = function(x) mean(x, na.rm = TRUE))
  m  <- mean(vc, na.rm = TRUE); s <- sd(vc, na.rm = TRUE)
  keep <- !is.na(vc) & vc >= m - 3 * s & vc <= m + 3 * s
  v[!keep] <- NA
  v
}

# ---- 3. Pairwise block correlations -------------------------------------------
pair_cor <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  n  <- sum(ok)
  if (n < 4) return(c(n = n, r = NA, ci_low = NA, ci_high = NA, p = NA))
  ct <- cor.test(x[ok], y[ok], method = "pearson")
  c(n = n, r = unname(ct$estimate), ci_low = ct$conf.int[1],
    ci_high = ct$conf.int[2], p = ct$p.value)
}

res  <- list()
sd3_log <- list()
for (tr in TRAITS) {
  vals <- list(raw = d[[tr]], sd3 = sd3_clean(d, tr))
  n_out <- sum(!is.na(d[[tr]])) - sum(!is.na(vals$sd3))
  sd3_log[[tr]] <- data.frame(trait = TRAIT_LABEL[[tr]],
                              n_nonNA = sum(!is.na(d[[tr]])),
                              n_removed_sd3 = n_out)
  for (an in names(vals)) {
    for (se in sort(unique(d$season))) {
      ds  <- d[d$season == se, ]
      v   <- vals[[an]][d$season == se]
      wide <- tapply(v, list(ds$short_Tag, ds$Block_no), identity)  # accession x block
      blocks <- sort(as.integer(colnames(wide)))
      for (i in seq_along(blocks)) for (j in seq_along(blocks)) if (i < j) {
        a <- blocks[i]; b <- blocks[j]
        pc <- pair_cor(wide[, as.character(a)], wide[, as.character(b)])
        res[[length(res) + 1]] <- data.frame(
          season = se, trait = TRAIT_LABEL[[tr]], analysis = an,
          block_a = a, block_b = b, n = as.integer(pc["n"]),
          r = pc["r"], ci_low = pc["ci_low"], ci_high = pc["ci_high"], p = pc["p"],
          row.names = NULL)
      }
    }
  }
}
all_cor <- do.call(rbind, res)
cat("\nValues removed by +/-3 SD cleaning (all seasons pooled):\n")
print(do.call(rbind, sd3_log), row.names = FALSE)

fmt <- function(df) {
  df$r <- round(df$r, 3); df$ci_low <- round(df$ci_low, 3)
  df$ci_high <- round(df$ci_high, 3); df$p <- signif(df$p, 3); df
}
write.table(fmt(all_cor), file.path(OUT_TAB, "interblock_correlations_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. Summary: block-4 pairs vs pairs among blocks 1-3 ------------------------
summ <- list()
for (se in sort(unique(all_cor$season))) for (tr in TRAIT_LABEL) for (an in c("raw", "sd3")) {
  s  <- all_cor[all_cor$season == se & all_cor$trait == tr & all_cor$analysis == an, ]
  b4 <- s$r[s$block_a == 4 | s$block_b == 4]
  b13 <- s$r[s$block_a <= 3 & s$block_b <= 3]
  summ[[length(summ) + 1]] <- data.frame(
    season = se, trait = tr, analysis = an,
    n_blocks = length(unique(c(s$block_a, s$block_b))),
    n_pairs_b4 = length(b4),
    b4_mean_r = mean(b4), b4_min_r = min(b4), b4_max_r = max(b4),
    n_pairs_b1to3 = length(b13),
    b1to3_mean_r = mean(b13), b1to3_min_r = min(b13), b1to3_max_r = max(b13),
    n_range = paste0(min(s$n), "-", max(s$n)))
}
summary_df <- do.call(rbind, summ)
num_cols <- grep("_r$", names(summary_df))
summary_df[num_cols] <- lapply(summary_df[num_cols], round, 3)
write.table(summary_df, file.path(OUT_TAB, "interblock_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSummary (block-4 pairs vs pairs among blocks 1-3):\n")
print(summary_df, row.names = FALSE)

# ---- 5. Online Resource table: 2019-20, raw ------------------------------------
or <- all_cor[all_cor$season == 2019 & all_cor$analysis == "raw", ]
or$pair <- paste0(or$block_a, "-", or$block_b)
or_tab <- data.frame(block_pair = unique(or$pair))
for (tr in TRAIT_LABEL) {
  s <- or[or$trait == tr, ]
  s <- s[match(or_tab$block_pair, s$pair), ]
  or_tab[[paste0(tr, "_r")]] <- sprintf("%.3f", s$r)
  or_tab[[paste0(tr, "_n")]] <- s$n
}
write.table(or_tab, file.path(OUT_TAB, "Table_OnlineResource_2019_block_correlations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nOnline Resource table (2019-20, raw):\n"); print(or_tab, row.names = FALSE)

# ---- 6. Plain-prose key numbers -------------------------------------------------
f3 <- function(x) sprintf("%.3f", x)
sm <- function(se, tr, an) summary_df[summary_df$season == se & summary_df$trait == tr &
                                        summary_df$analysis == an, ]
lines <- c(
  "Check 01 - inter-block correlations of the NIR grain traits (key numbers)",
  paste0("Generated by scripts/01_interblock_correlations.R on ", format(Sys.Date()),
         " (", R.version.string, ")."),
  "Pearson r between blocks, across wild accessions measured in both blocks (one plant per accession per block).",
  "Filters: cultivated checks and site 04 (Hachola) removed; 290 accessions per block; block 4 of 2019-20 retained.",
  "")
for (an in c("raw", "sd3")) {
  lines <- c(lines, if (an == "raw") "== 2019-20, raw values ==" else
    "== 2019-20, after +/-3 SD outlier removal (sensitivity) ==")
  for (tr in TRAIT_LABEL) {
    x <- sm(2019, tr, an)
    lines <- c(lines, sprintf(
      "%s: block 4 vs blocks 1-3, r = %s-%s (mean %s); among blocks 1-3, r = %s-%s (mean %s); n = %s accessions per pair.",
      tr, f3(x$b4_min_r), f3(x$b4_max_r), f3(x$b4_mean_r),
      f3(x$b1to3_min_r), f3(x$b1to3_max_r), f3(x$b1to3_mean_r), x$n_range))
  }
  lines <- c(lines, "")
}
for (se in c(2020, 2021)) {
  lines <- c(lines, sprintf("== Reference: %s, raw values (%d blocks) ==",
                            SEASON_LABEL[[as.character(se)]], sm(se, "Protein", "raw")$n_blocks))
  for (tr in TRAIT_LABEL) {
    x <- sm(se, tr, "raw")
    lines <- c(lines, sprintf(
      "%s: block-4 pairs r = %s-%s (mean %s); among blocks 1-3, r = %s-%s (mean %s); n = %s.",
      tr, f3(x$b4_min_r), f3(x$b4_max_r), f3(x$b4_mean_r),
      f3(x$b1to3_min_r), f3(x$b1to3_max_r), f3(x$b1to3_mean_r), x$n_range))
  }
  lines <- c(lines, "")
}
all_2019 <- all_cor[all_cor$season == 2019, ]
for (an in c("raw", "sd3")) {
  pp <- all_2019$p[all_2019$analysis == an]
  lines <- c(lines, sprintf("2019-20 (%s): two-sided P of the 24 pairwise correlations ranged from %s to %s.",
                            an, signif(min(pp), 3), signif(max(pp), 3)))
}
bg <- sm(2019, "Beta-glucan", "raw")
lines <- c(lines, "",
  "Sanity check vs older drafts (beta-glucan 2019-20: block 4 r = 0.328-0.391; blocks 1-3 r = 0.670-0.682):",
  sprintf("this analysis gives block 4 r = %s-%s and blocks 1-3 r = %s-%s (raw).",
          f3(bg$b4_min_r), f3(bg$b4_max_r), f3(bg$b1to3_min_r), f3(bg$b1to3_max_r)))
writeLines(lines, file.path(OUT_TAB, "results_numbers.txt"))
cat("\n"); cat(lines, sep = "\n")

# ---- 7. Figure: per-trait block x block heatmaps, three seasons ----------------
# Lower triangle, raw values. Rows = traits, columns = seasons; panels a-c = seasons.
FONT <- "Liberation Sans"   # metric-compatible with Arial (Arial not installed)
hm <- all_cor[all_cor$analysis == "raw", ]
hm$season_lab <- factor(SEASON_LABEL[as.character(hm$season)], levels = SEASON_LABEL)
hm$trait <- factor(hm$trait, levels = TRAIT_LABEL)
hm$bx <- factor(hm$block_a, levels = 1:5)   # x = lower-numbered block
hm$by <- factor(hm$block_b, levels = 1:5)   # y = higher-numbered block

make_panel <- function(se_lab) {
  s <- hm[hm$season_lab == se_lab, ]
  nb <- max(s$block_b)
  s$bx <- factor(s$block_a, levels = 1:(nb - 1))
  s$by <- factor(s$block_b, levels = nb:2)
  ggplot(s, aes(bx, by, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", r)), family = FONT, size = 8 / .pt,
              colour = ifelse(s$r > 0.55, "white", "black")) +
    facet_wrap(~ trait, nrow = 1) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08306B", limits = c(0, 1),
                        breaks = seq(0, 1, 0.25), name = "Pearson r") +
    labs(x = "Block", y = "Block") +
    coord_equal() +
    theme_minimal(base_size = 9, base_family = FONT) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8, colour = "black"),
          strip.text = element_text(size = 9, face = "bold"),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          plot.tag = element_text(size = 12, face = "bold"))
}
p <- (make_panel("2019-20") / make_panel("2020-21") / make_panel("2021-22")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right")
ggsave(file.path(OUT_FIG, "Fig_interblock_correlation_heatmaps.png"), p,
       width = 174, height = 170, units = "mm", dpi = 600, bg = "white",
       device = ragg::agg_png)
cat("\nFigure written:", file.path(OUT_FIG, "Fig_interblock_correlation_heatmaps.png"), "\n")
cat("Done.\n")
