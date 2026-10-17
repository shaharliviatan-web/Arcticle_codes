# =============================================================================
# 00_common.R -- shared setup of check 02 (raw phenotype checks for the supplement).
# Sourced by 01_boxplots.R, 02_correlations.R, 03_distributions.R.
#
# Input (read-only):
#   00_THIN_Generate_Plots_For_Publication/all years barley.csv          raw plant-level data
#   00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables/BLUP_10_traits.csv   genotype BLUPs
# Panel = the analysis panel of step 00 (01_phenotypic_analysis_no_GxE_v3_streamlined.R, section 2):
#   cultivated checks (short_Tag starting with c / C / M) removed; site 04 (Hachola) removed; block 4 of
#   2019-20 removed (low inter-block concordance; check 01) -- EXCEPT in the per-block boxplots, which show it,
#   marked as excluded. Values are the raw values as the analyses use them, before centering and outlier removal:
#   the 2021-22 starch values carry the +8.85 of step 00 (user decision ESM0 Q12a, 2026-10-01: "fix the data with
#   the +8.85 and use this version for the supplement"; the correction is NOT reported in any supplementary file or
#   in the manuscript; it cancels under within-season centering).
# The +/-3 SD outlier rule of step 00 (clean_trait_local(): values centered within season, pooled over seasons,
#   NA removed, kept if within mean +/- 3 SD) is shown, not applied: on the raw scale its limits are
#   season mean + mean(centered) +/- 3 SD(centered), per trait and season.
# Figures follow TAG (10_USED_Paper_writing/TAG_requirements.md): 174 mm wide, Liberation Sans (Arial metrics)
#   9 pt, panels a, b, c (bold lower case), no titles inside the image; PDF (cairo, vector) + PNG 600 dpi.
# Created 2026-10-01 (user request; supplementary session). Not cited in the manuscript: the outputs go to the
#   unreferenced supplement (10_USED_Paper_writing/new_publishing_paper/supplementary/unreferenced/).
# =============================================================================
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(dplyr); library(tidyr) })

PROJ      <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
INPUT_CSV <- file.path(PROJ, "00_THIN_Generate_Plots_For_Publication", "all years barley.csv")
BLUP_CSV  <- file.path(PROJ, "00_THIN_Generate_Plots_For_Publication", "outputs", "subsection_1", "tables", "BLUP_10_traits.csv")
CHECK_DIR <- file.path(PROJ, "12_USED_adittional_checks_for_supplementary", "02_raw_phenotype_checks")
OUT_TAB   <- file.path(CHECK_DIR, "results", "tables")
OUT_FIG   <- file.path(CHECK_DIR, "results", "figures")

FONT <- "Liberation Sans"; PT <- 9; DPI <- 600; W_MM <- 174

TRAITS <- data.frame(
  col   = c("ProteinAsis.", "StarchAsis.", "BetaglucansAsis.", "FiberAsis.", "flowering_time", "Tiller.Length",
            "Spike.Length", "Number.Of.Tillers", "Grain_weight", "number_of_grains"),
  short = c("Protein", "Starch", "β-glucan", "Fiber", "Flowering time", "Plant height", "Spike length",
            "Number of tillers", "Grain weight", "Number of grains"),
  unit  = c("%", "%", "%", "%", "days", "", "cm", "count", "", "count"),   # units as in Potapenko et al. 2026a (ESM0 Q11b)
  blup  = c("ProteinCenterd", "StarchCenterd", "BetaglucansCenterd", "FiberCenterd", "flowering_timeCenterd",
            "Tiller.LengthCenterd", "Spike.LengthCenterd", "Number.Of.TillersCenterd", "Grain_weightCenterd",
            "number_of_grainsCenterd"),
  stringsAsFactors = FALSE)
TRAITS$label <- ifelse(TRAITS$unit %in% c("count", ""), TRAITS$short, sprintf("%s (%s)", TRAITS$short, TRAITS$unit))
NUTRI <- TRAITS[1:4, ]
SEASON_LABEL <- c("2019" = "2019–20", "2020" = "2020–21", "2021" = "2021–22")
SEASON_COL   <- c("2019–20" = "#1B9E77", "2020–21" = "#D95F02", "2021–22" = "#7570B3")   # + linetype (accessibility)
SEASON_LTY   <- c("2019–20" = "solid", "2020–21" = "dashed", "2021–22" = "dotted")

load_raw <- function(keep_block4_2019 = FALSE) {
  d <- read.csv(INPUT_CSV, stringsAsFactors = FALSE, check.names = TRUE)
  d$StarchAsis.[d$season == 2021] <- d$StarchAsis.[d$season == 2021] + 8.85   # as step 00, Filter 1 (see header)
  d <- d[!grepl("^[cCM]", d$short_Tag), ]                 # cultivated checks (Morex, Clipper)
  d <- d[!grepl("^HS04", d$short_Tag), ]                  # site 04 (Hachola)
  d$block4_2019 <- d$season == 2019 & d$Block_no == 4
  if (!keep_block4_2019) d <- d[!d$block4_2019, ]
  d$acc <- paste0(substr(d$short_Tag, 3, 4), "_", substr(d$short_Tag, 5, 6))
  d$season_lab <- factor(SEASON_LABEL[as.character(d$season)], levels = SEASON_LABEL)
  d
}

# +/-3 SD limits of step 00 on the raw scale (computed on the analysis panel, block 4 of 2019-20 removed)
sd3_limits <- function(d) {
  do.call(rbind, lapply(seq_len(nrow(TRAITS)), function(i) {
    v <- TRAITS$col[i]
    sm <- tapply(d[[v]], d$season_lab, mean, na.rm = TRUE)
    cen <- d[[v]] - sm[as.character(d$season_lab)]
    cen <- cen[!is.na(cen)]
    m <- mean(cen); s <- sd(cen)
    out <- data.frame(trait = TRAITS$short[i], season = names(sm), season_mean = as.numeric(sm),
                      lower = as.numeric(sm) + m - 3 * s, upper = as.numeric(sm) + m + 3 * s)
    out[!is.na(out$season_mean), ]
  }))
}

theme_tag <- function(base_size = PT) {
  theme_bw(base_size = base_size, base_family = FONT) +
    theme(strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.4),
          strip.text = element_text(face = "bold", size = base_size),
          axis.text = element_text(size = base_size - 1, colour = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
          legend.position = "bottom", legend.text = element_text(size = base_size - 1),
          legend.title = element_text(size = base_size - 1, face = "bold"),
          plot.tag = element_text(size = base_size + 1, face = "bold"),
          plot.margin = margin(2, 3, 1, 2))
}

save_fig <- function(p, name, h_mm) {
  stopifnot(h_mm <= 234)                                   # TAG maximum height
  cairo_pdf(file.path(OUT_FIG, paste0(name, ".pdf")), width = W_MM / 25.4, height = h_mm / 25.4, family = FONT)
  print(p); invisible(dev.off())
  ragg::agg_png(file.path(OUT_FIG, paste0(name, ".png")), width = W_MM, height = h_mm, units = "mm", res = DPI,
                background = "white")
  print(p); invisible(dev.off())
  cat(sprintf("  figure %s (%d x %d mm; PDF + PNG %d dpi)\n", name, W_MM, h_mm, DPI))
}

write_tsv <- function(x, name) {
  write.table(x, file.path(OUT_TAB, paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  cat(sprintf("  table  %s (%d rows)\n", name, nrow(x)))
}
