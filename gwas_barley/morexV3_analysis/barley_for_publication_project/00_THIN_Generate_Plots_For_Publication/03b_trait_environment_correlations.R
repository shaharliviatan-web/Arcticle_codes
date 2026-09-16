# =============================================================================
# Script C2 — STAGE 2 : Trait <-> environment correlations (site level)
# =============================================================================
# Standalone SECOND stage of the revised environmental analysis.
#
# Scope (deliberately limited):
#   * Pearson correlations between the four site-mean nutritional-trait BLUPs
#     and the EIGHT environmental predictors retained after the independent
#     stage-1 collinearity assessment.
#   * Significance from RAW two-sided P values (NO multiple-comparison / FDR
#     correction), across all 4 x 8 = 32 tests.
#   * Spearman correlations as a sensitivity check.
#   * Quality-control scatterplots for raw-P-significant associations.
#   * NO multiple regression, NO LMM, NO predictor selection, NO VIF, and NO
#     causal interpretation in this stage.
#
# Stage-1 result (independent of any trait data): two |r|>0.8 pairs ->
#   March temperature ~ Elevation (r = -0.951)   ->  keep March temperature
#   Sand ~ Silt              (r = -0.934)         ->  keep Sand
# (Biological-interpretability choices supplied by the user.)
# Therefore Elevation and Silt are EXCLUDED here, leaving 8 predictors.
#
# Site = unit of analysis (29 sites; HS04/Hachola excluded -- no phenotype data).
#
# Inputs:
#   BLUP_10_traits.csv
#   merged_Evgeny_env_fitness_environmental_30_sites.csv
#
# Output root: outputs/subsection_1/C2_trait_environment_correlations/
#   tables/   figures/   reports/
#
# Author: Shahar Liviatan
# =============================================================================

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# ---- Inputs -----------------------------------------------------------------

BLUP_FILE <- "BLUP_10_traits.csv"
ENV_FILE  <- "merged_Evgeny_env_fitness_environmental_30_sites.csv"
if (!file.exists(BLUP_FILE)) {
  fb <- file.path("outputs", "subsection_1", "tables", "BLUP_10_traits.csv")
  if (file.exists(fb)) BLUP_FILE <- fb
}
stopifnot(file.exists(BLUP_FILE), file.exists(ENV_FILE))

# ---- Output structure -------------------------------------------------------

OUT_ROOT <- file.path("outputs", "subsection_1", "C2_trait_environment_correlations")
DIR_TAB  <- file.path(OUT_ROOT, "tables")
DIR_FIG  <- file.path(OUT_ROOT, "figures")
DIR_REP  <- file.path(OUT_ROOT, "reports")
DIR_SCAT <- file.path(DIR_FIG, "significant_scatterplots")
for (d in c(DIR_TAB, DIR_FIG, DIR_REP, DIR_SCAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Publication theme + save helper (identical to scripts 01/02) -----------

theme_pub <- function(base_size = 16) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey40", linewidth = 0.4),
      strip.text       = element_text(face = "bold", size = base_size),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1)
    )
}

save_plot <- function(plot, name, w = 7, h = 5, dir = DIR_FIG) {
  ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), plot, width = w, height = h,
         dpi = 600, bg = "white")
  invisible(plot)
}

global_sig <- function(q) dplyr::case_when(
  is.na(q)  ~ "",
  q < 0.001 ~ "***",
  q < 0.01  ~ "**",
  q < 0.05  ~ "*",
  TRUE      ~ "")

# ---- Trait + retained-predictor definitions ---------------------------------

TRAITS_NUTRI <- c("ProteinCenterd", "StarchCenterd",
                  "BetaglucansCenterd", "FiberCenterd")
TRAIT_LABELS <- c(ProteinCenterd = "Protein", StarchCenterd = "Starch",
                  BetaglucansCenterd = "β-glucan", FiberCenterd = "Fiber")
TRAIT_ORDER  <- c("Protein", "Starch", "β-glucan", "Fiber")

# 8 retained predictors (Silt + Elevation excluded by stage-1 collinearity).
FEATURE_MAP <- c(
  "pH"              = "pH",
  "Electrical_Cond" = "Electrical.Conductivity..mS.cm.",
  "Clay_pct"        = "Clay....",
  "Sand_pct"        = "Sand....",
  "March_Temp"      = "March..temp",
  "Precipitation"   = "Precipitation",
  "Organic_Carbon"  = "Organic.Carbon..mg.Kg.",
  "Total_N"         = "Total.N..mg.Kg."
)
PREDICTOR_LABELS <- c(
  "pH"              = "pH",
  "Electrical_Cond" = "Electrical conductivity (mS/cm)",
  "Clay_pct"        = "Clay (%)",
  "Sand_pct"        = "Sand (%)",
  "March_Temp"      = "March temperature (°C)",
  "Precipitation"   = "Precipitation (mm)",
  "Organic_Carbon"  = "Organic carbon (mg/kg)",
  "Total_N"         = "Total N (mg/kg)"
)
PREDICTOR_UNITS <- c(
  "pH" = "pH units (unitless)", "Electrical_Cond" = "mS/cm",
  "Clay_pct" = "%", "Sand_pct" = "%", "March_Temp" = "°C",
  "Precipitation" = "mm", "Organic_Carbon" = "mg/kg", "Total_N" = "mg/kg")
PREDICTOR_REASON <- c(
  "pH"              = "No strong collinearity in stage 1 (max |r| <= 0.8); retained.",
  "Electrical_Cond" = "No strong collinearity in stage 1 (max |r| <= 0.8); retained.",
  "Clay_pct"        = "Max |r| with another predictor = 0.75 (Sand), below 0.8; retained.",
  "Sand_pct"        = "Retained over Silt (Sand-Silt r = -0.93) on biological interpretability.",
  "March_Temp"      = "Retained over Elevation (March temp-Elevation r = -0.95) on biological interpretability.",
  "Precipitation"   = "No strong collinearity in stage 1 (max |r| <= 0.8); retained.",
  "Organic_Carbon"  = "No strong collinearity in stage 1 (max |r| <= 0.8); retained.",
  "Total_N"         = "No strong collinearity in stage 1 (max |r| <= 0.8); retained.")
PREDICTORS    <- names(FEATURE_MAP)
EXCLUDE_SITES <- c("HS04")     # Hachola/Hula: no phenotype data

# -----------------------------------------------------------------------------
# 1. BUILD SITE-LEVEL DATASET
# -----------------------------------------------------------------------------

message("--- Building site-level dataset ---")

blup_raw <- read.csv(BLUP_FILE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("Genotype" %in% colnames(blup_raw))
miss_tr <- setdiff(TRAITS_NUTRI, colnames(blup_raw))
if (length(miss_tr)) stop("BLUP file missing traits: ", paste(miss_tr, collapse = ", "))

n_geno_rows <- nrow(blup_raw)

blup_df <- blup_raw %>%
  dplyr::select(Genotype, all_of(TRAITS_NUTRI)) %>%
  mutate(site_id = paste0("HS", substr(Genotype, 1, 2))) %>%
  filter(!site_id %in% EXCLUDE_SITES)

site_means <- blup_df %>%
  group_by(site_id) %>%
  summarise(n_genotypes = dplyr::n(),
            across(all_of(TRAITS_NUTRI), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

env_raw <- read.csv(ENV_FILE, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8") %>% mutate(site_id = as.character(site_id))
n_env_sites_before <- length(unique(env_raw$site_id))

miss_cols <- setdiff(unname(FEATURE_MAP), colnames(env_raw))
if (length(miss_cols)) stop("Env file missing columns: ", paste(miss_cols, collapse = ", "))

env_site <- env_raw[, c("site_id", unname(FEATURE_MAP))]
names(env_site) <- c("site_id", PREDICTORS)
env_site <- env_site %>%
  mutate(across(all_of(PREDICTORS), ~ suppressWarnings(as.numeric(.)))) %>%
  filter(!site_id %in% EXCLUDE_SITES) %>%
  distinct(site_id, .keep_all = TRUE)

site_dat <- site_means %>% inner_join(env_site, by = "site_id") %>% arrange(site_id)
n_sites  <- nrow(site_dat)

# Missing-value audit.
miss_trait <- sapply(TRAITS_NUTRI, function(t) sum(!is.finite(site_dat[[t]])))
miss_pred  <- sapply(PREDICTORS,   function(p) sum(!is.finite(site_dat[[p]])))
any_missing <- sum(miss_trait) + sum(miss_pred) > 0

write.csv(site_dat, file.path(DIR_TAB, "C2corr_site_level_dataset.csv"),
          row.names = FALSE)

message(sprintf("    Genotype rows: %d | env sites before excl: %d | excluded: %s | matched sites: %d",
                n_geno_rows, n_env_sites_before, paste(EXCLUDE_SITES, collapse = ","), n_sites))
message(sprintf("    Genotypes/site: min %d, max %d, mean %.1f",
                min(site_dat$n_genotypes), max(site_dat$n_genotypes),
                mean(site_dat$n_genotypes)))
message(sprintf("    Missing values: traits=%d, predictors=%d (%s)",
                sum(miss_trait), sum(miss_pred),
                if (any_missing) "PRESENT" else "none"))

# ---- Predictor-definition table ---------------------------------------------
defs_df <- data.frame(
  internal_name     = PREDICTORS,
  source_csv_column = unname(FEATURE_MAP[PREDICTORS]),
  display_label     = unname(PREDICTOR_LABELS[PREDICTORS]),
  units             = unname(PREDICTOR_UNITS[PREDICTORS]),
  non_missing_sites = sapply(PREDICTORS, function(p) sum(is.finite(site_dat[[p]]))),
  reason_retained   = unname(PREDICTOR_REASON[PREDICTORS]),
  stringsAsFactors  = FALSE)
write.csv(defs_df, file.path(DIR_TAB, "C2corr_predictor_definitions.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 2. PEARSON + SPEARMAN CORRELATIONS  (4 traits x 8 predictors = 32 tests)
# -----------------------------------------------------------------------------

message("--- Trait x environment correlations (32 tests) ---")

corr_one <- function(trait, pred) {
  y <- site_dat[[trait]]; x <- site_dat[[pred]]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  cp <- suppressWarnings(cor.test(y, x, method = "pearson"))           # two-sided
  cs <- suppressWarnings(cor.test(y, x, method = "spearman", exact = FALSE))
  data.frame(Trait = unname(TRAIT_LABELS[trait]),
             trait_col = trait,
             Predictor = pred,
             PredLabel = unname(PREDICTOR_LABELS[pred]),
             n_sites   = length(x),
             pearson_r = unname(cp$estimate),
             pearson_p = cp$p.value,
             spearman_rho = unname(cs$estimate),
             spearman_p   = cs$p.value,
             stringsAsFactors = FALSE)
}

all32 <- do.call(rbind, lapply(TRAITS_NUTRI, function(tr)
  do.call(rbind, lapply(PREDICTORS, function(p) corr_one(tr, p)))))

stopifnot(nrow(all32) == 32)

# Significance based on RAW two-sided p-values -- NO multiple-comparison /
# FDR correction (per analysis decision). Stars: * p<0.05 ** p<0.01 *** p<0.001.
all32$sig_symbol   <- global_sig(all32$pearson_p)
all32$TraitLabel   <- factor(all32$Trait, levels = TRAIT_ORDER)
all32 <- all32 %>% arrange(TraitLabel, desc(abs(pearson_r)))

all32_out <- all32 %>%
  dplyr::select(Trait, Predictor, PredLabel, n_sites,
                pearson_r, pearson_p, sig_symbol)
write.csv(all32_out, file.path(DIR_TAB, "C2corr_trait_environment_all32.csv"),
          row.names = FALSE)

sig_only <- all32_out %>% filter(pearson_p < 0.05) %>% arrange(desc(abs(pearson_r)))
write.csv(sig_only, file.path(DIR_TAB, "C2corr_trait_environment_significant_only.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. SPEARMAN SENSITIVITY + single-site-influence (jackknife) QC
# -----------------------------------------------------------------------------

message("--- Spearman sensitivity + jackknife influence ---")

jackknife_r <- function(trait, pred) {
  y <- site_dat[[trait]]; x <- site_dat[[pred]]
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x)
  loo <- vapply(seq_len(n), function(i) cor(y[-i], x[-i]), numeric(1))
  list(min = min(loo), max = max(loo), site = site_dat$site_id[ok])
}

sens <- all32 %>%
  mutate(direction_pearson  = ifelse(pearson_r  >= 0, "positive", "negative"),
         direction_spearman = ifelse(spearman_rho >= 0, "positive", "negative"),
         direction_agreement = direction_pearson == direction_spearman,
         pearson_sig_raw  = pearson_p < 0.05,
         spearman_sig_raw = spearman_p < 0.05,
         conclusion_differs = (sign(pearson_r) != sign(spearman_rho)) |
                              (pearson_sig_raw != spearman_sig_raw),
         nonlinearity_flag = (abs(spearman_rho) - abs(pearson_r)) > 0.15)

# jackknife range for every test (used to flag single-site-driven results)
jk <- t(mapply(function(tr, pr) {
  z <- jackknife_r(tr, pr); c(loo_r_min = z$min, loo_r_max = z$max)
}, sens$trait_col, sens$Predictor))
sens$loo_r_min <- jk[, "loo_r_min"]
sens$loo_r_max <- jk[, "loo_r_max"]
# largest single-site swing in r away from the full-sample r
sens$max_single_site_dr <- pmax(abs(sens$pearson_r - sens$loo_r_min),
                                abs(sens$pearson_r - sens$loo_r_max))
# a significant pair is "fragile" if dropping one site can push |r| below 0.30
sens$fragile_one_site <- sens$pearson_sig_raw &
  (pmin(abs(sens$loo_r_min), abs(sens$loo_r_max)) < 0.30)

sens_out <- sens %>%
  dplyr::select(Trait, Predictor, PredLabel, n_sites,
                pearson_r, pearson_p,
                spearman_rho, spearman_p,
                direction_agreement, conclusion_differs,
                nonlinearity_flag, loo_r_min, loo_r_max,
                max_single_site_dr, fragile_one_site) %>%
  arrange(factor(Trait, levels = TRAIT_ORDER), desc(abs(pearson_r)))
write.csv(sens_out, file.path(DIR_TAB, "C2corr_pearson_spearman_sensitivity.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. MAIN FIGURE : significant-only dot plot (style of A12 Panel B)
# -----------------------------------------------------------------------------

message("--- Main dot plot (significant only) ---")

sig_plot <- all32 %>% filter(pearson_p < 0.05) %>%
  mutate(nutri = factor(Trait, levels = TRAIT_ORDER),
         xlab  = paste(PredLabel, as.integer(nutri), sep = "___"))

empties <- setdiff(TRAIT_ORDER, as.character(unique(sig_plot$Trait)))
if (length(empties)) {
  ph <- data.frame(Trait = empties, PredLabel = "", pearson_r = NA_real_,
                   sig_symbol = "", stringsAsFactors = FALSE)
  ph$nutri <- factor(ph$Trait, levels = TRAIT_ORDER)
  ph$xlab  <- paste("", as.integer(ph$nutri), sep = "___")  # "___k" -> strips to ""
  plotdat <- dplyr::bind_rows(sig_plot, ph)
} else plotdat <- sig_plot

# order x within each facet by r (placeholders sort first; harmless)
lev <- plotdat %>% arrange(nutri, pearson_r) %>% pull(xlab) %>% unique()
plotdat$xlab <- factor(plotdat$xlab, levels = lev)

p_dot <- ggplot(plotdat, aes(x = xlab, y = pearson_r)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70",
             linewidth = 0.4) +
  geom_point(data = subset(plotdat, !is.na(pearson_r)),
             size = 4, colour = "grey45") +
  geom_text(data = subset(plotdat, !is.na(pearson_r)),
            aes(label = sig_symbol), vjust = -0.9, size = 6.0,
            fontface = "bold", colour = "grey25") +
  facet_wrap(~ nutri, scales = "free_x", ncol = 2, drop = FALSE) +
  scale_x_discrete(labels = function(x) sub("___.*", "", x)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Pearson's r") +
  theme_pub() +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 15),
        axis.title.y     = element_text(size = 16),
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold", size = 16))

if (length(empties)) {
  p_dot <- p_dot +
    geom_text(data = subset(plotdat, is.na(pearson_r)),
              aes(x = xlab, label = "No significant correlations"),
              y = 0, size = 4.6, fontface = "italic", colour = "grey45")
}

save_plot(p_dot, "C2corr_significant_environmental_correlations_dotplot",
          w = 8, h = 7.6)

# -----------------------------------------------------------------------------
# 5. SUPPORTING FIGURE : full 32-test heatmap
# -----------------------------------------------------------------------------

message("--- Supporting heatmap (all 32) ---")

heat <- all32 %>%
  mutate(TraitLabel = factor(Trait, levels = rev(TRAIT_ORDER)),
         PredLabel  = factor(PredLabel,
                             levels = unname(PREDICTOR_LABELS[PREDICTORS])))
# Narrower (legend moved to the bottom, no forced-square coord_equal) with
# larger cell numbers/stars and bolder axis labels, so the panel reads well
# when placed at reduced width in the Figure-2 composite.
p_heat <- ggplot(heat, aes(x = PredLabel, y = TraitLabel, fill = pearson_r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  # emphasize cells significant at raw P < 0.05 with a black border (drawn
  # before the text layer so numbers/stars stay on top and are never covered)
  geom_tile(data = subset(heat, pearson_p < 0.05),
            fill = NA, colour = "black", linewidth = 1.3) +
  # value and stars on two lines so a larger font still fits inside each cell
  geom_text(aes(label = ifelse(sig_symbol == "", sprintf("%.2f", pearson_r),
                               sprintf("%.2f\n%s", pearson_r, sig_symbol))),
            size = 5.6, lineheight = 0.78, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson r") +
  scale_x_discrete(labels = function(x) sub("\\s*\\(.*\\)\\s*$", "", x)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = 13, barheight = 1.0,
                                title.position = "left", title.vjust = 0.9)) +
  theme_pub(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 16),
        panel.grid  = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 14),
        legend.text  = element_text(size = 12))
save_plot(p_heat, "C2corr_all32_heatmap", w = 7.6, h = 5.2)

# -----------------------------------------------------------------------------
# 5b. REDUCED heatmap : only Sand, Precipitation, Temperature
# -----------------------------------------------------------------------------
# Same style as the full 32-test heatmap above, but restricted to three
# environmental predictors on the x-axis. "March temperature" is relabelled to
# "Temperature" for this reduced panel.

message("--- Reduced heatmap (Sand, Precipitation, Temperature) ---")

SUBSET_PREDS  <- c("Sand_pct", "Precipitation", "March_Temp")
SUBSET_LABELS <- c(
  "Sand_pct"      = "Sand (%)",
  "Precipitation" = "Precipitation (mm)",
  "March_Temp"    = "Temperature (°C)"   # relabelled from "March temperature"
)

heat_sub <- all32 %>%
  filter(Predictor %in% SUBSET_PREDS) %>%
  mutate(TraitLabel = factor(Trait, levels = rev(TRAIT_ORDER)),
         PredLabel  = factor(unname(SUBSET_LABELS[Predictor]),
                             levels = unname(SUBSET_LABELS[SUBSET_PREDS])))

p_heat_sub <- ggplot(heat_sub, aes(x = PredLabel, y = TraitLabel, fill = pearson_r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = subset(heat_sub, pearson_p < 0.05),
            fill = NA, colour = "black", linewidth = 1.3) +
  geom_text(aes(label = ifelse(sig_symbol == "", sprintf("%.2f", pearson_r),
                               sprintf("%.2f\n%s", pearson_r, sig_symbol))),
            size = 5.6, lineheight = 0.78, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson r") +
  scale_x_discrete(labels = function(x) sub("\\s*\\(.*\\)\\s*$", "", x),
                   position = "top") +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 7,
                                title.position = "top", title.hjust = 0)) +
  theme_pub(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 16),
        panel.grid  = element_blank(),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 11),
        legend.text  = element_text(size = 10))
save_plot(p_heat_sub, "C2corr_sand_precip_temp_heatmap", w = 5.2, h = 5.2)

# -----------------------------------------------------------------------------
# 6. QC SCATTERPLOTS for raw-P-significant associations
# -----------------------------------------------------------------------------

message("--- QC scatterplots for significant associations ---")

make_scatter <- function(row) {
  tr <- row$trait_col; pr <- row$Predictor
  d  <- site_dat[, c("site_id", tr, pr)]
  names(d) <- c("site_id", "y", "x")
  d <- d[is.finite(d$x) & is.finite(d$y), ]
  m <- lm(y ~ x, data = d)
  d$resid_std <- rstandard(m)
  d$cook      <- cooks.distance(m)
  d$flag      <- abs(d$resid_std) > 2 | d$cook > (4 / nrow(d))
  lab <- sprintf("r = %.2f,  P = %s", row$pearson_r,
                 format.pval(row$pearson_p, digits = 2, eps = 1e-4))
  ggplot(d, aes(x, y)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey25", fill = "grey80", linewidth = 0.7) +
    geom_point(size = 2.6, colour = "grey35") +
    ggrepel::geom_text_repel(data = subset(d, flag), aes(label = site_id),
                             size = 3.4, colour = "grey25", min.segment.length = 0,
                             max.overlaps = 20) +
    labs(x = unname(PREDICTOR_LABELS[pr]),
         y = sprintf("Site-mean %s (BLUP)", unname(TRAIT_LABELS[tr])),
         caption = lab) +
    theme_pub(base_size = 14) +
    theme(plot.caption = element_text(hjust = 0, size = 12, face = "italic",
                                      colour = "grey20"))
}

scatter_paths <- character(0)
if (nrow(sig_only)) {
  sig_rows <- all32 %>% filter(pearson_p < 0.05)
  for (i in seq_len(nrow(sig_rows))) {
    g <- make_scatter(sig_rows[i, ])
    nm <- sprintf("scatter_%s__%s", sig_rows$trait_col[i], sig_rows$Predictor[i])
    save_plot(g, nm, w = 5.4, h = 4.6, dir = DIR_SCAT)
    scatter_paths <- c(scatter_paths, file.path(DIR_SCAT, paste0(nm, ".png")))
  }
}

# -----------------------------------------------------------------------------
# 7. TEXT REPORT
# -----------------------------------------------------------------------------

rep_path <- file.path(DIR_REP, "C2corr_analysis_report.txt")
sink(rep_path)
cat("Trait x environment correlations (STAGE 2)\n", strrep("=", 64), "\n", sep = "")
cat("\n[A] DATA AUDIT\n")
cat("  BLUP file:            ", BLUP_FILE, "\n", sep = "")
cat("  Env file:             ", ENV_FILE, "\n", sep = "")
cat("  Trait columns:        ", paste(TRAITS_NUTRI, collapse = ", "), "\n", sep = "")
cat("  Site-ID derivation:   site_id = paste0('HS', substr(Genotype, 1, 2))\n")
cat("  Excluded site:        ", paste(EXCLUDE_SITES, collapse = ", "),
    " (Hachola/Hula; no phenotype/BLUP data)\n", sep = "")
cat("  Genotype rows:        ", n_geno_rows, "\n", sep = "")
cat("  Env sites before excl:", n_env_sites_before, "\n", sep = "")
cat("  Matched sites:        ", n_sites, "\n", sep = "")
cat("  Genotypes/site:       min ", min(site_dat$n_genotypes), ", max ",
    max(site_dat$n_genotypes), ", mean ", round(mean(site_dat$n_genotypes), 2), "\n", sep = "")
cat("  Missing (traits/preds):", sum(miss_trait), "/", sum(miss_pred),
    if (any_missing) " (PRESENT)\n" else " (none -> all tests use n = ", "")
if (!any_missing) cat(n_sites, ")\n", sep = "")
cat("  Same n for all tests: ", all(all32$n_sites == n_sites), "\n", sep = "")

cat("\n[B] STATISTICAL AUDIT\n")
cat("  Pearson:   stats::cor.test(method='pearson'), two-sided (default)\n")
cat("  Spearman:  stats::cor.test(method='spearman', exact=FALSE)\n")
cat("  Multiple-comparison correction: NONE -- significance is RAW two-sided P.\n")
cat("  Thresholds: * P<0.05  ** P<0.01  *** P<0.001\n")
cat("  ", R.version.string, "\n", sep = "")
for (pk in c("dplyr","tidyr","tibble","ggplot2","ggrepel"))
  cat(sprintf("  pkg %-9s %s\n", pk, as.character(packageVersion(pk))))

cat("\n[C] SIGNIFICANT RESULTS (raw P < 0.05)\n")
if (nrow(sig_only)) {
  for (i in seq_len(nrow(sig_only)))
    cat(sprintf("  %-10s ~ %-28s n=%d  r=%+.3f  P=%.2e  %s\n",
                sig_only$Trait[i], sig_only$PredLabel[i], sig_only$n_sites[i],
                sig_only$pearson_r[i], sig_only$pearson_p[i],
                sig_only$sig_symbol[i]))
} else cat("  (none)\n")
no_sig_traits <- setdiff(TRAIT_ORDER, unique(sig_only$Trait))
cat("  Traits with NO significant env correlation: ",
    if (length(no_sig_traits)) paste(no_sig_traits, collapse = ", ") else "none", "\n", sep = "")

cat("\n[D] SENSITIVITY / QC\n")
dd <- sens %>% filter(conclusion_differs | nonlinearity_flag | fragile_one_site)
if (nrow(dd)) {
  for (i in seq_len(nrow(dd)))
    cat(sprintf("  %-10s ~ %-28s r=%+.2f rho=%+.2f | dir.agree=%s concl.differs=%s nonlin=%s fragile1site=%s (LOO r in [%.2f, %.2f])\n",
                dd$Trait[i], dd$PredLabel[i], dd$pearson_r[i], dd$spearman_rho[i],
                dd$direction_agreement[i], dd$conclusion_differs[i],
                dd$nonlinearity_flag[i], dd$fragile_one_site[i],
                dd$loo_r_min[i], dd$loo_r_max[i]))
} else cat("  No meaningful disagreements, nonlinearity, or single-site-driven results flagged.\n")

cat("\nNOTE: Exploratory, correlational analysis. No regression, LMM, predictor\n")
cat("selection, VIF, or causal interpretation were performed at this stage.\n")
sink()

# -----------------------------------------------------------------------------
# 8. CONSOLE SUMMARY
# -----------------------------------------------------------------------------

message("\n", strrep("=", 60))
message("STAGE 2 finished. Output root: ", OUT_ROOT)
message(strrep("=", 60))
message(sprintf("Sites: %d | tests: %d | significant (raw P<0.05): %d",
                n_sites, nrow(all32), nrow(sig_only)))
message("Significant associations:")
if (nrow(sig_only)) for (i in seq_len(nrow(sig_only)))
  message(sprintf("  %s ~ %s : r=%+.3f, P=%.3g %s",
                  sig_only$Trait[i], sig_only$PredLabel[i],
                  sig_only$pearson_r[i], sig_only$pearson_p[i],
                  sig_only$sig_symbol[i]))
message("Traits with no significant env correlation: ",
        if (length(no_sig_traits)) paste(no_sig_traits, collapse = ", ") else "none")
message("Scatterplots written: ", length(scatter_paths))
