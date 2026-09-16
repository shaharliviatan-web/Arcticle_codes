# =============================================================================
# Script B v2 (STREAMLINED) : Phenotypic analysis WITH G x E interaction
# =============================================================================
# Streamlined version of `02_phenotypic_analysis_GxE.R`. Retains every filter,
# the same per-trait QC, and the same LMM model with G x E term as v1, but
# produces ONLY the analyses required for Subsection 1 of the Results chapter:
#
#   1. Broad-sense heritability under the G x E model (CSV only, for
#      comparison with the no-GxE H^2 reported in Table 1 from Script A)
#   2. Variance partitioning, 8-trait stacked-bar plot
#   3. Reaction norms, per-trait highlights (4 representative accessions
#      picked independently in each panel)
#
# Model used   :  trait ~ (1 | short_Tag) + (1 | season_Bl) + (1 | short_Tag:season)
# Implemented  :  lme4 v1.1.37
#
# Filters applied ONCE, globally, before anything else runs:
#   1. Starch + 8.85 correction for season 2021    (post-NIR adjustment)
#   2. Drop Block 4 of season 2019                 (low inter-block concordance)
#   3. Drop control lines Morex / Clipper          (^[cCM])
#   4. Drop Site 04 (Hachola) accessions           (^04_)
#
# Per-trait local cleaning (inside every analysis function):
#   * NA exclusion for the trait under analysis
#   * +/- 3 SD outlier removal on the trait under analysis
#
# Output folder layout is preserved exactly as in v1:
#   outputs/subsection_1/pdf/
#   outputs/subsection_1/png/
#   outputs/subsection_1/tables/
#
# Author : Shahar Liviatan
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETUP -------------------------------------------------------------------
# -----------------------------------------------------------------------------

# ---- Auto-install any missing packages ------------------------------------
.required <- c("ggplot2", "lme4", "dplyr", "tidyr", "stringr")
.missing  <- .required[!vapply(.required, requireNamespace,
                               logical(1), quietly = TRUE)]
if (length(.missing)) {
  message("Installing missing packages: ", paste(.missing, collapse = ", "))
  install.packages(.missing)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- Configuration --------------------------------------------------------
INPUT_CSV <- "all years barley.csv"
OUT_ROOT  <- file.path("outputs", "subsection_1")

dir.create(file.path(OUT_ROOT, "pdf"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "png"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---- Trait definitions (8 traits used here: 4 nutri + 4 morpho) ----------
# NOTE: This 8-trait set differs from the 10-trait set in Script A.
# It omits Plant height (Tiller.Length) and Grain number, and includes
# Tillers (Number.Of.Tillers) and Spike length -- as in the original
# 02_phenotypic_analysis_GxE.R (v1).
TRAITS_RAW_NUTRI <- c("ProteinAsis.", "StarchAsis.", "BetaglucansAsis.", "FiberAsis.")
TRAITS_RAW_AGRO  <- c("flowering_time", "Number.Of.Tillers", "Grain_weight", "Spike.Length")
TRAITS_RAW_ALL   <- c(TRAITS_RAW_NUTRI, TRAITS_RAW_AGRO)
TRAITS_C_NUTRI   <- c("ProteinCenterd", "StarchCenterd",
                      "BetaglucansCenterd", "FiberCenterd")
TRAITS_C_AGRO    <- c("flowering_timeCenterd", "Number.Of.TillersCenterd",
                      "Grain_weightCenterd", "Spike.LengthCenterd")
TRAITS_C_ALL     <- c(TRAITS_C_NUTRI, TRAITS_C_AGRO)

LABEL_MAP <- c(
  ProteinAsis.            = "Protein",          ProteinCenterd            = "Protein",
  StarchAsis.             = "Starch",           StarchCenterd             = "Starch",
  BetaglucansAsis.        = "\u03B2-glucan",    BetaglucansCenterd        = "\u03B2-glucan",
  FiberAsis.              = "Fiber",            FiberCenterd              = "Fiber",
  flowering_time          = "Flowering time",   flowering_timeCenterd     = "Flowering time",
  Number.Of.Tillers       = "Tillers",          Number.Of.TillersCenterd  = "Tillers",
  Grain_weight            = "Grain weight",     Grain_weightCenterd       = "Grain weight",
  Spike.Length            = "Spike length",     Spike.LengthCenterd       = "Spike length"
)

TRAIT_CLASS_MAP <- c(
  "Protein"          = "Nutritional",
  "Starch"           = "Nutritional",
  "\u03B2-glucan"    = "Nutritional",
  "Fiber"            = "Nutritional",
  "Flowering time"   = "Agromorphological",
  "Tillers"          = "Agromorphological",
  "Grain weight"     = "Agromorphological",
  "Spike length"     = "Agromorphological"
)

# Trait factor order used for B1 (8-trait variance partitioning).
trait_order_8 <- c("Protein", "Starch", "\u03B2-glucan", "Fiber",
                   "Flowering time", "Tillers", "Grain weight", "Spike length")

# Variance-component palette (colourblind-friendly; residual is muted grey).
PAL_VARCOMP <- c(
  "Genetics"       = "#2CA02C",
  "Season & Block" = "#1F77B4",
  "G x E"          = "#FF7F0E",
  "Residual"       = "#BDBDBD"
)

# Highlight palette (Okabe-Ito-style, colourblind-friendly).
PAL_HIGHLIGHT <- c(
  "Most stable"        = "#009E73",
  "Strongest increase" = "#D55E00",
  "Strongest decrease" = "#0072B2",
  "Strongest crossover"= "#CC79A7"
)

# ---- Publication theme ----------------------------------------------------
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

# ---- Save helper: writes both PDF (vector) and PNG @ 600 dpi -------------
save_plot <- function(plot, name, w = 7, h = 5) {
  ggsave(file.path(OUT_ROOT, "pdf", paste0(name, ".pdf")),
         plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT_ROOT, "png", paste0(name, ".png")),
         plot, width = w, height = h, dpi = 600, bg = "white")
  invisible(plot)
}


# -----------------------------------------------------------------------------
# 2. PREPROCESSING (identical to v1; identical to Script A plus season->factor)
# -----------------------------------------------------------------------------

message("--- Loading and preprocessing the field-experiment data ---")

df_all <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

# Filter 1: Starch + 8.85 correction for 2021 only.
df_all$StarchAsis.[df_all$season == 2021] <-
  df_all$StarchAsis.[df_all$season == 2021] + 8.85

# Filter 2: drop Block 4 of season 2019.
df_all <- df_all[!(df_all$season == 2019 & df_all$Block_no == 4), ]

# Filter 3: drop control lines (Morex, Clipper).
df_all <- df_all[!grepl("^[cCM]", df_all$short_Tag), ]

# ID synchronisation with the VCF format ("HS0123" -> "01_23").
df_all$short_Tag <- str_replace(df_all$short_Tag, "HS", "")
df_all$short_Tag <- paste0(substr(df_all$short_Tag, 1, 2), "_",
                           substr(df_all$short_Tag, 3, 4))

# Filter 4: drop Site 04 (Hachola) accessions.
df_all <- df_all[!grepl("^04_", df_all$short_Tag), ]

# Build the season-by-block factor used as the environmental random effect.
df_all$season_Bl <- paste0(df_all$season, ".", df_all$Block_no)

df_all$short_Tag <- as.factor(df_all$short_Tag)
df_all$season_Bl <- as.factor(df_all$season_Bl)
df_all$season    <- as.factor(df_all$season)   # required for the G x E term

# Per-season mean-centering.
center_by_season <- function(df, raw_col, centered_col) {
  m <- aggregate(df[[raw_col]], list(df$season), FUN = mean, na.rm = TRUE)
  df[[centered_col]] <- df[[raw_col]]
  for (i in seq_len(nrow(m))) {
    yr <- m$Group.1[i]
    df[[centered_col]][df$season == yr] <-
      df[[raw_col]][df$season == yr] - m$x[i]
  }
  df
}
for (k in seq_along(TRAITS_RAW_ALL)) {
  df_all <- center_by_season(df_all, TRAITS_RAW_ALL[k], TRAITS_C_ALL[k])
}

message(sprintf("    Records retained: %d   |   Unique accessions: %d",
                nrow(df_all), length(unique(df_all$short_Tag))))


# -----------------------------------------------------------------------------
# 3. CORE FUNCTIONS ----------------------------------------------------------
# -----------------------------------------------------------------------------

# Per-trait local cleaning: drop NAs in `trait`, then drop +/- 3 SD outliers.
clean_trait_local <- function(df, trait) {
  d <- df[!is.na(df[[trait]]), ]
  v <- d[[trait]]
  m <- mean(v); s <- sd(v)
  d[v >= m - 3 * s & v <= m + 3 * s, , drop = FALSE]
}

# Variance partitioning under the G x E model.
# Returns: list(H2 = data.frame(pheno, heritability),
#               VarComp = data.frame(Trait, Variance, Component))
herit_calc_GxE <- function(df, traits) {
  vc_all <- list()
  h2     <- numeric(length(traits))
  for (i in seq_along(traits)) {
    tr <- traits[i]
    d  <- clean_trait_local(df, tr)
    fit <- lmer(d[[tr]] ~ (1 | short_Tag) + (1 | season_Bl) +
                  (1 | short_Tag:season),
                data = d)
    vc <- as.data.frame(VarCorr(fit))
    g   <- vc$vcov[vc$grp == "short_Tag"]
    e   <- vc$vcov[vc$grp == "season_Bl"]
    gxe <- vc$vcov[vc$grp == "short_Tag:season"]
    res <- vc$vcov[vc$grp == "Residual"]
    h2[i] <- g / (g + e + gxe + res)
    vc_all[[length(vc_all) + 1]] <- data.frame(
      Trait     = tr,
      Variance  = c(g, e, gxe, res),
      Component = c("Genetics", "Season & Block", "G x E", "Residual")
    )
    message(sprintf("    H2 (%s) = %.3f", tr, h2[i]))
  }
  list(H2       = data.frame(pheno = traits, heritability = round(h2, 3)),
       VarComp  = bind_rows(vc_all))
}


# -----------------------------------------------------------------------------
# 4. FIT G x E MODEL FOR ALL EIGHT TRAITS ------------------------------------
# -----------------------------------------------------------------------------

message("--- Fitting G x E mixed model for the 8 traits ---")
gxe_res <- herit_calc_GxE(df_all, TRAITS_C_ALL)


# -----------------------------------------------------------------------------
# 5. SAVE H^2 (GxE MODEL) AS CSV ---------------------------------------------
# -----------------------------------------------------------------------------
# CSV only -- no publication-style HTML/docx/png. The main-text Table 1 of
# the manuscript uses the no-GxE H^2 from Script A. This CSV is kept here
# purely as a comparison anchor (e.g. for a supplementary note showing how
# H^2 drops when G x E is included in the denominator).
# -----------------------------------------------------------------------------

message("--- Saving G x E H^2 values to CSV ---")

h2_gxe <- gxe_res$H2 %>%
  mutate(Trait      = unname(LABEL_MAP[pheno]),
         TraitClass = unname(TRAIT_CLASS_MAP[Trait]),
         H2_GxE     = round(heritability, 3)) %>%
  arrange(factor(TraitClass, levels = c("Nutritional", "Agromorphological")),
          desc(H2_GxE)) %>%
  select(Trait, TraitClass, H2_GxE)

write.csv(h2_gxe,
          file = file.path(OUT_ROOT, "tables", "H2_GxE.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 6. VARIANCE COMPONENTS : LONG + WIDE CSVs ----------------------------------
# -----------------------------------------------------------------------------
# Long  : one row per (Trait, Component) -- Variance + Percentage of total.
#         Direct backing table for the B1 stacked-bar plot.
# Wide  : one row per Trait, columns = % Genetics / % Season&Block / % G x E
#         / % Residual / Total variance. Convenient for manuscript tables and
#         for quoting specific numbers in the results text.
# -----------------------------------------------------------------------------

message("--- Saving variance-components tables (long + wide) ---")

var_components_long <- gxe_res$VarComp %>%
  group_by(Trait) %>%
  mutate(TotalVar   = sum(Variance),
         Percentage = round(Variance / TotalVar * 100, 2)) %>%
  ungroup() %>%
  mutate(TraitLabel = factor(unname(LABEL_MAP[Trait]),
                             levels = trait_order_8)) %>%
  arrange(TraitLabel,
          factor(Component, levels = c("Genetics", "Season & Block",
                                       "G x E", "Residual"))) %>%
  select(Trait, TraitLabel, Component, Variance, Percentage, TotalVar)

write.csv(var_components_long,
          file = file.path(OUT_ROOT, "tables",
                           "Variance_components_GxE.csv"),
          row.names = FALSE)

# Wide version: rows = trait, columns = % per component + total variance.
# Column names are sanitised (no spaces / no "x") to keep CSV consumers happy.
var_components_wide <- var_components_long %>%
  mutate(Component_clean = recode(Component,
                                  "Genetics"       = "Genetics_pct",
                                  "Season & Block" = "Season_Block_pct",
                                  "G x E"          = "GxE_pct",
                                  "Residual"       = "Residual_pct")) %>%
  select(TraitLabel, Component_clean, Percentage, TotalVar) %>%
  pivot_wider(names_from  = Component_clean,
              values_from = Percentage) %>%
  rename(Total_Variance = TotalVar) %>%
  select(TraitLabel, Genetics_pct, Season_Block_pct, GxE_pct, Residual_pct,
         Total_Variance) %>%
  arrange(factor(TraitLabel, levels = trait_order_8))

write.csv(var_components_wide,
          file = file.path(OUT_ROOT, "tables",
                           "Variance_components_GxE_WIDE.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 7. PLOT B1 : VARIANCE-PARTITIONING STACKED BAR (8 TRAITS) ------------------
# -----------------------------------------------------------------------------

message("--- Building B1 variance-partitioning plot (8 traits) ---")

var_perc <- var_components_long %>%
  mutate(TraitLabel = factor(TraitLabel, levels = trait_order_8),
         Component  = factor(Component,
                             levels = c("Residual", "G x E",
                                        "Season & Block", "Genetics")))

build_var_plot <- function(data) {
  ggplot(data, aes(x = TraitLabel, y = Percentage, fill = Component)) +
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = PAL_VARCOMP, name = NULL,
                      breaks = c("Genetics", "Season & Block",
                                 "G x E", "Residual")) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.02)),
      breaks = seq(0, 100, 25)
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(x = NULL, y = "Variance explained (%)") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))
}

p_var_8 <- build_var_plot(var_perc)
save_plot(p_var_8, "B1_variance_partitioning_8_traits", w = 8, h = 5)


# -----------------------------------------------------------------------------
# 8. REACTION-NORM DATA + PER-(TRAIT, ACCESSION) METRICS ---------------------
# -----------------------------------------------------------------------------
# Builds the underlying long table of genotype x season trait means and the
# per-(Trait, accession) metric panel used both to rank accessions and to
# pick the four representative highlights in B2.
#
# Both tables are persisted as CSVs so the results-writing section has
# access to the exact numbers behind every line in the figure.
# -----------------------------------------------------------------------------

message("--- Building reaction-norm long table + per-trait metrics ---")

# Per-trait locally cleaned long table of genotype x season means
# (nutritional traits only -- agromorpho traits aren't shown in B2).
gxe_long <- bind_rows(lapply(TRAITS_C_NUTRI, function(tr) {
  d <- clean_trait_local(df_all, tr)
  d %>%
    group_by(short_Tag, season) %>%
    summarise(Mean_Value = mean(.data[[tr]], na.rm = TRUE), .groups = "drop") %>%
    mutate(Trait = unname(LABEL_MAP[tr]))
})) %>%
  mutate(Trait = factor(Trait, levels = c("Protein", "Starch",
                                          "\u03B2-glucan", "Fiber")))

write.csv(gxe_long,
          file = file.path(OUT_ROOT, "tables", "gxe_long_data.csv"),
          row.names = FALSE)

# Per (trait, genotype) summary metrics for ranking.
gxe_metrics <- gxe_long %>%
  pivot_wider(names_from = season, values_from = Mean_Value,
              names_prefix = "y") %>%
  rowwise() %>%
  mutate(
    n_obs        = sum(!is.na(c_across(starts_with("y")))),
    change_19_21 = y2021 - y2019,
    stability    = sum(abs(diff(c(y2019, y2020, y2021))), na.rm = TRUE),
    range_all    = max(c(y2019, y2020, y2021), na.rm = TRUE) -
                   min(c(y2019, y2020, y2021), na.rm = TRUE),
    slope1       = y2020 - y2019,
    slope2       = y2021 - y2020,
    crossover    = ifelse(!is.na(slope1) & !is.na(slope2) &
                            sign(slope1) != sign(slope2), range_all, 0)
  ) %>%
  ungroup() %>%
  filter(n_obs == 3)   # require all three seasons for ranking

write.csv(gxe_metrics,
          file = file.path(OUT_ROOT, "tables",
                           "gxe_metrics_per_trait.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 9. PER-TRAIT HIGHLIGHTS : SELECTION + TABLE --------------------------------
# -----------------------------------------------------------------------------
# In each trait panel, four representative accessions are picked sequentially
# (already-picked accessions are excluded from later picks, so all four are
# guaranteed to be distinct within a panel):
#   * Most stable          = lowest `stability`
#   * Strongest increase   = highest `change_19_21`
#   * Strongest decrease   = lowest  `change_19_21`
#   * Strongest crossover  = highest `crossover`
# -----------------------------------------------------------------------------

message("--- Selecting per-trait highlights and exporting selection table ---")

per_trait_highlights <- gxe_metrics %>%
  group_by(Trait) %>%
  group_modify(~ {
    d <- .x
    ms <- d %>% arrange(stability)              %>% slice(1) %>% pull(short_Tag)
    si <- d %>% filter(short_Tag != ms)         %>% arrange(desc(change_19_21)) %>% slice(1) %>% pull(short_Tag)
    sd_<- d %>% filter(!short_Tag %in% c(ms, si)) %>% arrange(change_19_21)      %>% slice(1) %>% pull(short_Tag)
    sc <- d %>% filter(!short_Tag %in% c(ms, si, sd_)) %>%
            arrange(desc(crossover))    %>% slice(1) %>% pull(short_Tag)
    data.frame(short_Tag = c(ms, si, sd_, sc),
               Pattern   = c("Most stable", "Strongest increase",
                             "Strongest decrease", "Strongest crossover"))
  }) %>%
  ungroup() %>%
  mutate(Pattern = factor(Pattern, levels = names(PAL_HIGHLIGHT)))

# Export the selected highlights with their full metric values + season-by-
# season trait values -- so results-text sentences can quote any number
# without recomputation.
highlights_with_values <- per_trait_highlights %>%
  left_join(gxe_metrics, by = c("Trait", "short_Tag")) %>%
  mutate(metric_driving_selection = case_when(
    Pattern == "Most stable"         ~ stability,
    Pattern == "Strongest increase"  ~ change_19_21,
    Pattern == "Strongest decrease"  ~ change_19_21,
    Pattern == "Strongest crossover" ~ crossover
  )) %>%
  select(Trait, Pattern, short_Tag,
         y2019, y2020, y2021,
         change_19_21, stability, range_all,
         slope1, slope2, crossover,
         metric_driving_selection) %>%
  arrange(Trait, Pattern)

write.csv(highlights_with_values,
          file = file.path(OUT_ROOT, "tables",
                           "B2_reaction_norm_highlights_selected.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 10. PLOT B2 : REACTION NORM (PER-TRAIT HIGHLIGHTS) -------------------------
# -----------------------------------------------------------------------------
# Background = ALL ~290 genotypes drawn as faint grey LINES, no points.
# 4 highlighted accessions per trait panel (picked independently per panel).
# -----------------------------------------------------------------------------

message("--- Building B2 reaction-norm plot (per-trait highlights) ---")

plot_data_pertrait <- gxe_long %>%
  left_join(per_trait_highlights, by = c("Trait", "short_Tag")) %>%
  mutate(Group  = ifelse(is.na(Pattern), "Background", "Highlight"),
         season = factor(season, levels = c("2019", "2020", "2021")))

p_reaction_pertrait <- ggplot() +
  geom_line(data = plot_data_pertrait %>% filter(Group == "Background"),
            aes(x = season, y = Mean_Value, group = short_Tag),
            colour = "grey80", linewidth = 0.3, alpha = 0.45) +
  geom_line(data = plot_data_pertrait %>% filter(Group == "Highlight"),
            aes(x = season, y = Mean_Value, group = short_Tag,
                colour = Pattern),
            linewidth = 1.15) +
  geom_point(data = plot_data_pertrait %>% filter(Group == "Highlight"),
             aes(x = season, y = Mean_Value, colour = Pattern),
             size = 2.1) +
  scale_colour_manual(values = PAL_HIGHLIGHT, name = "Representative pattern",
                      labels = c("Most stable"         = "stable",
                                 "Strongest increase"  = "increase",
                                 "Strongest decrease"  = "decrease",
                                 "Strongest crossover" = "crossover")) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 2) +
  labs(x = "Growing season", y = "Centered phenotypic value") +
  theme_pub() +
  guides(colour = guide_legend(nrow = 1, title.position = "top",
                               title.hjust = 0.5)) +
  theme(legend.title = element_text(hjust = 0.5))

save_plot(p_reaction_pertrait, "B2_reaction_norm_per_trait_highlights",
          w = 8, h = 6.5)


# -----------------------------------------------------------------------------
# 11. FINAL MESSAGE ----------------------------------------------------------
# -----------------------------------------------------------------------------
message("\n--- Script B v2 (streamlined) finished. ---")
message("    Outputs root:  ", OUT_ROOT)
message("    PDFs:          ", file.path(OUT_ROOT, "pdf"))
message("    PNGs:          ", file.path(OUT_ROOT, "png"))
message("    CSV / tables:  ", file.path(OUT_ROOT, "tables"))
