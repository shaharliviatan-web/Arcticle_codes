# =============================================================================
# Script A v3 (STREAMLINED) : Phenotypic analysis WITHOUT G x E interaction
# =============================================================================
# Streamlined version of `01_phenotypic_analysis_no_GxE_v2_10traits.R`.
# Retains every filter, every QC step, and the same LMM model as v2, but
# produces ONLY the analyses required for Subsection 1 of the Results chapter:
#
#   1. BLUPs for all 10 traits (GWAS phenotype input + descriptive stats)
#   2. Broad-sense heritability (H^2) table for all 10 traits
#   3. Section A4: sampling-site boxplot of the 4 nutritional traits
#                  on the centered BLUP scale + per-region summary tables
#                  (centered AND absolute % scale) for results-text use
#   4. Section A12 Panel A: 4 x 4 nutritional Pearson-correlation heatmap
#                  TWO variants:
#                     (i)  raw p-value asterisks
#                     (ii) BH-FDR-corrected asterisks (local; 6 nutri-nutri pairs)
#   5. Section A12 Panel B: 4 nutri x 4 morpho Pearson-correlation dotplot
#                  with BH-FDR-corrected asterisks (local; 16 nutri-morpho pairs)
#
# All filters, QC, ID re-formatting, per-season centering, and the random-
# effects model are IDENTICAL to v2:
#
#   Model used :  trait ~ (1 | short_Tag) + (1 | season_Bl)
#   Package    :  lme4 v1.1.37
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
# Output folder layout is preserved exactly as in v2:
#   outputs/subsection_1/pdf/
#   outputs/subsection_1/png/
#   outputs/subsection_1/tables/
#
# Author : Shahar Liviatan
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETUP -------------------------------------------------------------------
# -----------------------------------------------------------------------------

# ---- Auto-install any missing packages (helpful on a fresh Posit Cloud) ---
.required <- c("ggplot2", "lme4", "insight", "dplyr", "tidyr", "stringr",
               "scales", "gt", "flextable", "officer", "ggpubr", "tibble")
.missing  <- .required[!vapply(.required, requireNamespace,
                               logical(1), quietly = TRUE)]
if (length(.missing)) {
  message("Installing missing packages: ", paste(.missing, collapse = ", "))
  install.packages(.missing)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(lme4)
  library(insight)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
  library(gt)
  library(flextable)
  library(officer)
  library(ggpubr)
  library(tibble)
})

# ---- Configuration --------------------------------------------------------
INPUT_CSV <- "all years barley.csv"   # change this if your file is elsewhere
OUT_ROOT  <- file.path("outputs", "subsection_1")

dir.create(file.path(OUT_ROOT, "pdf"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "png"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "tables"), recursive = TRUE, showWarnings = FALSE)

# ---- Trait definitions ----------------------------------------------------
# Raw (As-is, percent / absolute scale) column names exactly as in the CSV.
TRAITS_RAW_NUTRI <- c("ProteinAsis.", "StarchAsis.", "BetaglucansAsis.", "FiberAsis.")
TRAITS_RAW_AGRO  <- c("flowering_time", "Tiller.Length", "Grain_weight",
                      "number_of_grains", "Number.Of.Tillers", "Spike.Length")
TRAITS_RAW_ALL   <- c(TRAITS_RAW_NUTRI, TRAITS_RAW_AGRO)

# Centered column names (created during preprocessing).
TRAITS_C_NUTRI   <- c("ProteinCenterd", "StarchCenterd",
                      "BetaglucansCenterd", "FiberCenterd")
TRAITS_C_AGRO    <- c("flowering_timeCenterd", "Tiller.LengthCenterd",
                      "Grain_weightCenterd", "number_of_grainsCenterd",
                      "Number.Of.TillersCenterd", "Spike.LengthCenterd")
TRAITS_C_ALL     <- c(TRAITS_C_NUTRI, TRAITS_C_AGRO)

# Display labels (publication style: Greek beta, American English).
LABEL_MAP <- c(
  ProteinAsis.            = "Protein",          ProteinCenterd            = "Protein",
  StarchAsis.             = "Starch",           StarchCenterd             = "Starch",
  BetaglucansAsis.        = "\u03B2-glucan",    BetaglucansCenterd        = "\u03B2-glucan",
  FiberAsis.              = "Fiber",            FiberCenterd              = "Fiber",
  flowering_time          = "Flowering time",   flowering_timeCenterd     = "Flowering time",
  Tiller.Length           = "Plant height",     Tiller.LengthCenterd      = "Plant height",
  Grain_weight            = "Grain weight",     Grain_weightCenterd       = "Grain weight",
  number_of_grains        = "Grain number",     number_of_grainsCenterd   = "Grain number",
  Number.Of.Tillers       = "Tillers",          Number.Of.TillersCenterd  = "Tillers",
  Spike.Length            = "Spike length",     Spike.LengthCenterd       = "Spike length"
)

TRAIT_CLASS_MAP <- c(
  "Protein" = "Nutritional", "Starch" = "Nutritional",
  "\u03B2-glucan" = "Nutritional", "Fiber" = "Nutritional",
  "Flowering time" = "Agromorphological", "Plant height" = "Agromorphological",
  "Grain weight" = "Agromorphological", "Grain number" = "Agromorphological",
  "Tillers" = "Agromorphological", "Spike length" = "Agromorphological"
)

# Trait factor orders used downstream.
trait_order      <- c("Protein", "Starch", "\u03B2-glucan", "Fiber")
trait_order_all  <- c("Protein", "Starch", "\u03B2-glucan", "Fiber",
                      "Flowering time", "Plant height", "Grain weight",
                      "Grain number", "Tillers", "Spike length")

# ---- Colour palette (kept consistent across the whole subsection) --------
PAL_BLUP   <- "#9E9E9E"               # neutral grey (kept for back-compat)

# ---- Sampling-site -> ecological region mapping --------------------------
# Site code = first two digits of short_Tag (after the HS-strip reformat).
# Site 04 (Hachola) is filtered out upstream; kept here for completeness.
SITE_NAME_MAP <- c(
  "01" = "Meizar",          "02" = "Beit Zeida",     "03" = "Mount Hermon",
  "05" = "Misgav River",
  "06" = "Zichron Yaakov",  "07" = "Mount Tayasim",  "08" = "Rohama",
  "09" = "Nahala",          "10" = "Zikim",
  "11" = "Masciot",         "12" = "Mount Ramon",    "13" = "Masada",
  "14" = "Maale Akrabim",   "15" = "Arad",
  "16" = "Cabri",           "17" = "Meron",          "18" = "Alon Hagalil",
  "19" = "Afek",            "20" = "Rechasim",
  "21" = "Meirav",          "22" = "Sde Ilan",       "23" = "Gesher",
  "24" = "Givat Hamoreh",   "25" = "Tirat Zvi",
  "26" = "Megido Forest",   "27" = "Revivim",        "28" = "Orim",
  "29" = "Hagosh Junction", "30" = "Amona"
)

SITE_REGION_MAP <- c(
  "01" = "North",  "02" = "North",  "03" = "North",
  "05" = "North",
  "06" = "Coast",  "07" = "Coast",  "08" = "Coast",
  "09" = "Coast",  "10" = "Coast",
  "11" = "Desert", "12" = "Desert", "13" = "Desert",
  "14" = "Desert", "15" = "Desert",
  "16" = "HZ1 (North-Coast)",  "17" = "HZ1 (North-Coast)",
  "18" = "HZ1 (North-Coast)",  "19" = "HZ1 (North-Coast)",
  "20" = "HZ1 (North-Coast)",
  "21" = "HZ2 (North-Desert)", "22" = "HZ2 (North-Desert)",
  "23" = "HZ2 (North-Desert)", "24" = "HZ2 (North-Desert)",
  "25" = "HZ2 (North-Desert)",
  "26" = "HZ3 (Coast-Desert)", "27" = "HZ3 (Coast-Desert)",
  "28" = "HZ3 (Coast-Desert)", "29" = "HZ3 (Coast-Desert)",
  "30" = "HZ3 (Coast-Desert)"
)

REGION_ORDER <- c("North", "Coast", "Desert",
                  "HZ1 (North-Coast)",
                  "HZ2 (North-Desert)",
                  "HZ3 (Coast-Desert)")

# Region palette: three pure ecoregions use strong, semantically meaningful
# hues (cool blue = North, green = Coast, warm red = Desert); the three
# hybrid zones use a blend tone derived from their parent regions.
PAL_REGION <- c(
  "North"               = "#2166AC",   # blue
  "Coast"               = "#1B7837",   # green
  "Desert"              = "#B2182B",   # red
  "HZ1 (North-Coast)"   = "#67A9CF",   # light blue   (N + C)
  "HZ2 (North-Desert)"  = "#762A83",   # purple       (N + D)
  "HZ3 (Coast-Desert)"  = "#E08214"    # orange       (C + D)
)

# ---- Publication theme (used by every plot in this script) ---------------
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      plot.title       = element_blank(),         # titles added at figure assembly
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

# ---- Save helper: writes both PDF (vector) and PNG @ 600 dpi ------------
save_plot <- function(plot, name, w = 7, h = 5) {
  ggsave(file.path(OUT_ROOT, "pdf", paste0(name, ".pdf")),
         plot, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT_ROOT, "png", paste0(name, ".png")),
         plot, width = w, height = h, dpi = 600, bg = "white")
  invisible(plot)
}


# -----------------------------------------------------------------------------
# 2. PREPROCESSING -----------------------------------------------------------
# -----------------------------------------------------------------------------

message("--- Loading and preprocessing the field-experiment data ---")

df_all <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)

# Filter 1: Starch + 8.85 correction for 2021 only.
# (Post-NIR adjustment retained from Tal Amar's pipeline; see project notes.)
df_all$StarchAsis.[df_all$season == 2021] <-
  df_all$StarchAsis.[df_all$season == 2021] + 8.85

# Filter 2: drop Block 4 of season 2019 (low inter-block concordance).
df_all <- df_all[!(df_all$season == 2019 & df_all$Block_no == 4), ]

# Filter 3: drop control lines (Morex, Clipper) -- tags starting with c/C/M.
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
# season is kept numeric here; Script B coerces it to a factor for the
# G x E interaction term.

# Per-season mean-centering (all ten traits).
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

# H^2 calculator (no G x E).
herit_calc <- function(df, traits) {
  out <- vector("numeric", length(traits))
  for (i in seq_along(traits)) {
    tr <- traits[i]
    d  <- clean_trait_local(df, tr)
    fit <- lmer(d[[tr]] ~ (1 | short_Tag) + (1 | season_Bl), data = d)
    vc  <- get_variance(fit)
    out[i] <- vc$var.intercept[1] / (vc$var.random + vc$var.residual)
    message(sprintf("    H2 (%s) = %.3f", tr, out[i]))
  }
  data.frame(pheno = traits, heritability = round(out, 3))
}

# BLUP extractor (no G x E).
blup_calc <- function(df, trait) {
  d <- clean_trait_local(df, trait)
  fit <- lmer(d[[trait]] ~ (1 | short_Tag) + (1 | season_Bl), data = d)
  bl  <- coef(fit)$short_Tag
  data.frame(Genotype = rownames(bl), value = bl[, 1]) %>%
    rename(!!sym(trait) := value)
}

# Significance-asterisk helper based on RAW p-values.
star_fn <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", ""))))
}

# Significance-asterisk helper based on (BH-FDR-corrected) q-values.
local_sig <- function(q) {
  case_when(
    is.na(q)   ~ "",
    q < 0.001  ~ "***",
    q < 0.01   ~ "**",
    q < 0.05   ~ "*",
    TRUE       ~ ""
  )
}


# -----------------------------------------------------------------------------
# 4. CALCULATE H^2 AND BLUPs  (10 traits) ------------------------------------
# -----------------------------------------------------------------------------

message("--- Estimating broad-sense heritability for the 10 traits ---")
herit_df <- herit_calc(df_all, TRAITS_C_ALL)

message("--- Estimating BLUPs for the 10 traits ---")
blup_list <- lapply(TRAITS_C_ALL, function(tr) blup_calc(df_all, tr))
blup_df   <- Reduce(function(a, b) merge(a, b, by = "Genotype", all = TRUE), blup_list)

# Persist the BLUP table (this is the GWAS phenotype file -- must match exactly).
write.csv(blup_df,
          file = file.path(OUT_ROOT, "tables", "BLUP_10_traits.csv"),
          row.names = FALSE)
write.csv(herit_df,
          file = file.path(OUT_ROOT, "tables", "H2_10_traits_raw.csv"),
          row.names = FALSE)

# ---- BLUP descriptive statistics (NEW: convenience table for results text) ---
# Per-trait mean, SD, median, min, max, range, n_accessions across all genotypes.
# Useful for headline sentences such as "Fiber BLUPs ranged from X to Y
# (mean +/- SD = Z +/- W; n = N accessions)".
blup_long <- blup_df %>%
  pivot_longer(-Genotype, names_to = "trait_raw", values_to = "blup") %>%
  mutate(Trait      = factor(unname(LABEL_MAP[trait_raw]), levels = trait_order_all),
         TraitClass = unname(TRAIT_CLASS_MAP[as.character(Trait)]))

blup_descr <- blup_long %>%
  filter(!is.na(blup)) %>%
  group_by(Trait, TraitClass) %>%
  summarise(
    n_accessions = n(),
    mean_blup    = round(mean(blup),   4),
    sd_blup      = round(sd(blup),     4),
    median_blup  = round(median(blup), 4),
    min_blup     = round(min(blup),    4),
    max_blup     = round(max(blup),    4),
    range_blup   = round(max(blup) - min(blup), 4),
    .groups      = "drop"
  ) %>%
  arrange(factor(TraitClass, levels = c("Nutritional", "Agromorphological")),
          Trait)

write.csv(blup_descr,
          file = file.path(OUT_ROOT, "tables", "BLUP_descriptive_stats.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. TABLE 1 : BROAD-SENSE HERITABILITY  (10 traits) -------------------------
# -----------------------------------------------------------------------------
# Saved as HTML (via gt), .docx (via flextable), and PNG.
# -----------------------------------------------------------------------------

message("--- Building Table 1 (H^2 of 10 traits) ---")

h2_table <- herit_df %>%
  mutate(Trait      = unname(LABEL_MAP[pheno]),
         TraitClass = unname(TRAIT_CLASS_MAP[Trait]),
         H2         = sprintf("%.3f", heritability)) %>%
  arrange(factor(TraitClass, levels = c("Nutritional", "Agromorphological")),
          desc(heritability)) %>%
  select(Trait, TraitClass, H2)

write.csv(h2_table,
          file = file.path(OUT_ROOT, "tables", "H2_10_traits_clean.csv"),
          row.names = FALSE)

# ---- HTML via gt ------------------------------------------------------------
gt_table <- h2_table %>%
  gt() %>%
  tab_header(
    title = md("**Table 1.** Broad-sense heritability (*H*\u00B2) of grain nutritional and agromorphological traits in wild barley")
  ) %>%
  cols_label(
    Trait      = "Trait",
    TraitClass = "Trait class",
    H2         = md("*H*\u00B2")
  ) %>%
  cols_align(align = "left",   columns = c(Trait, TraitClass)) %>%
  cols_align(align = "center", columns = H2) %>%
  tab_row_group(
    label = "Agromorphological",
    rows  = TraitClass == "Agromorphological"
  ) %>%
  tab_row_group(
    label = "Nutritional",
    rows  = TraitClass == "Nutritional"
  ) %>%
  cols_hide(columns = TraitClass) %>%
  tab_source_note(
    source_note = md("*H*\u00B2 was estimated from a linear mixed model fit per trait with genotype and season-by-block as crossed random effects (lme4 v1.1.37). Per-trait local quality control was applied: missing values were excluded and observations exceeding +/- 3 SD from the trait mean were removed prior to model fitting. *H*\u00B2 = \u03C3\u00B2_G / (\u03C3\u00B2_G + \u03C3\u00B2_E + \u03C3\u00B2_R).")
  ) %>%
  tab_options(
    heading.title.font.size      = 14,
    column_labels.font.weight    = "bold",
    table.font.names             = "Arial",
    table.border.top.color       = "black",
    table.border.top.width       = px(2),
    table.border.bottom.color    = "black",
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(1.5),
    row_group.border.top.color   = "grey70",
    row_group.border.bottom.color = "grey70",
    table_body.hlines.color      = "white",
    table.background.color       = "white"
  )

gtsave(gt_table, file.path(OUT_ROOT, "tables", "Table1_H2.html"))

# ---- .docx + PNG via flextable ---------------------------------------------
ft <- flextable(h2_table) %>%
  set_header_labels(
    Trait      = "Trait",
    TraitClass = "Trait class",
    H2         = "H\u00B2"
  ) %>%
  merge_v(j = "TraitClass") %>%
  align(j = "H2", align = "center", part = "all") %>%
  align(j = c("Trait", "TraitClass"), align = "left", part = "all") %>%
  bold(part = "header") %>%
  flextable::font(fontname = "Arial", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  border_remove() %>%
  hline_top(border = officer::fp_border(color = "black", width = 1.5),
            part = "header") %>%
  hline_bottom(border = officer::fp_border(color = "black", width = 1.5),
               part = "header") %>%
  hline_bottom(border = officer::fp_border(color = "black", width = 1.5),
               part = "body") %>%
  add_footer_lines(values = "H\u00B2 estimated from an LMM fit per trait with genotype and season-by-block as crossed random effects (lme4 v1.1.37). Per-trait local QC: NA exclusion + +/- 3 SD outlier removal. H\u00B2 = sigma2_G / (sigma2_G + sigma2_E + sigma2_R).") %>%
  italic(part = "footer") %>%
  flextable::fontsize(size = 8, part = "footer") %>%
  autofit()

save_as_docx(ft,
             path = file.path(OUT_ROOT, "tables", "Table1_H2.docx"))

png_path <- file.path(OUT_ROOT, "tables", "Table1_H2.png")
png_ok <- tryCatch({
  save_as_image(ft, path = png_path, webshot = "webshot2")
  TRUE
}, error = function(e) {
  message("    flextable PNG export failed (", conditionMessage(e),
          "); falling back to ggpubr.")
  FALSE
})

if (!png_ok) {
  h2_display <- h2_table
  names(h2_display) <- c("Trait", "Trait class", "H\u00B2")
  ggp <- ggpubr::ggtexttable(
    as.data.frame(h2_display),
    rows  = NULL,
    theme = ggpubr::ttheme(
      base_style = "blank",
      colnames.style = ggpubr::colnames_style(
        face = "bold", size = 11, fill = "white",
        linewidth = 1, linecolor = "black"
      ),
      tbody.style = ggpubr::tbody_style(
        size = 10, fill = "white",
        linewidth = 0.4, linecolor = "grey70"
      )
    )
  )
  ggsave(png_path, ggp, width = 5.5, height = 3.5, dpi = 600, bg = "white")
}


# -----------------------------------------------------------------------------
# 6. SECTION A4 : SAMPLING-SITE BOXPLOTS (centered BLUP) ---------------------
# -----------------------------------------------------------------------------
# Each box represents one sampling site (first two digits of short_Tag, after
# the HS-strip / site_id_within reformatting). Each point is the BLUP of one
# accession sampled at that site. Four nutritional traits are panelled.
#
# Only the CENTERED-BLUP plot is produced in this streamlined script. The
# absolute-%-scale summary tables (mean +/- SD per region in % units) are
# still written to CSV because they are useful for the results-text
# wording of "X% protein, Y% starch" etc.
# -----------------------------------------------------------------------------

message("--- Building Section A4: sampling-site boxplot (centered BLUP) ---")

# Trait grand means -- used to shift BLUPs from the centered scale to the
# absolute % scale in the per-region summary tables.
grand_means <- sapply(TRAITS_RAW_NUTRI, function(tr)
  mean(clean_trait_local(df_all, tr)[[tr]], na.rm = TRUE))
names(grand_means) <- unname(LABEL_MAP[TRAITS_RAW_NUTRI])

# Tidy long table of BLUPs for the 4 nutritional traits, with the ecological
# region and the human-readable location name attached.
site_long_centered <- blup_df %>%
  select(Genotype, all_of(TRAITS_C_NUTRI)) %>%
  pivot_longer(-Genotype, names_to = "trait_raw", values_to = "value") %>%
  mutate(
    Trait    = factor(unname(LABEL_MAP[trait_raw]), levels = trait_order),
    site     = substr(as.character(Genotype), 1, 2),
    location = unname(SITE_NAME_MAP[site]),
    region   = factor(unname(SITE_REGION_MAP[site]), levels = REGION_ORDER)
  ) %>%
  filter(!is.na(value)) %>%
  select(region, site, location, short_Tag = Genotype, Trait, value)

# Absolute-scale version (used ONLY for the summary tables, not for the plot).
site_long_absolute <- site_long_centered %>%
  mutate(value = value + grand_means[as.character(Trait)])

# Persist a single tidy CSV with BOTH scales (long format) -- accession-level.
site_csv <- site_long_centered %>%
  rename(value_centered = value) %>%
  left_join(site_long_absolute %>%
              rename(value_absolute = value),
            by = c("region", "site", "location", "short_Tag", "Trait")) %>%
  arrange(Trait, region, site, short_Tag)

write.csv(site_csv,
          file = file.path(OUT_ROOT, "tables", "A4_site_boxplot_data.csv"),
          row.names = FALSE)

# Order sites along the x-axis BY REGION (North -> Coast -> Desert ->
# HZ1 -> HZ2 -> HZ3), then numerically within each region.
site_levels <- site_csv %>%
  distinct(region, site) %>%
  arrange(region, site) %>%
  pull(site)
site_long_centered$site <- factor(site_long_centered$site, levels = site_levels)
site_long_absolute$site <- factor(site_long_absolute$site, levels = site_levels)

# ---- Plot: centered BLUP ---------------------------------------------------
p_site_centered <- ggplot(site_long_centered,
                          aes(x = site, y = value, fill = region)) +
  geom_boxplot(colour = "grey30", outlier.shape = NA,
               linewidth = 0.35, alpha = 0.75) +
  geom_jitter(width = 0.18, height = 0, size = 0.9,
              alpha = 0.55, colour = "grey20") +
  facet_wrap(~ Trait, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = PAL_REGION, name = "Ecological region",
                    drop = FALSE) +
  labs(x = "Sampling site", y = "BLUP (centered)") +
  theme_pub() +
  theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

save_plot(p_site_centered, "A4_site_boxplots_centered", w = 11, h = 7.5)

# ---- A4 SUPPLEMENT : Per-region summary statistics (for manuscript text) ---
# Computes mean +/- SD and accession count for every (trait x region)
# combination, on BOTH the centered BLUP scale AND the absolute % scale.
# A "wide" version (traits as columns) is also saved to make it easy to
# copy numbers directly into the manuscript.
# ----------------------------------------------------------------------------

message("--- Building Section A4 supplement: per-region summary stats ---")

summarise_by_region <- function(data, value_col_label) {
  data %>%
    group_by(region, Trait) %>%
    summarise(
      n_accessions = n(),
      mean_blup    = round(mean(value, na.rm = TRUE), 4),
      sd_blup      = round(sd(value,   na.rm = TRUE), 4),
      median_blup  = round(median(value, na.rm = TRUE), 4),
      min_blup     = round(min(value,  na.rm = TRUE), 4),
      max_blup     = round(max(value,  na.rm = TRUE), 4),
      .groups      = "drop"
    ) %>%
    arrange(Trait, region) %>%
    rename_with(~ sub("blup", value_col_label, .x),
                starts_with("mean_") | starts_with("sd_") |
                  starts_with("median_") | starts_with("min_") | starts_with("max_"))
}

region_summary_centered <- summarise_by_region(site_long_centered, "centered")
region_summary_absolute <- summarise_by_region(site_long_absolute, "absolute_pct")

write.csv(region_summary_centered,
          file = file.path(OUT_ROOT, "tables", "A4_region_summary_centered.csv"),
          row.names = FALSE)
write.csv(region_summary_absolute,
          file = file.path(OUT_ROOT, "tables", "A4_region_summary_absolute.csv"),
          row.names = FALSE)

# Wide version (traits as columns, mean +/- SD formatted as "X.XXX +/- X.XXX")
# - convenient for pasting directly into the manuscript table.
make_wide_summary <- function(summary_df, val_col, sd_col) {
  summary_df %>%
    mutate(mean_sd = paste0(
      formatC(.data[[val_col]], format = "f", digits = 3), " \u00B1 ",
      formatC(.data[[sd_col]],  format = "f", digits = 3)
    )) %>%
    select(region, Trait, n_accessions, mean_sd) %>%
    tidyr::pivot_wider(names_from  = Trait,
                       values_from = c(n_accessions, mean_sd),
                       names_glue  = "{Trait}_{.value}") %>%
    arrange(factor(region, levels = REGION_ORDER))
}

region_wide_centered <- make_wide_summary(region_summary_centered,
                                          "mean_centered", "sd_centered")
region_wide_absolute <- make_wide_summary(region_summary_absolute,
                                          "mean_absolute_pct", "sd_absolute_pct")

write.csv(region_wide_centered,
          file = file.path(OUT_ROOT, "tables", "A4_region_summary_centered_WIDE.csv"),
          row.names = FALSE)
write.csv(region_wide_absolute,
          file = file.path(OUT_ROOT, "tables", "A4_region_summary_absolute_WIDE.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 7. CORRELATION COMPUTATIONS (feeds A12 Panels A and B) ---------------------
# -----------------------------------------------------------------------------
# Pairwise Pearson correlations between BLUPs of the 10 traits. The raw
# correlation matrix and raw p-value matrix are saved as reference tables;
# no plot is produced from them in this streamlined script. The 45-test
# global BH-FDR correction is intentionally NOT computed -- Panels A and B
# each apply their OWN local FDR correction over their own pair set.
# -----------------------------------------------------------------------------

message("--- Computing pairwise Pearson correlations between BLUPs (10 traits) ---")

blup_matrix <- blup_df %>%
  select(Genotype, all_of(TRAITS_C_ALL)) %>%
  rename_with(~ unname(LABEL_MAP[.x]), all_of(TRAITS_C_ALL)) %>%
  tibble::column_to_rownames("Genotype")
blup_matrix <- blup_matrix[, trait_order_all]

cor_matrix_with_p <- function(x) {
  n_t   <- ncol(x)
  cn    <- colnames(x)
  r_mat <- matrix(NA_real_, n_t, n_t, dimnames = list(cn, cn))
  p_mat <- matrix(NA_real_, n_t, n_t, dimnames = list(cn, cn))
  n_mat <- matrix(0L,       n_t, n_t, dimnames = list(cn, cn))
  for (i in seq_len(n_t)) for (j in seq_len(n_t)) {
    if (i == j) {
      r_mat[i, j] <- 1
      p_mat[i, j] <- 0
      n_mat[i, j] <- sum(!is.na(x[, i]))
    } else {
      ok <- !is.na(x[, i]) & !is.na(x[, j])
      n_mat[i, j] <- sum(ok)
      if (n_mat[i, j] >= 3) {
        ct <- cor.test(x[ok, i], x[ok, j])
        r_mat[i, j] <- unname(ct$estimate)
        p_mat[i, j] <- ct$p.value
      }
    }
  }
  list(r = r_mat, p = p_mat, n = n_mat)
}
cm <- cor_matrix_with_p(as.matrix(blup_matrix))

# Wide CSV exports (raw r and raw p matrices -- master reference tables).
write.csv(round(cm$r, 4),
          file = file.path(OUT_ROOT, "tables", "A5_r_matrix_10x10.csv"))
write.csv(signif(cm$p, 4),
          file = file.path(OUT_ROOT, "tables", "A5_p_matrix_10x10.csv"))

# Long edge list with raw p-value significance flags (no FDR).
edge_long <- expand.grid(Trait1 = rownames(cm$r),
                         Trait2 = colnames(cm$r),
                         stringsAsFactors = FALSE) %>%
  mutate(
    i = match(Trait1, trait_order_all),
    j = match(Trait2, trait_order_all),
    r = mapply(function(a, b) cm$r[a, b], Trait1, Trait2),
    p = mapply(function(a, b) cm$p[a, b], Trait1, Trait2),
    n = mapply(function(a, b) cm$n[a, b], Trait1, Trait2),
    sig = star_fn(p)
  ) %>%
  filter(i < j) %>%                # upper triangle only (45 unique pairs)
  mutate(Class1 = unname(TRAIT_CLASS_MAP[Trait1]),
         Class2 = unname(TRAIT_CLASS_MAP[Trait2]),
         EdgeType = case_when(
           Class1 == "Nutritional" & Class2 == "Nutritional" ~ "nutri-nutri",
           Class1 == "Nutritional" | Class2 == "Nutritional" ~ "nutri-morpho",
           TRUE                                              ~ "morpho-morpho"
         )) %>%
  select(Trait1, Trait2, Class1, Class2, EdgeType, r, p, n, sig)

write.csv(edge_long,
          file = file.path(OUT_ROOT, "tables", "A5_edges_long_10x10.csv"),
          row.names = FALSE)

# Pull all nutri-morpho pairs (used by Panel B) and persist them.
nutri_morpho_pairs <- edge_long %>%
  filter(EdgeType == "nutri-morpho") %>%
  mutate(nutri  = ifelse(Class1 == "Nutritional", Trait1, Trait2),
         morpho = ifelse(Class1 == "Nutritional", Trait2, Trait1)) %>%
  select(nutri, morpho, r, p, n, sig)

write.csv(nutri_morpho_pairs,
          file = file.path(OUT_ROOT, "tables", "A8_nutri_morpho_pairs.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 8. SECTION A12 PANEL A : 4 x 4 NUTRITIONAL CORRELATION HEATMAP -------------
# -----------------------------------------------------------------------------
# TWO variants are produced:
#   (i)  raw-p version    : asterisks based on raw p-values from cor.test
#   (ii) local-FDR version : asterisks based on BH-FDR computed locally
#                            within the 6 nutri-nutri pairs only
# Raw r values are identical in both variants -- only the asterisks change.
# -----------------------------------------------------------------------------

message("--- Building Section A12 Panel A: nutri-nutri heatmap (2 variants) ---")

# Compact triangle: drop the "Protein" column (no tiles) and the "Fiber" row
# (no tiles), reducing the 4 x 4 grid down to a 3 x 3 triangular block.
x_levels <- c("Starch", "\u03B2-glucan", "Fiber")    # no Protein column
y_levels <- c("Protein", "Starch", "\u03B2-glucan")  # no Fiber row

# ---- (i) raw-p variant -----------------------------------------------------
nutri_rawP <- edge_long %>%
  filter(EdgeType == "nutri-nutri") %>%             # the 6 nutri-nutri pairs
  mutate(sig_raw = star_fn(p))

write.csv(nutri_rawP,
          file = file.path(OUT_ROOT, "tables",
                           "A12_nutri_pairs_uncorrectedP.csv"),
          row.names = FALSE)

heat_nutri_rawP <- nutri_rawP %>%
  mutate(Trait1 = factor(Trait1, levels = y_levels),
         Trait2 = factor(Trait2, levels = x_levels),
         label  = sprintf("%.2f%s", r, sig_raw))

p_panel_A_rawP <- ggplot(heat_nutri_rawP,
                         aes(x = Trait2, y = Trait1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label  = label,
                colour = ifelse(abs(r) > 0.45, "white", "black")),
            size = 4.0, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       name = "Pearson r") +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(y_levels), drop = FALSE) +
  scale_x_discrete(limits = x_levels,      drop = FALSE,
                   position = "top") +
  labs(x = NULL, y = NULL,
       caption = "* p < 0.05   ** p < 0.01   *** p < 0.001 (uncorrected)") +
  coord_equal() +
  theme_pub() +
  theme(panel.grid   = element_blank(),
        axis.text.x  = element_text(angle = 0, hjust = 0.5),
        plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"))

save_plot(p_panel_A_rawP,
          "A12_panelA_nutri_heatmap_rawP", w = 4.2, h = 4.0)

# ---- (ii) local-FDR variant (BH-FDR over the 6 nutri-nutri pairs) ---------
nutri_localFDR <- edge_long %>%
  filter(EdgeType == "nutri-nutri") %>%
  mutate(q_local   = p.adjust(p, method = "BH"),
         sig_local = local_sig(q_local))

write.csv(nutri_localFDR,
          file = file.path(OUT_ROOT, "tables",
                           "A12_nutri_pairs_localFDR.csv"),
          row.names = FALSE)

heat_nutri_localFDR <- nutri_localFDR %>%
  mutate(Trait1 = factor(Trait1, levels = y_levels),
         Trait2 = factor(Trait2, levels = x_levels),
         label  = sprintf("%.2f%s", r, sig_local))

p_panel_A_localFDR <- ggplot(heat_nutri_localFDR,
                             aes(x = Trait2, y = Trait1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label  = label,
                colour = ifelse(abs(r) > 0.45, "white", "black")),
            size = 4.0, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1),
                       name = "Pearson r") +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(y_levels), drop = FALSE) +
  scale_x_discrete(limits = x_levels,      drop = FALSE,
                   position = "top") +
  labs(x = NULL, y = NULL,
       caption = paste0("BH-FDR over the 6 nutri-nutri pairs.  ",
                        "* q < 0.05   ** q < 0.01   *** q < 0.001")) +
  coord_equal() +
  theme_pub() +
  theme(panel.grid   = element_blank(),
        axis.text.x  = element_text(angle = 0, hjust = 0.5),
        plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"))

save_plot(p_panel_A_localFDR,
          "A12_panelA_nutri_heatmap_localFDR", w = 4.2, h = 4.0)


# -----------------------------------------------------------------------------
# 9. SECTION A12 PANEL B : 4 NUTRI x 4 MORPHO DOTPLOT (local FDR) ------------
# -----------------------------------------------------------------------------
# Significance asterisks are based on BH-FDR computed LOCALLY within the
# 16 nutri x morpho pairs displayed in this plot. Raw r and raw p are
# unchanged (they are pairwise quantities and don't depend on the corpus).
# -----------------------------------------------------------------------------

message("--- Building Section A12 Panel B: nutri-morpho dotplot (local FDR) ---")

NETWORK_MORPHO_4 <- c("Flowering time", "Plant height",
                      "Grain weight",   "Grain number")
NETWORK_NUTRI    <- c("Protein", "Starch", "\u03B2-glucan", "Fiber")

nutri_morpho_local <- nutri_morpho_pairs %>%
  filter(morpho %in% NETWORK_MORPHO_4) %>%          # 16 nutri x morpho pairs
  mutate(q_local   = p.adjust(p, method = "BH"),
         sig_local = local_sig(q_local),
         nutri     = factor(nutri, levels = NETWORK_NUTRI))

write.csv(nutri_morpho_local,
          file = file.path(OUT_ROOT, "tables",
                           "A12_nutri_morpho_pairs_localFDR.csv"),
          row.names = FALSE)

# Sort morpho traits BY r value WITHIN EACH nutri panel.
# Trick: tag each x-level with the nutri index, build a global factor with
# the correct global ordering, then strip the tag in the axis labels.
d_ord <- nutri_morpho_local %>%
  mutate(morpho_lab = paste(morpho, as.integer(nutri), sep = "___"))
level_order      <- d_ord %>% arrange(nutri, r) %>% pull(morpho_lab)
d_ord$morpho_lab <- factor(d_ord$morpho_lab, levels = level_order)

p_panel_B <- ggplot(d_ord, aes(x = morpho_lab, y = r)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey70", linewidth = 0.4) +
  geom_point(size = 4, colour = "grey45") +
  geom_text(aes(label = sig_local, y = r),
            vjust = -0.9, size = 4.2, fontface = "bold",
            colour = "grey25") +
  facet_wrap(~ nutri, scales = "free", ncol = 2) +
  scale_x_discrete(labels = function(x) sub("___.*", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(x = NULL, y = "Pearson's r",
       caption = paste0("BH-FDR over the 16 displayed pairs.  ",
                        "* q < 0.05    ** q < 0.01    *** q < 0.001")) +
  theme_pub() +
  theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
        axis.title.y     = element_text(size = 11),
        panel.grid.minor = element_blank(),
        strip.text       = element_text(face = "bold", size = 12),
        plot.caption     = element_text(size = 8, hjust = 0.5,
                                        face = "italic"))

save_plot(p_panel_B,
          "A12_panelB_nutri_morpho_dotplot_localFDR", w = 8, h = 6.5)


# -----------------------------------------------------------------------------
# 10. FINAL MESSAGE ----------------------------------------------------------
# -----------------------------------------------------------------------------
message("\n--- Script A v3 (streamlined) finished. ---")
message("    Outputs root:  ", OUT_ROOT)
message("    PDFs:          ", file.path(OUT_ROOT, "pdf"))
message("    PNGs:          ", file.path(OUT_ROOT, "png"))
message("    CSV / tables:  ", file.path(OUT_ROOT, "tables"))
