# =============================================================================
# Script C2 v2 (STREAMLINED) : Environmental drivers of nutritional BLUPs
# =============================================================================
# Streamlined version of `03_env_LMM_drivers_r08.R`. Same pipeline, same data,
# same collinearity filter (|r| > 0.8), same backward elimination via
# lmerTest::step(); only the effect-size plot styling is simplified and the
# table outputs are extended with four results-writing convenience CSVs.
#
# Purpose (unchanged)
# -------------------
# Identify environmental predictors of nutritional BLUPs in wild barley using
# a linear mixed-effects model with sampling site as a random intercept.
# Workflow:
#   1. Load 10 a priori predictors from Evgeny's direct field/lab measurements.
#   2. Iteratively drop predictors with pairwise |r| > 0.8 (logged).
#   3. Standardize retained predictors (mean = 0, SD = 1).
#   4. Fit full LMM per trait, then lmerTest::step() for backward elimination
#      of fixed main effects (no interactions; random intercept fixed).
#   5. Diagnostics: VIF on final model; marginal + conditional R^2 (MuMIn).
#
# Changes vs. v1 (r08)
# --------------------
#   * Effect-size plots: dropped the red-vs-grey p-value < 0.05 colour split
#     and the matching alpha encoding. Because lmerTest::step() drops
#     non-significant predictors by design, every coefficient surviving into
#     the final model is significant, so the colour split was vacuous. The
#     asterisks (* p < 0.05  ** p < 0.01  *** p < 0.001) still encode the
#     gradation of significance. A single neutral grey is used throughout.
#   * Four NEW CSV tables added for results-text / supplementary writing:
#       - C2r08_predictor_descriptive_stats.csv
#       - C2r08_coefficients_standardized_WIDE.csv
#       - C2r08_step_elimination_trace.csv
#       - C2r08_results_summary_per_trait.csv
#   * No outputs removed; no model behaviour changed.
#
# Inputs (expected in project root):
#   * BLUP_10_traits.csv
#   * merged_Evgeny_env_fitness_environmental_30_sites.csv
#
# Output root: outputs/subsection_1/C2_env_LMM_r08/
#
# Author: Shahar Liviatan
# =============================================================================


# -----------------------------------------------------------------------------
# 1. SETUP
# -----------------------------------------------------------------------------

.required <- c("dplyr", "readr", "stringr", "tidyr", "ggplot2",
               "tibble", "purrr", "lme4", "lmerTest", "MuMIn",
               "broom.mixed", "car", "performance")
.missing  <- .required[!vapply(.required, requireNamespace,
                               logical(1), quietly = TRUE)]
if (length(.missing)) {
  message("Installing missing packages: ", paste(.missing, collapse = ", "))
  install.packages(.missing)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(lme4)
  library(lmerTest)
  library(MuMIn)
  library(broom.mixed)
  library(car)
  library(performance)
})

# ---- Input files ------------------------------------------------------------

BLUP_FILE <- "BLUP_10_traits.csv"
ENV_FILE  <- "merged_Evgeny_env_fitness_environmental_30_sites.csv"

if (!file.exists(BLUP_FILE)) {
  fallback_blup <- file.path("outputs", "subsection_1", "tables", "BLUP_10_traits.csv")
  if (file.exists(fallback_blup)) BLUP_FILE <- fallback_blup
}
if (!file.exists(BLUP_FILE))
  stop("Could not find BLUP_10_traits.csv.")
if (!file.exists(ENV_FILE))
  stop("Could not find ", ENV_FILE)

# ---- Output structure -------------------------------------------------------

OUT_ROOT <- file.path("outputs", "subsection_1", "C2_env_LMM_r08")
dir.create(file.path(OUT_ROOT, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_ROOT, "models"),  recursive = TRUE, showWarnings = FALSE)

# ---- Publication theme + save helper ----------------------------------------

theme_pub <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text       = element_text(face = "bold"),
          legend.position  = "right")
}

save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT_ROOT, "figures", paste0(name, ".png")),
         p, width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(OUT_ROOT, "figures", paste0(name, ".pdf")),
         p, width = w, height = h)
}

# ---- Trait definitions ------------------------------------------------------

TRAITS_NUTRI <- c("ProteinCenterd", "StarchCenterd",
                  "BetaglucansCenterd", "FiberCenterd")
TRAIT_LABELS <- c(ProteinCenterd     = "Protein",
                  StarchCenterd      = "Starch",
                  BetaglucansCenterd = "beta-glucan",
                  FiberCenterd       = "Fiber")
TRAIT_ORDER  <- c("Protein", "Starch", "beta-glucan", "Fiber")


# -----------------------------------------------------------------------------
# 2. LOAD BLUP DATA
# -----------------------------------------------------------------------------

message("--- Loading BLUP table ---")

blup_raw <- read.csv(BLUP_FILE, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("Genotype" %in% colnames(blup_raw))
missing_traits <- setdiff(TRAITS_NUTRI, colnames(blup_raw))
if (length(missing_traits))
  stop("BLUP file is missing trait columns: ", paste(missing_traits, collapse = ", "))

blup_df <- blup_raw %>%
  dplyr::select(Genotype, all_of(TRAITS_NUTRI)) %>%
  mutate(site_id = paste0("HS", substr(Genotype, 1, 2)))

message(sprintf("    %d genotype rows; %d unique sites.",
                nrow(blup_df), length(unique(blup_df$site_id))))


# -----------------------------------------------------------------------------
# 3. LOAD ENV CSV + FEATURE SELECTION
# -----------------------------------------------------------------------------

message("--- Loading environmental table and selecting 10 predictors ---")

env_raw <- read.csv(ENV_FILE, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8")

env_raw <- env_raw %>% mutate(site_id = as.character(site_id))

# ---- Feature map: 10 biologically motivated predictors ---------------------
#   pH                   : soil pH        -> nutrient availability
#   Electrical_Cond      : salinity proxy -> ion balance, osmotic stress
#   Clay_pct             : texture        -> water retention, CEC
#   Silt_pct             : texture        -> nutrient availability
#   Sand_pct             : texture        -> drainage
#   March_Temp           : early-season temperature -> phenology
#   Precipitation        : water supply   -> grain filling
#   Organic_Carbon       : SOM indicator  -> N mineralisation
#   Total_N              : available N    -> protein synthesis
#   Elevation            : climate proxy  -> season length / temperature

FEATURE_MAP <- c(
  "pH"                      = "pH",
  "Electrical_Cond"         = "Electrical.Conductivity..mS.cm.",
  "Clay_pct"                = "Clay....",
  "Silt_pct"                = "Silt....",
  "Sand_pct"                = "Sand....",
  "March_Temp"              = "March..temp",
  "Precipitation"           = "Precipitation",
  "Organic_Carbon"          = "Organic.Carbon..mg.Kg.",
  "Total_N"                 = "Total.N..mg.Kg.",
  "Elevation"               = "Elevation"
)

# Publication-ready display labels (used in all figures and tables).
PREDICTOR_LABELS <- c(
  "pH"                      = "pH",
  "Electrical_Cond"         = "Electrical conductivity (mS/cm)",
  "Clay_pct"                = "Clay (%)",
  "Silt_pct"                = "Silt (%)",
  "Sand_pct"                = "Sand (%)",
  "March_Temp"              = "March temperature (\u00B0C)",
  "Precipitation"           = "Precipitation (mm)",
  "Organic_Carbon"          = "Organic carbon (mg/kg)",
  "Total_N"                 = "Total N (mg/kg)",
  "Elevation"               = "Elevation (m)"
)

# Verify all source columns exist.
missing_cols <- setdiff(unname(FEATURE_MAP), colnames(env_raw))
if (length(missing_cols))
  stop("Missing columns in env file: ", paste(missing_cols, collapse = ", "))

env_filtered <- env_raw[, c("site_id", unname(FEATURE_MAP))]
names(env_filtered) <- c("site_id", names(FEATURE_MAP))

# Save feature-selection log.
write.csv(
  data.frame(readable_name = names(FEATURE_MAP),
             csv_column    = unname(FEATURE_MAP),
             source        = "Evgeny direct field/lab measurement"),
  file.path(OUT_ROOT, "tables", "C2r08_features_selected.csv"),
  row.names = FALSE
)

# Force numeric + drop any non-numeric.
candidate_cols <- names(FEATURE_MAP)
env_num <- env_filtered %>%
  mutate(across(all_of(candidate_cols), ~ suppressWarnings(as.numeric(.))))

bad_cols <- vapply(env_num[candidate_cols],
                   function(v) !is.numeric(v) || all(is.na(v)), logical(1))
if (any(bad_cols)) {
  warning("Dropping non-numeric columns: ", paste(names(bad_cols)[bad_cols], collapse = ", "))
  candidate_cols <- setdiff(candidate_cols, names(bad_cols)[bad_cols])
  env_num <- env_num %>% dplyr::select(all_of(c("site_id", candidate_cols)))
}

message(sprintf("    %d candidate predictors entering |r| > 0.8 filter.",
                length(candidate_cols)))


# -----------------------------------------------------------------------------
# 4. ITERATIVE COLLINEARITY FILTER  |r| > 0.8
# -----------------------------------------------------------------------------

message("--- Iterative |r| > 0.8 collinearity filter ---")

remove_correlated_predictors <- function(df, predictors, cutoff = 0.8) {
  X <- df %>%
    dplyr::select(all_of(predictors)) %>%
    drop_na()
  if (nrow(X) < 3) stop("Too few complete site rows for correlation filter.")

  cm           <- cor(X, use = "pairwise.complete.obs")
  removal_log  <- list()
  kept         <- predictors

  repeat {
    sub_cm    <- cm[kept, kept, drop = FALSE]
    diag(sub_cm) <- 0
    if (all(abs(sub_cm) <= cutoff, na.rm = TRUE)) break

    max_idx   <- which(abs(sub_cm) == max(abs(sub_cm), na.rm = TRUE),
                       arr.ind = TRUE)[1, ]
    a         <- rownames(sub_cm)[max_idx[1]]
    b         <- colnames(sub_cm)[max_idx[2]]
    others    <- setdiff(kept, c(a, b))
    mean_abs_a <- mean(abs(sub_cm[a, others]), na.rm = TRUE)
    mean_abs_b <- mean(abs(sub_cm[b, others]), na.rm = TRUE)
    drop      <- if (is.na(mean_abs_a) || is.na(mean_abs_b)) a
                 else if (mean_abs_a >= mean_abs_b) a else b
    keep_partner <- if (drop == a) b else a

    removal_log[[length(removal_log) + 1]] <- data.frame(
      removed      = drop,
      kept_partner = keep_partner,
      r            = round(sub_cm[a, b], 3)
    )
    kept <- setdiff(kept, drop)
  }

  list(
    kept               = kept,
    removed            = setdiff(predictors, kept),
    removal_log        = if (length(removal_log)) bind_rows(removal_log)
                         else data.frame(removed = character(),
                                         kept_partner = character(),
                                         r = numeric()),
    correlation_matrix = cm
  )
}

site_only  <- env_num %>%
  dplyr::select(all_of(c("site_id", candidate_cols))) %>%
  distinct(site_id, .keep_all = TRUE)

cor_filter <- remove_correlated_predictors(
  df         = site_only,
  predictors = candidate_cols,
  cutoff     = 0.8
)

env_kept <- cor_filter$kept

# Save filter results.
write.csv(
  as.data.frame(round(cor_filter$correlation_matrix, 4)) %>%
    tibble::rownames_to_column("predictor"),
  file.path(OUT_ROOT, "tables", "C2r08_predictor_correlation_matrix.csv"),
  row.names = FALSE
)
write.csv(cor_filter$removal_log,
          file.path(OUT_ROOT, "tables", "C2r08_predictors_removed_collinearity.csv"),
          row.names = FALSE)
write.csv(data.frame(kept_predictor = env_kept,
                     display_label  = unname(PREDICTOR_LABELS[env_kept])),
          file.path(OUT_ROOT, "tables", "C2r08_predictors_kept.csv"),
          row.names = FALSE)

message(sprintf("    Removed by |r|>0.8 filter: %d  (%s)",
                length(cor_filter$removed),
                paste(cor_filter$removed, collapse = ", ")))
message(sprintf("    Retained for LMM: %d  (%s)",
                length(env_kept),
                paste(env_kept, collapse = ", ")))

# ---- NEW : Site-level descriptive stats of the 10 input predictors --------
# Useful for results-text gradient description: "Sites spanned March temps
# from X to Y degC, elevations from Z to W m, ..." -- without re-computing.

predictor_descr <- site_only %>%
  pivot_longer(-site_id, names_to = "predictor", values_to = "value") %>%
  filter(predictor %in% candidate_cols) %>%
  group_by(predictor) %>%
  summarise(
    n_sites = sum(!is.na(value)),
    mean    = round(mean(value,   na.rm = TRUE), 3),
    sd      = round(sd(value,     na.rm = TRUE), 3),
    median  = round(median(value, na.rm = TRUE), 3),
    min     = round(min(value,    na.rm = TRUE), 3),
    max     = round(max(value,    na.rm = TRUE), 3),
    range   = round(max(value,    na.rm = TRUE) -
                    min(value,    na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(
    display_label             = unname(PREDICTOR_LABELS[predictor]),
    retained_after_r08_filter = predictor %in% env_kept
  ) %>%
  dplyr::select(predictor, display_label, retained_after_r08_filter,
                n_sites, mean, sd, median, min, max, range) %>%
  arrange(desc(retained_after_r08_filter), predictor)

write.csv(predictor_descr,
          file.path(OUT_ROOT, "tables",
                    "C2r08_predictor_descriptive_stats.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. JOIN BLUP + ENV, STANDARDIZE PREDICTORS
# -----------------------------------------------------------------------------

message("--- Joining BLUP + env and standardizing retained predictors ---")

env_kept_df <- env_num %>%
  dplyr::select(all_of(c("site_id", env_kept))) %>%
  distinct(site_id, .keep_all = TRUE)

dat <- blup_df %>% inner_join(env_kept_df, by = "site_id")

n_dropped <- nrow(blup_df) - nrow(dat)
if (n_dropped > 0)
  message(sprintf("    %d genotypes dropped (no env match).", n_dropped))

message(sprintf("    Analysis frame: %d genotypes, %d sites.",
                nrow(dat), length(unique(dat$site_id))))

dat_std <- dat
dat_std[env_kept] <- scale(dat_std[env_kept], center = TRUE, scale = TRUE)

write.csv(dat,
          file.path(OUT_ROOT, "tables", "C2r08_joined_BLUP_env_rawunits.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 6. FIT FULL LMM + lmerTest::step() PER TRAIT
# -----------------------------------------------------------------------------

message("--- Fitting LMMs and running backward elimination per trait ---")

bt <- function(x) paste0("`", x, "`")

fit_one_trait <- function(trait, data_std, data_raw, predictors) {
  fml_full <- as.formula(paste(bt(trait), "~",
                               paste(bt(predictors), collapse = " + "),
                               "+ (1 | site_id)"))

  m_full <- lmer(fml_full, data = data_std, REML = TRUE)

  # Backward elimination of fixed effects; random intercept is kept throughout.
  s          <- lmerTest::step(m_full, reduce.fixed = TRUE, reduce.random = FALSE)
  m_final_std <- get_model(s)
  fe_terms    <- attr(terms(m_final_std), "term.labels")

  # Refit final formula on raw-unit data.
  if (length(fe_terms)) {
    fml_raw <- as.formula(paste(bt(trait), "~",
                                paste(bt(fe_terms), collapse = " + "),
                                "+ (1 | site_id)"))
  } else {
    fml_raw <- as.formula(paste(bt(trait), "~ 1 + (1 | site_id)"))
  }
  m_final_raw <- lmer(fml_raw, data = data_raw, REML = TRUE)

  list(trait       = trait,
       full_model  = m_full,
       step_obj    = s,
       final_std   = m_final_std,
       final_raw   = m_final_raw,
       fe_terms    = fe_terms,
       n_obs       = nrow(data_std),
       n_groups    = length(unique(data_std$site_id)))
}

model_results <- lapply(TRAITS_NUTRI, function(tr) {
  message(sprintf("    Fitting trait: %s", TRAIT_LABELS[tr]))
  fit_one_trait(tr, dat_std, dat, env_kept)
})
names(model_results) <- TRAITS_NUTRI


# -----------------------------------------------------------------------------
# 7. DIAGNOSTICS: VIF + R^2
# -----------------------------------------------------------------------------

message("--- Computing VIF and R^2 diagnostics ---")

vif_safe <- function(m) {
  fe <- attr(terms(m), "term.labels")
  if (length(fe) < 2) return(NULL)
  out <- tryCatch(car::vif(m), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  data.frame(Predictor = names(out), VIF = as.numeric(out), row.names = NULL)
}

diag_rows <- list()
vif_rows  <- list()

for (tr in TRAITS_NUTRI) {
  res <- model_results[[tr]]
  r2  <- suppressWarnings(MuMIn::r.squaredGLMM(res$final_std))
  r2m <- as.numeric(r2[1, "R2m"])
  r2c <- as.numeric(r2[1, "R2c"])

  diag_rows[[tr]] <- data.frame(
    Trait                = unname(TRAIT_LABELS[tr]),
    n_genotypes          = res$n_obs,
    n_sites              = res$n_groups,
    n_predictors_entered = length(env_kept),
    n_predictors_final   = length(res$fe_terms),
    R2_marginal          = round(r2m, 3),
    R2_conditional       = round(r2c, 3)
  )

  vf <- vif_safe(res$final_std)
  if (!is.null(vf)) {
    vf$Trait    <- unname(TRAIT_LABELS[tr])
    vf$VIF_flag <- dplyr::case_when(
      vf$VIF >= 10 ~ "VIF \u2265 10 (serious)",
      vf$VIF >  5  ~ "VIF > 5 (review)",
      TRUE         ~ "OK"
    )
    vif_rows[[tr]] <- vf
  }
}

diag_df <- bind_rows(diag_rows)
vif_df  <- if (length(vif_rows)) bind_rows(vif_rows) else data.frame()

write.csv(diag_df,
          file.path(OUT_ROOT, "tables", "C2r08_model_fit_diagnostics.csv"),
          row.names = FALSE)
write.csv(vif_df,
          file.path(OUT_ROOT, "tables", "C2r08_VIF_final_models.csv"),
          row.names = FALSE)

message("    R\u00B2 diagnostics:")
for (i in seq_len(nrow(diag_df))) {
  message(sprintf("      %-12s  predictors: %d/%d  R2m=%.3f  R2c=%.3f",
                  diag_df$Trait[i],
                  diag_df$n_predictors_final[i],
                  diag_df$n_predictors_entered[i],
                  diag_df$R2_marginal[i],
                  diag_df$R2_conditional[i]))
}

if (nrow(vif_df)) {
  bad_vif <- vif_df %>% dplyr::filter(VIF_flag != "OK")
  if (nrow(bad_vif)) {
    message("    *** VIF WARNINGS ***")
    for (i in seq_len(nrow(bad_vif)))
      message(sprintf("      %s | %s | VIF=%.2f",
                      bad_vif$Trait[i], bad_vif$Predictor[i], bad_vif$VIF[i]))
  } else {
    message("    All retained predictors: VIF \u2264 5.")
  }
}


# -----------------------------------------------------------------------------
# 8. EXTRACT COEFFICIENTS (long format)
# -----------------------------------------------------------------------------

message("--- Extracting standardized and raw-unit coefficients ---")

tidy_fixed <- function(m, trait, kind) {
  td <- broom.mixed::tidy(m, effects = "fixed", conf.int = TRUE,
                          conf.method = "Wald")
  td %>%
    dplyr::filter(term != "(Intercept)") %>%
    dplyr::mutate(
      Trait      = unname(TRAIT_LABELS[trait]),
      TraitLabel = factor(Trait, levels = TRAIT_ORDER),
      Kind       = kind,
      PredLabel  = ifelse(term %in% names(PREDICTOR_LABELS),
                          PREDICTOR_LABELS[term], term),
      p_stars    = dplyr::case_when(
        is.na(p.value)  ~ "",
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE            ~ ""
      )
    ) %>%
    dplyr::select(Trait, TraitLabel, Kind,
                  Predictor     = term,
                  PredLabel,
                  estimate, std.error,
                  ci_low  = conf.low,
                  ci_high = conf.high,
                  df,
                  p_value = p.value,
                  p_stars)
}

coef_std <- bind_rows(lapply(TRAITS_NUTRI, function(tr) {
  res <- model_results[[tr]]
  if (length(res$fe_terms)) tidy_fixed(res$final_std, tr, "standardized")
  else data.frame()
}))

coef_raw <- bind_rows(lapply(TRAITS_NUTRI, function(tr) {
  res <- model_results[[tr]]
  if (length(res$fe_terms)) tidy_fixed(res$final_raw, tr, "raw-unit")
  else data.frame()
}))

write.csv(coef_std,
          file.path(OUT_ROOT, "tables", "C2r08_coefficients_standardized.csv"),
          row.names = FALSE)
write.csv(coef_raw,
          file.path(OUT_ROOT, "tables", "C2r08_coefficients_raw_units.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 8b. NEW : COEFFICIENTS WIDE (predictors x traits) -- manuscript paste-ready
# -----------------------------------------------------------------------------
# Rows  = every retained predictor (env_kept); columns = the 4 traits.
# Cell  = "beta [CI_low, CI_high] stars".  NA  when the predictor was
# dropped during backward elimination for that trait.

message("--- Building wide standardized-coefficients table (predictors x traits) ---")

format_cell <- function(beta, lo, hi, stars) {
  if (any(is.na(c(beta, lo, hi)))) return(NA_character_)
  sprintf("%.3f [%.3f, %.3f]%s", beta, lo, hi, ifelse(is.na(stars), "", stars))
}

cells_by_trait <- coef_std %>%
  rowwise() %>%
  mutate(cell = format_cell(estimate, ci_low, ci_high, p_stars)) %>%
  ungroup() %>%
  dplyr::select(Predictor, Trait, cell) %>%
  pivot_wider(names_from = Trait, values_from = cell)

all_pred_rows <- data.frame(
  Predictor     = env_kept,
  display_label = unname(PREDICTOR_LABELS[env_kept]),
  stringsAsFactors = FALSE
)

coef_wide_std <- all_pred_rows %>%
  left_join(cells_by_trait, by = "Predictor")

# Guarantee all four trait columns are present (even if a trait kept nothing).
for (tc in TRAIT_ORDER) {
  if (!tc %in% colnames(coef_wide_std)) {
    coef_wide_std[[tc]] <- NA_character_
  }
}
coef_wide_std <- coef_wide_std[, c("Predictor", "display_label", TRAIT_ORDER)]

write.csv(coef_wide_std,
          file.path(OUT_ROOT, "tables",
                    "C2r08_coefficients_standardized_WIDE.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 8c. NEW : STEP() ELIMINATION TRACE -- documents the model-selection process
# -----------------------------------------------------------------------------
# One row per (Trait, Predictor) showing whether the predictor was eliminated
# during backward selection and at which step. Extracted from lmerTest's
# step()$fixed.

message("--- Extracting backward-elimination trace per trait ---")

extract_step_trace <- function(step_obj, trait_label) {
  s_fixed <- tryCatch(as.data.frame(step_obj$fixed),
                      error = function(e) NULL)
  if (is.null(s_fixed) || nrow(s_fixed) == 0) return(NULL)

  s_fixed$Predictor <- rownames(s_fixed)

  out <- data.frame(
    Trait                = trait_label,
    Predictor            = s_fixed$Predictor,
    step_when_eliminated = if ("Eliminated" %in% colnames(s_fixed))
                             as.numeric(s_fixed[["Eliminated"]]) else NA_real_,
    F_statistic          = if ("F value"    %in% colnames(s_fixed))
                             as.numeric(s_fixed[["F value"]])    else NA_real_,
    DenDF                = if ("DenDF"      %in% colnames(s_fixed))
                             as.numeric(s_fixed[["DenDF"]])      else NA_real_,
    p_value              = if ("Pr(>F)"     %in% colnames(s_fixed))
                             as.numeric(s_fixed[["Pr(>F)"]])     else NA_real_,
    stringsAsFactors     = FALSE
  )
  # Eliminated == 0 means "retained in the final model".
  out$retained_in_final <- !is.na(out$step_when_eliminated) &
                           out$step_when_eliminated == 0
  out
}

step_trace_df <- bind_rows(lapply(TRAITS_NUTRI, function(tr) {
  res <- model_results[[tr]]
  extract_step_trace(res$step_obj, unname(TRAIT_LABELS[tr]))
}))

write.csv(step_trace_df,
          file.path(OUT_ROOT, "tables",
                    "C2r08_step_elimination_trace.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 8d. NEW : RESULTS SUMMARY PER TRAIT -- headline numbers for the manuscript
# -----------------------------------------------------------------------------
# One row per trait. Top positive and top negative predictor in the final
# model (by signed beta), plus R^2 and the mean |beta|.

message("--- Building per-trait results summary ---")

results_summary_rows <- lapply(TRAITS_NUTRI, function(tr) {
  res <- model_results[[tr]]
  trait_lbl <- unname(TRAIT_LABELS[tr])
  d   <- diag_df %>% dplyr::filter(Trait == trait_lbl)

  trait_coef <- coef_std %>% dplyr::filter(Trait == trait_lbl)

  if (nrow(trait_coef)) {
    top_pos <- trait_coef %>% dplyr::filter(estimate > 0) %>%
                              dplyr::arrange(desc(estimate)) %>% dplyr::slice(1)
    top_neg <- trait_coef %>% dplyr::filter(estimate < 0) %>%
                              dplyr::arrange(estimate)      %>% dplyr::slice(1)
    mean_abs <- mean(abs(trait_coef$estimate))
  } else {
    top_pos <- data.frame(PredLabel = NA_character_, estimate = NA_real_,
                          p_value = NA_real_, p_stars = NA_character_)
    top_neg <- top_pos
    mean_abs <- NA_real_
  }

  data.frame(
    Trait                  = trait_lbl,
    n_predictors_final     = d$n_predictors_final[1],
    R2_marginal            = d$R2_marginal[1],
    R2_conditional         = d$R2_conditional[1],
    top_positive_predictor = if (nrow(top_pos)) top_pos$PredLabel[1] else NA_character_,
    top_positive_beta      = if (nrow(top_pos)) round(top_pos$estimate[1], 3) else NA_real_,
    top_positive_p         = if (nrow(top_pos)) signif(top_pos$p_value[1], 3) else NA_real_,
    top_positive_stars     = if (nrow(top_pos)) top_pos$p_stars[1] else NA_character_,
    top_negative_predictor = if (nrow(top_neg)) top_neg$PredLabel[1] else NA_character_,
    top_negative_beta      = if (nrow(top_neg)) round(top_neg$estimate[1], 3) else NA_real_,
    top_negative_p         = if (nrow(top_neg)) signif(top_neg$p_value[1], 3) else NA_real_,
    top_negative_stars     = if (nrow(top_neg)) top_neg$p_stars[1] else NA_character_,
    mean_abs_beta          = round(mean_abs, 3),
    stringsAsFactors       = FALSE
  )
})

results_summary_df <- bind_rows(results_summary_rows)
write.csv(results_summary_df,
          file.path(OUT_ROOT, "tables",
                    "C2r08_results_summary_per_trait.csv"),
          row.names = FALSE)


# -----------------------------------------------------------------------------
# 9. EFFECT PLOTS  (colour-cleaned : single neutral grey, asterisks encode sig.)
# -----------------------------------------------------------------------------

message("--- Building effect-size plots ---")

EFFECT_PLOT_COLOUR <- "grey25"   # single neutral dark tone

build_effect_plot <- function(coef_df, x_axis_label) {
  if (!nrow(coef_df)) {
    message("    (no retained predictors -- skipping plot)")
    return(NULL)
  }

  # Order predictors within each facet by signed effect size.
  coef_df <- coef_df %>%
    dplyr::group_by(TraitLabel) %>%
    dplyr::mutate(PredOrder = reorder(PredLabel, estimate)) %>%
    dplyr::ungroup()

  ggplot(coef_df,
         aes(x = estimate, y = PredOrder)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               linewidth = 0.4, colour = "grey45") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                   height = 0.26, linewidth = 0.85,
                   colour = EFFECT_PLOT_COLOUR) +
    geom_point(size = 3.4, colour = EFFECT_PLOT_COLOUR) +
    geom_text(aes(label = p_stars,
                  hjust = ifelse(estimate >= 0, -0.4, 1.4)),
              size = 7.5, vjust = 0.35, colour = "black") +
    facet_wrap(~ TraitLabel, scales = "free", ncol = 1) +
    scale_x_continuous(expand = expansion(mult = 0.13)) +  # room so * stars don't clip the panel edge
    # Shorten y-axis labels: drop trailing units in parentheses (e.g. "Silt (%)" -> "Silt")
    scale_y_discrete(labels = function(x) sub("\\s*\\(.*\\)\\s*$", "", x)) +
    labs(x = x_axis_label, y = NULL) +
    theme_pub(base_size = 19) +
    theme(panel.grid.major.y = element_blank(),
          legend.position    = "none",
          strip.text         = element_text(face = "bold", size = 22),
          axis.title.x       = element_text(size = 22),
          axis.text.y        = element_text(size = 20, face = "bold"),
          axis.text.x        = element_text(size = 22))
}

# Use a plotmath expression for β so it renders as the Greek symbol in both the
# cairo PDF and PNG (a literal "β" string can drop to ".." after PDF conversion).
p_std <- build_effect_plot(coef_std, expression("Standardized effect (" * beta * ")"))
p_raw <- build_effect_plot(coef_raw, "Estimated effect (raw units)")

if (!is.null(p_std)) save_plot(p_std, "C2r08_effects_standardized", w = 7.5, h = 9)
if (!is.null(p_raw)) save_plot(p_raw, "C2r08_effects_raw_units",    w = 7.5, h = 9)


# -----------------------------------------------------------------------------
# 10. FINAL-FORMULA TABLE + MODEL SUMMARIES
# -----------------------------------------------------------------------------

message("--- Writing final-formula table and model summaries ---")

formula_rows <- lapply(TRAITS_NUTRI, function(tr) {
  res <- model_results[[tr]]
  data.frame(
    Trait             = unname(TRAIT_LABELS[tr]),
    Final_Formula     = if (length(res$fe_terms))
                          paste0(tr, " ~ ",
                                 paste(res$fe_terms, collapse = " + "),
                                 " + (1 | site_id)")
                        else paste0(tr, " ~ 1 + (1 | site_id)"),
    n_predictors_final = length(res$fe_terms)
  )
})
write.csv(bind_rows(formula_rows),
          file.path(OUT_ROOT, "tables", "C2r08_final_formulas.csv"),
          row.names = FALSE)

sink(file.path(OUT_ROOT, "models", "C2r08_LMM_model_summaries.txt"))
cat("LMM environmental drivers of nutritional BLUPs (Script C2 v2 streamlined)\n")
cat(strrep("=", 64), "\n\n", sep = "")
cat("BLUP file:        ", BLUP_FILE, "\n", sep = "")
cat("Env  file:        ", ENV_FILE,  "\n", sep = "")
cat("Random structure: (1 | site_id), held fixed throughout step()\n")
cat("Selection:        lmerTest::step(reduce.fixed=TRUE, reduce.random=FALSE)\n")
cat("Collinearity:     iterative |r| > 0.8 filter (", length(cor_filter$removed),
    " removed)\n\n", sep = "")
cat("10 input predictors: ", paste(candidate_cols, collapse = ", "), "\n\n")
cat("Removed by filter (", length(cor_filter$removed), "): ",
    paste(cor_filter$removed, collapse = ", "), "\n\n", sep = "")
cat("Retained for LMM (", length(env_kept), "): ",
    paste(env_kept, collapse = ", "), "\n\n", sep = "")

for (tr in TRAITS_NUTRI) {
  res <- model_results[[tr]]
  cat("\n", strrep("=", 64), "\n", sep = "")
  cat("Trait: ", unname(TRAIT_LABELS[tr]), "  (", tr, ")\n", sep = "")
  cat("Final predictors: ",
      if (length(res$fe_terms)) paste(res$fe_terms, collapse = ", ")
      else "intercept only", "\n", sep = "")
  cat("\n--- Standardized model summary ---\n")
  print(summary(res$final_std))
  cat("\n--- Raw-unit model summary ---\n")
  print(summary(res$final_raw))
  cat("\n--- step() decision log ---\n")
  print(res$step_obj)
}
sink()


# -----------------------------------------------------------------------------
# 11. PREDICTOR CORRELATION HEATMAP (kept = green labels, removed = grey)
# -----------------------------------------------------------------------------

message("--- Building predictor correlation heatmap ---")

cm_mat  <- cor_filter$correlation_matrix
cm_long <- expand.grid(Var1 = rownames(cm_mat),
                       Var2 = colnames(cm_mat),
                       stringsAsFactors = FALSE) %>%
  dplyr::mutate(r = mapply(function(a, b) cm_mat[a, b], Var1, Var2))

ax_cols   <- ifelse(rownames(cm_mat) %in% env_kept, "#1B7837", "grey50")

cm_long <- cm_long %>%
  dplyr::mutate(
    Var1_lab = ifelse(Var1 %in% names(PREDICTOR_LABELS),
                      PREDICTOR_LABELS[Var1], Var1),
    Var2_lab = ifelse(Var2 %in% names(PREDICTOR_LABELS),
                      PREDICTOR_LABELS[Var2], Var2)
  )

ax_order_lab <- ifelse(rownames(cm_mat) %in% names(PREDICTOR_LABELS),
                       PREDICTOR_LABELS[rownames(cm_mat)],
                       rownames(cm_mat))

p_cor <- ggplot(cm_long, aes(x = Var2_lab, y = Var1_lab, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4.0, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  scale_x_discrete(limits = ax_order_lab) +
  scale_y_discrete(limits = rev(ax_order_lab)) +
  coord_equal() +
  labs(x = NULL, y = NULL,
       caption = "Green labels = retained after |r|>0.8 filter; grey = removed.") +
  theme_pub(base_size = 12) +
  theme(panel.grid   = element_blank(),
        axis.text.x  = element_text(angle = 55, hjust = 1, colour = ax_cols),
        axis.text.y  = element_text(colour = rev(ax_cols)),
        plot.caption = element_text(hjust = 0.5, size = 12, face = "italic"))

save_plot(p_cor, "C2r08_predictor_correlation_heatmap",
          w = max(7, 0.32 * ncol(cm_mat) + 2.5),
          h = max(7, 0.32 * nrow(cm_mat) + 2.5))


# -----------------------------------------------------------------------------
# 12. FINAL MESSAGE
# -----------------------------------------------------------------------------

message("\n==============================================================")
message("Script C2 v2 (streamlined) finished. Outputs under: ", OUT_ROOT)
message("==============================================================")
message("Predictors removed by |r|>0.8: ", paste(cor_filter$removed, collapse = ", "))
message("Predictors in LMM:             ", paste(env_kept, collapse = ", "))
message("\nFigures:")
message("  * C2r08_effects_standardized        (main)")
message("  * C2r08_effects_raw_units           (supplementary)")
message("  * C2r08_predictor_correlation_heatmap")
message("\nTables (existing):")
message("  * C2r08_features_selected.csv")
message("  * C2r08_predictor_correlation_matrix.csv")
message("  * C2r08_predictors_removed_collinearity.csv")
message("  * C2r08_predictors_kept.csv")
message("  * C2r08_joined_BLUP_env_rawunits.csv")
message("  * C2r08_final_formulas.csv")
message("  * C2r08_coefficients_standardized.csv")
message("  * C2r08_coefficients_raw_units.csv")
message("  * C2r08_VIF_final_models.csv")
message("  * C2r08_model_fit_diagnostics.csv")
message("\nTables (NEW for results-writing):")
message("  * C2r08_predictor_descriptive_stats.csv")
message("  * C2r08_coefficients_standardized_WIDE.csv")
message("  * C2r08_step_elimination_trace.csv")
message("  * C2r08_results_summary_per_trait.csv")
message("\nText: ", file.path(OUT_ROOT, "models", "C2r08_LMM_model_summaries.txt"))
