# =============================================================================
# Section A4 (add-on) : Statistical test of regional differences
# =============================================================================
# Tests whether the four grain nutritional traits differ across the six
# ecological regions, using the centered BLUP scale. This complements the
# already-existing descriptive summaries (A4_region_summary_*.csv) WITHOUT
# regenerating any of them.
#
# RATIONALE (replaces the previous Kruskal-Wallis / Dunn analysis):
#   Accessions are clustered within sampling Sites (diagnostic ICC ~0.25-0.61),
#   so accession-level observations are NOT independent. A one-way test
#   (Kruskal-Wallis or plain ANOVA) ignores this clustering and pseudoreplicates.
#   We therefore model Region as a fixed effect and Site as a random intercept:
#
#       value_centered ~ Region + (1 | Site)
#
#   Diagnostics (non-singular fits, approximately normal residuals, homogeneous
#   residual variance across Regions once Site is modelled) support a standard
#   homoscedastic linear mixed model for all four traits.
#
# Per nutritional trait (Protein, Starch, beta-glucan, Fiber):
#   * Linear mixed model  value_centered ~ Region + (1 | Site)  [lmerTest]
#   * Overall Region effect : Type III F-test with Satterthwaite df -> F, df, p
#   * Site variance + ICC   : clustering of accessions within sampling sites
#   * Estimated marginal means (EMMs) per Region [emmeans]
#   * THREE pre-specified post-hoc contrasts only (ecotone groups excluded):
#       North vs Desert, North vs Coast, Coast vs Desert
#     Benjamini-Hochberg (FDR) corrected across these three contrasts, within
#     each trait.
#       (North vs Desert flagged for the results text)
#
# Inputs  (read-only, NOT modified) :
#   outputs/subsection_1/tables/A4_site_boxplot_data.csv
# Outputs (new files only) :
#   outputs/subsection_1/tables/A4_region_LMM_overall.csv
#   outputs/subsection_1/tables/A4_region_LMM_emmeans.csv
#   outputs/subsection_1/tables/A4_region_LMM_pairwise.csv
#
# Author : Shahar Liviatan
# =============================================================================

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

# emmeans was installed into a project-local library (avoids writing under ~).
.libPaths(c("/mnt/data/shahar/Rlibs", .libPaths()))

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)   # adds Satterthwaite F-tests / df to lmer
  library(emmeans)
})
emm_options(lmerTest.limit = 5000, pbkrtest.limit = 5000)

# ---- Configuration (mirrors 01_phenotypic_analysis_no_GxE_v3_streamlined.R) -
TABLES_DIR  <- file.path("outputs", "subsection_1", "tables")
BOXPLOT_CSV <- file.path(TABLES_DIR, "A4_site_boxplot_data.csv")

# Region factor levels and trait order, identical to the generating script.
REGION_ORDER <- c("North", "Coast", "Desert",
                  "HZ1 (North-Coast)", "HZ2 (North-Desert)", "HZ3 (Coast-Desert)")
trait_order  <- c("Protein", "Starch", "β-glucan", "Fiber")

# ---- Load accession-level data ---------------------------------------------
dat <- read.csv(BOXPLOT_CSV, check.names = FALSE, stringsAsFactors = FALSE)
dat$region <- factor(dat$region, levels = REGION_ORDER)
dat$site   <- factor(dat$site)
dat$Trait  <- factor(dat$Trait,  levels = trait_order)

stopifnot(!any(is.na(dat$region)), !any(is.na(dat$Trait)))

# Significance-star helper (same thresholds as the old Dunn table).
star_fn <- function(p) {
  ifelse(is.na(p), "ns",
  ifelse(p < 0.0001, "****",
  ifelse(p < 0.001,  "***",
  ifelse(p < 0.01,   "**",
  ifelse(p < 0.05,   "*", "ns")))))
}

# ---- Per-trait LMM : overall Region test + EMMs + pairwise contrasts -------
overall_rows  <- list()
emm_rows       <- list()
pairwise_rows  <- list()

for (tr in trait_order) {
  d <- dat %>% filter(Trait == tr) %>%
    select(region, site, value_centered) %>%
    filter(!is.na(value_centered))
  d$region <- droplevels(d$region)
  d$site   <- droplevels(d$site)

  N <- nrow(d)

  # ---- Fit:  value_centered ~ Region + (1 | Site) --------------------------
  fit <- lmer(value_centered ~ region + (1 | site), data = d, REML = TRUE)

  # Overall Region effect: Type III F-test, Satterthwaite denominator df.
  aov_tab <- anova(fit, type = 3)
  F_val   <- aov_tab[["F value"]][1]
  num_df  <- aov_tab[["NumDF"]][1]
  den_df  <- aov_tab[["DenDF"]][1]
  p_overall <- aov_tab[["Pr(>F)"]][1]

  # Variance components -> Site ICC (clustering of accessions within sites).
  vc      <- as.data.frame(VarCorr(fit))
  v_site  <- vc$vcov[vc$grp == "site"]
  v_resid <- vc$vcov[vc$grp == "Residual"]
  icc     <- v_site / (v_site + v_resid)
  singular <- isSingular(fit, tol = 1e-4)

  overall_rows[[tr]] <- data.frame(
    trait        = tr,
    n_accessions = N,
    n_sites      = nlevels(d$site),
    F_value      = round(F_val, 4),
    NumDF        = round(num_df, 0),
    DenDF        = round(den_df, 2),
    p_region     = p_overall,
    p_signif     = star_fn(p_overall),
    var_site     = round(v_site, 4),
    var_resid    = round(v_resid, 4),
    ICC_site     = round(icc, 4),
    singular     = singular,
    stringsAsFactors = FALSE
  )

  # ---- Estimated marginal means per Region ---------------------------------
  emm <- emmeans(fit, ~ region)
  emm_df <- as.data.frame(emm)
  emm_rows[[tr]] <- data.frame(
    trait    = tr,
    region   = as.character(emm_df$region),
    emmean   = round(emm_df$emmean, 4),
    SE       = round(emm_df$SE, 4),
    df       = round(emm_df$df, 2),
    lower_CL = round(emm_df$lower.CL, 4),
    upper_CL = round(emm_df$upper.CL, 4),
    stringsAsFactors = FALSE
  )

  # ---- Three pre-specified contrasts only, BH (FDR)-corrected --------------
  # Ecotone groups (HZ1/HZ2/HZ3) are intentionally excluded from the post-hoc.
  # Raw p-values come from the model contrasts (adjust = "none"); Benjamini-
  # Hochberg FDR is then applied across exactly these three contrasts, per trait.
  emm_levels <- levels(d$region)
  ct_list <- list(
    "North vs Desert" = as.integer(emm_levels == "North") -
                        as.integer(emm_levels == "Desert"),
    "North vs Coast"  = as.integer(emm_levels == "North") -
                        as.integer(emm_levels == "Coast"),
    "Coast vs Desert" = as.integer(emm_levels == "Coast") -
                        as.integer(emm_levels == "Desert")
  )
  pw <- as.data.frame(contrast(emm, method = ct_list, adjust = "none"))

  g1 <- sub(" vs .*$", "", as.character(pw$contrast))
  g2 <- sub("^.* vs ", "", as.character(pw$contrast))

  p_bh <- p.adjust(pw$p.value, method = "BH")   # Benjamini-Hochberg across 3 contrasts

  pw_out <- data.frame(
    trait          = tr,
    group1         = g1,
    group2         = g2,
    estimate       = round(pw$estimate, 4),  # EMM difference (group1 - group2)
    SE             = round(pw$SE, 4),
    df             = round(pw$df, 2),
    t_ratio        = round(pw$t.ratio, 4),
    p_raw          = pw$p.value,
    p_adj_BH       = p_bh,
    p_adj_signif   = star_fn(p_bh),
    stringsAsFactors = FALSE
  )
  pw_out$highlight <- ifelse(pw_out$group1 == "North" & pw_out$group2 == "Desert",
                             "North_vs_Desert", "")

  pairwise_rows[[tr]] <- pw_out
}

overall_tab  <- do.call(rbind, overall_rows);  rownames(overall_tab)  <- NULL
emm_tab      <- do.call(rbind, emm_rows);       rownames(emm_tab)      <- NULL
pairwise_tab <- do.call(rbind, pairwise_rows);  rownames(pairwise_tab) <- NULL

# ---- Write outputs (new files only) ----------------------------------------
overall_path  <- file.path(TABLES_DIR, "A4_region_LMM_overall.csv")
emm_path      <- file.path(TABLES_DIR, "A4_region_LMM_emmeans.csv")
pairwise_path <- file.path(TABLES_DIR, "A4_region_LMM_pairwise.csv")

write.csv(overall_tab,  overall_path,  row.names = FALSE)
write.csv(emm_tab,      emm_path,      row.names = FALSE)
write.csv(pairwise_tab, pairwise_path, row.names = FALSE)

# ---- Console report --------------------------------------------------------
cat("\n=== A4_region_LMM_overall.csv  (overall Region effect per trait) ===\n")
print(overall_tab, row.names = FALSE)

cat("\n=== A4_region_LMM_emmeans.csv  (estimated marginal means per Region) ===\n")
print(emm_tab, row.names = FALSE)

cat("\n=== A4_region_LMM_pairwise.csv  (3 pre-specified contrasts, BH/FDR) ===\n")
print(pairwise_tab, row.names = FALSE)

cat("\n--- North vs Desert (focal contrast) ---\n")
print(pairwise_tab[pairwise_tab$highlight == "North_vs_Desert",
                   c("trait", "group1", "group2", "estimate",
                     "t_ratio", "p_raw", "p_adj_BH", "p_adj_signif")],
      row.names = FALSE)

cat("\nWrote:\n  ", overall_path,
    "\n  ", emm_path,
    "\n  ", pairwise_path, "\n", sep = "")
