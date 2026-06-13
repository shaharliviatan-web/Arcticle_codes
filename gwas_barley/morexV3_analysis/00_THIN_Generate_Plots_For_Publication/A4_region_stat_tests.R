# =============================================================================
# Section A4 (add-on) : Statistical test of regional differences
# =============================================================================
# Tests whether the four grain nutritional traits differ across the six
# ecological regions, using the centered BLUP scale. This complements the
# already-existing descriptive summaries (A4_region_summary_*.csv) WITHOUT
# regenerating any of them.
#
# Per nutritional trait (Protein, Starch, beta-glucan, Fiber):
#   * Kruskal-Wallis test of value_centered across the six regions
#       -> H (chi-squared), df, p
#   * Epsilon-squared effect size  ->  eps^2 = H / (N - 1)
#   * Dunn's post-hoc, Benjamini-Hochberg corrected (all pairwise contrasts;
#       North vs Desert flagged for the results text)
#
# Inputs  (read-only, NOT modified) :
#   outputs/subsection_1/tables/A4_site_boxplot_data.csv
# Outputs (new files only) :
#   outputs/subsection_1/tables/A4_region_KW_tests.csv
#   outputs/subsection_1/tables/A4_region_dunn_posthoc.csv
#
# Author : Shahar Liviatan
# =============================================================================

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

suppressPackageStartupMessages({
  library(dplyr)
  library(rstatix)
})

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
dat$Trait  <- factor(dat$Trait,  levels = trait_order)

stopifnot(!any(is.na(dat$region)), !any(is.na(dat$Trait)))

# ---- Per-trait Kruskal-Wallis + epsilon-squared + Dunn post-hoc ------------
kw_rows   <- list()
dunn_rows <- list()

for (tr in trait_order) {
  d <- dat %>% filter(Trait == tr) %>%
    select(region, value_centered) %>%
    filter(!is.na(value_centered))
  d$region <- droplevels(d$region)

  N <- nrow(d)
  k <- nlevels(d$region)

  # Kruskal-Wallis
  kw <- kruskal.test(value_centered ~ region, data = d)
  H  <- unname(kw$statistic)
  df <- unname(kw$parameter)
  p  <- kw$p.value

  # Epsilon-squared effect size:  eps^2 = H / (N - 1)
  eps2 <- H / (N - 1)

  kw_rows[[tr]] <- data.frame(
    trait           = tr,
    H               = round(H, 4),
    df              = df,
    p               = p,
    epsilon_squared = round(eps2, 4),
    n               = N,
    stringsAsFactors = FALSE
  )

  # Dunn's post-hoc, Benjamini-Hochberg corrected (all pairwise contrasts).
  dn <- d %>%
    dunn_test(value_centered ~ region, p.adjust.method = "BH") %>%
    as.data.frame()

  dn_out <- data.frame(
    trait        = tr,
    group1       = dn$group1,
    group2       = dn$group2,
    n1           = dn$n1,
    n2           = dn$n2,
    Z            = round(dn$statistic, 4),
    p            = dn$p,
    p_adj_BH     = dn$p.adj,
    p_adj_signif = dn$p.adj.signif,
    stringsAsFactors = FALSE
  )
  # Flag the focal North vs Desert contrast.
  dn_out$highlight <- ifelse(
    (dn_out$group1 == "North" & dn_out$group2 == "Desert") |
    (dn_out$group1 == "Desert" & dn_out$group2 == "North"),
    "North_vs_Desert", "")

  dunn_rows[[tr]] <- dn_out
}

kw_tab   <- do.call(rbind, kw_rows);   rownames(kw_tab)   <- NULL
dunn_tab <- do.call(rbind, dunn_rows); rownames(dunn_tab) <- NULL

# ---- Write outputs (new files only) ----------------------------------------
kw_path   <- file.path(TABLES_DIR, "A4_region_KW_tests.csv")
dunn_path <- file.path(TABLES_DIR, "A4_region_dunn_posthoc.csv")

write.csv(kw_tab,   kw_path,   row.names = FALSE)
write.csv(dunn_tab, dunn_path, row.names = FALSE)

# ---- Console report --------------------------------------------------------
cat("\n=== A4_region_KW_tests.csv ===\n")
print(kw_tab, row.names = FALSE)

cat("\n=== A4_region_dunn_posthoc.csv ===\n")
print(dunn_tab, row.names = FALSE)

cat("\n--- North vs Desert (focal contrast) ---\n")
print(dunn_tab[dunn_tab$highlight == "North_vs_Desert",
               c("trait", "group1", "group2", "Z", "p", "p_adj_BH", "p_adj_signif")],
      row.names = FALSE)

cat("\nWrote:\n  ", kw_path, "\n  ", dunn_path, "\n", sep = "")
