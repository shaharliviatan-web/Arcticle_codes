# =============================================================================
# Script C2 — STAGE 1 : Collinearity among the ten environmental predictors
# =============================================================================
# Standalone, self-contained FIRST stage of the revised environmental analysis.
#
# Scope (deliberately limited):
#   * Assess Pearson collinearity among the 10 a priori environmental predictors
#     using site-level observations only.
#   * NO nutritional-trait / BLUP data are loaded or used.
#   * NO FDR correction, NO regression, NO predictor selection, NO automatic
#     variable removal. Which collinear variable to keep is a later, manual,
#     biology-driven decision and is NOT made here.
#
# Input:
#   merged_Evgeny_env_fitness_environmental_30_sites.csv   (one row per site)
#
# Output root: outputs/subsection_1/C2_environmental_predictor_collinearity/
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
})

# ---- Input ------------------------------------------------------------------

ENV_FILE <- "merged_Evgeny_env_fitness_environmental_30_sites.csv"
if (!file.exists(ENV_FILE)) stop("Could not find ", ENV_FILE)

# ---- Output structure -------------------------------------------------------

OUT_ROOT <- file.path("outputs", "subsection_1",
                      "C2_environmental_predictor_collinearity")
DIR_TAB  <- file.path(OUT_ROOT, "tables")
DIR_FIG  <- file.path(OUT_ROOT, "figures")
DIR_REP  <- file.path(OUT_ROOT, "reports")
for (d in c(DIR_TAB, DIR_FIG, DIR_REP))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Publication theme + save helper (matches scripts 01/02) ----------------

theme_pub <- function(base_size = 14) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "black"),
      panel.grid       = element_blank(),
      panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1)
    )
}

save_plot <- function(plot, name, w, h) {
  ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), plot, width = w, height = h,
         device = cairo_pdf)
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), plot, width = w, height = h,
         dpi = 600, bg = "white")
  invisible(plot)
}

# ---- The ten predictors: internal name -> source column / label / units -----
#   (identical source columns + display labels to the existing env scripts)

FEATURE_MAP <- c(
  "pH"              = "pH",
  "Electrical_Cond" = "Electrical.Conductivity..mS.cm.",
  "Clay_pct"        = "Clay....",
  "Silt_pct"        = "Silt....",
  "Sand_pct"        = "Sand....",
  "March_Temp"      = "March..temp",
  "Precipitation"   = "Precipitation",
  "Organic_Carbon"  = "Organic.Carbon..mg.Kg.",
  "Total_N"         = "Total.N..mg.Kg.",
  "Elevation"       = "Elevation"
)
PREDICTOR_LABELS <- c(
  "pH"              = "pH",
  "Electrical_Cond" = "Electrical conductivity (mS/cm)",
  "Clay_pct"        = "Clay (%)",
  "Silt_pct"        = "Silt (%)",
  "Sand_pct"        = "Sand (%)",
  "March_Temp"      = "March temperature (°C)",
  "Precipitation"   = "Precipitation (mm)",
  "Organic_Carbon"  = "Organic carbon (mg/kg)",
  "Total_N"         = "Total N (mg/kg)",
  "Elevation"       = "Elevation (m)"
)
PREDICTOR_UNITS <- c(
  "pH"              = "pH units (unitless)",
  "Electrical_Cond" = "mS/cm",
  "Clay_pct"        = "%",
  "Silt_pct"        = "%",
  "Sand_pct"        = "%",
  "March_Temp"      = "°C",
  "Precipitation"   = "mm",
  "Organic_Carbon"  = "mg/kg",
  "Total_N"         = "mg/kg",
  "Elevation"       = "m"
)
PREDICTORS <- names(FEATURE_MAP)
COLLIN_CUT <- 0.8

# Site 04 = "Hachola" (Hula): its accessions were missing from the phenotyping
# experiments and were excluded from ALL analyses (see Methods). The env file
# still lists it (30 sites), so we drop it here to use the SAME 29 sites that
# carry BLUP data and that every other analysis uses.
EXCLUDE_SITES <- c("HS04")

# -----------------------------------------------------------------------------
# 1. LOAD + VERIFY ONE ROW PER SITE
# -----------------------------------------------------------------------------

message("--- Loading environmental table ---")

env_raw <- read.csv(ENV_FILE, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8") %>%
  mutate(site_id = as.character(site_id))

dups <- env_raw$site_id[duplicated(env_raw$site_id)]
site_unique <- length(dups) == 0
if (!site_unique)
  warning("Duplicate site_id values found: ", paste(unique(dups), collapse = ", "))

n_sites_total <- length(unique(env_raw$site_id))
message(sprintf("    Sites in file: %d rows, %d unique site_id (%s).",
                nrow(env_raw), n_sites_total,
                if (site_unique) "each appears once" else "DUPLICATES PRESENT"))

# Verify all source columns exist, then extract + coerce to numeric.
miss_cols <- setdiff(unname(FEATURE_MAP), colnames(env_raw))
if (length(miss_cols)) stop("Missing source columns: ",
                            paste(miss_cols, collapse = ", "))

env <- env_raw[, c("site_id", unname(FEATURE_MAP))]
names(env) <- c("site_id", PREDICTORS)
env <- env %>%
  distinct(site_id, .keep_all = TRUE) %>%
  mutate(across(all_of(PREDICTORS), ~ suppressWarnings(as.numeric(.))))

# Drop excluded site(s) so collinearity uses the same site set as all analyses.
present_excluded <- intersect(EXCLUDE_SITES, env$site_id)
env <- env %>% filter(!site_id %in% EXCLUDE_SITES)
if (length(present_excluded))
  message(sprintf("    Excluded site(s): %s (e.g. site 04 = Hachola/Hula, no phenotype data).",
                  paste(present_excluded, collapse = ", ")))

n_sites <- nrow(env)
message(sprintf("    Sites used for collinearity: %d", n_sites))

# -----------------------------------------------------------------------------
# 2. MISSINGNESS PER PREDICTOR
# -----------------------------------------------------------------------------

message("--- Missingness per predictor ---")

missing_df <- data.frame(
  predictor    = PREDICTORS,
  display_label= unname(PREDICTOR_LABELS[PREDICTORS]),
  total_sites  = n_sites,
  non_missing_sites = vapply(PREDICTORS, function(p) sum(!is.na(env[[p]])), integer(1)),
  stringsAsFactors = FALSE
)
missing_df$missing_sites <- missing_df$total_sites - missing_df$non_missing_sites
write.csv(missing_df[, c("predictor", "total_sites",
                         "non_missing_sites", "missing_sites")],
          file.path(DIR_TAB, "C2col_predictor_missingness.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. PREDICTOR DEFINITIONS TABLE
# -----------------------------------------------------------------------------

defs_df <- data.frame(
  internal_name     = PREDICTORS,
  source_csv_column = unname(FEATURE_MAP[PREDICTORS]),
  display_label     = unname(PREDICTOR_LABELS[PREDICTORS]),
  units             = unname(PREDICTOR_UNITS[PREDICTORS]),
  non_missing_sites = missing_df$non_missing_sites,
  stringsAsFactors  = FALSE
)
write.csv(defs_df, file.path(DIR_TAB, "C2col_predictor_definitions.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. PEARSON CORRELATION MATRIX (pairwise complete observations)
# -----------------------------------------------------------------------------

message("--- Pearson correlation matrix (pairwise complete obs) ---")

X <- as.matrix(env[, PREDICTORS])
cor_mat <- cor(X, use = "pairwise.complete.obs", method = "pearson")

write.csv(as.data.frame(round(cor_mat, 6)) %>% rownames_to_column("predictor"),
          file.path(DIR_TAB, "C2col_environment_correlation_matrix.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. PAIRWISE TABLE + n complete sites per pair + |r| > 0.8 flag
# -----------------------------------------------------------------------------

n_complete_pair <- function(a, b) sum(stats::complete.cases(env[[a]], env[[b]]))

pair_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
pairs_df <- data.frame(
  predictor_1 = rownames(cor_mat)[pair_idx[, 1]],
  predictor_2 = colnames(cor_mat)[pair_idx[, 2]],
  stringsAsFactors = FALSE
)
pairs_df$label_1       <- unname(PREDICTOR_LABELS[pairs_df$predictor_1])
pairs_df$label_2       <- unname(PREDICTOR_LABELS[pairs_df$predictor_2])
pairs_df$n_complete_sites <- mapply(n_complete_pair, pairs_df$predictor_1, pairs_df$predictor_2)
pairs_df$pearson_r     <- round(cor_mat[pair_idx], 6)
pairs_df$abs_r         <- round(abs(cor_mat[pair_idx]), 6)
pairs_df$abs_r_gt_0.8  <- pairs_df$abs_r > COLLIN_CUT
pairs_df$direction     <- ifelse(pairs_df$pearson_r >= 0, "positive", "negative")
pairs_df <- pairs_df[order(-pairs_df$abs_r), ]

write.csv(pairs_df[, c("predictor_1", "predictor_2", "label_1", "label_2",
                       "n_complete_sites", "pearson_r", "abs_r",
                       "abs_r_gt_0.8", "direction")],
          file.path(DIR_TAB, "C2col_environment_pairwise_correlations.csv"),
          row.names = FALSE)

high_df <- pairs_df[pairs_df$abs_r_gt_0.8, ]
write.csv(high_df[, c("predictor_1", "predictor_2", "label_1", "label_2",
                      "n_complete_sites", "pearson_r", "abs_r", "direction")],
          file.path(DIR_TAB, "C2col_pairs_abs_r_gt_0.8.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. HEATMAP (10 x 10), r in each cell, |r|>0.8 cells outlined in black
# -----------------------------------------------------------------------------

message("--- Building correlation heatmap ---")

lev   <- unname(PREDICTOR_LABELS[PREDICTORS])
long  <- expand.grid(p1 = PREDICTORS, p2 = PREDICTORS, stringsAsFactors = FALSE)
long$r        <- mapply(function(a, b) cor_mat[a, b], long$p1, long$p2)
long$lab1     <- factor(unname(PREDICTOR_LABELS[long$p1]), levels = lev)
long$lab2     <- factor(unname(PREDICTOR_LABELS[long$p2]), levels = rev(lev))
long$is_diag  <- long$p1 == long$p2
long$flag_hi  <- abs(long$r) > COLLIN_CUT & !long$is_diag   # mark strong off-diagonal pairs

short_axis <- function(x) sub("\\s*\\(.*\\)\\s*$", "", x)

p <- ggplot(long, aes(lab1, lab2, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  # bold black outline marks |r| > 0.8 pairs
  geom_tile(data = subset(long, flag_hi),
            fill = NA, colour = "black", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%.2f", r),
                fontface = ifelse(flag_hi, "bold", "plain")),
            size = 3.7, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.5, 0, 0.5, 1), name = "Pearson r") +
  scale_x_discrete(labels = short_axis, position = "top") +
  scale_y_discrete(labels = short_axis) +
  coord_equal() +
  labs(x = NULL, y = NULL,
       caption = "Black outline / bold value: |r| > 0.8") +
  theme_pub(base_size = 13) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 0, face = "bold"),
        axis.text.y  = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0.5, size = 11, face = "italic"),
        legend.position = "right")

save_plot(p, "C2col_environment_correlation_heatmap", w = 8.5, h = 8)

# -----------------------------------------------------------------------------
# 7. PLAIN-TEXT REPORT
# -----------------------------------------------------------------------------

rep_path <- file.path(DIR_REP, "C2col_collinearity_report.txt")
sink(rep_path)
cat("Environmental-predictor collinearity (STAGE 1)\n")
cat(strrep("=", 64), "\n", sep = "")
cat("Input file:          ", ENV_FILE, "\n", sep = "")
cat("Sites in file:       ", n_sites_total, "\n", sep = "")
cat("Excluded site(s):    ", if (length(present_excluded))
      paste(present_excluded, collapse = ", ") else "none",
    "  (site 04 = Hachola/Hula: no phenotype data, excluded from all analyses)\n", sep = "")
cat("Sites used:          ", n_sites, "\n", sep = "")
cat("Each site once:      ", if (site_unique) "YES" else "NO (duplicates!)", "\n", sep = "")
cat("Predictors:          ", length(PREDICTORS), "\n\n", sep = "")
cat("Missingness per predictor (predictor: non-missing / total):\n")
for (i in seq_len(nrow(missing_df)))
  cat(sprintf("  %-16s %d / %d (missing %d)\n",
              missing_df$predictor[i], missing_df$non_missing_sites[i],
              missing_df$total_sites[i], missing_df$missing_sites[i]))
cat("\nPairs with |r| > ", COLLIN_CUT, " (highest |r| first):\n", sep = "")
if (nrow(high_df)) {
  for (i in seq_len(nrow(high_df)))
    cat(sprintf("  %-16s ~ %-16s  r = %+.3f (%s)  | n = %d complete sites\n",
                high_df$predictor_1[i], high_df$predictor_2[i],
                high_df$pearson_r[i], high_df$direction[i],
                high_df$n_complete_sites[i]))
} else cat("  (none)\n")
cat("\nNOTE: No nutritional-trait/BLUP data were used. No FDR, no regression,\n")
cat("no predictor selection, and no automatic variable removal were performed.\n")
sink()

# -----------------------------------------------------------------------------
# 8. CONSOLE SUMMARY
# -----------------------------------------------------------------------------

message("\n", strrep("=", 60))
message("STAGE 1 (collinearity only) finished. Output root: ", OUT_ROOT)
message(strrep("=", 60))
message(sprintf("Sites: %d (unique: %s)", n_sites, site_unique))
message("Total missing values across all predictors: ", sum(missing_df$missing_sites))
message("Pairs with |r| > ", COLLIN_CUT, ":")
if (nrow(high_df)) {
  for (i in seq_len(nrow(high_df)))
    message(sprintf("  %s ~ %s : r=%+.3f (%s), n=%d",
                    high_df$predictor_1[i], high_df$predictor_2[i],
                    high_df$pearson_r[i], high_df$direction[i],
                    high_df$n_complete_sites[i]))
} else message("  none")
message("\nTables:  ", DIR_TAB)
message("Figures: ", DIR_FIG)
message("Report:  ", rep_path)
