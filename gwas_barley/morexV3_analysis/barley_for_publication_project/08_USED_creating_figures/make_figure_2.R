#!/usr/bin/env Rscript
# make_figure_2.R — Fig. 2 of the TAG manuscript (Results ch. 1): ecological architecture of the
# grain nutritional traits.
#
#   a  centered trait BLUPs across the 29 sampling sites, ordered and coloured by ecological region
#   b  Pearson correlations among the four nutritional traits
#   c  Pearson correlations between nutritional and agro-morphological traits
#   d  Pearson correlations between site-mean nutritional BLUPs and eight environmental variables
#
# Layout (changed 2026-09-27 on the user's request: panel a was too dense to read at half width):
# a spans the full width on top, b | c share the middle row, d spans the full width at the bottom.
# The mini paper packed the four panels into a 2 x 2 block with a and d on top, which left a with
# 87 mm for 29 sites x 4 traits. Each letter still carries exactly the content the Ch. 1 caption
# gives it, and the letters now also run in reading order.
#
# WHY THIS SCRIPT EXISTS, and what is NOT changing: see the header of make_figure_1.R. Same rule
# here (user instruction 2026-09-27) — the panels are redrawn from the step-00 output tables with
# the same data, geoms, palettes, factor orders, labels and significance symbols as
# `00_THIN_.../01_phenotypic_analysis_no_GxE_v3_streamlined.R` (A4, A12 panels A and B) and
# `00_THIN_.../03b_trait_environment_correlations.R` (C2 all-32 heatmap). Step 00 is neither
# re-run nor modified. Only the rendering differs: 174 mm instead of 160 mm, a true 9 pt for all
# lettering instead of per-panel rescales of a base_size-16 theme, panel letters a-d instead of
# A-D, 600 dpi, and TIFF alongside PNG.
#
# The step-00 panels carry hand-enlarged fonts (axis text 24-30 pt, in-tile labels 4.8-5.6 mm)
# because each was drawn large and then shrunk by `assemble_figures.py`. Drawing at final size
# makes that unnecessary: every size here is the size on the page.
#
# Inputs (read-only, 00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/):
#   tables/A4_site_boxplot_data.csv                            a
#   tables/A12_nutri_pairs_uncorrectedP.csv                    b
#   tables/A12_nutri_morpho_pairs_localFDR.csv                 c (raw-p stars, column `sig`)
#   C2_trait_environment_correlations/tables/C2corr_trait_environment_all32.csv   d
# Outputs: Figure_2/Fig2.png (docx build), Figure_2/Fig2.tif (submission, LZW)
#
# Created 2026-09-27. Run: Rscript make_figure_2.R   (seconds)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
S1   <- file.path(ROOT, "00_THIN_Generate_Plots_For_Publication/outputs/subsection_1")
TAB  <- file.path(S1, "tables")
C2   <- file.path(S1, "C2_trait_environment_correlations/tables")
OUT  <- file.path(ROOT, "08_USED_creating_figures", "Figure_2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── TAG style ─────────────────────────────────────────────────────────────────
FONT <- "Liberation Sans"
PT   <- 9                        # TRUE size on the page (TAG 8-12)
DPI  <- 600
W_MM <- 174
H_MM <- 230
mm_text <- function(pt) pt / .pt  # geom_text size is in mm

theme_tag <- function(base_size = PT) {          # theme_pub() of step 00 at base_size = PT
  theme_bw(base_size = base_size, base_family = FONT) +
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
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(3.2, "mm"),
      legend.margin    = margin(t = 0, b = 0),
      plot.margin      = margin(2, 3, 1, 2)
    )
}

TRAIT_ORDER <- c("Protein", "Starch", "β-glucan", "Fiber")
REGION_ORDER <- c("North", "Coast", "Desert", "HZ1 (North-Coast)",
                  "HZ2 (North-Desert)", "HZ3 (Coast-Desert)")
PAL_REGION <- c("North" = "#2166AC", "Coast" = "#1B7837", "Desert" = "#B2182B",
                "HZ1 (North-Coast)" = "#67A9CF", "HZ2 (North-Desert)" = "#762A83",
                "HZ3 (Coast-Desert)" = "#E08214")
FILL_R <- function(...) scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                             midpoint = 0, limits = c(-1, 1),
                                             breaks = c(-1, -0.5, 0, 0.5, 1),
                                             name = "Pearson r", ...)

# ── a: site boxplots ──────────────────────────────────────────────────────────
site <- read.csv(file.path(TAB, "A4_site_boxplot_data.csv"), fileEncoding = "UTF-8",
                 colClasses = c(site = "character")) %>%
  mutate(Trait  = factor(Trait, levels = TRAIT_ORDER),
         region = factor(region, levels = REGION_ORDER))
site$site <- factor(site$site, levels = site %>% distinct(region, site) %>%
                      arrange(region, site) %>% pull(site))
stopifnot(!anyNA(site$Trait), !anyNA(site$region), nlevels(site$site) == 29)

p_a <- ggplot(site, aes(x = site, y = value_centered, fill = region)) +
  geom_boxplot(colour = "grey30", outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
  geom_jitter(width = 0.18, height = 0, size = 0.25, alpha = 0.55, colour = "grey20") +
  facet_wrap(~ Trait, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = PAL_REGION, name = "Ecological region", drop = FALSE) +
  scale_x_discrete(labels = function(x)                       # label every 5th site
    ifelse(!is.na(suppressWarnings(as.integer(x))) & suppressWarnings(as.integer(x)) %% 5 == 0,
           suppressWarnings(as.integer(x)), "")) +
  labs(x = "Sampling site", y = "BLUP (centered)") +
  theme_tag() +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, title.position = "top", title.hjust = 0.5)) +
  theme(legend.title = element_text(hjust = 0.5))

# ── b: nutritional trait correlations ─────────────────────────────────────────
x_levels <- c("Starch", "β-glucan", "Fiber")     # no Protein column
y_levels <- c("Protein", "Starch", "β-glucan")   # no Fiber row
heat_nutri <- read.csv(file.path(TAB, "A12_nutri_pairs_uncorrectedP.csv"), fileEncoding = "UTF-8") %>%
  mutate(Trait1 = factor(Trait1, levels = y_levels),
         Trait2 = factor(Trait2, levels = x_levels),
         label  = sprintf("%.2f%s", r, ifelse(is.na(sig_raw), "", sig_raw)))
stopifnot(nrow(heat_nutri) == 6)

p_b <- ggplot(heat_nutri, aes(x = Trait2, y = Trait1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label, colour = ifelse(abs(r) > 0.45, "white", "black")),
            size = mm_text(7), fontface = "bold") +
  FILL_R() + scale_colour_identity() +
  scale_y_discrete(limits = rev(y_levels), drop = FALSE) +
  scale_x_discrete(limits = x_levels, drop = FALSE, position = "top") +
  labs(x = NULL, y = NULL) +
  coord_equal() +
  guides(fill = guide_colourbar(barwidth = unit(0.25, "cm"), barheight = unit(2.0, "cm"))) +
  theme_tag() +
  theme(panel.grid = element_blank(), legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# ── c: nutritional x agro-morphological correlations ──────────────────────────
NUTRI <- TRAIT_ORDER
d_ord <- read.csv(file.path(TAB, "A12_nutri_morpho_pairs_localFDR.csv"), fileEncoding = "UTF-8") %>%
  mutate(nutri = factor(nutri, levels = NUTRI),
         sig   = ifelse(is.na(sig), "", sig),
         morpho_lab = paste(morpho, as.integer(nutri), sep = "___"))
stopifnot(nrow(d_ord) == 16, !anyNA(d_ord$nutri))
d_ord$morpho_lab <- factor(d_ord$morpho_lab,
                           levels = d_ord %>% arrange(nutri, r) %>% pull(morpho_lab))

p_c <- ggplot(d_ord, aes(x = morpho_lab, y = r)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_point(size = 1.4, colour = "grey45") +
  geom_text(aes(label = sig, y = r), vjust = -0.7, size = mm_text(8),
            fontface = "bold", colour = "grey25") +
  facet_wrap(~ nutri, scales = "free", ncol = 2) +
  scale_x_discrete(labels = function(x) sub("___.*", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0.14, 0.50)), n.breaks = 3) +   # fewer ticks + headroom for the stars at 9 pt
  labs(x = NULL, y = "Pearson's r") +
  theme_tag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())

# ── d: trait x environment correlations ───────────────────────────────────────
PRED_LEVELS <- c("pH", "Electrical conductivity (mS/cm)", "Clay (%)", "Sand (%)",
                 "March temperature (°C)", "Precipitation (mm)",
                 "Organic carbon (mg/kg)", "Total N (mg/kg)")
heat_env <- read.csv(file.path(C2, "C2corr_trait_environment_all32.csv"), fileEncoding = "UTF-8") %>%
  mutate(TraitLabel = factor(Trait, levels = rev(TRAIT_ORDER)),
         PredLabel  = factor(PredLabel, levels = PRED_LEVELS),
         sig_symbol = ifelse(is.na(sig_symbol), "", sig_symbol))
stopifnot(nrow(heat_env) == 32, !anyNA(heat_env$TraitLabel), !anyNA(heat_env$PredLabel))

p_d <- ggplot(heat_env, aes(x = PredLabel, y = TraitLabel, fill = pearson_r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = subset(heat_env, pearson_p < 0.05), fill = NA, colour = "black",
            linewidth = 0.7) +
  geom_text(aes(label = ifelse(sig_symbol == "", sprintf("%.2f", pearson_r),
                               sprintf("%.2f\n%s", pearson_r, sig_symbol))),
            size = mm_text(6.5), lineheight = 0.78, colour = "black") +
  FILL_R() +
  scale_x_discrete(labels = function(x) sub("\\s*\\(.*\\)\\s*$", "", x)) +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(barwidth = unit(38, "mm"), barheight = unit(2.6, "mm"),
                                title.position = "left", title.vjust = 1)) +
  theme_tag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(), legend.position = "bottom")

# ── assemble: a | d over b | c ────────────────────────────────────────────────
fig <- p_a / (p_b + p_c + plot_layout(widths = c(0.8, 1.2))) / p_d +
  plot_layout(heights = c(1.12, 1.40, 0.80)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = FONT, face = "bold", size = PT, hjust = 0, vjust = 1))

png(file.path(OUT, "Fig2.png"), width = W_MM, height = H_MM, units = "mm", res = DPI,
    type = "cairo", family = FONT, bg = "white")
print(fig); invisible(dev.off())
stopifnot(system2("python3", c("-c", shQuote(sprintf(
  "from PIL import Image; Image.open('%s').convert('RGB').save('%s', compression='tiff_lzw', dpi=(%d,%d))",
  file.path(OUT, "Fig2.png"), file.path(OUT, "Fig2.tif"), DPI, DPI)))) == 0)
cat(sprintf("[fig2] OK -> %s  (%d x %d mm, %d dpi, %d pt)\n", OUT, W_MM, H_MM, DPI, PT))
