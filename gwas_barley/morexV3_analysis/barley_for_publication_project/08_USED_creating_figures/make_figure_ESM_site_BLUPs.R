#!/usr/bin/env Rscript
# make_figure_ESM_site_BLUPs.R — supplementary figure (Online Resource, ESM number fixed later) of
# the TAG manuscript: centered trait BLUPs of the four nutritional traits across the 29 sampling
# sites, ordered and coloured by ecological region. This is OLD Fig. 2a.
#
# WHY (S. Hübner's comments on the figures, 2026-09-30, user decision): "2a can go to the supmat".
# Old Fig. 2 was dissolved: 2b + 2c merged and 2d moved into the new Fig. 1 (make_figure_1.R), and
# 2a moved here. TAG wants supplementary figures as PDF (TAG_requirements.md, ESM).
#
# NOTHING ABOUT THE ANALYSIS CHANGES. Drawn exactly as panel a of the former make_figure_2.R
# (2026-09-27; itself panel A4 of 00_THIN_.../01_phenotypic_analysis_no_GxE_v3_streamlined.R): same
# data, geoms, palette, region order, site order, axis labels and legend wording. Differences, all
# rendering only: its own page (174 x 120 mm instead of 174 x ~78 mm inside Fig. 2), PDF (vector,
# cairo_pdf, fonts embedded) plus a PNG preview, and a fixed seed for the jittered points
# (set.seed(1); make_figure_2.R had none, so the points moved slightly at every run).
#
# Input (read-only): 00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables/A4_site_boxplot_data.csv
# Outputs: Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf (the Online Resource file), ESM_site_BLUPs.png (preview)
# The ESM title page TAG asks for (article title, journal, authors, corresponding author) is added
# when the Online Resources are numbered and assembled, not here.
#
# Created 2026-09-30 (S. Hübner's structural comments). Run: Rscript make_figure_ESM_site_BLUPs.R

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
TAB  <- file.path(ROOT, "00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables")
OUT  <- file.path(ROOT, "08_USED_creating_figures", "Figure_ESM_site_BLUPs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FONT <- "Liberation Sans"; PT <- 9; DPI <- 600; W_MM <- 174; H_MM <- 120

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

TRAIT_ORDER  <- c("Protein", "Starch", "β-glucan", "Fiber")
REGION_ORDER <- c("North", "Coast", "Desert", "HZ1 (North-Coast)",
                  "HZ2 (North-Desert)", "HZ3 (Coast-Desert)")
PAL_REGION <- c("North" = "#2166AC", "Coast" = "#1B7837", "Desert" = "#B2182B",
                "HZ1 (North-Coast)" = "#67A9CF", "HZ2 (North-Desert)" = "#762A83",
                "HZ3 (Coast-Desert)" = "#E08214")

site <- read.csv(file.path(TAB, "A4_site_boxplot_data.csv"), fileEncoding = "UTF-8",
                 colClasses = c(site = "character")) %>%
  mutate(Trait  = factor(Trait, levels = TRAIT_ORDER),
         region = factor(region, levels = REGION_ORDER))
site$site <- factor(site$site, levels = site %>% distinct(region, site) %>%
                      arrange(region, site) %>% pull(site))
stopifnot(!anyNA(site$Trait), !anyNA(site$region), nlevels(site$site) == 29,
          length(unique(site$short_Tag)) == 290)

set.seed(1)
p <- ggplot(site, aes(x = site, y = value_centered, fill = region)) +
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

cairo_pdf(file.path(OUT, "ESM_site_BLUPs.pdf"), width = W_MM / 25.4, height = H_MM / 25.4,
          family = FONT)
set.seed(1); print(p); invisible(dev.off())
png(file.path(OUT, "ESM_site_BLUPs.png"), width = W_MM, height = H_MM, units = "mm", res = DPI,
    type = "cairo", family = FONT, bg = "white")
set.seed(1); print(p); invisible(dev.off())
cat(sprintf("[ESM site BLUPs] OK -> %s  (%d x %d mm; PDF + PNG preview)\n", OUT, W_MM, H_MM))
