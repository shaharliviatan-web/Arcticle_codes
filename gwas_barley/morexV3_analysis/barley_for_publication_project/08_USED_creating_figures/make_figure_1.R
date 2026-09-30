#!/usr/bin/env Rscript
# make_figure_1.R — Fig. 1 of the TAG manuscript (Results ch. 1): genetic architecture of the
# grain nutritional and agro-morphological traits.
#
#   a  variance partitioning of 8 traits (genetics / season + block / G x E / residual)
#   b  per-genotype reaction norms across the three seasons, 4 nutritional traits,
#      with four representative accessions highlighted per trait
#
# WHY THIS SCRIPT EXISTS. The content is exactly the mini paper's Figure 1, which was assembled by
# `assemble_figures.py 1`: it rescaled two panel PDFs from step 00 to a 160 mm-wide canvas at
# 300 dpi. That route cannot meet TAG: the width is not 84 or 174 mm, plots with lettering need
# >= 600 dpi, there is no TIFF, the stamped panel letters were upper case at ~15 pt (TAG allows
# 8-12 pt, lower case a/b), and every panel's lettering ended up at whatever size its own rescale
# factor produced. Here the panels are redrawn at final size instead, so the type size is exact.
#
# NOTHING ABOUT THE ANALYSIS CHANGES (user instruction 2026-09-27). The panels are rebuilt from the
# step-00 output tables, with the same data, geoms, palettes, factor orders, axis labels and legend
# wording as `00_THIN_Generate_Plots_For_Publication/02_phenotypic_analysis_GxE_v2_streamlined.R`
# (plots B1 and B2). Step 00 is neither re-run nor modified; only the rendering differs:
#   * one figure drawn at 174 mm, not two panels rescaled to 160 mm;
#   * a true 9 pt for all lettering (theme_pub's base_size 16 at ~0.45x rescale became ~7 pt);
#   * panel letters a / b, bold, 9 pt (were A / B at ~15 pt);
#   * 600 dpi, RGB, PNG (docx build) + TIFF/LZW (submission).
#
# TAG figure spec: 10_USED_Paper_writing/TAG_requirements.md.
#
# Inputs (read-only, 00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables/):
#   Variance_components_GxE.csv              a: one row per trait x variance component
#   gxe_long_data.csv                        b: genotype x season means, 4 nutritional traits
#   B2_reaction_norm_highlights_selected.csv b: the 4 highlighted accessions per trait
# Outputs: Figure_1/Fig1.png (docx build), Figure_1/Fig1.tif (submission, LZW)
#
# Created 2026-09-27. Run: Rscript make_figure_1.R   (seconds)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
TAB  <- file.path(ROOT, "00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables")
OUT  <- file.path(ROOT, "08_USED_creating_figures", "Figure_1")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── TAG style ─────────────────────────────────────────────────────────────────
FONT <- "Liberation Sans"        # metric-compatible with Arial
PT   <- 9                        # TRUE size on the page (TAG 8-12)
DPI  <- 600
W_MM <- 174                      # TAG full width
H_MM <- 74

# theme_pub() of step 00, at base_size = PT instead of 16. Element-for-element identical.
theme_tag <- function(base_size = PT) {
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

# Factor orders and palettes, copied from step-00 script 02
TRAIT_ORDER_8 <- c("Protein", "Starch", "β-glucan", "Fiber",
                   "Flowering time", "Tillers", "Grain weight", "Spike length")
TRAIT_ORDER_4 <- TRAIT_ORDER_8[1:4]
PAL_VARCOMP   <- c("Genetics" = "#2CA02C", "Season & Block" = "#1F77B4",
                   "G x E" = "#FF7F0E", "Residual" = "#BDBDBD")
PAL_HIGHLIGHT <- c("Most stable" = "#009E73", "Strongest increase" = "#D55E00",
                   "Strongest decrease" = "#0072B2", "Strongest crossover" = "#CC79A7")

# ── a: variance partitioning ──────────────────────────────────────────────────
var_perc <- read.csv(file.path(TAB, "Variance_components_GxE.csv"), fileEncoding = "UTF-8") %>%
  mutate(TraitLabel = factor(TraitLabel, levels = TRAIT_ORDER_8),
         Component  = factor(Component,
                             levels = c("Residual", "G x E", "Season & Block", "Genetics")))
stopifnot(nlevels(droplevels(var_perc$TraitLabel)) == 8, !anyNA(var_perc$TraitLabel),
          all(abs(tapply(var_perc$Percentage, var_perc$TraitLabel, sum) - 100) < 0.2))

p_a <- ggplot(var_perc, aes(x = TraitLabel, y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = PAL_VARCOMP, name = NULL,
                    breaks = c("Genetics", "Season & Block", "G x E", "Residual")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), breaks = seq(0, 100, 25)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "Variance explained (%)") +
  theme_tag() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))

# ── b: reaction norms ─────────────────────────────────────────────────────────
gxe <- read.csv(file.path(TAB, "gxe_long_data.csv"), fileEncoding = "UTF-8")
hl  <- read.csv(file.path(TAB, "B2_reaction_norm_highlights_selected.csv"), fileEncoding = "UTF-8")
stopifnot(nrow(hl) == 16, length(unique(gxe$short_Tag)) == 290)

pd <- gxe %>%
  left_join(hl[, c("Trait", "Pattern", "short_Tag")], by = c("Trait", "short_Tag")) %>%
  mutate(Group   = ifelse(is.na(Pattern), "Background", "Highlight"),
         Pattern = factor(Pattern, levels = names(PAL_HIGHLIGHT)),
         season  = factor(season, levels = c("2019", "2020", "2021")),
         Trait   = factor(Trait, levels = TRAIT_ORDER_4))
stopifnot(!anyNA(pd$Trait), sum(pd$Group == "Highlight") == nrow(hl) * 3)

p_b <- ggplot() +
  geom_line(data = filter(pd, Group == "Background"),
            aes(x = season, y = Mean_Value, group = short_Tag),
            colour = "grey80", linewidth = 0.3, alpha = 0.45) +
  geom_line(data = filter(pd, Group == "Highlight"),
            aes(x = season, y = Mean_Value, group = short_Tag, colour = Pattern),
            linewidth = 0.7) +
  geom_point(data = filter(pd, Group == "Highlight"),
             aes(x = season, y = Mean_Value, colour = Pattern), size = 1.0) +
  scale_colour_manual(values = PAL_HIGHLIGHT, name = "Representative pattern",
                      labels = c("Most stable" = "stable", "Strongest increase" = "increase",
                                 "Strongest decrease" = "decrease",
                                 "Strongest crossover" = "crossover")) +
  facet_wrap(~ Trait, scales = "free_y", ncol = 2) +
  labs(x = "Growing season", y = "Centered phenotypic value") +
  theme_tag() +
  guides(colour = guide_legend(nrow = 1, title.position = "top", title.hjust = 0.5)) +
  theme(legend.title = element_text(hjust = 0.5))

# ── assemble ──────────────────────────────────────────────────────────────────
fig <- (p_a | p_b) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = FONT, face = "bold", size = PT, hjust = 0, vjust = 1))

png(file.path(OUT, "Fig1.png"), width = W_MM, height = H_MM, units = "mm", res = DPI,
    type = "cairo", family = FONT, bg = "white")
print(fig); invisible(dev.off())
stopifnot(system2("python3", c("-c", shQuote(sprintf(
  "from PIL import Image; Image.open('%s').convert('RGB').save('%s', compression='tiff_lzw', dpi=(%d,%d))",
  file.path(OUT, "Fig1.png"), file.path(OUT, "Fig1.tif"), DPI, DPI)))) == 0)
cat(sprintf("[fig1] OK -> %s  (%d x %d mm, %d dpi, %d pt)\n", OUT, W_MM, H_MM, DPI, PT))
