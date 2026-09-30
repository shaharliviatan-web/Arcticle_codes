#!/usr/bin/env Rscript
# make_figure_1.R — Fig. 1 of the TAG manuscript (Results ch. 1): genetic and ecological architecture
# of the grain nutritional traits.
#
#   a  variance partitioning of 8 traits (genetics / season + block / G x E / residual)   = old Fig. 1a
#   b  per-genotype reaction norms across the three seasons, 4 nutritional traits,
#      four highlighted accessions per trait                                              = old Fig. 1b
#   c  ONE correlation heatmap, the 4 nutritional traits as rows: their correlations with one
#      another as a lower triangle (columns Protein, Starch, beta-glucan), then with the
#      4 morphological traits (columns Flowering time, Plant height, Grain weight, Grain number) = old Fig. 2b + 2c
#   d  site-mean nutritional BLUPs x eight environmental variables                        = old Fig. 2d
#
# HISTORY. Rewritten 2026-09-30 for S. Hübner's comments on the figures (user decisions; hand-over
# 10_USED_Paper_writing/new_publishing_paper/build/STRUCTURE_CHANGES_2026-10.md): old Fig. 2 was
# dissolved — 2a went to the Online Resources (make_figure_ESM_site_BLUPs.R), 2b and 2c were merged
# into one heatmap in the style of 2b/2d, 2d joined this figure — and old Figs. 3-5 became Figs. 2-4.
# The previous make_figure_1.R (a and b only, 174 x 74 mm, 2026-09-27) and make_figure_2.R (old
# Fig. 2, 2026-09-27/28) were replaced by this script (both remain in git). Layout and caption
# approved by the user 2026-09-30 (build/proposed/S0_QUESTIONS.md, Q1 "triangle", made a lower
# triangle; Q2 caption A).
#
# NOTHING ABOUT THE ANALYSIS CHANGES. Step 00 is neither re-run nor modified. Panels a and b are
# drawn exactly as in the previous make_figure_1.R; panel d exactly as panel d of the previous
# make_figure_2.R, with the rendering changes below; panel c is new only in form.
#   c  values: A12_nutri_pairs_uncorrectedP.csv (r, stars `sig_raw`; old 2b) and
#      A12_nutri_morpho_pairs_localFDR.csv (r, stars `sig`, the raw-P column; old 2c, never
#      `sig_local`). SAME SIGNIFICANCE BASIS in both blocks, checked here (script stops otherwise):
#      both star columns are star_fn() of step-00 script 01 applied to the raw two-sided P of
#      cor.test() on the same 290-accession BLUP matrix (n = 290 for every pair), thresholds
#      0.05 / 0.01 / 0.001, no multiple-testing correction. Panel d uses the same thresholds on
#      raw P (C2corr_trait_environment_all32.csv, 29 sites). The caption says so.
#      Each nutritional pair is drawn once, below the diagonal (row = the later trait of
#      Protein, Starch, beta-glucan, Fiber); the Fiber column would be empty and is not drawn, so c
#      has 7 columns (d has 8; their columns are not aligned). A dark line separates the two blocks.
#   c, d  tile text 8 pt (old 2b/2d had 6.5-7 pt, below TAG's 8 pt minimum), black on every tile
#      (contrast >= 4.5:1 at |r| <= 0.78, the largest |r| drawn; white text fails there); one
#      shared colour bar for c and d.
#
# TAG figure spec: 10_USED_Paper_writing/TAG_requirements.md (174 mm, <= 234 mm high, Arial 8-12 pt,
# 600 dpi, RGB, panels a-d, no titles inside the image).
#
# Inputs (read-only, 00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/):
#   tables/Variance_components_GxE.csv                          a
#   tables/gxe_long_data.csv, tables/B2_reaction_norm_highlights_selected.csv   b
#   tables/A12_nutri_pairs_uncorrectedP.csv, tables/A12_nutri_morpho_pairs_localFDR.csv   c
#   C2_trait_environment_correlations/tables/C2corr_trait_environment_all32.csv   d
# Outputs: Figure_1/Fig1.png (docx build), Figure_1/Fig1.tif (submission, LZW, converted from the PNG)
#
# Run: Rscript make_figure_1.R   (seconds)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
S1   <- file.path(ROOT, "00_THIN_Generate_Plots_For_Publication/outputs/subsection_1")
TAB  <- file.path(S1, "tables")
C2   <- file.path(S1, "C2_trait_environment_correlations/tables")
OUT  <- file.path(ROOT, "08_USED_creating_figures", "Figure_1")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── TAG style ─────────────────────────────────────────────────────────────────
FONT <- "Liberation Sans"        # metric-compatible with Arial
PT   <- 9                        # TRUE size on the page (TAG 8-12)
TILE_PT <- 8                     # in-tile numbers (TAG minimum)
DPI  <- 600
W_MM <- 174                      # TAG full width
H_MM <- 205                      # TAG limit 234
mm_text <- function(pt) pt / .pt # geom_text size is in mm

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

TRAIT_ORDER_8 <- c("Protein", "Starch", "β-glucan", "Fiber",
                   "Flowering time", "Tillers", "Grain weight", "Spike length")
TRAIT_ORDER_4 <- TRAIT_ORDER_8[1:4]
MORPHO_4      <- c("Flowering time", "Plant height", "Grain weight", "Grain number")
PAL_VARCOMP   <- c("Genetics" = "#2CA02C", "Season & Block" = "#1F77B4",
                   "G x E" = "#FF7F0E", "Residual" = "#BDBDBD")
PAL_HIGHLIGHT <- c("Most stable" = "#009E73", "Strongest increase" = "#D55E00",
                   "Strongest decrease" = "#0072B2", "Strongest crossover" = "#CC79A7")
FILL_R <- function(...) scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                             midpoint = 0, limits = c(-1, 1),
                                             breaks = c(-1, -0.5, 0, 0.5, 1),
                                             name = "Pearson r", ...)
star_fn <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",   # step-00 script 01, star_fn()
                              ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# ── a: variance partitioning (as make_figure_1.R) ─────────────────────────────
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

# ── b: reaction norms (as make_figure_1.R) ────────────────────────────────────
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

# ── c: merged trait-correlation heatmap (old 2b + 2c) ─────────────────────────
nn <- read.csv(file.path(TAB, "A12_nutri_pairs_uncorrectedP.csv"), fileEncoding = "UTF-8")
nm <- read.csv(file.path(TAB, "A12_nutri_morpho_pairs_localFDR.csv"), fileEncoding = "UTF-8")
nn$sig_raw[is.na(nn$sig_raw)] <- ""; nm$sig[is.na(nm$sig)] <- ""
# same significance basis in both blocks (see header): raw two-sided P, n = 290, same thresholds
stopifnot(nrow(nn) == 6, nrow(nm) == 16, all(nn$n == 290), all(nm$n == 290),
          identical(nn$sig_raw, star_fn(nn$p)), identical(nm$sig, star_fn(nm$p)),
          setequal(nm$morpho, MORPHO_4), setequal(c(nn$Trait1, nn$Trait2), TRAIT_ORDER_4))

# lower triangle: each pair once, row = the later trait in TRAIT_ORDER_4, column = the earlier one
ord4 <- function(x) match(x, TRAIT_ORDER_4)
nutri_cells <- data.frame(
  row = ifelse(ord4(nn$Trait1) > ord4(nn$Trait2), nn$Trait1, nn$Trait2),
  col = ifelse(ord4(nn$Trait1) > ord4(nn$Trait2), nn$Trait2, nn$Trait1),
  r = nn$r, sig = nn$sig_raw)
morpho_cells <- data.frame(row = nm$nutri, col = nm$morpho, r = nm$r, sig = nm$sig)
COL_C <- c(TRAIT_ORDER_4[1:3], MORPHO_4)              # the Fiber column would be empty
heat_c <- bind_rows(nutri_cells, morpho_cells) %>%
  mutate(row = factor(row, levels = rev(TRAIT_ORDER_4)), col = factor(col, levels = COL_C))
stopifnot(!anyNA(heat_c$row), !anyNA(heat_c$col), nrow(heat_c) == 6 + 16,
          !anyDuplicated(heat_c[, c("row", "col")]),
          all(ord4(nutri_cells$row) > ord4(nutri_cells$col)))
N_NUTRI_COLS <- 3

tile_label <- function(r, s) ifelse(s == "", sprintf("%.2f", r), sprintf("%.2f\n%s", r, s))

p_c <- ggplot(heat_c, aes(x = col, y = row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = tile_label(r, sig)), size = mm_text(TILE_PT), lineheight = 0.78,
            colour = "black") +
  geom_vline(xintercept = N_NUTRI_COLS + 0.5, colour = "grey25", linewidth = 0.6) +
  FILL_R() +
  scale_x_discrete(limits = COL_C, drop = FALSE) +
  scale_y_discrete(limits = rev(TRAIT_ORDER_4), drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_tag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(), panel.grid.major = element_blank())

# ── d: trait x environment correlations (as make_figure_2.R d) ────────────────
PRED_LEVELS <- c("pH", "Electrical conductivity (mS/cm)", "Clay (%)", "Sand (%)",
                 "March temperature (°C)", "Precipitation (mm)",
                 "Organic carbon (mg/kg)", "Total N (mg/kg)")
heat_env <- read.csv(file.path(C2, "C2corr_trait_environment_all32.csv"), fileEncoding = "UTF-8") %>%
  mutate(TraitLabel = factor(Trait, levels = rev(TRAIT_ORDER_4)),
         PredLabel  = factor(PredLabel, levels = PRED_LEVELS),
         sig_symbol = ifelse(is.na(sig_symbol), "", sig_symbol))
stopifnot(nrow(heat_env) == 32, !anyNA(heat_env$TraitLabel), !anyNA(heat_env$PredLabel),
          identical(heat_env$sig_symbol, star_fn(heat_env$pearson_p)))

p_d <- ggplot(heat_env, aes(x = PredLabel, y = TraitLabel, fill = pearson_r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = subset(heat_env, pearson_p < 0.05), fill = NA, colour = "black",
            linewidth = 0.7) +
  geom_text(aes(label = tile_label(pearson_r, sig_symbol)),
            size = mm_text(TILE_PT), lineheight = 0.78, colour = "black") +
  FILL_R() +
  scale_x_discrete(labels = function(x) sub("\\s*\\(.*\\)\\s*$", "", x)) +
  labs(x = NULL, y = NULL) +
  theme_tag() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(), panel.grid.major = element_blank())

# ── assemble: a | b over c over d; one colour bar for c and d ─────────────────
CBAR <- guides(fill = guide_colourbar(barwidth = unit(38, "mm"), barheight = unit(2.6, "mm"),
                                      title.position = "left", title.vjust = 1))
cd  <- (p_c + CBAR) / (p_d + CBAR) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig <- (p_a | p_b) / cd +
  plot_layout(heights = c(74, H_MM - 74)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = FONT, face = "bold", size = PT, hjust = 0, vjust = 1))

png(file.path(OUT, "Fig1.png"), width = W_MM, height = H_MM, units = "mm", res = DPI,
    type = "cairo", family = FONT, bg = "white")
print(fig); invisible(dev.off())
stopifnot(system2("python3", c("-c", shQuote(sprintf(
  "from PIL import Image; Image.open('%s').convert('RGB').save('%s', compression='tiff_lzw', dpi=(%d,%d))",
  file.path(OUT, "Fig1.png"), file.path(OUT, "Fig1.tif"), DPI, DPI)))) == 0)
cat(sprintf("[fig1] OK -> %s  (%d x %d mm, %d dpi, %d pt)\n", OUT, W_MM, H_MM, DPI, PT))
