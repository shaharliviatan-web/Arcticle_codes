# 09_figure_desert_enrichment.R — Figure: the Desert and six-region enrichment tests.
#
# Same rows and labels as the carrier-origin figure (script 08), so the two read side by side.
#   a  per lead SNP: % of minor-allele carriers (filled) and of major-allele carriers (open)
#      from the Desert
#   b  per haplotype group: % of the group's members (filled) and of the other assigned
#      accessions of the same gene (open) from the Desert
# Hairline: the Desert share of the panel (50 of 290). Right: BH q of the two site-permutation
# tests (T02, T04): Desert (two-sided) and six regions; for the haplotype groups the six-region
# q is the gene-level test, printed on the gene's header row. q <= 0.05 in bold.
#
# Style as script 08 (TAG: 174 mm, 9 pt, 600 dpi, PNG + TIFF). The two point types differ in
# fill, not colour.
#
# Inputs:  results/tables/T02_lead_enrichment_tests.tsv, T04_haplotype_group_enrichment_tests.tsv
# Outputs: results/figures/Fig_desert_enrichment.{png,tif}

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_figure_helpers.R"))

t02 <- read_tsv(file.path(TAB, "T02_lead_enrichment_tests.tsv"), show_col_types = FALSE)
t04 <- read_tsv(file.path(TAB, "T04_haplotype_group_enrichment_tests.tsv"), show_col_types = FALSE)
rows <- build_rows(t02, t04)
PANEL_PCT <- 100 * 50 / 290
X_Q1 <- 118; X_Q2 <- 140; X_MAX <- 150

# 3 significant digits, so a q just above 0.05 never prints as "0.050" (bold marks q <= 0.05)
fmt_q <- function(q) ifelse(is.na(q), "", ifelse(q < 0.001, gsub("e-0?", "e\u2212", formatC(q, format = "e", digits = 1)),
                                                  formatC(q, digits = 3, format = "fg", flag = "#")))
pts <- bind_rows(
  t02 %>% transmute(key = lead_SNP, focal = desert_pct_minor, other = desert_pct_major,
                    q1 = q_desert, q2 = q_regions),
  t04 %>% transmute(key = paste(gene, group, sep = "|"), focal = desert_pct_group,
                    other = desert_pct_other_groups, q1 = q_desert, q2 = NA_real_),
  t04 %>% distinct(gene, q_regions_gene) %>%
    transmute(key = paste0("hdr_", gene), focal = NA_real_, other = NA_real_, q1 = NA_real_, q2 = q_regions_gene)) %>%
  left_join(select(rows, key, panel, y), by = "key")
stopifnot(!anyNA(pts$y))

dot_panel <- function(p) {
  r <- rows %>% filter(panel == p); d <- pts %>% filter(panel == p)
  lab <- r %>% filter(type != "header"); dd <- d %>% filter(!is.na(focal))
  q <- bind_rows(d %>% transmute(y, x = X_Q1, q = q1), d %>% transmute(y, x = X_Q2, q = q2)) %>%
    filter(!is.na(q)) %>% mutate(lab = fmt_q(q), face = ifelse(q <= 0.05, "bold", "plain"))
  top <- max(r$y) + 1
  ggplot() +
    geom_vline(xintercept = c(0, 25, 50, 75, 100), colour = GRID, linewidth = 0.25) +
    geom_vline(xintercept = PANEL_PCT, colour = INK2, linewidth = 0.3) +
    geom_segment(data = dd, aes(x = other, xend = focal, y = y, yend = y), colour = GREY, linewidth = 0.5) +
    geom_point(data = dd, aes(x = other, y = y, shape = "other"), size = 1.3, stroke = 0.45, colour = INK, fill = "white") +
    geom_point(data = dd, aes(x = focal, y = y, shape = "focal"), size = 1.4, colour = INK, fill = INK) +
    hdr_layer(r, 0) +
    geom_text(data = q, aes(x = x, y = y, label = lab, fontface = face), family = FONT, size = PT / .pt,
              colour = INK, hjust = 1) +
    annotate("text", x = c(X_Q1, X_Q2), y = top, label = c("q Desert", "q regions"), family = FONT,
             size = PT / .pt, colour = INK, hjust = 1, fontface = "bold") +
    scale_shape_manual(values = c(focal = 21, other = 21),
                       labels = c(focal = "Minor-allele carriers (a) / group members (b)",
                                  other = "Major-allele carriers (a) / other groups (b)"),
                       guide = guide_legend(nrow = 2, override.aes = list(fill = c(INK, "white"), size = 1.6))) +
    scale_y_continuous(breaks = lab$y, labels = lab$label, limits = c(min(r$y) - 0.6, top + 0.6), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, X_MAX), breaks = c(0, 25, 50, 75, 100), expand = c(0.01, 0)) +
    labs(x = NULL, y = NULL, tag = p) + theme_fig()
}
pa <- dot_panel("a") + theme(axis.text.x = element_blank(), legend.position = "none")
pb <- dot_panel("b") + labs(x = "Accessions from the Desert (%)") +
  theme(axis.title.x = element_text(hjust = 100 / X_MAX / 2 * 0.95), legend.position = "bottom",
        legend.justification = "center", legend.margin = margin(0, 0, 0, 0))

na <- sum(rows$panel == "a") + 1; nb <- sum(rows$panel == "b") + 1
ROW_MM <- 3.0
fig <- pa / pb + plot_layout(heights = unit(c(na * ROW_MM, nb * ROW_MM), "mm"))
render_fig(fig, "Fig_desert_enrichment", h_mm = (na + nb) * ROW_MM + 27)
