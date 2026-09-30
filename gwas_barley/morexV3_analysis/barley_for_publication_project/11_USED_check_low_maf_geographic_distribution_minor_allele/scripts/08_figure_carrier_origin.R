# 08_figure_carrier_origin.R — Figure: where the accessions come from.
#
#   a  one row per lead SNP (36), grouped by trait: minor-allele carrier / major-allele
#      carrier / no call, for each of the 290 accessions
#   b  one row per haplotype group of the four presented genes: member / member of another
#      group of the same gene / unassigned
# Columns: accessions ordered by region (Ch. 1 order), then site, then accession; region band
# above (named, Ch. 1 colours), site names below. Replaces the 2026-09-25 three-panel figure
# (its panels b and c were dropped for readability, user 2026-09-30; the Desert shares are in
# the enrichment figure, script 09).
#
# TAG spec, as 08_USED_creating_figures: 174 mm wide, <= 234 mm, Liberation Sans (Arial
# metric) 9 pt (8 pt until 2026-09-30), 600 dpi, RGB, PNG + LZW TIFF. Cell colors: CELL_* in 00_config.R. Region is never shown by colour alone: the band
# names every region (the Ch. 1 palette fails a colour-vision check on its red/green pair).
#
# Inputs:  intermediates/lead_genotype_long.tsv, haplotype_group_long.tsv, T02, T04
# Outputs: results/figures/Fig_carrier_origin.{png,tif}

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_config.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_figure_helpers.R"))

g   <- read_tsv(file.path(INTER, "lead_genotype_long.tsv"), show_col_types = FALSE)
h   <- read_tsv(file.path(INTER, "haplotype_group_long.tsv"), show_col_types = FALSE)
t02 <- read_tsv(file.path(TAB, "T02_lead_enrichment_tests.tsv"), show_col_types = FALSE)
t04 <- read_tsv(file.path(TAB, "T04_haplotype_group_enrichment_tests.tsv"), show_col_types = FALSE)

rows <- build_rows(t02, t04)
AX   <- build_axis(load_accessions())
STATE <- c("dark", "grey", "pale")
LEG <- c(dark = "Minor-allele carrier (a) / group member (b)",
         grey = "Major-allele carrier (a) / other group (b)",
         pale = "No call (a) / unassigned (b)")

# ---- cells ---------------------------------------------------------------------------------------
cell_a <- g %>% transmute(IID, key = lead_SNP,
                          state = recode(as.character(allele), minor = "dark", major = "grey", missing = "pale"))
cell_b <- bind_rows(lapply(split(h, h$gene), function(d) {
  bind_rows(lapply(setdiff(sort(unique(d$hap)), "0"), function(k)
    d %>% transmute(IID, key = paste(gene, k, sep = "|"),
                    state = ifelse(hap == k, "dark", ifelse(hap == "0", "pale", "grey")))))
}))
cells <- bind_rows(cell_a, cell_b) %>%
  left_join(select(AX$acc, IID, x), by = "IID") %>% left_join(select(rows, key, panel, y), by = "key") %>%
  mutate(state = factor(state, levels = STATE))
stopifnot(!anyNA(cells$x), !anyNA(cells$y))

bar_panel <- function(p, site_labels = FALSE) {
  r <- rows %>% filter(panel == p); d <- cells %>% filter(panel == p)
  lab <- r %>% filter(type != "header")
  ggplot(d, aes(x = x, y = y, fill = state)) +
    geom_tile(width = 1, height = 0.78) +
    hdr_layer(r, AX$xlim[1] + 0.5) +
    scale_fill_manual(values = c(dark = CELL_MINOR, grey = CELL_MAJOR, pale = CELL_NOCALL), labels = LEG, drop = FALSE,
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_y_continuous(breaks = lab$y, labels = lab$label, expand = expansion(add = c(0.35, 0.8))) +   # top room for the 9-pt header
    scale_x_continuous(limits = AX$xlim, expand = c(0, 0),
                       breaks = if (site_labels) AX$sites$xmid else NULL,
                       labels = if (site_labels) AX$sites$location else NULL) +
    labs(x = NULL, y = NULL, tag = p) + theme_fig()
}
band <- ggplot(AX$regions) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 1), fill = AX$regions$fill) +
  geom_text(aes(x = xmid, y = 0.5, label = lab), colour = AX$regions$txt, family = FONT,
            fontface = "bold", size = PT / .pt) +
  scale_x_continuous(limits = AX$xlim, expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(1, 1, 0.8, 1, "mm"))
pa <- bar_panel("a") + theme(legend.position = "none")
pb <- bar_panel("b", site_labels = TRUE) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "bottom", legend.justification = "center",
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-1, 0, 0, 0, "mm"))

na <- sum(rows$panel == "a"); nb <- sum(rows$panel == "b")
ROW_MM <- 3.0
fig <- band / pa / pb + plot_layout(heights = unit(c(5, na * ROW_MM, nb * ROW_MM), "mm"))
render_fig(fig, "Fig_carrier_origin", h_mm = 5 + (na + nb) * ROW_MM + 53)
