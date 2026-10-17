# =============================================================================
# 03_distributions.R -- check 02, item 3 (user request 2026-10-01): distributions of the four nutritional traits
# and the six morphological and phenological traits of Results Ch. 1 (flowering time, plant height, spike length,
# number of tillers, grain weight, number of grains):
#   Fig. E6  raw values by season (frequency polygons, one line per season: colour + line type);
#            x axis = the window of the step-00 +/-3 SD limits of the trait (all seasons); values outside it are
#            counted on the panel (they are the values the outlier rule removes; T1 gives them per season);
#   Fig. E7  the genotype BLUPs used in all analyses (290 accessions), one histogram per trait.
# Table T5: summary of the BLUPs.  Setup, panel and filters: 00_common.R.
# (Fig. E5 is check 01's block-correlation figure; numbering = the order of the user's list: boxplots,
#  correlations between seasons / traits / blocks, distributions.)
# Run: /usr/bin/Rscript scripts/03_distributions.R (from the check dir)
# =============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
cat("[03] distributions\n")

d   <- load_raw(keep_block4_2019 = FALSE)
lim <- sd3_limits(d)
long <- d %>% select(acc, season_lab, all_of(TRAITS$col)) %>%
  pivot_longer(all_of(TRAITS$col), names_to = "col", values_to = "value") %>% filter(!is.na(value)) %>%
  mutate(trait = factor(TRAITS$label[match(col, TRAITS$col)], levels = TRAITS$label),
         short = TRAITS$short[match(col, TRAITS$col)])
win <- lim %>% group_by(trait) %>% summarise(lo = min(lower), hi = max(upper), .groups = "drop") %>%
  mutate(trait = factor(TRAITS$label[match(trait, TRAITS$short)], levels = TRAITS$label))
long <- long %>% left_join(win, by = "trait") %>% mutate(inside = value >= lo & value <= hi)
outside <- long %>% group_by(trait) %>% summarise(n_out = sum(!inside), .groups = "drop")
inside  <- long %>% filter(inside)
pE5 <- ggplot(inside, aes(value, colour = season_lab, linetype = season_lab)) +
  geom_freqpoly(bins = 40, linewidth = 0.4) +
  geom_text(data = outside, aes(x = Inf, y = Inf, label = ifelse(n_out > 0, paste0(n_out, " outside"), "")),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, family = FONT, size = (PT - 1) / .pt) +
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  scale_colour_manual(values = SEASON_COL, name = "Season") + scale_linetype_manual(values = SEASON_LTY, name = "Season") +
  labs(x = "Raw value", y = "Plants") + theme_tag()
save_fig(pE5, "Fig_E6_raw_distributions_by_season", 215)

b <- read.csv(BLUP_CSV, check.names = FALSE)
stopifnot(nrow(b) == 290)
bl <- b %>% pivot_longer(all_of(TRAITS$blup), names_to = "col", values_to = "value") %>%
  mutate(trait = factor(TRAITS$label[match(col, TRAITS$blup)], levels = TRAITS$label))
T5 <- bl %>% group_by(trait) %>% summarise(n = sum(!is.na(value)), mean = mean(value), sd = sd(value), min = min(value),
                                           median = median(value), max = max(value), .groups = "drop")
write_tsv(as.data.frame(T5), "T5_BLUP_summary")
pE6 <- ggplot(bl, aes(value)) + geom_histogram(bins = 30, fill = "grey75", colour = "grey30", linewidth = 0.2) +
  facet_wrap(~ trait, ncol = 5, scales = "free", labeller = label_wrap_gen(width = 15)) +
  labs(x = "Genotype BLUP (centered; units of the trait)", y = "Accessions") +
  theme_tag() + theme(strip.text = element_text(size = PT - 1))
save_fig(pE6, "Fig_E7_BLUP_distributions", 90)
cat("[03] done\n")
