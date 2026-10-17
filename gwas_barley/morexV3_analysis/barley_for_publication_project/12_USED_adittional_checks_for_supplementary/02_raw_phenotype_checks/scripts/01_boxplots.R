# =============================================================================
# 01_boxplots.R -- check 02, item 1 (user request 2026-10-01): boxplots of the raw values of the four
# nutritional traits (a) by season, all seasons together; (b) by block within each season.
# Also the summary tables behind them (T1 by season, all ten traits; T2 by block, nutritional traits).
# Setup, panel and filters: 00_common.R.  Run: /usr/bin/Rscript scripts/01_boxplots.R (from the check dir)
# =============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
cat("[01] boxplots\n")

d   <- load_raw(keep_block4_2019 = FALSE)                 # analysis panel
lim <- sd3_limits(d)

# ---- checks against the manuscript (M&M, Phenotypic data analysis) ----------------------------------------------
stopifnot(length(unique(d$acc)) == 290)
stopifnot(identical(as.integer(table(d$season_lab)), c(870L, 1160L, 1450L)))             # plants in the panel
nir <- !is.na(d$ProteinAsis.)
stopifnot(identical(as.integer(table(d$season_lab[nir])), c(862L, 1139L, 1163L)))         # plants with NIR values
stopifnot(all(!is.na(d$StarchAsis.[nir]) & !is.na(d$BetaglucansAsis.[nir]) & !is.na(d$FiberAsis.[nir])))

# ---- T1: summary by trait and season (all ten traits) -----------------------------------------------------------
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
T1 <- do.call(rbind, lapply(seq_len(nrow(TRAITS)), function(i) {
  v <- TRAITS$col[i]
  do.call(rbind, lapply(levels(d$season_lab), function(s) {
    x <- d[[v]][d$season_lab == s]; a <- d$acc[d$season_lab == s][!is.na(x)]; x <- x[!is.na(x)]
    if (!length(x)) return(NULL)
    L <- lim[lim$trait == TRAITS$short[i] & lim$season == s, ]
    data.frame(trait = TRAITS$short[i], unit = TRAITS$unit[i], season = s, n_plants = length(x),
               n_accessions = length(unique(a)), mean = mean(x), sd = sd(x), cv_pct = 100 * sd(x) / mean(x),
               min = min(x), q1 = q(x, .25), median = median(x), q3 = q(x, .75), max = max(x),
               sd3_lower = L$lower, sd3_upper = L$upper, n_below_sd3 = sum(x < L$lower), n_above_sd3 = sum(x > L$upper))
  }))
}))
write_tsv(T1, "T1_summary_by_season")

# ---- Figure E1: nutritional traits by season ----------------------------------------------------------------------
long <- d[nir, ] %>% select(acc, season_lab, all_of(NUTRI$col)) %>%
  pivot_longer(all_of(NUTRI$col), names_to = "col", values_to = "value") %>%
  mutate(trait = factor(NUTRI$label[match(col, NUTRI$col)], levels = NUTRI$label))
limN <- lim %>% filter(trait %in% NUTRI$short) %>%
  mutate(trait = factor(NUTRI$label[match(trait, NUTRI$short)], levels = NUTRI$label),
         season_lab = factor(season, levels = SEASON_LABEL))
nlab <- long %>% count(trait, season_lab)
panel_e1 <- function(tr) {
  nn <- filter(nlab, trait == tr)
  L  <- filter(limN, trait == tr) %>% mutate(x = as.numeric(season_lab))
  ggplot(filter(long, trait == tr), aes(season_lab, value)) +
    geom_boxplot(fill = "grey90", colour = "grey20", linewidth = 0.3, outlier.size = 0.4, outlier.alpha = 0.6, width = 0.6) +
    geom_segment(data = L, aes(x = x - 0.38, xend = x + 0.38, y = lower, yend = lower), inherit.aes = FALSE,
                 linetype = "dashed", colour = "#B2182B", linewidth = 0.35) +
    geom_segment(data = L, aes(x = x - 0.38, xend = x + 0.38, y = upper, yend = upper), inherit.aes = FALSE,
                 linetype = "dashed", colour = "#B2182B", linewidth = 0.35) +
    scale_x_discrete(labels = setNames(paste0(nn$season_lab, "\nn = ", nn$n), as.character(nn$season_lab))) +
    labs(x = NULL, y = tr) + theme_tag()
}
pE1 <- wrap_plots(lapply(NUTRI$label, panel_e1), ncol = 2) + plot_annotation(tag_levels = "a")
save_fig(pE1, "Fig_E1_nutritional_traits_by_season", 130)

# ---- T2 + Figure E2: nutritional traits by block within season (block 4 of 2019-20 shown) -------------------------
db <- load_raw(keep_block4_2019 = TRUE)
db <- db[!is.na(db$ProteinAsis.), ]
longb <- db %>% select(acc, season_lab, Block_no, block4_2019, all_of(NUTRI$col)) %>%
  pivot_longer(all_of(NUTRI$col), names_to = "col", values_to = "value") %>%
  mutate(trait = factor(NUTRI$label[match(col, NUTRI$col)], levels = NUTRI$label),
         Block = factor(Block_no), status = ifelse(block4_2019, "Excluded from the analyses", "Analyzed"))
T2 <- longb %>% group_by(trait, season = season_lab, block = Block_no, status) %>%
  summarise(n_plants = sum(!is.na(value)), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE), min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE),
            .groups = "drop") %>% mutate(trait = sub(" \\(%\\)", "", as.character(trait)))
write_tsv(as.data.frame(T2), "T2_summary_by_block")
pE2 <- ggplot(longb, aes(Block, value, fill = status)) +
  geom_boxplot(colour = "grey20", linewidth = 0.3, outlier.size = 0.3, outlier.alpha = 0.6, width = 0.65) +
  facet_grid(trait ~ season_lab, scales = "free", space = "free_x", switch = "y") +
  scale_fill_manual(values = c("Analyzed" = "grey90", "Excluded from the analyses" = "#F4A582"), name = NULL) +
  labs(x = "Block", y = NULL) + theme_tag() +
  theme(strip.placement = "outside", strip.background.y = element_blank(), strip.text.y = element_text(face = "plain"))
save_fig(pE2, "Fig_E2_nutritional_traits_by_block", 200)
cat("[01] done\n")
