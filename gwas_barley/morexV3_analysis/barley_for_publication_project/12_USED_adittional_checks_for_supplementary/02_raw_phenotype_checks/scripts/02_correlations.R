# =============================================================================
# 02_correlations.R -- check 02, item 2 (user request 2026-10-01): Pearson correlations of the raw values of the
# four nutritional traits, from accession means:
#   (i)  between seasons, per trait: each season's accession means against the other seasons' and against the
#        all-season accession means (mean of the values centered within season, so that season differences and
#        missing seasons do not shift an accession);
#   (ii) between traits: over all seasons (accession means of the season-centered values) and within each season
#        (accession means of the raw values).
#   (iii) between blocks within each season: already computed by check 01 (01_block4_2019_interblock_correlations,
#        interblock_correlations_all.tsv and Fig_interblock_correlation_heatmaps.png); not recomputed here.
# Two-sided cor.test, complete accessions per pair; descriptive, no correction for multiple testing.
# Setup, panel and filters: 00_common.R.  Run: /usr/bin/Rscript scripts/02_correlations.R (from the check dir)
# =============================================================================
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "00_common.R"))
cat("[02] correlations\n")

d <- load_raw(keep_block4_2019 = FALSE)
d <- d[!is.na(d$ProteinAsis.), ]
for (v in NUTRI$col) {                                         # season-centered values (as step 00, raw scale)
  d[[paste0(v, "_c")]] <- d[[v]] - ave(d[[v]], d$season_lab, FUN = function(x) mean(x, na.rm = TRUE))
}
acc_season <- d %>% group_by(acc, season_lab) %>% summarise(across(all_of(NUTRI$col), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
acc_all    <- d %>% group_by(acc) %>% summarise(across(all_of(paste0(NUTRI$col, "_c")), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
names(acc_all) <- c("acc", NUTRI$col)
stopifnot(nrow(acc_all) == 290)
SCOPES <- c(SEASON_LABEL, "All seasons")
SCOPES <- unname(SCOPES)

ct <- function(x, y) { ok <- complete.cases(x, y); t <- cor.test(x[ok], y[ok]); c(r = unname(t$estimate), n = sum(ok), p = t$p.value) }

# ---- (i) between seasons, per trait ---------------------------------------------------------------------------
wide_by_scope <- function(v) {
  w <- acc_season %>% select(acc, season_lab, all_of(v)) %>% pivot_wider(names_from = season_lab, values_from = all_of(v))
  w <- full_join(w, acc_all %>% select(acc, `All seasons` = all_of(v)), by = "acc")
  w[, c("acc", SCOPES)]
}
T3 <- do.call(rbind, lapply(seq_len(nrow(NUTRI)), function(i) {
  w <- wide_by_scope(NUTRI$col[i])
  do.call(rbind, combn(SCOPES, 2, simplify = FALSE, FUN = function(pr) {
    s <- ct(w[[pr[1]]], w[[pr[2]]])
    data.frame(trait = NUTRI$short[i], scope_1 = pr[1], scope_2 = pr[2], r = s["r"], n = s["n"], p = s["p"])
  }))
}))
rownames(T3) <- NULL
write_tsv(T3, "T3_correlations_between_seasons")

# ---- (ii) between traits, over all seasons and per season -----------------------------------------------------
T4 <- do.call(rbind, lapply(SCOPES, function(sc) {
  w <- if (sc == "All seasons") acc_all else acc_season[acc_season$season_lab == sc, ]
  do.call(rbind, combn(seq_len(nrow(NUTRI)), 2, simplify = FALSE, FUN = function(ij) {
    s <- ct(w[[NUTRI$col[ij[1]]]], w[[NUTRI$col[ij[2]]]])
    data.frame(scope = sc, trait_1 = NUTRI$short[ij[1]], trait_2 = NUTRI$short[ij[2]], r = s["r"], n = s["n"], p = s["p"])
  }))
}))
rownames(T4) <- NULL
write_tsv(T4, "T4_correlations_between_traits")

# ---- figures: lower-triangle heatmaps ---------------------------------------------------------------------------
heat <- function(df, levels, lab_x = NULL, lab_y = NULL) {
  # lower triangle without its empty first row and last column
  df$a <- factor(df$a, levels = levels[-length(levels)]); df$b <- factor(df$b, levels = rev(levels[-1]))
  ggplot(df, aes(a, b, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > 0.6), family = FONT, size = (PT - 1) / .pt, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
    scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", limits = c(-1, 1), name = "Pearson r") +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    coord_equal() + labs(x = lab_x, y = lab_y) + theme_tag() +
    theme(panel.grid.major = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1))
}
pE3 <- wrap_plots(lapply(NUTRI$short, function(tr) {
  s <- T3[T3$trait == tr, ]
  heat(data.frame(a = s$scope_1, b = s$scope_2, r = s$r), SCOPES) + labs(subtitle = tr) +
    theme(plot.subtitle = element_text(face = "bold", hjust = 0.5, size = PT))
}), ncol = 2) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
save_fig(pE3, "Fig_E3_correlations_between_seasons", 150)
pE4 <- wrap_plots(lapply(SCOPES[c(4, 1, 2, 3)], function(sc) {
  s <- T4[T4$scope == sc, ]
  heat(data.frame(a = s$trait_1, b = s$trait_2, r = s$r), NUTRI$short) + labs(subtitle = sc) +
    theme(plot.subtitle = element_text(face = "bold", hjust = 0.5, size = PT))
}), ncol = 2) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
save_fig(pE4, "Fig_E4_correlations_between_traits", 150)
cat("[02] done\n")
