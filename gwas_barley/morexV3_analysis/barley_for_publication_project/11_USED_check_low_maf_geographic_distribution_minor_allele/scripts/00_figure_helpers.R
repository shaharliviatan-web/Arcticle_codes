# 00_figure_helpers.R — row layout, accession axis, theme and export shared by the two
# figures (08, 09), so both draw the same rows with the same labels in the same order.
# Sourced after 00_config.R.

suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

MINUS <- "−"
fmt_signed <- function(x, d = 2) gsub("-", MINUS, sprintf(paste0("%+.", d, "f"), x))

# ---- rows ------------------------------------------------------------------------------------
# Panel a: one header row per trait, then its leads in genomic order.
# Panel b: one header row per gene, then its haplotype groups (A, B, ...).
# y runs downwards (row 1 = y -1) within each panel.
build_rows <- function(t02, t04) {
  a <- bind_rows(lapply(TRAIT_ORDER, function(t) {
    l <- t02 %>% filter(trait == t) %>% mutate(gene = coalesce(gene, ""))   # read_tsv turns "" into NA
    bind_rows(
      tibble(type = "header", trait = t, key = paste0("hdr_", t),
             label = paste0(TRAIT_LAB[t], " (", nrow(l), if (nrow(l) == 1) " locus)" else " loci)")),
      tibble(type = "lead", trait = t, key = l$lead_SNP,
             label = sprintf("%s:%.2f  %.2f  %s%s", sub(":.*", "", l$lead_SNP),
                             as.numeric(sub(".*:", "", l$lead_SNP)) / 1e6, l$MAF,
                             ifelse(l$minor_effect == "raises", "+", MINUS),
                             ifelse(l$gene != "", paste0("  ", l$gene), ""))))
  })) %>% mutate(panel = "a", y = -seq_len(n()))
  b <- bind_rows(lapply(unique(t04$gene), function(gn) {
    x <- t04 %>% filter(gene == gn)
    tr <- if (x$trait[1] == "fiber; starch") "fiber and starch" else x$trait[1]
    mean_lab <- if (gn == "GDSL") paste0(fmt_signed(x$mean_BLUP_fiber), " / ", fmt_signed(x$mean_BLUP_starch))
                else fmt_signed(ifelse(is.na(x$mean_BLUP_fiber), x$mean_BLUP_starch, x$mean_BLUP_fiber))
    bind_rows(
      tibble(type = "header", trait = x$trait[1], key = paste0("hdr_", gn),
             label = paste0(if (gn == "GDSL") "GDSL (7H)" else gn, ", ", tr)),
      tibble(type = "group", trait = x$trait[1], key = paste(gn, x$group, sep = "|"),
             label = sprintf("%s  n = %d  %s", x$group, x$n, mean_lab)))
  })) %>% mutate(panel = "b", y = -seq_len(n()))
  bind_rows(a, b)
}

# ---- accession axis: region > site > accession, with gaps ---------------------------------------
build_axis <- function(acc) {
  acc <- acc %>% mutate(region = factor(region, levels = REGION_ORDER)) %>% arrange(region, site, IID)
  x <- numeric(nrow(acc)); pos <- 0
  for (i in seq_len(nrow(acc))) {
    if (i > 1) pos <- pos + 1 + (acc$site[i] != acc$site[i - 1]) * 1.4 + (acc$region[i] != acc$region[i - 1]) * 4
    x[i] <- pos
  }
  acc$x <- x
  list(acc = acc,
       sites = acc %>% group_by(region, site, location) %>%
         summarise(xmid = mean(x), .groups = "drop"),
       regions = acc %>% group_by(region) %>%
         summarise(xmid = mean(range(x)), xmin = min(x) - 0.5, xmax = max(x) + 0.5, .groups = "drop") %>%
         mutate(lab = REGION_SHORT[as.integer(region)], fill = PAL_REGION[as.integer(region)],
                txt = REGION_TXT[as.integer(region)]),
       xlim = c(min(x) - 1, max(x) + 1))
}

# ---- theme -------------------------------------------------------------------------------------
theme_fig <- function() {
  theme_minimal(base_size = PT, base_family = FONT) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = INK, size = PT), axis.title = element_text(colour = INK, size = PT),
          axis.ticks = element_blank(), legend.text = element_text(colour = INK, size = PT),
          legend.title = element_blank(), legend.key.size = unit(3, "mm"),
          plot.tag = element_text(family = FONT, face = "bold", size = PT + 2, colour = INK),
          plot.tag.position = c(0, 1), plot.margin = margin(1, 1, 1, 1, "mm"))
}
hdr_layer <- function(rows, x) {
  h <- rows %>% filter(type == "header")
  geom_label(data = h, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, vjust = 0.5,
             family = FONT, fontface = "bold", size = PT / .pt, colour = INK, fill = "white",
             linewidth = 0, label.padding = unit(0.4, "mm"), label.r = unit(0, "mm"))
}

# ---- export: PNG (docx build) + TIFF (submission), pixel-identical, as the 08 figures -----------
render_fig <- function(plot, name, h_mm) {
  stopifnot(h_mm <= H_MAX)
  png_f <- file.path(FIG, paste0(name, ".png")); tif_f <- file.path(FIG, paste0(name, ".tif"))
  ragg::agg_png(png_f, width = W_MM, height = h_mm, units = "mm", res = DPI, background = "white")
  print(plot); invisible(dev.off())
  stopifnot(system2("python3", c("-c", shQuote(sprintf(
    "from PIL import Image; Image.open('%s').convert('RGB').save('%s', compression='tiff_lzw', dpi=(%d,%d))",
    png_f, tif_f, DPI, DPI)))) == 0)
  cat(sprintf("[fig] %s.{png,tif}  %d x %.0f mm, %d dpi\n", name, W_MM, h_mm, DPI))
}
