# =============================================================================
# 04_figures.R -- one figure per gene per SNP-matching version
# =============================================================================
# Layout, as specified:
#
#   TOP     violin plot of the wild haplotype groups (phenotype BLUP ~ group),
#           boxplot inside, x-axis label "A" over "(n=85)". Group n only --
#           no mean marker, no significance brackets.
#   BOTTOM  aligned genotype barcodes, one row per entity, drawn top to bottom:
#             - one row per wild haplotype group (per-SNP majority consensus)
#             - a gap
#             - one row per elite line
#           Row labels name the group or the cultivar. Columns are SNPs in
#           GENOMIC ORDER, shared between panels so the two read as one figure.
#
# Pairwise significance brackets were added on 2026-09-18, on request, so these
# figures carry the same test as step 04's combined PDFs. The method is ported
# verbatim from 04_.../01_scripts/R/plot_combined_pdf.R so the two agree:
#   * Wilcoxon rank-sum of EVERY group against the LARGEST group (deterministic
#     reference; not all pairs, which would multiply the comparisons)
#   * Holm correction across those k-1 comparisons, within the gene
#   * symbols  **** <= 1e-4, *** <= 1e-3, ** <= 0.01, * <= 0.05, else ns
#   * ns is shown, not hidden, as in step 04
# The brackets are drawn with plain geom_segment/geom_text rather than
# ggpubr::stat_pvalue_manual: ggpubr 0.6.3 errors against ggplot2 4.0.1 here
# (invalid NULL fontsize passed to grid::gpar). The statistics are unchanged.
# The Kruskal-Wallis p in the subtitle is the overall test; the brackets are the
# follow-up. Values are also written to Table_pairwise_group_tests.tsv.
#
# Deliberately NOT drawn (decided 2026-09-10): gene-model strip, crosshap marker
# group annotation bar, group means, per-elite "% identity to group" column.
#
# Elite lines are NOT assigned to a haplotype group anywhere in this figure --
# crosshap never saw them. The comparison is left to the reader's eye, which is
# the whole point of aligning the barcodes.
#
# Colours are inherited from step 04's plot_heatmaps.R so these figures read the
# same as the per-gene heatmaps already in the paper, plus three new states:
#   HET        elite only; the wild VCFs contain no heterozygous call (verified)
#   no record  version `filled_marked` only; a wild SNP with no elite record
#   Triallelic filled versions only; an elite line carrying an allele that is
#              neither the wild REF nor the wild ALT (added 2026-09-14)
# The legend lists only the states actually drawn in each figure, so e.g. the
# triallelic colour appears only if some line really carries the third allele.
#
# Rows are separated by a white gap -- tiles are drawn at TILE_HEIGHT of the row
# (added 2026-09-14) -- so each row reads as its own barcode, not one heatmap.
#
# Output -> results/figures/<version>/<trait>__<short_name>__<gene_id>.{pdf,png}
#
# Author : Shahar Liviatan
# Created: 2026-09-10
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
})
source(file.path(script_dir, "_load_params.R"))
P <- load_params(file.path(script_dir, "..", "config", "params.sh"))
Sys.setenv(TMPDIR = P$TMPDIR)

genes <- target_genes_df(P)
elite <- read_elite_lines(P)
msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
pw_tbl <- list()

# functional names for the titles, from step 05's paper table
gene_titles <- local({
  if (!file.exists(P$STEP05_PAPER_TABLE)) return(setNames(genes$short_name, genes$gene_id))
  t5 <- read.delim(P$STEP05_PAPER_TABLE, stringsAsFactors = FALSE)
  cn <- intersect(c("final_call", "annotation", "description"), names(t5))
  if (!length(cn)) return(setNames(genes$short_name, genes$gene_id))
  setNames(t5[[cn[1]]], t5$gene_id)
})

STATE_LEVELS <- c("Reference", "Alternate", "Heterozygous",
                  "Triallelic (elite-only allele)", "Missing", "No elite record")
state_cols <- c(Reference        = P$COL_REF,
                Alternate        = P$COL_ALT,
                Heterozygous     = P$COL_HET,
                `Triallelic (elite-only allele)` = P$COL_TRIALLELIC,
                Missing          = P$COL_MISS,
                `No elite record`= P$COL_NORECORD)

code_to_state <- function(x) {
  s <- rep("Missing", length(x))
  s[!is.na(x) & x == 0] <- "Reference"
  s[!is.na(x) & x == 1] <- "Alternate"
  s[!is.na(x) & x == 2] <- "Heterozygous"
  s[!is.na(x) & x == 3] <- "No elite record"
  s[!is.na(x) & x == 4] <- "Triallelic (elite-only allele)"
  factor(s, levels = STATE_LEVELS)
}

# --- pairwise test, ported from step 04's plot_combined_pdf.R -----------------
p_to_signif_symbol <- function(p) {
  if (is.na(p)) return("ns")
  if (p <= 0.0001) return("****")
  if (p <= 0.001)  return("***")
  if (p <= 0.01)   return("**")
  if (p <= 0.05)   return("*")
  "ns"
}

# Wilcoxon of every group vs the LARGEST group, Holm-corrected across k-1 tests.
pairwise_vs_largest <- function(d) {
  groups <- sort(unique(d$hap))
  if (length(groups) <= 1) return(NULL)
  sizes <- table(d$hap)
  ref <- names(sizes)[which.max(sizes)]          # largest group, deterministic
  others <- setdiff(groups, ref)
  if (!length(others)) return(NULL)
  raw_p <- vapply(others, function(g)
    tryCatch(wilcox.test(d$Pheno[d$hap == ref], d$Pheno[d$hap == g],
                         exact = FALSE)$p.value, error = function(e) NA_real_),
    numeric(1))
  p_adj <- p.adjust(raw_p, method = "holm")
  y_max <- max(d$Pheno, na.rm = TRUE); y_min <- min(d$Pheno, na.rm = TRUE)
  y_span <- max(1e-8, y_max - y_min)
  data.frame(group1 = ref, group2 = others, n_ref = as.integer(sizes[[ref]]),
             n_other = as.integer(sizes[others]), p = raw_p, p.adj = p_adj,
             p.adj.signif = vapply(p_adj, p_to_signif_symbol, character(1)),
             y.position = y_max + y_span * (0.10 + seq_along(others) * 0.08),
             stringsAsFactors = FALSE)
}

theme_pub <- function(base = 11) {
  theme_bw(base_size = base, base_family = "sans") +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
          axis.text = element_text(colour = "black"),
          plot.title = element_blank(), plot.subtitle = element_blank())
}

for (ver in P$VERSIONS) {
  outdir <- file.path(P$DIR_FIGURES, ver)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(genes))) {
    gid <- genes$gene_id[i]; short <- genes$short_name[i]
    f <- file.path(P$DIR_MATRICES, paste0(gid, "__", ver, ".rds"))
    if (!file.exists(f)) { msg("missing matrix: ", basename(f)); next }
    D <- readRDS(f)

    grp_lvls <- rownames(D$consensus)                        # A, B, C, ...
    gs <- D$group_stats
    gs <- gs[match(grp_lvls, gs$hap), ]

    # ---------------- TOP: violin ------------------------------------------
    pheno <- D$indfile
    pheno$hap <- factor(pheno$hap, levels = grp_lvls)
    pheno <- pheno[!is.na(pheno$hap) & !is.na(pheno$Pheno), ]

    xlabs <- setNames(sprintf("%s\n(n=%d)", gs$hap, gs$n), gs$hap)
    grp_fill <- setNames(
      grDevices::hcl.colors(length(grp_lvls), palette = "Dark 3")[seq_along(grp_lvls)],
      grp_lvls)

    p_violin <- ggplot(pheno, aes(x = hap, y = Pheno, fill = hap)) +
      geom_violin(trim = FALSE, colour = "black", linewidth = 0.35, width = 0.9) +
      geom_boxplot(width = 0.14, fill = "white", colour = "black",
                   outlier.size = 0.6, linewidth = 0.35) +
      scale_fill_manual(values = grp_fill, guide = "none") +
      scale_x_discrete(labels = xlabs) +
      labs(x = "Wild haplotype group", y = paste0(D$trait, " BLUP")) +
      theme_pub()

    # pairwise brackets (same test as step 04; identical across the three versions)
    pw <- pairwise_vs_largest(pheno)
    if (!is.null(pw) && nrow(pw)) {
      tip <- diff(range(pheno$Pheno, na.rm = TRUE)) * 0.012
      br <- data.frame(x = match(pw$group1, grp_lvls), xend = match(pw$group2, grp_lvls),
                       y = pw$y.position, lab = pw$p.adj.signif, stringsAsFactors = FALSE)
      p_violin <- p_violin +
        geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y),
                     inherit.aes = FALSE, linewidth = 0.3) +
        geom_segment(data = br, aes(x = x, xend = x, y = y - tip, yend = y),
                     inherit.aes = FALSE, linewidth = 0.3) +
        geom_segment(data = br, aes(x = xend, xend = xend, y = y - tip, yend = y),
                     inherit.aes = FALSE, linewidth = 0.3) +
        geom_text(data = br, aes(x = (x + xend) / 2, y = y + tip * 0.6, label = lab),
                  inherit.aes = FALSE, size = 3.2, vjust = 0) +
        expand_limits(y = max(pw$y.position) + tip * 4)
      if (identical(ver, P$VERSIONS[1]))
        pw_tbl[[gid]] <- cbind(gene_id = gid, short_name = short, trait = D$trait,
                               test = "Wilcoxon vs largest group, Holm", pw)
    }

    # ---------------- BOTTOM: aligned barcodes ------------------------------
    cons_df <- as.data.frame(as.table(D$consensus), stringsAsFactors = FALSE)
    names(cons_df) <- c("row_id", "site_key", "code")
    cons_df$block <- "Wild haplotype groups"
    cons_df$label <- paste0("Group ", cons_df$row_id)

    el_df <- as.data.frame(as.table(D$elite), stringsAsFactors = FALSE)
    names(el_df) <- c("row_id", "site_key", "code")
    el_df$block <- "Elite lines"
    el_df$label <- el_df$row_id

    bar <- rbind(cons_df, el_df)
    bar$state <- code_to_state(bar$code)
    bar$pos <- as.integer(sub("^[^:]+:([0-9]+):.*$", "\\1", bar$site_key))

    # column order = genomic position; row order = groups then elites, with a gap
    site_order <- unique(bar$site_key[order(bar$pos)])
    bar$site_key <- factor(bar$site_key, levels = site_order)

    grp_labels   <- paste0("Group ", grp_lvls)
    elite_labels <- rownames(D$elite)
    SPACER <- "​"                                   # zero-width, draws blank
    row_order <- c(rev(elite_labels), SPACER, rev(grp_labels))
    bar$label <- factor(bar$label, levels = row_order)

    # Tile borders are light grey, not white as in step 04's heatmaps: the
    # "No elite record" state of the `filled_marked` version is white, and on a
    # white ground with white borders it is invisible -- both in the panel and in
    # the legend key. A grey border makes the empty state readable as a state.
    p_bar <- ggplot(bar, aes(x = site_key, y = label, fill = state)) +
      geom_tile(colour = "grey78", linewidth = 0.25,
                height = as.numeric(P$TILE_HEIGHT)) +      # white gap between rows
      scale_fill_manual(values = state_cols, drop = TRUE, name = NULL) +  # present states only
      guides(fill = guide_legend(override.aes = list(colour = "grey40"))) +
      scale_y_discrete(drop = FALSE) +
      labs(x = sprintf("%s SNPs in %s:%s-%s (gene %s-%s, ±%d bp), genomic order",
                       length(site_order), D$chr,
                       format(D$win_start, big.mark = ","), format(D$win_end, big.mark = ","),
                       format(D$gene_start, big.mark = ","), format(D$gene_end, big.mark = ","),
                       P$WINDOW_BP),
           y = NULL) +
      theme_pub() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position = "bottom", legend.key.size = unit(0.4, "cm"),
            axis.text.y = element_text(size = 9))

    ttl <- gene_titles[[gid]]
    if (is.null(ttl) || is.na(ttl)) ttl <- short
    header <- sprintf("%s  |  %s  |  %s  (%s)", short, ttl, gid, D$trait)
    sub <- sprintf(
      "wild %d SNPs, elite %d records in window; shared %d, triallelic %d, wild-only %d, elite-only %d (not drawn); %d columns drawn  |  version: %s",
      D$n_shared + D$n_triallelic + D$n_wild_only, D$n_shared + D$n_triallelic + D$n_elite_only,
      D$n_shared, D$n_triallelic, D$n_wild_only, D$n_elite_only, ncol(D$consensus), ver)

    n_rows <- length(grp_lvls) + length(elite_labels) + 1
    fig <- (p_violin / p_bar) +
      plot_layout(heights = c(2.1, max(1, n_rows * 0.22))) +
      plot_annotation(title = header, subtitle = sub,
                      theme = theme(plot.title = element_text(face = "bold", size = 12),
                                    plot.subtitle = element_text(size = 8.5, colour = "grey25")))

    base <- file.path(outdir, sprintf("%s__%s__%s", D$trait, short, gid))
    h <- as.numeric(P$FIG_HEIGHT_IN) + max(0, n_rows - 8) * 0.16
    ggsave(paste0(base, ".pdf"), fig, width = as.numeric(P$FIG_WIDTH_IN), height = h,
           device = grDevices::cairo_pdf)
    ggsave(paste0(base, ".png"), fig, width = as.numeric(P$FIG_WIDTH_IN), height = h,
           dpi = as.numeric(P$FIG_DPI))
    msg("[", ver, "] ", short, " -> ", basename(base), ".{pdf,png}")
  }
}

if (length(pw_tbl)) {
  write.table(do.call(rbind, pw_tbl),
              file.path(P$DIR_TABLES, "Table_pairwise_group_tests.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  msg("pairwise tests -> ", file.path(P$DIR_TABLES, "Table_pairwise_group_tests.tsv"))
}

msg("DONE -> ", P$DIR_FIGURES)
