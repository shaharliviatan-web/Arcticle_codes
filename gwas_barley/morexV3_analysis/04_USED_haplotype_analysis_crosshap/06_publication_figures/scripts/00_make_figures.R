#!/usr/bin/env Rscript
# 00_make_figures.R
# Publication haplotype figures for the 4 final candidate genes (wild-only, 290
# accessions; NO elite/cultivated lines). For each gene: a stacked figure with a
# violin (phenotype per haplotype group, Kruskal-Wallis + Wilcoxon/Holm vs the
# largest group) on top and a one-representative-row-per-group genotype heatmap
# below. Each gene also gets its own comprehensive stats CSVs so the results
# chapter can be written entirely from the tables. Outputs PNG (600 dpi) + vector
# PDF (cairo, embedded fonts) per gene, plus a 2x2 composite with (a)-(d) tags.
#
# Inputs : config/figure_genes.tsv  (serial_no, panel_tag, trait, gene_id,
#            gene_symbol, display_name, y_label, MGmin, epsilon)
#          step-04 cache HapObject.rds per gene (haplotyping + raw VCF path + IDs)
#          step-04 Stats/gene_shortlist.csv (gene-level FDR/Bonferroni + window info)
# Outputs: results/figures/  results/tables/

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(data.table); library(dplyr); library(patchwork) })

BASE   <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap"
HERE   <- file.path(BASE, "06_publication_figures")
source(file.path(HERE, "scripts", "R", "figure_utils.R"))

RUN_ID   <- "candidate_genes_1000bp_V1"
CACHE    <- file.path(BASE, "04_runs", RUN_ID, "Cache")
SHORTLIST<- file.path(BASE, "04_runs", RUN_ID, "Stats", "gene_shortlist.csv")
FIG_DIR  <- file.path(HERE, "results", "figures")
TAB_DIR  <- file.path(HERE, "results", "tables")
ensure_dir(FIG_DIR); ensure_dir(TAB_DIR)

cfg <- fread(file.path(HERE, "config", "figure_genes.tsv"))
sl  <- if (file.exists(SHORTLIST)) fread(SHORTLIST) else data.table()

PNG_W <- 13; PNG_H <- 10; PNG_DPI <- 600   # inches / dpi for standalone figures
HEATMAP_REL_H <- 1.0                        # heatmap height relative to violin = 3 (compact strip)

combined_plots <- list()
all_summary <- list(); all_groups <- list()

for (i in seq_len(nrow(cfg))) {
  g  <- cfg[i]
  lab <- paste0("Haplotypes_MGmin", g$MGmin, "_E", g$epsilon)
  cache_file <- file.path(CACHE, g$trait, g$gene_id, paste0("MGmin_", g$MGmin), "HapObject.rds")
  if (!file.exists(cache_file)) stop("Missing cache: ", cache_file)
  o  <- readRDS(cache_file)
  ho <- o$HapObject[[lab]]
  if (is.null(ho)) stop("Label not in cache: ", lab, " (", g$gene_id, ")")

  indfile <- ho$Indfile
  hap_map <- as.data.frame(indfile)[, c("Ind", "hap")]
  hap_map$hap <- as.character(hap_map$hap)

  # ---- Violin ----
  vio <- build_violin(indfile, title = g$display_name, y_label = g$y_label)

  # ---- Heatmap (reconstruct genotype matrix exactly as haplotyped) ----
  m01 <- build_geno_matrix(o$raw_path, o$common_ids, o$vcf_raw_ids, g$MGmin)
  het <- build_heatmap(m01, hap_map, vio$hap_levels)

  # ---- Stacked figure ----
  combined <- vio$plot / het$plot + plot_layout(heights = c(3.0, HEATMAP_REL_H))
  combined_plots[[i]] <- combined

  sym_safe <- gsub("[^A-Za-z0-9]+", "", g$gene_symbol)   # filesystem-safe (e.g. AP2/ERF -> AP2ERF)
  stub <- sprintf("%d_%s_%s", g$serial_no, sym_safe, sub("HORVU.MOREX.r3.", "", g$gene_id))
  ggplot2::ggsave(file.path(FIG_DIR, paste0(stub, ".png")), combined,
                  width = PNG_W, height = PNG_H, dpi = PNG_DPI, bg = "white", type = "cairo")
  ggplot2::ggsave(file.path(FIG_DIR, paste0(stub, ".pdf")), combined,
                  width = PNG_W, height = PNG_H, device = grDevices::cairo_pdf, bg = "white")

  # ---- Gene-level info from the step-04 shortlist (FDR/Bonferroni + window) ----
  s <- if (nrow(sl)) sl[trait == g$trait & gene_id == g$gene_id] else data.table()
  getv <- function(col) if (nrow(s) && col %in% names(s)) s[[col]][1] else NA

  summ <- data.frame(
    serial_no = g$serial_no, panel_tag = g$panel_tag, gene_symbol = g$gene_symbol,
    trait = g$trait, gene_id = g$gene_id, display_name = g$display_name,
    chr = getv("chr"), gene_start = getv("gene_start"), gene_end = getv("gene_end"),
    n_snps_window = ncol(m01), MGmin = g$MGmin, epsilon = g$epsilon,
    n_total = vio$n_total, n_groups = vio$n_groups, reference_group = vio$ref_hap,
    kw_statistic = vio$kw_stat, kw_df = vio$kw_df, kw_p_raw = vio$kw_p,
    epsilon_squared = vio$eps2,
    lead_SNP = getv("lead_SNP"), class = getv("class"), dist_to_lead_bp = getv("dist_to_lead_bp"),
    within_gene_holm_p = getv("best_kw_p_holm_within_gene"),
    trait_fdr_p = getv("trait_fdr_p"), significant_fdr = getv("significant_fdr"),
    trait_bonferroni_p = getv("trait_bonferroni_p"), significant_bonferroni = getv("significant_bonferroni"),
    stringsAsFactors = FALSE)

  # ---- Per-group table: distribution + pairwise test + heatmap representative ----
  grp <- vio$gstats
  grp$is_reference <- grp$hap == vio$ref_hap
  pw <- vio$pairwise
  grp$wilcox_vs_ref_p_raw  <- if (nrow(pw)) pw$p_raw[match(grp$hap, pw$hap)] else NA_real_
  grp$wilcox_vs_ref_p_holm <- if (nrow(pw)) pw$p_holm[match(grp$hap, pw$hap)] else NA_real_
  grp$signif_stars         <- if (nrow(pw)) pw$p.signif[match(grp$hap, pw$hap)] else NA_character_
  grp <- merge(grp, het$rep_info, by = "hap", all.x = TRUE, sort = FALSE)
  grp <- grp[order(match(grp$hap, vio$hap_levels)), ]
  grp <- cbind(serial_no = g$serial_no, gene_symbol = g$gene_symbol, trait = g$trait,
               gene_id = g$gene_id, grp, stringsAsFactors = FALSE)

  fwrite(summ, file.path(TAB_DIR, paste0(stub, "_summary.csv")))
  fwrite(grp,  file.path(TAB_DIR, paste0(stub, "_groups.csv")))
  all_summary[[i]] <- summ; all_groups[[i]] <- grp

  cat(sprintf("[%d/%d] %-5s %s  KW P=%.2e  eps2=%.3f  groups=%d (ref=%s)  SNPs=%d\n",
              i, nrow(cfg), g$gene_symbol, g$gene_id, vio$kw_p, vio$eps2,
              vio$n_groups, vio$ref_hap, ncol(m01)))
}

# ---- Combined masters ----
fwrite(rbindlist(all_summary, fill = TRUE), file.path(TAB_DIR, "all_genes_summary.csv"))
fwrite(rbindlist(all_groups,  fill = TRUE), file.path(TAB_DIR, "all_genes_groups.csv"))

# ---- 2x2 composite with (a)-(d) tags ----
wrapped <- lapply(combined_plots, patchwork::wrap_elements)
composite <- patchwork::wrap_plots(wrapped, ncol = 2) +
  patchwork::plot_annotation(tag_levels = "a") &
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 26, face = "bold"))
ggplot2::ggsave(file.path(FIG_DIR, "composite_4panel.png"), composite,
                width = 24, height = 22, dpi = 400, bg = "white", type = "cairo", limitsize = FALSE)
ggplot2::ggsave(file.path(FIG_DIR, "composite_4panel.pdf"), composite,
                width = 24, height = 22, device = grDevices::cairo_pdf, bg = "white", limitsize = FALSE)

cat("\nDone. Figures + tables in:\n  ", FIG_DIR, "\n  ", TAB_DIR, "\n")
