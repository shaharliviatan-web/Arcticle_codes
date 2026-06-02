#!/usr/bin/env Rscript
# 05_publication_plots.R
# Publication-quality plots (PNG 300 dpi + vector PDF) for a CURATED final gene list,
# after the user's biological filtering. Reuses the cached HapObjects from the run
# (no re-haplotyping).
#
# Input: a curated TSV (default 00_config/final_genes.tsv) with columns:
#   trait  gene_id  MGmin  epsilon   [title_override optional]
# Output (per gene) under 05_results/publication_plots/<trait>/<gene_id>/:
#   <gene>__tree.(png|pdf)
#   <gene>__crosshap_E<eps>.(png|pdf)
#   <gene>__violin_E<eps>.(png|pdf)
#   <gene>__heatmap_E<eps>.pdf            (one page per haplotype group; vector)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({
  library(yaml); library(data.table); library(ggplot2); library(crosshap)
})

argv <- commandArgs(trailingOnly = FALSE)
file_arg <- argv[grepl("^--file=", argv)]
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "plot_combined_pdf.R"))
source(file.path(script_dir, "R", "plot_heatmaps.R"))

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)
final_path <- if (length(args) >= 2) args[2] else file.path(cfg$output_root, "00_config", "final_genes.tsv")

if (!file.exists(final_path))
  stop("Curated final-gene list not found: ", final_path,
       "\nCreate it with columns: trait  gene_id  MGmin  epsilon  [title_override]")

run_root <- file.path(cfg$output_root, "04_runs", cfg$run_id)
out_root <- file.path(cfg$output_root, "05_results", "publication_plots")
gw <- as.data.table(read_gene_windows(file.path(cfg$output_root, "00_config", "gene_windows.tsv")))

final <- fread(final_path)
stopifnot(all(c("trait","gene_id","MGmin","epsilon") %in% names(final)))

save_gg <- function(plot, base, w, h) {
  ok <- tryCatch({
    ggplot2::ggsave(paste0(base, ".png"), plot = plot, width = w, height = h, dpi = 300, bg = "white")
    ggplot2::ggsave(paste0(base, ".pdf"), plot = plot, width = w, height = h, device = cairo_pdf)
    TRUE
  }, error = function(e) { message("  save failed for ", base, ": ", conditionMessage(e)); FALSE })
  ok
}

for (k in seq_len(nrow(final))) {
  trait <- final$trait[k]; gene_id <- final$gene_id[k]
  MGmin <- as.integer(final$MGmin[k]); eps <- final$epsilon[k]
  gene_file <- paste0(gene_id, ".vcf.gz")
  ttl <- if ("title_override" %in% names(final) && !is.na(final$title_override[k]) &&
             nzchar(final$title_override[k])) final$title_override[k] else {
    gr <- gw[trait == trait & gene_id == gene_id]
    if (nrow(gr) >= 1) gr$title[1] else paste0(gene_id, " | ", trait)
  }
  label <- paste0("Haplotypes_MGmin", MGmin, "_E", eps)

  cache_rds <- file.path(run_root, "Cache", trait, gene_id, paste0("MGmin_", MGmin), "HapObject.rds")
  if (!file.exists(cache_rds)) { message("MISSING cache: ", cache_rds); next }
  res <- readRDS(cache_rds); HapObject <- res$HapObject
  if (is.null(HapObject[[label]])) { message("No HapObject for ", label, " (", gene_id, ")"); next }

  outdir <- file.path(out_root, trait, gene_id); ensure_dir(outdir)
  base <- file.path(outdir, gene_id)
  message(sprintf("[%d/%d] %s  MGmin=%d eps=%s", k, nrow(final), gene_id, MGmin, eps))

  # tree
  tryCatch({
    if ("clustree_viz" %in% getNamespaceExports("crosshap")) {
      tp <- crosshap::clustree_viz(HapObject, type = "hap") + ggtitle(paste0(ttl, " | CrossHap tree"))
      save_gg(tp, paste0(base, "__tree"), 8, 8)
    }
  }, error = function(e) message("  tree failed: ", conditionMessage(e)))

  # crosshap viz
  tryCatch({
    vz <- suppressWarnings(suppressMessages(crosshap::crosshap_viz(HapObject, epsilon = as.numeric(eps))))
    save_gg(vz, paste0(base, "__crosshap_E", eps_tag(eps)), 12, 8)
  }, error = function(e) message("  crosshap_viz failed: ", conditionMessage(e)))

  # violin (largest-group-vs-others, Holm)
  tryCatch({
    vio <- plot_violin_with_holm(HapObject, label, title_text = paste0(ttl, " | eps=", eps))
    save_gg(vio, paste0(base, "__violin_E", eps_tag(eps)), 8, 6)
  }, error = function(e) message("  violin failed: ", conditionMessage(e)))

  # heatmap (vector PDF; one page per haplotype group)
  tryCatch({
    write_heatmaps_for_eps(HapObject = HapObject, label = label, gene_file = gene_file,
      title = ttl, trait = trait, eps = eps, MGmin = MGmin, raw_path = res$raw_path,
      common_ids = res$common_ids, vcf_raw_ids = res$vcf_raw_ids,
      out_pdf = paste0(base, "__heatmap_E", eps_tag(eps), ".pdf"))
  }, error = function(e) message("  heatmap failed: ", conditionMessage(e)))
}

cat(sprintf("Publication plots written under %s for %d gene(s).\n", out_root, nrow(final)))
