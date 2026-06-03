#!/usr/bin/env Rscript
# 07_build_publication_table.R
# Build the master "publishing genes" table for the final biologically-filtered
# candidate genes. One row per curated gene, joining everything generated across
# steps 04 (haplotype stats + plot paths), 05 (Swiss-Prot evidence), 06 (InterPro/
# Pfam/GO). The per-gene MGmin/epsilon selected for the figures drive the displayed
# haplotype grouping; the FDR/Bonferroni significance is the gene-level value.
#
# Input : 00_config/final_genes.tsv  (serial_no, trait, gene_id, gene_symbol,
#          final_annotation, MGmin, epsilon)
# Output: 05_results/tables/publication_genes_table.csv

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages({ library(yaml); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
cfg_path <- if (length(args) >= 1) args[1] else
  "/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap/00_config/config.yaml"
cfg <- yaml::read_yaml(cfg_path)

base04   <- cfg$output_root
analysis <- dirname(base04)                                  # morexV3_analysis/
run_root <- file.path(base04, "04_runs", cfg$run_id)
cache_root <- file.path(run_root, "Cache")

final_path <- file.path(base04, "00_config", "final_genes.tsv")
sl_path    <- file.path(run_root, "Stats", "gene_shortlist.csv")
sp_path    <- file.path(analysis, "05_USED_gene_annotation", "results", "tables", "fdr_gene_annotation_swissprot.tsv")
ip_path    <- file.path(analysis, "06_USED_interpro_domains", "results", "tables", "fdr_gene_interpro.tsv")

stopifnot(file.exists(final_path), file.exists(sl_path))
fg <- fread(final_path)
sl <- fread(sl_path)
sp <- if (file.exists(sp_path)) fread(sp_path) else data.table()
ip <- if (file.exists(ip_path)) fread(ip_path) else data.table()

## ---- selected-epsilon haplotype stats from the cached HapObject ----
sel_stat <- function(trait, gene_id, MGmin, eps) {
  f <- file.path(cache_root, trait, gene_id, paste0("MGmin_", MGmin), "HapObject.rds")
  if (!file.exists(f)) return(list(n_groups = NA_integer_, sizes = NA_character_, kw = NA_real_))
  ho <- readRDS(f)$HapObject
  lab <- paste0("Haplotypes_MGmin", MGmin, "_E", eps)
  if (is.null(ho[[lab]])) return(list(n_groups = NA_integer_, sizes = NA_character_, kw = NA_real_))
  d <- as.data.table(ho[[lab]]$Indfile)
  d[, hap := as.character(hap)][, Pheno := as.numeric(Pheno)]
  d <- d[hap != "0" & !is.na(Pheno)]
  szs <- sort(as.integer(table(d$hap)))
  kw <- tryCatch(kruskal.test(Pheno ~ hap, data = d)$p.value, error = function(e) NA_real_)
  list(n_groups = length(szs), sizes = paste(szs, collapse = "|"), kw = kw)
}

rows <- vector("list", nrow(fg))
for (i in seq_len(nrow(fg))) {
  g <- fg$gene_id[i]; tr <- fg$trait[i]; M <- as.integer(fg$MGmin[i]); e <- fg$epsilon[i]
  s <- sl[trait == tr & gene_id == g]
  if (nrow(s) == 0) stop("Gene not found in shortlist: ", tr, " / ", g)
  st <- sel_stat(tr, g, M, e)
  spr <- if (nrow(sp)) sp[gene_id == g & trait == tr] else data.table()
  if (nrow(spr) == 0 && nrow(sp)) spr <- sp[gene_id == g][1]   # SP hit is gene-level
  ipr <- if (nrow(ip)) ip[gene_id == g][1] else data.table()

  rows[[i]] <- data.table(
    serial_no   = fg$serial_no[i],
    gene_symbol = fg$gene_symbol[i],
    trait       = tr,
    gene_id     = g,
    final_annotation = fg$final_annotation[i],
    chr         = s$chr, gene_start = s$gene_start, gene_end = s$gene_end, strand = s$strand,
    locus_id    = s$locus_id,
    lead_SNP    = s$lead_SNP, lead_pos = as.integer(sub(".*:", "", s$lead_SNP)),
    lead_neg_log10p = s$lead_neg_log10p, class = s$class, dist_to_lead_bp = s$dist_to_lead_bp,
    n_snps_window = s$n_snps_window,
    MGmin = M, epsilon = e,
    n_groups = st$n_groups, group_sizes = st$sizes, kw_p_raw = st$kw,
    holm_within_gene_p = s$best_kw_p_holm_within_gene,
    trait_fdr_p = s$trait_fdr_p, significant_fdr = s$significant_fdr,
    trait_bonferroni_p = s$trait_bonferroni_p, significant_bonferroni = s$significant_bonferroni,
    swissprot_name = if (nrow(spr)) spr$sp_name else NA_character_,
    swissprot_acc  = if (nrow(spr)) spr$sp_accession else NA_character_,
    swissprot_organism = if (nrow(spr)) spr$sp_organism else NA_character_,
    swissprot_pident = if (nrow(spr)) spr$pident else NA_real_,
    swissprot_qcov   = if (nrow(spr)) spr$qcovhsp else NA_real_,
    swissprot_evalue = if (nrow(spr)) spr$evalue else NA_real_,
    swissprot_status = if (nrow(spr)) spr$status else NA_character_,
    pfam_ids       = if (nrow(ipr)) ipr$pfam_ids else NA_character_,
    pfam_descs     = if (nrow(ipr)) ipr$pfam_descs else NA_character_,
    interpro_ids   = if (nrow(ipr)) ipr$interpro_ids else NA_character_,
    interpro_descs = if (nrow(ipr)) ipr$interpro_descs else NA_character_,
    go_terms       = if (nrow(ipr)) ipr$go_terms else NA_character_,
    combined_pdf_path = s$combined_pdf_path, heatmap_dir = s$heatmap_dir
  )
}
out <- rbindlist(rows, fill = TRUE)
setorder(out, serial_no)

out_path <- file.path(base04, "05_results", "tables", "publication_genes_table.csv")
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
fwrite(out, out_path)

cat(sprintf("Wrote %s (%d genes, %d columns)\n", out_path, nrow(out), ncol(out)))
print(out[, .(serial_no, gene_symbol, trait, gene_id=sub("HORVU.MOREX.r3.","",gene_id),
              class, MGmin, epsilon, n_groups, kw_p_raw=signif(kw_p_raw,3),
              FDR=signif(trait_fdr_p,3), Bonf=signif(trait_bonferroni_p,3),
              sigF=significant_fdr, sigB=significant_bonferroni)])
