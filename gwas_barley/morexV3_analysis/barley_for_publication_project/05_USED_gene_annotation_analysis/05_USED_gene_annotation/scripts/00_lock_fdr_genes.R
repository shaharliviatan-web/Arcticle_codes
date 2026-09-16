#!/usr/bin/env Rscript
# 00_lock_fdr_genes.R
# Step 00 of 05_USED_gene_annotation.
#
# Lock the input set: the genes that passed BH-FDR in step 04.
#
# REPOINTED 2026-09-09 to the rebuilt step 04. The v1 version is in
#   ../_ARCHIVE_v1_45genes_2026-09-09/00_lock_fdr_genes.R
#
# WHAT CHANGED
#   old source : 04_.../05_results/tables/gene_annotation_review.csv  (v1, 48 rows / 45 genes)
#   new source : 04_.../04_runs/<run>/Significant_genes/significant_genes.tsv (21 genes)
#   Column renames forced by the step-04 rewrite:
#     trait_fdr_p            -> fdr_q
#     significant_bonferroni -> REMOVED (v1 reported BH and Bonferroni in parallel;
#                               step 04 now reports one correction only)
#     annotation             -> REMOVED (v1's non-reproducible manual "legacy_annotation";
#                               it was comparison-only and never the functional call)
#   Effect sizes (eta_squared, delta_top_bottom_sd) are new and are carried through,
#   because they, not the q-value, are what ranks these genes biologically.
#
#   The new source is ALREADY filtered to significant_fdr == TRUE, so no filtering
#   happens here; the row count is asserted against the source instead.
#
# Outputs:
#   inputs/fdr_genes.txt        unique gene_id, one per line
#   inputs/fdr_genes_table.tsv  one row per (trait, gene_id) with step-04 statistics

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

base_dir <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/05_USED_gene_annotation_analysis/05_USED_gene_annotation"
src_tsv  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Significant_genes/significant_genes.tsv"

inputs_dir <- file.path(base_dir, "inputs")
dir.create(inputs_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(src_tsv))
df <- read.delim(src_tsv, stringsAsFactors = FALSE, check.names = FALSE)

need <- c("trait","gene_id","locus_id","class","lead_SNP","dist_to_lead_bp",
          "kw_p_raw","fdr_q","eta_squared","delta_top_bottom_sd","n_groups")
missing <- setdiff(need, names(df))
if (length(missing)) stop("Source is missing columns: ", paste(missing, collapse=", "),
                          "\nHas step 04 changed again? See this script's header.")

out <- data.frame(
  trait               = df$trait,
  gene_id             = df$gene_id,
  locus_id            = df$locus_id,
  locus_class         = df$class,
  lead_SNP            = df$lead_SNP,
  dist_to_lead_bp     = df$dist_to_lead_bp,
  n_haplotype_groups  = df$n_groups,
  kw_p_raw            = df$kw_p_raw,
  fdr_q               = df$fdr_q,
  eta_squared         = df$eta_squared,
  delta_top_bottom_sd = df$delta_top_bottom_sd,
  stringsAsFactors    = FALSE)
out <- out[order(out$trait, out$gene_id), , drop = FALSE]

write.table(out, file.path(inputs_dir, "fdr_genes_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

genes <- sort(unique(out$gene_id))
writeLines(genes, file.path(inputs_dir, "fdr_genes.txt"))

cat(sprintf("Locked %d (trait, gene) rows = %d unique genes\n", nrow(out), length(genes)))
cat("Per trait:\n"); print(table(out$trait))
if (any(duplicated(out$gene_id)))
  cat("\nGenes appearing under >1 trait: ",
      paste(unique(out$gene_id[duplicated(out$gene_id)]), collapse=", "), "\n")
cat("\nWrote inputs/fdr_genes.txt and inputs/fdr_genes_table.tsv\n")
