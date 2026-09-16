#!/usr/bin/env Rscript
# =============================================================================
# 01_build_paper_table.R
#
# WHAT   One publication-ready table for the significant haplotype genes: identity,
#        position, GWAS origin, haplotype statistics, effect sizes and the functional
#        call with its evidence, ordered by significance. This is the table the
#        manuscript is written from.
#
# READS  results/tables/fdr_annotation_master.tsv           (step 08 merge)
#        ../../04_.../<run>/Significant_genes/significant_genes.tsv  (figure paths)
# WRITES results/tables/Table_significant_genes_annotated.tsv   full, all evidence
#        results/tables/Table_significant_genes_paper.tsv       trimmed, paper columns
#
# serial_no is assigned HERE, by ascending BH q, and is the stable public identifier
# for a gene in the manuscript and in the figure filenames.
#
# trait_candidate_strength is written EMPTY on purpose: it is a manual judgement made
# on the reported evidence (literature support for a link to the trait), not something
# the pipeline can decide. Fill it in review; re-running will not overwrite a value
# already present in results/tables/trait_candidate_calls.tsv if that file exists.
# =============================================================================
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))

ROOT <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
ANNO <- file.path(ROOT, "05_USED_gene_annotation_analysis")
OUT  <- file.path(ANNO, "08_USED_annotation_master/results/tables")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

m <- fread(file.path(OUT, "fdr_annotation_master.tsv"))

# figure paths + group sizes from the step-04 collection
run <- Sys.getenv("STEP04_RUN", "loci_LDspan_eps06_V4")
sg_path <- file.path(ROOT, "04_USED_haplotype_analysis_crosshap", "04_runs", run,
                     "Significant_genes", "significant_genes.tsv")
if (file.exists(sg_path)) {
  sg <- fread(sg_path)
  keep <- c("gene_id","trait","chr","gene_start","gene_end","strand","group_sizes",
            "n_ind","n_unassigned","top_group","top_n","top_mean",
            "bottom_group","bottom_n","bottom_mean","figure_stem")
  m <- merge(m, sg[, ..keep], by = c("gene_id","trait"), all.x = TRUE)
}

setorder(m, fdr_q, kw_p_raw)
m[, serial_no := .I]
m[, gene_short := sub("^HORVU\\.MOREX\\.r3\\.", "", gene_id)]

# manual trait-candidacy calls, if a review file exists
cc <- file.path(OUT, "trait_candidate_calls.tsv")
if (file.exists(cc)) {
  cur <- fread(cc)
  m <- merge(m, cur[, .(gene_id, trait, trait_candidate_strength, candidate_rationale)],
             by = c("gene_id","trait"), all.x = TRUE)
  setorder(m, serial_no)
} else {
  m[, trait_candidate_strength := NA_character_]
  m[, candidate_rationale := NA_character_]
}

full_cols <- c("serial_no","trait","gene_short","gene_id","chr","gene_start","gene_end","strand",
               "locus_id","lead_SNP","lead_SNP_class","dist_to_lead_bp",
               "n_haplotype_groups","group_sizes","n_ind","n_unassigned",
               "kw_p_raw","fdr_q","eta_squared","delta_top_bottom_sd",
               "top_group","top_n","top_mean","bottom_group","bottom_n","bottom_mean",
               "final_call","annotation_source","annotation_from_check",
               "n_sources_with_call","needs_review_two_calls",
               "call_check1_swissprot","call_check2_interpro","call_check3_nr",
               "sp_vs_interpro_agreement",
               "trait_candidate_strength","candidate_rationale",
               "sp_name","sp_accession","sp_organism","sp_pident","sp_qcovhsp","sp_evalue","sp_status",
               "interpro_status","n_interpro","interpro_ids","interpro_descs",
               "n_pfam","pfam_ids","pfam_descs","n_go","go_terms",
               "nr_title","nr_accession","nr_pident","nr_qcovs","nr_evalue","nr_status",
               "figure_stem")
full_cols <- intersect(full_cols, names(m))
fwrite(m[, ..full_cols], file.path(OUT, "Table_significant_genes_annotated.tsv"), sep = "\t", na = "NA")

paper_cols <- intersect(c("serial_no","trait","gene_short","gene_id","chr","gene_start","gene_end",
                          "locus_id","lead_SNP","dist_to_lead_bp",
                          "n_haplotype_groups","group_sizes","n_ind",
                          "kw_p_raw","fdr_q","eta_squared","delta_top_bottom_sd",
                          "final_call","annotation_source","annotation_from_check",
                          "n_sources_with_call","needs_review_two_calls",
                          "call_check1_swissprot","call_check2_interpro","call_check3_nr",
                          "trait_candidate_strength","candidate_rationale"), names(m))
fwrite(m[, ..paper_cols], file.path(OUT, "Table_significant_genes_paper.tsv"), sep = "\t", na = "NA")

cat(sprintf("Wrote Table_significant_genes_annotated.tsv (%d rows x %d cols)\n", nrow(m), length(full_cols)))
cat(sprintf("Wrote Table_significant_genes_paper.tsv     (%d rows x %d cols)\n", nrow(m), length(paper_cols)))
cat("\nannotation_source:\n"); print(table(m$annotation_source))
cat("\nper trait:\n"); print(m[, .N, by = trait])
if (all(is.na(m$trait_candidate_strength)))
  cat("\ntrait_candidate_strength is empty - fill it by review, see the script header.\n")
