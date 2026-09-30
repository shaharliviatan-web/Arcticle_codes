#!/usr/bin/env Rscript
# 07_distances_current_loci.R
# Distance from each canonical (1,3;1,4)-beta-glucan gene to the nearest GWAS lead SNP,
# measured against the CURRENT 36-locus set.
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# `06_build_table.R` answered the same question in June 2026, but against the retired
# candidate-gene pipeline: a +/-200 kb window around each lead and a 20-locus set. Both are
# gone (step 01 now defines loci by iterative LD clumping -> 36 loci; step 03 searches the
# LD span with no flank), so every distance in
# `results/tables/canonical_betaglucan_gene_distances.tsv` is STALE and must not be quoted.
# Some of the lead SNPs in that table are no longer significant at all.
#
# What is NOT stale is the expensive part: the sequence-based mapping of the 10 canonical
# genes onto Morex V3 (steps 01-05). That mapping is assembly-independent and unchanged, so
# this script reuses the gene coordinates from the old table and recomputes only the
# distances, against the live locus table. It does not re-run any BLAST.
#
# METHOD
#   * distance is measured from the nearest edge of the gene body to the lead SNP, and is 0
#     if the lead falls inside the gene;
#   * only leads on the SAME chromosome are considered -- a cross-chromosome "distance" is
#     meaningless. A gene on a chromosome that carries no locus of that trait gets NA, which
#     is itself the answer (e.g. HvCslF6 is on 7H and no beta-glucan locus is on 7H);
#   * one row per (canonical gene x trait), as in the old table, so the two are comparable.
#
# Inputs (read-only):
#   results/tables/canonical_betaglucan_gene_distances.tsv   gene -> Morex V3 coordinates
#   01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv
# Output:
#   results/tables/canonical_gene_distances_current_loci.tsv
#
# Created 2026-09-23 for Results ch. 3 (the "no canonical beta-glucan gene lies near any
# beta-glucan lead" statement). Run: Rscript scripts/07_distances_current_loci.R

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))

PROJ <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
B09  <- file.path(PROJ, "05_USED_gene_annotation_analysis", "09_USED_canonical_betaglucan_gene_check")
LOCI <- file.path(PROJ, "01_USED_GWAS_V2_pipeline", "results", "00_FINAL_BLUP_3PC",
                  "02_loci_FINAL", "tables", "Table_loci_master.tsv")
OLD  <- file.path(B09, "results", "tables", "canonical_betaglucan_gene_distances.tsv")
OUT  <- file.path(B09, "results", "tables", "canonical_gene_distances_current_loci.tsv")

G <- unique(fread(OLD)[, .(canonical_gene, source_accession, horvu_id, chr, gene_start,
                           gene_end, strand, pct_identity, pct_coverage, reciprocal_best_hit)])
L <- fread(LOCI)
stopifnot(nrow(G) == 10L, nrow(L) == 36L,
          all(c("trait", "chr", "lead_SNP", "lead_bp", "lead_neg_log10_p") %in% names(L)))

TRAITS <- c("betaglucan", "fiber", "starch")   # protein is not a cell-wall trait; kept out, as in step 06

res <- rbindlist(lapply(TRAITS, function(tr) {
  G[, {
    s <- L[trait == tr & chr == .BY$chr]
    if (!nrow(s)) {
      .(trait = tr, nearest_lead_SNP = NA_character_, nearest_locus_id = NA_character_,
        lead_neg_log10p = NA_real_, distance_bp = NA_real_,
        note = "no locus of this trait on this chromosome")
    } else {
      d <- pmax(0, pmin(abs(s$lead_bp - gene_start), abs(s$lead_bp - gene_end)))
      d[s$lead_bp >= gene_start & s$lead_bp <= gene_end] <- 0
      i <- which.min(d)
      .(trait = tr, nearest_lead_SNP = s$lead_SNP[i], nearest_locus_id = s$locus_id[i],
        lead_neg_log10p = s$lead_neg_log10_p[i], distance_bp = d[i], note = "")
    }
  }, by = .(canonical_gene, horvu_id, chr, gene_start, gene_end)]
}))
res[, distance_Mb := round(distance_bp / 1e6, 1)]
setorder(res, trait, distance_bp, na.last = TRUE)
fwrite(res, OUT, sep = "\t")

bg <- res[trait == "betaglucan"]
cat(sprintf("[09-07] %d canonical genes x %d traits -> %s\n", nrow(G), length(TRAITS), OUT))
cat(sprintf("[09-07] beta-glucan: nearest canonical gene is %s at %.1f Mb from %s (%s)\n",
            bg[which.min(distance_bp), canonical_gene], bg[which.min(distance_bp), distance_Mb],
            bg[which.min(distance_bp), nearest_lead_SNP], bg[which.min(distance_bp), nearest_locus_id]))
cat(sprintf("[09-07] on a chromosome with no beta-glucan locus: %s\n",
            paste(bg[is.na(distance_bp), paste0(canonical_gene, " (", chr, ")")], collapse = ", ")))
cat("[09-07] every canonical gene, beta-glucan:\n")
print(bg[, .(canonical_gene, chr, nearest_lead_SNP, distance_Mb)])
