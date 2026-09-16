# ============================================================================
# config/params.sh -- single source of truth for step 03.
#
# Sourced by every script (bash sources it directly; R reads it via
# scripts/_load_params.R). Change a value HERE and nowhere else, then re-run
# scripts/run_all.sh.
# ============================================================================

# ---- Paths -----------------------------------------------------------------
export STEP03_BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/03_USED_candidate_genes_around_leading_snps
export PROJECT_ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project

# The locus table handed over by step 01. THIS IS THE ONLY SNP/LOCUS INPUT.
export LOCI_HANDOFF_TSV="$PROJECT_ROOT/01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_for_gene_search.tsv"
# Step-01 parameter table, copied into results/ so the methods section is self-contained.
export GWAS_PARAMS_TSV="$PROJECT_ROOT/01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_analysis_parameters.tsv"

# Morex V3 annotation.
export GFF_GZ=/mnt/data/shahar/gwas_barley/morexV3_analysis/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3.gz

# ---- Tools (absolute; bare names resolve to broken versions on this host) ---
export BEDTOOLS=/usr/bin/bedtools

# ---- Search rule -----------------------------------------------------------
# FLANK_BP: bp added to EACH side of the step-01 LD locus span.
#
#   0 = the search interval IS the LD locus span, nothing added.
#
# Chosen 2026-09-08, unchanged 2026-09-09. Rationale: step 01 already defines each
# locus as a real LD block -- PLINK --clump at r2=0.5, then severed at the first
# internal gap larger than the gap rule. The block is therefore the unit in LD with
# the lead SNP and is the honest search space; adding an arbitrary flank on top would
# re-import the distance-only assumption that the LD clumping was introduced to replace.
#
# Consequence, accepted deliberately: 18 of the 37 loci contain no annotated gene
# and are reported in results/tables/loci_without_genes.tsv rather than silently
# dropped. scripts/04_flank_sensitivity.R quantifies what a flank would add.
#
# NOTE: step 01's clumping parameters are NOT set here and have been revised twice.
# History, so a stale result set can be recognised:
#   2026-09-08  clump_kb 1000, max_span 2000 kb, gap 50 kb, +2 sub-threshold protein
#               peaks   -> 34 loci,  9.40 Mb, 64 genes
#   2026-09-09a clump_kb 2000, max_span 4000 kb, gap 60 kb, +2 peaks
#                       -> 37 loci, 14.72 Mb, 90 genes
#   2026-09-09b clump_kb 2000, max_span 4000 kb, gap 50 kb, iterative, NO extra peaks
#                       -> 36 loci,  7.91 Mb, 55 genes   <- CURRENT
# This step reads whatever step 01 hands over and stores a stamped copy of it in
# inputs/loci_handoff_snapshot.tsv; the authoritative parameter list is mirrored
# into results/tables/analysis_parameters.tsv on every run.

export FLANK_BP=0

# Flank values (bp) reported side by side by 04_flank_sensitivity.R.
export FLANK_SENSITIVITY_SET="0 25000 50000 100000 200000"

# Include the 2 sub-threshold protein peaks that step 01 kept separate from the
# 32 significant loci? They carry include_in_gene_search=TRUE in the handoff
# table. "yes" honours that flag; "no" keeps only class=significant_locus.
export INCLUDE_SUBTHRESHOLD_PEAKS=yes

# ---- Environment -----------------------------------------------------------
export TMPDIR=/mnt/data/shahar/.tmp
