#!/usr/bin/env bash
# config.sh - every parameter of the REPLACEMENT analysis. Single source of truth.
export TMPDIR=/mnt/data/shahar/.tmp
export PROJECT_ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project
export BRANCH="$PROJECT_ROOT/03_01_7H_branch_Starch_Fiber_shared_signal_explore"
export TEMP_ROOT="$BRANCH"     # 2026-09-24: this analysis became the branch root (was $BRANCH/REPLACEMENT_ANALYSIS_MGmin3_eps0.9)
export ARCHIVE_V1="$BRANCH/archive_v1_MGmin2_eps0.6"   # the previous version (MGmin 2, eps 0.6), archived
export GENE_ID="HORVU.MOREX.r3.7HG0729030"; export SHORT_NAME="GDSL"
export GENE_CHR="7H"; export GENE_START=573604051; export GENE_END=573605728
export GENE_STRAND="-"; export GENE_TSS=573605728
export WIN_START=573603051; export WIN_END=573606728; export WINDOW_BP=1000
export FIBER_LEAD=573606306; export STARCH_LEAD=573606460
export SIGNAL_SNPS="573606282 573606306 573606460 573606491"
# --- THE TWO PARAMETERS THAT DIFFER FROM THE PIPELINE ------------------------
export MGMIN=3          # pipeline uses 2; at 2 this gene has NO usable grouping
export EPSILON=0.9      # pipeline uses 0.6
export MINHAP=9         # unchanged
# -----------------------------------------------------------------------------
export STEP04_ROOT="$PROJECT_ROOT/04_USED_haplotype_analysis_crosshap"
export STEP07_ROOT="$PROJECT_ROOT/07_USED_elite_lines_compariosn_to_wild_lines"
export PHENO_ROOT=/mnt/data/shahar/gwas_barley/data/inputs
export PHENO_SUFFIX="_corrected_V3.pheno"
export SOURCE_VCF=/mnt/data/shahar/gwas_barley/data/inputs/morexV3_with_ids.vcf.gz
export STAGED_IMPUTED="$STEP04_ROOT/03_per_gene_vcfs/staging/morexV3_final_imputed.vcf.gz"
export SAMPLE_RENAME="$STEP04_ROOT/03_per_gene_vcfs/staging/sample_rename_doubled_to_single.txt"
export KEEP_290="$STEP04_ROOT/00_config/samples_keep_290.txt"
export DIVBROWSE_BASE="https://divbrowse.ipk-gatersleben.de/barley_pangenome_v2"
export BCFTOOLS=/usr/bin/bcftools; export BGZIP=/usr/bin/bgzip
export PLINK=/usr/local/bin/plink   # bare `plink` is a broken v0.76
