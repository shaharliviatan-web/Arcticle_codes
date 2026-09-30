#!/usr/bin/env bash
# config.sh - every parameter for the 7H shared-signal branch, single source of truth.
# Created 2026-09-17. Exploratory branch: nothing here feeds the main pipeline.
export TMPDIR=/mnt/data/shahar/.tmp
export PROJECT_ROOT=/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project
export BRANCH="$PROJECT_ROOT/03_01_7H_branch_Starch_Fiber_shared_signal_explore/archive_v1_MGmin2_eps0.6"   # archived 2026-09-24 (was the branch root)

# --- the gene under investigation -------------------------------------------
export GENE_ID="HORVU.MOREX.r3.7HG0729030"
export GENE_CHR="7H"
export GENE_START=573604051          # Ensembl r62 GFF3, gene feature
export GENE_END=573605728
export GENE_STRAND="-"               # minus strand => TSS is at GENE_END
export GENE_TSS=573605728

# --- the two GWAS signals (step 01, run 00_FINAL_BLUP_3PC) -------------------
export FIBER_LEAD=573606306          # fiber_L17, -log10p 7.2451
export STARCH_LEAD=573606460         # starch_L05, -log10p 6.4367
export BONF_THRESHOLD=6.0454         # Bonferroni alpha=0.10 over 111,017 pruned SNPs

# --- window / haplotyping ----------------------------------------------------
export WINDOW_BP=1000                # same rule as step 04 (gene +/- 1000 bp)
export MGMIN=2
export EPS_FIXED=0.6                 # the pipeline's locked value
export EPS_EXPLORE="0.4 0.6 0.8 1.0 1.2 1.5 2.0"

# --- external data -----------------------------------------------------------
export SOURCE_VCF=/mnt/data/shahar/gwas_barley/data/inputs/morexV3_with_ids.vcf.gz
export STAGED_IMPUTED="$PROJECT_ROOT/04_USED_haplotype_analysis_crosshap/03_per_gene_vcfs/staging/morexV3_final_imputed.vcf.gz"
export SAMPLE_RENAME="$PROJECT_ROOT/04_USED_haplotype_analysis_crosshap/03_per_gene_vcfs/staging/sample_rename_doubled_to_single.txt"
export KEEP_290="$PROJECT_ROOT/04_USED_haplotype_analysis_crosshap/00_config/samples_keep_290.txt"
export PHENO_ROOT=/mnt/data/shahar/gwas_barley/data/inputs
export PHENO_SUFFIX="_corrected_V3.pheno"
export ASSOC_DIR="$PROJECT_ROOT/01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/01_assoc"
export FREQ_FRQ="$PROJECT_ROOT/01_USED_GWAS_V2_pipeline/intermediates/morexV3_290_freq.frq"
export LOCI_TABLE="$PROJECT_ROOT/01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_for_gene_search.tsv"
export PGSB_GFF=/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.gff3
export PGSB_PROT=/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.aa.fa
export ENSEMBL_GFF=/mnt/data/shahar/gwas_barley/morexV3_analysis/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3.gz
export SWISSPROT_DMND="$PROJECT_ROOT/05_USED_gene_annotation_analysis/05_USED_gene_annotation/intermediates/swissprot.dmnd"

# --- binaries ----------------------------------------------------------------
export BCFTOOLS=/usr/bin/bcftools
export PLINK=/usr/local/bin/plink     # bare `plink` is a broken v0.76 - use the absolute path
export DIAMOND=/usr/bin/diamond
export BEDTOOLS=/usr/bin/bedtools
