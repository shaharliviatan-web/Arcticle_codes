#!/bin/bash
set -e

# ====================================================================
# MASTER PIPELINE SCRIPT
# This script runs the entire Imputation workflow from start to finish.
# It calls the R script for the core logic and then performs
# Bash-based data reshaping to create a valid, PLINK-compatible VCF.
#
# Run this script from the 'imputation' directory.
# ====================================================================

# --- Define File Names ---
ORIGINAL_VCF="/mnt/data/shahar/gwas_barley/data/inputs/morexV3.vcf.gz"
SAMPLES_TO_REMOVE="/mnt/data/shahar/gwas_barley/data/inputs/samples_to_remove_V3.txt"
CLEAN_VCF="morexV3_samples_removed.vcf"
IMPUTED_LFMM="morexV3_samples_removed.lfmm_imputed.lfmm"
TRANSPOSED_FILE="imputed_transposed.txt"
GENOTYPES_FILE="imputed_genotypes_final.txt"
SKELETON_FILE="skeleton_data_only.txt"
FINAL_VCF_OUT="morexV3_final_imputed.vcf"


echo "--- STEP 1: Filtering VCF (PLINK) ---"
plink --vcf $ORIGINAL_VCF \
      --remove $SAMPLES_TO_REMOVE \
      --recode vcf \
      --allow-extra-chr \
      --out morexV3_samples_removed

echo "--- STEP 2: Running Core Imputation (R) ---"
# This will run the R script which handles conversion, sNMF, and imputation
# This step will take many hours.
Rscript imputation_core.R

echo "--- STEP 3: Transposing Imputed Data (AWK) ---"
# Transpose the (Samples x SNPs) file to (SNPs x Samples)
# Use TAB ('\t') as the delimiter. This is a very memory-heavy step.
awk '{for (i=1; i<=NF; i++) {a[i,NR] = $i}} END {for (i=1; i<=NF; i++) {for (j=1; j<=NR; j++) {printf "%s%s", a[i,j], (j==NR ? "" : "\t")} printf "\n"}}' $IMPUTED_LFMM > $TRANSPOSED_FILE

echo "--- STEP 4: Converting Genotypes (AWK) ---"
# Convert numeric (0, 2) format to VCF format (0/0, 1/1)
# Force the output field separator (OFS) to be a TAB.
awk 'BEGIN{OFS="\t"} {for(i=1;i<=NF;i++) {if($i==0) $i="0/0"; else if($i==2) $i="1/1"} print}' $TRANSPOSED_FILE > $GENOTYPES_FILE

echo "--- STEP 5: Assembling Final VCF (GREP/CUT/PASTE) ---"
# 1. Add the '##' header lines
grep "^##" $CLEAN_VCF > $FINAL_VCF_OUT

# 2. Add the '#CHROM' data header line
grep "^#CHROM" $CLEAN_VCF >> $FINAL_VCF_OUT

# 3. Create a temporary "skeleton" file (just the 9 info columns)
grep -v "^#" $CLEAN_VCF | cut -f 1-9 > $SKELETON_FILE

# 4. Paste the skeleton and the new genotypes together
paste $SKELETON_FILE $GENOTYPES_FILE >> $FINAL_VCF_OUT

echo "--- STEP 6: Cleaning up temporary files ---"
rm $TRANSPOSED_FILE
rm $GENOTYPES_FILE
rm $SKELETON_FILE

echo "--- PIPELINE COMPLETE! ---"
echo "Final imputed VCF is ready at: $FINAL_VCF_OUT"