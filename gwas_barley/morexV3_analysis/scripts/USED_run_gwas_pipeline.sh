#!/bin/bash

# ====================================================================

# GWAS Pipeline Script - v1.0 - DEDICATED to morexV3 ANALYSIS

#

# This script is refactored for efficiency based on the V2 analysis.

# Phenotype-independent steps are performed once, before the loop.

# It is fully parameterized, including the number of PCs.

# ====================================================================



# --- Step 0: Script Setup ---

set -e



# --- Step 1: CONTROL PANEL & Main Variables ---

echo "--- Initializing GWAS Batch Run for morexV3 ---"



# --- Path Definitions ---

SCRIPT_DIR=$(dirname $(readlink -f "$0"))

ANALYSIS_DIR=$(dirname "$SCRIPT_DIR")

PROJECT_DIR=$(dirname "$ANALYSIS_DIR")

INPUT_DIR="${PROJECT_DIR}/data/inputs"



# --- Define the static list of phenotypes for V3 analysis ---

PHENOTYPES_TO_RUN=(

    "protein_corrected_V3"

    "starch_corrected_V3"

    "fiber_corrected_V3"

    "betaglucan_corrected_V3"

)

# ====================================================================



# ====================================================================

# --- CONTROL PANEL (Interactive) ---

# ====================================================================

read -p "Enter MAF value [default: 0.05]: " MAF

MAF=${MAF:-0.05}



read -p "Enter GENO value (max missingness) [default: 0.1]: " GENO

GENO=${GENO:-0.1}



read -p "Enter thinning distance in bp (use 0 for none) [default: 0]: " THIN

THIN=${THIN:-0}



read -p "Enter number of PCs to calculate [default: 20]: " PCS

PCS=${PCS:-20}

# ====================================================================



# --- Static & Dynamic Paths ---

TOOLS_DIR="${PROJECT_DIR}/tools"

SAMPLES_TO_REMOVE="${INPUT_DIR}/samples_to_remove_V3.txt" # Adapted for V3

INTERMEDIATES_DIR="${ANALYSIS_DIR}/data/intermediates_maf${MAF}_geno${GENO}_thin${THIN}_pc${PCS}" # Now includes PCs

RESULTS_DIR="${ANALYSIS_DIR}/data/results_maf${MAF}_geno${GENO}_thin${THIN}_pc${PCS}" # Now includes PCs

VCF_NICKNAME="morexV3" # Adapted for V3

VCF_IN="${INPUT_DIR}/morexV3.vcf.gz" # Adapted for V3



# --- All Major File Prefixes ---

PLINK_BASE_PREFIX="${INTERMEDIATES_DIR}/${VCF_NICKNAME}_barley_full_filtered"

PCA_PREFIX="${RESULTS_DIR}/${VCF_NICKNAME}_pca"

KINSHIP_FILE="${RESULTS_DIR}/${VCF_NICKNAME}_kinship_matrix.aIBS.kinf"

GWAS_TPED_PREFIX="${RESULTS_DIR}/${VCF_NICKNAME}_gwas_genotypes"



echo ""

echo "--- Strategy Summary ---"

echo "Target VCF: ${VCF_NICKNAME}"

echo "Phenotypes: ${PHENOTYPES_TO_RUN[@]}"

echo "Parameters: MAF=${MAF}, GENO=${GENO}, THIN=${THIN}bp, PCs=${PCS}"

echo "Output will be in: ${RESULTS_DIR}"

echo "---------------------------------"



mkdir -p $INTERMEDIATES_DIR $RESULTS_DIR



# ====================================================================

# --- PRE-LOOP: PHENOTYPE-INDEPENDENT STEPS ---

# ====================================================================



# --- Step 2.1: VCF to PLINK conversion ---

if [ ! -f "${PLINK_BASE_PREFIX}.bed" ]; then

    echo "--- Step 2.1: Creating filtered PLINK binary file ---"

    plink --vcf $VCF_IN --maf ${MAF} --geno ${GENO} --allow-extra-chr --make-bed --out $PLINK_BASE_PREFIX

else

    echo "--- Step 2.1: Found existing PLINK file. Skipping. ---"

fi



# --- Step 2.2: Prepare files for PCA and Kinship ---

if [ ! -f "${PCA_PREFIX}.eigenvec" ] || [ ! -f "$KINSHIP_FILE" ]; then

    echo "--- Step 2.2: Preparing Covariate Files (PCA & Kinship) ---"

    PRUNED_PREFIX="${INTERMEDIATES_DIR}/${VCF_NICKNAME}_pruned_for_covs"

    plink --bfile $PLINK_BASE_PREFIX --remove $SAMPLES_TO_REMOVE --indep-pairwise 50 5 0.2 --allow-extra-chr --out $PRUNED_PREFIX

   

    PRUNED_COV_PREFIX="${INTERMEDIATES_DIR}/${VCF_NICKNAME}_pruned_covariates_data"

    plink --bfile $PLINK_BASE_PREFIX --remove $SAMPLES_TO_REMOVE --extract "${PRUNED_PREFIX}.prune.in" --recode12 transpose --allow-extra-chr --out $PRUNED_COV_PREFIX

   

    echo "Calculating PCA..."

    plink --tfile $PRUNED_COV_PREFIX --pca ${PCS} --allow-extra-chr --out $PCA_PREFIX

   

    echo "Calculating Kinship Matrix..."

    ${TOOLS_DIR}/emmax-kin-intel64 -v -s -d 10 -x $PRUNED_COV_PREFIX

    mv "${PRUNED_COV_PREFIX}.aIBS.kinf" "$KINSHIP_FILE"

else

    echo "--- Step 2.2: Found existing PCA/Kinship files. Skipping. ---"

fi



# --- Step 2.3: Prepare main genotype file for GWAS ---

if [ ! -f "${GWAS_TPED_PREFIX}.tped" ]; then

    echo "--- Step 2.3: Preparing Main Genotype File ---"

   

    FINAL_PLINK_PREFIX=$PLINK_BASE_PREFIX

   

    if [ "$THIN" -ne 0 ]; then

        echo "Thinning strategy detected (thin = ${THIN}bp)."

        VCF_FILTERED_PREFIX="${INTERMEDIATES_DIR}/${VCF_NICKNAME}_filtered_samples"

        plink --bfile $PLINK_BASE_PREFIX --remove $SAMPLES_TO_REMOVE --recode vcf --allow-extra-chr --out $VCF_FILTERED_PREFIX

       

        THINNED_VCF_OUT="${INTERMEDIATES_DIR}/${VCF_NICKNAME}_thinned_${THIN}bp.vcf.gz"

        vcftools --vcf ${VCF_FILTERED_PREFIX}.vcf --thin $THIN --recode --recode-INFO-all --stdout | gzip > $THINNED_VCF_OUT

       

        THINNED_PLINK_PREFIX_TEMP="${INTERMEDIATES_DIR}/${VCF_NICKNAME}_final_for_gwas"

        plink --vcf $THINNED_VCF_OUT --allow-extra-chr --make-bed --out $THINNED_PLINK_PREFIX_TEMP

        FINAL_PLINK_PREFIX=$THINNED_PLINK_PREFIX_TEMP

    fi



    # Create the final TPED file for EMMAX, consistently removing samples

    plink --bfile $FINAL_PLINK_PREFIX --remove $SAMPLES_TO_REMOVE --recode12 transpose --allow-extra-chr --out $GWAS_TPED_PREFIX

   

    echo "Copying final SNP map..."

    cp "${FINAL_PLINK_PREFIX}.bim" "${RESULTS_DIR}/snp_map.bim"

else

    echo "--- Step 2.3: Found existing GWAS TPED file. Skipping. ---"

    if [ ! -f "${RESULTS_DIR}/snp_map.bim" ]; then

      echo "WARNING: GWAS TPED exists but snp_map.bim is missing! You may need to delete intermediates to regenerate."

    fi

fi

echo "---------------------------------"



# ====================================================================

# --- PHENOTYPE LOOP: RUN GWAS FOR EACH PHENOTYPE ---

# ====================================================================

echo "--- Starting GWAS runs for each phenotype ---"

for PHENOTYPE_NAME in "${PHENOTYPES_TO_RUN[@]}"; do

   

    PHENOTYPE_FILE="${INPUT_DIR}/${PHENOTYPE_NAME}.pheno"

    EMMAX_OUT_PREFIX="${RESULTS_DIR}/${VCF_NICKNAME}_${PHENOTYPE_NAME}_gwas"



    if [ ! -f "$PHENOTYPE_FILE" ]; then

        echo "WARNING: Phenotype file for ${PHENOTYPE_NAME} not found. Skipping."

        continue

    fi



    if [ ! -f "${EMMAX_OUT_PREFIX}.ps" ]; then

        echo "--- Running EMMAX for ${PHENOTYPE_NAME} ---"

        ${TOOLS_DIR}/emmax-intel64 -v -d 10 \

            -t $GWAS_TPED_PREFIX \

            -p $PHENOTYPE_FILE \

            -k $KINSHIP_FILE \

            -c "${PCA_PREFIX}.eigenvec" \

            -o $EMMAX_OUT_PREFIX

    else

        echo "--- GWAS result for ${PHENOTYPE_NAME} already exists. Skipping. ---"

    fi

    echo "---------------------------------"

done



echo "#####################################################################"

echo "### All V3 GWAS batch runs completed successfully! ###"

echo "#####################################################################"