#!/bin/bash
# ====================================================================
# GCTA Heritability Analysis Script - V3
# ====================================================================

set -e

echo "--- Initializing GCTA Heritability Analysis for V3 ---"

# --- 1. PATH DEFINITIONS ---
SCRIPT_DIR=$(dirname $(readlink -f "$0"))
BASE_DIR=$(dirname "$SCRIPT_DIR")
TOOLS_DIR="${BASE_DIR}/tools/gcta-1.95.0-linux-kernel-3-x86_64"
GCTA_INPUTS_DIR="${BASE_DIR}/inputs"

# --- 2. ANALYSIS-SPECIFIC DIRECTORIES ---
V3_DIR="${BASE_DIR}/analysis_V3"
V3_INTERMEDIATES="${V3_DIR}/intermediates"
V3_RESULTS="${V3_DIR}/results"
mkdir -p $V3_INTERMEDIATES
mkdir -p $V3_RESULTS
echo "Project directories are set up in: ${V3_DIR}"

# --- 3. INPUT FILE DEFINITIONS (V3 Specific) ---
PLINK_PREFIX_IN="${V3_INTERMEDIATES}/barley_V3_filtered"

# --- 4. PHENOTYPE DEFINITIONS ---
PHENOTYPES=( "Protein" "Starch" "Betaglucans" "Fiber" )

echo "--- Setup complete. Starting main analysis steps. ---"
echo "======================================================="

# ====================================================================
# --- STAGE 1: Create Genetic Relationship Matrix (GRM) ---
# ====================================================================
echo "--- Stage 1: Creating the Genetic Relationship Matrix (GRM) ---"
GRM_PREFIX_OUT="${V3_INTERMEDIATES}/barley_V3_grm"
if [ ! -f "${GRM_PREFIX_OUT}.grm.bin" ]; then
    ${TOOLS_DIR}/gcta64 --bfile $PLINK_PREFIX_IN \
        --autosome-num 7 \
        --make-grm \
        --out $GRM_PREFIX_OUT
    echo "GRM created successfully."
else
    echo "Found existing GRM file. Skipping creation."
fi
echo "--- Stage 1 Complete ---"
echo "======================================================="

# ====================================================================
# --- STAGE 2: Run Heritability Analysis (REML) for Each Trait ---
# ====================================================================
echo "--- Stage 2: Starting REML analysis for all traits ---"
for TRAIT in "${PHENOTYPES[@]}"; do
    echo ""
    echo "--- Analyzing trait: ${TRAIT} ---"
    PHENOTYPE_FILE="${GCTA_INPUTS_DIR}/${TRAIT}_for_gcta.txt"
    RESULT_PREFIX="${V3_RESULTS}/result_${TRAIT}_h2"
    if [ ! -f "${RESULT_PREFIX}.hsq" ]; then
        if [ -f "$PHENOTYPE_FILE" ]; then
            ${TOOLS_DIR}/gcta64 --grm $GRM_PREFIX_OUT \
                --pheno $PHENOTYPE_FILE \
                --reml \
                --out $RESULT_PREFIX
            echo "REML analysis for ${TRAIT} completed."
        else
            echo "WARNING: Phenotype file not found for ${TRAIT}. Skipping."
        fi
    else
        echo "Result file for ${TRAIT} already exists. Skipping analysis."
    fi
done
echo ""
echo "======================================================="
echo "### All GCTA analyses completed successfully! ###"
echo "### Results are located in: ${V3_RESULTS} ###"
echo "======================================================="
