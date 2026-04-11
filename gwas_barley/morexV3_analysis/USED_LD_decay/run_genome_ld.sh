#!/bin/bash
# =============================================================
# run_genome_ld.sh
#
# This script runs the full-genome PLINK LD decay calculation.
# v1.1 - Changed --ld-window-bp to --ld-window-kb
#        to support older PLINK 1.9 builds.
# =============================================================

# Stop the script if any command fails
set -e

echo "--- Script started in directory: $(pwd) ---"
echo "--- This will take SEVERAL HOURS ---"

# --- 1. Define file paths ---
PLINK_PREFIX="../data/intermediates_maf0.000001_geno1_thin0_pc5/plink_with_chr_pos_ids/morexV3_full_filtered_chr_pos_ids"
SAMPLES_TO_REMOVE="../../data/inputs/samples_to_remove_V3.txt"
OUTPUT_FILE_PREFIX="genome_ld_decay_data"

# --- 2. Run PLINK for Full Genome LD Decay (Corrected) ---
plink --bfile "$PLINK_PREFIX" \
      --remove "$SAMPLES_TO_REMOVE" \
      --allow-extra-chr \
      --chr 1H,2H,3H,4H,5H,6H,7H \
      --r2 gz \
      --ld-window-kb 5000 \
      --ld-window 999999 \
      --ld-window-r2 0.01 \
      --out "$OUTPUT_FILE_PREFIX"

echo "--- 🥳 PLINK calculation finished! ---"
echo "Raw LD data is in $(pwd)/${OUTPUT_FILE_PREFIX}.ld.gz"