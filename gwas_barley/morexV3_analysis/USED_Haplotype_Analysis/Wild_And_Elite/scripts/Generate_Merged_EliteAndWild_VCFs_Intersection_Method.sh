#!/bin/bash

# ==============================================================================
# SCRIPT: Merge_Wild_Elite_8Genes_Intersection.sh
# DESCRIPTION: Creates a STRICT INTERSECTION merge for the selected 8 genes.
#              Final VCF contains ONLY the exact SNPs present in the harmonized
#              Wild VCF. Elite VCF is sanitized, subsetted to 6 selected lines,
#              chromosome-renamed, sample-renamed, merged with -0, then filtered
#              to the Wild SNP list.
# ==============================================================================

set -euo pipefail

# ==============================================================================
# 1. DIRECTORY CONFIGURATION
# ==============================================================================
ROOT_DIR="/mnt/data/shahar/gwas_barley/morexV3_analysis"
WILD_ROOT="$ROOT_DIR/USED_haplotype_vcfs_per_gene"
ELITE_ROOT="$ROOT_DIR/USED_Haplotype_Analysis/Wild_And_Elite/Elite_VCFs"
SCRIPT_DIR="$ROOT_DIR/USED_Haplotype_Analysis/Wild_And_Elite/scripts"
OUT_DIR="$ROOT_DIR/USED_Haplotype_Analysis/Wild_And_Elite/Wild_And_Elite_VCFs__IntersectionMethod"

mkdir -p "$OUT_DIR"

REPORT_FILE="$OUT_DIR/merge_intersection_report.txt"
> "$REPORT_FILE"

log() {
    echo -e "[$(date '+%H:%M:%S')] $1" | tee -a "$REPORT_FILE"
}

log "=== STARTING 8-GENE STRICT INTERSECTION MERGE PIPELINE ==="

# ==============================================================================
# 2. ELITE LINES TO KEEP
# ==============================================================================
# SAMEA -> final sample name
SELECTED_LINES=(
    "SAMEA6619031:Corvette"
    "SAMEA6618923:Triumph"
    "SAMEA6619023:Iris"
    "SAMEA6619038:Gunnar"
    "SAMEA110189323:Manchuria_A"
    "SAMEA110189222:Rapid"
)

# ==============================================================================
# 3. TARGET GENES
# Format: "Trait:Gene:TargetChrom"
# ==============================================================================
TARGETS=(
    "betaglucan:HORVU.MOREX.r3.2HG0112910:2H"
    "betaglucan:HORVU.MOREX.r3.3HG0299450:3H"
    "betaglucan:HORVU.MOREX.r3.3HG0322580:3H"
    "fiber:HORVU.MOREX.r3.7HG0642350:7H"
    "fiber:HORVU.MOREX.r3.7HG0729010:7H"
    "starch:HORVU.MOREX.r3.3HG0301710:3H"
    "starch:HORVU.MOREX.r3.3HG0301750:3H"
    "starch:HORVU.MOREX.r3.3HG0301770:3H"
)

# ==============================================================================
# 4. MAIN LOOP
# ==============================================================================
for item in "${TARGETS[@]}"; do
    TRAIT=$(echo "$item" | cut -d':' -f1)
    GENE=$(echo "$item" | cut -d':' -f2)
    TARGET_CHROM=$(echo "$item" | cut -d':' -f3)
    RAW_CHROM="chr${TARGET_CHROM}"

    log "\n------------------------------------------------------------"
    log "Processing: $GENE ($TRAIT)"

    WILD_VCF="$WILD_ROOT/${TRAIT}_corrected_V3_targets_gene_220000bp_HARMONIZED/${GENE}.vcf.gz"
    ELITE_RAW="$ELITE_ROOT/${TRAIT}/${GENE}_custom_export.vcf"

    if [ ! -f "$WILD_VCF" ]; then
        log "   [ERROR] Missing wild VCF: $WILD_VCF"
        exit 1
    fi
    if [ ! -f "$ELITE_RAW" ]; then
        log "   [ERROR] Missing elite raw VCF: $ELITE_RAW"
        exit 1
    fi

    WORK_TMP="$OUT_DIR/tmp_${GENE}"
    mkdir -p "$WORK_TMP"

    SAMPLES_LIST="$WORK_TMP/samples_to_keep.txt"
    RENAME_MAP="$WORK_TMP/rename_map.txt"

    : > "$SAMPLES_LIST"
    : > "$RENAME_MAP"
    for pair in "${SELECTED_LINES[@]}"; do
        sample_id="${pair%%:*}"
        sample_name="${pair#*:}"
        echo "$sample_id" >> "$SAMPLES_LIST"
        echo "$sample_id $sample_name" >> "$RENAME_MAP"
    done

    # --- A. Verify Wild VCF & Extract Exact Wild Sites ---
    if [ ! -f "${WILD_VCF}.csi" ] && [ ! -f "${WILD_VCF}.tbi" ]; then
        log "   Indexing harmonized wild VCF..."
        bcftools index "$WILD_VCF"
    fi

    WILD_SITES="$WORK_TMP/wild_exact_sites.tsv"
    bcftools query -f '%CHROM\t%POS\n' "$WILD_VCF" > "$WILD_SITES"

    WILD_START=$(head -n 1 "$WILD_SITES" | cut -f2)
    WILD_END=$(tail -n 1 "$WILD_SITES" | cut -f2)

    # --- B. Clean and Prepare Elite VCF ---
    log "   Step 1: Sanitizing elite raw VCF..."
    grep "^##" "$ELITE_RAW" > "$WORK_TMP/header_meta.vcf" || true
    echo "##contig=<ID=$RAW_CHROM>" >> "$WORK_TMP/header_meta.vcf"
    grep "^#CHROM" "$ELITE_RAW" > "$WORK_TMP/header_columns.vcf"
    grep -v "^#" "$ELITE_RAW" | awk 'BEGIN{FS="\t"; OFS="\t"} {$7="."; print}' > "$WORK_TMP/body.vcf"

    cat "$WORK_TMP/header_meta.vcf" "$WORK_TMP/header_columns.vcf" "$WORK_TMP/body.vcf" \
        | bcftools view -Oz -o "$WORK_TMP/Elite_Import.vcf.gz"
    bcftools index "$WORK_TMP/Elite_Import.vcf.gz"

    log "   Step 2: Subsetting elite lines, renaming chromosomes and samples..."
    echo "$RAW_CHROM $TARGET_CHROM" > "$WORK_TMP/chr_map.txt"

    bcftools view -S "$SAMPLES_LIST" "$WORK_TMP/Elite_Import.vcf.gz" \
        | bcftools annotate --rename-chrs "$WORK_TMP/chr_map.txt" \
        | bcftools reheader -s "$RENAME_MAP" \
        | bcftools view -Oz -o "$WORK_TMP/Elite_Prep.vcf.gz"

    bcftools index "$WORK_TMP/Elite_Prep.vcf.gz"

    # --- C. Merge & Strict Intersection ---
    log "   Step 3: Merging and filtering to strict wild-site intersection..."
    MERGED_VCF="$OUT_DIR/${GENE}_Merged_Intersection.vcf.gz"

    bcftools merge -0 "$WILD_VCF" "$WORK_TMP/Elite_Prep.vcf.gz" -Ou \
        | bcftools view -T "$WILD_SITES" -Oz -o "$MERGED_VCF"

    bcftools index "$MERGED_VCF"

    # --- D. Sanity Stats ---
    WILD_SNPS=$(wc -l < "$WILD_SITES")
    MERGED_SNPS=$(bcftools query -f '%POS\n' "$MERGED_VCF" | wc -l)

    ELITE_TOTAL=$(bcftools query -f '%POS\n' "$WORK_TMP/Elite_Prep.vcf.gz" | wc -l)

    ELITE_IN_BOUNDS=$(bcftools view -r "${TARGET_CHROM}:${WILD_START}-${WILD_END}" "$WORK_TMP/Elite_Prep.vcf.gz" 2>/dev/null \
        | bcftools query -f '%POS\n' | wc -l)

    ELITE_MATCHING=$(bcftools view -T "$WILD_SITES" "$WORK_TMP/Elite_Prep.vcf.gz" 2>/dev/null \
        | bcftools query -f '%POS\n' | wc -l)

    LOST_BOUNDARIES=$((ELITE_TOTAL - ELITE_IN_BOUNDS))
    LOST_NOT_IN_WILD=$((ELITE_IN_BOUNDS - ELITE_MATCHING))

    LOST_MONO=$(bcftools view -r "${TARGET_CHROM}:${WILD_START}-${WILD_END}" -T ^"$WILD_SITES" "$WORK_TMP/Elite_Prep.vcf.gz" 2>/dev/null \
        | bcftools view -C 0 \
        | bcftools query -f '%POS\n' | wc -l)

    LOST_NOVEL=$((LOST_NOT_IN_WILD - LOST_MONO))

    ADDED_AS_REF=$((WILD_SNPS - ELITE_MATCHING))

    log "   -> Final Stats:"
    log "      Wild SNPs Target: $WILD_SNPS"
    log "      Final Merged SNPs: $MERGED_SNPS"
    log "      -----------------------------------"
    log "      Elite SNPs Found (Total): $ELITE_TOTAL"
    log "      Elite SNPs in Wild Bounds: $ELITE_IN_BOUNDS"
    log "      Elite SNPs perfectly matched: $ELITE_MATCHING"
    log "      Added as REF during merge: $ADDED_AS_REF"
    log "      -----------------------------------"
    log "      Lost due to boundaries:   $LOST_BOUNDARIES"
    log "      Lost within boundaries:   $LOST_NOT_IN_WILD"
    log "        ├─> Monomorphic Noise:  $LOST_MONO"
    log "        └─> True Novel Elite:   $LOST_NOVEL"
    log "      -----------------------------------"

    if [ "$WILD_SNPS" -eq "$MERGED_SNPS" ]; then
        log "      [SUCCESS] Perfect 1:1 wild-to-merged SNP count."
    else
        log "      [WARNING] Wild and merged SNP counts do not match."
    fi

    rm -rf "$WORK_TMP"
done

log "\n============================================================"
log "PIPELINE COMPLETE. All files saved to: $OUT_DIR"