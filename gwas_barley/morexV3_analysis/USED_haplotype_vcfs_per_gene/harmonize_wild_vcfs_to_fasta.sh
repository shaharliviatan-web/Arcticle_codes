#!/bin/bash

# ==============================================================================
# SCRIPT: harmonize_wild_vcfs_to_fasta.sh (v3 - Full Scan)
# AUTHOR: Shahar / Gemini Analysis Partner
# DATE:   2026-03-03
#
# DESCRIPTION:
# This script harmonizes VCF files against a reference FASTA genome. It ensures
# that the REF allele in the VCF matches the base in the FASTA at that position.
#
# v3 CHANGES:
# - Removed the hard-coded GENE_FILES list.
# - The script now automatically finds and processes ALL *.vcf.gz files
#   within each specified trait directory.
# ==============================================================================

set -e
set -o pipefail

# --- USER CONFIGURATION ---

# Path to the indexed reference genome
FASTA_REF="/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_haplotype_vcfs_per_gene/Final_publication_genes/generate_elite_barley_VCFs/MorexV3.fa"

# Base directory containing the trait folders
BASE_VCF_DIR="/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_haplotype_vcfs_per_gene"

# List of trait folders to process
TRAIT_FOLDERS=(
    "starch_corrected_V3_targets_gene_220000bp"
    "starch_corrected_V3_imputed_VCFs_for_LD_targets_gene_220000bp"
    "fiber_corrected_V3_targets_gene_220000bp"
    "fiber_corrected_V3_imputed_VCFs_for_LD_targets_gene_220000bp"
    "betaglucan_corrected_V3_targets_gene_220000bp"
    "betaglucan_corrected_V3_imputed_VCFs_for_LD_targets_gene_220000bp"
    "protein_corrected_V3_targets_gene_220000bp"
    "protein_corrected_V3_imputed_VCFs_for_LD_targets_gene_220000bp"

)

# --- END OF CONFIGURATION ---


# --- SCRIPT LOGIC ---

log() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Check for tools
for tool in bcftools samtools bgzip; do
    if ! command -v $tool &> /dev/null; then
        log "ERROR: Required tool '$tool' is not installed or not in PATH."
        exit 1
    fi
done

if [ ! -f "${FASTA_REF}.fai" ]; then
    log "ERROR: FASTA index '${FASTA_REF}.fai' not found. Please index the reference."
    exit 1
fi

log "=== STARTING VCF HARMONIZATION (v3 - Full Scan) ==="

for trait_dir_name in "${TRAIT_FOLDERS[@]}"; do
    INPUT_DIR="${BASE_VCF_DIR}/${trait_dir_name}"
    OUTPUT_DIR="${INPUT_DIR}_HARMONIZED"
    
    if [ ! -d "$INPUT_DIR" ]; then
        log "WARNING: Input directory not found, skipping: $INPUT_DIR"
        continue
    fi
    
    mkdir -p "$OUTPUT_DIR"
    log "Processing Trait Directory: $trait_dir_name"
    log "Output will be saved to: $OUTPUT_DIR"

    # Find all .vcf.gz files in the input directory and loop through them
    find "$INPUT_DIR" -maxdepth 1 -type f -name "*.vcf.gz" | while read -r INPUT_VCF; do
        
        gene_file=$(basename "$INPUT_VCF")
        OUTPUT_VCF="${OUTPUT_DIR}/${gene_file}"
        
        log "  -> Harmonizing Gene: $gene_file"
        
        TEMP_DIR=$(mktemp -d)
        # Setup cleanup for temporary directory on script exit
        trap 'rm -rf -- "$TEMP_DIR"' EXIT

        POSITIONS_FILE="$TEMP_DIR/positions.txt"
        FASTA_MAP_FILE="$TEMP_DIR/fasta_map.txt"
        MISMATCH_LOG="$OUTPUT_DIR/${gene_file}.mismatch.log"
        > "$MISMATCH_LOG"

        # 1. Create a map of true reference bases from the FASTA
        bcftools query -f '%CHROM\t%POS\n' "$INPUT_VCF" > "$POSITIONS_FILE"
        
        while IFS=$'\t' read -r chrom pos; do
            ref_base=$(samtools faidx "$FASTA_REF" "${chrom}:${pos}-${pos}" | sed -n '2p' | tr 'a-z' 'A-Z')
            echo -e "${chrom}:${pos}\t${ref_base}"
        done < "$POSITIONS_FILE" > "$FASTA_MAP_FILE"

        # 2. Process the VCF with awk for harmonization
        bcftools view "$INPUT_VCF" | awk -v mismatch_log="$MISMATCH_LOG" '
            BEGIN {
                FS=OFS="\t";
                while ((getline < "'"$FASTA_MAP_FILE"'") > 0) {
                    fasta_map[$1] = $2;
                }
                close("'"$FASTA_MAP_FILE"'");
            }

            /^#/ { print; next; }
            
            {
                key = $1 ":" $2;
                vcf_ref = $4;
                vcf_alt = $5;
                
                if (vcf_alt ~ /,/ || length(vcf_ref) != 1 || length(vcf_alt) != 1) {
                    print "WARNING: Skipping complex site " key " (multi-allelic or indel)" > "/dev/stderr";
                    next;
                }

                if (key in fasta_map) {
                    fasta_ref = fasta_map[key];

                    if (vcf_ref == fasta_ref) {
                        print;
                    }
                    else if (vcf_alt == fasta_ref) {
                        $4 = vcf_alt;
                        $5 = vcf_ref;
                        
                        for (i=10; i<=NF; i++) {
                            split($i, parts, ":");
                            gt = parts[1];
                            
                            separator = "/";
                            if (index(gt, "|")) {
                                separator = "|";
                            }

                            if (gt == "0" separator "0") {
                                gt = "1" separator "1";
                            } else if (gt == "1" separator "1") {
                                gt = "0" separator "0";
                            }
                            
                            parts[1] = gt;
                            $i = parts[1];
                            for(j=2; j<=length(parts); j++) $i = $i ":" parts[j];
                        }
                        print;
                    }
                    else {
                        print "MISMATCH at " key ": VCF_REF=" vcf_ref ", VCF_ALT=" vcf_alt ", FASTA_REF=" fasta_ref >> mismatch_log;
                    }
                } else {
                    print "WARNING: Position " key " not found in FASTA map. Skipping." > "/dev/stderr";
                }
            }
        ' | bgzip -c > "$OUTPUT_VCF"
        
        bcftools index -f "$OUTPUT_VCF"
        log "    -> Done. Harmonized VCF created at $OUTPUT_VCF"
        if [ -s "$MISMATCH_LOG" ]; then
            log "    -> WARNING: Mismatches were found and logged to ${MISMATCH_LOG}"
        else
            rm "$MISMATCH_LOG"
        fi
        
        trap - EXIT
    done
done

log "=== HARMONIZATION COMPLETE ==="
