#!/bin/bash
# ============================================================================
# 02_extract_genes.sh
#
# WHAT   Intersect the search intervals with the Morex V3 gene annotation.
#
# READS  intermediates/loci_intervals.bed   (from 01_build_intervals.R)
#        $GFF_GZ                            Morex V3 GFF3, Ensembl Plants r62
#
# WRITES intermediates/genes_1to7H.gff      gene features on 1H-7H only
#        intermediates/intersect_raw.tsv    bedtools intersect -wa -wb
#        intermediates/genes_per_locus.tsv  gene COUNT per locus, zeros included
#
# COORDINATES  loci_intervals.bed is 0-based half-open; the GFF is 1-based
#        inclusive. bedtools infers this from the .gff extension and converts
#        internally, so the intersection is coordinate-correct. Do not rename
#        genes_1to7H.gff to .bed or .txt -- that silently breaks the offset.
#
# SCOPE  Only col3 == "gene" (not mRNA/CDS/exon, which would multiply-count a
#        gene by its transcripts). Unplaced CAJHDD* scaffolds are excluded:
#        they carry no GWAS SNPs, so no locus can fall on them.
# ============================================================================

set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/../config/params.sh"

BED="$STEP03_BASE/intermediates/loci_intervals.bed"
GENES="$STEP03_BASE/intermediates/genes_1to7H.gff"
RAW="$STEP03_BASE/intermediates/intersect_raw.tsv"
COUNTS="$STEP03_BASE/intermediates/genes_per_locus.tsv"

echo "=== 02_extract_genes.sh ==="

[ -s "$BED" ] || { echo "ERROR: $BED missing or empty. Run 01_build_intervals.R first." >&2; exit 1; }
[ -s "$GFF_GZ" ] || { echo "ERROR: GFF not found: $GFF_GZ" >&2; exit 1; }
command -v "$BEDTOOLS" >/dev/null || { echo "ERROR: bedtools not found at $BEDTOOLS" >&2; exit 1; }

echo "GFF      : $GFF_GZ"
echo "bedtools : $("$BEDTOOLS" --version)"

# --- Step 1: gene features on 1H..7H -----------------------------------------
zcat "$GFF_GZ" | awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && $1 ~ /^[1-7]H$/' > "$GENES"
echo "Gene features on 1H-7H : $(wc -l < "$GENES")"

if grep -q "CAJHDD" "$GENES"; then
  echo "ERROR: an unplaced CAJHDD scaffold leaked into the filtered gene set." >&2
  exit 1
fi

# --- Step 2: intersect --------------------------------------------------------
"$BEDTOOLS" intersect -a "$BED" -b "$GENES" -wa -wb > "$RAW"
echo "Intersect rows (gene x locus) : $(wc -l < "$RAW")"

# --- Step 3: per-locus counts, INCLUDING loci that hit nothing -----------------
# -c reports 0 for empty intervals; -wa -wb drops them entirely. Both are needed:
# the zeros are a reportable result here, not an absence of data.
"$BEDTOOLS" intersect -a "$BED" -b "$GENES" -c \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "locus_id","n_genes"} {print $4, $5}' > "$COUNTS"
echo "Loci with zero genes : $(awk -F'\t' 'NR>1 && $2==0' "$COUNTS" | wc -l) of $(( $(wc -l < "$COUNTS") - 1 ))"

echo "Wrote intermediates/genes_1to7H.gff, intersect_raw.tsv, genes_per_locus.tsv"
