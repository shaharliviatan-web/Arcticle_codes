#!/bin/bash
# 02_extract_genes.sh
# Intersect candidate-interval BED (from 01_build_loci.R) with Morex V3 GFF3 genes.
# Genes: col3 == "gene", chromosomes 1H..7H only (CAJHDD* scaffolds excluded).
#
# Outputs:
#   intermediates/genes_1to7H.gff      filtered gene records (1-based GFF)
#   intermediates/intersect_raw.tsv    bedtools intersect -wa -wb output
#
# Coordinates: loci_intervals.bed is 0-based half-open; the GFF is 1-based.
# bedtools recognises the .gff extension and converts internally, so the
# intersection is coordinate-correct.

set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp

BEDTOOLS=/usr/bin/bedtools
BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/03_USED_candidate_genes_around_leading_snps
GFF=/mnt/data/shahar/gwas_barley/morexV3_analysis/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3.gz

BED="$BASE/intermediates/loci_intervals.bed"
GENES="$BASE/intermediates/genes_1to7H.gff"
RAW="$BASE/intermediates/intersect_raw.tsv"

if [ ! -s "$BED" ]; then
  echo "ERROR: $BED missing or empty. Run 01_build_loci.R first." >&2
  exit 1
fi

# Step 1: filter GFF to gene features on chromosomes 1H..7H only.
zcat "$GFF" | awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" && $1 ~ /^[1-7]H$/' > "$GENES"
echo "Filtered genes (1H-7H): $(wc -l < "$GENES")" >&2

# Sanity: no scaffolds leaked in.
if grep -q "CAJHDD" "$GENES"; then
  echo "ERROR: CAJHDD scaffold leaked into filtered genes." >&2
  exit 1
fi

# Step 2: intersect intervals (-a, BED) with genes (-b, GFF), keep both sides.
"$BEDTOOLS" intersect -a "$BED" -b "$GENES" -wa -wb > "$RAW"
echo "Intersect rows (gene x locus): $(wc -l < "$RAW")" >&2
echo "Wrote $GENES and $RAW" >&2
