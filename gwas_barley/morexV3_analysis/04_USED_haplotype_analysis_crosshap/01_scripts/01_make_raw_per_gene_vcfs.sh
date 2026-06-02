#!/bin/bash
# 01_make_raw_per_gene_vcfs.sh
# Cut per-gene RAW VCFs (gene_start-window .. gene_end+window) from the original
# wild-barley VCF, subset to the 290 analysis samples.
#
# Uses bcftools (NOT plink): preserves REF/ALT exactly (no allele flip) and keeps
# single-ID sample names (HS0103). Biallelic SNPs only.
#
# Input : 00_config/gene_windows.tsv, 00_config/samples_keep_290.txt
# Output: 03_per_gene_vcfs/raw_1000bp/{trait}/{gene_id}.vcf.gz (+ .csi)
#         03_per_gene_vcfs/raw_1000bp_manifest.tsv

set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap
CFG="$BASE/00_config/config.yaml"
cfgval() { grep -E "^[[:space:]]*$1:" "$CFG" | head -1 | sed -E "s/^[[:space:]]*$1:[[:space:]]*\"?([^\"]*)\"?[[:space:]]*$/\1/"; }

BCFTOOLS=$(cfgval bcftools_bin)
SOURCE_VCF=$(cfgval source_vcf)
OUTDIR="$BASE/03_per_gene_vcfs/raw_1000bp"
GW="$BASE/00_config/gene_windows.tsv"
KEEP="$BASE/00_config/samples_keep_290.txt"
MANIFEST="$BASE/03_per_gene_vcfs/raw_1000bp_manifest.tsv"

[ -s "$GW" ] || { echo "ERROR: missing $GW (run 00_build_gene_windows.R first)" >&2; exit 1; }
[ -s "$KEEP" ] || { echo "ERROR: missing $KEEP" >&2; exit 1; }
[ -s "$SOURCE_VCF" ] || { echo "ERROR: missing source VCF $SOURCE_VCF" >&2; exit 1; }

printf "trait\tgene_id\tregion\tn_snps\tn_samples\n" > "$MANIFEST"

n=0
# columns: 1 trait, 2 gene_id, 3 chr, 4 gene_start, 5 gene_end, 6 strand, 7 win_start, 8 win_end
tail -n +2 "$GW" | while IFS=$'\t' read -r trait gene_id chr gstart gend strand wstart wend rest; do
  n=$((n+1))
  mkdir -p "$OUTDIR/$trait"
  out="$OUTDIR/$trait/${gene_id}.vcf.gz"
  region="${chr}:${wstart}-${wend}"

  "$BCFTOOLS" view -r "$region" -S "$KEEP" -m2 -M2 -v snps \
      -Oz -o "$out" "$SOURCE_VCF"
  "$BCFTOOLS" index -f "$out"

  nsnp=$("$BCFTOOLS" view -H "$out" | wc -l)
  nsamp=$("$BCFTOOLS" query -l "$out" | wc -l)
  printf "%s\t%s\t%s\t%s\t%s\n" "$trait" "$gene_id" "$region" "$nsnp" "$nsamp" >> "$MANIFEST"
  printf "[%3d] %-12s %-28s %-26s snps=%-5s samples=%s\n" "$n" "$trait" "$gene_id" "$region" "$nsnp" "$nsamp"
done

echo "----"
echo "Wrote raw per-gene VCFs to: $OUTDIR"
echo "Manifest: $MANIFEST"
echo "Total genes: $(($(wc -l < "$MANIFEST")-1))"
