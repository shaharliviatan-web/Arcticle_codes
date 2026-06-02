#!/bin/bash
# 02_make_imputed_per_gene_vcfs.sh
# Cut per-gene IMPUTED VCFs (same gene+/-window regions) from the genome-wide LEA
# sNMF K=3 imputed VCF. These are used ONLY for PLINK LD (r2 is invariant to the
# REF/ALT flip that plink introduced in the imputed VCF, so it is harmless here).
#
# Samples are reheadered from doubled (HS0103_HS0103) to single (HS0103) IDs so they
# match the raw per-gene VCFs from script 01.
#
# Input : 00_config/gene_windows.tsv, 00_config/samples_keep_290.txt, imputed_vcf (config)
# Output: 03_per_gene_vcfs/imputed_1000bp/{trait}/{gene_id}.vcf.gz (+ .csi)
#         03_per_gene_vcfs/imputed_1000bp_manifest.tsv
#         03_per_gene_vcfs/staging/morexV3_final_imputed.vcf.gz (+ .csi)  [one-time bgzip]

set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/04_USED_haplotype_analysis_crosshap
CFG="$BASE/00_config/config.yaml"
cfgval() { grep -E "^[[:space:]]*$1:" "$CFG" | head -1 | sed -E "s/^[[:space:]]*$1:[[:space:]]*\"?([^\"]*)\"?[[:space:]]*$/\1/"; }

BCFTOOLS=$(cfgval bcftools_bin)
IMP_VCF=$(cfgval imputed_vcf)
OUTDIR="$BASE/03_per_gene_vcfs/imputed_1000bp"
STAGE_DIR="$BASE/03_per_gene_vcfs/staging"
STAGE="$STAGE_DIR/morexV3_final_imputed.vcf.gz"
GW="$BASE/00_config/gene_windows.tsv"
KEEP="$BASE/00_config/samples_keep_290.txt"
MAP="$STAGE_DIR/sample_rename_doubled_to_single.txt"
MANIFEST="$BASE/03_per_gene_vcfs/imputed_1000bp_manifest.tsv"

mkdir -p "$STAGE_DIR"
[ -s "$GW" ] || { echo "ERROR: missing $GW" >&2; exit 1; }
[ -s "$IMP_VCF" ] || { echo "ERROR: missing imputed VCF $IMP_VCF" >&2; exit 1; }

# --- One-time: bgzip + index the genome-wide imputed VCF (enables region queries) ---
if [ ! -s "$STAGE" ] || [ ! -s "$STAGE.csi" ]; then
  echo "[stage] bgzip imputed VCF (8.4G, one-time)..."
  bgzip -@ 8 -c "$IMP_VCF" > "$STAGE"
  echo "[stage] index..."
  "$BCFTOOLS" index -f "$STAGE"
else
  echo "[stage] using existing $STAGE"
fi

# --- Reheader map: doubled (HS0103_HS0103) -> single (HS0103) ---
"$BCFTOOLS" query -l "$STAGE" | awk '{s=$1; sub(/_.*/,"",s); print $1"\t"s}' > "$MAP"
echo "[stage] reheader map: $(wc -l < "$MAP") samples"

printf "trait\tgene_id\tregion\tn_snps\tn_samples\n" > "$MANIFEST"

n=0
# columns: 1 trait, 2 gene_id, 3 chr, 4 gene_start, 5 gene_end, 6 strand, 7 win_start, 8 win_end
tail -n +2 "$GW" | while IFS=$'\t' read -r trait gene_id chr gstart gend strand wstart wend rest; do
  n=$((n+1))
  mkdir -p "$OUTDIR/$trait"
  out="$OUTDIR/$trait/${gene_id}.vcf.gz"
  tmp="$TMPDIR/imp_${trait}_${gene_id}.vcf.gz"
  region="${chr}:${wstart}-${wend}"

  "$BCFTOOLS" view -r "$region" -m2 -M2 -v snps -Oz -o "$tmp" "$STAGE"
  "$BCFTOOLS" reheader -s "$MAP" "$tmp" | "$BCFTOOLS" view -S "$KEEP" -Oz -o "$out"
  "$BCFTOOLS" index -f "$out"
  rm -f "$tmp"

  nsnp=$("$BCFTOOLS" view -H "$out" | wc -l)
  nsamp=$("$BCFTOOLS" query -l "$out" | wc -l)
  printf "%s\t%s\t%s\t%s\t%s\n" "$trait" "$gene_id" "$region" "$nsnp" "$nsamp" >> "$MANIFEST"
  printf "[%3d] %-12s %-28s %-26s snps=%-5s samples=%s\n" "$n" "$trait" "$gene_id" "$region" "$nsnp" "$nsamp"
done

echo "----"
echo "Wrote imputed per-gene VCFs to: $OUTDIR"
echo "Manifest: $MANIFEST"
echo "Total genes: $(($(wc -l < "$MANIFEST")-1))"
