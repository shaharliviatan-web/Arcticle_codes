#!/usr/bin/env bash
# 02_fetch_elite_vcf.sh - elite panel genotypes for this gene's window.
# Mirrors 07_.../scripts/01_fetch_elite_vcfs.sh: same DivBrowse endpoint, same
# 136-line spring elite pool (reused from step 07's cached inputs, read-only).
# NOTE the pad is 2 kb, not step 07's 50 kb: this region is variant-dense and the
# larger request streams for many minutes. Skips if the output already exists.
set -euo pipefail
S="$(dirname "$(readlink -f "$0")")"; source "$S/config.sh"
D="$TEMP_ROOT/work/elite"; mkdir -p "$D"
OUT="$D/7HG0729030.elite.vcf.gz"
if [ -s "$OUT" ]; then echo "elite VCF already present: $OUT"; exit 0; fi
POOL="$STEP07_ROOT/inputs/elite_pool_samples.tsv"
[ -s "$POOL" ] || { echo "ERROR: $POOL missing (run step 07's 00_panel_metadata.sh)" >&2; exit 1; }
SAMPLES="[$(awk -F'\t' 'NR>1{printf "%s\"%s\"",(n++?",":""),$1}' "$POOL")]"
echo "requesting chr${GENE_CHR}:$((WIN_START-2000))-$((WIN_END+2000)) for $(( $(wc -l < "$POOL") - 1 )) samples"
curl -s -m 1800 -X POST "${DIVBROWSE_BASE}/vcf_export" \
  --data-urlencode "chrom=chr${GENE_CHR}" --data-urlencode "startpos=$((WIN_START-2000))" \
  --data-urlencode "endpos=$((WIN_END+2000))" --data-urlencode "samples=${SAMPLES}" \
  -o "$D/raw.vcf"
# DivBrowse leaves INFO empty, which vcfR cannot parse; fill it, then trim to the window
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print; next} {if($8=="") $8="."; print}' "$D/raw.vcf" | "$BGZIP" -c > "$D/raw.vcf.gz"
"$BCFTOOLS" index -f "$D/raw.vcf.gz"
"$BCFTOOLS" view -m2 -M2 -v snps -r "chr${GENE_CHR}:${WIN_START}-${WIN_END}" -Oz -o "$OUT" "$D/raw.vcf.gz"
"$BCFTOOLS" index -f "$OUT"; rm -f "$D/raw.vcf" "$D/raw.vcf.gz"*
echo "elite variants in window: $("$BCFTOOLS" view -H "$OUT" | wc -l)"
