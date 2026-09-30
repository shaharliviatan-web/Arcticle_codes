#!/usr/bin/env bash
# 01_extract_lead_genotypes.sh
# Extract the genotypes of the 36 GWAS lead SNPs for the 290 accessions.
#
# Inputs (read-only):
#   01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv
#   01_USED_GWAS_V2_pipeline/intermediates/morexV3_290.{bed,bim,fam}   (the GWAS genotype set)
# Output:
#   intermediates/lead_snps.txt              36 lead SNP IDs
#   intermediates/lead_genotypes.raw         PLINK --recode A, counted allele = PLINK A1 (minor)
#   intermediates/lead_genotypes.frq         PLINK allele frequencies of the 36 leads (MAF check)
set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp TEMP=/mnt/data/shahar/.tmp TMP=/mnt/data/shahar/.tmp

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PROJ="$(cd "$HERE/.." && pwd)"
PLINK=/usr/local/bin/plink            # bare `plink` resolves to a broken v0.76
LOCI="$PROJ/01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv"
BFILE="$PROJ/01_USED_GWAS_V2_pipeline/intermediates/morexV3_290"
OUT="$HERE/intermediates"

tail -n +2 "$LOCI" | cut -f4 > "$OUT/lead_snps.txt"
[[ $(wc -l < "$OUT/lead_snps.txt") -eq 36 ]] || { echo "expected 36 lead SNPs"; exit 1; }

# --recode A counts copies of A1. The bim/Table_loci_master A1 is the minor allele,
# and --keep-allele-order stops PLINK from re-deciding it.
"$PLINK" --bfile "$BFILE" --extract "$OUT/lead_snps.txt" --keep-allele-order --allow-extra-chr --memory 4000 \
         --recode A --out "$OUT/lead_genotypes" > "$HERE/logs/01_plink_recode.log" 2>&1
"$PLINK" --bfile "$BFILE" --extract "$OUT/lead_snps.txt" --keep-allele-order --allow-extra-chr --memory 4000 \
         --freq --out "$OUT/lead_genotypes" > "$HERE/logs/01_plink_freq.log" 2>&1
rm -f "$OUT"/lead_genotypes.nosex; mv "$OUT"/lead_genotypes.log "$HERE/logs/01_plink_last.log"
echo "[01] done: $(( $(head -1 "$OUT/lead_genotypes.raw" | wc -w) - 6 )) SNPs x $(( $(wc -l < "$OUT/lead_genotypes.raw") - 1 )) accessions"
