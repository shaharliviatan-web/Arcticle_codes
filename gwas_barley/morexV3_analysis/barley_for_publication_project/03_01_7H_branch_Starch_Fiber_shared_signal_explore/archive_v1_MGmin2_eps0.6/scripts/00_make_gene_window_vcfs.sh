#!/usr/bin/env bash
# 00_make_gene_window_vcfs.sh - raw + imputed per-gene VCFs for 7HG0729030 (gene +/- 1000 bp).
# Mirrors step 04's 01_/02_ scripts exactly: biallelic SNPs only, the same 290-sample
# keep-list, and the same doubled->single sample rename for the imputed VCF.
# Reuses step 04's already-staged bgzipped+indexed imputed VCF (read-only).
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"
INT="$BRANCH/intermediates"; mkdir -p "$INT"
REGION="${GENE_CHR}:$((GENE_START-WINDOW_BP))-$((GENE_END+WINDOW_BP))"
echo "region: $REGION"
"$BCFTOOLS" view -r "$REGION" -S "$KEEP_290" -m2 -M2 -v snps -Oz -o "$INT/gene_window_raw.vcf.gz" "$SOURCE_VCF"
"$BCFTOOLS" index -f "$INT/gene_window_raw.vcf.gz"
"$BCFTOOLS" view -r "$REGION" -m2 -M2 -v snps -Oz -o "$TMPDIR/g30_imp_tmp.vcf.gz" "$STAGED_IMPUTED"
"$BCFTOOLS" reheader -s "$SAMPLE_RENAME" "$TMPDIR/g30_imp_tmp.vcf.gz" \
  | "$BCFTOOLS" view -S "$KEEP_290" -Oz -o "$INT/gene_window_imputed.vcf.gz"
"$BCFTOOLS" index -f "$INT/gene_window_imputed.vcf.gz"; rm -f "$TMPDIR/g30_imp_tmp.vcf.gz"
"$PLINK" --vcf "$INT/gene_window_imputed.vcf.gz" --r2 square --keep-allele-order \
  --allow-extra-chr --double-id --silent --out "$INT/gene_window_ld_r2_square" 
mv -f "$INT/gene_window_ld_r2_square.ld" "$INT/gene_window_ld_r2_square.ld" 2>/dev/null || true
echo "raw SNPs     : $("$BCFTOOLS" view -H "$INT/gene_window_raw.vcf.gz" | wc -l)"
echo "imputed SNPs : $("$BCFTOOLS" view -H "$INT/gene_window_imputed.vcf.gz" | wc -l)"
echo "00 OK"
