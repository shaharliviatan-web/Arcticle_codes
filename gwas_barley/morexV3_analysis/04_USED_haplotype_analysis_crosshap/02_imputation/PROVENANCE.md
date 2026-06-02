# Imputation provenance (step 04 reuses an existing genome-wide imputed VCF)

The imputed VCF used for LD in this pipeline was NOT regenerated here. It was produced earlier in
`/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_imputation/` and is reused as-is. The two source
scripts are copied here verbatim for provenance:

- `imputation_core.R`            LEA sNMF imputation core (K=3, 10 reps, best run by cross-entropy, method='mode')
- `run_imputation_pipeline.sh`   master wrapper (sample removal, sNMF, lfmm -> VCF reconstruction)

## What was produced
- Output: `USED_imputation/morexV3_final_imputed.vcf` (8.4 GB, 7,110,996 sites, 290 samples, GT only:
  `0/0` / `1/1`, no heterozygotes, no missing).
- Method: LEA `snmf(K=3, repetitions=10, entropy=TRUE)`; best run = `which.min(cross.entropy(...))`;
  `impute(method='mode', K=3)`. K=3 is the accepted ancestry count for this Southern-Levant wild
  barley panel.
- Samples: the original 300-sample VCF minus the 10 in `data/inputs/samples_to_remove_V3.txt`
  (HS0401-HS0411) = the 290 analysis accessions.

## Known issue (documented, harmless for this step)
`run_imputation_pipeline.sh` STEP 1 runs `plink --vcf ... --recode vcf` WITHOUT
`--keep-allele-order`, so plink re-set REF = major allele. As a result the imputed VCF has REF/ALT
FLIPPED relative to the original at the ~half of sites where ALT is the major allele
(e.g. 1H:49058 original C/T -> imputed T/C). It also doubled the sample IDs (HS0103 -> HS0103_HS0103)
because of `--recode vcf` defaults.

Why this is harmless in step 04:
- The imputed VCF is used ONLY to compute the PLINK LD matrix (`--r2`). r2 is invariant to per-site
  REF/ALT swaps, so the LD values are identical regardless of the flip.
- Haplotyping and the genotype heatmap use the RAW per-gene VCFs (cut with bcftools from the original
  `morexV3_with_ids.vcf.gz`), which keep the correct, unflipped REF/ALT and single-ID samples.
- Step 02 (`02_make_imputed_per_gene_vcfs.sh`) reheaders the per-gene imputed samples back to single
  IDs so they match the raw per-gene VCFs.

## If you ever regenerate the imputation cleanly (e.g. for the wild+elite comparison)
Re-run with allele order preserved so the imputed VCF matches the original at every site:

    plink --vcf morexV3.vcf.gz --remove samples_to_remove_V3.txt \
          --keep-allele-order --recode vcf --allow-extra-chr --out morexV3_samples_removed

Keep the rest of `run_imputation_pipeline.sh` / `imputation_core.R` unchanged. This removes the flip
and the need for any downstream allele harmonization.
