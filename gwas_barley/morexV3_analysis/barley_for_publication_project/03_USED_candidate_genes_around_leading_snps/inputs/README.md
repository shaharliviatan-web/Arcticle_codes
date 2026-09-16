# inputs/

**Nothing here is hand-maintained.** This directory holds one auto-generated
provenance file.

## `loci_handoff_snapshot.tsv`

A frozen copy of the step-01 handoff table, written by `01_build_intervals.R` at
every run. The first four lines are `#` comments recording the source path, the
source file's mtime, and when the copy was taken.

Its purpose is provenance: the pipeline always *reads* the live table at
`01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_for_gene_search.tsv`,
so step 01 stays the single source of truth — but if step 01 is ever re-run, this
snapshot still records exactly which locus set produced the tables in `results/`.

## What used to be here

v1 kept a hand-curated `marginal_snps.tsv` (9 SNPs just below the old 6.7712
threshold) and a generated `curated_snps.tsv`.

They were obsolete because the v3 LD pruning moved the Bonferroni threshold from
6.7712 to **6.0454**, which makes most of those "marginal" SNPs simply significant.
The category no longer exists, and step 01's `class` column replaces it. Deleted
2026-09-10 with the rest of the v1 material.
