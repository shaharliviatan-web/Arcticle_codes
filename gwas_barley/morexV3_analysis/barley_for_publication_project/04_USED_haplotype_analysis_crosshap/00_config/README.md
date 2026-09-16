# 00_config/

Everything that parameterises step 04. **`config.yaml` is the single source of truth** —
both the shell scripts and the R scripts read it.

| file | role |
|---|---|
| **`config.yaml`** | all parameters, with the reasoning for each choice in comments |
| `gene_windows.tsv` | **generated** by `01_scripts/00_build_gene_windows.R` from step 03 — one row per gene with its ± `window_bp` window and its GWAS origin. Do not hand-edit |
| `samples_keep_290.txt` | **generated** — the 290 accessions (source VCF minus `samples_remove`) |
| `final_genes.tsv` | **hand-curated** publication gene list — ⚠ **stale**, see below |
| `final_genes.tsv.template` | template for rebuilding that list |

## The two parameters that matter

`mgmin_values: [2]` and `epsilon_vector: [0.6]` — **exactly one value each**. The
pipeline refuses to start otherwise. Chosen on genotype-only grounds; the full reasoning
is in `config.yaml`'s comments and in the top-level [`../README.md`](../README.md).

Note that **epsilon is not an r² cutoff** — it is a DBSCAN Euclidean radius over r²
profiles, so its stringency varies with SNP density in the window.

## Current state (2026-09-08)

`gene_windows.tsv` holds **55 genes** from step 03 (run `loci_LDspan_eps06_V4`):

| trait | genes | testable in step 04 |
|---|---|---|
| betaglucan | 7 | 6 |
| fiber | 32 | 19 |
| protein | 6 | **0** |
| starch | 10 | 5 |

Protein has no testable gene: 3 of its 6 candidates have zero SNPs in window, 1 has a
single SNP, and 2 fail to cluster. Verified not to be a parameter artefact — see
`../04_runs/loci_LDspan_eps06_V4/Diagnostics/protein_epsilon_rescue_scan.tsv`.

## ⚠ `final_genes.tsv` is stale

Its four curated genes were selected from v1's 108-gene candidate set. Step 03 now uses
the step-01 LD locus span with **no flanking window**, so three of them are no longer
candidates at all:

| gene | symbol | trait | in the current 64? |
|---|---|---|---|
| HORVU.MOREX.r3.3HG0301710 | PHT | starch | **yes** |
| HORVU.MOREX.r3.3HG0301750 | Pho | starch | no — 102,704 bp outside its locus |
| HORVU.MOREX.r3.3HG0299440 | AP2/ERF | betaglucan | no — 69,298 bp outside |
| HORVU.MOREX.r3.7HG0642350 | BAHD | fiber | no — 34,865 bp outside |

This was accepted deliberately (see step 03's README). The list must be re-curated from
`04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv` before the publication table and
figures are rebuilt.
