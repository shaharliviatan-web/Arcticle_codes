# 05_USED_gene_annotation_analysis

Functional annotation of the **23 haplotype-significant genes** from step 04
(run `loci_LDspan_eps06_V4`). Takes a gene accession and returns what the protein is,
with the evidence behind that call.

**Start here → [`08_USED_annotation_master/results/tables/Table_significant_genes_paper.tsv`](08_USED_annotation_master/results/tables/Table_significant_genes_paper.tsv)**

Last run 2026-09-10. Each sub-step keeps its previous 45-gene results in its own
`_ARCHIVE_v1_45genes_2026-09-09/`.

---

## The chain

Three independent evidence sources, checked in a fixed order, then merged. The order
**is** the check number reported in the final table.

| # | directory | method | local/remote | runtime |
|---|---|---|---|---|
| **1** | `05_USED_gene_annotation/` | DIAMOND BLASTP vs **UniProt Swiss-Prot** 2026_03 (corrected 2026-09-30 from 2026_01) | local (280 MB DB built) | ~25 s |
| **2** | `06_USED_interpro_domains/` | **EBI InterProScan 5** — Pfam/InterPro domains + **GO** | remote, 1 job per protein | ~15 min |
| **3** | `07_USED_BLASTP_genes_with_no_annotation_left/` | **NCBI nr** BLASTP, Viridiplantae — rescue only | remote | minutes to hours |
| — | `08_USED_annotation_master/` | merge → one auditable table per gene | local | s |
| — | `09_USED_canonical_betaglucan_gene_check/` | were the canonical CslF/CslH genes missed? **⚠ stale** | local | — |

```bash
cd 05_USED_gene_annotation                  && bash scripts/run_all.sh
cd ../06_USED_interpro_domains              && bash scripts/run_all.sh   # use screen
cd ../07_USED_BLASTP_genes_with_no_annotation_left && bash scripts/run_all.sh
cd ../08_USED_annotation_master             && bash scripts/run_all.sh && Rscript scripts/01_build_paper_table.R
```

Steps 2 and 3 are **resumable**: InterProScan caches one JSON per sequence, nr writes a
`.done` marker. Re-running only fills gaps.

## Result

| | |
|---|---|
| genes in | **23** (beta-glucan 6, fiber 12, starch 5) |
| named by check 1 (Swiss-Prot) | **16** |
| named by check 2 (InterPro) | **5** |
| no call from any source | **2** |
| carrying GO terms | 17 |
| carrying Pfam domains | 20 |

Protein contributes no genes — it has no testable gene in step-04 run V4.

## No confidence tier

Removed 2026-09-10. It conflated *how sure we are what the protein is* with *how
plausible the gene is for the trait*, which are independent, and hid the evidence.
The table now reports plain provenance — `annotation_source`, `annotation_from_check`
(1/2/3) — plus **every source's own call** in `call_check1_swissprot`,
`call_check2_interpro`, `call_check3_nr`, so when two sources differ both texts are
visible and can be judged together. `needs_review_two_calls` flags such pairs; it
decides nothing.

Trait candidacy is a separate manual judgement, recorded in
`08_USED_annotation_master/results/tables/trait_candidate_calls.tsv` (see that
directory's README).

## The gene list is not hard-coded

Every sub-step derives its input. v1 hard-coded a fixed count of 45, a specific gene's
peptide length, a literal 5-gene residual list, and a hand-made per-sequence FASTA
split — each of which silently broke or aborted when the gene set changed. All four are
now derived and asserted. See the "set-size independence" note in each sub-README.

## Where the input comes from

`04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Significant_genes/significant_genes.tsv`

Read by `05_USED_gene_annotation/scripts/00_lock_fdr_genes.R`, which also carries the
step-04 statistics (`fdr_q`, `eta_squared`, `delta_top_bottom_sd`, group counts) through
to the final table — a functional call is only interpretable beside the effect it explains.

**To re-point at a different step-04 run:** change the path in that one script.

## ⚠ `09_USED_canonical_betaglucan_gene_check/` is stale

It asks whether the canonical (1,3;1,4)-β-glucan genes (CslF6, CslH1, …) were missed by
the candidate-gene window, and is written around the retired **±200 kb** window. Step 03
now searches the LD locus span with **no flanking window**, so its conclusion no longer
describes the current pipeline. It was not re-run and is not part of the chain above.

## Conventions

- All paths absolute. `TMPDIR=/mnt/data/shahar/.tmp`; nothing under `/tmp` or `$HOME`.
- DIAMOND 2.0.14 (`/usr/bin/diamond`); proteome = PGSB Morex V3 r3 HC, longest isoform
  per gene.
- Long remote jobs run under GNU `screen`.
