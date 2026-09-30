# loci_LDspan_eps06_V4 — the live run

crosshap haplotype analysis, **MGmin = 2, epsilon = 0.6**, both fixed. Created 2026-09-09.

**Start here → [`Stats/results_chapter_numbers.txt`](Stats/results_chapter_numbers.txt)**
**Just the hits → [`Significant_genes/`](Significant_genes/)**

| | |
|---|---|
| genes in (step 03) | 55 |
| genes tested | 30 |
| genes not testable | 25 |
| significant, BH q ≤ 0.05 | **23 genes in 12 loci** |

| trait | genes in | tested | significant |
|---|---|---|---|
| beta-glucan | 7 | 6 | 6 |
| fiber | 32 | 19 | 12 |
| protein | 6 | **0** | **0** |
| starch | 10 | 5 | 5 |

## Why epsilon = 0.6

Fixed in advance, applied to every gene, chosen on **genotype-only** grounds — crosshap
builds haplotype groups from genotypes alone, so no p-value entered the choice.

Assignment rate and gene coverage trade off **monotonically against each other**
(measured on this run's 55 candidate genes,
[`Diagnostics/epsilon_choice_supplementary.tsv`](Diagnostics/epsilon_choice_supplementary.tsv);
per-gene data in `Diagnostics/epsilon_choice_per_gene.tsv`; both written by
`01_scripts/diagnostics/parameter_choice_tables.R`, 2026-09-30):

| ε | genes testable | accessions assigned | median groups |
|---|---|---|---|
| 0.2 | 23 / 55 | 70.0% | 2 |
| 0.4 | 27 / 55 | 69.3% | 3 |
| **0.6** | **30 / 55** | **65.0%** | **3** |
| 0.8 | 34 / 55 | 60.9% | 4 |
| 1.0 | 33 / 55 | 60.3% | 4 |

A low ε assigns a high proportion of accessions, but only in the minority of genes whose
variants cluster tightly — 23 of 55 at ε=0.2, and with a median of just 2 groups. A high
ε makes more genes testable but leaves more accessions unassigned in each. Neither
extreme is preferable: an unassigned accession contributes nothing to a gene's test, and
an untestable gene contributes nothing at all.

**0.6 sits at the midpoint** — 65% of accessions assigned while 55% of candidate genes
become testable, at a median of 3 haplotype groups. It is **not the maximum of either
metric**, and the methods text should not claim it is; the argument is balance, uniform
application, and pre-declaration.

> Superseded, recorded so it is not reintroduced: an earlier version claimed 0.6 had the
> highest assignment rate (~70%). That was measured on a retired 64-gene set and on a
> grid omitting ε < 0.6. On the real gene set, 0.2 leads on assignment and 0.8 on
> coverage.

**ε = 0.85 was rejected** because it maximised significance for PHT4;3 off its stable
plateau. Re-checked 2026-09-30 against three genotype-only criteria — coverage/assignment,
group evenness, smallest-group size — **none selects 0.85**, and it is worse than 0.6 on
evenness. Choosing it would be selection on the outcome.

## What this run is

Step 01's final locus definition: `--clump-kb 2000` (±2 Mb), r² ≥ 0.5, members on LD
alone, severed iteratively at the first internal gap > **50 kb**, and **no sub-threshold
peaks** — the 2 extra protein peaks earlier runs carried were dropped by decision.
Step 03 then searched the LD span itself with **no flanking window** (`FLANK_BP = 0`).

## Protein has no result in this run

Its 2 significant loci give 6 candidate genes, none testable: 3 have **zero** SNPs in
window, 1 has a single SNP, and 2 fail to cluster at ε = 0.6. Protein's only previous hit
came from a sub-threshold peak that is no longer included. This is a reportable outcome,
not a pipeline failure — see `Stats/genes_not_tested.tsv`.

## Layout

| dir | contents |
|---|---|
| `Stats/` | all result tables — `gene_results.tsv` is the main one |
| `Significant_genes/` | **the 23 hits**, ranked, with both PDFs each |
| `CombinedPDF/`, `Heatmaps/` | per-gene figures for all 30 tested genes |
| `Cache/<trait>/<gene>/MGmin_2/HapObject.rds` | crosshap objects — delete to force a re-run |
| `Logs/` | timestamped event log |

## Interpreting the counts

**23 significant genes sit in 12 loci.** Genes within a locus are in LD and are not
independent discoveries — report loci alongside genes. BH stays valid under that
correlation (positive regression dependence); the caveat is one of reporting.

**Read the effect sizes, not just q.** These genes sit at loci already selected for
association with the trait, so a small p is close to guaranteed. `eta_squared` and
`delta_top_bottom_sd` carry the information. In this run they split the list sharply:
ranks 1–12 have η² 0.08–0.28, ranks 18–23 have η² 0.019–0.035.
