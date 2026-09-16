# 01_USED_GWAS_V2_pipeline

GWAS of 290 wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*, Southern
Levant) against MorexV3, 7,110,996 SNPs, four grain-quality traits: beta-glucan,
fiber, protein, starch.

**Start here → [`results/00_FINAL_BLUP_3PC/`](results/00_FINAL_BLUP_3PC/)** — the chosen
configuration and the 36 loci that the candidate-gene search runs on.

Reorganised 2026-09-08.

---

## The result, in one table

| | |
|---|---|
| Configuration | BLUP phenotypes × 3 PCs, EMMAX aIBS kinship (290 × 290) |
| Significance | Bonferroni α = 0.10 over 111,017 LD-pruned SNPs → **−log10p ≥ 6.0454** (p < 9.008e-07) |
| Significant SNPs | **52** |
| **Loci** | **36** |
| Median locus span | **93.5 kb** (mean 219.8, max 1,727.8) |
| Member SNPs across all loci | 3,423 |
| Significant SNPs inside a locus | **52 / 52** (enforced by assertion) |
| Total gene-search space | **7.9 Mb** |

---

## Pipeline flow

| stage | scripts | produces |
|---|---|---|
| **1. Setup & QC** | `scripts/01_setup/` | VCF → 290-sample PLINK trio; LD pruning (`--indep-pairwise 1000kb 1 0.2`) → 111,017 SNPs |
| **2. Population structure** | `scripts/02_structure/` | PCA + aIBS kinship on the pruned set; PC-count justification |
| **3. GWAS** | `scripts/03_gwas/` | EMMAX tped; the 24-cell grid (4 traits × {BLUP,BLUE} × {3,5,10} PCs) |
| **4. Diagnostics** | `scripts/04_diagnostics/` | QQ + λ_GC, Manhattans, summary tables for all 24 cells |
| **5. v2 ↔ v3 comparison** | `scripts/05_comparison_v2_v3/` | evidence that the 1 Mb-window LD pruning improved the correction — methods justification |
| **6. Chosen config & loci** | `scripts/06_loci_FINAL/` | **the live result**: SNP-level tables, iterative LD clumping, the 50 kb gap rule, painted Manhattans, paper tables |

Why two population-structure corrections exist: v2 pruned with a **50-SNP** window,
which spans a wildly variable physical distance at this SNP density. v3 uses a
**1 Mb physical** window → 590,462 → 111,017 markers, and the Bonferroni threshold
moved 6.7712 → 6.0454. Stage 5 documents that λ_GC improved in 21 of 24 cells.

---

## Directory map

| path | contents |
|---|---|
| `scripts/` | all code, grouped by stage — see `scripts/README.md` |
| `intermediates/` | PLINK binaries, pruned marker set, PCA/kinship, covariates (11 GB) |
| `results/00_FINAL_BLUP_3PC/` | **the chosen configuration and the loci** |
| `results/emmax_ps/` | raw EMMAX output, 24 cells |
| `results/qq/`, `manhattan/`, `tables/`, `pc_selection/` | 24-cell diagnostics |
| `results/comparison_v2_vs_v3/` | v2 ↔ v3 comparison figures and tables |
| `results/_archive/_archive_v2_win50snp_2026-08-16/` | the frozen **v2** result set. Despite the `_archive` name this is **NOT dead weight** — it is a LIVE INPUT to stage 5, the v2↔v3 comparison that justifies the 1 Mb pruning window. Six scripts in `scripts/05_comparison_v2_v3/` read it by hard-coded path, as does `intermediates/_archive_v2_win50snp_2026-08-16/`. Deleting either would make the methods justification unreproducible. (The v1 set, `_archive_v1_FDR_Bonf005_2026-05-28`, was genuinely dead and was deleted 2026-09-10.) |
| `logs/` | run logs |

---

## Reproduce the locus set

```bash
MAX_GAP=50000 KB=2000 Rscript scripts/06_loci_FINAL/31_loci_clump_iterative.R
Rscript scripts/06_loci_FINAL/33_manhattan_loci_painted.R
Rscript scripts/06_loci_FINAL/36_paper_tables_loci.R
```

## Handoff to the gene search

`results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_for_gene_search.tsv`
carries the 36 search intervals with `chr`,
`start`, `end`, `span_kb`, `lead_SNP`, `class` and an `include_in_gene_search` flag.
That table is the pipeline's output contract.

## Environment notes

- `plink` must be `/usr/local/bin/plink` — the bare name resolves to a broken v0.76.
- Set `TMPDIR=/mnt/data/shahar/.tmp`; never write under `/tmp` or `$HOME`.
- Jobs over ~15 min run in a named GNU `screen` that stays open afterwards.
