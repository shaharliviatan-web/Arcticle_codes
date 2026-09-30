# 02_loci_FINAL — the locus definition used for the candidate-gene search

**36 loci.** Locked 2026-09-09. This is the set the gene search runs on.

## How a locus is defined

| step | rule |
|---|---|
| 1. **Lead SNPs** | `--clump-p1 9.008e-07` — only Bonferroni-significant SNPs (−log10p ≥ 6.0454) may lead |
| 2. **Members** | `--clump-p2 1` — **membership on LD alone, no p-value condition**. A causal variant need not itself be significant, so filtering members by p-value would bias the window |
| 3. **LD** | `--clump-r2 0.5` |
| 4. **Reach** | `--clump-kb 2000` → ±2 Mb, 4 Mb maximum span |
| 5. **Contiguity** | members sorted by position; the locus is **severed at the first gap > 50 kb** walking outward from the lead in each direction |
| 6. **Iteration** | any significant SNP still outside every locus becomes a candidate lead in the next pass; repeat until none remain |

Every parameter is written to `tables/locus_definition_params.tsv` by the script that
applies them, and all downstream figures and tables read their parameter strings from
that file — so a figure caption can never state a value the run did not use.

### Why step 5 matters

`--clump` reports the min→max position of every member it absorbs, so a single member
sitting across an empty stretch inflates the span to the full clump reach. Severing at
the first > 50 kb gap keeps a locus to the region that is actually contiguous in LD.

### Why step 6 exists

`--clump` is winner-take-all: each SNP joins exactly one clump. A significant SNP that
is not a lead, and that sits beyond a gap from the lead which owns it, is severed and
then belongs to nothing. Without the iteration such a SNP is simply absent from the
tables and the painted Manhattans, even though it passed the genome-wide threshold.

The iteration converged in **3 passes**. Five loci come from recovery passes; the
`pass` column records which. The script **asserts** that all 52 significant SNPs fall
inside a locus and fails loudly if not.

## Results

| trait | loci | median span | max span | from recovery passes |
|---|---|---|---|---|
| beta-glucan | 11 | 57.6 kb | 537.4 kb | 1 |
| fiber | 18 | 74.1 kb | **1,727.8 kb** | 4 |
| protein | 2 | 667.8 kb | 1,290.4 kb | 0 |
| starch | 5 | 169.0 kb | 745.2 kb | 0 |

**All 36:** median **93.5 kb**, mean 219.8 kb, range 0 – 1,727.8 kb. 4 loci are
lead-only (no member survived the gap rule). **Total search space 7.9 Mb.**

The 50 kb gap tolerance sets the total search space, so quote it explicitly in the
methods alongside the span statistics.

Protein carries its **2 Bonferroni-significant loci only**. The two sub-threshold 3H
peaks that were previously carried alongside them are not part of this set; their
earlier outputs and scripts are kept under
`../_ARCHIVE_superseded_2026-09-08/results/protein_subthreshold_extras/`.

## Files

| file | contents |
|---|---|
| `tables/Table_loci_master.tsv` | **the main table** — one row per locus: coordinates, span, members, lead A1/A2/MAF/beta/SE/p, severing diagnostics. **`lead_beta` = effect of `lead_A1` (minor allele)**; sign flipped from EMMAX on 2026-09-22 (see `../README.md`) |
| **`tables/Table_loci_for_gene_search.tsv`** | **the output contract** — 36 rows with `chr`, `start`, `end`, `span_kb`, `lead_SNP`, `class`, and an `include_in_gene_search` flag |
| `tables/Table_loci_members_full.tsv` | every member SNP: position, MAF, p-value, lead flag |
| `tables/Table_loci_per_trait.tsv` | per-trait counts and span distribution |
| `tables/Table_analysis_parameters.tsv` | every constant, one row each — for the methods section |
| `tables/locus_definition_params.tsv` | **single source of truth** for the parameters; read by all downstream scripts |
| `tables/loci_summary.tsv`, `loci_members.tsv` | raw output of script 31 (the `Table_*` files derive from these) |
| `figures/manhattan_<trait>__loci_r05.{png,pdf}` | Manhattans, every member painted in its locus colour; legend carries each locus's span coordinates, span kb, member count and lead position |
| `figures/locus_colour_key__<trait>.tsv` | locus → colour mapping |
| `plink/pass<k>/` | raw PLINK output per iteration pass |

## Reproduce

```bash
cd ../../..                      # 01_USED_GWAS_V2_pipeline
MAX_GAP=50000 KB=2000 Rscript scripts/06_loci_FINAL/31_loci_clump_iterative.R
Rscript scripts/06_loci_FINAL/36_paper_tables_loci.R
Rscript scripts/06_loci_FINAL/33_manhattan_loci_painted.R
```

Loci ~8 min (3 passes), figures ~25 min. Input is `../01_assoc/<trait>.assoc`.

A companion view with the **same clumping but no gap rule** is in
`../02b_REVIEW_raw_no_gap/` — review only, not part of the analysis.
