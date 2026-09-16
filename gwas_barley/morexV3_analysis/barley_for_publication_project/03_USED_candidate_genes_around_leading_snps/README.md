# 03_USED_candidate_genes_around_leading_snps

Annotated genes at the GWAS loci, for the wild-barley grain-quality paper.
290 wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*, Southern Levant),
Morex V3, 4 traits (beta-glucan, fiber, protein, starch).

**Rebuilt 2026-09-08; last re-run 2026-09-10** against step 01's final locus
definition: `--clump-kb 2000` (±2 Mb), max span 4000 kb, **50 kb gap rule**, iterative,
and **no sub-threshold protein peaks** — those were dropped by decision, so protein now
contributes only its 2 genome-wide-significant loci. The retired 200 kb single-linkage version was deleted on 2026-09-10.

**Start here → [`results/tables/results_chapter_numbers.txt`](results/tables/results_chapter_numbers.txt)** — every
number in this step, in prose, ready to paste into the manuscript.

---

## What this step does, in one line

Takes the **36 finished loci** from step 01, intersects their genomic intervals with
the Morex V3 gene annotation, and writes the candidate-gene tables that step 04
(crosshap haplotype analysis) consumes. **Tables only — no plots.**

It does **not** define loci. Step 01 does that, by LD clumping. This step only
converts loci into intervals and looks up what is inside them.

## The result, in one table

| | |
|---|---|
| Loci searched | **36** (all genome-wide significant; no sub-threshold peaks) |
| Search rule | interval = **the step-01 LD locus span**, `FLANK_BP = 0` |
| Total search space | **7.91 Mb** |
| Loci containing a gene | **15** |
| Loci with no annotated gene | **21** |
| **Candidate genes** | **55** (all protein-coding) |
| With a functional description | **3** |

| trait | loci | with genes | empty | genes | annotated |
|---|---|---|---|---|---|
| beta-glucan | 11 | 4 | 7 | 7 | 0 |
| fiber | 18 | 6 | 12 | 32 | 1 |
| protein | 2 | 2 | 0 | 6 | 0 |
| starch | 5 | 3 | 2 | 10 | 2 |

> **Locus IDs are not stable across step-01 re-runs.** When step 01 was rewidened,
> **19 of the 34 previous IDs came to mean a different locus** — new loci were
> inserted earlier in the sort order and pushed the `L##` numbering down (e.g. the
> locus led by `3H:543690146` went from `fiber_L04` to `fiber_L07`). Always join on
> `lead_SNP`, never on `locus_id`, when comparing across versions.

---

## ⚠ Read this before using the output

`FLANK_BP = 0` — no flanking window is added to the LD span. This was a deliberate
choice (see below), and it has two consequences you must know about:

**1. 21 of 36 loci return no gene.** They are listed in
[`results/tables/loci_without_genes.tsv`](results/tables/loci_without_genes.tsv), not
silently dropped.

**2. Three of the four genes the paper's haplotype chapter is built on are NOT in this
candidate list.** Step 04's curated `final_genes.tsv` names Pho, PHT, AP2/ERF and BAHD;
only **PHT** falls inside its LD span under the 50 kb gap rule. (Pho was briefly recovered
under the 60 kb rule, which produced wider spans; tightening back to 50 kb dropped it again.) Checked automatically every run by
`scripts/05_check_publication_genes.R` →
[`results/tables/publication_gene_recovery.tsv`](results/tables/publication_gene_recovery.tsv):

| gene | trait | in list? | nearest locus | gap to locus edge |
|---|---|---|---|---|
| PHT | starch | **yes** | starch_L03 | 0 |
| BAHD | fiber | no | fiber_L15 | 34,865 bp |
| AP2/ERF | beta-glucan | no | betaglucan_L04 | 69,298 bp |
| Pho | starch | no | starch_L03 | 102,704 bp |

Smallest flank that recovers all four: **102,704 bp**.

**If that is not what you want, it is a one-line change:** set `FLANK_BP` in
[`config/params.sh`](config/params.sh) (e.g. `200000`, the genome-wide LD-decay
distance) and re-run `bash scripts/run_all.sh`. The whole step takes ~10 seconds.
[`results/tables/flank_sensitivity.tsv`](results/tables/flank_sensitivity.tsv) shows
exactly what each option returns:

| flank | genes | empty loci | searched |
|---|---|---|---|
| **0 (current)** | **55** | **21/36** | **7.9 Mb** |
| 25 kb | 78 | 14/36 | 9.7 Mb |
| 50 kb | 93 | 10/36 | 11.5 Mb |
| 100 kb | 122 | 7/36 | 15.1 Mb |
| 200 kb | 198 | 3/36 | 22.3 Mb |

---

## Why `FLANK_BP = 0`

Step 01 already defines each locus as a real LD block: PLINK `--clump` with lead
p < 9.008e-07, **members admitted on LD alone** (r² ≥ 0.5, no p-value condition),
±2 Mb reach, then severed at the first internal gap > 50 kb, applied iteratively. That block is the
region genuinely in LD with the lead SNP, so it is the honest search space. Adding
a fixed flank on top would re-import the distance-only assumption that the LD
clumping was introduced to replace.

The cost is the 21 empty loci and the 3 missing curated genes above. The empty loci
are not biologically barren — their median span is **16.0 kb** against **207 kb** for
loci that do contain a gene; the interval is simply too small to hold one. **4 of them
are lead-only loci** (a single SNP, span 0), whose search interval is 1 bp — such a locus
yields a gene only if the SNP falls literally inside a gene body, and none do. The
alternative rule — the one used in the archived v1 and still defensible — is a
**±200 kb** flank, justified by the genome-wide LD decay crossing r² = 0.2 at
188 kb (`02_USED_LD_decay_V2_wholegenome`).

---

## Pipeline

```bash
bash scripts/run_all.sh          # ~10 seconds, end to end
```

| # | script | does |
|---|---|---|
| 1 | `01_build_intervals.R` | step-01 locus table → search intervals; freezes a snapshot of the input |
| 2 | `02_extract_genes.sh` | filters the GFF to genes on 1H–7H, `bedtools intersect` |
| 3 | `03_build_tables.R` | parses annotation, computes distances, writes every output table |
| 4 | `04_flank_sensitivity.R` | what other flank sizes would have returned |
| 5 | `05_check_publication_genes.R` | regression check against step 04's curated gene list |

All parameters live in **[`config/params.sh`](config/params.sh)** and nowhere else.
Bash sources it directly; R reads the same file via `scripts/_load_params.R`, so the
two can never drift.

## Directory map

| path | contents |
|---|---|
| `config/params.sh` | **every parameter**, single source of truth |
| `scripts/` | the five scripts + the shared R param loader — see `scripts/README.md` |
| `inputs/` | frozen snapshot of the step-01 handoff table (provenance) |
| `intermediates/` | intervals BED, filtered GFF, raw intersect, per-locus counts |
| `results/tables/` | **all output** — see `results/README.md` |
| `logs/` | one log per script, plus `run_all.log` |

## Where the input comes from

`01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_for_gene_search.tsv`

36 rows, each with `chr / start / end / lead_SNP / class / include_in_gene_search`.
A copy is frozen into `inputs/loci_handoff_snapshot.tsv` at every run, stamped with
the source path and mtime, so a result stays reproducible even if step 01 is re-run.

**This is the only SNP/locus input.** There is no hand-curated SNP list any more —
the v1 `marginal_snps.tsv` is obsolete, because the v3 Bonferroni threshold (6.0454,
down from 6.7712) already makes most of those SNPs significant.

Because step 01 is re-run independently, `01_build_intervals.R` freezes a stamped copy
of the handoff table into `inputs/loci_handoff_snapshot.tsv` at every run — that
snapshot, not the live step-01 table, records which locus set produced the current
`results/`.

## Handoff to step 04

Step 04 reads two files from `results/tables/` and requires these exact column names:

| file | required columns |
|---|---|
| `candidate_genes.tsv` | `trait, locus_id, lead_SNP, class, gene_id, chr, gene_start, gene_end, strand, dist_to_lead_bp, description` |
| `lead_loci.tsv` | `trait, locus_id, lead_pos, lead_neg_log10p` |

Both filenames and both spellings are load-bearing — `lead_pos` / `lead_neg_log10p`
keep **step 04's** spelling, not step 01's `lead_bp` / `lead_neg_log10_p`. Do not
"tidy" them. `03_build_tables.R` asserts this contract and stops if it breaks.

Verified 2026-09-08: step 04's `00_build_gene_windows.R` consumes the new tables
unchanged. `class` values are `significant_locus` (and, in earlier locus definitions,
`subthreshold_peak`); step 04 only carries `class` through as a label and never filters
on it, so no step-04 edit is needed.

Step 04 was itself rebuilt on 2026-09-08/09 and is current — its live run is
`loci_LDspan_eps06_V4`, built from this step's 55 candidate genes.

## Conventions

- All paths absolute. `bedtools` = `/usr/bin/bedtools` (v2.30.0). R 4.1.2.
- `TMPDIR=/mnt/data/shahar/.tmp`; nothing written under `/tmp` or `$HOME`.
- Scripts write only named files under `intermediates/` and `results/`.
- Genes = GFF3 col3 `gene` only (not mRNA/CDS, which would multiply-count by transcript).
- Chromosomes 1H–7H; unplaced `CAJHDD*` scaffolds excluded (they carry no GWAS SNPs).
