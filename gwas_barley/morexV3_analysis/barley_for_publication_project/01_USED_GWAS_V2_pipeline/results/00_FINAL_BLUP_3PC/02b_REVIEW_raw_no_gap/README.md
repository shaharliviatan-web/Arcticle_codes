# 02b_REVIEW_raw_no_gap — review only, NOT the analysis set

The same clumping as the live definition, with **no contiguity rule applied**. Each
locus here is the full min→max span of every SNP `--clump` absorbed, however large the
internal gaps.

**The analysis set is `../02_loci_FINAL/`. Nothing here feeds the gene search.**

## What this is

| parameter | value |
|---|---|
| lead SNPs | `--clump-p1 9.008e-07` |
| members | `--clump-p2 1` — LD alone, no p-value condition |
| LD | `--clump-r2 0.5` |
| reach | `--clump-kb 2000` (±2 Mb, 4 Mb max) |
| **contiguity** | **none — no severing** |

Clumps are read directly from `../02_loci_FINAL/plink/pass1/*.clumped`, so this view
and the analysis set come from exactly the same clumping run.

No iteration is needed: without severing, every SNP `--clump` touched stays in its
clump, so no significant SNP can be orphaned.

## Result vs the analysis set

| | raw, no gap rule | with the 60 kb gap rule |
|---|---|---|
| clumps / loci | 31 | 35 |
| median span | **2,049.3 kb** | 114.2 kb |
| mean span | **2,087.4 kb** | 336.8 kb |
| total span | **64.7 Mb** | 14.7 Mb |

Per trait (raw): beta-glucan 10 clumps, median 1,420.5 kb, max 3,822.7 ·
fiber 14, median 2,529.0, max 3,990.8 · protein 2, median 2,605.8, max 3,185.8 ·
starch 5, median 1,019.1, max 2,049.3.

Including the two sub-threshold protein peaks (raw): `3H:173351655` spans 1,996.5 kb
with 1,147 members and a largest internal gap of 167.1 kb; `3H:198076308` spans
2,008.1 kb with 424 members and a 70.9 kb largest gap. Both are close to the ±1 Mb
half-reach on each side, i.e. bounded by the search radius rather than by the data.

Several clumps sit within a few kb of the 4 Mb ceiling, i.e. their span is set by the
`--clump-kb` reach rather than by anything in the data. That is the behaviour the
contiguity rule exists to remove.

## Files

| file | contents |
|---|---|
| `tables/loci_raw_summary.tsv` | per clump: lead, coordinates, span, member count, largest internal gap |
| `figures/manhattan_<trait>__raw_nogap.{png,pdf}` | Manhattans, every member painted in its clump colour; legend gives span coordinates, span kb, member count and lead position |
| `figures/manhattan_protein__raw_nogap_with_extra.{png,pdf}` | protein, additionally including the 2 sub-threshold 3H peaks as their own index SNPs (diamonds, green / amber) |
| `tables/protein_with_extra_raw_summary.tsv` | the 4 protein entries with the extras, raw spans |

`max_internal_gap_kb` in the table is the useful column: it shows how far apart the
absorbed members actually are inside each span.

## Reproduce

```bash
Rscript ../../../scripts/06_loci_FINAL/37_review_raw_no_gap.R
Rscript ../../../scripts/06_loci_FINAL/37b_review_protein_with_extra_raw.R
```
