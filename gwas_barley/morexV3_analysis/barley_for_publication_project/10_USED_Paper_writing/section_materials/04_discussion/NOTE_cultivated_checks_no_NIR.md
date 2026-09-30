# Cultivated checks carry no NIR data — note for the Discussion

Checked 2026-09-24. **Affects one sentence in the Discussion; read before writing it.**

## What was checked

Morex and Clipper were grown as cultivated checks in the same common garden — 17 rows in
`00_THIN_Generate_Plots_For_Publication/all years barley.csv`:

| label | season | rows |
|---|---|---|
| `Morex`, `Clipper` | 2020 | 4 each |
| `MOREX_01`–`MOREX_05` | 2021 | 5 |
| `Clipper_01`–`Clipper_05` (no `_04`) | 2021 | 4 |

**Every one of them is `NA` at all four NIR columns** (`ProteinAsis.`, `StarchAsis.`,
`BetaglucansAsis.`, `FiberAsis.`). The checks were phenotyped for the agro-morphological traits but
not for grain composition.

## What follows

- **We cannot compare wild and cultivated grain composition from this experiment.** The
  Introduction's motivating claim — wild barley has roughly 50% more protein, 38% less starch and
  65% more fiber (Friedman & Atsmon 1988) — must be carried by that citation and must never be
  asserted from our own results.
- The elite-cultivar comparison in Results Ch. 3 and Ch. 4 is **genotypic only**. It says which
  haplotypes five malting cultivars carry. It says nothing about their grain composition.

## Open — worth one check before the Discussion is written

The NIR calibration was developed on the cultivated package, so cultivated grain was measured at
some point. It is possible that check values exist outside this table:

- [ ] the raw NIR exports or the calibration set, rather than the merged phenotype table
- [ ] the collaborators' records (the checks were grown for a reason)
- [ ] N. Pintel's thesis, p. 22, which reports summary statistics with Morex and Clipper as
      cultivated checks — **but that is MorexV2-era and two seasons only; it is not citable
      (unpublished) and its numbers must not be quoted.** It does indicate the values existed.

**If cultivated-check values are found**, the wild-vs-cultivated statement can be made from this
experiment rather than from a 1988 citation, which would be considerably stronger and is worth the
search. **Until then, insert a TODO comment at that sentence in the Discussion draft** recording the
question, so it is visible in the Word file rather than lost here.

<!-- src: checked 2026-09-24 against 00_THIN_Generate_Plots_For_Publication/all years barley.csv, columns 35-38, all 17 check rows; wild rows n = 3,527 with values -->
