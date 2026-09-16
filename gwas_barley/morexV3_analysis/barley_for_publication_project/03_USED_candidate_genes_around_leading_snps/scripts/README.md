# scripts/

Five scripts, run in order by `run_all.sh` (~10 s total). Every script carries a
header block stating what it does, what it reads, what it writes, and why any
non-obvious decision was made.

**No parameter is hard-coded in a script.** Everything lives in
[`../config/params.sh`](../config/params.sh). Bash sources it; R reads the same
file through `_load_params.R`, which parses the `export KEY=VALUE` lines. That way
a value cannot drift between the two languages.

```bash
bash run_all.sh                  # everything
Rscript 01_build_intervals.R     # or a single step, in order
```

| script | reads | writes |
|---|---|---|
| `_load_params.R` | `config/params.sh` | *(sourced, not run)* |
| `01_build_intervals.R` | step-01 handoff table | `inputs/loci_handoff_snapshot.tsv`, `intermediates/loci.tsv`, `loci_intervals.bed` |
| `02_extract_genes.sh` | that BED + Morex V3 GFF3 | `intermediates/genes_1to7H.gff`, `intersect_raw.tsv`, `genes_per_locus.tsv` |
| `03_build_tables.R` | those intermediates | **all 8 output tables** in `results/tables/` |
| `04_flank_sensitivity.R` | `loci.tsv`, `genes_1to7H.gff` | `flank_sensitivity.tsv`, `flank_sensitivity_by_trait.tsv` |
| `05_check_publication_genes.R` | step 04's `final_genes.tsv` | `publication_gene_recovery.tsv` |

## What each one is for

**`01_build_intervals.R`** — the thin adapter. It does *not* define loci; step 01
already did that by LD clumping. It freezes a stamped copy of the input, applies
`FLANK_BP`, and writes a BED.
*Note:* step 01 names the two sub-threshold peaks by their SNP ID (`3H:173351655`).
A `:` in an ID is unsafe once it reaches filenames and plot titles in step 04, so
they are renamed `protein_X01` / `protein_X02`; the original is kept in
`source_search_id`. The script asserts that no locus ID contains an unsafe character.

**`02_extract_genes.sh`** — the lookup. Two intersects, deliberately:
`-wa -wb` gives the gene rows but **silently drops intervals that hit nothing**,
while `-c` reports a 0 for them. Under `FLANK_BP=0`, 18 loci hit nothing, and that
is a reportable result — hence both.
*Coordinates:* the BED is 0-based half-open, the GFF is 1-based inclusive. bedtools
infers this from the `.gff` extension and converts internally. **Do not rename
`genes_1to7H.gff`** to `.txt` or `.bed` — that silently shifts every coordinate by one.

**`03_build_tables.R`** — everything paper-facing. It also *asserts the step-04
contract* (the required column names) and stops if a future edit breaks it, so a
break surfaces here rather than three steps downstream.

**`04_flank_sensitivity.R`** — changes nothing; re-runs the intersection at
0/25/50/100/200 kb so the `FLANK_BP` choice is defensible with numbers, and so its
cost stays visible instead of hidden.

**`05_check_publication_genes.R`** — a regression test. Step 04 has a curated list
of the genes the paper is actually built on, chosen under the *old* search rule.
This asks every run whether the current rule still recovers them, and if not, by
how many bp it misses. **It currently reports 1 of 4 recovered** — see the warning
in the top-level README.

## To change the search rule

Edit `FLANK_BP` in `../config/params.sh`, then `bash run_all.sh`. Nothing else.
