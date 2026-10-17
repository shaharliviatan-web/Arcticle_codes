# 12_USED_adittional_checks_for_supplementary

**Renamed 2026-10-01** (user request) from `12_USED_adittional_checks`. Only the path changed: the one
hard-coded path (`CHECK_DIR` in check 01's script) was updated, and a re-run in scratch reproduced every output
byte for byte (only the run date in `results_numbers.txt` differs). References were updated in the manuscript's
hidden notes, `10_USED_Paper_writing/CLAUDE.md`, `new_publishing_paper/supplementary/supplementary.md` and the ESM
review files; dated records (the run log of check 01, `build/proposed/MM9_edits.json`) keep the old name.

This folder holds small stand-alone checks for the **supplementary material (Online Resources)** of the TAG
manuscript: checks that back specific statements in the text, and descriptive checks of the raw phenotypic data
that go to the supplement only (not cited in the text; the user decides later whether they enter the main
supplement). Each check has its own subfolder with `scripts/`, `results/`,
`logs/` and a README. Checks only read inputs from other project steps and
never write to them.

| # | check | question | inputs | key output | date run |
|---|---|---|---|---|---|
| 01 | [`01_block4_2019_interblock_correlations`](01_block4_2019_interblock_correlations/README.md) | Is block 4 of 2019–20 less concordant with the other blocks in the four NIR grain traits? (This is the basis for excluding it.) | `00_THIN_Generate_Plots_For_Publication/all years barley.csv` | [`results_numbers.txt`](01_block4_2019_interblock_correlations/results/tables/results_numbers.txt) | 2026-09-29 |
| 02 | [`02_raw_phenotype_checks`](02_raw_phenotype_checks/README.md) | Descriptive checks of the raw phenotypes, for the supplement only (not cited): boxplots by season and by block, correlations between seasons, traits and blocks (blocks = check 01), distributions of the raw values and of the BLUPs. | `00_THIN_Generate_Plots_For_Publication/all years barley.csv`, `…/BLUP_10_traits.csv` | [`results_numbers.txt`](02_raw_phenotype_checks/results/tables/results_numbers.txt) | 2026-10-01 |
