# 12_additional_checks

This folder holds small stand-alone checks that back specific statements in the
TAG manuscript. Each check has its own subfolder with `scripts/`, `results/`,
`logs/` and a README. Checks only read inputs from other project steps and
never write to them.

| # | check | question | inputs | key output | date run |
|---|---|---|---|---|---|
| 01 | [`01_block4_2019_interblock_correlations`](01_block4_2019_interblock_correlations/README.md) | Is block 4 of 2019–20 less concordant with the other blocks in the four NIR grain traits? (This is the basis for excluding it.) | `00_THIN_Generate_Plots_For_Publication/all years barley.csv` | [`results_numbers.txt`](01_block4_2019_interblock_correlations/results/tables/results_numbers.txt) | 2026-09-29 |
