# 02_raw_phenotype_checks

Descriptive checks of the **raw** phenotypic data, for the supplement only: **not cited in the manuscript** (user
request, 2026-10-01). Their assembled version is the unnumbered supplement
`10_USED_Paper_writing/new_publishing_paper/supplementary/unreferenced/ESM_extra_phenotype_checks.{pdf,xlsx}`; the user
decides later whether any of it joins the numbered Online Resources (it would then need a citation in the text).

Run 2026-10-01, R 4.1.2 (`/usr/bin/Rscript`). **Start here → [`results/tables/results_numbers.txt`](results/tables/results_numbers.txt).**

```bash
bash scripts/run_all.sh     # about a minute; logs in logs/
```

## Data and panel

Input (read-only): `00_THIN_Generate_Plots_For_Publication/all years barley.csv` (raw plant-level values) and
`…/outputs/subsection_1/tables/BLUP_10_traits.csv` (the genotype BLUPs of all analyses).
**Panel = the analysis panel of step 00** (its script 01, section 2): cultivated checks (Morex, Clipper) removed; site 04
(Hachola) removed; block 4 of 2019–20 removed — except in Fig. E2, which shows it, marked "excluded". 290 accessions;
870 / 1,160 / 1,450 plants (2019–20 / 2020–21 / 2021–22), 862 / 1,139 / 1,163 with NIR values (asserted; = M&M).
**Values are the raw values as the analyses use them**, before centering and outlier removal: the 2021–22 starch values
carry the **+8.85 of step 00** (user decision, ESM0 Q12a, 2026-10-01: the supplement uses the corrected values; **the
correction is not reported** in any supplementary file or in the manuscript; it cancels under within-season centering).
Units as in Potapenko et al. (2026a): none for plant height and grain weight (ESM0 Q11b). The ±3 SD outlier rule of step 00
(values centered within season, pooled over seasons, kept within mean ± 3 SD) is **shown, not applied**: on the raw
scale its limits are season mean ± 3 SD of the pooled centered values (T1; dashed lines in Fig. E1).
*Decided by the user (ESM0 Q12, 2026-10-01): (a) starch 2021–22 with the +8.85, not reported; (b) both correlation readings; (c) raw and BLUP distributions.*

## Outputs (TAG figure spec: 174 mm wide, Arial-metric 9 pt, panels a, b, …, no titles in the image; PDF + PNG 600 dpi)

| user's item | figure | table | script |
|---|---|---|---|
| 1a boxplots, the four nutritional traits by season | `Fig_E1_nutritional_traits_by_season` (n per season; ±3 SD limits dashed) | `T1_summary_by_season.tsv` (all ten traits: n plants, n accessions, mean, SD, CV, min, quartiles, max, ±3 SD limits, values beyond) | `01_boxplots.R` |
| 1b boxplots by block within each season | `Fig_E2_nutritional_traits_by_block` (block 4 of 2019–20 shaded, excluded) | `T2_summary_by_block.tsv` | `01_boxplots.R` |
| 2a, 2b correlations between seasons, per trait (each season and the all-season mean) | `Fig_E3_correlations_between_seasons` | `T3_correlations_between_seasons.tsv` | `02_correlations.R` |
| 2a, 2b correlations between the traits, all seasons and per season | `Fig_E4_correlations_between_traits` | `T4_correlations_between_traits.tsv` | `02_correlations.R` |
| 2c correlations between blocks, per season | **Fig. E5 = check 01's** `Fig_interblock_correlation_heatmaps.png` (not recomputed) | check 01's `interblock_correlations_all.tsv` | `../01_block4_…` |
| 3 distributions: the four nutritional and the six morphological and phenological traits of Results Ch. 1 | `Fig_E6_raw_distributions_by_season` (raw, one line per season; x = the ±3 SD window, values outside counted) | (T1) | `03_distributions.R` |
| 3 distributions of the genotype BLUPs | `Fig_E7_BLUP_distributions` | `T5_BLUP_summary.tsv` | `03_distributions.R` |
| key numbers | — | `results_numbers.txt` | `04_results_numbers.py` |

Correlations: Pearson, two-sided `cor.test`, complete accessions per pair, no correction for multiple testing
(descriptive). Accession means: per season, the mean of the raw values over blocks; over all seasons, the mean of the
season-centered values (so that season differences and missing seasons do not shift an accession).
**Two readings of the user's items 2a/2b were open (ESM0 Q12b); both are computed** (between seasons, Fig. E3; between
traits, Fig. E4).

## Findings worth knowing (for the user; none enters the manuscript text)

- **Starch 2021–22 has negative values (internal; not for the supplement text):** as recorded, 25 of 1,163 plants below
  0 % (minimum −21.6 %); with the +8.85, 4 stay below 0 (minimum −12.8 %), and the ±3 SD rule removes only the lowest.
  NIR predictions outside the physical range: raise with N. Pintel together with TODO 213.
- Season means (with the correction) differ: starch 33.5 / 43.0 / 29.6 %; within-season centering absorbs this.
- Raw extremes outside the ±3 SD limits: e.g. 35 tiller counts, 12 spike lengths, 5 plant heights (up to 664 cm).
- Between-season correlations of accession means: starch 0.66–0.77, protein 0.39–0.46, fiber 0.41–0.54, β-glucan 0.26–0.50.
- Trait correlations from raw means agree with the BLUP-based Fig. 1c (starch–fiber −0.77 vs −0.78; starch–β-glucan
  −0.54 vs −0.53; β-glucan–fiber 0.38 vs 0.35; protein–starch −0.13 vs −0.12).

## History

- 2026-10-01: created (user request, supplementary session; questions in `10_USED_Paper_writing/new_publishing_paper/build/proposed/ESM0_QUESTIONS.md`, Q12).
