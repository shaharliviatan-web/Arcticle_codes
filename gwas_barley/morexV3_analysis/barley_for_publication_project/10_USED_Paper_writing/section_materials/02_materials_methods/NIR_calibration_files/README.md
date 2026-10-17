# NIR calibration files — upload here

Created 2026-10-01 for the supplementary session (ESM0 Q5). **Put the Excel sheets on the NIR calibration in this
folder**, as they are (do not rename or edit them); add a line per file to the table below if you can.

They are needed for:
- **Online Resource 2** (`new_publishing_paper/supplementary/ESM_2.xlsx`): a sheet with the validation of the NIR
  calibration (Pearson r of the NIR predictions against the laboratory values), cited in the M&M, NIR ¶2 (TODO 214).
- **TODO 213** (whether one calibration was applied to all three seasons) and **TODO 212** (Kjeldahl protein).

What helps most (any of it):
1. The **laboratory values** of the 30 calibration accessions (one per site): starch and β-glucan (Megazyme kits,
   mean of the subsamples), protein (Kjeldahl); with the accession IDs.
2. The **NIR predictions** for the same accessions: the calibration scans, and the scans of the season(s) compared
   (Pintel thesis Table 1 used the 2019–20 season, called "Season 2020" there).
3. Any **correlation table** already computed (the thesis Table 1 values: kit vs calibration scans, kit vs season scans,
   season scans vs calibration scans), with the **number of accessions** per comparison.
4. Which **calibration** (instrument package, date, version) was applied to which season.

Only what you upload here is used; the values will be copied, not re-estimated, unless you ask otherwise.

| file | what it holds | from whom / date |
|---|---|---|
| `NIR.zip` | N. Pintel's NIR folder: lab results of the 30 calibration accessions (Megazyme, Kjeldahl incl. the calculation sheet, N × 6.25), DA 7250 exports (calibration scans 7 Jan 2020 with the factory calibration names; season scans 2020 and 2021), the paired files behind thesis Table 1 (`NIR2020.vs.NIRcali/`) and Noam's R scripts. Read 2026-10-01; findings in `new_publishing_paper/build/proposed/ESM0_QUESTIONS.md`, Q17 | user, 2026-10-01 |
