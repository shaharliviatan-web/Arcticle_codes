# NIR calibration — source files of Online Resource 2 (sheets "NIR calibration", "Calibration correlations")

Copied unchanged on 2026-10-01 from N. Pintel's folder, uploaded by the user as
`10_USED_Paper_writing/section_materials/02_materials_methods/NIR_calibration_files/NIR.zip` (paths inside the zip below).
Read by `../../scripts/build_ESM_xlsx.py` (ESM 2). Findings: `../../../build/proposed/ESM0_QUESTIONS.md`, Q17.

| file | zip path | content |
|---|---|---|
| `kits.results.xlsx` | `NIR/NIR2020.vs.NIRcali/` | lab values of the 30 calibration accessions (one per site): protein (Kjeldahl), starch and β-glucan (Megazyme), % |
| `calibration.scans.xlsx` | same | NIR predictions of the same grain: mean of 5 DA 7250 scans, 7 Jan 2020, factory cultivated-barley calibrations; its β-glucan column is a copy of the kit values (no β-glucan calibration then), so it is not used |
| `season2020.xlsx` | same | the same accessions in the 2019–20 season ("Season 2020" in the thesis): one plant each (block 1), identical to `00_THIN_…/all years barley.csv` (checked) |
| `kits.cali.2020.R` | same | N. Pintel's script that computed thesis Table 1 from these three files (for the record) |
| `DA7250_Extensive Report_…_avg.xlsx` | `NIR/NIR device data/` | the 150 calibration scans (30 × 5): predictions, calibration file names (`Barley_*_Asis_7July2011_Micro.41M`, bias 0, slope 1); used to check the 5-scan means |
| `Results for NIR calibration.xlsx` | `NIR/Calibration/` | lab subsamples: β-glucan and starch a–d (4; starch of one accession 2), protein a–c (3) |
| `PROTEIN CALC NIR.xlsx` | `NIR/Calibration/` | Kjeldahl calculation: 0.1 N HCl titrant, N × 6.25, subsample weights 0.08–0.26 g |
