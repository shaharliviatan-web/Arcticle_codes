# results/

| path | contents |
|---|---|
| **`00_FINAL_BLUP_3PC/`** | **the chosen configuration — start here.** SNP-level tables, the 36 loci, figures, paper tables, and the superseded-analysis archive |
| `emmax_ps/` | raw EMMAX output for all 24 grid cells (`morexV3__<trait>__<BLUP\|BLUE>__pc<3\|5\|10>.ps`) |
| `qq/` | QQ plots, all 24 cells |
| `manhattan/` | Manhattans, all 24 cells × 2 threshold variants |
| `tables/` | 24-cell summary, `lambda_table.tsv`, top-SNP tables |
| `pc_selection/` | PC variance / scree data |
| `comparison_v2_vs_v3/` | v2 ↔ v3 comparison: PC selection, QQ/λ, Manhattans, SNP concordance |
| `_archive/` | frozen earlier result sets — `_archive_v1_FDR_Bonf005_2026-05-28` and `_archive_v2_win50snp_2026-08-16`. Read-only; the comparison scripts read v2 from here |

`emmax_ps/`, `qq/`, `manhattan/`, `tables/`, `pc_selection/` cover **all 24 grid cells**
and exist to justify the configuration choice. Only **BLUP × 3 PCs** was carried
forward; those results are in `00_FINAL_BLUP_3PC/`.
