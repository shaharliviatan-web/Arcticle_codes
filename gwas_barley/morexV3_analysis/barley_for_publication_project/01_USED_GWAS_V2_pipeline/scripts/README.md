# scripts/

Grouped by pipeline stage. Every script carries a header comment stating what it
does, what it reads, what it writes, and — where a parameter was changed — why.

`plink` must be `/usr/local/bin/plink` (bare `plink` is a broken v0.76 here).
Set `TMPDIR=/mnt/data/shahar/.tmp`.

| dir | stage |
|---|---|
| `01_setup/` | VCF → PLINK, LD pruning |
| `02_structure/` | PCA + kinship, PC-count justification |
| `03_gwas/` | EMMAX tped, 24-cell grid |
| `04_diagnostics/` | QQ/λ, Manhattans, summary tables (all 24 cells) |
| `05_comparison_v2_v3/` | v2 ↔ v3 comparison — methods justification |
| **`06_loci_FINAL/`** | **the live result: chosen config, loci, paper tables** |
| `helpers/`, `_archive_v1/` | shared helpers; retired v1 code |

## 01_setup
| script | does |
|---|---|
| `01_vcf_to_plink.sh` | VCF → canonical 290-sample bed/bim/fam |
| `02_ld_prune.sh` | `--indep-pairwise 1000kb 1 0.2` → 111,017 SNPs. The **1 Mb physical** window replaced a 50-SNP window, which spanned a wildly variable distance at 1,701 SNPs/Mb |

## 02_structure
| script | does |
|---|---|
| `03_pca_kinship.sh` | PCA + EMMAX aIBS kinship on the pruned set; writes the PC covariate files |
| `08c_pc_selection_diagnostics.R` | evidence for choosing 3 PCs |
| `10b_extra_scree_plots.R` | scree / elbow plots |

## 03_gwas
| script | does |
|---|---|
| `04_make_emmax_tped.sh` | recode the full 7.1M-SNP set to tped (independent of pruning) |
| `05_run_emmax_grid.sh` | 24 cells: 4 traits × {BLUP,BLUE} × {3,5,10} PCs |

## 04_diagnostics
`06_qq_lambda.R`, `06c_qq_lambda_v3_grid.R`, `07_manhattan.R`,
`08_summary_tables.R`, `08b_top10_snps.R`, `09_comparison_views.R` — QQ + λ_GC,
Manhattans and summary tables across all 24 cells.

## 05_comparison_v2_v3
`06b`, `07b`, `07c`, `07d`, `07e`, `08d_pc_selection_compare`, `09b`, `09c` —
side-by-side QQ/λ, mirrored Manhattans, contact sheets, SNP concordance.
These read the frozen v2 set in `results/_archive/_archive_v2_win50snp_2026-08-16/`.

## 06_loci_FINAL — the live flow
| order | script | does |
|---|---|---|
| 1 | `11_chosen_config_tables_v3.R` | SNP-level package for BLUP × 3 PCs → `00_snp_level/` |
| 2 | `31_loci_clump_iterative.R` | LD clumping + the **50 kb gap rule** + **iteration until no significant SNP is orphaned** → the 36 loci. Also writes `locus_definition_params.tsv`, the single source of truth for every parameter string used downstream |
| 3 | `36_paper_tables_loci.R` | `Table_*.tsv` — publication tables + the gene-search handoff |
| 4 | `33_manhattan_loci_painted.R` | Manhattans, every member painted in its locus colour |
| — | `37_review_raw_no_gap.R`, `37b_review_protein_with_extra_raw.R` | review-only figures showing the same clumping with **no** gap rule → `02b_REVIEW_raw_no_gap/` |

`31` asserts that all 52 Bonferroni-significant SNPs land inside a locus and fails
loudly otherwise.

Superseded analyses and their scripts were archived under
`results/00_FINAL_BLUP_3PC/_ARCHIVE_superseded_2026-09-08/` and **deleted on 2026-09-10**
during a project-wide cleanup. Nothing live referenced them.
