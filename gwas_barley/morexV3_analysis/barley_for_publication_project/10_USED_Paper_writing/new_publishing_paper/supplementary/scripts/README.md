# supplementary/scripts — builders of the supplementary tables (Online Resources)

Created 2026-09-30 for S. Hübner's structural comments on the figures and tables (hand-over:
`../../build/STRUCTURE_CHANGES_2026-10.md`). **No tables remain in the main text**; the two former tables are
Online Resources, built here from the project outputs. Nothing is re-computed: every value is copied from a
read-only project file, and each script stops if a check fails.

| script | output (`../tables/`) | replaces | sources (read-only) | checks |
|---|---|---|---|---|
| `make_ESM_variance_components.py` | `ESM_variance_components.xlsx` (+ `.tsv`) | **Table 1** (H² only) | `00_THIN_…/outputs/subsection_1/tables/Variance_components_GxE.csv` (+ `_WIDE.csv`), `H2_GxE.csv` | every % equals Fig. 1a's input and the WIDE table; % sum to 100; H² = V~G~ / total equals `H2_GxE.csv` and the former Table 1 (tillers "< 0.001") |
| `make_ESM_candidate_genes.py` | `ESM_candidate_genes.xlsx` (+ `.tsv`) | **Table 2**, the 23-gene annotation table (formerly "Online Resource 2") and the planned 55-gene table | `03_…/results/tables/candidate_genes.tsv`; `04_…/04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv`, `genes_not_tested.tsv`; `05_…/08_USED_annotation_master/results/tables/Table_significant_genes_{paper,annotated}.tsv`; `06_…/results/tables/Table_2_genes_carried_forward.tsv` | joins on `gene_id` (+ `lead_SNP`), never `locus_id`; 55 = 30 tested + 25 untestable; 23 significant in 12 loci; per-trait counts of Results Ch. 2–3; the 3 carried-forward genes reproduce every cell of the former Table 2; the relevance-screen categories (strong / plausible / unlikely / no annotation) appear nowhere in the output (user rule, 2026-09-30) |

```bash
python3 make_ESM_variance_components.py    # seconds (openpyxl)
python3 make_ESM_candidate_genes.py        # seconds (openpyxl)
```

Each `.xlsx` has the table on sheet 1 and, on sheet 2 ("Notes"), a draft caption, column definitions and the sources.
**Still to add when the Online Resources are numbered** (not part of the 2026-09-30 task): the file names
`ESM_N.xlsx`, and the TAG title block in each file (article title, journal, authors, corresponding author's
affiliation and e-mail).

The supplementary **figure** from the same task (former Fig. 2a) is built with the manuscript figures:
`08_USED_creating_figures/make_figure_ESM_site_BLUPs.R` → `Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf`.
