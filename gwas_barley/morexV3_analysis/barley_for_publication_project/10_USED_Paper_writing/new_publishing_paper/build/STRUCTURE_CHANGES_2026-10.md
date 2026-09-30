# Structure changes of 2026-09-30 — hand-over for the text-review session

S. Hübner commented on the figures and tables ("supmat" = supplementary materials). The user settled the
decisions, and this session applied them. Items **S0–S5** are logged at the end of
`build/REVIEW_2026-09-27_Results_Discussion.md`, with questions and answers in `build/proposed/S0_QUESTIONS.md`
and edits in `build/proposed/S*_edits.json`. The working file is `Methods_Results_Discussion_Conclusions.md`,
and its reading copy was regenerated. Nothing was committed.

## 1. Old → new map

| old | new | where it lives now |
|---|---|---|
| Fig. 1a (variance partitioning) | **Fig. 1a** | `08_USED_creating_figures/make_figure_1.R` → `Figure_1/Fig1.{png,tif}` |
| Fig. 1b (reaction norms, 4 sub-panels, 4 highlighted accessions) | **Fig. 1b** (same caption text) | same |
| Fig. 2a (site BLUPs by region) | **Online Resource** (PDF), cited under TODOs 321, 421, 422 | `make_figure_ESM_site_BLUPs.R` → `Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf` |
| Fig. 2b (nutritional × nutritional) + Fig. 2c (nutritional × morphological) | **Fig. 1c**, one heatmap: nutritional traits as rows; lower triangle (each pair once; columns Protein, Starch, β-glucan), then flowering time, plant height, grain weight, grain number | `make_figure_1.R` |
| Fig. 2d (site means × 8 environmental variables) | **Fig. 1d** | `make_figure_1.R` |
| Fig. 3 (Manhattan + QQ) | **Fig. 2** (pixel-identical) | `make_figure_2.R` → `Figure_2/Fig2.*` |
| Fig. 4 (haplotypes + elite barcodes) | **Fig. 3** (pixel-identical) | `make_figure_3.R` → `Figure_3/Fig3.*` |
| Fig. 5 (7H GDSL gene) | **Fig. 4** (pixel-identical) | `make_figure_4.R` → `Figure_4/Fig4.*` |
| Table 1 (H²) | **Online Resource**: variance components + H² | `supplementary/tables/ESM_variance_components.xlsx` |
| Table 2 (the 3 genes) | **Online Resource**: all 55 candidate genes | `supplementary/tables/ESM_candidate_genes.xlsx` |

**No tables remain in the main text.** In the old → new reference changes, panel letters were mapped as follows:
"2b" and "2c" → "1c", "2d" → "1d", "2a" → "Online Resource N", "Fig. 3x" → "Fig. 2x", "Fig. 4x" → "Fig. 3x", "Fig. 5x"
→ "Fig. 4x", and "Figs. 4, 5" → "Figs. 3, 4". Two Discussion sentences cited the whole old "(Fig. 2)". Each now
cites the panels the sentence rests on (a judgement call; see §3, S4). **Hidden `<!-- -->` notes written before
2026-09-30 keep the old numbers**, and a note at the top of the working file gives this map.

## 2. New Online Resources (numbering still open; placeholders "Online Resource N")

| Online Resource | file (built) | built by | sources (read-only) | cited under TODO |
|---|---|---|---|---|
| Variance components and H² of the 8 traits: V~G~, V~E~ (season & block), V~G×E~, V~R~, each as variance and % of total, plus total variance and H² (replaces Table 1) | `supplementary/tables/ESM_variance_components.xlsx` (+ `.tsv`; sheet 2 has the draft caption) | `supplementary/scripts/make_ESM_variance_components.py` | `00_THIN_…/outputs/subsection_1/tables/Variance_components_GxE.csv` (+ `_WIDE.csv`), `H2_GxE.csv` | M&M **238**; Results **320** (Ch. 1 ¶1), **322** (Ch. 2 ¶2) |
| All 55 candidate genes, one row each, nothing filtered: trait, locus, lead SNP, position, distance, MorexV3 description, SNPs in window, testable + reason, accessions grouped, groups, group sizes, KW *P*, BH *q*, η², SD difference, significant; for the 23 significant genes: annotation call and source, the call of each source, best Swiss-Prot hit + identity / coverage + accepted, Pfam, best nr hit; carried forward. It replaces Table 2, the 23-gene annotation table (formerly "Online Resource 2") and the planned 55-gene table | `supplementary/tables/ESM_candidate_genes.xlsx` (+ `.tsv`; sheet 2 has the caption and column definitions) | `supplementary/scripts/make_ESM_candidate_genes.py` | `03_…/candidate_genes.tsv`; `04_…/Stats/gene_results.tsv`, `genes_not_tested.tsv`; `05_…/08_USED_annotation_master/results/tables/Table_significant_genes_{paper,annotated}.tsv`; `06_…/results/tables/Table_2_genes_carried_forward.tsv` | M&M **222, 224, 227**; Results **303** (Ch. 3 ¶2), **323–325** (GPAT6, GH17, PHT4;3 paragraphs) |
| Site figure: centered trait BLUPs across the 29 sites by ecological region (old Fig. 2a) | `08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf` (+ PNG preview) | `08_USED_creating_figures/make_figure_ESM_site_BLUPs.R` | `00_THIN_…/tables/A4_site_boxplot_data.csv` | Results **321** (Ch. 1 ¶2); Discussion **421** (§1), **422** (§2 ¶1) |

Checks built into the scripts (they stop on failure):
- **Variance components:** every % matches Fig. 1a's input, and every H² matches the former Table 1 (tillers "< 0.001").
- **Candidate genes:** the three carried-forward rows reproduce the former Table 2 cell by cell. They also matched the transposed block in the manuscript. Joins use `gene_id`/`lead_SNP`, never `locus_id`, and the per-trait counts match Results Ch. 2–3.
- **Screen categories:** strong / plausible / unlikely / no annotation must not appear in any supplementary file (user rule, S0 Q6). The gene script checks this.

All are registered in `supplementary/supplementary.md`.

## 3. Text passages changed (all through `review_tool.py`)

- **S1: Table 1 out.**
  - Removed the table and its caption; the hidden note keeps both.
  - "Table 1" → "Online Resource N" in three places: M&M Phenotypic ¶3 now reads "(Fig. 1a; Online Resource N)", plus Results Ch. 1 ¶1 and Ch. 2 ¶2.
- **S2: Table 2 out.**
  - Removed the caption, table and footnotes; a hidden note keeps them.
  - Results Ch. 3 ¶2: "(Table 2; Online Resource 2)" → "(Online Resource N)". TODO 303 was rewritten; the user's comment 318 stays.
  - "Table 2" → "Online Resource N" in the GPAT6, GH17 and PHT4;3 paragraphs.
  - The texts of TODOs 222, 224 and 227 were updated. No wording or numbers were added.
- **S3: Figs. 3, 4, 5 → 2, 3, 4.** 25 references and captions changed: M&M ×6, Results Ch. 2–4, Discussion §3 and §5. This includes "figure 3" → "figure 2" in the user's Word comment 311.
- **S4: old Fig. 2 references.**
  - M&M ¶4: "(Fig. 2b) … (Fig. 2c)" became one "(Fig. 1c)" at the end of the sentence.
  - Old 2d → 1d in the M&M, Results and Discussion. Old 2b/2c → 1c in Results Ch. 1 ¶2, Ch. 2 ¶4 and Ch. 4 ¶3.
  - Old 2a → Online Resource N (TODO 321).
  - Discussion §1 "(Fig. 2)" → "(Fig. 1d; Online Resource N)" and §2 ¶1 "(Fig. 2)" → "(Fig. 1c, d; Online Resource N)" (TODOs 421, 422).
- **S5: new Fig. 1 caption.**
  - Title: "Genetic and ecological architecture of grain nutritional traits". Panels a and b keep their text word for word.
  - "**c** Pearson correlations of the nutritional-trait BLUPs (rows) with one another (lower triangle) and with four morphological-trait BLUPs".
  - Panel d keeps the old 2d wording. New significance line: "from the two-sided test of each correlation, without correction for multiple testing". Checked: all stars in c and d come from the same raw P and thresholds.
  - Removed the old Fig. 2 image and caption (a hidden note keeps the caption), and updated the image links and src notes. Added the numbering note at the top of the file.
- **No reply to Sariel** (user, S0 Q4/Q5).

The M&M wording was touched only in figure and table references and in the Online Resource citations. docx build
check: 75 Word comments (66 + 9 new TODOs: 238, 320–325, 421, 422), 4 figures, 0 tables. Next free IDs: **239**
(M&M), **326** (Results), **423** (Discussion).

## 4. Files and scripts created or changed

- **Created:**
  - Scripts: `08_USED_creating_figures/make_figure_ESM_site_BLUPs.R`; `supplementary/scripts/make_ESM_variance_components.py`, `make_ESM_candidate_genes.py`, `README.md`.
  - Outputs: `Figure_1/Fig1.{png,tif}` (new content); `Figure_ESM_site_BLUPs/`; `supplementary/tables/*.{xlsx,tsv}`.
  - Review files: `build/proposed/S0_QUESTIONS.md`, `S1–S5_edits.json` + `_PROPOSED.md`, and this file.
- **Rewritten:** `make_figure_1.R` (new Fig. 1; the header records its history).
- **Renamed with plain `mv`** (nothing staged): `make_figure_{3,4,5}.R` → `make_figure_{2,3,4}.R`, and `Figure_{3,4,5}/Fig{3,4,5}.*` → `Figure_{2,3,4}/Fig{2,3,4}.*`. Only the header, output names and log tags changed; each script was re-run and is pixel-identical to the files it replaces.
- **Removed** (user-approved, S0 Q3):
  - The previous `make_figure_1.R` and `make_figure_2.R`, and `Figure_1/Fig1.*` and `Figure_2/Fig2.*` of old Figs. 1–2. The scripts and TIFFs are in git.
  - `Figure_5/` (now empty) and the draft `Figure_1_NEW/`.
  - Copies of all old images and scripts are in this session's scratchpad (`…/scratchpad/old_figs/`, temporary).
- **Docs and comments updated:**
  - `10_USED_Paper_writing/CLAUDE.md`: module map, What changed, sources of the numbers, Figure sources, TODO IDs, layout, pandoc test counts.
  - `08_USED_creating_figures/README.md`, `supplementary/supplementary.md`, `03_01_7H_branch_…/README.md` (Fig. 5 → 4), `06_USED_genes_selected_to_present/README.md` (Table 2 note + history row).
  - `11_USED_…/scripts/00_config.R` (one comment), `_pdf_index.html` (its four dead `06_USED_figures/Figure_N` entries replaced by the ESM PDF).
  - The review log.

## 5. Noticed outside the scope, not fixed

- **Results Ch. 3 ¶1** says the 25 untestable genes "carried no SNPs, too few SNPs, or SNPs too weakly linked". One of them, `3HG0301250` (47 SNPs), failed with a crosshap internal error. The gene Online Resource now shows that reason, so the reader will see it (the Ch. 3 note records "crosshap software failure not reported" as a decision).
- **User Word comments now answered by the new Online Resources:** 317 (Ch. 3 ¶1, "Add supplementary table of the genes?") and 318 (reply to 303, "Add supplementary table"). Ch. 3 ¶1 itself does not cite the gene Online Resource (the text was kept as is). Comments 314/315 (per-SNP table) and 316 (co-localized loci table) are still open.
- **`_pdf_index.html`** is stale beyond the figures: generated 2026-08-17, 959 of its 1,105 links are dead, and no generator script was found.
- **`07_USED_…/README.md`** history row of 2026-09-24 names "Fig. 4 (`make_figure_4.R`)". It is a dated record and was left as written.
- **Fig. 1:** the columns of c (7) and d (8) are not aligned, a consequence of the triangle (accepted in S0 Q1).

## 6. Left open

- **Online Resource numbering** (another session): numbers in citation order, files `ESM_N.xlsx` / `ESM_N.pdf`, the TAG title block in each file, and a concise caption for each in the manuscript. The draft captions are on sheet 2 of each `.xlsx` and in `supplementary.md`. The placeholders already in the text ("Online Resource 1", "3", "X") are unchanged.
- **Site figure vs sampling-site table (TODO 218):** same sites and regions, so they could be one Online Resource.
- **Variance-components caption:** it says "squared units of each trait". Grain weight and spike length have no stated unit in the project (cf. TODO 209).
- **Sariel's comment on TODO 301:** no reply in the file (user decision). If he raises it again, the answer is that in TAG an Online Resource *is* a supplementary file (ESM), cited as "Online Resource N". Nothing is uploaded elsewhere.
