# S0 — questions on Sariel's structural comments (2026-09-30)

> **ANSWERED 2026-09-30 (user):** Q1 triangle, drawn as a **lower** triangle (caption: "lower triangle"); Q2 A,
> first sentence OK; Q3 approved (also update the image links, `10_USED_Paper_writing/CLAUDE.md`, the 08 README,
> the 7H-branch README, the comment in `11_…/scripts/00_config.R` and `_pdf_index.html`; pixel-check Figs. 2–4);
> Q4, Q5 **no reply to Sariel**; Q6 **no** — the screen categories must not appear in any supplementary file. The new
> Fig. 1 is written as Fig1.png + Fig1.tif. Applied as items S3–S5; the draft images linked below
> (`Figure_1_NEW/`) were removed — the final figure is `08_USED_creating_figures/Figure_1/Fig1.png`.

Everything else is going ahead without you (tables built, text items S1–S2 applied). These five need your answer.
Reply with e.g. "Q1 full, Q2 A, Q3 yes, Q4 OK, Q5 a, Q6 no".

---

## Q1. New Fig. 1 layout (draft)

![draft new Fig. 1](../../../../08_USED_creating_figures/Figure_1_NEW/Fig1.png)

174 × 205 mm (TAG limit 234), 600 dpi, true 9 pt. **a | b** exactly as the current Fig. 1 (b keeps its 4 sub-panels and
the four highlighted accessions). **c** = old 2b + 2c as one heatmap: nutritional traits as rows; the 4 nutritional +
4 morphological traits as columns, a dark line between the two blocks. **d** = old 2d. One colour bar for c and d;
c and d columns are aligned.

Rendering changes (no value changed): in-tile numbers 8 pt (old 2b/2d had 6.5–7 pt, below TAG's 8 pt minimum);
black text on every tile (white text on the darkest tiles fails TAG's 4.5:1 contrast; black passes).

**Variant of c:** the draft shows the nutritional block **full** (each pair twice, diagonal blank), so every row reads as
one trait's complete profile, and c lines up column-for-column with d. The alternative, each pair once (old 2b
triangle, 7 columns, not aligned with d):
[Fig1_c_triangle.png](../../../../08_USED_creating_figures/Figure_1_NEW/Fig1_c_triangle.png).
**Recommend: full.**

## Q2. New Fig. 1 caption

Panels a and b are the current caption word for word; c and d are the old Fig. 2 b/c/d captions merged.

> **Fig. 1** Genetic and ecological architecture of grain nutritional traits. **a** Variance partitioning for four
> nutritional and four morphological traits into genetic, season + block, genotype × environment, and residual
> components. **b** Per-genotype reaction norms across the three seasons for the four nutritional traits; gray lines,
> all accessions; colored lines, four highlighted accessions per trait with values in all three seasons: the most
> stable, the strongest increase and the strongest decrease from 2019–20 to 2021–22, and the strongest crossover
> between seasons. **c** Pearson correlations between accession BLUPs of the four nutritional traits (rows) and of the
> four nutritional and four morphological traits (columns). **d** Pearson correlations between site-mean
> nutritional-trait BLUPs and eight environmental variables; cells outlined in black are significant at *P* < 0.05.
> Significance (**c**, **d**), from the two-sided test of each correlation, **without correction for multiple
> testing**: \**P* < 0.05, \*\**P* < 0.01, \*\*\**P* < 0.001

Checked: every star in c (both blocks) and d is the same test — raw two-sided `cor.test` P, same thresholds
(0.05 / 0.01 / 0.001), n = 290 accessions in c, 29 sites in d (0 mismatches over 54 cells). The old 2c stars were
the raw-P column `sig`, not the FDR column. The caption says so in one clause.

- **A (recommended):** as above.
- **B:** drop "without correction for multiple testing" (the M&M note says raw P was "not stressed"; the old caption
  did not say it). The caption would still say one test covers all stars.

The first sentence merges the two old titles ("Genetic architecture of grain nutritional and reference traits" +
"Ecological architecture of grain nutritional traits"). Change it if you prefer.

## Q3. Renaming scripts and outputs (nothing renamed or removed yet)

| now | after | note |
|---|---|---|
| `make_figure_1.R` → `Figure_1/Fig1.{png,tif}` (old Fig. 1) | replaced by `make_figure_1_NEW.R`, renamed `make_figure_1.R`, writing `Figure_1/Fig1.{png,tif}` | old Fig1.tif is in git; `Figure_1_NEW/` removed |
| `make_figure_2.R` → `Figure_2/Fig2.*` (old Fig. 2) | **removed**; its panel a → `make_figure_ESM_site_BLUPs.R` (already written), b–d → new `make_figure_1.R` | old script and Fig2.tif are in git |
| `make_figure_3.R` → `Figure_3/Fig3.*` (Manhattan) | `make_figure_2.R` → `Figure_2/Fig2.*` | image unchanged |
| `make_figure_4.R` → `Figure_4/Fig4.*` (haplotypes) | `make_figure_3.R` → `Figure_3/Fig3.*` | image unchanged |
| `make_figure_5.R` → `Figure_5/Fig5.*` (7H GDSL) | `make_figure_4.R` → `Figure_4/Fig4.*`; `Figure_5/` removed | image unchanged |
| — | `make_figure_ESM_site_BLUPs.R` → `Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf` (+ PNG preview) | old Fig. 2a, **done** |

Done with plain `mv` (no `git mv`, nothing staged). Inside each renamed script only the output folder, file names,
log tags and header change; I then re-run all four and check that Figs. 2–4 are pixel-identical to the files they
replace (Fig. 2, the Manhattan, takes ~4 min / 10 GB). The manuscript image links and the figure src notes change in
the same step. Old PNGs are not in git, but each has a pixel-identical TIFF that is. **Approve?**

## Q4. Reply to Sariel (his comment on TODO 301: "what's that? don't upload tables, put them as supmat.")

> In TAG, supplementary files are called Online Resources (Electronic Supplementary Material). They are submitted
> with the manuscript and cited in the text as "Online Resource 1", "Online Resource 2", etc. So this one is a
> supplementary table (.xlsx); nothing is uploaded elsewhere. The numbers will be fixed once the list of
> supplementary files is final.

OK, or your wording?

## Q5. Where the reply goes

Sariel's comment is **not in the working file** (nor in any file of the project), so there is no comment ID to reply
to.

- **a (recommended):** add his comment, word for word, as comment 326 (author "Sariel Hübner", `parent="301"`) on the
  TODO 301 placeholder, and the reply as comment 327 (`parent="326"`, author "שחר לויתן" like your other comments).
  In Word both appear on "Online Resource 1".
- **b:** only the reply, as comment 326 with `parent="301"`.

## Q6. (optional) One more column in the all-genes Online Resource?

The Results say the relevance screen "kept only those with a documented or plausible role… The three genes with the
clearest such role… are presented" (8 kept: 3 strong + 5 plausible, `05_…/08_USED_annotation_master/TRAIT_CANDIDACY.md`).
The table shows the annotation and "carried forward", but not the screen's category (strong / plausible / unlikely /
no annotation) for the 23 genes. Add it? **Recommend: no** for now — it would publish the screen's internal grading,
which the text does not describe; say "yes" if a reviewer should see which 8 genes passed.

The built table: [ESM_candidate_genes.xlsx](../../supplementary/tables/ESM_candidate_genes.xlsx)
(TSV twin: [ESM_candidate_genes.tsv](../../supplementary/tables/ESM_candidate_genes.tsv)); variance components:
[ESM_variance_components.xlsx](../../supplementary/tables/ESM_variance_components.xlsx).
