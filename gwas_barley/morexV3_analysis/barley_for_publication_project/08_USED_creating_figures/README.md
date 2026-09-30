# 08_USED_creating_figures — manuscript figures (TAG)

Renamed 2026-09-22 from `08_USED_figures**NOT_UPDATED**`. Builds the assembled figures of the TAG
manuscript. Figure specs: `10_USED_Paper_writing/TAG_requirements.md`.

| figure | status | built by | output |
|---|---|---|---|
| **Fig. 1 — variance partition, reaction norms** | **rebuilt to TAG spec and adopted 2026-09-27** | `make_figure_1.R` | `Figure_1/Fig1.png` (docx build), `Figure_1/Fig1.tif` (submission) |
| **Fig. 2 — ecology, trait correlations** | **rebuilt to TAG spec and adopted 2026-09-27** | `make_figure_2.R` | `Figure_2/Fig2.png` (docx build), `Figure_2/Fig2.tif` (submission) |
| **Fig. 3 — GWAS Manhattan + QQ** | **rebuilt and approved 2026-09-22; re-rendered 2026-09-27** (true 9 pt) | `make_figure_3.R` | `Figure_3/Fig3.png` (docx build), `Figure_3/Fig3.tif` (submission) |
| **Fig. 4 — haplotype structure + elite cultivars** | **rebuilt and approved 2026-09-23; re-rendered 2026-09-24** (true 9 pt, ties → REF) | `make_figure_4.R` | `Figure_4/Fig4.png` (docx build), `Figure_4/Fig4.tif` (submission) |
| **Fig. 5 — GDSL esterase at the shared 7H signal** (Results ch. 4) | **built and approved 2026-09-24** (status corrected here 2026-09-27; Ch. 4 approved to stay in the paper 2026-09-25) | `make_figure_5.R` | `Figure_5/Fig5.png` (docx build), `Figure_5/Fig5.tif` (submission) |

**One script per figure.** Each `make_figure_N.R` draws its figure at final size — the only way to
hold TAG's 8–12 pt lettering and 600 dpi — and writes `FigN.png` (docx build) + `FigN.tif` (LZW,
submission), 174 mm wide, true 9 pt, RGB.

*History:* the mini paper's figures were assembled by `assemble_figures.py`, which rescaled panel
PDFs onto a 160 mm canvas at 300 dpi. It was **deleted 2026-09-27** (user approval) once every
figure had its own script: targets 1 and 2 went with the files they produced, target 4 was the
obsolete mini-paper gene figure (Pho, PHT, AP2/ERF, BAHD — three of them not in the current
candidate set; its outputs were deleted 2026-09-23), and target 3 only called `make_figure_3.R`.

---

## Figs. 1 and 2 — `make_figure_1.R`, `make_figure_2.R` (rebuilt 2026-09-27, awaiting approval)

```bash
Rscript make_figure_1.R        # seconds
Rscript make_figure_2.R        # seconds
```

**Why they were rebuilt.** Both were still the mini paper's assemblies: `assemble_figures.py` rescaled
panel PDFs from step 00 onto a 160 mm canvas at 300 dpi. That cannot meet TAG — the width is neither
84 nor 174 mm, plots with lettering need ≥ 600 dpi, there was no TIFF, the stamped panel letters were
upper case at ~15 pt (TAG: 8–12 pt, lower case), and each panel's lettering ended up at whatever its
own rescale factor produced. The old files also clipped the top y-axis label of two panels behind the
stamped letter.

**What did NOT change (user instruction 2026-09-27).** The analysis is untouched: step 00 is neither
re-run nor modified, and the panels are redrawn from its saved output tables with the same data,
geoms, palettes, factor orders, axis labels, legend wording and significance symbols as
`01_phenotypic_analysis_no_GxE_v3_streamlined.R` (A4, A12 A/B), `02_phenotypic_analysis_GxE_v2_streamlined.R`
(B1, B2) and `03b_trait_environment_correlations.R` (C2 all-32 heatmap). Only the rendering differs.

**Fig. 1** (174 × 74 mm): **a** variance partitioning of 8 traits, **b** per-genotype reaction norms
over the three seasons for the 4 nutritional traits, four highlighted accessions per trait.
**Fig. 2** (174 × 230 mm, TAG limit 234): **a** site boxplots by ecological region, **d** trait × environment
heatmap (top row), **b** nutritional trait correlations, **c** nutritional × agro-morphological
correlations (bottom row). Panel positions and letter-to-content mapping are the mini paper's, so the
letters are deliberately not in reading order — the Ch. 1 caption, which is locked text, names them
that way.

**TAG spec (both).** 174 mm wide, Liberation Sans (Arial-metric) at a true 9 pt for every label
(matching Figs. 3–5), 600 dpi, RGB, PNG (docx build) + TIFF/LZW (submission), panel letters a–d in
bold 9 pt. Drawing at final size removes the step-00 panels' hand-enlarged fonts (axis text 24–30 pt,
in-tile labels 4.8–5.6 mm), which existed only to survive the old down-scaling.

**Rendering-only adjustments, forced by the smaller type** (no value changes): Fig. 2c shows 3 y-axis
breaks instead of ggplot's default and reserves 50 % headroom for the significance stars, which at
9 pt would otherwise collide with the facet strip.

**Star clipping fixed 2026-09-28** (user: the stars over Grain weight in the Starch facet of **c**, and
others, were cut off). The stars are drawn a fixed number of text lines above their point, so in a
short panel that offset ate more of the y-range than the headroom allowed. Fixed by making the figure
taller rather than any panel shorter: 208 → **230 mm** (TAG allows 234), with the extra height going to
rows 1 and 2 — panels **a**, **b** and **c** all gained height — plus the star offset trimmed
(vjust −0.9 → −0.7) and the headroom raised to 50 %. No value, colour or position changed.

**Built-in checks** (the scripts stop if any fails): Fig. 1 — 8 traits present, variance components
sum to 100 % per trait, 16 highlighted accessions over 3 seasons, 290 genotypes. Fig. 2 — 29 sites,
all region and trait levels matched, 6 nutritional pairs, 16 nutritional × morphological pairs,
32 trait × environment tests.

**Inputs** (read-only, `00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/`):
`tables/Variance_components_GxE.csv`, `tables/gxe_long_data.csv`,
`tables/B2_reaction_norm_highlights_selected.csv`, `tables/A4_site_boxplot_data.csv`,
`tables/A12_nutri_pairs_uncorrectedP.csv`, `tables/A12_nutri_morpho_pairs_localFDR.csv`,
`C2_trait_environment_correlations/tables/C2corr_trait_environment_all32.csv`.

**Ch. 1 captions were converted to TAG style on 2026-09-27** (user approval): "Fig. 1a" in the text,
bold `Fig. 1` / `Table 1` with no punctuation after the number, bold lower-case panel letters, no end
punctuation. **Adopted into the manuscript 2026-09-27** (user green light): `Results_Discussion.md` (the working
file since 2026-09-27) and the frozen `03_results.md` now embed `Figure_1/Fig1.png` and
`Figure_2/Fig2.png`. The stray Table-1 image artefact of the docx conversion was removed at the same
time, so the order now reads table → **Table 1** caption → Fig. 1 image → **Fig. 1** caption. A build
check confirms all five figures embed at 600 dpi.

**Superseded files deleted 2026-09-27** (user approval): `Figure_1/Figure_1.{png,pdf}`,
`Figure_2/Figure_2.{png,pdf}`, and the `assemble_figures.py` targets 1 and 2 with them. Nothing
referenced those files any more — the manuscript embeds `Fig1.png` / `Fig2.png`.

**Panel a of Fig. 2 was widened on 2026-09-27** (user: too dense to read). The mini paper's 2 x 2 block
gave it 87 mm for 29 sites x 4 traits; it now spans the full 174 mm on its own row, with b | c in the
middle row and d full width at the bottom. Every letter still carries the content the Ch. 1 caption
gives it, and the letters now run in reading order.

---

## Fig. 3 — `make_figure_3.R` (approved 2026-09-22)

```bash
Rscript make_figure_3.R        # ~4 min, ~10 GB RAM
```

**Layout.** One row per trait, as in the mini paper: **a** β-glucan, **b** fiber, **c** protein,
**d** starch; Manhattan left (128 mm), QQ right (46 mm). Panel letter + trait name above each row.

**Manhattan.** Every SNP (7,110,996 per trait) one dot, all dots the same size. The members of each
of the 36 LD-clumped loci (`01_.../02_loci_FINAL`, 50 kb gap rule) are painted in one colour per
locus over an alternating-grey background. No lead-SNP enlargement, no special marker for loci left
with only their lead by the gap rule, no title, legend or locus labels; locus identities go to the
Online Resource locus table. Colours only separate neighbouring loci (8-colour CVD-validated palette,
reused along the genome). One line: Bonferroni threshold, −log10p = 6.0454 (α = 0.10 / 111,017
LD-pruned SNPs). "Chromosome" axis title on the bottom row only.

**QQ.** Observed vs expected −log10(p), y = x line, λGC (3 decimals) only. All 50,000 smallest
p-values plus a random 150,000 of the rest (seed 1). "Expected" axis title on the bottom row only.

**TAG spec.** 174 × 183 mm; Liberation Sans (Arial-metric; Arial is not installed) at a **true 9 pt** for all
lettering (re-rendered 2026-09-27: `layout()` had silently set `par(cex = 0.66)`, so the version approved
2026-09-22 as "12 pt" measured ~7.9 pt on the page, below TAG's 8 pt minimum; `par(cex = 1)` is now set
after `layout()`, `PT = 9` matches Fig. 4, and `PT_CEX` was rescaled 0.40 → 0.35 so the SNP dots keep the
approved size. Same data, same colours, same layout); lines 0.75 pt; RGB; 600 dpi. Raster only: 4 × 7.1 M points make a vector file impractical.
The TIFF (LZW) is converted from the PNG, so the two are pixel-identical.

**Built-in checks** (the script stops if any fails): 7,110,996 SNPs in the same order in every trait
file; λGC equals `01_.../results/tables/lambda_table.tsv` (BLUP, 3 PCs); every genome-wide-significant
SNP lies in a painted locus. Last run: 36 loci (11/18/2/5), 52 significant SNPs (20/24/2/6),
λGC 0.985 / 0.975 / 1.013 / 0.985 (β-glucan / fiber / protein / starch).

**Inputs** (read-only, `01_USED_GWAS_V2_pipeline/`): `results/00_FINAL_BLUP_3PC/01_assoc/<trait>.assoc`,
`results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/loci_{summary,members}.tsv`,
`intermediates/morexV3_pruned_for_covs.prune.in`, `results/tables/lambda_table.tsv`.

**Differs from the step-01 working figures** (`02_loci_FINAL/figures/manhattan_*__loci_r05`,
`04_diagnostics/qq/`), which remain the analysis-review versions: those carry titles, locus keys,
enlarged/open-circle leads, and a palette including a grey that blended into the background.

---

## Fig. 4 — `make_figure_4.R` (approved 2026-09-23)

```bash
Rscript make_figure_4.R        # seconds
```

**Layout.** One row per gene, violins left (70 mm) and the aligned genotype barcodes right (104 mm):
**a**/**d** GPAT6 (fiber), **b**/**e** GH17 (fiber), **c**/**f** PHT4;3 (starch) — the three genes
carried forward by step 06, in the order Results ch. 3 discusses them. Panel letters run down the
violin column (a–c) and then down the barcode column (d–f). One genotype legend, under the barcode
column and centred on the tiles.

**Violins.** Trait BLUP by wild haplotype group, boxplot inside, group letter and *n* below. Densities
are computed **untrimmed** (`cut = 2`) and the y-range is taken from the drawn shapes, so no violin is
cut off flat at the extreme observations. Above each panel: the BH *q* and η² of step 04, and the
Wilcoxon-vs-largest-group brackets of step 07 (Holm-corrected within the gene, `****` ≤ 1e-4 … `ns`).

**Barcodes.** One row per wild haplotype group (per-SNP majority consensus over the group's members,
ties drawn as missing), a white gap, then one row per elite cultivar in order of release. Columns are
SNPs in genomic order. **Only the `shared_sites` version is drawn** (decided in
`10_USED_Paper_writing/new_publishing_paper/build/BLUEPRINT.md`): a column exists only where both call
sets have a record with identical REF and ALT, so nothing is assumed anywhere in the figure. Elite
lines are **not assigned to a haplotype group** — crosshap never saw them; the rows are aligned and
the comparison is left to the reader. Colours are step 07's, themselves inherited from step 04's
heatmaps: Reference `#FFFACD`, Alternate `#2F4F4F`, Heterozygous `#C46210`, Missing `grey70`. The
legend lists only the states actually drawn.

**TAG spec.** 174 × 187.7 mm (limit 234); Liberation Sans (Arial-metric) at a **true 9 pt** for all lettering;
lines 0.75 pt; RGB; 600 dpi. No title inside the image — the gene names are panel labels, and the
windows, SNP counts and colour definitions belong to the caption. The TIFF (LZW) is converted from
the PNG, so the two are pixel-identical.

**Built-in checks** (the script stops if any fails): no shared position is allele-swapped or has a
differing REF (re-asserted from `Table_allele_concordance.tsv`, so a stale matrix cannot produce a
silently inverted figure); the group count and group sizes match step 04's `gene_results.tsv`; the
group means are recomputed from the individuals; the wild consensus contains only REF/ALT/missing;
the elite rows are the five configured cultivars in release order. Last run: GPAT6 3 groups / 43
shared SNPs / q 5.26e-06 / η² 0.208, GH17 5 / 33 / 7.47e-05 / 0.222, PHT4;3 4 / 12 / 3.83e-04 / 0.078 (last run 2026-09-24).

**Fixed 2026-09-24 (user decisions).**
- **Text size.** `layout()` with three or more rows silently sets `par(cex = 0.66)` and nothing reset
  it, so the version approved 2026-09-23 as "11 pt" measured **~7.3 pt** on the page, below TAG's
  8 pt minimum; the margin arithmetic (`LINE_IN`) was wrong for the same reason. `par(cex = 1)` is now
  set after every `layout()` call and `PT = 9`: a true 9 pt (TAG 8–12), ~25% larger than the approved look.
- **Consequences of the larger text, fixed in the drawing only (no statistic changed):** brackets are
  stacked one text line apart in inches above the violins (step 07's `y.position` now sets only their
  order; its data-unit spacing let the stars sit on the next bracket); y-axis ticks stop at the top of
  the data; every `(n=..)` label is drawn (`gap.axis = -1`; R had dropped two touching labels in
  panel b); the legend is drawn by `hlegend()` with per-item widths (R 4.1's `legend()` gives every
  item the widest label's width and no longer fitted), centred on the whole barcode column.
- **Consensus ties → REF** (step 07's `consensus_row()`, see its README): turns the two grey GH17 wild
  cells in panel e (group A `5H:462,728,332`, group E `5H:462,729,123`) into REF.
- **Applied to Fig. 3 on 2026-09-27** (user decision: match Fig. 4 at a true 9 pt). See the Fig. 3 section.

**Inputs** (read-only): `07_USED_elite_lines_compariosn_to_wild_lines/intermediates/matrices/<gene>__shared_sites.rds`,
`07_.../results/tables/{Table_pairwise_group_tests,Table_allele_concordance}.tsv`,
`07_.../config/elite_lines.tsv`,
`04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv`.

**Differs from the step-07 working figures** (`07_.../results/figures/shared_sites/`), which remain the
analysis-review versions: those are one file per gene at 11 × 8.5 in with a title and a subtitle inside
the image, and they also exist in the `filled_marked` and `filled_silent` variants.

**A stacked layout was drawn first and rejected** (2026-09-22/23): violins above the barcodes, one
full-width panel per gene, as in the mini paper's Fig. 4. At 174 mm it reached 226 mm of the 234 mm
limit and still left the violins cramped. `draw_stacked()` is kept in the script so it can be
regenerated, but it is no longer rendered.

---

## Fig. 5 — `make_figure_5.R` (approved 2026-09-24)

```bash
Rscript make_figure_5.R        # seconds
```

Results ch. 4: the GDSL esterase/lipase `HORVU.MOREX.r3.7HG0729030` at the shared 7H
fiber/starch signal. Source analysis:
`03_01_7H_branch_Starch_Fiber_shared_signal_explore/` (the branch root since 2026-09-24; it was the subfolder `REPLACEMENT_ANALYSIS_MGmin3_eps0.9`)
(crosshap MGmin = 3, ε = 0.9, gene ± 1 kb). Built on `make_figure_4.R`: same drawing code, palette,
true 9 pt lettering, line widths, row height, 70 | 104 mm columns, shared-sites barcodes, step 07's
consensus rule (ties → REF; this gene has none).

**Differences from Fig. 4 (user decisions 2026-09-24).**
- **Layout:** violins **a** fiber and **b** starch stacked on the left, **one** barcode **c** spanning
  both rows on the right. The two traits give identical groups (the same 212 | 34 accessions; the
  script stops if not), so their barcodes would be identical.
- **Statistic:** raw Kruskal–Wallis *P* with η² above each violin (single pre-specified gene, outside
  the step-04 BH family), not a BH *q*. Brackets: Wilcoxon A vs B, Holm (one comparison).
- **Red triangles** under the **three genome-wide significant SNPs** (7H:573,606,306, 573,606,460,
  573,606,491), legend entry "Genome-wide significant SNP". The fourth haplotype-defining SNP
  (7H:573,606,282; −log10P 2.35 / 3.13) is not marked.
- **Legend on two lines** (genotype states; the triangle): one line does not fit 104 mm at 9 pt. It sits
  **inside column c, directly under the barcode** (moved 2026-09-24, user request), the two centred as
  one block; there is no separate legend row.

**TAG spec.** 174 × 119.4 mm; true 9 pt; lines 0.75 pt; RGB; 600 dpi; TIFF (LZW) converted from the PNG.

**Built-in checks** (the script stops if any fails): no shared position swapped or with a differing
REF; shared-site count equals `Table_site_overlap.tsv`; elite states re-derived from the VCF equal
`Table_elite_genotypes_wide__*.tsv`; group count and sizes equal `mgmin3_gene_results.tsv`; fiber and
starch groupings identical; wild consensus REF/ALT/missing only; the three significant SNPs are shared
sites and members of the haplotype-defining marker group. Last run 2026-09-24: 2 groups (212 | 34),
22 shared SNPs, triangles at columns 20–22, fiber *P* 2.09e-03 / η² 0.035, starch *P* 1.07e-03 / η² 0.040.

