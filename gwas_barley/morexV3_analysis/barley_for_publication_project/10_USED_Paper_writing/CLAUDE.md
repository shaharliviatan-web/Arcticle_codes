# 10_USED_Paper_writing — rules for writing the TAG manuscript

Wild-barley grain-quality GWAS paper for **Theoretical and Applied Genetics (TAG)**.
290 *H. spontaneum* accessions, **MorexV3**, 4 traits (protein, starch, β-glucan, fiber). Created 2026-09-15.

---

## How we work

1. **One section per session.** Each future session works on one section only. Sections already approved
   are context for writing the others: **read the approved sections in `new_publishing_paper/` before
   starting.**
2. **Section rules files.** A session may get an extra rules file for its section,
   `section_materials/<NN_section>/SECTION_RULES.md`. It **adds to** this file and does not replace it.
   If the two conflict, ask the user.
3. **English only.** All replies and all manuscript text are in English, even when the user writes in Hebrew.
4. **Short replies.** Chat answers are brief, focused and to the point, with no long reports for the user to read.
5. **Review only, unless approved.** Reading and reviewing code, analyses and project files is always fine.
   **Never change them without the user's explicit approval.** Without approval, write only to
   `new_publishing_paper/`. If a number looks wrong, report it. After an approved analysis change, document it
   in the script comments and READMEs (project `CLAUDE.md`).
6. **The user approves each section.** Commit only when asked.
7. **Results and Discussion sessions include brainstorming.** Before drafting, explore the outputs with the
   user to find the results worth writing about. Suggest candidates with their evidence (file and numbers);
   the user decides what goes in.

---

## The study in brief

- **This study uses MorexV3.** The earlier work on the same collection, Noam Pintel's thesis and Evgenii
  Potapenko's papers, used **MorexV2**. Never carry V2 coordinates, gene IDs or SNP counts into this paper.
- **All data for this study** are in `/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/`.
- **Some bioinformatics methods changed** since the mini paper, so some results and conclusions in this
  paper differ from it.
- **The mini paper's analysis folder** (`barley_project_assignment/`) is for comparison. Before writing a
  section, check what changed there, from the bioinformatics methods to the results, conclusions and
  discussion (module map below).

## Locations

| role | path | access |
|---|---|---|
| **New analysis (source of truth)** | `/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/` | read |
| **Old analysis** (behind the old mini paper) | `/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_project_assignment/`: 390 GB, never copy or symlink it | read only |
| **Old mini paper** | `old_mini_paper/md/`: `00_FULL_PAPER.md`, then `01_…`–`08_…` in section order, figures in `media_*/`. Originals in `old_mini_paper/docx/` | read only |
| **Section materials** (extra text cut from the mini paper, per section) | `section_materials/<NN_section>/` | read (user adds) |
| **Reference sources** | `reference_sources/<source>/`: original (pdf/docx) and readable `.md` side by side | read only |
| **Journal (TAG)** | `TAG_requirements.md` (the rules); `journal_TAG/` (saved guideline pages, 2026-09-15) | read |
| **New manuscript** | `new_publishing_paper/` | write |

## Input texts and how far to trust them

**The old mini paper is the strongest text reference.** It is the most recent write-up of this research,
and it ranks below only the project files.
- **Style and level:** much of it was written by **Sariel Hübner** (supervisor). That is the writing level and style to match: voice, sentence length, paragraph flow, density. Rewrite text from every other source (section materials, thesis, papers) into this style.
- **Research process:** it shows how the study was done and argued. Reuse its framing, arguments and wording, but check every topic against the project files first: some analyses changed, and with them the results and conclusions. Keep the style, update the facts.
- **It is condensed.** It was cut to the ~13–14-page limit of a BSc bioinformatics project, so whole parts are shortened or missing (e.g. the genotyping pipeline). The user has much of the cut text and adds it per section (below). Something absent from the mini paper was not necessarily absent from the study: check the project files, or ask the user.

**Section materials.** When a section is written, the user adds that section's extra material (text cut
from the mini paper, collaborators' text, notes) to `section_materials/<NN_section>/`, numbered like the
manuscript files, and sends it to that session.
- Before writing a section, **read that section's folder.**
- Check its factual claims against the project files before using them.

**Reference sources** (`reference_sources/`) are earlier work on the same collection. They are accurate but
**MorexV2**-based, only partly relevant, and each is good for different things:

| source | files | use it for | do not use for |
|---|---|---|---|
| **Noam Pintel, MSc thesis** *The Genomic Basis of Nutritional Traits in Wild Barley* (Tel Hai; supervisor S. Hübner) | `reference_sources/Pintel_MSc_thesis/`: `.docx`, `.pdf`. **Text:** `Pintel_MSc_thesis.md`. **Figures/tables:** `pdf_pages/thesis_page-NN.png`, pp. 18–28 (the docx images are EMF and unreadable). Page index: p18 Table 1 NIR calibration correlations; p20 Fig 1 season-to-season correlations; p22 summary statistics + H² table; p23 Fig 2 trait correlations; p24 Fig 3 environment effects; p26 Fig 4 fitness correlations; p27 Fig 5 Miami plots; p28 Table 3 MorexV2 candidate genes | **the experiment and phenotyping:** net-house setup, plant material, phenotypic measurements, **NIR grain phenotyping and its calibration** (Megazyme β-glucan and starch kits, Kjeldahl protein; **fiber was calibrated with the cultivated-barley package only; no wild-barley lab measurements** (Table 1D has no kit values)); Morex and Clipper as cultivated checks (p22). **Central ideas and ecological explanations**; background literature | numbers: the thesis covers **2 seasons** (2019–20, 2020–21), and a third was added since; GWAS methods (GEMMA, 4 PCs, SimpleM); MorexV2 coordinates and candidate genes. **Season labels are shifted by a year:** the thesis's "Season 2020" = this project's 2019, "Season 2021" = this project's 2020 (its flowering-time GDD ranges match `Z11_GDD` exactly). **Not citable** (unpublished) |
| **Potapenko et al. 2026 *Mol Ecol*** "Haplotype blocks are associated with rapid local adaptation to environmental shifts in wild barley", https://doi.org/10.1111/mec.70540 | `reference_sources/Potapenko_2026_MolEcol/`: `.pdf`; `.md` converted from the PMC full-text XML | **bioinformatics reference and data source.** Evgeny published the 300-accession WGS data here (ENA PRJEB79623). It is also **the same common-garden experiment as this study** (confirmed by the user 2026-09-15: seasons 2019–2021, Hula Valley net-house, randomized blocks). **The 2021 season had 5 blocks**, not the paper's 4 replicates; state this in Methods. **Citable to shorten Methods, because these steps do not depend on the reference genome:** plant material and sampling design; DNA extraction, library preparation and sequencing (the raw reads are the same, ENA PRJEB79623); the common-garden experiment. **Morphological traits: same experiment and same measurement seasons, but not the same set.** This study uses 6 of their 8 (flowering time, plant height, spike length, tillers, grain weight, grain number) and drops flag-leaf length and above-ground biomass. One definition still to verify before citing: their "number of grains per spike" vs this study's `number_of_grains`. **Flowering time is settled**: same measured trait, days from sowing to heading, no GDD (`section_materials/02_materials_methods/NOTE_flowering_time.md`). Data availability wording. **Decision (2026-09-16): describe the traits ourselves and cite this paper once for the common garden** (same trial, location, seasons, block design; say that 2021 had 5 blocks). Then list this study's traits with their definitions and the seasons each was recorded in, and describe the NIR grain traits in full, since they are new here. Read QC and trimming (FastQC, fastp) are probably the same trimmed reads too, but verify in the Methods session. The grain nutrient traits (NIR) are new in this study | **everything from the alignment onward: it is MorexV2 and a different run.** Their alignment (BWA-MEM2 v2.2.1 → MorexV2), duplicate marking, variant calling and filters, and their 33 M SNP set. This study's MorexV3 re-call must be described in full (different reference, different tool versions, different filters, 7,110,996 SNPs). Also not for: haplotype-block results. The research questions differ, so expect little text to reuse |
| **Potapenko et al. 2026 *Mol Biol Evol*** 43:msag051, https://doi.org/10.1093/molbev/msag051 (published version of the bioRxiv genome-size preprint) | `reference_sources/Potapenko_2026_MBE/`: `.pdf`; `.md` converted from the PMC full-text XML | collection sampling design (30 sites × 10 plants, ≥ 10 m apart, two rounds of single-seed descent); environmental data sources (Israeli Meteorological Service, ISRIC soil); net-house location | MorexV2 SNP sets and filters; genome-size and *H. bulbosum* results |

Rules for these sources:
- **Cite published versions only.** Never cite the preprint, the thesis or other unpublished work; TAG's reference list takes published or accepted works only.
- **Paraphrase; do not copy wording.** The journal screens for plagiarism and requires transparency on text recycling.
- **Numbers, settings, versions and counts** taken from them must match the current project files, or be clearly marked as belonging to the earlier dataset.

**Precedence when sources disagree:** project files and logs > **old mini paper** (for topics the analyses did not
change) > section materials (after checking) > published papers > thesis.

## Hard rules

1. **Every number, gene, figure and method detail comes from the NEW project files.** The old
   mini paper is the model for structure, argument flow, wording and style. Even where its numbers
   look unchanged, quote the project files, not the paper.
2. **Module READMEs outrank this file.** This file is a map written 2026-09-15; the
   analyses may be re-run. Check each README's "last run" date before quoting it.
3. **Join on `lead_SNP` or `gene_id`, never on `locus_id`**: locus IDs were renumbered
   between runs (e.g. `fiber_L04` → `fiber_L07`).
4. Next to each number, add a source comment in the md: `<!-- src: 04_.../results_chapter_numbers.txt -->`.
   Pandoc drops HTML comments on export, so they never reach Word.
5. **Missing details: mark them, never guess.** Use the TODO convention below.

## Writing style

- **Results: results only.** Report what was found, at the same level as the mini paper's Results. Broad
  biological interpretation belongs in the Discussion, not here.
- **Discussion: the interpretation.** This is where results are expanded, explained and set against other work.
- **Do not invent a new style.** Keep the mini paper's structure and depth. The change is that this is a full
  paper built on the new analyses, not a condensed one: expand and add where the condensed version had to cut,
  without drifting from its voice or its level of detail.
- **Write for the reader, never against the previous analysis.** (Rule added 2026-09-23 after a Ch. 3 draft.)
  The reader has never seen the earlier runs, the retired parameters or the v1 pipeline, so a sentence that
  exists only to answer a criticism of them is meaningless to them and must not be written. The test: strike
  the clause and ask whether the reader loses anything. Clauses like *"so that neither the grouping nor the
  test was tuned per gene"*, *"one fixed set of parameters for every gene"*, *"used exactly as defined"*,
  *"the correction is applied only once"* are defences of a method choice — they belong in **Materials and
  methods**, or nowhere. Exception: a guarantee the reader needs in order to trust the result (for example
  that the elite cultivars took no part in defining the haplotype groups, which rules out circularity) is a
  fact about the analysis, not a defence, and stays.
- **Report; do not declare what you report.** (Added 2026-09-23.) Sentences of the form *"loci are reported
  alongside genes throughout"*, *"effect sizes are given for every gene"*, *"this is stated for each locus"*
  describe the write-up rather than the study. Drop them and let the text do the thing: if loci belong
  beside genes, put them there. The reader never needs to be told what they are about to be told.
- **Method detail belongs in Materials and methods.** Results names the test only as far as its statistics
  need (e.g. "Kruskal-Wallis, Benjamini-Hochberg within trait, q <= 0.05"); windows, parameter values and
  their justification go to M&M.
- **Materials and methods: E. Potapenko's voice** (user decision 2026-09-30, M&M voice review MM1–MM11; log in
  `new_publishing_paper/build/REVIEW_2026-09-27_Results_Discussion.md`). Model: the M&M of Potapenko et al. 2026
  *Mol Ecol* and *Mol Biol Evol* (`reference_sources/`). Applies to any later change to the M&M:
  - **Aim first, then the tool:** "To define loci, we grouped…", "To test whether…, we…". Use **"we"** at paragraph
    openings (user-approved); passive is fine elsewhere.
  - **Cite, do not repeat, what is already published:** "as described in detail by Potapenko et al. (2026a).
    Briefly, …". **Every detail kept next to such a citation must be verified in the cited paper** (user rule).
  - **Name standard procedures, do not explain them** (tool, version, key setting; no formulas for λGC, BLUP
    arithmetic, imputation internals).
  - **No figure-construction details in the M&M:** display rules go to the figure caption (e.g. the four
    highlighted reaction norms are named in the Fig. 1b caption).
  - **No results in the M&M** (values go to the Results or the Online Resource caption), and a method choice is
    justified in one clause at most.
  - **Plain sentences:** one action per sentence; no lists set off by commas mid-sentence (use a colon or
    parentheses); define a term before using it; the trait set is "nutritional traits" ("grain composition" only
    for the concept).
- **Results, Discussion and Conclusions: S. Hübner's voice first** (user decision 2026-10-01, voice review R1–R5,
  F1–F2, D1–D3, C1; same log). Keep his mini-paper wording, framing and closing sentences wherever the facts hold
  (the user kept, e.g., "indicating the respective genetic control…", his correlation sentences, "Together, these
  associations…", "QTLs" in the Ch. 2 heading, the causal wording of Discussion §2 ¶1). Use E. Potapenko's voice
  (aim first, plain sentences) only for text the mini paper lacks (allelic direction, co-localization, elite lines,
  7H, carrier origin, canonical β-glucan genes). Fix only clarity, logic and facts; "we" sparingly. The cultivars
  were **not assigned haplotypes**: write "genotypes" for them (text and captions). Captions are checked against the
  figure image.

## Marking missing information (becomes a Word comment)

When a detail is missing and writing should go on, insert:

```markdown
[TODO: what is missing and where to find it]{.comment-start id="201" author="TODO"}[[short placeholder]]{.mark}[]{.comment-end id="201"}
```

- In the docx this becomes a **Word comment** holding the TODO note, anchored to the placeholder text **highlighted in yellow** (tested with pandoc 3.11).
- **IDs must be numeric** (Word requires it) **and unique across the whole manuscript**: section number × 100 + n. Title/abstract 1–99; Introduction 101, 102…; Methods 201…; Results 301…; Discussion 401…; Declarations 601….
- List open items: `grep -n 'TODO:' new_publishing_paper/Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` (the working file, renamed 2026-10-01; see Manuscript layout). The Introduction uses 101–103 (item AB1). The M&M uses 201–238 (238 added 2026-09-30, structure item S1).
- **The user's own Word comments** (converted from the Word file on 2026-09-28) keep their author and follow the same
  numbering: 311–319 in the Results, 414–416 in the Discussion (412 and 413 removed 2026-09-30, review item 11). The
  structure items of 2026-09-30 (S. Hübner's figure/table comments) added TODOs 238, 320–325, 421–422. The next free IDs
  are **239** (M&M), **326** (Results) and **423** (Discussion). Retired, never reuse: 403, 404, 405, 417, 418 (418 drafted, never applied). The
  dir-11 carrier-origin Online Resource is cited under TODOs 419 and 420. Pandoc cannot write Word reply threads, so a reply carries `parent="<id>"` in the md and
  appears in Word as a separate comment on the same text.
- `<!-- src: … -->` comments are invisible in Word; TODO comments are visible. Use each for its own purpose.

---

## Module map: old → new

The folder numbers are **not** a one-to-one map. Old `06` = figures; new `06` = gene
selection; new `07` has no old counterpart; new `08` holds the old figures.

| topic | old project | new project | status |
|---|---|---|---|
| Phenotypes, variance partition, H², G×E, ecology/environment correlations (old Fig 1, Fig 2, Table 1) | `00_THIN_...` | `00_THIN_...` | **UNCHANGED**: scripts and all 86 output files byte-identical |
| LD decay | `02_USED_LD_decay_V2_wholegenome` | same | **UNCHANGED outputs, CHANGED role**: see below |
| GWAS + locus definition | `01_USED_GWAS_V2_pipeline` (flat `scripts/`, `results/publication_BonfOnly_BLUP_3PC`) | `01_...` (staged `scripts/`, `results/00_FINAL_BLUP_3PC`) | **CHANGED method + results** |
| Loci → candidate genes | `03_USED_candidate_genes_around_leading_snps` | same | **CHANGED** |
| Haplotype analysis (crosshap) | `04_.../04_runs/candidate_genes_1000bp_V1` | `04_.../04_runs/loci_LDspan_eps06_V4` | **CHANGED method + results** |
| Gene annotation (Swiss-Prot → InterPro → nr) | `05_...` (45 genes) | `05_...` (23 genes, no confidence tier) | **CHANGED** |
| Final gene selection | inside old `04` (`06_publication_figures`) | `06_USED_genes_selected_to_present` + `05_.../08_USED_annotation_master/TRAIT_CANDIDACY.md` | **CHANGED logic** |
| Fiber/starch 7H trade-off direction | `04_.../07_fiber_starch_tradeoff_direction` (on V1) | none | **OLD-ONLY**, not re-run |
| Elite lines vs wild haplotypes | none | `07_USED_elite_lines_compariosn_to_wild_lines` | **NEW** |
| v2 ↔ v3 pruning comparison | none | `01_.../scripts/05_comparison_v2_v3` | **NEW** (methods justification) |
| Assembled figures | `06_USED_figures` | `08_USED_creating_figures` (renamed from `08_USED_figures**NOT_UPDATED**` 2026-09-22) | all rebuilt to TAG spec. **Renumbered 2026-09-30** (S. Hübner's comments): Fig. 1 = old Fig. 1 + old 2b/2c (one heatmap) + old 2d (`make_figure_1.R`); Figs. 2–4 = old Figs. 3–5 (`make_figure_2.R`–`make_figure_4.R`, pixel-identical); old 2a → Online Resource (`make_figure_ESM_site_BLUPs.R`, PDF). See `08_USED_creating_figures/README.md` |

## What changed: key differences to carry into Methods, Results and Discussion

| | old mini paper | new manuscript |
|---|---|---|
| LD pruning for PCA/threshold | 50-SNP window, 5-SNP step → 590,462 SNPs | 1 Mb physical window (`--indep-pairwise 1000kb 1 0.2`) → **111,017** SNPs |
| Bonferroni threshold (α = 0.10) | −log10p 6.7712 | **6.0454** |
| Significance classes | 15 significant + 9 "marginal" SNPs | **52 significant**; no marginal/sub-threshold peaks (dropped 2026-09-09) |
| Locus definition | single-linkage at 200 kb → 20 loci | iterative LD clumping (±2 Mb, 50 kb gap rule, max span 4 Mb) → **36 loci** |
| Gene-search window | lead ±200 kb (LD-decay distance) → 108 genes | LD locus span, no flank → 7.9 Mb, **55 genes**; 21/36 loci contain no gene |
| crosshap parameters | 7 ε × MGmin {2,3}, best result per gene | **fixed ε = 0.6, MGmin = 2** (chosen on genotype-only criteria) |
| Haplotype statistics | within-gene Holm → BH + Bonferroni across genes | **one Kruskal–Wallis per gene, BH on raw p within trait**, q ≤ 0.05 |
| Haplotype result | 76 tested / 48 significant | **30 tested / 23 significant** (12 loci) |
| Protein | 1 locus, 1 gene | 2 loci, 6 genes, **0 testable → no protein gene** |
| Final genes | Pho `3HG0301750`, PHT `3HG0301710` (starch); AP2/ERF (β-glucan); BAHD (fiber) | **GPAT6 `3HG0301300`, GH17 `5HG0487060` (fiber); PHT4;3 `3HG0301710` (starch)**. Only PHT carries over, and β-glucan has no presented gene |
| Gene selection logic | biological plausibility after the stats | **still two-stage**: significance is the gate (only the 23 haplotype-significant genes are considered), then protein function + literature choose among them; rank by q / effect size is **not** used to choose among the 23 (`TRAIT_CANDIDACY.md`) |
| 7H shared fiber/starch locus | 3 dual-trait genes (`7HG0729020/090/100`), inverse haplotype direction; central to the Discussion "carbon-allocation" argument | both loci still exist (`fiber_L17` 7H:573606306; `starch_L05` 7H:573606460, now significant) **but those genes are not in the new step-03 candidates or V4 significant genes**. Do not reuse that evidence; the trade-off argument must stand on the 00_THIN phenotype/ecology results or be re-analysed |
| LD decay | set the ±200 kb gene window | descriptive context only; gene window no longer depends on it |
| Elite lines | not in the paper | new Results subsection (step 07) |
| Figures and tables | Figs. 1–4 + Table 1 (H²) | **4 figures, no tables** (S. Hübner, 2026-09-30): Fig. 1 variance, reaction norms, trait correlations (one heatmap), environment; Fig. 2 Manhattan + QQ; Fig. 3 haplotypes + elite cultivars; Fig. 4 the 7H GDSL gene. Table 1 → variance-components Online Resource; Table 2 → all-candidate-genes Online Resource; old Fig. 2a → site-figure Online Resource |

## Where the numbers come from (read these, in this order, per section)

| section content | source |
|---|---|
| Plant material, field experiment, sequencing, genotype calling (M&M, Data availability) | Field experiment: the same 3-season common garden as Potapenko et al. 2026 *Mol Ecol*; phenotype table `00_THIN_Generate_Plots_For_Publication/all years barley.csv`. Raw reads: **ENA PRJEB79623**, published in Potapenko et al. 2026 *Mol Ecol* (https://doi.org/10.1111/mec.70540; WGS, originally called against MorexV2). Our **MorexV3 re-call**: the GATK/bcftools commands in the header of `/mnt/data/shahar/gwas_barley/data/inputs/morexV3_with_ids.vcf.gz` (300 samples; 10 site-04 samples removed via `samples_to_remove_V3.txt` → 290). Sample IDs `HSsspp` = ENA alias `ss_pp` |
| Phenotypes, H², G×E, ecology | `00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/{tables,png,C2_*}` |
| LD decay | `02_.../results/ld_decay_summary_for_results_chapter.txt` |
| GWAS, PCs, loci | `01_.../README.md` → `results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/` (`Table_loci_master.tsv`, `Table_loci_per_trait.tsv`, `Table_analysis_parameters.tsv`), `03_pc_selection/Table_S2_pc_variance.tsv` |
| Candidate genes | `03_.../results/tables/results_chapter_numbers.txt` |
| Haplotypes | `04_.../04_runs/loci_LDspan_eps06_V4/Stats/results_chapter_numbers.txt`, `Significant_genes/significant_genes.tsv` |
| Annotation | `05_.../08_USED_annotation_master/results/tables/Table_significant_genes_paper.tsv` |
| Final 3 genes | `06_USED_genes_selected_to_present/README.md`, `TRAIT_CANDIDACY.md` |
| **Supplementary tables (Online Resources) built 2026-09-30** (they replace Tables 1 and 2; no tables in the main text) | `new_publishing_paper/supplementary/tables/ESM_variance_components.xlsx` (8 traits: V~G~, V~E~, V~G×E~, V~R~ with %, total, H²; from 00_THIN `Variance_components_GxE.csv`, `H2_GxE.csv`) and `ESM_candidate_genes.xlsx` (all 55 genes, one row each, testability, haplotype test, annotation, carried forward), built by `new_publishing_paper/supplementary/scripts/` (README there). The former **Table 2** source, `06_USED_genes_selected_to_present/results/tables/Table_2_genes_carried_forward.tsv` (`06_.../scripts/01_make_table2_paper.R`), is kept: the gene ESM reads it for "carried forward" and checks its three rows cell by cell. **The relevance-screen categories (strong / plausible / unlikely / no annotation) must not appear in any supplementary file** (user, 2026-09-30) |
| Elite lines | `07_.../results/tables/results_chapter_numbers.txt` + README |

## Figure sources (new)

Renumbered 2026-09-30 (S. Hübner's comments; hand-over `new_publishing_paper/build/STRUCTURE_CHANGES_2026-10.md`):
old Fig. 2 was dissolved and old Figs. 3–5 are now Figs. 2–4. Hidden notes written before that date use the old numbers.

| content | files |
|---|---|
| **Fig. 1** variance partitioning, reaction norms, trait-correlation heatmap (lower triangle + morphological traits), site means × environment; **rebuilt 2026-09-30**, layout and caption approved | `08_USED_creating_figures/Figure_1/Fig1.png` (docx) / `Fig1.tif` (submission), `make_figure_1.R`. Inputs: `00_THIN_.../outputs/subsection_1/tables/` (+ `C2_trait_environment_correlations/tables/`) |
| **Online Resource: site BLUPs by region** (old Fig. 2a) | `08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf`, `make_figure_ESM_site_BLUPs.R` |
| Phenotype / ecology, step-00 working versions | `00_THIN_.../outputs/subsection_1/png/`, `C2_*/figures/` |
| **Fig. 2** Manhattan (loci painted) + QQ, **approved 2026-09-22** (as Fig. 3) | `08_USED_creating_figures/Figure_2/Fig2.png` (docx) / `Fig2.tif` (submission), `make_figure_2.R`. Step-01 working versions: `01_.../02_loci_FINAL/figures/manhattan_*__loci_r05.png`, `04_diagnostics/qq/` |
| LD decay | `02_.../results/figure_ld_decay_{genomewide,perchromosome}.png` |
| **Fig. 3** Haplotype violins + elite barcodes, **approved 2026-09-23** (as Fig. 4) | `08_USED_creating_figures/Figure_3/Fig3.png` (docx) / `Fig3.tif` (submission), `make_figure_3.R`. Working versions: `07_.../results/figures/shared_sites/`, `04_.../04_runs/loci_LDspan_eps06_V4/Significant_genes/{by_trait,CombinedPDF,Heatmaps}`, `06_.../{fiber,starch}/*.pdf` |
| **Fig. 4** GDSL esterase/lipase at the shared 7H signal, **approved 2026-09-24** (as Fig. 5) | `08_USED_creating_figures/Figure_4/Fig4.png` (docx) / `Fig4.tif` (submission), `make_figure_4.R`; analysis `03_01_7H_branch_Starch_Fiber_shared_signal_explore/` |
| Elite vs wild | `07_.../results/figures/{shared_sites,filled_marked,filled_silent}/` |

## Statements the READMEs require in the manuscript

- ~~**ε is not an r² cutoff** — state in the methods~~ **Withdrawn 2026-09-23 (user decision):** not written in the manuscript. Methods give ε = 0.6 and MGmin = 2 as fixed, with the genotype-only grounds for the choice (`04_.../README.md` § The two fixed parameters), and no parameter-space explanation. The point stays recorded in the step-04 README for future runs.
- BH remains valid under LD (positive regression dependence), and genes that could not be tested are excluded from the BH denominator but reported.
- **Elite lines:** Morex (the reference genome) is itself an elite cultivar. For GH17 and PHT4;3 the elite "match" is reference identity only; **only GPAT6 supports a shared-haplotype claim.** The figure shows which haplotypes elite lines carry, not that breeding selected them.

---

## Manuscript layout (`new_publishing_paper/`)

```
Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md               <- THE working file (title … Conclusions)
Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley_reading_copy.md  <- same text without comments / hidden notes (generated)
05_references.md  06_statements_declarations.md   (not written yet; title page/declarations may go into the working file too)
figures/  tables/  supplementary/  build/
```
**The whole manuscript text is one file (user decisions 2026-09-30 and 2026-10-01): `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` is the working file. Edit only that file.** Since 2026-10-01 (item AB1) it opens with the title, Key message, Abstract, Keywords and Introduction, and it was **renamed after the title** (it was `Methods_Results_Discussion_Conclusions.md`; older notes and logs use that name). The M&M + R + D + C part It was built by concatenating the
approved M&M (`02_materials_methods.md`, written 2026-09-29/30 from `build/MM_BLUEPRINT.md`) and the R+D+C working file
`build/UPDATED_Results_Discussion_Conclusions.md` (itself converted 2026-09-28 from the user's Word-edited
`build/Results_Discussion_Conclusions.docx`, all Word comments kept, hidden notes restored). Bodies unchanged; only the
figure paths were re-rooted from `build/` to `new_publishing_paper/` (verified). Its first-line comment records this.
- **Reading copy:** `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley_reading_copy.md` has every Word comment and hidden `<!-- -->`
  note removed, for reading. **Never edit it; regenerate it after every change** to the working file:
  `python3 build/make_reading_copy.py Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley_reading_copy.md`
  (run from `new_publishing_paper/`).
- **Removed 2026-09-30 (user):** `02_materials_methods.md` (now inside the working file), the frozen copies
  `03_results.md`, `04_discussion.md`, `Results_Discussion.md`, `Results_Discussion_Conclusions.md`, and the M&M-only
  reading copy `build/MM_reading_copy.md`. Also removed the same day (user): the frozen pre-review copies
  `Methods_Results_Discussion_Conclusions_BEFORE_REVIEWING_AND_SHORTENING.md` (+ its reading copy), after the M&M voice
  review was approved.
- **M&M voice review, 2026-09-30 (done):** items MM1–MM11 (`build/proposed/MM*_edits.json`, previews `MM*_PROPOSED.md`)
  rewrote the M&M in E. Potapenko's voice (rules under "Writing style"), about 3,530 → 3,190 visible words. Structure
  changes: the 7H subsection is now the last paragraph of "Functional annotation and choice of candidate genes"; the
  M&M ends with two new subsections, **Data and material availability** (ENA accession, material-sharing placeholder,
  TODOs 234–235) and **Use of generative artificial intelligence** (TODO 236). All M&M TODOs 201–237 are still open.
- **Title, Key message, Abstract, Keywords, Introduction, 2026-10-01 (item AB1, approved):** drafted from the mini paper (S. Hübner's text) and updated to this study; S. Hübner's title kept; TODOs 101–103. Title page (authors, affiliations, ORCID, Acknowledgments) and author contributions: Declarations session.
- **Results + Discussion + Conclusions voice review, 2026-10-01 (done):** items R1–R4 (R5 discarded), F1–F2 (Fig. 1–3
  captions; Fig. 4 caption in R4), D1–D3 (+ D1b), C1 (`build/proposed/`; log at the end of the review file). Rules
  under "Writing style". Results about 2,690 and Discussion + Conclusions about 2,320 visible words (with captions).
- **Structure changes, 2026-09-30 (done; S. Hübner's comments on the figures and tables):** items S1–S5 (`build/proposed/S*_edits.json`;
  questions and answers `build/proposed/S0_QUESTIONS.md`; hand-over **`build/STRUCTURE_CHANGES_2026-10.md`**). Tables 1 and 2 left
  the main text (→ Online Resources); new Fig. 1 (old 1 + old 2b/2c/2d); old Fig. 2a → Online Resource; Figs. 3–5 → 2–4. No reply to
  Sariel's Online Resource comment (user). Online Resource numbering is still open (placeholders "Online Resource N", each with a TODO).
- **`build/DECLARATIONS_HANDOVER.md`:** text moved out of the M&M for the Declarations session (the Zenodo VCF /
  BLUP Online Resource / GitHub code sentence of the Data availability statement), word for word.
- **Superseded, kept for the record:** `build/UPDATED_Results_Discussion_Conclusions.md` (top note says so; never edit
  it). Historical notes in `build/` (`BLUEPRINT.md`, `MM_BLUEPRINT.md`, `DISCUSSION_CANDIDATES.md`, `CHANGED.md`,
  `MM_HANDOVER_from_review.md`, the review log) may name the old file names inside their records; their headers point
  to the working file.

**Changing the working file: propose first, apply after approval (since 2026-09-27).** `build/proposed/review_tool.py`
targets `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` since 2026-10-01 (`WORK` in the script; before, `Methods_Results_Discussion_Conclusions.md`) and regenerates the reading
copy after every `apply`. It matches the exact old text, not line numbers, so the merge does not affect it.
- Each change is an item with an edit file, `build/proposed/<ID>_edits.json`, holding exact old text and new
  text marked `{-deleted-}{+inserted+}`.
- `python3 build/proposed/review_tool.py preview <ID>` writes `<ID>_PROPOSED.md`. The user reads it in the VS
  Code markdown preview: red strikethrough for deletions, green for insertions, each paragraph in full.
- `review_tool.py apply <ID>` writes the change only after approval, and refuses if anything else in the file
  would change.
- **The user reads only the proposal file, not the chat.** Put every question in the file's "Questions for
  you" block (`{"questions": [...], "edits": [...]}`).
- Decisions are logged at the end of `build/REVIEW_2026-09-27_Results_Discussion.md`.

File order follows TAG (declarations after references); formatting rules are in `TAG_requirements.md`.
- **Figures:** in the md, link the source output directly (relative path). TAG wants figures embedded in the text, so the docx build embeds them. For submission, export `figures/Fig1.eps|tif` (Arial 8–12 pt, 84 or 174 mm wide, panels a/b/c).
- **Tables:** **none in the main text since 2026-09-30** (S. Hübner; user decision). Tables go to the Online Resources as `.xlsx`, built by `supplementary/scripts/` into `supplementary/tables/` (named `ESM_N.xlsx` when numbered), cited as "Online Resource N". In TAG an Online Resource is a supplementary file (ESM); nothing is uploaded elsewhere.
- **References:** TAG author-year style, `(Godfray et al. 2010)` with no comma. The old paper's citations and reference list do **not** follow TAG (commas, `&`, full journal names, no DOIs). Reuse its content, not its format.

## Tools

Pandoc 3.11 (not on PATH): `/mnt/data/shahar/gwas_barley/tools/pandoc-3.11/bin/pandoc`
```bash
P=/mnt/data/shahar/gwas_barley/tools/pandoc-3.11/bin/pandoc
# md → docx (all sections, in order; styles from a Word template if present). Run from new_publishing_paper/,
# because the figure paths in the md are relative to it. The working file holds title … Conclusions (renamed 2026-10-01);
# add 05_references.md / 06_statements_declarations.md when they exist.
cd new_publishing_paper
W=Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md
$P "$W" [05_*.md 06_*.md] -o build/manuscript.docx --resource-path=. [--reference-doc=template.docx]
# reading copy (after every change to the working file)
python3 build/make_reading_copy.py "$W" "${W%.md}_reading_copy.md"
# docx → md (e.g. supervisor's edited version); --track-changes=all keeps Word comments and tracked changes
$P in.docx -t markdown --wrap=none --track-changes=all --extract-media=media_in -o in.md
# open TODO items (the working file only; the reading copy has none)
grep -n 'TODO:' new_publishing_paper/Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md
```
Word conversions leave `[..]{dir="rtl"}` spans (Hebrew keyboard artefacts). Strip them.

## TAG journal requirements

Kept separately in [`TAG_requirements.md`](TAG_requirements.md). Read it before drafting any section,
and before building the docx.
