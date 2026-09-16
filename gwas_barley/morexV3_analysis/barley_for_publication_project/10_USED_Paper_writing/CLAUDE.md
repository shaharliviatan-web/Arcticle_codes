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

## Marking missing information (becomes a Word comment)

When a detail is missing and writing should go on, insert:

```markdown
[TODO: what is missing and where to find it]{.comment-start id="201" author="TODO"}[[short placeholder]]{.mark}[]{.comment-end id="201"}
```

- In the docx this becomes a **Word comment** holding the TODO note, anchored to the placeholder text **highlighted in yellow** (tested with pandoc 3.11).
- **IDs must be numeric** (Word requires it) **and unique across the whole manuscript**: section number × 100 + n. Title/abstract 1–99; Introduction 101, 102…; Methods 201…; Results 301…; Discussion 401…; Declarations 601….
- List open items: `grep -rn 'TODO:' new_publishing_paper/*.md`.
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
| Assembled figures | `06_USED_figures` | `08_USED_figures**NOT_UPDATED**` | byte-identical to old: **not a source**; old Fig 3–4 content is obsolete |

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
| Gene selection logic | biological plausibility after the stats | protein function + literature only; stats not used for selection (`TRAIT_CANDIDACY.md`) |
| 7H shared fiber/starch locus | 3 dual-trait genes (`7HG0729020/090/100`), inverse haplotype direction; central to the Discussion "carbon-allocation" argument | both loci still exist (`fiber_L17` 7H:573606306; `starch_L05` 7H:573606460, now significant) **but those genes are not in the new step-03 candidates or V4 significant genes**. Do not reuse that evidence; the trade-off argument must stand on the 00_THIN phenotype/ecology results or be re-analysed |
| LD decay | set the ±200 kb gene window | descriptive context only; gene window no longer depends on it |
| Elite lines | not in the paper | new Results subsection (step 07) |

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
| Elite lines | `07_.../results/tables/results_chapter_numbers.txt` + README |

## Figure sources (new)

| content | files |
|---|---|
| Phenotype / ecology | `00_THIN_.../outputs/subsection_1/png/`, `C2_*/figures/` |
| Manhattan (loci painted), QQ | `01_.../results/00_FINAL_BLUP_3PC/02_loci_FINAL/figures/manhattan_*__loci_r05.png`, `04_diagnostics/qq/` |
| LD decay | `02_.../results/figure_ld_decay_{genomewide,perchromosome}.png` |
| Haplotypes | `04_.../04_runs/loci_LDspan_eps06_V4/Significant_genes/{by_trait,CombinedPDF,Heatmaps}`; `06_.../{fiber,starch}/*.pdf` |
| Elite vs wild | `07_.../results/figures/{shared_sites,filled_marked,filled_silent}/` |

## Statements the READMEs require in the manuscript

- **ε is not an r² cutoff.** It is a Euclidean radius in r²-profile space; its stringency depends on the SNP count (04 README: "State this in the methods").
- BH remains valid under LD (positive regression dependence), and genes that could not be tested are excluded from the BH denominator but reported.
- **Elite lines:** Morex (the reference genome) is itself an elite cultivar. For GH17 and PHT4;3 the elite "match" is reference identity only; **only GPAT6 supports a shared-haplotype claim.** The figure shows which haplotypes elite lines carry, not that breeding selected them.

---

## Manuscript layout (`new_publishing_paper/`)

```
00_title_abstract_keymessage.md  01_introduction.md  02_materials_methods.md
03_results.md  04_discussion.md  05_references.md  06_statements_declarations.md
figures/  tables/  supplementary/  build/
```
File order follows TAG (declarations after references); formatting rules are in `TAG_requirements.md`.
- **Figures:** in the md, link the source output directly (relative path). TAG wants figures embedded in the text, so the docx build embeds them. For submission, export `figures/Fig1.eps|tif` (Arial 8–12 pt, 84 or 174 mm wide, panels a/b/c).
- **Tables:** the script TSV/CSV is the source. Small tables go in as pipe tables (they become real Word tables); large ones go to `supplementary/ESM_N.xlsx`, cited as "Online Resource N".
- **References:** TAG author-year style, `(Godfray et al. 2010)` with no comma. The old paper's citations and reference list do **not** follow TAG (commas, `&`, full journal names, no DOIs). Reuse its content, not its format.

## Tools

Pandoc 3.11 (not on PATH): `/mnt/data/shahar/gwas_barley/tools/pandoc-3.11/bin/pandoc`
```bash
P=/mnt/data/shahar/gwas_barley/tools/pandoc-3.11/bin/pandoc
# md → docx (all sections, in order; styles from a Word template if present)
$P new_publishing_paper/0*.md -o new_publishing_paper/build/manuscript.docx --resource-path=.:.. [--reference-doc=template.docx]
# docx → md (e.g. supervisor's edited version); --track-changes=all keeps Word comments and tracked changes
$P in.docx -t markdown --wrap=none --track-changes=all --extract-media=media_in -o in.md
# open TODO items
grep -rn 'TODO:' new_publishing_paper/*.md
```
Word conversions leave `[..]{dir="rtl"}` spans (Hebrew keyboard artefacts). Strip them.

## TAG journal requirements

Kept separately in [`TAG_requirements.md`](TAG_requirements.md). Read it before drafting any section,
and before building the docx.
