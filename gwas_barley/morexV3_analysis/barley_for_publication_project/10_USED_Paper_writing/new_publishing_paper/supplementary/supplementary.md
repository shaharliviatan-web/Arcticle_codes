# Supplementary material

**Where everything is: `README.md` in this folder** (folder map, upload folder, rebuild command). This file is the list.

## Current numbering (2026-10-17) — BUILT

S. Hübner's comment on the review copy: TAG papers call the supplement "Supplementary" and number figures and tables in two running series (Fig. S1…, Table S1…), cited that way in the text. User decision 2026-10-17: two files, **`ESM_1.xlsx` = Supplementary Tables S1–S22** (one sheet per table) and **`ESM_2.pdf` = Supplementary Figures S1–S5** (each with its legend); one number per table or figure, numbered by first citation (checked by `scripts/build_ESM_tables.py` on every build). The content is unchanged (the 22 tables checked cell by cell against the former Online Resource sheets); only the numbering and the assembly changed. Manuscript item SUPP2.

| new | from | caption (start) |
|---|---|---|
| Table S1 | part 1 (former Online Resource 1), sheet “Sites” | Sampling sites: the 30 sites with their ecological region, coordinates and ten environmental variables |
| Table S2 | part 2 (former Online Resource 2), sheet “NIR calibration” | Near-infrared (NIR) calibration: for the 30 calibration accessions (one per sampling site), the laboratory values, the NIR predictions of the same gra… |
| Table S3 | part 2 (former Online Resource 2), sheet “Calibration correlations” | Validation of the NIR calibration: Pearson correlations among the laboratory values, the NIR predictions and the 2019–20 NIR values of the 30 calibrat… |
| Table S4 | part 2 (former Online Resource 2), sheet “Blocks 2019–20” | Pearson correlations between the blocks of the 2019–20 season for the four nutritional traits, from the raw NIR values of the accessions grown in both… |
| Table S5 | part 1 (former Online Resource 1), sheet “Accessions” | The 290 analyzed accessions: identifiers in the European Nucleotide Archive (project PRJEB79623), site, ecological region and the genotype BLUPs of th… |
| Table S6 | part 3 (former Online Resource 3), sheet “Variance components” | Variance components and broad-sense heritability (H²) of the four nutritional and four morphological and phenological traits, from a linear mixed mode… |
| Table S7 | part 1 (former Online Resource 1), sheet “Environment correlations” | Pearson correlations among the ten environmental variables across the 29 sites of the analyzed accessions |
| Table S8 | part 5 (former Online Resource 5), sheet “Loci” | The 36 loci of the genome-wide association analysis of the four nutritional traits (EMMAX; 290 accessions; MorexV3 coordinates): position, span, lead … |
| Table S9 | part 5 (former Online Resource 5), sheet “Significant SNPs” | The 52 SNPs that exceeded the genome-wide threshold (Bonferroni, α = 0.10; −log10(P) = 6.0454), with their locus, allele frequencies and effects |
| Table S10 | part 6 (former Online Resource 6), sheet “Candidate genes” | The 55 candidate genes of the 36 loci, one row per gene: locus and position, whether the gene's haplotypes could be tested and why not, the Kruskal–Wa… |
| Table S11 | part 7 (former Online Resource 7), sheet “Summary” | Choice of the crosshap clustering radius (ε), summary: for five values of ε (MGmin = 2), the number of the 55 candidate genes whose haplotypes could b… |
| Table S12 | part 7 (former Online Resource 7), sheet “Per gene” | Choice of the crosshap clustering radius (ε), per gene: the same quantities for each of the 55 candidate genes at each of the five values of ε (MGmin … |
| Table S13 | part 8 (former Online Resource 8), sheet “Canonical genes” | The ten canonical (1,3;1,4)-β-glucan synthase and endohydrolase genes located on the MorexV3 genome from their reference protein sequences, and the di… |
| Table S14 | part 10 (former Online Resource 10), sheet “Site accounting” | Comparison with elite cultivars: per gene, the SNP sites of the wild and elite MorexV3 call sets in the gene window (gene ± 1 kb), the sites compared … |
| Table S15 | part 10 (former Online Resource 10), sheet “Excluded sites” | Comparison with elite cultivars: every SNP site excluded from the comparison, with the reason |
| Table S16 | part 10 (former Online Resource 10), sheet “Screen” | Call-rate screen of the 136 spring elite lines of the barley pangenome panel, by pooled call rate over the compared sites of GPAT6, GH17 and PHT4;3… |
| Table S17 | part 10 (former Online Resource 10), sheet “Cultivars” | The five spring malting cultivars compared with the wild haplotype groups, and their call rates at the four presented genes |
| Table S18 | part 11 (former Online Resource 11), sheet “Lead SNPs – tests” | Geographic origin of the carriers of the lead-SNP alleles of the 36 loci: carriers per ecological region, and the site-permutation tests of desert enr… |
| Table S19 | part 11 (former Online Resource 11), sheet “Lead SNPs – by region” | Minor- and major-allele carriers and missing calls of each lead SNP, per ecological region |
| Table S20 | part 11 (former Online Resource 11), sheet “Haplotype groups – tests” | Geographic origin of the haplotype groups of the four presented genes: members per ecological region, and the site-permutation tests of desert enrichm… |
| Table S21 | part 11 (former Online Resource 11), sheet “Haplotype groups – by region” | Members of each haplotype group of the four presented genes per ecological region, and the unassigned accessions |
| Table S22 | part 11 (former Online Resource 11), sheet “Carriers vs site-mates” | Comparison of the minor-allele carriers of each lead SNP with the major-allele carriers of their own site |
| Fig. S1 | former Online Resource 4 | Decay of linkage disequilibrium with physical distance across the genome (a) and per chromosome (b) |
| Fig. S2 | former Online Resource 9 | crosshap haplotype plots of GPAT6, GH17, PHT4;3 and the GDSL esterase/lipase gene |
| Fig. S3 | former Online Resource 12 | Geographic origin of the accessions carrying the lead-SNP alleles and the haplotypes of the four candidate genes |
| Fig. S4 | former Online Resource 12 | Desert enrichment of the lead-SNP alleles and the haplotypes of the four candidate genes |
| Fig. S5 | former Online Resource 13 | Genotype BLUPs of the four nutritional traits across the 29 sampling sites, by ecological region |

**Upload:** `for_submission/` (`ESM_1.xlsx`, `ESM_2.pdf`, `ESM_all_supplementary.zip`, `CAPTIONS.md`). **Review copy for S. Hübner:** `for_review_Hubner/` (`Supplementary_Tables.xlsx`, `Supplementary_Figures.docx`, plus the two uncited raw-phenotype files and a zip).

## Former numbering (2026-10-01 to 2026-10-17; superseded by the table above, kept for the record)

Numbered by first citation in the manuscript (checked: 1 → 13 in reading order). Files in this folder, built by
`scripts/build_all.sh` (README there: sources, checks, layout). Questions and answers: `../build/proposed/ESM0_QUESTIONS.md`;
hand-over: `../build/SUPPLEMENTARY_DONE_2026-10.md`. Every citation is resolved; the "supmat" Word comments that tracked them were deleted on 2026-10-01.

| # | file | content | status |
|---|---|---|---|
| 1 | `ESM_1.xlsx` | sampling sites, accessions (ENA identifiers, genotype BLUPs of 10 traits), site environments, environmental correlations | built; no unit for plant height and grain weight (as Mol Ecol) |
| 2 | `ESM_2.xlsx` | validation of the NIR calibration (30 calibration accessions; N. Pintel's files) and inter-block correlations 2019–20 | built |
| 3 | `ESM_3.xlsx` | variance components and H² | built |
| 4 | `ESM_4.pdf` | LD decay, a genome-wide, b per chromosome | built |
| 5 | `ESM_5.xlsx` | 36 loci; 52 significant SNPs | built |
| 6 | `ESM_6.xlsx` | 55 candidate genes | built |
| 7 | `ESM_7.xlsx` | choice of ε | built |
| 8 | `ESM_8.xlsx` | canonical β-glucan genes on MorexV3 | built |
| 9 | `ESM_9.pdf` | crosshap plots of the four presented genes, GDSL for fiber and starch (replaces the 7H MGmin × ε grid, dropped by the user) | built |
| 10 | `ESM_10.xlsx` | elite cultivars: lines, screen, site accounting, excluded sites | built |
| 11 | `ESM_11.xlsx` | carrier origin, tables | built |
| 12 | `ESM_12.pdf` | carrier origin, Figs. S1–S2 | built |
| 13 | `ESM_13.pdf` | site figure (former Fig. 2a) | built |

**For download / upload:** `for_submission/` — `xlsx/` (9), `pdf/` (4), `ESM_all_Online_Resources.zip` (all 13) and `CAPTIONS.md` (file → caption for the submission system); refreshed by `scripts/build_all.sh`. One file per Online Resource, as TAG asks (no merged PDF).

Title blocks: author names, the corresponding author's affiliation and e-mail are **placeholders** (`scripts/esm_common.py`;
manuscript checklist item 1). Captions in the manuscript: last section, "Supplementary Information" (item CAP1, applied).
**Not cited in the text** (user, 2026-10-01): `unreferenced/` — the raw-phenotype checks; the user decides later whether
they join the numbered set.

## Planning list (2026-09-29 to 2026-09-30; kept for the record, superseded by the table above)

Collected while writing. Started 2026-09-29.

**TAG rules:** cite as "Online Resource N"; files `ESM_N.xlsx` / `ESM_N.pdf`, numbered in citation order; each file carries
the article title, journal, authors and the corresponding author's affiliation and e-mail; a short caption for each in the
manuscript. Tables → .xlsx/.csv; figure or text collections → PDF. A large supplement is welcome.

**2026-09-30, S. Hübner's structural comments (hand-over `build/STRUCTURE_CHANGES_2026-10.md`):** no tables remain in the
main text. Table 1 → the variance-components Online Resource; Table 2 → the all-candidate-genes Online Resource (which
also absorbs the 55-gene table and the 23-gene annotation table, formerly "Online Resource 2"); old Fig. 2a → the
site-figure Online Resource. The three are **built** (files below). In TAG, "Online Resource N" = a supplementary
file (ESM); nothing is uploaded elsewhere. Numbering was still open then (fixed 2026-10-01, above).
Built ESM files: tables in `supplementary/tables/` (scripts in `supplementary/scripts/`, README there); figures in
`08_USED_creating_figures/Figure_ESM_*/`.

**Rule (user, 2026-09-30): the categories of the biological-relevance screen (strong / plausible / unlikely / no
annotation, `TRAIT_CANDIDACY.md`, `trait_candidate_calls.tsv`) must not appear in any supplementary file** — not as a
column, not in a caption, and not in file names (e.g. the `T1__…` prefix of `06_USED_genes_selected_to_present/*/`).
`scripts/make_ESM_candidate_genes.py` stops if they appear in its output.

| # | content | why | cited from | source |
|---|---|---|---|---|
| 1 | 36-locus table (+ member SNP counts, spans, co-localized regions' coordinates) | coordinates removed from the text | Results Ch. 2 (TODO 301, 309, 310; comment 316); M&M §6 ¶1 (TODO 220) | `01_…/02_loci_FINAL/tables/Table_loci_master.tsv` |
| – | 52 significant SNPs: MAF, β, *P* | per-SNP detail | Results Ch. 2 (comments 314, 315) | `01_…/00_FINAL_BLUP_3PC/00_snp_level/significant_snps.tsv` |
| – | **All candidate genes (BUILT 2026-09-30)**: the 55 genes, one row each, nothing filtered — trait, locus, lead SNP, position, distance to lead, MorexV3 description, SNPs in window, testable Y/N + reason, accessions grouped, groups, group sizes, KW *P*, BH *q*, η², SD difference, significant Y/N; for the 23 significant genes the annotation call, its source, the call of each source (Swiss-Prot / InterPro / nr), best Swiss-Prot hit + identity / coverage + accepted, Pfam, best nr hit; carried forward Y/N. **Replaces Table 2, the 23-gene table (formerly "Online Resource 2") and the 55-gene table** (S. Hübner, 2026-09-30) | gene search, testability, haplotype tests, annotation, choice of the 3 genes; the 3 carried-forward rows reproduce the former Table 2 cell by cell (script check) | Results Ch. 3 ¶2 (TODO 303; comments 317, 318), ¶3–¶5 (TODOs 323–325, were "Table 2"); M&M §6 ¶2 (TODO 222), §7 ¶3 (TODO 224), §8 ¶1 (TODO 227) | **`supplementary/tables/ESM_candidate_genes.xlsx`** (+ `.tsv`), built by `supplementary/scripts/make_ESM_candidate_genes.py` from `03_…/candidate_genes.tsv`, `04_…/Stats/gene_results.tsv` + `genes_not_tested.tsv`, `05_…/08_USED_annotation_master/results/tables/Table_significant_genes_{paper,annotated}.tsv`, `06_…/Table_2_genes_carried_forward.tsv`. Relevance-screen category **not** included (user, Q6; the script checks) |
| – | **Variance components and H² (BUILT 2026-09-30)**: 8 traits × variance and % of total for V~G~, V~E~ (season & block), V~G×E~, V~R~, total variance, H². **Replaces Table 1** (S. Hübner: "a full table with all variation components") | the full variance partition behind Fig. 1a and the H² values in the text; % and H² checked against Fig. 1a and the former Table 1 (script check; tillers H² "< 0.001") | M&M §4 ¶3 (TODO 238); Results Ch. 1 ¶1 (TODO 320), Ch. 2 ¶2 (TODO 322) | **`supplementary/tables/ESM_variance_components.xlsx`** (+ `.tsv`), built by `supplementary/scripts/make_ESM_variance_components.py` from `00_THIN_…/outputs/subsection_1/tables/Variance_components_GxE.csv` (+ `_WIDE.csv`), `H2_GxE.csv` |
| – | **Site figure, former Fig. 2a (BUILT 2026-09-30; PDF)**: centered trait BLUPs of the four nutritional traits across the 29 sampling sites, ordered and colored by ecological region. Caption draft (the old Fig. 2a caption): "Centered trait BLUPs of the four nutritional traits across 29 sampling sites, ordered and colored by ecological region" | the regional pattern (north vs desert) described in Results Ch. 1 ¶2 and the Discussion | Results Ch. 1 ¶2 (TODO 321); Discussion §1 (TODO 421), §2 ¶1 (TODO 422) (applied 2026-09-30, items S4–S5) | **`08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf`** (+ PNG preview), `make_figure_ESM_site_BLUPs.R`, from `00_THIN_…/tables/A4_site_boxplot_data.csv`. Could be merged with the sampling-site table below (same sites and regions) |
| 3 | Elite comparison: site accounting and removed sites | only shared sites compared | Results Ch. 3 (TODO 306; comment 319); M&M §10 ¶1 (TODO 232; moved from ¶2 by item MM4, 2026-09-30) | `07_…/results/tables/Table_site_overlap.tsv`, `Table_removed_sites.tsv` |
| – | Elite lines: 5 cultivars, call rates per gene (GPAT6, GH17, PHT4;3, pooled, **GDSL**) + 136-line screen | line choice; the 85% rule is pooled over 3 genes (136 lines scored, 91 pass), and at the GDSL gene only 2 lines pass (this result was removed from the M&M text by item MM4, 2026-09-30, so the table must show the GDSL call rates) | M&M §10 ¶1 (TODO 231) | `07_…/results/tables/Table_elite_lines.tsv`, `Table_elite_line_screen.tsv`; `03_01_…/results/tables/Table_elite_line_screen.tsv` |
| – | 7H GDSL gene: MGmin × ε grid, site overlap, elite genotypes | parameter choice | Results Ch. 4; M&M §9 ¶2 (TODO 229) | `03_01_…/results/tables/param_grid_eps0.05-1.5_MGmin2-3.tsv` (written by `scripts/06_param_grid.R`, 2026-09-30; decide whether to show its kw_p / eta2 / delta columns) + `Table_*` |
| – | Carrier origin of the lead-SNP alleles and the haplotype groups: **T02** (+ T01 as sheet 2), **T04** (+ T03), **Figs S1–S2**; supporting T05, T06, T08 (T07, T09–T11 not used: carrier sharing and kinship are not in the text) | backs the Discussion; analysis finished 2026-09-30; not in the Results (user) | Discussion §3 ¶2 (TODO 419), "Wild alleles" (TODO 420); M&M §11 (TODO 237; methods cover T02, T04 and T08 only) | `11_USED_check_low_maf_geographic_distribution_minor_allele/results/` — captions ready in `results/CAPTIONS.md` |
| – | Canonical β-glucan genes: UniProt accession → MorexV3 gene ID, position, identity / coverage, reciprocal best hit, tBLASTn agreement, distance to the nearest β-glucan lead | backs the Discussion; shows how the genes were located on MorexV3 | Discussion (TODO 410; comment 414); M&M §8 ¶3 (TODO 228) | mapping columns of `05_…/09_…/results/tables/canonical_betaglucan_gene_distances.tsv` (its distances are stale) + `canonical_gene_distances_current_loci.tsv` |
| – | Imputation summary (optional): missing rate before imputation, sNMF cross-entropy of 10 runs at K = 3 | imputation used for LD in the haplotype analysis | M&M §7 ¶1 | `/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_imputation/*.snmfProject` (**check**) |
| – | Variant-calling tools, versions and filter thresholds (optional; Evgeny suggests it) | the calling pipeline in one place | M&M §1 ¶3–4 | VCF header + Evgeny's reply (`build/MM_BLUEPRINT.md`, Appendix A) |
| – | LD decay, genome-wide + per chromosome | descriptive. **Caption must carry** (moved out of the M&M text, item MM4, 2026-09-30): pairs with r² ≥ 0.001 within 2.5 Mb, 1-kb bins, LOESS span 0.10 weighted by pair count; the smoothed r² falls to 0.2 at about 188 kb | M&M §5 ¶3 only (TODO 219, `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md`) | `02_USED_LD_decay_V2_wholegenome/results/figure_ld_decay_*.png`; value from `02_…/results/ld_decay_summary_for_results_chapter.txt` |
| – | crosshap ε comparison on the 55 candidate genes (ε 0.2–1.0: genes testable, accessions assigned, median groups) | why ε = 0.6 (balance of testable genes and assigned accessions) | M&M §7 ¶2 (TODO 223, `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md`) | `04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Diagnostics/epsilon_choice_supplementary.tsv` (summary) + `epsilon_choice_per_gene.tsv` (55 genes × 5 ε), from `01_scripts/diagnostics/parameter_choice_tables.R` (2026-09-30) |
| – | PC variance (scree), 3 PCs | structure correction | **not cited** (user 2026-09-29: the M&M cites Potapenko et al. 2026b (Mol Ecol; 2026a before the REF1 letter swap) for 3 PCs, no justification in the text); keep only if wanted | `01_…/00_FINAL_BLUP_3PC/03_pc_selection/Table_S2_pc_variance.tsv` |
| – | NIR calibration validation (Pearson r) | calibration quality | M&M §3 ¶2 (TODO 214, `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md`) | old M&M 2.5.2026 / thesis Table 1 — **Noam to verify** |
| – | Block 4 of 2019–20: inter-block correlations of the 4 NIR traits | why it was removed (β-glucan) | M&M §4 ¶1 (TODO 215, `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md`) | `12_USED_adittional_checks_for_supplementary/01_block4_2019_interblock_correlations/results/tables/Table_OnlineResource_2019_block_correlations.tsv` (generated 2026-09-29; optional figure `results/figures/Fig_interblock_correlation_heatmaps.png`) |
| – | Genotype BLUPs, 290 accessions × traits | TAG expects phenotype data for QTL studies (as in Potapenko et al.: phenotypes as supplementary, no separate deposit) | M&M closing paragraph (TODO 234, `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md`); Data availability | `00_THIN_…/outputs/subsection_1/tables/BLUP_10_traits.csv` |
| – | Sampling sites and accessions: site number, name, ecotype region or intermediate zone, coordinates; per accession its ENA sample alias (VCF `HS0103` = ENA `01_03`) | the site-to-region division behind the site figure (former Fig. 2a); links accessions to the raw reads | M&M §1 ¶1 (TODO 218; moved there from §4 ¶5 by item MM3, 2026-09-30) and closing paragraph (TODO 234), `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` | `00_THIN_…/01_phenotypic_analysis_no_GxE_v3_streamlined.R` (SITE_NAME_MAP, SITE_REGION_MAP) + longitude/latitude in `00_THIN_…/all years barley.csv`; could be merged with the environmental-variables table below |
| – | Site environmental variables + collinearity matrix | Fig. 1d inputs (former Fig. 2d); the two variables dropped for collinearity | M&M §4 ¶5 (TODO 217, `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md`) | `00_THIN_…/outputs/subsection_1/C2_trait_environment_correlations/tables/C2corr_site_level_dataset.csv` + `C2_environmental_predictor_collinearity/tables/C2col_environment_correlation_matrix.csv` |

**Optional ideas:** per-gene crosshap figures for the 23 genes (PDF; without the screen-tier file names, see the rule above); Spearman results for Fig. 1d (former Fig. 2d); per-season
phenotype summaries.
