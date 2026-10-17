# ESM0 — Online Resources: proposed final list and open decisions

2026-10-01. Task: arrange the supplementary material (checklist item 5).

> **Status (latest first)**
> - **All questions answered (Q1–Q17) and applied. All 37 supmat comments ✅.**
> - **Open, not questions:** author names, affiliation and e-mail for the title blocks (checklist item 1); TODO 213 (S. Hübner: one NIR calibration for all seasons?).
> - **Also open (not questions):** author names, affiliation and e-mail for the title blocks (checklist item 1).
>
> Q1–Q16 were answered by you in chat (recorded below in *italics*). Write under each open question, or "OK" to take the
> recommendation.

**Your broad comments (2026-10-01), and what I do with them:**
1. *"supmat" Word comment on every supplementary reference* → item **SUPMAT1** adds one to each citation (and to the
   two places that still need one); format and list in **Q13**. Updated to ✅ as each one is resolved.
2. *Captions — required, and where?* → **TAG requires a concise caption for each ESM file in the manuscript** ("The
   manuscript contains a descriptive caption for each supplementary material"; guidelines p. 39). TAG also requires
   the title block *inside* each file. It does not require explanations inside the file, but a caption plus column
   definitions inside the file is good practice (files are published as received, so the reader may open one without
   the paper). So: caption in the manuscript (**Q9**: where), and the same caption plus column definitions on the
   first sheet / first page of each file (**Q8**).
3. *Rename `12_USED_adittional_checks` → `12_USED_adittional_checks_for_supplementary`* → **done**. The one
   hard-coded path (check 01's script) was updated, and a re-run in scratch reproduced every output byte for byte
   (only the run date in `results_numbers.txt` differs). Paths updated in the folder READMEs, `supplementary.md`,
   this file, `10_USED_Paper_writing/CLAUDE.md`; the manuscript's TODO 215 and its hidden note get the new path in
   item SUPMAT1. Dated records (check 01's run log, `MM9_edits.json`) keep the old name.
4. *Raw-phenotype checks for the supplement (boxplots, correlations, distributions), not cited in the text* →
   planned as check **02** in that folder; three points to settle first in **Q12**. Their outputs go to a separate,
   unnumbered supplement (`supplementary/unreferenced/`), so TAG's citation rule is not broken.

---

## Q1. The final list and its order (the main decision)

Numbered by **first citation in the manuscript, from the title down**. Every Online Resource is cited at least once
(TAG: "the text must make specific mention of the material as a citation"). The TODO IDs are the placeholders that
get the number; "first cited" is the section of the first one.

| # | file | content (sheets / pages) | first cited → also cited | sources (read-only) |
|---|---|---|---|---|
| **1** | `ESM_1.xlsx` | **Sampling sites, accessions and site environments.** Sheet *Sites*: 30 sites (site 04 flagged "not analyzed"): number, name, ecotype region or intermediate zone, latitude, longitude, accessions analyzed, and the 10 environmental variables. Sheet *Accessions*: the 290 accessions: accession ID (= ENA sample alias, e.g. `01_03`), site, region, and the **genotype BLUPs of the 10 traits**. Sheet *Environment correlations*: the 10 × 10 Pearson matrix over the 29 sites (elevation–March temperature and silt–sand marked). | M&M §1 ¶1 (TODO 218) → M&M §4 ¶5 (217); M&M Data availability (234); [M&M §4 ¶2, BLUPs — see Q2] | `00_THIN_…/01_…R` (SITE_NAME_MAP, SITE_REGION_MAP); `all years barley.csv` (longitude, latitude); `outputs/subsection_1/tables/BLUP_10_traits.csv`; `C2_trait_environment_correlations/tables/C2corr_site_level_dataset.csv`; `merged_Evgeny_env_fitness_environmental_30_sites.csv` (elevation, silt); `C2_environmental_predictor_collinearity/tables/C2col_environment_correlation_matrix.csv` |
| **2** | `ESM_2.xlsx` | **Quality of the NIR data.** Sheet *NIR calibration*: Pearson r of NIR predictions vs laboratory values (see Q5). Sheet *Blocks 2019–20*: r and n for the 4 traits × 6 block pairs, raw values. | M&M §3 (214) → M&M §4 ¶1 (215) | Pintel MSc thesis Table 1 (p. 18) — see Q5; `12_USED_adittional_checks_for_supplementary/01_block4_…/results/tables/Table_OnlineResource_2019_block_correlations.tsv` |
| **3** | `ESM_3.xlsx` | **Variance components and H²** of the 8 traits (already built; replaces Table 1). | M&M §4 ¶3 (238) → Results Ch. 1 ¶1 (320), Ch. 2 ¶2 (322) | `supplementary/tables/ESM_variance_components.*` (script exists) |
| **4** | `ESM_4.pdf` | **LD decay**: **a** genome-wide, **b** per chromosome (the existing vector PDFs, unchanged). Caption carries r² ≥ 0.001 within 2.5 Mb, 1-kb bins, LOESS span 0.10 weighted by pair count, r² = 0.2 at ~188 kb. **Per-chromosome panel checked (your Q1 point): legitimate and current** — same run as the genome-wide curve (PLINK v1.90b6.4, 7,110,996 SNPs = the GWAS set, 290 accessions; `logs/plink_chr*.log`), and both figures are drawn by the same script from the same per-chromosome tables (`03_make_publication_figure.R`, 2026-06-10); the genome-wide profile is the sum of the per-chromosome bins. So both panels go in. | M&M §5 ¶3 (219) only | `02_USED_LD_decay_V2_wholegenome/results/figure_ld_decay_{genomewide,perchromosome}.pdf`, `ld_decay_summary_for_results_chapter.txt` |
| **5** | `ESM_5.xlsx` | **The 36 loci** (sheet 1: trait, lead SNP, span, members, lead MAF, β, SE, *P*, …) and **the 52 significant SNPs** (sheet 2: locus, MAF, β, SE, *P*). Sheet 2 answers your comments 314/315 ("supplementary table of each SNP and MAF"); the locus coordinates answer 316 (co-localized regions). | M&M §6 ¶1 (220) → Results Ch. 2 (301; TODOs 309, 310 check its content) | `01_…/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv`, `Table_loci_members_full.tsv`; `00_snp_level/significant_snps.tsv` (β = minor allele, both files agree) |
| **6** | `ESM_6.xlsx` | **All 55 candidate genes** (already built; replaces Table 2 and the 23-gene table). | M&M §6 ¶2 (222) → M&M §7 (224), §8 (227); Results Ch. 3 (303, 323–325) | `supplementary/tables/ESM_candidate_genes.*` (script exists) |
| **7** | `ESM_7.xlsx` | **Choice of ε**: summary for ε 0.2–1.0 (genes testable, accessions assigned, median groups) and the per-gene table (55 genes × 5 ε). | M&M §7 ¶2 (223) only | `04_…/loci_LDspan_eps06_V4/Diagnostics/epsilon_choice_supplementary.tsv`, `epsilon_choice_per_gene.tsv` |
| **8** | `ESM_8.xlsx` | **Canonical (1,3;1,4)-β-glucan genes** on MorexV3: reference accession, MorexV3 gene, position, identity / coverage, reciprocal best hit, tBLASTn agreement, and distance to the nearest β-glucan lead SNP (current 36 loci only). | M&M §8 ¶3 (228) → [Discussion §3, TODO 410 — see Q3] | `05_…/09_USED_canonical_betaglucan_gene_check/results/tables/canonical_betaglucan_gene_distances.tsv` (mapping columns only; its distances are stale) + `canonical_gene_distances_current_loci.tsv` |
| **9** | `ESM_9.pdf` | ~~7H MGmin × ε grid~~ **dropped (your Q1)**. **Replaced by: crosshap plots of the four presented genes**, at the setting used: *GPAT6*, *GH17*, *PHT4;3* (MGmin 2, ε 0.6) and the GDSL gene (MGmin 3, ε 0.9). One page each = page 5 of each combined PDF (the full crosshap plot; page 6 is the violin). See **Q10**. | M&M §8 ¶4, 7H (229) → places in Q10 | `04_…/loci_LDspan_eps06_V4/CombinedPDF/{fiber,starch}/<gene>/CrosshapTree+Violin__<gene>__MGmin2__Eps0p6.pdf` p. 5; `03_01_7H_…/results/figures/crosshap/CrosshapTree+Violin__7HG0729030__<trait>__MGmin3_eps0.9.pdf` p. 5 |
| **10** | `ESM_10.xlsx` | **Elite cultivars.** Sheet *Cultivars*: the 5 lines (BioSamples accession, year, breeder, country) with call rates per gene (GPAT6, GH17, PHT4;3, pooled, and the 7H gene). Sheet *Screen*: the 136 spring elite lines. Sheet *Site accounting*: per gene (4 genes) wild SNPs, elite records, shared, wild-only, elite-only, triallelic, … Sheet *Removed sites*: site by site. | M&M §10 ¶1 (232, then 231) → Results Ch. 3 (306, was "Online Resource 3") | `07_…/results/tables/Table_elite_lines.tsv`, `Table_elite_line_screen.tsv`, `Table_site_overlap.tsv`, `Table_removed_sites.tsv`; `03_01_7H_…/results/tables/` same four |
| **11** | `ESM_11.xlsx` | **Geographic origin of allele and haplotype carriers — tables**: T02 lead-SNP tests (+ T01 per region), T04 haplotype-group tests (+ T03 per region), T08 carriers vs their own site-mates. Caption carries the BH groups, the χ² region test and the region fall-back (TODO 237). See Q7. | M&M §11 (237) → Discussion §3 ¶2 (419, was "X"), §5 (420, was "X") | `11_USED_…/results/tables/T01–T04, T08`; captions `results/CAPTIONS.md` |
| **12** | `ESM_12.pdf` | **Same analysis — figures**: Fig. S1 carrier origin, Fig. S2 desert enrichment (existing 600-dpi images, unchanged). Cited together with 11: "(Online Resources 11, 12)". | same places as 11 | `11_USED_…/results/figures/Fig_carrier_origin.png`, `Fig_desert_enrichment.png` |
| **13** | `ESM_13.pdf` | **Site figure** (former Fig. 2a; already built as PDF). | Results Ch. 1 ¶3 (321) → Discussion §1 (421), §2 ¶1 (422) | `08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf` |

**Merges I made, and why:**
- **Sites + ENA aliases + genotype BLUPs + environment (1):** all are per-site or per-accession descriptors of the
  same panel; the BLUPs need a home that is cited in the text, and the site table is it.
- **NIR calibration + block 4 (2):** both are data-quality checks of the NIR traits. If the calibration values are
  later withdrawn (Q5), one sheet goes and nothing is renumbered.
- **Elite lines + site accounting (10):** same step (07 + the 7H branch), cited in the same M&M paragraph.
- **Tables vs figures of the carrier-origin analysis (11, 12):** one analysis, two files, because TAG wants tables
  as .xlsx and figures as PDF. *Alternative:* drop the two figures (the Discussion rests on the tables).

**Not merged:** the site figure (13) and the site table (1): different formats (PDF vs xlsx), and the figure is first
cited in the Results. ε (7) stays apart from the gene table (6): a different question (parameter choice).
**Not included** (not cited anywhere, and TAG requires a citation): imputation summary, variant-calling tool list,
PC variance table, per-gene crosshap figures, T05–T07 and T09–T11 of step 11.

**Answer Q1:** *(you, chat 2026-10-01)* ESM 1: check that all the information is there first → **checked, see Q11**. ESM 2: copying Noam's thesis table is fine; you may find more calibration details / Excel tables → see Q5. ESM 4: genome-wide is enough unless the per-chromosome figure is 100% legitimate and current → **it is (row 4), so both panels go in**. ESM 9: drop the grid; instead the crosshap plots of the four genes (row 9, Q10). Everything else: agreed.


---

## Q2. Cite the genotype BLUPs where they are defined?

ESM 1 is already cited (sites, TODO 218), so TAG is satisfied either way. But a reader looking for the phenotype
data will look where the BLUPs are defined. **Recommendation:** add "(Online Resource 1)" to M&M §4 ¶2:
"The genotype best linear unbiased predictions (BLUPs) of this model were used as the phenotypes in all downstream
analyses (Online Resource 1)." The Data availability draft (`build/DECLARATIONS_HANDOVER.md`) also gets "Online
Resource 1".

**Answer Q2:** *(you)* Agreed: after the files are built, add "(Online Resource 1)" to M&M §4 ¶2 (and the Data availability draft).


---

## Q3. Discussion: canonical β-glucan genes (TODO 410 and your comment 414)

TODO 410 asks whether the 83.5 Mb distance goes to the Results or to an Online Resource; comment 414 (yours) asks
"where to add it?". The M&M now cites ESM 8 for these distances. **Recommendation:** cite it in the Discussion
sentence ("…*HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP (Online Resource 8), and…") and delete TODO 410 and
comment 414 (answered).

**Answer Q3:** *(you)* Agreed, with the "supmat" comment rule (Q13): the Discussion sentence gets "(Online Resource 8)" and a supmat comment; TODO 410 is resolved.


---

## Q4. Your own Word comments that the Online Resources answer

| comment | text (short) | answered by | recommendation |
|---|---|---|---|
| 314, 315 | "add supplementary table of each SNP and MAF?" | ESM 5 sheet 2 (52 SNPs) | delete |
| 316 (reply to TODO 310) | "Add supplementary table the loci co-localized regions?" | ESM 5 (all locus coordinates) | delete |
| 317 | "Add supplementary table of the genes?" | ESM 6 | delete |
| 318 (reply to TODO 303) | "Add supplementary table" | ESM 6 | delete |
| 319 (reply to TODO 306) | "…table that shows the SNPs that got into the analysis…" | ESM 10 (site accounting + removed sites) | delete |
| 414 (reply to TODO 410) | "where to add it?" | ESM 8 (Q3) | delete |
| 311, 312, 313 | how to report the 4H / 1H signals | not an ESM question (answered in the text by R2; you kept them) | **keep**; their TODO parent 309 is deleted (resolved: ESM 5 lists both 4H loci), so they stay as stand-alone comments |

**Answer Q4:** *(you)* Do not delete them: mark each with the supmat rule, "resolved ✅" / "not resolved ❌" (Q13).


---

## Q5. NIR calibration values (TODO 214): the only source is the Pintel MSc thesis

The M&M says the calibration "was validated by Pearson correlation … (Online Resource N)". The only values are in
the thesis, Table 1 (p. 18; the thesis's "Season 2020" = our 2019–20): kit vs calibration scans (starch 0.785,
protein 0.915, β-glucan 1.000), kit vs 2019–20 scans (0.466, 0.556, 0.514), 2019–20 scans vs calibration scans
(0.548, 0.649, 0.514; fiber 0.526). No *n*, no *P* thresholds; nothing in the project files.
- **A (recommended):** put these values in ESM 2 (sheet *NIR calibration*), marked in the sheet as taken from the
  calibration report of this study; **keep TODO 214 open**, reworded to what is missing: N. Pintel to confirm the
  values, the number of accessions per comparison and the significance levels. The β-glucan r = 1.000 gets a note
  (the calibration was built on those kit values).
- **B:** leave the sheet empty until Noam answers (ESM 2 then holds only the block-4 sheet).
- *Your note (Q1): you may find more calibration details or Excel tables.* If they give the number of accessions per
  comparison (and ideally the paired NIR and laboratory values), I build the sheet from them instead, and TODO 214
  closes. Please put them in `10_USED_Paper_writing/section_materials/02_materials_methods/` and tell me the file name.

**Answer Q5:** *(you, chat)* You will bring more information on the NIR calibration in a later message. **Pending:** ESM 2 is built with the block-4 sheet now; the calibration sheet waits for your files.


---

## Q6. ~~7H grid: show the haplotype-test columns?~~ — no longer needed

The grid was dropped (your Q1). TODO 229's "decide whether to show kw_p, eta2, delta" is moot; TODO 229 now points
to the crosshap plots (ESM 9, Q10).

**Answer Q6:** *(you)* A (moot: the grid was dropped).

---

## Q7. Carrier-origin tables (ESM 11): columns the text does not use

T02 and T04 carry an accession-level desert odds ratio (`desert_OR`), and T08 a check "without the two most extreme
accessions". The M&M (item MM7) left both out of the text. **Recommendation:** drop these columns from the ESM, so
the file matches the method described. Everything else as in the source tables (definitions on the caption sheet).

**Answer Q7:** *(you)* Recommendation accepted: `desert_OR` and the extreme-accession check are dropped from ESM 11.


---

## Q8. Title block of every file: authors and corresponding author

Not in any project file. I will use placeholders until you give them:
- author names, in order: `[AUTHOR NAMES]`
- corresponding author: name, affiliation (institution, department, city, country), e-mail: `[AFFILIATION]`, `[E-MAIL]`

Layout: each **.xlsx** opens with a sheet *About* (title block, "Online Resource N", the caption, the list of sheets,
the column definitions and the sources); the data sheets start with the header row, so they stay machine-readable.
Each **PDF** opens with the title block and the caption on page 1, the figure(s) after.

**Answer Q8:** *(you)* Placeholders for now, and **add it to the open tasks at the top of the manuscript** → item CHK1 (checklist item 1 names the ESM title blocks).


---

## Q9. Where the ESM captions go in the manuscript

TAG: "a concise caption for each supplementary material" in the manuscript; no position is prescribed.
**Recommendation:** a new last heading of the working file, **"Supplementary Information"**, after the Conclusions
(the References and Declarations will follow in their own session), with one entry per file:
"**Online Resource 1** Sampling sites, … (ESM_1.xlsx)". I show you the 13 captions before applying them.

**Answer Q9:** *(you)* "Show them to me. Are you sure the captions go in the manuscript and not in the supmat files? That's odd, and it will be long text for the main manuscript."
  **Reply:** yes, TAG asks for both, but the manuscript part is short:
  - TAG guidelines, "Supplementary Information (SI)" (`journal_TAG/TAG_submission_guidelines_2026-09-15.pdf`, pp. 38–39): "For each supplementary material, please supply a **concise caption** describing the content of the file", and under Accessibility: "**The manuscript contains a descriptive caption** for each supplementary material". Springer prints these captions at the end of the published article, next to each file's download link.
  - **Concise** = one or two lines per file, e.g. "**Online Resource 4** Linkage disequilibrium decay across the genome (a) and per chromosome (b) (PDF)". 13 files ≈ 250–350 words at the very end of the manuscript, not in the body text.
  - **The long explanation** (column definitions, settings, sources, figure legends) goes **inside each file** (the *About* sheet / the PDF's first page), with the title block.
  I show you all 13 short captions (and the in-file texts) before anything goes into the manuscript.


---

## Q10. ESM 9 — the crosshap plots of the four presented genes

Checked: page 5 of each combined PDF is the full crosshap plot (bar plot of haplotype sizes, marker-group matrix,
allele-count bars, phenotype by haplotype, the MG table); page 6 is the violin. Group sizes on page 5 match the text
(*PHT4;3*: A 137, B 41, C 17, D 10; GDSL: A 212, B 34). Page 5 is copied unchanged (vector PDF, 10 × 7 in).
- **(a) GDSL gene:** its crosshap PDF exists twice, for fiber and for starch. The haplotypes are identical; only the
  phenotype panels differ. **Recommendation, as you said ("total 4 plots"):** the fiber page (the stronger GWAS
  signal). *Alternative:* both pages (5 in total).
- **(b) Where it is cited** (it needs at least one citation; first citation keeps it number 9):
  1. M&M 7H paragraph, TODO 229's place: "…at which the three significant SNPs were included in the marker group
     that defined the haplotypes (Online Resource 9)" — required, replaces the grid;
  2. **Recommendation, also:** Results Ch. 3, the three gene paragraphs: "(*q* = …, η² = …; Online Resources 6, 9;
     Fig. 3a)", and Results Ch. 4 ¶3: "…with the same grouping for both traits (Fig. 4; Online Resource 9)".

**Answer Q10:** *(you, 2026-10-01)* Both GDSL plots; cite in the Results and M&M where needed, citations only. → **Done:** ESM 9 has 5 plots (GDSL fiber d, starch e). Item **E2** added only citations: Results Ch. 3 GPAT6, GH17, PHT4;3 "(… Online Resources 6, 9; Fig. 3a–c)" and Ch. 4 ¶3 "(Fig. 4; Online Resource 9)". M&M: already cited at the 7H paragraph; an earlier M&M citation would renumber 7–9, so none added. No wording changed.


---

## Q11. ESM 1 — check of the information (your Q1 point): what is there, and four points to confirm

**Everything needed is in the project files** (checked today):
- Sites: number, name (`SITE_NAME_MAP`; Hachola for site 04 from the 30-site file), region / zone (`SITE_REGION_MAP`),
  one latitude/longitude per site (`all years barley.csv`; each site has exactly one pair).
- Accessions: the 290 genotype BLUPs of the 10 traits (`BLUP_10_traits.csv` = the GWAS phenotype files); their site
  means equal the inputs of Fig. 1d (max difference 1 × 10⁻⁹).
- Environment: the 8 variables of Fig. 1d are identical in `C2corr_site_level_dataset.csv` and the 30-site file;
  elevation and silt (dropped for collinearity) come from the 30-site file; the 10 × 10 matrix from `C2col_…`.
- **ENA aliases verified today against ENA's public API** (read-only query, project PRJEB79623): all 290 accession
  IDs (`01_03` …) are ENA sample aliases, one run each; the other 10 ENA samples are site 04 (`04_01` … `04_11`).

To confirm:
- **(a) ENA accession numbers:** add the ENA sample accession (SAMEA…) and run accession (ERR…) next to the alias?
  **Recommendation: yes** — it links each accession to its reads unambiguously. The ENA report is saved as a dated
  source file in `supplementary/sources/`.
- **(b) Units not stated in the project:** plant height (tallest tiller length) and grain weight. From the values
  (means 85 and 2.1) they look like **cm** and **g**; Mol Ecol gives cm for spike length only. Confirm, or I write
  "units as recorded" and list it as open (TODO 209 already asks how grain weight was measured).
- **(c) BLUP scale:** the BLUPs are on the season-centered scale (deviation from the season mean, in trait units),
  as used in all analyses. The sheet says so. OK?
- **(d) Electrical conductivity:** *(you)* "just write one of the units; do not list it as open" → **mS/cm**, the unit
  of the source table and of Fig. 1d (µS/cm would contradict Fig. 1d). Not listed as open.

**Answer Q11:** *(you)* a OK (ENA accessions kept); b "as in Mol Ecol" → **done:** no unit for plant height and grain weight (ESM 1, ESM 3, extra supplement); d OK (mS/cm).
  **c — what I meant (please confirm):** the BLUPs in ESM 1 are **not grain contents**. Before the model was fitted, each
  season's mean was subtracted from that season's values (to remove the season differences), so each BLUP is the
  accession's genetic value as a **deviation from the average accession**, in the units of the trait. Example: starch BLUP
  +8.1 for accession 01_03 = about 8.1 percentage points more starch than the average accession, not "8.1 % starch". These
  are exactly the numbers used in every analysis (GWAS, haplotype tests, correlations) and plotted in ESM 13 ("BLUP
  (centered)"). The alternative would be to add the overall mean back so the numbers look like contents (e.g. starch around
  35 %), but those would not be the values analyzed, and that step is not in the M&M. **Recommendation: keep them as they
  are** (the sheet already says "deviation from the season mean"). → write "OK" or tell me otherwise.


---

## Q12. The raw-phenotype checks (new check 02 in `12_USED_adittional_checks_for_supplementary/`)

Plan: `02_raw_phenotype_checks/` with `scripts/`, `results/{figures,tables}/`, `logs/`, README; figures to TAG
specs (174 mm wide, 600 dpi, Arial-metric 9 pt, panels **a, b, c**, no titles in the image; PDF + PNG). Input
`all years barley.csv`, read-only. Three points to settle:
- **(a) Which plants.** **Recommendation:** the analysis panel (290 accessions; site 04 and the cultivated checks
  removed; block 4 of 2019–20 removed), raw values *before* the ±3 SD outlier step, with the ±3 SD limits drawn and
  the values beyond them counted. Exception: the per-block boxplots show block 4 of 2019–20, marked "excluded".
  **Starch 2021–22:** shown as recorded, *without* the +8.85 the scripts add (the M&M does not describe it; TODO 213).
  Note that the raw season means then differ strongly (starch 33.0 / 42.9 / 20.7%), and a reader will see it.
  *Alternative:* show starch 2021–22 with the +8.85, labelled.
- **(b) Your correlation items 2a and 2b — which correlations?** Two readings:
  - **(i) between seasons, per trait:** accession means of each season against each other season (2b), and each
    season against the mean over all seasons (2a). Together with 2c (between blocks) this is the repeatability
    check.
  - **(ii) between traits:** the 4 × 4 correlations of the nutritional traits, from accession means over all
    seasons (2a) and within each season (2b).
  **Recommendation: both** (each is one small figure + table, all descriptive).
- **(c) Distributions:** the 4 nutritional traits and the 6 morphological and phenological traits of Results Ch. 1
  (flowering time, plant height, spike length, number of tillers, grain weight, number of grains): raw values by
  season, and the accession BLUPs used in the analyses. **Recommendation: both.** (The raw values include extreme
  entries, e.g. plant height up to 664 and spike length up to 78; they show what the ±3 SD step removes.)

> **Built meanwhile with the recommended defaults (easy to change):** check 02 in
> `12_USED_adittional_checks_for_supplementary/02_raw_phenotype_checks/` (scripts, Figs. E1–E4, E6, E7 + tables T1–T5,
> README, `results_numbers.txt`); Fig. E5 = check 01's block correlations. Assembled as the unnumbered supplement
> `supplementary/unreferenced/ESM_extra_phenotype_checks.{pdf,xlsx}`.
> **New fact for (a): raw starch 2021–22 has 25 of 1,163 values below 0 % (minimum −21.6 %); 4 stay below 0 even with the
> +8.85, and the ±3 SD rule removes only those below about −6 %.** Showing raw starch makes these visible to a reader of
> the supplement. Your call: keep as recorded (current), show with the +8.85 (labelled), or leave starch 2021–22 out of
> the extra supplement. (Listed in the hand-over as something to check with N. Pintel / the NIR calibration.)

**Answer Q12:** *(you)* a: use the starch data with the +8.85 for everything in the supplement, **do not report it** → **done** (check 02 re-run; starch 2021–22 mean now 29.6 %; no file mentions it — checked); b both; c both (both were already built).


---

## Q13. The "supmat" Word comments (your broad comment 1)

Format (one Word comment, author **supmat**, beside the existing TODO/your comment on the same text):
- `supmat: not resolved ❌ — <what is missing>`
- `supmat: resolved ✅ — Online Resource N (ESM_N.xlsx)`

Placed on: every Online Resource citation (29 placeholders + TODO 234's sentence) and every place that still needs
one: M&M §4 ¶2 (BLUPs, Q2), Discussion §3 (*HvGlb1* distance, Q3), the crosshap places (Q10), and your comments
314/315 (SNP table), 316 (co-localized loci), 317 (gene table), 318, 319, 414 (Q4). IDs from 239 (M&M), 326
(Results), 423 (Discussion). All start as ❌ and turn ✅ when the file is built and the number is in place.
**Two citations I propose for the places your comments mark** (no citation there today):
Results Ch. 2 ¶3 (lead-SNP MAF and effects; comments 314/315) → "(Online Resource 5)";
Results Ch. 3 ¶1 (30 of 55 genes tested; comment 317) → "(Online Resource 6)".
They are deleted with the other comments before submission (checklist item 5 says so).
**Done: item SUPMAT1 applied (36 comments, all ❌; preview `SUPMAT1_PROPOSED.md`). The Results crosshap places (Q10b) and the two citations below get theirs when you answer Q10 and Q13.**

**Answer Q13:** *(you)* OK and OK → **done in E2:** Results Ch. 2 ¶3 "(median MAF = 0.068; Online Resource 5)"; Ch. 3 ¶1 "…distributed over 12 of the 36 loci (Online Resource 6)." Supmat 331 and 333 → ✅.


---

## Q14. The supplement that is not cited (your broad comment 4)

Everything not cited in the text goes to `supplementary/unreferenced/`, unnumbered, so the numbered set stays TAG
compliant: for now the raw-phenotype checks (Q12), as one PDF of figures + one .xlsx of tables, each with the title
block and captions, named `ESM_extra_phenotype_checks.{pdf,xlsx}`. If you later move one into the main set, it gets
a number and a citation. OK?

> **Done as proposed:** `supplementary/unreferenced/ESM_extra_phenotype_checks.pdf` (title page + Figs. E1–E7, each with its
> caption) and `.xlsx` (About + six tables), same title block; built by `supplementary/scripts/build_ESM_extra.py`.

**Answer Q14:** *(you)* OK.


---

## Q15. The 13 captions for the manuscript (item CAP1, not applied yet)

A new last section, **"Supplementary Information"**, after the Conclusions, with one short caption per file (about 270 words in all). The long legends and column definitions are already inside each file. Same text, with the change marked, in `CAP1_PROPOSED.md`.

> **Online Resource 1** Sampling sites (ecological region, coordinates and environmental variables), the 290 accessions (European Nucleotide Archive identifiers, site and genotype BLUPs of the ten traits) and the correlations among the environmental variables

> **Online Resource 2** Pearson correlations between the blocks of the 2019–20 season for the four nutritional traits

> **Online Resource 3** Variance components and broad-sense heritability of the four nutritional and four morphological and phenological traits

> **Online Resource 4** Decay of linkage disequilibrium with physical distance across the genome (**a**) and per chromosome (**b**)

> **Online Resource 5** The 36 loci of the genome-wide association analysis and the 52 genome-wide significant SNPs

> **Online Resource 6** The 55 candidate genes: position, haplotype tests, functional annotation and the three genes carried forward

> **Online Resource 7** Number of testable genes and proportion of accessions assigned to haplotype groups at five values of the crosshap clustering radius (ε)

> **Online Resource 8** The ten canonical (1,3;1,4)-β-glucan synthase and endohydrolase genes located on MorexV3, and their distances to the β-glucan lead SNPs

> **Online Resource 9** crosshap haplotype plots of *GPAT6*, *GH17*, *PHT4;3* and the GDSL esterase/lipase gene

> **Online Resource 10** Elite cultivars: selection and call rates, and the SNP sites compared between the wild and elite call sets

> **Online Resource 11** Geographic origin of the carriers of the lead-SNP alleles and of the haplotype groups of the four presented genes: tests of desert enrichment and regional distribution, and comparison with site-mates

> **Online Resource 12** Geographic origin (Fig. S1) and desert enrichment (Fig. S2) of the lead-SNP allele carriers and haplotype groups

> **Online Resource 13** Genotype BLUPs of the four nutritional traits across the 29 sampling sites, by ecological region

- **(a)** Approve the captions, or tell me what to change.
- **(b)** Place: after the Conclusions (recommended), or after the References / Declarations once they exist?

**Answer Q15:** *(you)* OK; asked whether we follow TAG (a published TAG paper as example). → **Checked** (the linked Springer page blocks scripted access, so three open-access TAG GWAS papers of 2026 were read in PubMed Central: PMC13616022, PMC13546169, PMC13219082). Each ends with a **"Supplementary Information"** section built by Springer from the files and their captions ("Below is the link to the electronic supplementary material", one line per file), and has a **Key message**; back matter: Acknowledgements, Author contributions, Funding, Data availability, Declarations (conflicts of interest). Our manuscript matches; the back matter is the Declarations session. Citation style in those papers varies ("Supplementary Table S1" or "Online Resource"); we follow the guidelines ("Online Resource N"). → **CAP1 applied** (+ CHK4); at submission, enter the same captions in the submission system.


---

## Q16. ESM 5 (loci) — two leftovers from your note of 2026-09-16

`section_materials/03_results/NOTE_loci_table_for_supplementary.md` left two points open:
- **(a) Carrier count:** a column with the number of accessions carrying the minor allele at each lead SNP (out of those
  with a call), "more intuitive than MAF alone". It is already in a project output, so no new analysis: step 11, T02
  (`n_minor`, `n_major`, `n_missing`). **Recommendation: add "Minor-allele carriers" and "Accessions called"** next to
  the lead MAF.
- **(b) The four clumping-diagnostic columns** (largest internal gap, clump span and members before the gap rule, members
  removed by it): "probably supplement-only, possibly dropped". They are in ESM 5 now, defined on the About sheet.
  **Recommendation: keep** — they show what the 50-kb gap rule did to each locus.

**Answer Q16:** *(you)* OK and OK → **done:** ESM 5 has "Minor-allele carriers" and "Accessions called" (from step 11, T02; they reproduce the lead MAF exactly); the four clumping-diagnostic columns stay.


---

## Q17. Your NIR.zip (Noam Pintel's files): what is in it, and what to put in ESM 2 (closes Q5)

**What I found (all checked against the raw files and our data):**
1. **The raw data behind thesis Table 1** (`NIR/NIR2020.vs.NIRcali/`: `kits.results.xlsx`, `calibration.scans.xlsx`,
   `season2020.xlsx`, Noam's script `kits.cali.2020.R`). These are **30 calibration accessions, one per site**
   (site 04 included). For each, there are three values:
   - **Lab values.** Protein by Kjeldahl, 3 subsamples of about 0.2 g, N × 6.25. Starch and β-glucan with the Megazyme
     kits, mean of 4 subsamples (one accession has 2 for starch).
   - **"Calibration scans".** Mean of 5 DA 7250 scans of the same grain, 7 Jan 2020. These were predicted with
     Perten's **factory cultivated-barley calibrations** (`Barley_Protein/Starch/Fiber_Asis_7July2011`, bias 0,
     slope 1, i.e. before any adjustment). That set has **no β-glucan calibration**.
   - **The 2019–20 value of the same accession**: one plant (block 1), identical to our data for all 30 accessions.
2. **Every thesis r reproduces exactly** (n = 30): protein 0.915 / 0.556 / 0.649; starch 0.785 / 0.466 / 0.548;
   β-glucan 1.000 / 0.514 / 0.514.
3. **Two thesis values are not validations:**
   - **β-glucan, kit vs calibration scans: r = 1.000.** The "calibration scan" column is a copy of the kit values,
     because there was no β-glucan calibration in January 2020. The r = 1 holds by construction.
   - **Fiber, r = 0.526** (2019–20 vs calibration scans). This comes from `Book3.csv`, in which 5 of the 30 rows differ
     from the source data. From the instrument scans and our 2019–20 data, **r = 0.25 (P = 0.19)**. Fiber has no lab
     value in any case.

**Proposal for ESM 2.** Add two sheets next to the block correlations, and the caption gains "and the validation of the
near-infrared calibration":
- **"NIR calibration".** The 30 accessions, with their lab values, calibration-scan predictions and 2019–20 values for
  protein, starch and β-glucan.
- **"Calibration correlations".** Pearson r, n and P:
  - protein and starch: lab vs calibration scans, lab vs 2019–20, and 2019–20 vs calibration scans;
  - β-glucan: lab vs 2019–20 only.

The options:
- **A (recommended):** as above. It leaves out the β-glucan r = 1 (by construction) and the fiber value (no lab
  reference; the thesis number does not reproduce). TODO 214 closes, and supmat 240 → ✅ (the last ❌).
- **B:** reproduce thesis Table 1 as it is, with notes. Not recommended, because the fiber 0.526 does not hold.

**Also found, for the M&M (outside the supplement; nothing changed):**
- **Kjeldahl:** protein = N × **6.25** (0.1 N HCl titration). TODO 212 says this factor is missing.
- **Lab subsamples:** protein 3, starch and β-glucan **4** (one accession 2). The M&M says "three subsamples of about
  200 mg".
- **Calibration scans:** 5 scans per accession, averaged; the names of the factory calibration files are above. This is
  the "manufacturer's calibration for cultivated barley" that the M&M says was then adjusted (TODO 213).
- **TODO 213 stays open.** The 2021 season scans (Oct 2021, `NIR season 2021.xlsx`) carry no calibration names, so
  whether one calibration served all three seasons is still unknown.

- **(a)** ESM 2: option A or B?
- **(b)** Should I add the M&M facts above to TODO comments 212 and 213, for the M&M session? This changes the comment
  text only, not the M&M.

**Answer Q17:** *(you, 2026-10-01)* Option A; fix the M&M facts in the text; (b) yes, add them to TODOs 212 and 213; TODO 213 stays open (S. Hübner to inspect). → **Done:** ESM 2 has the sheets "NIR calibration" and "Calibration correlations" (sources copied to `supplementary/sources/NIR_calibration/`; every thesis r reproduced, P checked against R). Item **E3**: M&M "scanned five times", "four subsamples", "Kjeldahl method (nitrogen × 6.25) in three subsamples"; TODO 212 (citation still open) and TODO 213 (open) carry the facts; TODO 214 removed; caption of Online Resource 2 extended; checklist ✅. **All 37 supmat comments are ✅.**


---

## Also noticed (no answer needed unless you disagree)

- Electrical conductivity units: moved to Q11(d).
- **LD-decay per-chromosome PDF** has a title inside the image ("LD decay by chromosome (wild barley, n = 290)").
  TAG forbids titles inside *main-text* figures; ESM is published as received. I use the existing PDF unchanged
  (step 02 is read-only); say if you want it redrawn.
- Placeholders "Online Resource 1" (TODO 301) and "Online Resource 3" (TODO 306) in the Results are old numbers;
  they become 5 and 10.
