# Materials and methods — blueprint

Plan for the M&M of the TAG manuscript. **Not the text.** Built 2026-09-29 with the user; tuned item by item before
the writing session. Paths are relative to `barley_for_publication_project/` unless they start with `new_publishing_paper/`.

**Status (2026-09-30): WRITTEN.** The M&M was drafted and approved section by section (as `02_materials_methods.md`),
then merged the same day with the Results, Discussion and Conclusions into the combined working file
`new_publishing_paper/Methods_Results_Discussion_Conclusions.md` (reading copy without comments:
`Methods_Results_Discussion_Conclusions_reading_copy.md`, regenerated with `build/make_reading_copy.py`).
Where the text differs from this plan, the text and its hidden `<!-- src -->` notes are right (decisions taken while
writing: e.g. §3 no as-is basis; §5 PLINK = Purcell et al. 2007 only; §7 ε justified on the current 55 genes;
§9 shortened; §11 compact, covering only what the Discussion uses; §12 closing paragraphs + LLM statement). Every open
item is now a TODO Word comment (IDs 201–237) in that file. This blueprint is kept for the record.
**2026-09-30, later:** the M&M was then rewritten in E. Potapenko's voice (review items MM1–MM11, log at the end of
`build/REVIEW_2026-09-27_Results_Discussion.md`); the text in the working file supersedes this plan's wording and order
(e.g. §9 is now the last paragraph of §8; §12 became two subsections).

*Status (2026-09-29, superseded):* reviewed with the user section by section, §1–§10 approved; §11 is a placeholder; open items below
are settled together with the user in the writing session.

**Open items**
- **User decisions:** ε numbers, (a) without numbers or (b) re-run on the current genes (§7 ¶2); a positive reason for
  MGmin = 2 (§7 ¶2); code and MorexV3 VCF availability (§12); LLM-use statement, **last** (§12).
- **To generate (user):** 2019 inter-block correlations of the NIR traits (§4 ¶1). (The carrier-origin analysis, §11, was finished 2026-09-30.)
- **Ask Evgeny:** MorexV3 citation (Mascher et al. 2021, not Monat et al. 2019) (§1 ¶3); missingness threshold and filter
  order; whether the V3 run re-trimmed the reads (§1 ¶2–4); which grain-number column he used (§2 ¶2).
- **Ask Noam (via the user):** single-seed descent rounds (§1 ¶1); Hachola wording (§1 ¶4); pre-sowing protocol (§2 ¶1);
  flowering-time start point, plant height, grain weight and grain number definitions (§2 ¶2); threshing in 2021 (§2 ¶3);
  NIR calibration values, incl. β-glucan r = 1.000 (§3 ¶2).
- **Citations to check:** GATK (Poplin et al. 2017 is a preprint) (§1 ¶3); Kjeldahl (Beljkaš et al. 2010) (§3 ¶2);
  DivBrowse and the barley pangenome (§10 ¶1).

---

## 0. Read first (writing session)

**This blueprint is the rule set for the M&M session.** Where it differs from `10_USED_Paper_writing/CLAUDE.md`
(e.g. how to rank the sources), follow this file.

### Files
- **Write to:** `new_publishing_paper/02_materials_methods.md` (new). TODO Word-comment IDs start at **201**. *(Done; merged
  2026-09-30 into `new_publishing_paper/Methods_Results_Discussion_Conclusions.md`, the working file from now on.)*
- **Rules:** `10_USED_Paper_writing/CLAUDE.md` (sources of truth, writing rules, TODO convention) · `10_USED_Paper_writing/TAG_requirements.md`.
- **What the M&M must support:** `new_publishing_paper/build/UPDATED_Results_Discussion_Conclusions.md` *(superseded 2026-09-30: now inside `new_publishing_paper/Methods_Results_Discussion_Conclusions.md`)* (the working
  Results + Discussion + Conclusions). Read it in full, including the hidden `<!-- -->` notes marked "M&M MUST".
- **Hand-over lists:** `new_publishing_paper/build/MM_HANDOVER_from_review.md` (text moved out of the Results) ·
  `new_publishing_paper/build/BLUEPRINT.md` (historical; use only its "⟶ M&M" notes).
- **Supplementary plan:** `new_publishing_paper/supplementary/supplementary.md`.
- **Old M&M texts** (`10_USED_Paper_writing/old_mini_paper/`):
  - `md/04_materials_methods.md`: the mini paper's M&M. **Model for style and voice (S. Hübner).**
  - `M&M/old_long_M&M_versions_from_old_mini_paper_that_been_shorten_for_mini_paper/`:
    `Materials and methods 2.5.2026.md` (long; 69 comments, mostly citations and "verify" flags; NIR calibration table) and
    `Materials_and_methods_minipaper 4.6.26.md` (condensed). Useful for the wet lab and genotyping; **their analysis parts are retired.**
- **Reference sources** (`10_USED_Paper_writing/reference_sources/`, md + pdf):
  - Potapenko et al. 2026 *Mol Ecol*: same collection and **same common garden**, *H. spontaneum* only, WGS (ENA PRJEB79623).
  - Potapenko et al. 2026 *MBE*: same collection, *H. spontaneum* **and *H. bulbosum*** (do not mix them up); sampling design, environmental data.
  - Pintel MSc thesis: wet lab and NIR detail. **Not citable.** Its "Season 2020" = our 2019.
- **Section materials:** `10_USED_Paper_writing/section_materials/02_materials_methods/NOTE_flowering_time.md`.

### How to use the sources
- **All sources are complementary.** Use them together to write the M&M: the project files, Evgeny's two papers, the old M&M
  drafts, the mini paper, the thesis, section materials, and the user's answers.
- **The project files win.** They are this study's research: whatever the project directories can answer (analysis
  steps, parameters, versions, numbers, which data were used), they answer.
- **Where the project files cannot answer** (the plants, the common garden, phenotyping, DNA and sequencing, which are not
  in these directories), **cite Evgeny's published papers** (same collection and data; *Mol Ecol* first, *MBE* for
  sampling and environmental data).
- **Contradictions:** write it as Evgeny's paper has it, add a TODO Word comment quoting the conflicting statement(s) and
  their source (Evgeny's other paper, the thesis, the old M&M drafts, the mini paper), and raise it with the user, who
  checks it with Noam Pintel, the source of truth for the wet experiment and phenotyping.
- **Style and voice** come from the mini paper.

### Rules
- **§4 onward (the analyses): verify in the code, not in this blueprint.** For every step, read the scripts, the
  algorithms and settings they actually run, and their outputs, and write what the code does, precisely. The
  blueprint only says what to cover and where to look; its summaries can be incomplete or out of date.
- **Every number, version and parameter from the project files**, with a `<!-- src: … -->` note. Missing → TODO Word comment.
- **Nothing from MorexV2**; no retired pipeline, no comparison with it. Justifications are fine here, stated positively.
- **Voice:** the mini paper's M&M, fuller. Paraphrase every source (TAG screens for text recycling). Cite published work only.
- **TAG:** enough detail to repeat the work; headings ≤ 3 levels; accession numbers at the end of M&M; state the graphics
  software; material-sharing restrictions, if any. **LLM-use statement: last, the user handles it at the end.**
- **Online Resources:** cite as "Online Resource N" with a TODO for the number; each one is logged in `supplementary.md`.

### Tags used below
**SRC** source file · **CITE** reference · **CHECK** verify in the project files · **NOAM** contradiction or open wet-lab fact: TODO comment + raise with the user (checked with Noam)
· **SUPP** Online Resource · **MUST** statement the Results/Discussion depend on · **DROP** do not write · **DECIDE** user decision

---

## 1. Plant material and genotyping

**¶1 Collection**
- 300 accessions, 30 sites × 10 plants ≥ 10 m apart; 3 ecotypes + 3 intermediate zones; design breaks the link between
  environment and geography; coordinates recorded. **CITE** *Mol Ecol* (sampling) / *MBE*.
- Single-seed descent and selfing. **NOAM** rounds: *Mol Ecol* "three times" vs *MBE* "two rounds" → cite one, TODO on the conflict.
- **DROP** "Hübner et al. 2009, 2012, 2013" for the design (ambiguous citation); **CHECK** "established in 2017" (thesis only).

**¶2 Sequencing** — cite *Mol Ecol* **only up to the raw reads** (reference-independent)
- DNA extraction, libraries, NovaSeq 6000, ≥ 5×; raw reads ENA **PRJEB79623**. **CITE** *Mol Ecol*. **DROP** the Shermeister thesis citation.
- Read QC and trimming (FastQC 0.11.5, fastp 0.20.0) are also before alignment: **CHECK** whether the MorexV3 run started
  from the same trimmed reads (then cite *Mol Ecol*) or re-trimmed (then describe as ours).
- **Never cite *Mol Ecol* from alignment onward:** its alignment, calling, filters and SNP set are MorexV2.

**¶3 Alignment and MorexV3 variant calling** — **SRC** Evgeny's reply (Appendix A; he ran the re-call) + VCF header
(`/mnt/data/shahar/gwas_barley/data/inputs/morexV3_with_ids.vcf.gz`). Write in our words from his paragraph.
- BWA-MEM2 v2.2.1, default parameters (**CITE** Vasimuddin et al. 2019) → MorexV3.
- Duplicates marked with Picard MarkDuplicates v3.4.0 (**CITE** Picard Toolkit, Broad Institute; TAG "online" format).
- GATK v4.6.2.0, best-practices workflow: HaplotypeCaller in GVCF mode, ploidy 2 → ReblockGVCF → GenomicsDBImport →
  GenotypeGVCFs, max 2 alternate alleles; run in non-overlapping genomic intervals (genome size). **VERIFIED** in the
  header: all four tools at 4.6.2.0, ploidy 2, max 2 alternate alleles, 60-Mb intervals.
- BWA-MEM2 and Picard versions: from Evgeny only (not in the VCF).
- **MorexV3 citation:** Evgeny suggested Monat et al. 2019 (TRITEX), but that paper is the **MorexV2** assembly
  (*Mol Ecol* cites it for V2). Use **Mascher et al. 2021** *Plant Cell* 33:1888–1906 (MorexV3) → TODO comment to confirm with Evgeny.
- **GATK citation:** Evgeny cites Poplin et al. 2017, a bioRxiv preprint; TAG lists published works only → check TAG
  on preprints, or cite McKenna et al. 2010 *Genome Res* (GATK) and/or Van der Auwera and O'Connor 2020 (GATK book).
- **CHECK** trimming: Evgeny's text starts at "Cleaned reads"; whether the V3 run re-trimmed is still open (see ¶2).

**¶4 Variant filtering → 290 × 7,110,996** — bcftools v1.13 (**CITE** Danecek et al. 2021) + VCFtools v0.1.15 (**CITE** Danecek et al. 2011)
- Biallelic SNPs only, INDELs removed; hard filters: QD < 5, MQ < 45, FS > 60, SOR > 3, MQRankSum < −2.5, QUAL < 140, site
  DP < 900 or > 1,800. **VERIFIED** header: exactly Evgeny's `bcftools filter` command. Our filters differ from *Mol Ecol*'s
  (MorexV2); write ours. Give the values in the text (or an Online Resource table, as Evgeny suggests).
- Heterozygous genotypes → missing (**VERIFIED** header, `setGT`); genotypes with DP < 3 → missing (**VERIFIED** in the
  data: 0 of 122 M called genotypes in the first 500,000 SNPs have DP < 3).
- MAF ≥ 0.05 and missingness ≤ ~30% (your email; Evgeny did not answer): **VERIFIED** in the data over the 300 accessions
  (first 500,000 SNPs: min MAF 0.0500, max missingness 30.7%). **CHECK** the exact missingness threshold and filter order with Evgeny.
- **MUST** the filters were applied to the **300-accession** call set; the 10 Hachola (site 04) accessions were removed
  afterwards and the SNPs were not re-filtered. In the 290 panel, 90,151 SNPs have MAF < 0.05 (min 0.017), 155 exceed 30%
  missingness, and one lead SNP sits just below 0.05 (fiber_L13, 5H:437,005,013, MAF 0.0496). **SRC**
  `01_…/intermediates/morexV3_290_freq.frq`, `02_loci_FINAL/tables/Table_loci_master.tsv`.
- Hachola removal: write "error in the sampling site" + **NOAM** TODO to fine-tune. **CHECK** script comments still say
  "missing from the phenotyping experiments" (e.g. `00_THIN/03a`).
- **DROP** the thesis's filter sentence (MorexV2 run; "MAF 10%" typo).

## 2. Common-garden experiment and phenotyping

**¶1 Common garden**
- 3 seasons, 2019–20, 2020–21, 2021–22 (= data labels 2019/2020/2021); insect-proof net-house, Hula Valley (33°09′08.4″ N
  35°37′15.8″ E); 5-L pots; randomized complete block design, one plant per accession per block. **CITE** *Mol Ecol* once.
- Blocks 4 / 4 / **5**: state that 2021 had 5 (*Mol Ecol* says 4 replicates). **SRC** `all years barley.csv` (Block_no).
- Pre-sowing: write as *Mol Ecol* + cite. **NOAM** TODO: thesis 3% NaOCl, 16 °C for 18 d · *Mol Ecol* 4%, 4 °C for 10 d ·
  *MBE* 3%, 4 °C for 10 d then 20 °C for 14 d.
- **DROP** the cultivated checks (Morex, Clipper): not used anywhere in the paper (user, 2026-09-29); they stay only in the
  Word comment of TODO 401.

**¶2 Traits recorded** — 10 traits used: 4 NIR grain traits (§3) + 6 morphological. **CITE** *Mol Ecol* for how the 6 were
measured (all 6 are in its trait list); keep our own definitions minimal. **SRC** `all years barley.csv`; seasons read from
the data match *Mol Ecol*.

| trait (column) | seasons | used in | definition / note |
|---|---|---|---|
| Flowering time (`flowering_time`) | 3 | Table 1, Fig. 1, 2c | days from sowing to heading (Zadoks; FieldBook); 2020 range 73–147 d = *Mol Ecol*. **DROP** GDD. **NOAM** thesis "from transplanting" |
| Plant height (`Tiller.Length`) | 3 | Fig. 2c | *Mol Ecol* "tallest tiller length" (= the column name). **NOAM** thesis "ground to base of the highest flag leaf" |
| Tillers (`Number.Of.Tillers`) | 2019, 2020 | Table 1, Fig. 1 | *Mol Ecol*; thesis: counted at the end of the season |
| Spike length (`Spike.Length`) | 3 | Table 1, Fig. 1 | defined nowhere → cite *Mol Ecol* only |
| Grain weight (`Grain_weight`) | 2019, 2021 | Table 1, Fig. 1, 2c | defined nowhere → cite *Mol Ecol*. **NOAM** weight of the grain of how many spikes? |
| Grain number (`number_of_grains`) | 2019 | Fig. 2c | write as *Mol Ecol*: number of grains per spike + **CITE**. **NOAM** TODO: our variable is the total grains on 3 spikes (2 spikes for 29 plants); per plant it correlates with grains per spike at r = 0.986 (845 plants, 2019, block 4 excluded; `Grain_per_spike` = grains / spikes). **Ask Evgeny** which column he used: his Methods say "grains per spike", but his Figs S9–S10 label the trait "Number of grains" (beside "Tiller length", "Flag length", "Plant biomass", i.e. the raw column names), so he may have used `number_of_grains` too. Decided with the user 2026-09-29 |

**¶3 Post-harvest**
- 5 spikes per plant, bagged, oven-dried 35 °C for 30 d; threshed (manually 2019; Haldrup LT-21 later); seed blower.
  **SRC** old drafts / thesis → **NOAM** (threshing in 2021; equipment).

## 3. Grain composition by NIR

**¶1 Measurement**
- NIR DA 7250 (Perten), whole grain, non-destructive; protein, starch, β-glucan, fiber. **SRC** NIR columns in `all years barley.csv`.
- 3,553 wild plants carry NIR values (2019: 1,191; 2020: 1,173; 2021: 1,189; before removing block 4 of 2019), **SRC** data; the old drafts' "~3,900" is all plants. Values are **% of grain weight on an as-is basis** (not corrected to dry matter; confirmed by the user 2026-09-29); these are the units of the BLUPs and of the allelic effects in Results Ch. 2.
- Instrument courtesy of Equinom → Acknowledgments.

**¶2 Calibration and validation**
- 30 accessions (one per site): starch Megazyme K-TSTA, β-glucan Megazyme mixed-linkage kit (**CHECK** catalog code),
  milling (IKA A11, 50 mesh), 3 × 200 mg subsamples; protein by Kjeldahl (Tel Hai lab; **CITE** Beljkaš et al. 2010 is weak).
- Protein and starch: cultivated package adjusted with the lab values; β-glucan calibrated from the kit values; **fiber:
  cultivated package, no wild-barley reference** (state as a fact).
- Validation, Pearson r: kit vs calibration scans, kit vs 2019 scans, 2019 vs calibration scans. **SUPP** calibration
  table. **NOAM** values come from the thesis (Table 1); β-glucan r = 1.000 looks circular.
- **DROP** the 2021 starch +8.85 adjustment (user).

## 4. Phenotypic data analysis

**¶1 Data preparation** — **SRC** `00_THIN_…/01_…no_GxE…R`, `02_…GxE…R` (headers)
- Removed: cultivated checks, site 04, **block 4 of 2019** (low inter-block concordance). **SUPP** block correlations —
  **TO GENERATE (user) before writing this paragraph and before the supplement:** the 2019 inter-block correlations of the
  four NIR traits. The old drafts' r (β-glucan 0.328–0.391 vs 0.670–0.682) have no source in the project; do not quote them.
- Values centred within season; per trait: missing values excluded, ± 3 SD outliers removed.
- **CHECK** R version for 00_THIN (old drafts: 4.4.3).

**¶2 BLUPs**
- lme4 1.1.35.3 in R 4.1.2 (corrected 2026-09-29, user: the 00_THIN scripts ran on fidel, not Posit Cloud; their headers were corrected the same day) (**CITE** Bates et al. 2015): trait ~ (1 | genotype) + (1 | season:block). BLUPs → GWAS, correlations, site analyses.

**¶3 Heritability and G×E**
- Same model + (1 | genotype:season); VG, VE (season:block), VG×E, VR; H² = VG / total (**CITE** Holland et al. 2003).
  Gives Table 1 and Fig. 1a (8 traits). Tillers H² < 0.001.
- Reaction norms (Fig. 1b): 4 nutritional traits; 4 highlighted accessions per trait, picked in turn: most stable, strongest
  increase, strongest decrease (2019 → 2021), strongest crossover. **SRC** `02_…GxE…R` §9.

**¶4 Trait correlations**
- Pearson among BLUPs: 4 × 4 nutritional (Fig. 2b); nutritional × 4 morphological (Fig. 2c). Stars = raw *P* (no need to stress it).

**¶5 Environment** — **SRC** `03a`, `03b` headers; `merged_Evgeny_env_fitness_environmental_30_sites.csv`
- 29 sites; site-mean BLUPs; 10 site variables (**CITE** *MBE* for the sources: Israeli Meteorological Service, ISRIC soil, field samples).
- Collinearity |r| > 0.8: elevation (vs March temperature, r = −0.951) and silt (vs sand, r = −0.934) dropped → 8 variables:
  clay, sand, precipitation, organic carbon, pH, EC, total N, March temperature.
- Two-sided Pearson (Fig. 2d); Spearman as a sensitivity check. Ecological regions for Fig. 2a from the sampling design.
- **DROP** the 4 June draft's mixed model with backward elimination.

## 5. Genome-wide association

**¶1 Model** — **SRC** `01_…/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_analysis_parameters.tsv`
- EMMAX (**CITE** Kang et al. 2010; **CHECK** version); BLUPs; all 7,110,996 SNPs.
- **One LD-pruned set, used three times** (PCA, kinship, Bonferroni count): PLINK 1.9 `--indep-pairwise 1000kb 1 0.2` → 111,017 SNPs;
  a 1 Mb physical window because SNP density varies widely along the genome (**SRC** `01_…/README.md`, `scripts/01_setup/02_ld_prune.sh`).
  **MUST** this is the only place the pruned set and its count appear (moved out of the Results, `MM_HANDOVER` 2.1a).
- PCA → first 3 PCs as fixed covariates (**CITE** *Mol Ecol*
  for 3 PCs; **SUPP** PC variance table); **aIBS kinship** (emmax-kin) as random effect: state it (*Mol Ecol* used Balding–Nichols).

**¶2 Significance** — **SRC** `MM_HANDOVER` item 2.1a
- Bonferroni α = 0.10 over the 111,017 pruned SNPs → −log10(*P*) = 6.0454; reason (independent tests; permissive α for
  n = 290), as in the mini paper. λGC per trait; Manhattan and QQ plots (Fig. 3).
- **MUST** effect direction: β = effect of the minor allele, in % of grain (Results Ch. 2 ¶3).

**¶3 LD decay** — **SRC** `02_USED_LD_decay_V2_wholegenome/results/ld_decay_summary_for_results_chapter.txt`
- PLINK r² for SNP pairs within 2.5 Mb; 1-kb bins; pair-weighted LOESS (span 0.10); r² = 0.2 at ~188 kb. **SUPP**
  genome-wide + per-chromosome plots. Cited only from the M&M.

## 6. Loci and candidate genes

**¶1 Locus definition** — **SRC** `01_…/02_loci_FINAL/tables/locus_definition_params.tsv`, `02_loci_FINAL/README.md`
- PLINK `--clump`: lead *P* < 9.008 × 10⁻⁷; members by LD alone (p2 = 1, r² ≥ 0.5); ± 2 Mb; locus cut at the first gap
  > 50 kb; repeated until every significant SNP lies in a locus → 36 loci. **SUPP** locus table (Online Resource 1).

**¶2 Gene search** — **SRC** `03_…/README.md`, `results/tables/results_chapter_numbers.txt`
- Interval = the locus span, no flank; bedtools intersect 2.30.0 (**CITE** Quinlan and Hall 2010) with the MorexV3 gene
  annotation (Ensembl Plants release 62), 1H–7H → 55 genes. **SUPP** 55-gene table.
- Only 3 of 55 carry a functional description → why the annotation chain (§8). **SRC** `MM_HANDOVER` item 5.4.
- **CHECK** high-confidence genes only? (low-confidence genes not searched).

## 7. Haplotype analysis

**¶1 Input** — **SRC** `04_…/README.md`, `00_config/config.yaml`, `02_imputation/PROVENANCE.md`
- Per-gene windows: gene ± 1 kb, raw genotypes (bcftools), 290 accessions.
- Imputed genotypes used **only for LD**: LEA sNMF, K = 3, 10 runs, lowest cross-entropy, mode (**CITE** Frichot and François
  2015; **CHECK** LEA version). K = 3 → **CITE** *Mol Ecol* (same collection; sNMF/LEA, K = 3–5 compared by PCA and structure
  plots, three genetic clusters, its Fig. 1b; the number of clusters does not depend on the reference genome); Hübner et al. 2009/2013
  optional as a second citation for the ecotypes (decided with the user 2026-09-29). Imputed file made outside this project:
  `/mnt/data/shahar/gwas_barley/morexV3_analysis/USED_imputation/` (scripts copied into `04_…/02_imputation/`).
- **SUPP** (optional) imputation summary: missing-genotype rate before imputation, cross-entropy of the 10 sNMF runs at
  K = 3. **CHECK** recoverable from the `.snmfProject` in `USED_imputation/`.

**¶2 crosshap** (**CITE** Marsh et al. 2023; **CHECK** version)
- LD matrix from imputed data (PLINK r²), haplotypes from raw genotypes; minHap = 9; unassigned accessions excluded.
- **ε = 0.6 and MGmin = 2, fixed for all genes: report the reason for the choice** (the 2026-09-23 decision dropped only the
  explanation of what ε is). Reason, **SRC** `04_…/README.md` § "The two fixed parameters" +
  `04_runs/loci_LDspan_eps06_V4/Diagnostics/epsilon_coverage_comparison.tsv`: fixed ε values compared on genotype-only criteria
  (no *P* values); 0.6 assigned the most accessions to haplotype groups (median 70% vs 61–65% at 0.8–1.0) at the cost of
  slightly fewer testable genes (30 vs 33–34); an unassigned accession contributes nothing to the test.
  - **DECIDE** those numbers come from an **earlier 64-gene candidate set** (V2 locus definition), not the current 55 genes:
    (a) state the criterion without numbers, or (b) re-run the ε comparison on the current set
    (`01_scripts/diagnostics/eps_coverage_scan.R`; analysis re-run, user approval) and quote it (**SUPP** possible).
  - **CHECK** MGmin = 2: the README's only reason ("the v1 publication genes used it") refers to the retired pipeline and
    cannot be written. Find a positive reason (e.g. the smallest marker group, which lets SNP-sparse genes be tested;
    check the crosshap documentation) or ask the user.
- **What ε and MGmin are: do not explain, cite crosshap** (Marsh et al. 2023, where the parameters are defined). **DROP** the
  explanation of ε as a radius in the r² profile space, not an r² cutoff (user, 2026-09-23; confirmed 2026-09-29).

**¶3 Tests**
- One Kruskal–Wallis test per gene; BH on raw *P* within trait; *q* ≤ 0.05.
- **MUST** BH valid under positive dependence (**CITE** Benjamini and Yekutieli 2001); untestable genes left out of the BH
  denominator but reported with the reason (**SUPP** 55-gene table).
- Effect sizes: η² = (H − k + 1)/(n − k); highest − lowest group difference in SD. Post hoc: Wilcoxon vs the largest group,
  Holm (Fig. 4 brackets).

## 8. Annotation and choice of candidate genes

**¶1 Annotation** — **SRC** `05_USED_gene_annotation_analysis/README.md` + sub-READMEs
- PGSB MorexV3 high-confidence proteome, longest isoform; DIAMOND 2.0.14 BLASTP vs UniProt Swiss-Prot (2026_03; corrected 2026-09-30 from 2026_01, per intermediates/versions.txt) → InterProScan 5
  (Pfam, InterPro, GO) → NCBI nr BLASTP (Viridiplantae) for the unresolved; fixed priority for the final call.
  **SUPP** Online Resource 2 (23 genes). **CHECK** hit thresholds; InterProScan version.

**¶2 Choice** — **SRC** `05_…/08_USED_annotation_master/TRAIT_CANDIDACY.md`
- Significance first (23 genes), then protein function and literature: a documented or plausible role in the biosynthesis,
  transport or regulation of the grain component. Rank by *q* or effect size not used.

**¶3 Canonical β-glucan genes** — **SRC** `05_…/09_…/README.md`, `results/tables/canonical_gene_distances_current_loci.tsv`
- 10 synthase and endohydrolase genes located on MorexV3 from their reference protein sequences (published positions are not
  MorexV3); distance to β-glucan leads. **MUST** (Discussion, TODO 410). **CHECK** the mapping method in the step-09 scripts
  01–05 (BLASTP reciprocal best hit + tBLASTn agreement).
- **SUPP** one table, mapping + distances: gene, UniProt accession, MorexV3 gene ID, position, % identity / coverage,
  reciprocal best hit, tBLASTn agreement (mapping columns of `canonical_betaglucan_gene_distances.tsv`, whose distance
  columns are stale) + distance to the nearest β-glucan lead (`canonical_gene_distances_current_loci.tsv`).

## 9. The shared 7H signal (Results Ch. 4)

**¶1 Gene and annotation** — **SRC** `03_01_7H_branch_…/archive_v1_MGmin2_eps0.6/results/tables/`
- Region examined directly (both loci too narrow to hold a gene; the gene lies 554 bp beyond the starch locus).
- InterProScan, SignalP, Phobius (**CHECK** versions); association profile relative to the annotated TSS; LD among the signal SNPs.

**¶2 Haplotypes** — **SRC** `03_01_…/README.md` §2 + §6; `MM_HANDOVER` item 3b
- **MUST** both settings: at MGmin 2 / ε 0.6 the signal SNPs were not grouped (*P* = 0.16 fiber, 0.36 starch); then MGmin = 3
  (smallest value grouping them) and ε = 0.9 (highest assignment; tied with 1.0, smaller taken). GWAS-informed; no
  haplotype-test *P* used.
- Raw *P* (one pre-specified gene); Kruskal–Wallis, η², SD; Wilcoxon A vs B. **SUPP** parameter grid.

## 10. Elite cultivars

**¶1 Source and lines** — **SRC** `07_…/README.md`, `results/tables/results_chapter_numbers.txt`
- IPK DivBrowse barley pangenome v2 (MorexV3, unimputed; **CITE** the pangenome and DivBrowse papers, verify).
- 136 spring elite lines → 5 spring malting cultivars released 2012–2018; names verified in EBI BioSamples.
- **Call-rate rule, stated precisely:** ≥ 85% of the shared sites **pooled over the three genes** (a genotype-quality
  criterion; per line 86–99%). Per gene it is lower: at *PHT4;3* four lines are at 83% (10 of 12 sites).
  **SRC** `07_…/results/tables/Table_elite_lines.tsv` (call_rate_* columns), `scripts/02_screen_elite_lines.R`.
  **SUPP** elite-line table with per-gene call rates (incl. the GDSL gene).

**¶2 Comparison**
- Same windows (gene ± 1 kb); only sites with identical CHROM:POS:REF:ALT in both call sets (allele concordance checked)
  (`MM_HANDOVER` 2.1c). **SUPP** site accounting (Online Resource 3).
- Wild group row = per-SNP majority, ties → REF.
- **MUST** crosshap was not re-run and the cultivars were never assigned to a group (rules out circularity).
- **MUST** Morex, the reference, is itself an elite cultivar (the Discussion relies on it).

**¶3 The 7H GDSL gene** (Fig. 5c) — same five cultivars and procedure, **not re-screened** for this gene
- **MUST** state the call rates here: over its 22 shared sites only Avalon and KWS Irina reach 85% (95.5%); Odyssey 82%,
  LG Diablo 77%, Laureate 55%. Fig. 5c draws the missing calls; the Results sentence ("…at every called site") stays as is.
  **SRC** `03_01_…/results/tables/Table_elite_line_screen.tsv`; `03_01_…/README.md` § 3 (call-rate caveat). Decided with the user 2026-09-29.

## 11. Geographic origin of minor-allele carriers — analysis finished 2026-09-30

**SRC** `11_USED_check_low_maf_geographic_distribution_minor_allele/README.md` § "Methods (facts for M&M §11)"; verify in
`scripts/00_config.R` and scripts 02–05. In the M&M and the Discussion only, **not in the Results** (user); the Discussion
cites it as an Online Resource in §3 ¶2 (TODO 419) and in "Wild alleles" (TODO 420); its 7H sentence (§2 ¶2) has no citation (user).
- Units: the lead SNP of each of the 36 loci (minor- vs major-allele carriers; missing calls excluded) and the haplotype
  groups of the four presented genes (as in Fig. 4 / Table 2 and Fig. 5; unassigned accessions excluded).
- Regions and sites: the Fig. 2a assignment, 29 sites. **Wording (user):** the three ecotypes and the transitional zones
  between them, never "six regions".
- Site-permutation tests: accessions kept in their site, region labels permuted among the 29 sites (10,000 permutations);
  Desert share (two-sided) and distribution over regions; BH within trait (leads), within gene (groups), across the 4 genes;
  conservative where carriers come from one site.
- Matched genome-wide background: 500 SNPs per lead (MAF ± 0.005, ± 10 called accessions, > 2 Mb from any lead) → which
  region the carriers mainly come from, observed vs expected.
- Within-site comparison (descriptive, no re-test): carriers vs the non-carrier accessions of their own site.
- **DROP** (user, 2026-09-30): the carrier-sharing analysis (T09, T10) and the kinship-by-region comparison (T11); neither is in the text.
- Software: R 4.1.2 and packages, PLINK 1.9 (listed in the README).
- **SUPP**: T02 (+ T01 as a second sheet), T04 (+ T03), Figs S1–S2 (captions ready in `results/CAPTIONS.md`); supporting T05, T06, T08.

## 12. Software, availability, statements

- Software: R (**CHECK** 4.1.2 for steps 01–07, 4.4.3 for 00_THIN?) and main packages; figures drawn in R at final size
  (TAG: name the graphics software).
- End of M&M: accession **ENA PRJEB79623**; MorexV3 VCF (Zenodo or on request → **DECIDE**); code (GitHub hubner-lab / Zenodo → **DECIDE**).
- **LLM-use statement** (TAG requires it): **last — the user handles it at the end.**

---

## Coverage check — Results/Discussion items the M&M must back

| item in the manuscript | M&M |
|---|---|
| BLUPs, H², G×E, reaction norms (Ch. 1, Table 1, Fig. 1) | §4 ¶2–3 |
| site/region display, trait and environment correlations (Fig. 2) | §4 ¶4–5 |
| GWAS, λGC, threshold, minor-allele β (Ch. 2, Fig. 3) | §5 |
| 36 loci, spans, co-localization (Ch. 2) | §6 ¶1 |
| 55 genes, LD span (Ch. 2) | §6 ¶2 |
| tested / not tested, *q*, η², SD, brackets (Ch. 3, Table 2, Fig. 4) | §7 |
| annotation, identity/coverage (Table 2), choice of 3 genes | §8 ¶1–2 |
| elite barcodes, consensus, shared sites (Ch. 3, Fig. 4d–f) | §10 |
| 7H gene, promoter, MGmin 3 / ε 0.9, raw *P* (Ch. 4, Fig. 5) | §9 |
| canonical β-glucan gene distances (Discussion) | §8 ¶3 |
| carrier origin (Discussion) | §11 |
| grouping setting can miss a gene (Discussion) | §9 ¶2 |


---

## Appendix A — Evgeny Potapenko's reply on the MorexV3 calling (received by the user, 2026)

Context: the user asked Evgeny, who ran the MorexV3 re-call, to update the thesis paragraph (MorexV2: BWA-MEM2, Picard
2.8.1, GATK 4.1.2.0, hard filters, heterozygotes → missing, INDELs removed, VCFtools MAF filter) and to confirm MAF 0.05 and
missingness ~30%. His reply, verbatim (his citation for MorexV3, Monat et al. 2019, is the V2 assembly; see §1 ¶3):

> Cleaned reads were aligned to the Morex V3 reference genome (Monat et al., 2019) using BWA-MEM2 v2.2.1 (Vasimuddin et al., 2019) with default parameters. Duplicate reads were identified and marked using Picard MarkDuplicates v3.4.0 (Picard Toolkit, Broad Institute <CITE PICARD TOOL>).
>
> Variant discovery was performed following the GATK4 best-practices workflow using GATK v4.6.2.0 (Poplin et al., 2017). For each sample, variants were first called in GVCF mode using HaplotypeCaller with a diploid ploidy setting. To accommodate the large genome size variant calling was conducted in non-overlapping genomic intervals. Resulting GVCFs were subsequently processed using ReblockGVCF to improve compression and standardize block structure prior to joint genotyping.
>
> Joint genotyping was performed per genomic interval using GenomicsDBImport, followed by GenotypeGVCFs to produce cohort-level VCF files. To reduce memory usage and spurious multi-allelic calls in low-coverage regions, genotyping was restricted to a maximum of two alternate alleles.
>
> Post-genotyping variant filtering was conducted using bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic SNPs were retained, and INDELs were removed. Variants were filtered using hard thresholds on quality metrics (QUAL, QD, MQ, FS, SOR, MQRankSum, and site depth) [Write here parameters or link to the supplementary table]. Heterozygous genotypes were set to missing data, and genotypes with read depth below 3 were excluded. Variants failing any filter were removed from downstream analyses.
>
> Filtering command: `bcftools filter -e 'QD<5 || MQ<45 || FS>60 || INFO/DP<900 || INFO/DP>1800 || QUAL<140 || SOR>3 || MQRankSum<-2.5'`
