# 11_USED_check_low_maf_geographic_distribution_minor_allele

**Where do the accessions come from that carry the GWAS lead-SNP alleles and the haplotype groups of
the four candidate genes: the Desert, or across the collection?** Built 2026-09-25; reorganized into a
finished, supplement-ready analysis 2026-09-30 (see *History*). 290 *H. spontaneum* accessions, 29 sites,
6 ecological regions, MorexV3.

It backs the Discussion of
`10_USED_Paper_writing/new_publishing_paper/build/UPDATED_Results_Discussion_Conclusions.md` (§2 ¶2, 7H sentence; §3 ¶2;
"Wild alleles"), cited there as an Online Resource (TODOs 419, 420; manuscript review items 11–11c, 2026-09-30), and
M&M §11 of `build/MM_BLUEPRINT.md`. It is not in the Results (user). Its tables and figures are Online Resource candidates; the user
decides what enters the supplement.

**Start here**
- [`results/CAPTIONS.md`](results/CAPTIONS.md): ready-to-paste TAG-style captions for both figures and every table, with column definitions
- [`results/results_numbers.txt`](results/results_numbers.txt): every number of this step, as prose (generated)
- [`results/figures/Fig_carrier_origin.png`](results/figures/Fig_carrier_origin.png) and [`Fig_desert_enrichment.png`](results/figures/Fig_desert_enrichment.png)
- *Verdict* and *Discussion claims → evidence* below

---

## The question

At every genome-wide significant locus the minor allele of the lead SNP moves the trait the same way: it
raises fiber (18 of 18 loci) and β-glucan (11 of 11) and lowers starch (5 of 5); 25 of the 36 leads have
MAF < 0.10. Two readings:
- **Biology.** Desert accessions have low starch and high fiber and β-glucan (Ch. 1). Locally adapted
  minor alleles would push grain composition toward the desert end of that axis.
- **Artefact.** PC + kinship correction is weakest for rare alleles concentrated in a few related
  populations, and λGC does not detect it. If the minor alleles sit mostly in desert accessions, which
  already have that phenotype, the uniform direction could be residual population structure.

The same question is asked of the haplotype groups of the four presented genes (GPAT6, GH17, PHT4;3,
GDSL): do the groups that separate the phenotype come from different regions or sites?

## Constraints (user, 2026-09-25)

- Read-only on the project; every output is inside this folder.
- **No GWAS re-run, no region covariate** (not approved). Every quantity is a count, a permutation test
  of origin, or a mean of the phenotypes the GWAS already used.
- The minor allele is never called "derived" (no outgroup).
- The call set has no heterozygous genotypes (set to missing upstream; asserted in scripts 02 and 04),
  so a carrier is a homozygote.

---

## Folder layout and pipeline

```
scripts/        00_config.R  00_figure_helpers.R  01 … 10  run_all.sh
intermediates/  lead genotypes, long tables, background SNP genotypes (regenerated)
results/        tables/T01–T11 · figures/Fig_*.{png,tif} · CAPTIONS.md · results_numbers.txt
logs/           one log per script
```

```bash
bash scripts/run_all.sh          # ~4 min, rebuilds everything in order
```
Seeds are fixed (20260925); a clean rebuild reproduced every table exactly (checked 2026-09-30).

| script | does | writes |
|---|---|---|
| `00_config.R` | paths, constants, loaders, **the site-permutation tests**, table writer | — |
| `00_figure_helpers.R` | shared rows/labels, accession axis, theme, PNG + TIFF export | — |
| `01_extract_lead_genotypes.sh` | the 36 lead SNPs × 290 accessions, PLINK 1.9 `--recode A --keep-allele-order` (counted allele = A1 = minor) | `intermediates/lead_genotypes.{raw,frq}` |
| `02_lead_carrier_origin.R` | minor / major / missing per lead; region and site of the carriers; **Desert and six-region tests** | **T01**, **T02** |
| `03_haplotype_group_origin.R` | haplotype groups of the 4 genes; region and site of the members; **group- and gene-level tests** | **T03**, **T04** |
| `04_matched_background.R` | 500 matched genome-wide SNPs per lead; where their carriers sit | T05, T06, T07 |
| `05_lead_within_site.R` | carriers vs their own non-carrier site-mates; dependence on extreme accessions | T08 |
| `06_lead_carrier_sharing.R` | shared carriers between loci; minor alleles per accession | T09, T10 |
| `07_kinship_by_region.R` | aIBS kinship within sites, between sites of a region, between regions | T11 |
| `08_figure_carrier_origin.R` | Fig. S1 (barcode: a leads, b haplotype groups) | `Fig_carrier_origin.{png,tif}` |
| `09_figure_desert_enrichment.R` | Fig. S2 (Desert share + q of both tests) | `Fig_desert_enrichment.{png,tif}` |
| `10_results_numbers.R` | every number, as prose | `results_numbers.txt` |

### Inputs (read-only)

| what | file |
|---|---|
| Lead SNPs, minor allele (`lead_A1`), MAF, beta (minor allele) | `01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv` |
| Genotypes (the GWAS set) and allele frequencies | `01_USED_GWAS_V2_pipeline/intermediates/morexV3_290.{bed,bim,fam}`, `morexV3_290_freq.frq` |
| Kinship (the GWAS's, EMMAX aIBS; rows = .fam order) | `01_USED_GWAS_V2_pipeline/intermediates/morexV3_kinship.aIBS.kinf` |
| Region, site, location per accession (the Fig. 2a assignment) | `00_THIN_Generate_Plots_For_Publication/outputs/subsection_1/tables/A4_site_boxplot_data.csv` |
| Phenotypes: the BLUPs EMMAX used | `/mnt/data/shahar/gwas_barley/data/inputs/<trait>_corrected_V3.pheno` |
| Haplotype groups, GPAT6 / GH17 / PHT4;3 (Fig. 4, Table 2) | `04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Cache/<trait>/<gene>/MGmin_2/HapObject.rds` (Indfile; MGmin 2, ε 0.6) |
| Haplotype groups, GDSL (Fig. 5) | `03_01_7H_branch_Starch_Fiber_shared_signal_explore/results/tables/mgmin3_haplotype_assignment.tsv` (MGmin 3, ε 0.9) |

Checks built into the scripts: all 290 accessions have a region (genotype ID `HS0103` = region-table
`01_03`); all 36 leads are in the bim with the loci table's A1/A2; the crosshap `Pheno` equals the GWAS
BLUPs; the GDSL grouping is identical for fiber and starch; no heterozygous calls.

---

## Methods (facts for M&M §11)

- **Units.** (i) The **lead SNP** of each of the 36 loci: accessions split into minor-allele carriers,
  major-allele carriers and missing calls (excluded from the tests). (ii) The **haplotype groups** of the
  four presented genes, exactly as defined for Fig. 4/Table 2 and Fig. 5 (crosshap not re-run);
  unassigned accessions excluded from the tests.
- **Regions and sites.** The six ecological regions of Fig. 2a (North 4 sites / 40 accessions; Coast,
  Desert, HZ1 North–Coast, HZ2 North–Desert, HZ3 Coast–Desert 5 sites / 50 accessions each), 29 sites.
- **Site-permutation tests.** Accessions stay in their site; the region labels are permuted among the 29
  sites (preserving 4 North and 5 per other region), 10,000 permutations, one shared permutation set.
  This keeps the clustering of accessions within sites, which an accession-level test would ignore; it
  is conservative for alleles or groups found in a single site, where region and site cannot be
  separated.
  - **Desert test**: statistic = desert share of the focal set − desert share of the other set (minor vs
    major carriers; a group vs the other groups of its gene); two-sided *P* = 2 × the smaller tail.
  - **Regions test**: chi-square statistic of set × six regions; upper-tail *P*. For the haplotype groups
    also at gene level (all groups of a gene × six regions).
  - Benjamini–Hochberg: within trait (leads); within gene (group tests); across the 4 genes (gene-level
    regions test).
  - Effect size: accession-level odds ratio of desert origin (Fisher conditional MLE), reported only as an
    effect size.
- **Matched genome-wide background** (per lead): 500 SNPs with |MAF − lead MAF| ≤ 0.005 and called
  chromosomes within ± 20 (± 10 accessions) of the lead, from the 7,110,996-SNP set, excluding ± 2 Mb
  around every lead (17,956 distinct SNPs). Per lead: position of its desert share and site count in its
  matched distribution; per trait: number of leads whose top carrier region is each region, observed vs
  the sum of the per-lead background probabilities (Poisson-binomial, 10⁵ simulations).
- **Within-site comparison** (descriptive; no re-test): observed difference between minor- and major-allele
  carriers (BLUPs, trait SD) vs the difference expected if each carrier had the mean of the major-allele
  carriers of its own site (sites without one fall back to the region); share of carriers exceeding their
  site-mates in the direction of the effect; the observed difference without the trait's two most extreme
  accessions.
- **Carrier sharing**: per pair of loci of one trait, shared minor-allele carriers vs chance
  (hypergeometric). **Kinship**: mean aIBS for same-site pairs, different-site pairs within a region, and
  between regions.
- **Software**: R 4.1.2 (dplyr 1.1.4, tidyr 1.3.1, readr 2.1.5, data.table 1.17.8; figures ggplot2 4.0.1,
  patchwork 1.3.2, ragg 1.3.0; TIFF via Pillow 12.2.0), PLINK 1.9 (v1.90b6.4).

---

## Results

Full numbers: `results/results_numbers.txt`.

**1. Lead SNPs: fiber and β-glucan minor alleles come mainly from the Desert and the coast–desert zone**
(T01, T02, T06; Fig. S1a, S2a)

| | fiber (18) | β-glucan (11) | starch (5) | protein (2) |
|---|---|---|---|---|
| minor-allele instances from Desert / Desert + HZ3 (panel 17% / 34%) | 59% / 74% | 50% / 78% | 19% / 36% | 0% / 11% |
| leads whose carriers come mainly from the Desert (expected from matched background) | **14** (3.9) | **8** (2.5) | 2 (1.2) | 0 |
| leads whose carriers come mainly from the North (expected) | **0** (2.8) | 0 (1.6) | 0 (0.6) | 0 |
| Desert test q ≤ 0.05, per lead (two-sided) | 8 enriched | 4 enriched | 1 enriched, 1 depleted | 1 depleted |
| six-region test q ≤ 0.05, per lead | 2 | 3 | 1 | 1 |

- Fiber + β-glucan: Desert is the main origin at **22 of 29** leads vs 6.5 expected from matched SNPs
  (*P* = 1e−5), North at none vs 4.5 (*P* = 0.008). Optimistic, because the loci share carriers (3).
- Per lead, significance is limited: many leads have 12–17 carriers from 3–6 sites, and the
  site-permutation test cannot separate region from site with that few sites. At 28 of the 29 fiber and
  β-glucan leads the minor-allele carriers are more often from the Desert than the major-allele carriers,
  and all 28 draw them from ≥ 2 desert sites.
- Rare alleles in this panel are regionally private as a rule (a matched SNP typically has over half its
  carriers in one region; Desert is the main origin of 22% of rare background SNPs). The leads are not
  more clustered than their matched SNPs; what is unusual is **which** region (T05, T07).
- **7H GDSL locus**: fiber 7H:573,606,306 — 41 carriers from 11 sites in 5 regions (0 North), 17% desert
  (major carriers 15%), q Desert 0.81, q regions 0.11; starch 7H:573,606,460 — 35 carriers, 10 sites,
  20% desert, q 0.82 / 0.41. No regional concentration.
- **Starch 6H:525,776,080** (common, MAF 0.42): carriers absent from the Desert (q = 0.024, depleted) and
  the only lead whose raw carrier difference is opposite to its GWAS beta (T08): the association exists
  only after the PC + kinship correction. **Protein 3H:106,623,911**: no desert carriers (q = 0.013, depleted).

**2. Haplotype groups: the high-fiber groups of GPAT6 and GH17 are desert groups; GDSL groups are not
structured** (T03, T04; Fig. S1b, S2b)

| gene | group (mean BLUP) | members N/C/D/HZ1/HZ2/HZ3, sites | desert % vs other groups | q Desert | gene-level regions q |
|---|---|---|---|---|---|
| GPAT6 (fiber) | **C** (+0.28, highest) | 0/0/12/0/0/6, 5 sites | 67 vs 7 | **0.034** | 0.097 |
| GH17 (fiber) | **C** (+0.20) | 0/0/15/1/0/0, 4 sites | 94 vs 11 | **0.043** | 0.097 |
| | E (+0.23, highest) | 1/2/5/1/0/1, 5 sites | 50 vs 21 | 0.34 | |
| | A, B, D (lower fiber) | | 0–8 vs 27–31 | 0.045 (depleted) | |
| PHT4;3 (starch) | **A** (+1.63, high starch) | 26/33/7/29/27/15, 25 sites | 5 vs 40 | **0.0024** (depleted) | 0.32 |
| | B, C, D (low starch) | | 34–60 vs 12–14 | 0.11–0.29 | |
| GDSL (fiber, starch) | B (high fiber, low starch) | 0/5/7/7/14/1, 10 sites | 21 vs 15 | 0.67 | 0.32 |

At GPAT6, GH17 and PHT4;3 the groups that carry the high-fiber or low-starch phenotype come mainly from
the Desert and HZ3, and the high-starch PHT4;3 group comes from everywhere but the Desert. The GDSL
haplotype B is spread over five regions with the panel's desert share, as its lead SNPs are.

**3. The loci share carriers** (T09, T10)
- Pairs of loci on different chromosomes sharing more carriers than chance (*P* < 0.001): fiber 41 of
  124, β-glucan 20 of 48, starch 3 of 9.
- HS1518 and HS1509 (Arad, Desert) carry the minor allele at 18 and 17 of the 18 fiber loci and 10 and 11
  of the 11 β-glucan loci. They are the #1 and #2 fiber accessions of the panel; HS1509 is #1 for β-glucan
  and the lowest in starch. Genome-wide they rank only 88th and 19th of 290 in rare-allele load.
- 31 of the 50 desert accessions carry fiber minor alleles at ≥ 3 loci, against 9 of the other 240.

**4. Carriers differ from their own site-mates** (T08)
- The carriers' sites predict a median 17% (β-glucan), 33% (fiber), 29% (starch) and 15% (protein) of the
  raw carrier difference; 81–94% of carriers exceed their site-mates in the direction of the effect;
  79–98% of the difference remains without the two most extreme accessions.
- Caveats: the loci were selected by a test that already conditioned on kinship (ascertainment); and at
  **9 loci fewer than 60% of carriers have any non-carrier site-mate** (whole sites fixed for the minor
  allele), including the **GH17** (7 of 16) and **PHT4;3** (12 of 31) loci, where allele and population
  cannot be separated.

**5. Relatedness** (T11). Mean aIBS between accessions of different sites of the same region: Coast
0.694, **Desert 0.688**, HZ1 0.683, HZ3 0.677, North 0.664, HZ2 0.660 (all pairs 0.671). Within sites,
Desert is the highest (0.782).

---

## Verdict

**Mixed. The uniform direction of the minor-allele effects is not an independent line of support for
the carbon-allocation reading.**
1. **The artefact premise holds.** Fiber and β-glucan minor alleles, and the high-fiber/low-starch
   haplotype groups of GPAT6, GH17 and PHT4;3, come mainly from desert and coast–desert accessions, the
   populations that already have that grain composition; none is concentrated in the North.
2. **The consistent directions are largely one observation.** The loci share carriers across chromosomes,
   and a small group of Arad, Masada, Mount Ramon, Maale Akrabim and Revivim accessions, including the
   panel's most extreme ones, carries most of the minor alleles. Rare alleles here are regionally
   private, the phenotypic extremes sit at the desert end, and rare variants are detected through extreme
   carriers; together these produce a uniform direction whether or not each effect is causal.
3. **A pure residual-structure artefact is not supported either.** Within sites, carriers are more extreme
   than their own site-mates; origin predicts a minority of the carrier difference; it survives removing
   the two most extreme accessions. This is weakened by ascertainment and unavailable at the 9 loci whose
   carriers sit in sites fixed for the minor allele (including GH17 and PHT4;3).
4. **The 7H GDSL signal is the clear exception**: its lead-SNP carriers and its haplotype B come from five
   regions at the panel's desert share, so this shared fiber–starch signal cannot be attributed to origin.

Separating "locally adapted desert alleles" from "desert-lineage tagging" would need a within-region
re-test or an independent population; neither is approved.

## Discussion claims → evidence

Current text, after manuscript review items 11, 11b and 11c (2026-09-30).

| Discussion statement (build/UPDATED_Results_Discussion_Conclusions.md) | evidence here | status |
|---|---|---|
| §2 ¶2: the 7H minor-allele carriers, "and of the haplotype carrying both effects, came from sites across the environmental gradient of the collection, so this shared signal cannot be explained as a geographic artifact of the accessions' origin" (no Online Resource citation, user) | T02 fiber_L17, starch_L05; T04 GDSL B; Fig. S1, S2 | supported: 10–11 sites, desert at the panel share (q 0.81 / 0.82; haplotype B q 0.67). Wording is the user's; none of the carriers is from the four North sites (0/41, 0/35), recorded in the hidden note |
| §3 ¶2: "at most fiber and β-glucan loci, they originated mainly from desert and coast–desert sites" (Online Resource, TODO 419) | T06 (Desert main origin at 22/29 vs 6.5 expected); instances 74–78% Desert + HZ3 | supported as a collective pattern; per lead the two-sided Desert test is q ≤ 0.05 at 12 of 29 |
| §3 ¶2: "The carriers also came from several desert sites rather than a single one" | T02 `n_desert_sites_minor` | supported: 28 of 29 fiber/β-glucan leads lean to the Desert, all from ≥ 2 desert sites |
| §3 ¶2: "at most loci they differed from non-carriers of their own site in the direction of the allelic effect" | T08 | supported (81–94% of carriers); 9 loci largely lack non-carrier site-mates |
| §3 ¶2: "the haplotypes associated with high fiber at *GPAT6* and *GH17* came mainly from several desert populations, and the high-starch haplotype of *PHT4;3* was nearly absent from them" | T04; Fig. S1b, S2b | supported: GPAT6 C q 0.034 (3 desert sites), GH17 C q 0.043 (C + E 20 of 26 desert, 4 desert sites; E alone q 0.34), PHT4;3 A 5% desert, q 0.0024 |
| §3 ¶2: "The 7H locus, whose carriers and haplotypes were not concentrated in the desert, therefore provides the clearest case" | T02, T04 | supported |
| "Wild alleles": "At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site" (Online Resource, TODO 420) | T04; `intermediates/haplotype_group_long.tsv` | supported: Masada 6, Mount Ramon 3, Arad 3; Revivim 5, Orim 1 |
| Dropped by the user (2026-09-30): "the desert populations are closely related to one another" | T11 | only partly supported (Coast higher than Desert); removed with TODO 417 |
| Not reported by the user's decision: the same few desert accessions carry the minor alleles at many loci | T09, T10 | kept here as background only |


---

## Outputs

| file | one row per | Online Resource? |
|---|---|---|
| `T01_lead_carriers_by_region.tsv` | lead × region (minor / major / missing) | sheet 2 of the T02 resource |
| **`T02_lead_enrichment_tests.tsv`** | lead: origin + Desert and regions tests | **main** |
| `T03_haplotype_group_members_by_region.tsv` | gene × group × region (incl. unassigned) | sheet 2 of the T04 resource |
| **`T04_haplotype_group_enrichment_tests.tsv`** | haplotype group: origin + tests (group and gene level) | **main** |
| `T05_lead_matched_background.tsv` | lead vs its 500 matched SNPs | supporting |
| `T06_top_region_vs_background.tsv` | trait × region: leads observed vs expected | supporting |
| `T07_background_top_region_by_MAF.tsv` | MAF class × region, all background SNPs | internal |
| `T08_lead_within_site.tsv` | lead: carriers vs site-mates, extreme accessions | supporting |
| `T09_lead_carrier_sharing.tsv` | pair of loci of one trait | supporting |
| `T10_accession_minor_allele_load.tsv` | accession | supporting |
| `T11_kinship_by_region.tsv` | region pair × pair type | supporting |
| `figures/Fig_carrier_origin.{png,tif}` | Fig. S1, 174 × 232 mm | main |
| `figures/Fig_desert_enrichment.{png,tif}` | Fig. S2, 174 × 207 mm | main |

Figures follow the TAG spec used in `08_USED_creating_figures` (174 mm wide, ≤ 234 mm, Liberation Sans/Arial
metric **9 pt** as Figs. 4–5 — TAG allows 8–12 pt at final size — 600 dpi, RGB, PNG + pixel-identical LZW TIFF).
Region is never shown by color alone (the Ch. 1 palette fails a color-vision check on its red/green pair; the
band names every region). Fig. S1 cells: minor/member `#1F2430` (black), major/other group `#F0D89E` (light tan),
no call/unassigned `#9AA1AB` (gray); OKLab lightness 0.26 / 0.89 / 0.71 vs white 1.00, so they separate in grayscale
and for color-vision deficiency; ΔE major vs no call 20.5, major vs white 13.7 (`CELL_*` in `00_config.R`).

---

## History

- **2026-09-25** built as a check for the Discussion: one-sided Desert test, three-panel figure
  (barcode + Desert share + observed vs origin-expected shift), compact table.
- **2026-09-30** reorganized to be cited from the M&M and the supplement (user request):
  - scripts split into `00`–`10` with one config; outputs renamed `T01`–`T11`; captions and a
    numbers file added;
  - **Desert test made two-sided** (user decision) and a **six-region test** added per lead; the
    per-lead Desert q ≤ 0.05 count fell from 13 to 8 (fiber) and 6 to 4 (β-glucan); the collective
    background result is unchanged;
  - **haplotype groups of the four genes added** (script 03, T03, T04, Fig. S1b, S2b);
  - **kinship by region added** (script 07, T11) for the "closely related" statement, previously computed
    ad hoc;
  - figure panels b and c dropped (user: readability); the barcode redrawn to the TAG spec and a separate
    Desert-enrichment figure added;
  - removed: `Table_COMPACT_carrier_origin.*` (→ T02 + T08), `Table_locus_origin_summary.tsv` (→ T02 +
    T08), `Table_locus_region_counts.tsv` (→ T01), `Table_trait_pooled_summary.tsv` (→ results_numbers.txt),
    `Table_accession_minor_load.tsv` + `Table_accession_genomewide_rare_allele_rate.tsv` (→ T10),
    `Table_locus_pair_overlap.tsv` (→ T09), `Table_matched_background_per_lead.tsv` (→ T05),
    `Table_top_region_vs_background.tsv` (→ T06), `Table_background_top_region.tsv` (→ T07),
    `Fig_carrier_origin.pdf`. TODO 403 in the manuscript named
    `results/tables/Table_COMPACT_carrier_origin.md`; **fixed 2026-09-30** in
    `build/UPDATED_Results_Discussion_Conclusions.md` (now T02 + T04; TODO 405 now points to T11; new
    TODO 417 on the "closely related" clause). The frozen copies (`04_discussion.md`, `Results_Discussion*.md`)
    still carry the old name.
  - later the same day (user): Fig. S1 cell colors changed — major and no-call were too similar
    (ΔE 12.3) and no-call too close to the white background (ΔE 4.6); now gray vs light tan (ΔE 20.5; 13.7
    from white). Lettering raised from 8 pt (the TAG minimum) to 9 pt in both figures, as in Figs. 4–5; the
    header rows of Fig. S1 given room so they are not clipped. Then, at the user's request, the two colors
    swapped: major = light tan, no call = gray. No number changed.

- **2026-09-30 (manuscript):** the Discussion was updated from this step with the user (review items 11, 11b, 11c in
  `build/REVIEW_2026-09-27_Results_Discussion.md`): 7H sentence reworded; Online Resource cited in §3 ¶2 and "Wild alleles"
  (TODOs 419, 420); "at most loci" added; the kinship clause dropped; the haplotype sentence rewritten on T04; a new GPAT6
  origin sentence. TODOs 403, 405, 417 and comments 412, 413 removed. The claims table above was rewritten to match.

## Not done (would need approval)

Re-testing the associations with region as a covariate or within region; a GWAS re-run; ancestral/derived
polarization (outgroup); a within-site test of the haplotype-group phenotypes.
