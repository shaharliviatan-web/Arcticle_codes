# Results and Discussion — structural blueprint

> **⚠ HISTORICAL — planning document, superseded (marked 2026-09-27, user decision).** This blueprint
> planned the Results and the Discussion before they were written. Both are now written, in
> **`../Methods_Results_Discussion_Conclusions.md`** (since 2026-09-30, the working file for M&M, Results, Discussion and
> Conclusions; the older `../Results_Discussion.md` named below was removed), and that file is the only source for what the manuscript says. Where this
> blueprint differs from it, the manuscript is right. Do not plan or write from this file.
> **Still useful:** the notes marked **⟶ M&M** or "M&M hand-over" (Ch. 2–4), which list what the Materials and
> methods section must state. They are a checklist for the M&M session, to check against the project files
> first. The Discussion's keep/drop decisions are in `DISCUSSION_CANDIDATES.md`.

**Wild-barley grain-quality GWAS, MorexV3 — manuscript for *Theoretical and Applied Genetics*.**
290 *Hordeum vulgare* ssp. *spontaneum* accessions, Southern Levant, three growing seasons,
four grain-quality traits (protein, starch, β-glucan, dietary fiber).

Prepared for review by S. Hübner · Results 2026-09-17, Discussion added 2026-09-24.
Numbers here are indicative; exact values are settled when each chapter is written.

**⟶ Discussion / ⟶ Ch. 3 / ⟶ Ch. 4** reported here, developed later ·
**[SUPPLEMENTARY]** Online Resource, not the printed text · **[OPTIONAL]** include-or-drop not decided.

What differs from the earlier mini paper is kept separately in **`CHANGED.md`**.

**File note (2026-09-27, user decision):** Results and Discussion are now one working file,
**`../Results_Discussion.md`**. It was built by joining `../03_results.md` + `../04_discussion.md` (verified
byte-identical). Those two are frozen on 2026-09-27 and kept for reference only. References to them below
are historical and now point to the matching part of `Results_Discussion.md`.

Status: **Ch. 1 locked · Ch. 2 written and approved 2026-09-22 (`../03_results.md`; 2 open TODOs: ESM number, Burton et al. citation) · Ch. 3 drafted 2026-09-23 and awaiting the user's review (`../03_results.md`; Fig. 4 approved; 4 open TODOs: 2 ESM numbers, GPAT6 and GH17 citations) · Ch. 4 optional, rebuilt 2026-09-24 on the MGmin = 3 haplotype analysis with Fig. 5, approved by the user 2026-09-24 (`../03_results.md`; 1 open TODO: GDSL citations, id 308), whether it stays in the paper is decided later · Discussion hand-over list, then the Discussion blueprint, at the end of this file.**

---

# Chapter 1 — Characterization of wild barley grain nutritional value traits

*Locked. The analysis behind it is unchanged — scripts and all 86 output files byte-identical — so
this chapter stands exactly as written, and is reproduced here only so the whole Results can be
judged as one piece.*

The four nutritional traits were quantified in 290 *H. spontaneum* accessions across three growing seasons, and for each trait genotype BLUPs were used to estimate accession-level genotypic values while accounting for season and block effects. Partitioning the phenotypic variance into genetic, season-and-block, genotype-by-season (G×E), and residual components (Figure 1A) showed that starch had the highest broad-sense heritability (H² = 0.474; Table 1), roughly twice as all other traits. The genotype-by-season interaction was small to moderate for all four traits (4.5% - 15.2%), where β-glucan combined the largest G×E component (15.2%) with the smallest genetic fraction indicating a strong genotype-specific response across the three growing seasons. Among morphological traits flowering time (H² = 0.747) had the highest heritability and tiller number the lowest (H^2^ = 0.031) indicating the respective genetic control for these traits. To further explore these trends, reaction norms were generated per genotype across seasons (Figure 1B). Overall, results were consistent with the heritability with substantial variation between specific genotypes indicating that the genotype X environment signature is notable for all traits.

  -----------------------------------------------
  **Trait**                    **H²**
  ---------------------------- ------------------
  *Starch*                     0.474

  *Fiber*                      0.248

  *Protein*                    0.242

  *β-glucan*                   0.234

  *Flowering time*             0.747

  *Spike length*               0.344

  *Grain weight*               0.223

  *Tillers*                    0.031
  -----------------------------------------------

![](../../old_mini_paper/md/media_Results/media/image1.png){width="6.427777777777778in" height="2.5708333333333333in"}**\
Table 1.** Broad-sense heritability (H²) of grain nutritional and morphological traits, estimated from the genotype × environment model as the genetic fraction of total phenotypic variance, for the four nutritional and four morphological traits, ordered by decreasing H².

**Figure 1.** Genetic architecture of grain nutritional and reference traits. (A) Variance partitioning for four nutritional and four morphological traits into genetic, season + block, genotype × environment, and residual components. (B) Per-genotype reaction norms across the three seasons for the four nutritional traits; grey lines, all accessions; colored lines, four representative trajectories.

To examine whether variation in grain nutritional traits followed the environmental gradient of the collection sites, trait BLUPs were first visualized across the 29 sites ordered by ecological region (Figure 2A). This descriptive pattern suggested that northern-derived genotypes generally had higher starch, whereas desert-derived genotypes tended to have lower starch and higher fiber and β-glucan; protein showed no clear directional trend. Among the nutritional traits, starch was strongly negatively correlated with fiber (r = −0.78; Figure 2B) and β-glucan (r = −0.53), while fiber and β-glucan were positively correlated (r = 0.35; all P \< 0.001). Thus, genotypes with higher starch tended to have lower fiber and β-glucan. Protein was poorly correlated with other traits, and significant negative correlation was obtained only for protein-starch (r = -0.12, P \< 0.05). Among the morphological and phenological traits, starch was strongly correlated with grain weight, plant height, and flowering time, while fiber showed the exact opposite pattern. Protein was also weakly negatively correlated with grain weight (r = −0.13, P \< 0.05). Overall, β-glucan and protein displayed weaker, mostly negative correlations, linking grain composition directly to whole-plant resource allocation and phenology (Figure 2C). To identify which environmental variables of the sites of origin were associated with each trait, site-mean nutritional-trait BLUPs were correlated with eight environmental variables (Figure 2D). Starch was positively correlated with clay, precipitation, and organic carbon, and negatively correlated with sand and pH. Fiber showed the opposite pattern for clay, sand, and precipitation, while β-glucan was positively correlated with sand and electrical conductivity and negatively correlated with clay, precipitation, and organic carbon. Protein was not significantly correlated with any of the measured environmental variables. Together, these associations linked the contrasting patterns of starch versus fiber and β-glucan to gradients in soil texture and precipitation, consistent with the broader geographic pattern observed across the collection sites.

![](../../old_mini_paper/md/media_Results/media/image2.png){width="5.954166666666667in" height="5.305555555555555in"}**Figure 2. Ecological architecture of grain nutritional traits.** (A) Centered trait BLUPs of the four nutritional traits across 29 sampling sites, ordered and coloured by ecological region. (B) Pearson correlations among nutritional-trait BLUPs. (C) Pearson correlations between nutritional- and morphological-trait BLUPs. (D) Pearson correlations between site-mean nutritional-trait BLUPs and eight environmental variables; cells outlined in black are significant at P \< 0.05. Significance (B--D): \*P \< 0.05, \*\*P \< 0.01, \*\*\*P \< 0.001.

---

# Chapter 2 — Identification of nutritional QTLs and candidate genes

*For approval.*

## Figures

* **Fig. 3** — four Manhattan plots beside four QQ plots, one row per trait. Locus members are
  painted in the locus colour over the grey background, so the figure shows the physical extent of
  each association, not only its peak.
* **Online Resource** — the 36 loci, one row each: position, span, member SNPs, allele frequency and
  carrier count, lead *P*-value. **[SUPPLEMENTARY]**

## What the chapter reports

**Mapping**

* BLUP phenotypes, three principal components, kinship, EMMAX, 7.1 M SNPs, 290 accessions.
* Population structure is adequately controlled — λ_GC close to 1 for all four traits.

**The loci**

* About fifty genome-wide-significant SNPs resolve into 36 loci.
* A locus is an LD block — members admitted on linkage disequilibrium with the lead, severed at
  internal gaps — not a fixed distance around a SNP.
* Fiber has the most loci, then β-glucan, then starch, then protein. Spans run from a single SNP to
  about 1.7 Mb, but the median is only about 90 kb — most loci are small and the widest is the
  exception.
* The strongest signal in the study is a β-glucan peak on 2H. Fiber's and starch's strongest peaks
  fall in the same 7H region, and protein's is on 3H.
* β-glucan's densest piece of architecture is a cluster of ten significant SNPs on 4H that resolves
  into two adjacent loci — its second (−log10p 7.33) and fourth (7.09) strongest signals; the third is
  3H:537,247,620 (7.19). [corrected 2026-09-22; was "second and third"]

**Architecture**

* Starch is by far the most heritable trait yet yields among the fewest loci; fiber is the opposite.
  ⟶ Discussion
* Most signals rest on rare alleles, carried by only a dozen or two accessions. The only two
  common-allele loci in the study are both starch. ⟶ Discussion

**Fiber and starch share genomic regions; β-glucan does not**

* Three regions carry fiber and starch loci together: **7H**, where the starch locus contains the
  fiber locus and the two leads sit a few hundred bases apart; **1H**, where a starch locus is
  nested between two fiber loci; and **3H**, where the widest fiber locus is followed by a starch
  locus under a megabase away.
* The 7H region — the tightest of the three — is one of the loci the gene search returned empty,
  because its LD span is only a few hundred bases. ⟶ Ch. 4
* β-glucan shares no region with any other trait — its nearest cross-trait locus is megabases away.
  ⟶ Discussion
* So the carbohydrate trade-off of chapter 1 has a genomic counterpart, but for starch and fiber
  only: β-glucan correlates with both phenotypically while sharing none of their loci.

**Two informative absences**

* **No β-glucan signal at the well-known *HvCslF6* locus on 7H.** Other barley GWAS report the same
  absence, and it is attributed to the gene carrying little variation — so it is not what drives
  β-glucan variation in this panel. ⟶ Ch. 3 and Discussion
* Fiber loci cluster on 1H, which carries a third of them on its own.

**Candidate genes**

* The search interval is the LD locus span itself, with no flanking window: the clump already is the
  region in LD with the lead.
* The 36 loci together cover about 7.9 Mb of search space, and return roughly 55 candidate genes at
  fewer than half of them.
* Fiber contributes most of the genes, then starch, while β-glucan returns only a handful from its
  eleven loci and protein a handful from its two — an ordering almost the inverse of the number of
  loci per trait.
* No lead SNP falls inside a gene body.
* Protein is the only trait with no gene-empty locus: both its loci contain genes. ⟶ Ch. 3
* 21 of the 36 loci contain no annotated gene, mostly because the interval is very small — four are
  single SNPs. 

**One locus to flag** ⟶ Ch. 3

* The widest locus, on 3H, spans about 1.7 Mb and holds a dozen candidate genes. Stated here so the
  constraint — those genes are one LD block and cannot be independent findings — is already on the
  table when chapter 3 reaches it.

## Open decisions

1. **Whether chapter 4 is written at all** — see that chapter.
2. **How hard to push the gene-desert point** — the empty-locus count depends on a gap parameter
   chosen as a reasonable value rather than from evidence. [update: drop it, dont write on it.]
3. **Whether the LD-decay result is kept and showen?**  +-200kb, and if so whether it moves to Materials and methods.
   [decided 2026-09-22: **LD decay is not mentioned in Results.**]
   [decided 2026-09-22: **no per-trait locus table in the body** — loci go to the Online Resource only.
   Protein's second locus sitting just above the threshold is not reported (unnecessary detail).
   Ch. 2 is written in TAG figure style ("Fig. 3a"); Ch. 1's "Figure 1A" is converted at the docx stage.]
   [decided 2026-09-22: **allele direction goes in** — one or two factual Results sentences, direction and
   lead-SNP beta only: the minor allele raises the trait at all 11 β-glucan and all 18 fiber loci, lowers
   it at all 5 starch loci, and protein is split (1 up, 1 down). Beta = effect of the minor allele (A1),
   from `Table_loci_master.tsv` (`lead_beta`). **The carriers' geographic origin / sampling sites are not
   used in the Results.** ⟶ Discussion: the rare-allele reading and the caveat that PC + kinship
   correction is weakest for rare, locally concentrated alleles (λGC does not detect it); to be settled in
   the Discussion session. **Not approved:** re-testing the effects with region as a covariate.]
4. [decided 2026-09-22: **Fig. 3 built and approved** — `08_USED_creating_figures/Figure_3/Fig3.{png,tif}`,
   by `make_figure_3.R`: per-locus colours, no key, no locus labels, uniform dots; QQ shows λGC only.]
   **Manhattan rendering** — the full per-locus colour key will not fit at print size with eight
   panels. Either a reduced two-colour scheme with locus identities carried by the Online Resource,
   or the full key as a supplementary figure.

---

# Chapter 3 — Haplotype analysis and candidate gene resolution

*Drafted 2026-09-23, awaiting review. Scope settled with the user before drafting:*
* *rank-vs-plausibility kept to two sentences (the two concrete examples only), not a section;*
* *median η² dropped; the crosshap software failure not reported; ε and MGmin values left to M&M;*
* *the β-glucan gene-desert fact cross-referred to Ch. 2, not repeated;*
* *"at most one gene can be causal" replaced by the plain LD-block statement — a locus in LD does not
  rule out more than one causal gene;*
* *added since this blueprint was written: the untestable genes fail for SNP density (21 of 25 with ≤ 5
  SNPs, 7 with none), the per-trait significant/tested counts (β-glucan 6/6, starch 5/5, fiber 12/19),
  the GPAT6–PHT4;3 1.47 Mb adjacency on 3H, and the elite lines being near-identical to one another.*

## Figures

* **Fig. 4 — built and approved 2026-09-23** — `08_USED_creating_figures/Figure_4/Fig4.{png,tif}`,
  by `make_figure_4.R`. One row per gene, violins left and aligned genotype barcodes right:
  **a**/**d** GPAT6 (fiber), **b**/**e** GH17 (fiber), **c**/**f** PHT4;3 (starch). Each violin
  carries its BH *q*, η² and the Wilcoxon-vs-largest-group brackets; each barcode carries one row
  per wild haplotype group, then the five elite cultivars under the same SNP columns.
  174 × 187.7 mm, 11 pt, 600 dpi. A stacked layout (violins above barcodes, as in the mini paper)
  was drawn first and rejected for height.
* Only positions called in **both** the wild and the elite data are drawn, so no genotype is assumed
  anywhere in the figure.
* **Re-rendered 2026-09-24 (user decisions):** true **9 pt** lettering (the approved file was ~7.3 pt,
  below TAG's 8 pt: `layout()` shrinks text to 0.66 and the script never reset it), and **consensus ties
  → REF** (step 07; only two GH17 cells in panel e change, from missing to REF). Brackets, y ticks, n
  labels and legend re-laid out to fit the larger text. No statistic or text number changed.
  **⟶ M&M:** the barcode consensus is the per-SNP majority over the group's members, ignoring missing
  calls, **with exact ties resolved to the reference allele**.
* **Fig. 3 has the same text-size bug** (~7.9 pt on the page instead of 12). Flagged 2026-09-24, not
  changed. Decide before submission.
* **Online Resource** — the significant genes with their haplotype statistics and annotation
  evidence. **[SUPPLEMENTARY]**

## What the chapter reports

**Why haplotypes**

* Grouping accessions by their haplotype across a whole gene asks a question single-SNP tests
  cannot: whether variation in the gene tracks the trait, and in which direction and by how much.

**The analysis**

* Fixed grouping parameters for every gene — one grouping, one test, one correction — so nothing is
  tuned per gene.
* Of the candidate genes, about half are testable; **23 are significant**, and they fall in **12
  loci**.
* Roughly half the untestable genes fail for lack of data in the window.
* Effect sizes are reported beside every p-value, because a p-value at a locus that was selected for
  association with that very phenotype is close to guaranteed. What a reader can judge is the share
  of variance explained and the gap between the highest and lowest haplotype group.

**Protein drops out entirely**

* Protein has candidate genes, but **not one of them is testable**. Protein therefore contributes no
  gene to this study. ⟶ Discussion
* shahar needs to check the method and quanitificaion chimical method of protein in the study, it may help explain ⟶ in the discussion ⟶ why no genes found for protein. 

**Genes are not independent findings**

* The 23 significant genes sit in only 12 loci, and one locus alone holds 7 of the 12 fiber genes.
  Genes within a locus are in one LD block, so **at most one can be causal** — loci are reported
  alongside genes throughout.

**What the significant genes are**

* Once annotated, most of the 23 genes offer no plausible route to their trait — they include
  transposable-element-derived, defence-related and housekeeping genes. Only a few can be connected
  to grain composition directly from their function.
* Statistical rank and biological plausibility point different ways: the top gene by significance
  has no functional annotation at all, and the second carries a transposase-derived domain.
* Selection was therefore two-stage: **significance was the gate** — only genes significant in the
  haplotype analysis were considered at all — and among those, the genes carried forward were chosen
  on biological relevance, from protein function and published literature, rather than on their
  rank by significance. ⟶ Discussion

**The three genes carried forward**

* Two fiber genes and one starch gene, out of 23. Each is reported with its mechanism and with the
  axis on which it is weakest:
  * **GPAT6** (fiber) — makes the committed precursor of cutin, which is recovered in the insoluble
    fibre fraction. Best significance of the three, but the signal is **one deviant haplotype above
    a flat background**, and under half the panel is grouped.
  * **GH17 glucanase** (fiber) — from the enzyme family that degrades barley mixed-linkage glucan,
    the dominant soluble fibre of the grain. Cleanest haplotype structure and largest effect, but it
    is an **uncharacterised paralog, not either known endohydrolase**, and it has the lowest
    assignment rate in the study.
  * **PHT4;3** (starch) — a plastid phosphate transporter, and phosphate is the allosteric brake on
    the committed step of starch synthesis. Much the best data quality, closest to its lead SNP and
    best protein identification, but the **weakest effect**, and the documented starch phenotype
    belongs to a **different family member**, so the mechanism is a family-level inference.


**β-glucan contributes no gene — a reportable negative** ⟶ Discussion

* **Revised 2026-09-26 (user):** the Ch. 3 paragraph now reports results only: no gene carried
  forward; all six testable genes significant, with small effects. The canonical-pathway distances
  moved to the Discussion (§5). "Informative rather than technical", the named proteins and the
  "weaker half of its signal" sentence were removed.

* Its three strongest loci return no annotated gene at all, so its genes come from the weaker half
  of its signal.
* No canonical β-glucan gene lies near any β-glucan lead.
* Every β-glucan effect size is small, well below the two fibre genes.
* β-glucan has real genetic architecture in this panel — more loci than starch, and stronger peaks —
  but it sits away from the annotated genes and away from the known pathway.

**Elite cultivars carry the wild haplotypes** ⟶ Discussion

* Five spring malting cultivars released 2012–2018, genotyped across the same windows.
* **⟶ M&M (decided 2026-09-23):** that the clustering never saw the elite lines and that no elite line is
  assigned to a group was written in the Results, then dropped as method. It is now in **neither the
  Results nor the Fig. 4 caption**, so the Methods session must state it — it is what rules out reading
  the elite/wild match as circular.
* **⟶ Discussion (decided 2026-09-23):** the five cultivars are **near-identical to one another** at all
  three genes — monomorphic at 42/43, 33/33 and 12/12 of the drawn columns
  (`07_.../results/tables/Table_site_overlap.tsv`). They are therefore **not five independent
  observations**, but one modern spring-malting haplotype per gene. Stated in the Results at first and
  dropped from it (visible in Fig. 4d–f); it must be picked up in the Discussion, together with the
  "carry, not selected" caveat above.
* **GPAT6 is the informative case**: all five cultivars match the low-fibre wild haplotype closely
  and the high-fibre haplotype hardly at all. **The high-fibre wild haplotype is essentially absent
  from modern malting germplasm.**
* For the other two genes the elite lines simply carry the **Morex reference haplotype** — and Morex
  is itself an elite cultivar, so the apparent match is reference identity.
  **Only GPAT6 supports a shared-haplotype claim.**
* **⟶ Discussion (decided 2026-09-23):** that whole reading is interpretation and was taken out of the
  Results. The Results now report only the bare percentages — elite reference alleles 12–13% at GPAT6
  against 14% / 12% / 65% for its three wild groups, and 100% at GH17 and PHT4;3. The Discussion must
  supply the argument: Morex is itself an elite cultivar, so reference identity is not shared ancestry,
  and **only GPAT6 supports a shared-haplotype claim**.

**Must be stated in the text**

* ~~The figure shows **which haplotypes elite cultivars carry — not that breeding selected them.**~~
  **Dropped from the Results 2026-09-23 (user decision).** It is an interpretation of what the comparison
  can and cannot support, not a result. **⟶ Discussion: pick it up there** — the caveat itself still holds
  and is recorded in `07_.../README.md` § "Caveat for the manuscript". The Results now close on the plain
  fact that the five cultivars are near-identical to one another at all three genes.
* Genes within a locus are in LD and are not independent discoveries. **Written 2026-09-23 as the plain
  fact only** — the sentence declaring that loci would be reported alongside genes was dropped; the text
  does it instead.

## Decided

* **All three presented genes get a figure panel** — GPAT6, GH17 and PHT4;3 — each with the
  elite-line barcodes beneath it.
* **The elite comparison uses the shared-sites version**: only positions called in both data sets
  are drawn. This costs columns but assumes nothing.
* **The genes not carried forward appear in the supplementary annotation table**, with their
  statistics and annotation. The reasoning for setting each aside is not published.

---

# Chapter 4 — A shared fiber/starch signal at 7H · **[OPTIONAL]**

**[OPTIONAL] — written and approved by the user 2026-09-23 (`../03_results.md`). Whether it stays in the paper is decided later.**
Everything below rests on an exploratory branch (first version, now archived in `03_01_7H_branch_Starch_Fiber_shared_signal_explore/archive_v1_MGmin2_eps0.6/`,
reproducible) that **feeds nothing**: no main table has been changed.
Full detail and caveats: `section_materials/03_results/NOTE_7H_shared_signal_GDSL_gene.md`.

*Why it is a separate chapter and not part of Ch. 3:* this gene did not come through the haplotype
pipeline and cannot be presented as a fourth member of that set. It is a targeted follow-up of an
observation made in Ch. 2, with its own method, and should be structurally separate so that is
obvious.

**Decided: the gene is not adopted into the main pipeline.** It stays outside the span-based
candidate set and outside every downstream step. The search rule is not changed for it, and nothing
is re-run. It is reported here as a single locus examined on its own, or not at all.

## 2026-09-24 — rebuild on the MGmin = 3 haplotype analysis (in progress)

* **New source analysis (user decision):**
  `03_01_7H_branch_Starch_Fiber_shared_signal_explore/` (the branch root since 2026-09-24; it was the subfolder
  `REPLACEMENT_ANALYSIS_MGmin3_eps0.9`, and the MGmin 2 version moved to `archive_v1_MGmin2_eps0.6/`), crosshap at
  MGmin = 3, ε = 0.9, which groups the gene into two haplotypes (212 | 34) and is significant for both
  traits. It replaces the MGmin = 2 negative result the 2026-09-23 draft reported. **Chapter text
  rewritten 2026-09-24** (`../03_results.md`), three paragraphs + Fig. 5, **approved by the user 2026-09-24**.
* **Fig. 5 drafted 2026-09-24, awaiting approval:** `08_USED_creating_figures/Figure_5/Fig5.{png,tif}`,
  by `make_figure_5.R`. Fig. 4's style (true 9 pt, 70 | 104 mm). Violins **a** fiber and **b** starch
  stacked left, **one** shared barcode **c** right (the two traits give identical groups). Label: raw
  KW *P* and η² (single gene, no BH family). Red triangles under the **three genome-wide significant
  SNPs** only (the fourth haplotype-defining SNP, 7H:573,606,282, is not significant and not marked).
  Legend moved under the barcode inside column c (174 × 119.4 mm). Approved by the user 2026-09-24.
* **Scope of the rewrite, decided with the user 2026-09-24** (all numbers from
  `03_01_7H_branch_Starch_Fiber_shared_signal_explore/results/tables/` unless noted; ε 0.9 only since 2026-09-24):
  * **Kept from the 2026-09-23 draft:** the gene and its position (578 / 732 bp, 554 bp beyond the
    starch locus, next gene 112 kb); the family-level secreted GDSL annotation; the three significant
    SNPs 578–763 bp upstream and the six null gene-body SNPs; the closing sentence, with "allele" →
    "haplotype".
  * **New, replacing the negative haplotype test and the genotype split:**
    * two haplotypes over 246 of 290 accessions, identical for both traits; minority haplotype B
      (n = 34) higher fiber (*P* = 2.1 × 10⁻³, η² = 0.035, 0.54 SD) and lower starch (*P* = 1.1 × 10⁻³,
      η² = 0.040, 0.56 SD); Fig. 5a, b;
    * the haplotypes are defined by the only SNP group formed: the three significant SNPs plus
      7H:573,606,282 (r² 0.44–0.54 with them);
    * haplotype B is exactly the 34 accessions carrying the minor allele at both leads (A = the 167
      major-at-both + 45 with discordant or missing lead calls). This replaces the genotype-split
      passage (201 of 205, ±0.60 SD, 85 missing), which is dropped;
    * ~~parameter disclosure in the Results, option a~~ **→ option b, 2026-09-24 (user): M&M only.** The
      Results report only the result, with no pointer to M&M (removed 2026-09-24, user). What M&M must say: at the Ch. 3
      settings the associated SNPs were not grouped (*P* = 0.16 / 0.36); MGmin raised to 3, the smallest
      value that groups them; procedure in M&M;
    * elite cultivars (Fig. 5c): at the three significant SNPs 11 of 15 elite genotypes called, all
      reference; no cultivar carries haplotype B's alleles.
  * **Family-function sentence added 2026-09-24 (user):** one sentence, as for each Ch. 3 gene —
    GDSL esterases act on cell-wall and cuticle esters (rice DARX1 deacetylates arabinoxylan; tomato
    GDSL1/CD1 polymerises cutin). Citations are TODO 308, to verify in the References session. The
    rest of DISCUSSION item 10 (BS1, barley hull attachment) and item 11 (starch link) stay in the Discussion.
  * **Optional items in:** effect size smallest of the four presented genes (η² 0.035–0.040 vs
    0.078–0.222); assignment rate highest of the four (246 vs 205 / 123 / 107; step 04 `gene_results.tsv`).
  * **Out (user decisions):** the ε 0.4–0.8 robustness plateau (not in the Results); the hypothetical
    in-context BH *q* (0.0042 / 0.0013); the GWAS lead betas (redundant with the haplotype direction);
    elite agreement % with A vs B (items above already carry the point).
  * **⟶ M&M:** the parameter-selection procedure, stated as GWAS-informed rather than phenotype-blind
    (step 1: smallest MGmin grouping the GWAS signal SNPs → 3; step 2: ε maximising assignment → 0.9;
    no haplotype-test *P* used; README § 2 of the replacement analysis); raw *P* (one pre-specified
    gene, no BH); the consensus tie rule (ties → REF). The window is gene ± 1 kb for wild and elite
    alike; the ±2 kb in `02_fetch_elite_vcf.sh` is only download padding, trimmed to the window, and
    is **not** a method detail to report.

## Scope decided with the user (2026-09-23), before drafting

* **Short, complementary to Ch. 2–3, no figure of its own.** The text cites Fig. 3b, d (the 7H
  peaks). A figure may be added later by the user. The exploratory crosshap and heatmap figures were
  rendered at a non-standard ε and must not be used.
* **In:**
  1. the gene: `HORVU.MOREX.r3.7HG0729030`, its position beside both loci, 554 bp outside the starch
     locus, and the 112 kb to the next gene;
  2. the annotation, a secreted GDSL esterase/lipase, **stated as a family-level call**, in one clause
     (domain-database detail cut 2026-09-23, user; evidence kept in the md src comment);
  3. the local association profile: the associated SNPs are upstream of the TSS and the gene-body
     SNPs are null, **as a fact only**;
  4. one allele with opposite effects: lead concordance (201 of 205), the effect of the minor allele
     in SD with η², **no p-values** (user decision; the genotype-split p is not kinship-corrected),
     and the GWAS lead betas for direction;
  5. the missing calls: 85 of 290 accessions lack a call at one or both leads;
  6. the Ch. 3 haplotype test is negative here, in one sentence, with the reason (the SNPs carrying
     the association were not grouped). **The ε sweep is never reported**. Its raw *P* = 0.16 (fiber) and
     0.36 (starch) **stay in the text** (user decision 2026-09-23): the no-p rule in item 4 covers the
     genotype split only;
  7. ~~a closing synthesis: three shared regions resolved at three levels~~ **Drafted, then moved
     to the Discussion 2026-09-23 (user):** it repeated Ch. 2–3, its only new fact (the 1H SH3-domain
     gene) is already in Online Resource 2, and the framing itself is interpretation. See DISCUSSION
     item 13. Ch. 4 now ends on the negative haplotype test.
* **Not reported (user decisions 2026-09-23):**
  * the GDSL family's routes to fiber (DARX1, BS1, CD1/GDSL1, barley hull attachment) as a Results
    sentence. ⟶ Discussion;
  * that Morex carries the major (low-fiber, high-starch) allele at both leads. Not reported
    anywhere: one reference cultivar is not enough to say anything about cultivated barley;
  * a GPAT6 × PHT4;3 carrier-overlap check at 3H: not run;
  * ~~the kinship-corrected effect / geographic clustering of minor-allele carriers~~: not written,
    by choice (see Open items 2).
* **The flank rule is not defended in the text** (CLAUDE.md, "write for the reader"). The Results
  state the fact instead: the 7H loci contain no gene, and the gene lies 554 bp beyond the starch
  locus.
* **⟶ M&M:** one sentence on the direct genotype split: the same Kruskal–Wallis effect-size measures
  as the haplotype analysis (η², difference in SD), applied to groups defined by the observed
  genotype at the two lead SNPs (accessions with discordant or missing calls excluded). The InterPro
  and signal-peptide checks (SignalP, Phobius) also need a line where the annotation chain is
  described.

## What the chapter reports

* Returning to the tightest of the three shared regions: a single gene sits immediately beside both
  leads, in a gene desert of more than a hundred kilobases. There was no other candidate to choose.
* It encodes a secreted GDSL esterase/lipase. Several databases agree on the family, and the protein
  is predicted to carry a signal peptide. **Family-level call, not a named ortholog.**
* The gene is on the minus strand, so its promoter runs into the signal: the associated SNPs sit
  upstream of the transcription start, while every SNP inside the gene body is null in both traits.
  The reading that this excludes a coding change and points to regulatory variation is
  ⟶ Discussion.
* One allele raises fiber and lowers starch by a similar amount, a near-symmetric inverse.
  ⟶ Discussion (one causal variant read out by two correlated traits).
* The haplotype test used in Ch. 3 is negative at this gene.
* ~~Three shared regions, three levels of resolution.~~ Moved to the Discussion (item 13).

## Open items before this can be written

1. **Direction of effect — settled 2026-09-22.** The minor allele at 7H raises fiber and lowers starch,
   in both the kinship-corrected GWAS (`lead_beta`, minor allele: fiber +0.119, starch −2.145) and the
   direct genotype split (+0.60 / −0.62 SD). `lead_beta` is reported for the minor allele (A1)
   (`01_.../results/00_FINAL_BLUP_3PC/README.md`, beta sign convention).
2. ~~**The effect sizes are not kinship-corrected.**~~ **Not written, by choice (user decision
   2026-09-23).** Not checked and not stated in the manuscript: no geographic-clustering check of the
   minor-allele carriers, no corrected effect estimates.
3. High missingness at the two leads, and the family-level annotation, both need stating.
   **Done in the 2026-09-23 draft.**

---
---

# DISCUSSION — items handed over from Results

*Added 2026-09-23 (user decision). One list of everything the Results sessions marked
"⟶ Discussion", so the Discussion session starts from it. Each item says where it came from; the
longer background stays in `section_materials/04_discussion/`, and each item points to it where it
exists. Items stay also where they were first marked, above.*

## From Ch. 2 (GWAS, loci, candidate genes)

1. **Mapping power does not follow heritability.** Starch is the most heritable trait but yields
   5 loci; fiber has about half its H² and yields 18. (Ch. 2, Architecture)
2. **Most signals rest on rare alleles** (MAF < 0.10 at 25 of 36 leads); the only two common-allele
   loci are both starch. Pick up with the caveat that PC + kinship correction is weakest for rare,
   locally concentrated alleles (λGC does not detect it). **Not approved:** re-testing with region
   as a covariate. The carriers' geographic origin is not used in the Results. (Ch. 2, Open
   decisions 3)
3. **β-glucan shares no genomic region with starch or fiber**, although it correlates with both
   phenotypically. (Ch. 2)
4. **No β-glucan signal at *HvCslF6* on 7H**, and no canonical β-glucan gene near any β-glucan lead.
   Background, literature and open checks: `section_materials/04_discussion/NOTE_cslf6_7H_absence.md`.
   (`NOTE_canonical_genes_not_recovered.md` in the same folder is background only, marked not for
   writing.) (Ch. 2, Ch. 3)

## From Ch. 3 (haplotypes, genes, elite lines)

5. **Protein contributes no gene** (no candidate gene testable). The user is to check the protein
   quantification method (NIR / Kjeldahl calibration), which may help explain it. (Ch. 3)
6. **Two-stage gene selection:** significance was the gate and biological relevance chose among the
   23. Statistical rank and plausibility diverge (the top gene has no functional call, the second is
   transposase-derived). (Ch. 3)
7. **β-glucan contributes no gene: a reportable negative.** Its strongest loci are gene-empty, its
   effects are small, and it sits away from the canonical pathway. (Ch. 3)
8. **The elite cultivars.**
   * The five cultivars are near-identical to one another at all three genes (monomorphic at 42/43,
     33/33, 12/12 drawn columns, `07_.../results/tables/Table_site_overlap.tsv`): one modern
     spring-malting haplotype per gene, not five independent observations.
   * Morex, the reference genome, is itself an elite cultivar, so the 100% reference identity at GH17
     and PHT4;3 is not shared ancestry. **Only GPAT6 supports a shared-haplotype claim.**
   * The comparison shows which haplotypes elite lines carry, **not that breeding selected them**
     (`07_.../README.md` § "Caveat for the manuscript").

## From Ch. 4 (the 7H shared signal) · only if Ch. 4 stays in

9. **Regulatory, not coding.** At `7HG0729030` all associated SNPs are upstream of the TSS and the
   six gene-body SNPs are null, which excludes a coding change in this gene and points to regulatory
   variation. Caveats: tag SNPs only (the causal variant may be an ungenotyped correlated one), and
   the promoter reading depends on the annotated TSS.
   (`03_01_7H_branch_Starch_Fiber_shared_signal_explore/archive_v1_MGmin2_eps0.6/results/tables/local_association_profile.tsv`)
10. **How a GDSL esterase/lipase can reach grain fiber.** Three routes in the family: arabinoxylan
    deacetylation (rice DARX1; also BS1), cutin polymerisation (tomato CD1/GDSL1, the same logic as
    GPAT6), and hull–caryopsis attachment in barley. The protein is secreted (signal peptide), i.e.
    targeted to where these processes happen. Family-level call only. Citations still to find and
    verify. (`section_materials/03_results/NOTE_7H_shared_signal_GDSL_gene.md`)
11. **The starch link is compositional, not mechanistic.** GDSL esterases are not starch enzymes; if
    the cell-wall/hull fraction shifts, starch as a share of the grain shifts the other way.
12. **One haplotype, two traits** (updated 2026-09-24 for the MGmin = 3 analysis). The same minority
    haplotype B (n = 34) carries higher fiber (0.54 SD) and lower starch (0.56 SD), a near-symmetric
    inverse, consistent with one causal variant read out by two correlated traits.
13. **The fiber/starch trade-off at the genome level: three shared regions, three levels of
    resolution.** Drafted as the closing paragraph of Ch. 4 and moved here 2026-09-23 (user). The
    facts, with sources:
    * **7H:** one gene (`7HG0729030`) beside both loci; one allele with higher fiber and lower starch
      (Ch. 4).
    * **3H:** neighbouring loci, 595 kb apart, each with one presented gene: GPAT6 (fiber) and PHT4;3
      (starch), 1.47 Mb apart (Ch. 2, Ch. 3, Fig. 4a, c).
    * **1H:** the three fiber loci contain no annotated gene
      (`03_.../results/tables/candidate_genes.tsv`). The starch locus holds two genes. `1HG0049620`
      is not testable (4 SNPs; `04_.../Stats/genes_not_tested.tsv`). `1HG0049610` is an SH3
      domain-containing protein, haplotype-significant (q = 0.026, η² = 0.026;
      `04_.../Stats/gene_results.tsv`, `05_.../Table_significant_genes_paper.tsv`), with no known
      role in either trait (SH3P2-like, endocytosis / autophagosome formation; `TRAIT_CANDIDACY.md`,
      starch "all too indirect").
    * The carbon-allocation / shared-control reading belongs here. If Ch. 4 is dropped, this rests on
      Ch. 2's co-localization and Ch. 3's GPAT6 / PHT4;3 adjacency only.
14. **Elite cultivars at the GDSL gene** (added 2026-09-24, only if Ch. 4 stays in). At the three
    significant SNPs every called elite genotype is reference (11 of 15 called), so the cultivars carry
    the low-fiber / high-starch alleles at the SNPs that carry the association. Elsewhere in the window
    they carry alternate alleles (53–75% reference over the 22 shared sites), so this is **not** plain
    reference identity across the gene. Read beside GPAT6 (item 8): at both fiber genes the high-fiber
    wild haplotype is absent from these malting cultivars.


---
---

# DISCUSSION — structural blueprint

*Drafted 2026-09-24 with the user. Built from the **written** Results chapters (`../03_results.md`)
and from the hand-over list above; the mini paper's Discussion
(`old_mini_paper/md/06_discussion.md`) is the model for structure, tone and length. Sariel's three
sections are kept where they survive — §1 and §2 intact, §3 with every example replaced — and §4–§6
are new. Bracketed numbers cite the hand-over items above.*

**Discipline for this section:** the Discussion interprets, it does not restate. Numbers appear only
where they carry an argument. Sariel's §3 is four sentences long; that is the target density.

## 1. Ecological differentiation of grain composition

* Extends the known ecological structuring of this collection — population structure, flowering
  time, growth (Hübner et al. 2009, 2013) — to grain composition.
* Because the comparisons rest on common-garden BLUPs, the regional pattern is heritable population
  difference, not plastic response. This is what licenses everything downstream.
* The direction has a physiological basis: heat and drought during grain filling shorten the filling
  period and suppress starch deposition (Savin & Nicolas 1996).
* β-glucan's unusually large G×E fits its known environmental lability, where the direction of
  response depends on genotype and stress timing (Anker-Nilssen et al. 2008).

## 2. A carbon-allocation axis, now visible in the genome

* The trade-off is no longer only phenotypic: it recurs in the direction of allelic effects, in locus
  position, and at a single haplotype — three independent lines of evidence from different data. [13]
* **Three shared regions, three levels of resolution** [13]: 7H resolves to one gene carrying both
  effects; 3H to two neighbouring genes, one per trait; 1H to nothing — its fiber loci contain no
  annotated gene and its starch locus offers none with a known role. The axis is visible at every
  level the data can reach, and where it is not, nothing contradicts it.
* Read as competition for one sucrose pool in the developing endosperm; *HvCslF6* overexpression
  raising β-glucan while lowering starch is the functional precedent (Lim et al. 2020).
* Grain composition is one expression of a whole-plant strategy tied to phenology and grain size,
  not an endosperm-autonomous property.
* **Protein is governed separately** — carbon–nitrogen balance rather than carbon partitioning
  (Corke et al. 1989; Bogard et al. 2010; Simmonds 1995). Practical consequence: raising fiber or
  β-glucan need not cost protein.
* Protein also contributes no gene [5]. Before interpreting that biologically, the protein
  quantification method (NIR calibration against Kjeldahl) is to be checked — it may be part of the
  explanation. **User action, still open.**

## 3. Genetic architecture

* Heritability did not predict mappability: three traits of near-identical heritability yielded 2, 11
  and 18 loci. Mapping success reflects allele frequency, effect size and LD rather than total
  genetic variance. [1]
* What this panel detects is largely low-frequency variation — the kind cultivated germplasm lacks,
  but also the kind whose effects a single panel estimates poorly. **Carry the caveat** [2]: PC and
  kinship correction is weakest for rare alleles concentrated in a few related accessions, and λ_GC
  does not detect that failure mode. Carrier geography is deliberately not used in the Results, and
  re-testing with region as a covariate was not approved.
* Statistical rank did not track biological plausibility, so significance alone is a poor guide to
  candidacy. [6]
* The genes carried forward point to four distinct routes into the grain: cutin deposition, glucan
  hydrolysis, plastid phosphate transport, and cell-wall ester modification. *(Replaces Sariel's
  closing four-pathway sentence, whose genes are gone.)*

## 4. What this design can and cannot resolve

*Both limitations are conservative: they lose genes, they do not invent them. Everything reported
passed; what is uncertain is what was missed.*

* The unit of discovery is the locus, not the gene: within an LD block the candidates are
  alternatives, not accumulating evidence.
* The gene list is a function of the locus definition — a wider window returns proportionally more
  genes — so the absence of a candidate is weak evidence of absence.
* The GDSL esterase is this paper's own worked example: a gene beside one of the strongest signals,
  excluded by the rule, recovered only by looking. [9]
* Haplotype testing needs variation the gene may not carry; one global clustering setting cannot suit
  genes of differing SNP density, and a gene with a single SNP is untestable by construction even if
  that SNP were causal.
* Where a gene was recovered, what it shows is bounded too: the 7H signal is regulatory rather than
  coding, but on tag SNPs and an annotated TSS [9], and its starch link is compositional rather than
  mechanistic [11].

## 5. β-glucan: strong signals, no genes

* β-glucan has the strongest associations in the study and yields no candidate gene. That contrast is
  itself the finding. [3, 7]
* **It is not the limitation of §4**: the canonical-gene distances are measured from lead SNPs and do
  not depend on the locus definition, and no β-glucan locus lies on 7H at all.
* **Why no signal at the canonical genes — carried by the literature, not asserted** [4]. *HvCslF6*
  is essential for mixed-linkage glucan synthesis and correspondingly conserved: across 1,336
  accessions, 288 of them wild, only three coding SNPs segregate and **none is associated with grain
  β-glucan content**, which the authors attribute to the gene's indispensable role (Garcia-Gimenez
  et al. 2019). A gene can be necessary for a trait and still carry no detectable association,
  because association requires segregating functional variation — and here there is almost none.
* **Our result is the common outcome, not an anomaly**: Geng et al. (2021) report the same failure
  and propose the same explanation, and Houston et al. (2014) found no *CslF6* association in 603
  cultivars; Shu & Rasmussen (2014) did detect the region in 254 European spring barleys, so
  detection is population-dependent. Background and open checks:
  `section_materials/04_discussion/NOTE_cslf6_7H_absence.md`.
* **Dropped deliberately:** any claim that the causal variation is regulatory or intergenic. With no
  β-glucan gene recovered there is no evidence for it.

## 6. Wild alleles for cultivated barley

* **Two genes, two traits, the same direction** [8, 14]: at GPAT6 the high-fiber haplotype and at the
  7H GDSL esterase the fiber-raising, starch-lowering haplotype are both essentially absent from the
  five malting cultivars. That convergence across independent genes is the substance of the claim.
* **Consistent with selection during malting breeding**, which targets high extract and low wort
  viscosity and so favours high starch and low fiber and β-glucan. A plausible cause, and one of this
  paper's central ideas.
* **Not demonstrated.** The comparison shows what cultivars carry, not what breeding chose; and the
  five cultivars are near-identical to one another at every gene [8] — one modern spring-malting
  haplotype, not five independent observations — so selection cannot be separated from drift or
  founder effects.
* At the glucanase and the phosphate transporter the cultivars carry the reference genome's alleles,
  and Morex is itself an elite cultivar, so those two genes are uninformative here. [8]
* Narrowing broad associations to a few plausible targets is what makes introgressed segments
  monitorable (Hübner & Kantar 2021).
* Next steps: validation in elite backgrounds, replication in independent populations, functional
  tests.

## Decided with the user 2026-09-25 (Discussion session)

* **Status 2026-09-26:** `../04_discussion.md` holds §1–§3 (approved by the user 2026-09-26) and
  §4–§6 (drafted 2026-09-26, reviewed and edited with the user 2026-09-27; the Discussion is finished for now). Open TODOs 401–406, 408 (407 removed with its paragraph 2026-09-27). The Conclusions are not drafted yet.

* **The item-by-item keep/drop decisions for the Discussion are in `DISCUSSION_CANDIDATES.md`**
  (2026-09-26). It overrides the section bullets below where they differ. It also holds the drafting
  rules: interpretation only, and a dir-11 Word comment on every inference from the carrier-origin check.

* **Ch. 4 is approved and stays in the paper.** Hand-over items 9–14 are usable.
* **Protein measurement method: not used in the text.** Insert only a TODO Word comment at the protein
  paragraph: check how protein was quantified; it may help discuss the protein GWAS and gene results.
* **Wild-vs-cultivated composition** rests on the literature only (Friedman and Atsmon 1988), with a
  TODO Word comment asking whether our own experiment holds such data (flag 1 below). Not discussed
  further in the text.
* **Conclusions = the last subsection of the Discussion** (TAG lists no Conclusions section; headings
  ≤ 3 levels). One paragraph, from the mini paper's `07_conclusions.md`; it takes over §6's closing
  breeding outlook (Hübner and Kantar, next steps). **Not drafted until the user gives the green light.**
  The mini paper's β-amylase example (a wild allele introgressed successfully) goes in with a TODO Word
  comment: citation to find.
* **"Domestication and breeding shifted cultivated barley to the starch-rich end"** (Friedman and
  Atsmon 1988) **stays in §2**, as in the mini paper.
* **Pending a separate check (user, new session):** minor alleles push every locus one way (fiber
  18/18 and β-glucan 11/11 up, starch 5/5 down). Carrier origin by region (North, Coast, Desert,
  HZ1–3) is being checked, to decide between the carbon-allocation reading and residual population
  structure. Items 1–2 wait for that result.

## Flags for the drafting session

1. **No wild-vs-cultivated nutritional comparison is possible from our data.** Morex and Clipper were
   grown in the common garden but carry **no NIR values** in
   `00_THIN_.../all years barley.csv` (checked 2026-09-24: `NA` at all four traits, both seasons, all
   17 check rows). Friedman & Atsmon (1988) must carry that claim. **Insert a TODO comment at that
   sentence** asking whether cultivated-check nutritional data exist elsewhere in the project or in
   the collaborators' records — if they do, the claim can be made from this experiment instead of a
   1988 citation, which would be considerably stronger. See
   `section_materials/04_discussion/NOTE_cultivated_checks_no_NIR.md`.
2. **Dead claims that must not resurface:** *Pho*, AP2/ERF and BAHD as candidates; the three 7H
   dual-trait genes (`7HG0729020/090/100`); "β-glucan distributed across all seven chromosomes";
   "starch required haplotype-level resolution" (checked and dropped — see the comment at
   `../03_results.md`, the haplotype-analysis paragraph of Ch. 3).
3. **§4 must come before §5.** The β-glucan argument has to be read as surviving the limitations, not
   as a casualty of them.
