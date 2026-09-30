# 06_USED_genes_selected_to_present

The **three candidate genes carried forward** out of the 23 haplotype-significant genes
of the wild-barley grain-quality GWAS — two for **fiber**, one for **starch**.

290 wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*, Southern Levant),
Morex V3. Step-04 run `loci_LDspan_eps06_V4`.

Created 2026-09-10; reduced to the three strong candidates 2026-09-10.
Figures are unmodified copies, verified byte-identical to source.

**Full reasoning for all 23 genes → [`../05_USED_gene_annotation_analysis/08_USED_annotation_master/TRAIT_CANDIDACY.md`](../05_USED_gene_annotation_analysis/08_USED_annotation_master/TRAIT_CANDIDACY.md)**

---

## The three, at a glance

| | **GPAT6** | **GH17 glucanase** | **PHT4;3 transporter** |
|---|---|---|---|
| gene | `3HG0301300` | `5HG0487060` | `3HG0301710` |
| trait | fiber | fiber | starch |
| locus | `fiber_L08` | `fiber_L14` | `starch_L03` |
| position | 3H:545.000–545.005 Mb | 5H:462.727–462.731 Mb | 3H:546.474–546.478 Mb |
| lead SNP (−log10p) | 3H:544084309 (6.24) | 5H:462564771 (6.05) | 3H:546433616 (6.34) |
| distance to lead | 916 kb | 162 kb | **40 kb** |
| nearest member SNP | **20 bp** | 1.5 kb | 311 bp |
| accessions grouped | 123 / 290 | 107 / 290 | **205 / 290** |
| haplotype groups | 3 | **5** | 4 |
| KW p (raw) | 1.38e-06 | 2.36e-05 | 3.06e-04 |
| **BH q** | **5.26e-06** | 7.47e-05 | 3.83e-04 |
| rank by q (of 23) | 7 | 9 | 11 |
| **η²** | 0.208 | **0.222** | 0.078 |
| Δ top–bottom | **1.51 SD** | 1.41 SD | 1.31 SD |
| identification | 69.3% id / 97.8% cov | 62.7% id / 90.9% cov | **84.2% id / 99.8% cov** |

Each PDF is the merged per-gene figure set (10–12 pages): the crosshap combined plots
followed by the LD/haplotype heatmaps.

---

## Why these three, and not the other twenty

Selection is two-stage. **Significance is the gate**: only the 23 genes that passed the
haplotype analysis (BH q ≤ 0.05) were considered at all. **Among those 23**, the genes were
judged on **protein function and published literature alone** — their rank by q and their effect
sizes played no part in choosing between them. That separation matters, because
the two rankings disagree: the top gene by q (`6HG0619810`, starch) has **no functional
call from Swiss-Prot, InterPro or NCBI nr**, and the second (`3HG0301210`, fiber) carries
a **domesticated-transposase** domain. Neither can be argued for from function.

These three were the only genes whose annotation gives a route to the trait that a reader
can accept **without a supporting experiment**. The remaining 20 fall into: mobile-element
derived (2), defence/immunity (3), housekeeping or broad regulation (5), indirect
developmental or signalling routes to the trait (5), and unannotated (2) — plus 3 more
that are plausible only at family level. All are documented in `TRAIT_CANDIDACY.md`.

---

## `3HG0301300` — GPAT6, glycerol-3-phosphate 2-O-acyltransferase → **fiber**

`fiber/T1__07__Glycerol-3-phosphate_2-O-acyltransferase_6__3HG0301300.pdf`

### The literature

GPAT4/6/8 are a land-plant-specific clade of **bifunctional acyltransferase/phosphatases**.
Unlike every other eukaryotic membrane-bound GPAT — which acylates the *sn*-1 position to
make lysophosphatidic acid — this clade prefers ***sn*-2 and carries its own phosphatase
domain**, so the product is *sn*-2 **monoacylglycerol**. That 2-MAG is the committed
precursor of **cutin**. Tomato GPAT6 is the central enzyme of fruit cutin biosynthesis,
and GPAT6 knockdown in tomato and *N. benthamiana* measurably alters **cell-wall
properties of the leaf epidermis**.

**Why that reaches fiber:** cutin and suberin of the pericarp and testa are recovered in
the **insoluble dietary-fibre fraction**. A gene that sets how much cutin the grain coat
deposits changes measured fibre directly, not through a long chain of inference.

The annotation carries the mechanism on its face — both halves of the bifunctional enzyme
are visible: `GO:0016746` acyltransferase **and** `GO:0016791` phosphatase, plus
`GO:0010143` **cutin biosynthetic process** as an explicit term. Pfam PF01553
(acyltransferase) + PF23270 (GPAT-specific).

### Location and locus

- Gene **3H:545,000,419–545,004,865**, + strand, 4.4 kb
- `fiber_L08` — 3H:543,987,340–545,715,157, **1,727.8 kb**, 899 member SNPs; the widest
  locus in the study
- Lead SNP 3H:544084309, −log10p 6.2422 — **916 kb away**

> **The 916 kb is not the problem it looks like.** The nearest clump member sits **20 bp
> from the gene body**, with 47 members within ±50 kb. `fiber_L08` is a contiguous LD
> block — step 01 severs a locus at the first internal gap >50 kb — so the gene is
> densely covered by SNPs at r² ≥ 0.5 with the lead, not stranded in a gap.
> Distance-to-lead is the wrong statistic in a 1.7 Mb LD block; **SNP coverage is the
> right one**, and this gene passes it best of the three.

### The crosshap result

54 SNPs in window · 123 of 290 accessions grouped (167 unassigned) · 3 groups

| group | n | mean fibre BLUP | SD |
|---|---|---|---|
| **C** | 18 | **+0.281** | 0.240 |
| B | 20 | −0.081 | 0.177 |
| A | 85 | −0.120 | 0.234 |

KW H = 26.98, df = 2, p = 1.38e-06 → **q = 5.26e-06**, η² = 0.208, Δ = 1.51 SD

**Read honestly:** the best q of the three, and the largest top-vs-bottom gap — but the
*weakest pattern*. Groups B (−0.081) and A (−0.120) are effectively identical; the entire
signal is **one distinct haplotype C on 18 accessions** sitting above a flat background.
A single-haplotype contrast, not a dose-response. Only 42% of the panel is grouped.

---

## `5HG0487060` — Glucan endo-1,3-β-glucosidase, GH17 + X8/CBM43 → **fiber**

`fiber/T1__09__Glucan_endo-1_3-beta-glucosidase_5__5HG0487060.pdf`

### The literature

The single most direct functional link in the whole gene set. Barley's own
**(1→3,1→4)-β-D-glucan endohydrolases EI and EII** — encoded by *HvGlb1* and *HvGlb2*,
the enzymes that hydrolyse high-molecular-weight mixed-linkage glucan during germination —
are classified in **glycoside hydrolase family 17**. This gene is a GH17 enzyme carrying
the **X8 (CBM43) carbohydrate-binding module** typical of the family.

Mixed-linkage β-glucan is the **dominant soluble dietary fibre of the barley grain**. An
enzyme from the family that degrades it, sitting under a fibre peak, needs no supporting
narrative.

Annotation: Pfam **PF00332** (GH17) + **PF07983** (X8); `GO:0004553` hydrolase acting on
O-glycosyl bonds, `GO:0005975` carbohydrate metabolic process.

> **Caveat, stated plainly.** Step 09's sequence-based mapping puts *HvGlb1* on 1H
> (396.0 Mb) and *HvGlb2* on 7H (624.4 Mb). This gene is on **5H (462.7 Mb)**, so it is an
> **uncharacterised GH17 paralog — not either known endohydrolase**. GH17 also contains
> many (1,3)-β-glucanases acting on callose rather than mixed-linkage glucan. Right
> family, unproven member. This should be stated in the manuscript, not glossed.

### Location and locus

- Gene **5H:462,726,793–462,731,160**, + strand, 4.4 kb
- `fiber_L14` — 5H:462,555,842–462,848,780, **292.9 kb**, 49 member SNPs
- Lead SNP 5H:462564771, −log10p 6.0506 — 162 kb away; nearest member SNP 1.5 kb, 13 within ±50 kb

Note this is the **weakest fibre lead in the study** (6.0506 against a threshold of
6.0454 — it clears Bonferroni by 0.005). The gene sits in the outer half of its locus.

### The crosshap result

41 SNPs in window · 107 of 290 accessions grouped (183 unassigned) · 5 groups

| group | n | mean fibre BLUP | SD |
|---|---|---|---|
| **E** | 10 | **+0.229** | 0.428 |
| C | 16 | +0.195 | 0.196 |
| B | 31 | +0.007 | 0.284 |
| D | 14 | −0.094 | 0.149 |
| **A** | 36 | **−0.169** | 0.191 |

KW H = 26.63, df = 4, p = 2.36e-05 → **q = 7.47e-05**, η² = 0.222, Δ = 1.41 SD

**Read honestly:** the **best-behaved haplotype structure of the three** — a clean
monotonic gradient across all five groups, which is what an allelic series at a genuinely
functional locus looks like, and far more persuasive than one deviant group. It also has
the **largest η² (0.222)** of the three. Against that: the weakest lead SNP, the **lowest
assignment rate in the entire study** (only 107 of 290 grouped — 63% unassigned, against
a study median of ~70% *assigned*), and a top group of n = 10.

---

## `3HG0301710` — Plastidic anion/phosphate transporter, PHT4 family → **starch**

`starch/T1__11__Probable_anion_transporter_3_chloroplastic__3HG0301710.pdf`

### The literature

**Phosphate is the allosteric brake on ADP-glucose pyrophosphorylase (AGPase)**, the
committed and rate-limiting step of starch synthesis. The Pi:triose-phosphate balance
across the plastid envelope therefore sets the rate at which a plastid makes starch, and
the transporters that move Pi across that membrane are direct controllers of it. In
Arabidopsis, loss of the plastidic Pi transporter **PHT4;2** alters starch accumulation
by exactly this route — defective Pi export raises stromal Pi and inhibits starch
synthesis.

The best identification in the set: Swiss-Prot **Q8W0H5 = rice PHT4;3**, 84.2% identity
over 99.8% coverage; UniProt places it in the chloroplast membrane (transit peptide
res. 1–76), major facilitator superfamily / sodium-anion cotransporter TC 2.A.1.14. The
Morex V3 GFF independently projects it from **AT3G46980 = Arabidopsis PHT4;3**. Pfam
PF07690 (MFS); `GO:0005315` inorganic phosphate transmembrane transporter, `GO:0009536`
plastid, `GO:0055085` transmembrane transport.

> **Be precise about which paralog.** This gene is a **PHT4;3** ortholog. The documented
> starch phenotype is for **PHT4;2** (AT2G38060) — a *different* family member. The
> mechanism is therefore a **family-level inference**, not an established result for this
> paralog, and PHT4 substrate assignments are not uniform: PHT4;4 was later shown to be a
> chloroplast **ascorbate** transporter. The defensible claim is *"a plastid-localised
> anion/phosphate transporter at a starch locus"* — not *"the starch transporter"*.

This is also the one gene of the four originally hand-curated (`00_config/final_genes.tsv`:
Pho, **PHT**, AP2/ERF, BAHD) that survived step 03's `FLANK_BP = 0` tightening — so the
pipeline and the earlier manual curation independently agree here.

### Location and locus

- Gene **3H:546,473,709–546,477,564**, + strand, 3.9 kb
- `starch_L03` — 3H:546,310,483–546,479,480, **169.0 kb**, 39 member SNPs
- Lead SNP 3H:546433616, −log10p 6.3441 — **40 kb away, the closest of the three**;
  nearest member SNP 311 bp, 18 within ±50 kb
- The gene ends **1.9 kb inside the locus edge** — it sits right at the margin

### The crosshap result

16 SNPs in window · **205 of 290 accessions grouped** (85 unassigned) · 4 groups

| group | n | mean starch BLUP | SD |
|---|---|---|---|
| **A** | 137 | **+1.632** | 4.125 |
| D | 10 | −0.174 | 4.298 |
| B | 41 | −1.601 | 6.984 |
| **C** | 17 | **−5.545** | 6.522 |

KW H = 18.76, df = 3, p = 3.06e-04 → **q = 3.83e-04**, η² = 0.078, Δ = 1.31 SD

**Read honestly:** the **weakest effect** (η² = 0.078, about a third of the two fibre
genes) — within-group SDs of 4.1–7.0 dwarf the differences between group means, which is
why a 1.31 SD top-to-bottom gap still yields a modest η². Against that, it has by far the
**best data quality**: the only one of the three tested on a *majority* of the panel, the
closest to its lead SNP, the best protein identification, and a **monotone 4-group
gradient** in which the common haplotype A (n = 137) is the high-starch class.

---

## How the three compare

None is clean on every axis, and they fail in different places:

| strongest on | gene |
|---|---|
| statistical significance (q) | **GPAT6** — 5.26e-06 |
| effect size (η²) and haplotype pattern | **GH17** — 0.222, monotonic across 5 groups |
| data quality, proximity to lead, protein ID | **PHT4;3** — 205/290 grouped, 40 kb, 84.2% id |

| weakest on | gene |
|---|---|
| haplotype pattern — one deviant group, flat background | GPAT6 |
| assignment rate (107/290) and lead strength (6.05) | GH17 |
| effect size (η² = 0.078) | PHT4;3 |

**A general caution.** A Kruskal–Wallis p-value at a locus that was itself selected for
association with this very phenotype is close to guaranteed; **η² and the top-vs-bottom
gap are what a reader can actually judge**, and are reported above for every gene. The
23 significant genes sit in only 12 loci, so genes within a locus are in LD and are not
independent discoveries — report loci alongside genes.

---

## The 3H long arm, ~544–546 Mb

Two of the three live **2.5 Mb apart on the same chromosome arm**:

- `3HG0301300` **GPAT6** — a cell-wall/cutin gene — in `fiber_L08` (544.0–545.7 Mb)
- `3HG0301710` **PHT4;3** — a starch-synthesis gene — in `starch_L03` (546.3–546.5 Mb)

A fibre gene and a starch gene as near neighbours, each under its own trait's peak. **If a
fibre/starch trade-off is to be argued, this region is where to look** — not the 7H
~573.6 Mb locus that the retired `04_.../07_fiber_starch_tradeoff_direction/` was built
on. That step's premise no longer holds under run V4: **no gene is tested for both fiber
and starch**, and the three 7HG0729xxx genes it analysed are not in the current candidate
set at all (see its `BREAKING_CHANGES.md`).

Note also that within `starch_L03` the functionally believable gene (**PHT4;3**,
q = 3.8e-04) is **not** the top-q gene (`3HG0301720`, a CRLK1-like kinase, q = 9.0e-06) —
the same divergence between statistical rank and biological plausibility seen in
`fiber_L08`.

---

## Beta-glucan contributes no gene here — deliberately

No beta-glucan gene reached tier 1, and the reason is worth recording:

1. **Its three strongest loci are gene deserts.** `betaglucan_L02` (−log10p **7.60**,
   the strongest beta-glucan signal), `_L05` (7.33) and `_L04` (7.19) return **zero**
   annotated genes. Mean −log10p of the gene-empty beta-glucan loci is **6.72** against
   **6.48** for gene-containing ones — the genes come from the *weaker* half of the signal.
2. **No canonical beta-glucan gene is near any beta-glucan lead.** Recomputed against the
   current 36-locus set: the nearest is *HvGlb1* at **83.5 Mb** from `betaglucan_L01`.
   **HvCslF6**, the main barley mixed-linkage-glucan synthase, is on 7H — and there is
   **no beta-glucan locus on 7H at all**. *HvCslH1* is 434 Mb away.
3. **Every beta-glucan effect size is small** — η² 0.019–0.078, against 0.21–0.28 for the
   two fibre genes here.

This is a **reportable negative**, not a gap: beta-glucan has real architecture in this
panel (11 loci, three stronger than anything in the starch set) but its causal variants
sit in intergenic space, away from the annotated genes and from the entire canonical
CslF/CslH/Glb pathway — pointing to regulatory or structural variation rather than coding
change in a known enzyme.

---

## Table 2 of the manuscript (script, added 2026-09-27)

> **2026-09-30: Table 2 left the main text** (S. Hübner: "put all genes in one supmat"; user decision: no tables in the
> main text). Its three rows are now part of the all-candidate-genes Online Resource,
> `10_USED_Paper_writing/new_publishing_paper/supplementary/tables/ESM_candidate_genes.xlsx`, built by
> `…/supplementary/scripts/make_ESM_candidate_genes.py`, which reads `results/tables/Table_2_genes_carried_forward.tsv`
> for the "carried forward" column and checks the three rows against it cell by cell. The script below and its
> outputs are kept unchanged as that source; its `.md` is no longer pasted into the manuscript.

`scripts/01_make_table2_paper.R` writes **Table 2 of the TAG manuscript**: the three genes carried forward,
one row each, with chromosome, gene position, lead SNP, distance to it, accessions grouped (groups), BH *q*,
η², the difference in SD between the highest and lowest haplotype group, and the Swiss-Prot identity/coverage.

```bash
Rscript 06_USED_genes_selected_to_present/scripts/01_make_table2_paper.R     # seconds
```

| output (`results/tables/`) | contents |
|---|---|
| `Table_2_genes_carried_forward.md` | the pipe table, caption and footnotes, **pasted verbatim** into `10_USED_Paper_writing/new_publishing_paper/Results_Discussion_Conclusions.md` |
| `Table_2_genes_carried_forward.tsv` | raw values beside each formatted cell, for checking |

Inputs, read-only: step 04 `Stats/gene_results.tsv`, and step 05 `Table_significant_genes_annotated.tsv`
(`sp_pident`, `sp_qcovhsp`; these two columns are **not** in `Table_significant_genes_paper.tsv`, the Online
Resource 2 source). If step 04 or 05 is re-run, re-run the script and re-paste the `.md`. The GDSL esterase at
7H (Results Ch. 4) is **not** in the table; the user decided on 2026-09-27 that its numbers stay in the text.

---

## Filename convention

```
<trait>/T1__<q_rank>__<functional_name>__<gene_id>.pdf
```

`T1` marks the candidacy tier (**strong**) among the 23; `<q_rank>` is the gene's rank by
BH q among all 23, carried over unchanged from the source filename so every figure traces
back. The two numbers disagree on purpose — **tier is biology, rank is statistics**.

## Source and provenance

```
04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Significant_genes/by_trait/<trait>/
```

Run `loci_LDspan_eps06_V4` — MGmin = 2, ε = 0.6, both fixed; Kruskal–Wallis per gene on
BLUP phenotype ~ haplotype group (unassigned `hap 0` and missing phenotypes dropped);
Benjamini–Hochberg on raw p within each trait, α = 0.05.

**These are copies.** If step 04 is re-run, regenerate them from the source above rather
than editing anything here. The selection is recorded machine-readably in
`../05_USED_gene_annotation_analysis/08_USED_annotation_master/results/tables/trait_candidate_calls.tsv`
(`trait_candidate_strength == "strong"`).

### History

| date | change |
|---|---|
| 2026-09-10 | created with 8 genes — 3 tier-1 (strong) + 5 tier-2 (plausible) |
| 2026-09-10 | **reduced to the 3 tier-1 genes**; the 5 tier-2 copies removed, originals intact in step 04. Empty `betaglucan/` directory removed |
| 2026-09-22 | **Wording fix, no change to the selection.** "The crosshap statistics played no part in the selection" read as if significance did not matter. Corrected here, in `TRAIT_CANDIDACY.md` and in `10_USED_Paper_writing/CLAUDE.md` to state the actual logic: significance is the gate, and function/literature choose among the 23 significant genes |
| 2026-09-27 | **`scripts/01_make_table2_paper.R` added** (user approval, manuscript review item 4): it writes Table 2 of the manuscript to `results/tables/Table_2_genes_carried_forward.{md,tsv}`. It only reads steps 04 and 05. The selection is unchanged |
| 2026-09-30 | **Table 2 removed from the manuscript** (S. Hübner's comments, item S2 of the manuscript review): its rows moved to the all-candidate-genes Online Resource (`10_USED_Paper_writing/new_publishing_paper/supplementary/`), which reads `Table_2_genes_carried_forward.tsv`. Nothing in this folder changed. The selection is unchanged |

## Sources

- [Glycoside Hydrolase Family 17 — CAZypedia](https://www.cazypedia.org/index.php/Glycoside_Hydrolase_Family_17)
- [Barley grain (1,3;1,4)-β-glucan content: transcript and sequence variation in synthase and endohydrolase genes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6872655/)
- [Enzymes in 3D: synthesis, remodelling and hydrolysis of cell wall (1,3;1,4)-β-glucans](https://academic.oup.com/plphys/article/194/1/33/7245803)
- [A distinct type of glycerol-3-phosphate acyltransferase with sn-2 preference and phosphatase activity producing 2-monoacylglycerol (PNAS)](https://www.pnas.org/doi/10.1073/pnas.0914149107)
- [Glycerol-3-phosphate acyltransferase GPAT6 from tomato plays a central role in fruit cutin biosynthesis](https://academic.oup.com/plphys/article/171/2/894/6115299)
- [GPAT6 controls filamentous pathogen interactions and cell wall properties of the leaf epidermis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6767537/)
- [Land-plant-specific GPAT family: substrate specificity, sn-2 preference, and evolution](https://academic.oup.com/plphys/article/160/2/638/6109449)
- [The sink-specific plastidic phosphate transporter PHT4;2 influences starch accumulation and leaf size in Arabidopsis](https://academic.oup.com/plphys/article/157/4/1765/6108878)
- [UniProt Q8W0H5 — Probable anion transporter 3, chloroplastic (PHT4;3), *Oryza sativa*](https://www.uniprot.org/uniprotkb/Q8W0H5/entry)
