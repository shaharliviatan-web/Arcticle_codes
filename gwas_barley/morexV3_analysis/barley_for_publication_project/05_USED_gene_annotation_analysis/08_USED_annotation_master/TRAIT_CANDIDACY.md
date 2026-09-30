# Trait candidacy — which of the 23 genes plausibly affect their trait

Literature-based judgement of the **23 haplotype-significant genes** from step-04 run
`loci_LDspan_eps06_V4`, annotated in `05_USED_gene_annotation_analysis/`.

Written 2026-09-10.

**Machine-readable version → [`results/tables/trait_candidate_calls.tsv`](results/tables/trait_candidate_calls.tsv)**
(`gene_id, trait, trait_candidate_strength, candidate_rationale`). `scripts/01_build_paper_table.R`
merges it into the paper tables automatically on the next run and will not overwrite it.

---

## What this document is, and what it is not

This is the **manual biological judgement** the annotation chain deliberately leaves
open. See `README.md` § *Identification is not biological candidacy*: knowing what a
protein **is** (Swiss-Prot / InterPro / nr) and believing it **affects the trait** are
independent questions, and the pipeline only answers the first.

Three deliberate constraints on what follows:

1. **Significance is the gate; function decides among the survivors.** Only the 23
   haplotype-significant genes (BH q ≤ 0.05) were eligible at all. *Among those 23*, genes
   were judged on protein function alone — the crosshap q-values, η² and effect sizes played
   no part in choosing between them. A gene is not more believable here because it ranked first.
2. **Believability, not proof.** The bar is "a reader can accept this connection from
   the gene's function without a supporting experiment" — not causality.
3. **No trait is favoured.** Protein contributes no genes (nothing testable in V4), so
   only beta-glucan, fiber and starch appear.

Where the two judgements disagree is itself informative, and is noted below.

---

## Tier 1 — genuinely believable (3 genes)

### `5HG0487060` — Glucan endo-1,3-β-glucosidase, GH17 + X8/CBM43 → **fiber**

The best of the set. Barley's own (1,3;1,4)-β-glucan endohydrolases **EI and EII
(HvGlb1, HvGlb2) are GH17 family enzymes** — this is precisely the enzyme family that
degrades barley mixed-linkage glucan, the dominant soluble fibre of the grain. A
cell-wall glucan hydrolase sitting under a fibre peak requires no supporting story.

> **Caveat, stated honestly.** Step 09's sequence-based mapping places HvGlb1 on 1H
> (396.0 Mb) and HvGlb2 on 7H (624.4 Mb). This gene is on 5H (462.6 Mb), so it is an
> **uncharacterised GH17 paralog, not either known endohydrolase**. Many GH17s are
> callose (1,3)-β-glucanases rather than (1,3;1,4). Right family, unproven member.

Evidence: Swiss-Prot 62.7% id / 90.9% cov; Pfam PF00332 (GH17) + PF07983 (X8);
GO:0004553 hydrolase acting on O-glycosyl bonds, GO:0005975 carbohydrate metabolism.

### `3HG0301710` — Plastidic phosphate transporter, PHT4 family → **starch**

Mechanistically direct. **Pi is the allosteric brake on AGPase**, so the rate at which
phosphate crosses the plastid envelope sets the rate of starch synthesis. In
Arabidopsis, loss of the plastidic Pi transporter `PHT4;2` alters starch accumulation by
exactly this route (defective Pi export → excess stromal Pi → inhibited starch
synthesis). A plastid Pi transporter under a starch peak is about as clean as candidate
genes get.

> **Corrected 2026-09-10 — be precise about which paralog.** The Swiss-Prot hit is
> **Q8W0H5 = rice PHT4;3** (84.2% id, 99.8% cov), and the Morex V3 GFF independently
> projects the gene from **AT3G46980 = Arabidopsis PHT4;3**. The documented starch
> phenotype above is for **PHT4;2** (AT2G38060), a *different* family member, so the
> mechanism is a **family-level inference, not an established result for this paralog** —
> and PHT4 substrate assignments are not uniform (PHT4;4 was later shown to be a
> chloroplast *ascorbate* transporter). The defensible claim is "a plastid-localised
> anion/phosphate transporter at a starch locus", not "the starch transporter".

This is also the **one gene of the four originally curated** (`00_config/final_genes.tsv`:
Pho, PHT, AP2/ERF, BAHD) that survived step 03's `FLANK_BP = 0` tightening — so the
pipeline and the earlier hand-curation agree here.

Evidence: Swiss-Prot 84.2% id / 99.8% cov; Pfam PF07690 (MFS); GO:0005315 inorganic
phosphate transmembrane transporter, GO:0009536 plastid, GO:0055085 transmembrane
transport.

### `3HG0301300` — GPAT6, glycerol-3-phosphate 2-O-acyltransferase → **fiber**

A **bifunctional acyltransferase/phosphatase producing sn-2 monoacylglycerol** — the
committed step of **cutin** biosynthesis, and its own GO annotation says so
(`GO:0010143 cutin biosynthetic process`). Cutin and suberin of the pericarp and testa
are recovered in the insoluble-fibre fraction, and tomato GPAT6 is documented to alter
epidermal cell-wall properties directly. Short, direct link.

Evidence: Swiss-Prot 69.3% id / 97.8% cov; Pfam PF01553 (acyltransferase) + PF23270;
GO:0010143, GO:0016746 acyltransferase, GO:0016791 phosphatase — the dual activity is
visible in the annotation itself.

---

## Tier 2 — plausible, but the argument needs hand-waving (5 genes)

| gene | trait | the route, and where it is weak |
|---|---|---|
| `3HG0301170` Myb/SANT TF | fiber | **R2R3-MYBs are the master regulators of secondary cell wall** biosynthesis (cellulose, xylan, lignin) — a strong family prior for a fibre locus. Weak on specificity: only **33% query coverage**, so this is a domain-level call, and MYB is a very large family with many unrelated roles |
| `3HG0301330` WOX9 homeobox | fiber | Developmental, not a cell-wall gene. Route is indirect: grain size and shape shift the **pericarp:endosperm ratio**, moving fibre as a *percentage* of grain without touching fibre synthesis. Identification is solid (82.9% id) |
| `4HG0341560` PEBP / FT-like | betaglucan | **Phenological, not biochemical**: flowering date sets the temperature and water regime during grain filling, and warm dry filling *raises* barley grain β-glucan while cool wet filling lowers it. A credible pleiotropic hit in a wild panel spanning a climate cline. Note the Swiss-Prot name "VERNALIZATION 3" only means *nearest curated PEBP* — true VRN3/HvFT1 is on **7H, not 4H** |
| `1HG0077790` Ca²⁺-ATPase (ACA) | betaglucan | CslF/CslH synthesise β-glucan in the Golgi and glycosyltransferases are divalent-cation dependent, so secretory-compartment ion homeostasis is a real mechanism — but an entirely **generic** one. Nothing specific to β-glucan |
| `5HG0487090` carboxylesterase, α/β hydrolase | fiber | The fold family does contain cell-wall esterases (cutinases, pectin acetylesterases, GDSL lipases). But the actual hit is a **tulip-specific** tuliposide-converting enzyme at 46% identity — effectively "some esterase" — and it sits in the **same locus** (`fiber_L14`) as the GH17 above, so LD hitch-hiking is the simpler explanation |

---

## Tier 3 — would not argue for these (13 genes)

**Mobile-element derived — no functional route to the trait**
- `3HG0301210` — ZMYM1/FAM200 RNase-like domain, a **domesticated transposase** domain. Flagged because it ranks **2nd of 23 by q**.
- `3HG0301260` — Retrovirus-related Pol polyprotein, TNT 1-94 **retroelement**.

**Defence / immunity — generic**
- `3HG0301110` HSPRO1 (nematode-resistance-like) · `3HG0301390` C2 domain, SRC2/BAP · `7HG0703630` L-type lectin receptor kinase IX.1 (LecRKs do touch cell-wall-integrity sensing, but generically, at 46.9% id).

**Housekeeping / broad regulation**
- `3HG0301400` PPR-DYW (organellar RNA editing) · `4HG0395770` SIL1/FES1 (ER BiP cofactor) · `4HG0395810` MA3 translation regulatory factor · `3HG0287270` RING-H2 ATL E3 ligase · `3HG0287240` LEA_2 / Lea14-like.

**Starch, all too indirect**
- `6HG0619820` **Kinesin KIN-7D** — a confident identification (82.8% id, full coverage) of a microtubule motor. The cleanest example in this set of *identification without candidacy*; the README already names it as such.
- `3HG0301720` CRLK1-like Ca²⁺/CaM cold-signalling kinase · `1HG0049610` SH3P2-like (endocytosis / autophagosome formation).

## Cannot be judged (2 genes)

`6HG0619810` (starch) and `4HG0339140` (fiber) returned **no call from any of the three
sources**, and `4HG0339140` has no Pfam or GO either. `6HG0619810` **ranks 1st of 23 by
q** — the strongest statistical signal in the study is a protein nobody can name.

---

## Two structural observations that matter for the manuscript

### 1. `fiber_L08` is a choice, not a collection

**7 of the 12 fiber genes sit in one locus** — `fiber_L08`, 3H, a **1.73 Mb** span, the
widest in the study. They are in one LD block, so **at most one is causal**: you are
choosing among them, not accumulating evidence.

| gene in `fiber_L08` | tier |
|---|---|
| `3HG0301300` GPAT6 | **strong** ← the one to back |
| `3HG0301170` Myb TF | plausible ← second choice |
| `3HG0301330` WOX9 | plausible |
| `3HG0301210` ZMYM1 | unlikely (transposase-derived) — **rank 2 by q** |
| `3HG0301260` TNT1-94 Pol | unlikely (retroelement) |
| `3HG0301390` C2/SRC2 | unlikely |
| `3HG0301400` PPR | unlikely |

The two top-ranked by q in this locus are the two mobile-element-derived ones. This is
the clearest case in the study where statistical rank and biological plausibility point
in different directions, and it is an argument for reporting **loci rather than gene
counts** — which step 04's `locus_summary.tsv` and its counting caveat already do.

### 2. The 3H long arm ~544–546 Mb carries both signals

`fiber_L08` (544.0–545.7 Mb) and `starch_L03` (546.3–546.5 Mb) are **neighbours**, and
the two strongest candidates of the whole set live there side by side:

- `3HG0301300` **GPAT6** — cutin, a cell-wall/fibre gene — in `fiber_L08`
- `3HG0301710` **PHT4 Pi transporter** — starch synthesis — in `starch_L03`

A cell-wall gene and a starch gene adjacent on the same chromosome arm. **If a
fiber/starch trade-off narrative is wanted, this region is where to look** — not the 7H
~573.6 Mb locus that the retired `04_.../07_fiber_starch_tradeoff_direction/` was built
on. That step's premise no longer holds under run V4: no gene is tested for both fiber
and starch, and the three 7HG0729xxx genes it analysed are not in the current candidate
set at all (see its `BREAKING_CHANGES.md`).

Note also that within `starch_L03` the functionally believable gene (**PHT4**,
q = 3.8e-04) is **not** the top-q gene (`3HG0301720` CRLK1, q = 9.0e-06) — the same
rank-vs-plausibility divergence as in `fiber_L08`.

---

## Sources

- [Glycoside Hydrolase Family 17 — CAZypedia](https://www.cazypedia.org/index.php/Glycoside_Hydrolase_Family_17)
- [Barley grain (1,3;1,4)-β-glucan content: transcript and sequence variation in synthase and endohydrolase genes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6872655/)
- [Enzymes in 3D: synthesis, remodelling and hydrolysis of cell wall (1,3;1,4)-β-glucans](https://academic.oup.com/plphys/article/194/1/33/7245803)
- [The sink-specific plastidic phosphate transporter PHT4;2 influences starch accumulation and leaf size in Arabidopsis](https://academic.oup.com/plphys/article/157/4/1765/6108878)
- [A distinct type of glycerol-3-phosphate acyltransferase with sn-2 preference and phosphatase activity producing 2-monoacylglycerol (PNAS)](https://www.pnas.org/doi/10.1073/pnas.0914149107)
- [GPAT6 controls filamentous pathogen interactions and cell wall properties of the leaf epidermis](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6767537/)
- [Glycerol-3-phosphate acyltransferase GPAT6 from tomato plays a central role in fruit cutin biosynthesis](https://academic.oup.com/plphys/article/171/2/894/6115299)
- [Molecular and functional characterization of PEBP genes in barley](https://academic.oup.com/plphys/article/149/3/1341/6107773)
- [Barley and wheat beta-glucan content influenced by weather, fertilization, and genotype](https://www.frontiersin.org/journals/sustainable-food-systems/articles/10.3389/fsufs.2023.1326716/full)
- [Changes caused by genotype and environmental conditions in beta-glucan content of spring barley](https://pubmed.ncbi.nlm.nih.gov/18551369/)

---

## Provenance

| | |
|---|---|
| genes assessed | 23 (beta-glucan 6, fiber 12, starch 5) |
| source table | `results/tables/Table_significant_genes_paper.tsv` (built 2026-09-10) |
| step-04 run | `loci_LDspan_eps06_V4` |
| tiers | strong 3 · plausible 5 · unlikely 13 · no_annotation 2 |
| basis | only the 23 haplotype-significant genes eligible (**significance is the gate**); among them, protein function and published literature only — **q-rank and effect sizes not used to choose between them** |
