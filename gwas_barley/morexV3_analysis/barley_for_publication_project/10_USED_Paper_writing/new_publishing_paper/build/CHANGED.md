# What changed since the mini paper

Companion to **`BLUEPRINT.md`**. The blueprint says what the Results will report; this file says where
that differs from the earlier BSc mini paper, so nothing is carried over by habit.

The earlier write-up used **MorexV2**; this study is a re-call and re-analysis against **MorexV3**,
and several steps were rebuilt. Some results changed, a few reversed.

Status: covers Results **Chapter 2** (GWAS, loci, candidate genes), **Chapter 3** (haplotypes,
genes, elite lines) and the **Discussion**. Written 2026-09-17, Discussion added 2026-09-24.

---

# Chapter 2 — GWAS, loci and candidate genes

## Statements that are now wrong and must not be reused

| the mini paper said | now |
|---|---|
| "Highly heritable starch yielded no significant signals" | Starch maps: it has genome-wide-significant SNPs in several loci. The argument that mapping power does not follow heritability **survives**, but this example must be rewritten |
| β-glucan is the most mappable trait | **Fiber** is, by a wide margin. β-glucan is second |
| 24 SNPs — 15 significant plus 9 "marginal" — in 20 QTLs | About fifty genome-wide-significant SNPs in **36 loci**. **The "marginal" class no longer exists**; every locus reported is genome-wide significant |
| Protein yielded a single locus | **Two** loci, both on 3H. The second only just clears the threshold |
| Bonferroni threshold −log₁₀*P* = 6.7712 | The threshold moved down after the LD pruning was redone on a physical rather than a marker-count window |
| A shared 7H window containing nine genes, three of them significant for both starch and fiber | Those genes are **not** in the current candidate set. The co-localization itself is now *stronger* — both leads are genome-wide significant — but the gene-level evidence behind the old carbon-allocation argument is gone. See the open question below |
| 108 candidate genes at 20 loci; β-glucan contributed the most (72) | Roughly **55** genes at fewer than half the loci. The per-trait ordering **inverted**: β-glucan went from the most genes to the fewest, fiber from few to most |
| The candidate-gene window was ±200 kb around each significant SNP, from the LD-decay distance | The window is the **LD locus span itself, with no flank** |

---

## Method changes behind those results

* **LD pruning** for the PCA and the significance threshold now uses a physical (1 Mb) window
  instead of a fixed marker count, which at this SNP density spanned a wildly variable distance.
  The threshold moved with it.
* **Locus definition** changed from single-linkage at a fixed distance to **iterative LD clumping**,
  with each locus severed at the first large internal gap. A locus is now a real LD block, and
  carries a span and a member count that did not exist before.
* **Sub-threshold peaks were dropped by decision.** Nothing below the genome-wide threshold is
  carried into the gene search.
* **The gene search no longer depends on LD decay.** LD decay is descriptive context only, which is
  why its panel leaves the figure.

---

## Individual signals that moved

* **Four of the earlier 24 SNPs are no longer significant** — three β-glucan (on 2H, 1H and 7H) and
  one starch (7H). The other twenty survived, and **every surviving "marginal" SNP is now
  genome-wide significant**.
* **The 4H β-glucan cluster grew and split.** What was described as four significant SNPs in one
  cluster is now ten SNPs over a wider interval, resolving into **two separate loci**.
* **The starch locus set changed composition**: one earlier peak (7H) fell below threshold, and two
  loci not previously reported (1H, 3H) came in.

---

## New in this manuscript

* **A 1H region where fiber and starch loci interleave** — a starch locus nested between two fiber
  loci, with a third fiber locus just upstream. Not reported before, and a second independent
  instance of the starch–fiber co-localization.
* **β-glucan shares no genomic region with any other trait**, despite correlating with both fiber
  and starch phenotypically. Its nearest cross-trait locus is megabases away.
* **Most signals rest on rare alleles** — carried by only a dozen or two accessions — while the only
  two common-allele loci in the study are both starch.
* **Genomic inflation is reported as actual values** (λ_GC close to 1 for all four traits) rather
  than as "λ > 0.9".

---

## Figures

* **Manhattan plots** now paint each locus over the grey background, showing the physical extent of
  each association rather than only its peak.
* **Final Fig. 3 (approved 2026-09-22):** Manhattan + QQ only, one row per trait; no locus key or
  labels (locus identities go to the Online Resource table); QQ carries λGC only.
* **The PC scree panel is dropped** — the correction is the one already used for this collection by
  Potapenko et al. (2026) and is cited rather than re-derived.
* **The LD-decay panel is dropped** — it no longer defines the candidate-gene window.
* Old Fig. 3 and Fig. 4 are obsolete in content. Fig. 1 and Fig. 2 (chapter 1) are unchanged and
  stand as they are.

---

## Open question carried forward

**The 7H story may not be dead.** Both the fiber and the starch lead SNP lie under a kilobase from
the *same* annotated gene. The current zero-flank rule excludes it, because both loci are only a few
hundred bases wide. A rule that admits a gene this close to a lead SNP regardless of locus span
would bring it back — and the downstream haplotype analysis would then have to be re-run for it.
This is the one place where the earlier manuscript's central argument might be recoverable in a new
form, and it should be settled before chapter 2 is written.

---
---

# Chapter 3 — Haplotypes, gene resolution and elite lines

## Statements that are now wrong and must not be reused

| the mini paper said | now |
|---|---|
| 45 of 108 candidate genes retained at FDR < 0.05 | 55 genes go in, **30 are testable**, **23 are significant**, and they sit in only **12 loci** |
| Five biologically plausible candidates; four presented as panels | **Three** carried forward — two fiber, one starch |
| *Pho* (α-glucan phosphorylase) and *PHT* for starch, AP2/ERF for β-glucan, BAHD for fiber | **Only *PHT* survives.** *Pho*, AP2/ERF and BAHD are **not in the current candidate set at all** and must not be mentioned as candidates |
| The β-glucan candidate is an AP2/ERF transcription factor regulating *HvCslF6* | **β-glucan contributes no presented gene.** This is a deliberate, reportable negative, not a gap |
| Three genes in the shared 7H window passed FDR for **both** starch and fiber, with reciprocal haplotype direction — the locus-level support for the trade-off | **No gene is tested for both traits.** Those three genes are not in the current candidate set. **This evidence cannot be reused in any form** |
| *PHT* showed "the strongest haplotype separation among the genes examined" (P = 1.7e-10, ε² = 0.23) | *PHT4;3* has the **weakest effect of the three** presented genes. Its strength is data quality, not effect size |
| Protein yielded no significant candidate genes | Protein has candidate genes but **none is testable** — a different and more specific statement |

## Method changes behind those results

* **The per-gene epsilon sweep is gone.** One grouping per gene at fixed parameters, one test, one
  correction. Sweeping and taking the smallest p was a forking path.
* **The double correction is gone.** The old version fed FWER-adjusted p-values into BH, which does
  not control FDR at any interpretable rate. BH now receives the raw p.
* **Gene selection is still two-stage, as before** — significance first, biological relevance
  second. What changed is that the biological judgement is now **recorded explicitly per gene**, with
  its rationale, separately from the annotation, and that rank by significance is not used to choose
  among the significant genes.
* **Annotation is now a three-source chain** with explicit provenance, and every source's own call is
  reported. The old HIGH/MEDIUM/LOW confidence tier was removed for conflating *what the protein is*
  with *whether it plausibly affects the trait*.
* **Effect sizes accompany every p-value** — η² and the top-vs-bottom gap in phenotype SD — because
  a p-value at a locus selected for that very phenotype is close to guaranteed.

## New in this manuscript

* **A comparison with elite cultivars** — five spring malting lines at the three presented genes.
  Not in the earlier version at all.
* **The divergence between statistical rank and biological plausibility.** The top gene by q has no
  functional call from any source, and the second carries a transposase-derived domain. Neither can
  be argued for from function.
* **An explicit counting caveat.** The 23 significant genes sit in 12 loci, and a single locus holds
  7 of the 12 fiber genes — they are one LD block, so at most one can be causal. Loci are reported
  alongside genes.
* **Honest per-gene weaknesses.** Each presented gene is reported with the axis on which it is
  weakest, rather than only its best statistic.

## Figure

* **Fig. 4 is rebuilt, not edited** (approved 2026-09-23, `08_USED_creating_figures/make_figure_4.R`).
  The old Fig. 4 showed four mini-paper genes — Pho, PHT, AP2/ERF, BAHD — of which only PHT survives
  in the current candidate set, so its content is obsolete and its outputs were deleted. The new
  figure shows the three genes carried forward, and adds the **elite-cultivar barcodes**, which did
  not exist in the earlier version at all. The old panels' Holm post-hoc brackets are kept; the
  omnibus statistic beside each panel is now the **BH q with η² next to it**, not a bare
  Kruskal-Wallis *P*.

---
---

# Discussion

The mini paper's Discussion has three sections. **Two survive; the third keeps its claim and loses
every example.** Three further sections are new.

## Section by section

| mini-paper section | status |
|---|---|
| **1. Ecological differentiation of grain composition** | **survives intact.** Ch. 1 is unchanged, so the argument, the citations and the physiological explanation all carry over |
| **2. A shared carbon-allocation axis** | **survives and is much stronger.** The mini paper could argue the trade-off from phenotype and one shared 7H locus. It now also appears in the direction of allelic effects at every locus, in three co-localized regions, and at a single haplotype carrying both effects |
| **3. Heritability and GWAS reveal different aspects of genetic architecture** | **claim survives, every example replaced** — see below |
| — | **new: what this design can and cannot resolve** (limitations of locus definition and haplotype testing) |
| — | **new: β-glucan — strong signals, no genes** |
| — | **new: wild alleles for cultivated barley**, expanded from the mini paper's Conclusions |

## What section 3 loses

| the mini paper said | now |
|---|---|
| "Despite exhibiting the highest broad-sense heritability, starch produced a sparse mapping profile" and yielded no significant signals | Starch maps: five loci. The claim that mapping power does not follow heritability **survives**, but on a different and cleaner demonstration — three traits of near-identical H² yield 2, 11 and 18 loci |
| β-glucan "yielded a robust, distributed signal comprising 11 loci across all seven chromosomes" | Still 11 loci, but **not across all seven chromosomes**: there is none on 7H, and that absence is now a finding in its own right |
| "the starch locus on chromosome 3H illustrates how marginal signal in the GWAS analysis produced a strong separation between haplotypes" | **There is no marginal class any more.** The example cannot be used in any form |
| The closing sentence: four genes identifying "direct starch metabolism, phosphate associated grain filling, glucan transcriptional regulation, and cell wall fiber structural modification" | **Three of those four genes no longer exist as candidates.** The replacement routes are cutin deposition (GPAT6), glucan hydrolysis (GH17), plastid phosphate transport (PHT4;3) and cell-wall ester modification (the 7H GDSL esterase) |

## Mini-paper content that carries over unchanged

Worth listing, because it is what keeps the Discussion in Sariel's voice rather than rebuilt from
nothing: the common-garden/BLUP argument that regional patterns are heritable rather than plastic;
Savin & Nicolas on heat and drought shortening grain filling; Anker-Nilssen on β-glucan's
environmental lability; the sucrose-competition reading of the carbohydrate axis with Lim et al. as
functional precedent; the whole-plant allocation framing; the protein C/N-balance argument with
Corke, Bogard and Simmonds; the practical conclusion that raising fiber need not cost protein; and
the linkage-drag argument with Hübner & Kantar.

## New constraints the Discussion must respect

* **No wild-vs-cultivated nutritional comparison is possible from our data** — the cultivated checks
  carry no NIR values (`section_materials/04_discussion/NOTE_cultivated_checks_no_NIR.md`). The mini
  paper's Introduction claim (Friedman & Atsmon 1988) stays a citation.
* **The five elite cultivars are one haplotype, not five observations**, so the absence of the
  high-fiber haplotype from them is consistent with selection but cannot demonstrate it.
* **The rare-allele result carries a correction caveat**: PC and kinship correction is weakest for
  rare alleles concentrated in a few related accessions, and λ_GC does not detect that failure mode.
