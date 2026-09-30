# Discussion — candidate list and the user's decisions

Built 2026-09-25; decisions recorded 2026-09-26. Every idea from the mini paper, the blueprint
(structural blueprint + hand-over items), the Results chapters 1–4, the project directories, the
dir-11 carrier-origin check, and new points from the session.

Sources: **MP** = mini paper · **BP** = BLUEPRINT.md · **R1–R4** = Results chapters · **D11** =
`11_USED_check_low_maf_geographic_distribution_minor_allele` · **NEW** = new in the 2026-09-25 session.

## Rules for drafting (user, 2026-09-26)

1. **Interpretation only.** Every kept item is written as interpretation and discussion, short and to
   the point. Results are not repeated unless a number carries the argument. Where the Results
   already made the point, the Discussion builds on it without restating it.
2. **dir-11 Word comment.** Every inference that rests on the dir-11 carrier-origin check gets a TODO
   Word comment (id 4xx at drafting), with this text:
   `[TODO: this inference rests on the dir-11 carrier-origin check (11_USED_check_low_maf_geographic_distribution_minor_allele), which is not in the Results. Decide whether to add its results (Online Resource + one M&M sentence) or keep the Discussion statement without them]`
   Applies to **2.12, 3.5** (4.8 dropped 2026-09-27).
3. Parameters (ε, MGmin, the 7H grouping choice) are reported in M&M only.

Status: **K** keep · **D** drop · **S** settled earlier · **?** open.

---

## §1 Ecological differentiation of grain composition (MP §1)

| # | candidate | source | status |
|---|---|---|---|
| 1.1 | Extends the known ecological structuring of this collection (Hübner et al. 2009, 2013) to grain composition | MP, BP | K |
| 1.2 | Common-garden BLUPs → regional differences are heritable, not plastic | MP, BP | K |
| 1.3 | North starch-rich vs Desert fiber/β-glucan-rich along soil texture and rainfall | MP, R1 | K — interpretation only, no restating of R1 |
| 1.4 | Heat and drought during grain filling cut starch deposition (Savin and Nicolas 1996) | MP, BP | K |
| 1.5 | Temperature correlates with no trait; weight on water availability | NEW | D |
| 1.6 | Heritable drought-escape syndrome | NEW | D |
| 1.7 | β-glucan's G×E fits its environmental lability (Anker-Nilssen et al. 2008) | MP, BP | K |
| 1.8 | Closing: the traits follow the environment | MP | K |
| 1.9 | Potapenko et al. 2026 *Mol Ecol* (same collection, local adaptation) | NEW | D |

## §2 A carbon-allocation axis (MP §2)

| # | candidate | source | status |
|---|---|---|---|
| 2.1–2.4 | Protein axis vs carbohydrate axis (Corke et al. 1989); sucrose competition; Lim et al. 2020; whole-plant strategy | MP, BP | K |
| 2.5 | Domestication shifted cultivated barley to the starch-rich end (Friedman and Atsmon 1988) + Word comment: own data? | MP | S — stays in §2 |
| 2.6–2.7 | Protein by C/N balance (Simmonds 1995; Bogard et al. 2010; Corke et al. 1989); raising fiber need not cost protein | MP, BP | K |
| 2.8 | Protein contributes no gene + Word comment: check the protein quantification method | BP item 5 | S — comment only |
| 2.9–2.10 | Three shared fiber/starch regions, three levels of resolution (7H one gene; 3H two genes; 1H none) | BP item 13 | K |
| 2.11 | One haplotype, two traits at 7H | BP item 12 | K |
| 2.12 | The 7H locus is the one shared signal geography does not explain (carriers from 10–11 coastal, desert and transitional sites, none northern; no desert excess) | D11 | K + **dir-11 Word comment** |
| 2.13 | Direction of allelic effects as evidence for the axis | BP §2, R2 | K — as evidence here, phrased "consistent with"; qualified in §3 (3.5) |
| 2.14 | Fiber tracks grain size; structural axis | NEW | D |
| 2.15 | Fiber genes are grain-coat genes | NEW | D |
| 2.16 | The 7H starch link is compositional, not mechanistic | BP item 11 | K |
| 2.17 | Closure caveat (traits are % of grain) | NEW | D |

## §3 Genetic architecture (MP §3: claim kept, examples replaced)

| # | candidate | source | status |
|---|---|---|---|
| 3.1–3.2 | Heritability did not predict mappability; mapping reflects allele frequency, effect size, LD | MP, BP item 1 | K |
| 3.3 | Starch's common-allele loci → polygenic | NEW | D |
| 3.4 | Most signals are rare alleles | BP item 2 | K |
| 3.5 | Minor-allele direction + carrier origin, with the user's four arguments and their limits (local adaptation; several Desert sites but a cohesive Desert lineage; within-site contrast; function and haplotypes) | BP item 2, D11, user | K + **dir-11 Word comment** |
| 3.6 | Wording rule: "minor", never "derived" | NEW | K (rule) |
| 3.7 | Rank vs plausibility; two-stage selection | BP item 6 | K |
| 3.8 | Closing sentence: four genes, four routes into the grain (cutin deposition, glucan hydrolysis, plastid phosphate transport, cell-wall ester modification), for future investigation — the counterpart of MP §3's last sentence | BP, MP | K — **in Sariel's wording, tone and brevity** (MP `06_discussion.md`, last sentence), only the four routes replaced |
| 3.9 | One short mechanism per gene, MP style (incl. Ma et al. 2021 for PHT) | MP, BP item 10 | K |
| 3.10 | 7H: regulatory rather than coding, with caveats | BP item 9 | K |
| 3.11 | GH17 maps to fiber, not β-glucan | NEW | D |
| 3.12 | PEBP/FT-like β-glucan gene → phenological route | NEW | D (no literature for the link; direction contradicts our own flowering correlation) |

## §4 What this design can and cannot resolve (new)

| # | candidate | source | status |
|---|---|---|---|
| 4.1–4.3 | Locus, not gene, is the unit of discovery; gene list depends on the locus definition; GDSL as the worked example | BP | K |
| 4.4 | Locus span vs allele frequency | NEW | D |
| 4.5 | Reference bias (pan-genome) | NEW | D |
| 4.6 | Haplotype testing needs SNP variation: an untested gene is not a gene without effect (protein's missing gene is a testing gap); one fixed setting cannot suit genes of different SNP density | BP | K — **merged into 4.9** |
| 4.7 | Low assignment rates | NEW | D |
| 4.8 | At some loci (incl. GH17, PHT4;3) whole sites carry the minor allele: allele and population cannot be separated | D11, user | D — drafted, then removed by the user 2026-09-27 (not mentioned) |
| 4.9 | **One paragraph (user, 2026-09-26), merging 4.2, 4.3, 4.6:** both steps that narrow the signal to genes can miss good candidates. **Locus step:** each locus is delimited by LD with its lead SNP and severed at large gaps between SNPs, so its span reflects these thresholds of the clumping procedure (values in M&M: r² 0.5, 50 kb gap), not only the physical extent of the causal variant; a gene just outside a narrow locus is never considered (the 7H GDSL gene was found only by examining the region directly). **Haplotype step:** 4.6. **Conclusion:** the presented genes are those that passed both filters, not a complete list. Say it smoothly, no defence of the method | BP, user | K |
| 4.10 | Fiber calibration (cultivated package only, no wild-barley lab reference) | NEW | **NOT in the text.** Reminder only: a TODO Word comment for the user, no sentence or paragraph |
| 4.11 | The 7H grouping choice | BP | D — M&M only |

## §5 β-glucan: strong signals, no genes (new)

| # | candidate | source | status |
|---|---|---|---|
| 5.1–5.3 | Strongest signals, no gene; small effects; no region shared with starch or fiber | BP items 3, 7 | K — **one opening sentence only**, no numbers |
| 5.3b | **Moved from Results Ch. 3 (user, 2026-09-26). Write the HvGlb1 distance, and restate the 7H / HvCslF6 absence briefly (it is also in Ch. 2) (user, 2026-09-26):** β-glucan loci do not lie near the canonical mixed-linkage glucan pathway; of the ten canonical synthases and endohydrolases the nearest to any β-glucan lead is *HvGlb1*, 83.5 Mb away; *HvCslF6* and *HvGlb2* lie on 7H, which carries no β-glucan locus. Src: `05_.../09_USED_canonical_betaglucan_gene_check/results/tables/canonical_gene_distances_current_loci.tsv`. M&M must state that the canonical genes were located on MorexV3 by protein sequence | R3 (moved) | K |
| 5.3c | Its strongest loci contain no annotated gene, so the genes tested came from its weaker loci (cut from Ch. 3 2026-09-26) | R2, R3 (cut) | D (§4 covers gene-empty loci) |
| 5.4 | "Not the §4 limitation" | BP | D |
| 5.5–5.6 | *HvCslF6* conserved (Garcia-Gimenez et al. 2019); common outcome (Geng et al. 2021; Houston et al. 2014) vs detection (Shu and Rasmussen 2014) | BP item 4 | K |
| 5.7 | Measure *HvCslF6* variation in our panel | NOTE | D |
| 5.8 | Rule: no "regulatory/intergenic" claim | BP | K (rule) |

## §6 Wild alleles for cultivated barley (new)

| # | candidate | source | status |
|---|---|---|---|
| 6.1 | The GPAT6 high-fiber haplotype and GDSL haplotype B are absent from the five malting cultivars | BP items 8, 14 | K — short |
| 6.2 | Consistent with malting selection | BP (2026-09-24) | D |
| 6.3 | The five cultivars share one haplotype at each gene and lack the high-fiber haplotype. **A plain fact, not a reservation**, with "as expected in modern germplasm" (**no "malting"**, user 2026-09-26) | BP item 8, user | K — reworded |
| 6.4 | GH17 and PHT4;3: cultivars carry the Morex alleles, and Morex is itself elite | BP item 8 | K — short |
| 6.5 | GDSL: reference at the three associated SNPs, alternate elsewhere → not plain reference identity | BP item 14 | K — short |
| 6.6 | What cultivars lack is the Desert-concentrated variation | D11 | D |
| 6.7 | Haplotype frequency across the 136-line elite pool | NEW | D — not approved |
| 6.8 | Landraces to separate domestication from breeding | NEW | D |

## Conclusions

**Written and approved 2026-09-28** (review item 9, `proposed/9_PROPOSED.md`): the last subsection of the Discussion in `../Results_Discussion.md`. Was: deferred until the Discussion is written (user, 2026-09-26). Candidates kept for then: carbon
allocation + protein axis; linkage drag (Hübner and Kantar 2021); β-amylase example + citation TODO;
next steps (validation, replication, functional tests); closing the Introduction's resolution loop.

## Do not use (dead claims)

*Pho*, AP2/ERF and BAHD as candidates · the three 7H genes `7HG0729020/090/100` · "β-glucan across all
seven chromosomes" · "starch required haplotype-level resolution" · the "marginal" SNP class · the
±200 kb LD-decay window · the starch 6H:525.78 sign flip.
