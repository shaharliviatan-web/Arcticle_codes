# ⚠ NOT FOR WRITING USE — JUST ADDITIONAL INFO ⚠

**Do not write anything from this file into the manuscript.** It is background information only,
kept in case it is useful. Nothing here has been reproduced from a project script, and none of it
has been approved for the paper.

---

Canonical barley grain-quality genes vs. the 36 GWAS loci.
Reported 2026-09-17 by a separate session; work done in a scratchpad, nothing in the project touched.
Filed verbatim, as written.

---

Done. Nothing in the project was touched — all work in the scratchpad.

## 1. Your list: mostly right, with gaps

I mapped 29 canonical genes onto Morex V3 by protein sequence (86–100% identity). Verdict on your table:

| your entry | verdict |
|---|---|
| β-glucan: CslF6 → 7H | ✓ correct, 7H:366.7 Mb |
| β-glucan: CslH1 → 7H | ✗ 2H:475.7 Mb |
| β-glucan: Glb1/Glb2 → 7H | ✗ split — Glb1 is 1H:396.0 Mb, Glb2 is 7H:624.4 Mb |
| starch: sex6/SSIIa 7HS | ✓ 7H:94.1 Mb |
| starch: Wax/GBSSI 7HS | ✓ 7H:17.6 Mb |
| starch: amo1/SSIIIa 1HS | ✓ 1H:84.4 Mb |
| protein: HvNAM-1 6H | ✓ 6H:51.4 Mb |
| protein: HvNAM-2 2H | ✓ 2H:160.8 Mb |
| protein: Hor1/Hor2 1H | ✓ 1H:1.8–2.6 Mb (tip of 1HS) |

Your protein row is fully correct — and notably, a web search will tell you hordeins are on
"chromosome 5". That's the old barley nomenclature; old 5 = modern 1H. You had it right.

Missing from your list: the 2H CslF cluster (F3/F4/F8/F10 at ~180 Mb) — and the largest-effect
β-glucan QTL in barley is on 2H, coincident with it. For starch you're missing SBEIIa/SBEIIb (2H),
ISA1/2/3 (7H/1H/5H), AGPase subunits (1H/5H/7H), SSI and SBEI (7H).

## 2. Distances: nothing is close. Not one.

Minimum distance between any canonical gene and any of our 36 significant loci: **37.9 Mb**.

| trait | nearest canonical gene | distance |
|---|---|---|
| fiber | HvGlb2 → fiber_L18 | 37.9 Mb |
| starch | SBEI → starch_L05 | 54.3 Mb |
| β-glucan | HvGlb1 → betaglucan_L01 | 83.5 Mb |
| protein | — | no canonical gene is even on the same chromosome (both protein loci are on 3H; all protein genes are 1H/2H/6H) |

Everything else runs 90–556 Mb. There is no proximity to interpret.

## 3. Sub-threshold scan — one real near-miss

Best SNP within ±1 Mb of each gene (threshold = 6.0454):

| gene | trait | best SNP | −log10p | gap to gene |
|---|---|---|---|---|
| γ-hordein-3 / hordein cluster | protein | 1H:1719723 | 5.66 | 123 kb |
| HvNAM-1 | protein | 6H:51874045 | 4.75 | 497 kb |
| GBSSI / Wax | starch | 7H:17200853 | 4.57 | 386 kb |
| SSI | starch | 7H:43240602 | 4.33 | 820 kb |
| all others | | | ≤ 4.0 | |

The hordein hit is the one worth knowing about. −log10p 5.66 against a 6.0454 threshold — short by
0.38 log units — and 123 kb from γ-hordein-3, inside our LD-decay distance (~188 kb at r²=0.2).
MAF 0.113, so not a rare-allele fluke.

But it's a singleton. In the surrounding 4 Mb (2,248 SNPs) exactly one SNP exceeds 4.0; the next best
is 3.82. A real causal variant at MAF 0.11 should drag neighbours up through LD. It doesn't. Call
rate is also lowish — 231 of 290 individuals.

In fairness, though: our own betaglucan_L02 — the strongest β-glucan locus we report — is nearly as
spiky (max 7.60, only 1 SNP >5, 4 member SNPs total). So spikiness alone doesn't disqualify it in
this dataset.

HvNAM-1 is the opposite shape: lower max (4.75) but a genuinely broader elevation — 30 SNPs >3 and
4 SNPs >4 in its window, versus 3 and 1 for the hordein region. It looks more like a real weak
signal, just far below threshold.

## One methodological finding worth your attention

GBSSI, SSIIa (sex6) and SSI are in the Morex V3 low-confidence gene set, not the high-confidence set.
My first BLAST against the HC proteome returned garbage (~30% identity) until I re-searched against
HC+LC, which gave 98.5–100%.

That means: if one of our loci had landed directly on Waxy or sex6, step 03 would have reported "no
gene" — those genes aren't in the annotation the candidate search uses. It didn't affect anything
here (we're 480–556 Mb away), but it's a real blind spot in the gene-search step.

## Bottom line

No canonical gene for any of the four traits is validated by our GWAS. The wild panel is segregating
at different loci than cultivated barley — which is a legitimate and interesting result, not a
failure, and it's consistent with what we already found for β-glucan.

The only thing I'd consider mentioning in the paper is the 1HS hordein near-miss for protein — with
the singleton caveat stated honestly.

Two genes I couldn't map: lys5 (no clean reference sequence retrievable) and C-hordein (the only
entry is a 105 aa fragment). Not forced, as you said.
