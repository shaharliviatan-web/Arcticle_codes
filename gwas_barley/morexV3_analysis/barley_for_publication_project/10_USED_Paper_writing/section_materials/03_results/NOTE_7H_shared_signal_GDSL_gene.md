# The 7H shared fiber/starch signal and `7HG0729030` — note for Results and Discussion

Written 2026-09-17. **Candidate for the manuscript — structure not yet approved.**

> **Moved 2026-09-24.** The analysis this note describes (MGmin 2, ε 0.6) is now archived in
> `03_01_7H_branch_Starch_Fiber_shared_signal_explore/archive_v1_MGmin2_eps0.6/`; all paths below are relative to it. The branch root now holds the MGmin = 3,
> ε = 0.9 haplotype analysis that Ch. 4 is built on (see its `README.md`). This note's MGmin 2
> haplotype section is superseded; its gene, annotation, promoter-profile and split facts still hold.

Source analysis: `03_01_7H_branch_Starch_Fiber_shared_signal_explore/archive_v1_MGmin2_eps0.6/` (created 2026-09-17,
reproducible, `bash scripts/run_all.sh`, ~2 min). That branch reads steps 01/04/05 read-only and
**feeds nothing** — no main table has been changed.
**Start there → `results/tables/genotype_group_tests.tsv`; full reasoning in its `README.md`.**

---

## The finding

At **7H ≈ 573.6 Mb** the fiber lead (`7H:573606306`, −log10p 7.2451) and the starch lead
(`7H:573606460`, 6.4367) lie **154 bp apart** — the tightest fiber/starch co-localization in the
study. Step 03 returned both loci **gene-empty**, because the LD spans are 185 bp and 1,563 bp and
`FLANK_BP = 0`.

Relaxing the flank **for this locus only** recovers a single gene:
**`HORVU.MOREX.r3.7HG0729030`** — 578 bp from the fiber lead, 732 bp from the starch lead, with
**no other annotated gene for 112 kb**. A flank of ≥ 1 kb would have caught it.

## Why it is a credible fiber gene

Secreted **GDSL esterase/lipase** (SGNH hydrolase). Four databases converge — CDD `cd01837`
plant-specific SGNH lipase at **E = 2.6e-120**, Pfam **PF00657**, Gene3D, SUPERFAMILY — and both
**Phobius and SignalP** predict a signal peptide, so the protein is **apoplast/cell-wall targeted**,
which is where arabinoxylan acetylation, cutin deposition and hull adhesion are controlled.
**Family-level call, not a named ortholog** (43–45 % Swiss-Prot identity).

Three independent routes from this family to grain fibre, one of them in barley itself:

| organism | gene | route to fibre |
|---|---|---|
| rice | **DARX1** | deacetylates **arabinoxylan** — the other major dietary-fibre polysaccharide of barley grain |
| rice | **BS1** | xylan deacetylase, secondary wall patterning |
| tomato | **CD1 / GDSL1** | polymerises **cutin**, which reports into the insoluble-fibre fraction — the same logic as GPAT6 |
| **barley** | GDSL esterase/lipase | wax/cutin deposition and **hull–caryopsis attachment**, which directly changes measured grain fibre |

**For starch the link is compositional, not mechanistic**: if the cell-wall/hull fraction shifts,
starch as a percentage of grain shifts inversely.

## The cis-regulatory signature

The gene is on the **minus strand**, so its TSS is the high coordinate and its promoter runs
**into** the signal. All four associated SNPs (+554 to +763 bp from the TSS, mutual r² 0.44–0.94,
MAF 0.14–0.19) sit in the putative promoter. The **six SNPs inside the gene body are null in both
traits** (−log10p ≤ 1.6). That pattern **excludes a coding change** and points to regulatory
variation.

## The effect

Splitting the panel on the two leads (no haplotype clustering): 34 minor-allele accessions vs 167
major-allele accessions —

| trait | difference | Wilcoxon p |
|---|---|---|
| fiber | **+0.60 SD** | 1.1e-03 |
| starch | **−0.62 SD** | 6.2e-04 |

A near-symmetric inverse, consistent with **one causal variant read out by two correlated traits**.

> These tests are **not kinship-corrected** and are not comparable to the GWAS p-values. They give
> direction and effect size only.

**Direction confirmed by the kinship-corrected GWAS (checked 2026-09-22).** The lead-SNP betas agree:
for the minor allele, fiber `lead_beta` = +0.119 (`7H:573606306`) and starch `lead_beta` = −2.145
(`7H:573606460`), in BLUP units. **The minor allele raises fiber and lowers starch.** `lead_beta` is
the effect of the minor allele (A1) (beta sign convention in `01_.../results/00_FINAL_BLUP_3PC/README.md`).
<!-- src: 01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv -->

## Why it is not in the haplotype chapter's gene set

At the pipeline's fixed **ε = 0.6, crosshap discards all four signal SNPs as DBSCAN noise** and
returns null (fiber p = 0.16, starch p = 0.36). The single marker group that forms is built from two
*null* SNPs, so the test evaluates a partition unrelated to the signal. This is the documented ε trap
in a low-SNP-density window (27 SNPs), **a methodological failure, not evidence against the gene**.

> The branch's ε sweep is **diagnostic only**. The ε = 1.5 row must never be reported as a finding —
> selecting ε by outcome is the forking path the main pipeline closed by fixing ε.

**Consequence for the manuscript:** this gene cannot be presented as a fourth member of the
"23 haplotype-significant genes → 3 presented" line. It did not come through that pipeline. It needs
its own framing as a **targeted follow-up of a pre-specified observation** (the 154 bp
co-localization), with its own method — the direct genotype split.

> **Decided 2026-09-17: the gene is NOT adopted into the main pipeline.** The `FLANK_BP = 0` rule
> stays as it is, step 03 is not re-run, and no downstream step changes. This gene remains a single
> locus examined on its own, outside the span-based candidate set. Either it is reported that way in
> its own chapter, or it is left out entirely.

## How it changes the shared-signal story

It does **not** make 7H "the only shared signal" — there are three fiber/starch regions, and they
resolve at three different levels. That is the story:

| region | what is shared | resolution |
|---|---|---|
| **7H ≈ 573.6 Mb** | the loci **overlap**; leads 154 bp apart | **one gene, both traits** — `7HG0729030` |
| **3H ≈ 543.5–546.5 Mb** | neighbouring loci, 595 kb apart | **two genes, one per trait** — GPAT6 (fiber) and PHT4;3 (starch), both presented in Ch. 3, GPAT6 with the elite-line comparison |
| **1H ≈ 330 Mb** | a starch locus nested between two fiber loci, ~58 kb | **no gene** |

Three regions, three resolutions — a much stronger basis for the carbon-allocation argument than any
one of them alone.

## Honesty requirements if this is written

1. **State that the flank rule was relaxed for this locus.** It is a targeted follow-up of a
   pre-specified co-localization, not a fishing expedition — but it must be said, and the cost of
   `FLANK_BP = 0` acknowledged as known and accepted.
2. **Not causal.** These are tag SNPs; the promoter reading depends on the annotated TSS; the causal
   variant may be an ungenotyped correlated one. What *is* excluded is a coding change.
3. **Family-level annotation**, not a named ortholog.
4. **The starch link is compositional**, not mechanistic.
5. **High missingness** — 85 of 290 accessions lack a call at one or both leads.
6. **No kinship correction** in the group tests.
7. **The haplotype analysis is negative**, and the ε closest to working was chosen after the fact.

## Supersedes

The retired `04_.../07_fiber_starch_tradeoff_direction/`, built on three genes 156–187 kb from the
same lead. This gene is ~270× closer. **Do not reuse the old three-gene evidence.**
