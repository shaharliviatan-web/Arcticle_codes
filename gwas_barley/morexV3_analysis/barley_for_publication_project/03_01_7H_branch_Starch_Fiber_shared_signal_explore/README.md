# 03_01_7H_branch_Starch_Fiber_shared_signal_explore — GDSL esterase at the shared 7H signal (MGmin = 3, ε = 0.9)

Haplotype analysis of **`HORVU.MOREX.r3.7HG0729030`** — the GDSL esterase at the shared
7H fiber/starch signal — at **MGmin = 3, ε = 0.9**, with the full elite-cultivar
comparison and paper table set.

**This is the analysis Results Chapter 4 is built on** (adopted and approved by the user
2026-09-24), and it draws **Fig. 5** (`08_USED_creating_figures/make_figure_5.R`). It replaced
the first version of this branch (MGmin = 2, ε = 0.6), under which the gene is **null**.

**Folder history (2026-09-24, user decision).**
- This analysis was built in the subfolder `REPLACEMENT_ANALYSIS_MGmin3_eps0.9/` and was
  moved up to be the branch root. Only `scripts/config.sh` (`TEMP_ROOT = BRANCH`) and one
  table label (`Table_gene_windows.tsv` → `crosshap_run`) changed; `run_all.sh` was re-run
  from here and reproduced every table (ε 0.9 rows identical).
- The first version (MGmin = 2, ε = 0.6) is archived, whole and still runnable, in
  **`archive_v1_MGmin2_eps0.6/`**. Its steps 00–04 do not depend on MGmin and **remain valid
  and cited by Ch. 4**: gene position and distances, InterPro/Swiss-Prot annotation, the local
  association profile (promoter vs gene body), signal-SNP LD, and the direct genotype split.
  Fig. 5 reads its `signal_snps.tsv`. Only its steps 05–06 (MGmin 2 crosshap, ε sweep, figures)
  are superseded by this analysis.
- **ε 0.6 outputs were removed** from this folder the same day (user): its caches, figures and
  table rows. Only the chosen ε = 0.9 is run. The comparison that chose 0.9 over 0.6 is kept in
  `results/tables/param_grid_eps0.05-1.5_MGmin2-3.tsv` (§2).

290 wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*, Southern Levant),
Morex V3. Built 2026-09-24.

> **Inputs.** Apart from its own files, it reads step 04, step 07 and step 01 **read-only**
> and writes nothing outside this folder (the archive subfolder included).

**Start here → [`results/tables/results_chapter_numbers.txt`](results/tables/results_chapter_numbers.txt)** — every number below, in prose, ready to paste.

---

## 1. Why this analysis exists

At 7H ~573.6 Mb the **fiber and starch lead SNPs are 154 bp apart** — the only shared
signal in the study. Both LD spans are near-zero (185 bp and 1,563 bp), so step 03, which
searches the LD span with `FLANK_BP = 0`, returned the locus as **gene-empty**.

The first version of this branch (`archive_v1_MGmin2_eps0.6/README.md`) found that a gene does sit there:
`7HG0729030`, **578 bp** from the fiber lead and **732 bp** from the starch lead, with no
other annotated gene for 112 kb. It is a **secreted GDSL esterase/lipase** — a family that
deacetylates arabinoxylan (rice DARX1, BS1), polymerises cutin (tomato CD1), and in barley
controls hull–caryopsis attachment. All four associated SNPs sit in its **promoter**
(554–763 bp upstream of the TSS; the gene is on the minus strand) while the six SNPs
**inside the gene body are null** — a cis-regulatory signature.

**But the haplotype analysis at the pipeline's MGmin = 2, ε = 0.6 returned nothing**
(fiber KW p = 0.16, starch p = 0.36). The cause is mechanical, not biological: DBSCAN
assigns **all four GWAS signal SNPs to marker group 0 (noise)** and builds the haplotypes
from two *null* SNPs at the far end of the window. Only 109 of 290 accessions are assigned.

At **MGmin = 3** the four correlated signal SNPs (mutual r² 0.44–0.94) form a proper
marker group, and the analysis works.

## 2. Why MGmin = 3, ε = 0.9 — the selection procedure

> **Scope.** This justifies the parameters **for this one gene only**. The 55-gene
> pipeline (step 04) keeps **MGmin = 2, ε = 0.6**, unchanged and untouched; its parameters
> and the reasoning behind them are documented in
> `04_USED_haplotype_analysis_crosshap/README.md` and belong there, not here.
> **The two analyses are separate, each with its own methods and parameters.**

The parameters were selected by a **single procedure in two steps**, evaluated on the
14 ε × 2 MGmin grid in `results/tables/param_grid_eps0.05-1.5_MGmin2-3.tsv`:

> **Step 1 — MGmin is the smallest value at which the four GWAS signal SNPs form a marker
> group.** That is **MGmin = 3**.
>
> **Step 2 — ε is then the value maximising the assignment rate.** That is **ε = 0.9**
> (246 of 290 accessions, tied with ε = 1.0; the smaller, more conservative ε is taken).

**Neither step uses the haplotype test's p-values.** The grid's `kw_p` column played no
part; the selection reads only `sig_kept` and `assigned`.

### Why step 1 is needed — assignment rate alone does not work

This is worth stating explicitly, because the obvious simpler rule fails. Ranked purely by
assignment, the grid looks like this:

| rank | MGmin | ε | signal SNPs carried | assigned / 290 | |
|---|---|---|---|---|---|
| 1 | 2 | 0.10 | **0 / 4** | **253 (87.2%)** | ← max assignment, but useless |
| **2** | **3** | **0.90** | **4 / 4** | **246 (84.8%)** | **← chosen** |
| 3 | 3 | 1.00 | 4 / 4 | 246 (84.8%) | tied |
| 4–8 | 3 | 0.4–0.8 | 3 / 4 | 239 (82.4%) | |
| 9 | 3 | 1.50 | 4 / 4 | 231 (79.7%) | |
| 10 | 3 | 1.20 | 4 / 4 | 197 (67.9%) | |
| 15 | 2 | 1.50 | 4 / 4 | 147 (50.7%) | only MGmin 2 setting keeping all 4 |
| 18 | 2 | 0.60 | **0 / 4** | 109 (37.6%) | the pipeline's setting |

**A pure maximum-assignment rule would select MGmin = 2, ε = 0.10** — which carries **none**
of the four GWAS signal SNPs and returns KW p = 0.92. Assignment rate measures how much of
the panel the test uses; it says nothing about whether the test is aimed at the right
variants. Step 1 supplies that, and only then is assignment rate a meaningful criterion.

Within MGmin = 2, the only setting that carries all four signal SNPs is ε = 1.5, and it
assigns just 147 of 290 — barely half the panel. **MGmin = 3 is therefore not a free
choice made to improve the result; it is the only value at which the locus can be tested
on a majority of the panel at all.**

### Is this phenotype-free?

Step 2 is purely genotypic. Step 1 refers to the **GWAS result**, which is prior to and
independent of the haplotype test — but it is fair to call it *GWAS-informed* rather than
blind to phenotype in the absolute sense. The distinction that matters:

- **No p-value from the haplotype test entered the choice.** Verified: the selection uses
  only `sig_kept` and `assigned`, never `kw_p`.
- Requiring a haplotype analysis to represent the variants the GWAS already identified is
  a **specification of what is being tested**, not a search for significance.

State it that way in the methods and it is defensible. Claiming the choice was entirely
phenotype-blind would not be.

### The plateau

At MGmin = 3 the outcome is **identical across ε 0.4–0.8** (3/4 SNPs, 239 assigned) and
again across **ε 0.9–1.0** (4/4, 246 assigned). These are flat regions, not spikes, so the
choice does not sit on a knife edge. The grid table is the record of this; the ε 0.6 run
itself (caches, figures) was removed 2026-09-24.

Recorded from the ε 0.6 cache before it was deleted: **at MGmin = 3, ε = 0.6 the single
marker group is exactly the three genome-wide significant SNPs** (7H:573,606,306, 573,606,460,
573,606,491), for both traits. So MGmin = 3 is the smallest value that groups the significant
SNPs even at the pipeline's ε; raising ε to 0.9 adds the fourth signal SNP (7H:573,606,282)
and 7 more accessions (239 → 246).

## 3. Results

### Haplotypes (MGmin = 3, ε = 0.9)

| trait | SNPs in window | assigned | groups | KW H (df 1) | **KW p (raw)** | η² | `delta_top_bottom_sd` |
|---|---|---|---|---|---|---|---|
| fiber | 27 | **246 / 290 (85%)** | 2 — 212 \| 34 | 9.46 | **0.00209** | 0.035 | 0.54 |
| starch | 27 | **246 / 290 (85%)** | 2 — 212 \| 34 | 10.70 | **0.00107** | 0.040 | 0.56 |

| trait | group | n | mean | median | SD |
|---|---|---|---|---|---|
| fiber | **B** | 34 | **+0.0981** | +0.0810 | 0.304 |
| fiber | A | 212 | −0.0528 | −0.0946 | 0.273 |
| starch | A | 212 | +0.8109 | +2.1436 | 5.636 |
| starch | **B** | 34 | **−2.3861** | −1.9350 | 5.500 |

**The same minority haplotype raises fiber and lowers starch** — the trade-off, in one
haplotype, at one gene.

Three checks that this is not a lucky parameter:

1. **All four GWAS signal SNPs are retained** in marker group MG1 (mean r² = 0.74).
2. **The ε = 0.6 result is identical across five consecutive ε values (0.4–0.8)** — a
   plateau, not a spike. See `results/tables/param_grid_eps0.05-1.5_MGmin2-3.tsv`.
3. **The n = 34 haplotype is exactly the genotype split.** It is the *identical* 34
   accessions carrying the minor allele at both lead SNPs in the first version's direct
   split — 34 shared, 0 discordant either way. Mean fiber BLUP matches to four decimals.

### Would it survive BH in context?

| trait | raw p | tests in trait | BH q |
|---|---|---|---|
| fiber | 0.00209 | 20 | **0.0042** |
| starch | 0.00107 | 6 | **0.0013** |

Significant for both traits.

### Elite cultivars

Same five spring malting cultivars and the same 136-line pool as step 07, re-fetched for
this window. 27 wild SNVs × 131 elite records → **22 shared**, including **all four signal
SNPs**; 22 of 22 shared positions allele-identical, **0 swapped, 0 triallelic**.

**No elite cultivar carries the minority haplotype.** 14 Reference calls, 6 no-calls,
**0 Alternate** across 5 lines × 4 signal SNPs.

| cultivar | 573606282 | 573606306 | 573606460 | 573606491 | call rate (22 sites) |
|---|---|---|---|---|---|
| Avalon | REF | REF | REF | REF | 95.5% ✓ |
| KWS Irina | REF | REF | REF | REF | 95.5% ✓ |
| Odyssey | — | — | REF | REF | 81.8% |
| Laureate | — | — | REF | REF | 54.5% |
| LG Diablo | REF | REF | — | — | 77.3% |
| **wild hap A (212)** | REF | REF | REF | REF | |
| **wild hap B (34)** | **ALT** | **ALT** | **ALT** | **ALT** | |

> **Call-rate caveat.** Step 07 requires ≥ 85%, applied *pooled over its three genes*. At
> this gene alone only **Avalon and KWS Irina** clear it — and they are the two lines
> called at all four signal SNPs. The conclusion rests on them.

### Sites removed by the shared-site rule

| file | removed as | sites | polymorphic | monomorphic |
|---|---|---|---|---|
| wild | wild-only | 5 | 5 | 0 |
| elite | elite-only | 109 | 33 | 76 |

Polymorphism is judged over the same reference sets step 07 uses — wild over the 246
**assigned** accessions, elite over the 5 **configured** lines. Per-site detail is in
`Table_removed_sites.tsv`. None of the 27 wild SNPs is monomorphic over the assigned set.

## 4. Figures

| folder | contents |
|---|---|
| **`results/figures/shared_sites/`** | **the main figure** — violin + aligned elite barcodes, step 07's published layout |
| `results/figures/crosshap/` | step 04's combined per-gene PDF (tree + violin) |
| `results/figures/heatmaps/` | step 04's LD / haplotype heatmaps |
| `results/figures/elite/` | barcode-only version, no violin |

`shared_sites/` reproduces `07_.../results/figures/shared_sites/` — same layout, palette
(`#FFFACD` reference, `#2F4F4F` alternate, `grey70` missing), title and caption format,
and the same statistic: **Wilcoxon of every group against the largest group,
Holm-corrected** (one comparison here, A vs B; `**` for both traits). Values in
`Table_pairwise_group_tests.tsv`.

**The one deliberate departure from step 07:** red triangles ▲ mark the four columns of
crosshap marker group **MG1** — the GWAS signal SNPs. Step 07 draws no marker-group
annotation; these were added on request and are explained in the figure caption.

Figures exist for **ε 0.9 only** (ε 0.6 removed 2026-09-24). The manuscript figure is **Fig. 5**,
built from this folder by `08_USED_creating_figures/make_figure_5.R` (TAG layout; triangles
under the three genome-wide significant SNPs only).

## 5. Layout

| path | contents |
|---|---|
| `scripts/config.sh` | every parameter — single source of truth |
| `archive_v1_MGmin2_eps0.6/` | the first version of this branch (MGmin 2), archived 2026-09-24; steps 00–04 still valid and cited |
| `scripts/01_run_mgmin3.R` | crosshap at MGmin = 3, ε 0.9 (ε 0.6 dropped 2026-09-24); caches, stats, step-04 figures |
| `scripts/02_fetch_elite_vcf.sh` | elite panel VCF for this window, from DivBrowse |
| `scripts/02_elite_comparison.R` | wild-vs-elite barcodes (barcode-only figure) |
| `scripts/03_violin_plus_barcode.R` | the `shared_sites` violin + barcode figure |
| `scripts/04_paper_tables.R` | the full step-07 table set |
| `scripts/05_results_chapter_numbers.R` | the prose summary |
| `scripts/run_all.sh` | all of the above, in order |
| `scripts/06_param_grid.R` | the MGmin × ε grid of §2 (added 2026-09-30, user-approved; **not** in `run_all.sh`, it changes no result). Regenerates `results/tables/param_grid_eps0.05-1.5_MGmin2-3.tsv` byte-identically (checked 2026-09-30); before it, no script in the project wrote that table |
| `Cache/<trait>/<gene>/MGmin_3/eps_<ε>/HapObject.rds` | crosshap objects, step 04's cache layout |
| `work/` | per-gene raw + imputed VCFs, elite VCF |
| `results/figures/`, `results/tables/`, `logs/` | outputs |

### Tables

| file | contents |
|---|---|
| **`results_chapter_numbers.txt`** | **everything below, in prose** |
| `analysis_parameters.tsv` | every constant, for the methods section |
| `param_grid_eps0.05-1.5_MGmin2-3.tsv` | 14 ε × 2 MGmin × 2 traits, **for this gene** — the evidence for §2; written by `scripts/06_param_grid.R`; an Online Resource of the manuscript. Its `kw_p`, `eta2` and `delta` columns are recorded but were not used for the choice |
| `mgmin3_gene_results.tsv` | headline stats per trait × ε |
| `mgmin3_haplotype_groups.tsv`, `mgmin3_haplotype_assignment.tsv` | group stats, per-accession haplotype |
| `Table_gene_windows.tsv` | step-07 schema: window, locus, lead, stats |
| `Table_haplotype_groups.tsv` | step-07 schema |
| `Table_site_overlap.tsv` | shared / wild-only / elite-only / monomorphic counts |
| `Table_removed_sites.tsv`, `Table_removed_sites_summary.tsv` | every dropped site, and whether polymorphic |
| `Table_allele_concordance.tsv`, `Table_triallelic_sites.tsv` | REF/ALT agreement at shared positions |
| `Table_elite_lines.tsv`, `Table_elite_line_screen.tsv` | the five cultivars and their call rates |
| `Table_elite_genotypes_wide__<gene>.tsv` | elite genotype per shared site, signal SNPs flagged |
| `Table_pairwise_group_tests.tsv` | Wilcoxon vs largest group, Holm |
| `elite_vcf_provenance.tsv` | what was requested from DivBrowse and when |

## 6. Methods — what differs from the rest of the project

Everything is the project's existing method, run unmodified, **except the two parameters**:

| step | method | source |
|---|---|---|
| per-gene VCFs | `bcftools view`, gene ± 1 kb, `-m2 -M2 -v snps`, 290 keep-list, doubled→single reheader | step 04 `01_`/`02_` |
| LD | `plink --r2 square --keep-allele-order` | inside step 04's `run_crosshap.R` |
| haplotyping | crosshap, step 04's `run_crosshap.R` **unmodified** | step 04 |
| statistics | Kruskal–Wallis, η², `delta_top_bottom_sd` | step 04 |
| figure statistics | Wilcoxon vs largest group, Holm | step 07 `04_figures.R` |
| elite comparison | consumed haplotypes, majority consensus (exact ties → REF when REF is tied, else missing), `CHROM:POS:REF:ALT` matching | step 07 `03_`/`04_`/`05_`/`06_` |
| figures | step 04's and step 07's renderers | steps 04, 07 |

**The differences, in full:**

1. **`MGmin = 3`** instead of 2 — the substantive change, justified in §2.
2. **`ε = 0.9`** instead of 0.6 — retains all four signal SNPs; the ε 0.6 result is significant
   too (grid table; its run was removed 2026-09-24).
3. **Red triangles** marking marker group MG1 in the `shared_sites` figures.
4. **Elite fetch pad is 2 kb**, not step 07's 50 kb — this region is variant-dense and the
   larger request streams for many minutes without covering the window.
5. **p-values are raw** — one pre-specified gene, so no BH is applied here (step 04
   corrects within trait *across genes*). The in-context BH q values are given in §3.
**Consensus tie rule aligned with step 07 (2026-09-27, user approval).** `02_elite_comparison.R`
and `03_violin_plus_barcode.R` still resolved every exact tie to missing, step 07's rule before
2026-09-24. They now use step 07's current rule, as `make_figure_5.R` already did: a tie resolves to
REF when REF is one of the tied states, otherwise missing. This gene has no tie (wild consensus
cells: 46 REF, 42 ALT, 0 missing), so both scripts were re-run and all 19 tables are byte-identical.
The method is therefore the same for Ch. 3 and Ch. 4, with nothing to separate in M&M.

6. **InterProScan input is stop-codon-stripped** (first version, `archive_v1_MGmin2_eps0.6/`) — the PGSB proteome
   carries a trailing `*` that the EBI service rejects.

## 7. Adopted for Chapter 4 (2026-09-24) — what that means

- Chapter 4's gene count and the statement that this locus has no candidate gene both change.
- `7HG0729030` is **not** in step 03's candidate list (`FLANK_BP = 0`; a flank of ≥ 1 kb
  would have included it). Adopting it means stating that it was recovered by a targeted
  follow-up, not by the standard search.
- **This changes nothing about the 55-gene pipeline.** Step 04 keeps MGmin = 2, ε = 0.6
  for every other gene, and its parameter choice is documented in its own README. The two
  analyses stand side by side, each with its own methods and parameters, and the
  manuscript should present them that way: the genome-wide haplotype screen at the
  pipeline's settings, and this single targeted follow-up at settings justified on its own
  terms (§2).
- The honest framing is: the GWAS association (p = 5.7e-08 fiber / 3.7e-07 starch), the
  promoter position, and the GDSL function are the primary evidence. The haplotype result
  supports it; it is not independent of it.

## 8. Reproduce

```bash
bash scripts/run_all.sh
source scripts/config.sh && Rscript scripts/06_param_grid.R   # the §2 grid only (~1-2 min)
```

Reads step 01's association files, step 04's helper code and staged imputed VCF, and step
07's panel metadata and elite-line config — **all read-only**. The elite fetch is skipped
if `work/elite/7HG0729030.elite.vcf.gz` already exists.

## Conventions

- All paths absolute; `TMPDIR=/mnt/data/shahar/.tmp`; nothing under `/tmp` or `$HOME`.
- `plink` must be `/usr/local/bin/plink` — the bare name resolves to a broken v0.76.
