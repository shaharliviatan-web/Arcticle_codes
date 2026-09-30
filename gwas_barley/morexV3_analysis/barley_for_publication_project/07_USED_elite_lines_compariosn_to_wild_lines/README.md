# 07_elite_lines_comparison_to_wild_lines

Which haplotypes of the three publication candidate genes do modern elite barley
cultivars carry? 290 wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*,
Southern Levant) against five European elite malting cultivars, Morex V3.

Created 2026-09-10. Updated 2026-09-14 — see **History** at the end.

**Start here → [`results/tables/results_chapter_numbers.txt`](results/tables/results_chapter_numbers.txt)** —
every number in this step, in prose, ready to paste into the manuscript.

---

## What this step does, in one line

Takes the three genes carried forward by step 06, downloads the same ±1 kb window
from the IPK elite panel, and draws the wild haplotype groups and the elite lines
as **aligned genotype barcodes under the phenotype violins** — so a reader can see
at a glance which wild haplotype each cultivar carries.

**crosshap is not re-run and never sees the elite lines.** The wild haplotype
groups are consumed exactly as step 04 defined them.

## The three genes

| gene | short | trait | window | wild SNPs |
|---|---|---|---|---|
| `HORVU.MOREX.r3.3HG0301300` | GPAT6 | fiber | 3H:544,999,419–545,005,865 | 54 |
| `HORVU.MOREX.r3.5HG0487060` | GH17 | fiber | 5H:462,725,793–462,732,160 | 41 |
| `HORVU.MOREX.r3.3HG0301710` | PHT4;3 | starch | 3H:546,472,709–546,478,564 | 16 |

Windows are gene span ± **1000 bp**, matching step 04's `window_bp` so the wild
and elite barcodes describe the same sequence. The flank is deliberate: it carries
the promoter-proximal region, not only the coding sequence. `01_fetch_elite_vcfs.sh`
asserts this against step 04 and refuses to run if the two ever diverge.

## The five elite lines

All five come from the DivBrowse panel's `accession type = elite lines` subset —
the only subset where cultivar name, breeder and year of release are all present
in the EBI BioSamples record, so **every label in the figures is machine-verifiable.**

| cultivar | SAMEA | released | breeder | called (GPAT6 / GH17 / PHT4;3) |
|---|---|---|---|---|
| **Avalon** | SAMEA110189742 | 2012 | J. Breun | **96.6%** (100 / 97 / 83) |
| **KWS Irina** | SAMEA110189768 | 2012 | KWS Lochow | **94.3%** (95 / 97 / 83) |
| **Odyssey** | SAMEA110189789 | 2012 | Limagrain | **94.3%** (95 / 97 / 83) |
| **Laureate** | SAMEA110189847 | 2016 | Syngenta | **86.4%** (88 / 85 / 83) |
| **LG Diablo** | SAMEA110189781 | 2018 | Limagrain | **98.9%** (100 / 100 / 92) |

Avalon and Odyssey replaced RGT Planet (53.4% called) and Propino (62.5%) on
2026-09-14 — see History.

Three choices behind that list, all auditable in
[`config/elite_lines.tsv`](config/elite_lines.tsv):

- **All five are spring malting types, deliberately.** Winter and spring barley are
  strongly genetically differentiated, so a mixed set would put a
  population-structure axis into the elite block and a barcode difference could
  reflect winter-vs-spring rather than anything about the gene.
- **They span 2012–2018 within one breeding target.** If all five carry the same
  haplotype, that says the modern malting germplasm is fixed for it — a stronger
  statement than five arbitrary lines agreeing.
- **Every line is called at ≥ 85% of the shared sites** (added 2026-09-14). The
  panel is unimputed low-coverage sequencing, so some lines have many no-calls at
  our sites, drawn as grey tiles. All 136 spring elite lines are downloaded once and
  scored (`Table_elite_line_screen.tsv`), the five shown lines are subset from that
  same download, and the run fails if any configured line is below the threshold.
  This filters on genotype quality only — it never looks at which haplotype a line
  carries, so it cannot steer the comparison.

`00_panel_metadata.sh` re-verifies name, accession type and release year for all
five against EBI BioSamples on every run, and **fails the run** if any disagrees.

> **These are not the six lines of the earlier, pre-publication analysis** in
> `USED_Haplotype_Analysis/Wild_And_Elite`. Four of those six (Corvette, Triumph,
> Iris, Gunnar) are samples carrying only an IPK genebank number — HOR 11013,
> HOR 2970, HOR 12070, HOR 10702 — and **no cultivar name in any public record**,
> so their names could not be verified; "Triumph" additionally collides with a
> differently-named sample in the same panel (`SAMEA110189957` = *LG Triumph*).
> The other two, Manchuria (HOR 10695) and Rapid (HOR 8662), are correctly named
> but belong to the `precision collection` — genebank landraces, not elite
> breeding material.

## Where the elite genotypes come from

IPK **DivBrowse barley pangenome v2**, `https://divbrowse.ipk-gatersleben.de/barley_pangenome_v2`
— unimputed SNP variants for the SHAPE2 Core1000 and elite line genotypes, called
against Morex V3. Panel: **1,315 genotypes** = 315 `elite lines` + 804
`precision collection` + 196 with no accession type. **No wild barley
(ssp. *spontaneum*) is in the panel**, so no line here can accidentally be a wild
accession.

The earlier analysis used manual browser downloads. This step calls the API
instead, so the window and the line set are recorded rather than remembered:

```
POST <base>/vcf_export     Content-Type: application/x-www-form-urlencoded
chrom=chr3H&startpos=..&endpos=..&samples=["SAMEA..",..]
```

A JSON body returns HTTP 500 — it must be form-encoded.

### ⚠ The endpoint returns a shifted window

`/vcf_export` returns the **correct number** of variants but reads them from a
**shifted position range**, and the drift is not constant:

| requested | returned |
|---|---|
| 546,472,709–546,478,564 | 158 records spanning **546,476,514–546,481,290** |
| 546,473,000–546,474,000 | 52 records spanning **546,478,016–546,480,350** |

A naive call therefore returns real data for the **wrong part of the genome**,
silently. The pipeline requests a padded window (±50 kb), **asserts that the
returned span brackets the target window on both sides** — widening and retrying
while it does not — and only then trims locally with bcftools. Verified against
the earlier manual export: 158/158 sites recovered, zero missing. Two further
guards exist because both failure modes actually occurred during development:
a truncated transfer is caught by a per-line field-count check, and a
non-bracketing export **aborts the run** rather than producing a quiet wrong answer.

The shift is not constant even across requests: the 136-line pool requests of
2026-09-14 came back with almost exactly the requested span (e.g. PHT4;3 requested
546,422,709–546,528,564, returned 546,422,699–546,528,562), while the earlier
5-line requests were shifted by 1–5 kb. The padding and the assertion stay in
place for exactly that reason; `elite_vcf_provenance.tsv` records requested vs
returned spans for every run.

## ⚠ REF/ALT must be identical, and is checked every run

Both call sets are called against Morex V3 and the wild set was produced without
PLINK allele-flipping, so at a shared position REF and ALT must be **identical**.
If they were ever swapped, `0/0` in one file and `0/0` in the other would denote
**opposite alleles** and every barcode in this step would be silently inverted.

`03_build_matrices.R` classifies every shared position as `identical` / `swapped` /
`ref_differs` / `alt_differs`, writes
[`results/tables/allele_concordance.tsv`](results/tables/allele_concordance.tsv),
and **aborts** on `swapped` or `ref_differs` (`ALLELE_MISMATCH_ACTION="stop"`).

Verified: of 90 shared positions, **88 are identical, 0 swapped, 0 with a differing
REF**. The other 2 are **triallelic** (`alt_differs`): the same reference base, but a
different alternate allele in each panel — GPAT6 `3H:545,003,306` wild G/A vs elite
G/C, and GH17 `5H:462,728,384` wild G/A vs elite G/T. That is a population
difference, not a coding error, so it does not abort. `shared_sites` does not draw
these two columns. The filled versions draw each line's **real genotype** there,
never a fill: `0/0` → Reference (the same base as the wild REF), a call containing
the elite-only allele → Triallelic. What each shown line carries is in
`Table_triallelic_sites.tsv`.

## The SNPs differ between the two files, in both directions

The wild set is a QC-filtered call set on 290 Levantine accessions; the elite
export is unfiltered across 1,315 diverse cultivars. So there are wild-only SNPs
*and* elite-only SNPs. **Three treatments are built and all three are kept** — the
choice between them is deliberately left open:

| version | columns (GPAT6 / GH17 / PHT4;3) | what it draws |
|---|---|---|
| `shared_sites` | 43 / 33 / 12 | only SNPs present in both files with identical REF/ALT. No assumption; costs columns |
| `filled_marked` | 54 / 41 / 16 | all wild SNPs; elite tiles with no elite record get their own **"No elite record"** colour rather than being assumed REF |
| `filled_silent` | 54 / 41 / 16 | all wild SNPs; missing elite records silently become REF — the earlier `bcftools merge -0` behaviour, kept only for comparison |

The two filled versions differ only at the wild-only columns (10 / 7 / 4). The
two triallelic columns are drawn identically in both — with the line's real
genotype — and all five lines are REF (or missing) there, so no purple appears.

**Elite-only SNPs are never drawn.** The wild haplotype groups have no data there,
so they would add columns of grey to every group row. They are counted and
characterised in `Table_site_overlap.tsv` instead.

### Monomorphic sites — what `shared_sites` actually costs (added 2026-09-23)

`Table_site_overlap.tsv` carries, per gene and per file, how many sites are
**monomorphic** over all sites, over the sites `shared_sites` **keeps**, and over the
sites it **removes** (all = kept + removed). "Monomorphic" is judged only over the rows
that enter the comparison, and the reference set differs per file on purpose:

- **wild** — only the accessions crosshap **assigned** to a haplotype group for that gene
  (123 / 107 / 205). Unassigned accessions are dropped from the test and from the figure,
  so a site that varies only among them cannot make a drawn column informative. The
  reference set therefore differs from gene to gene.
- **elite** — only the **5 configured lines**, not the 136-line pool, because those five
  are the only rows drawn.

| gene | wild: mono / all · kept · removed | elite: mono / all · kept · removed |
|---|---|---|
| GPAT6 | 10/54 · 9 of 43 kept · 1 of 11 removed | 235/244 · **42 of 43 kept** · 193 of 201 removed |
| GH17 | 8/41 · 8 of 33 kept · 0 of 8 removed | 286/306 · **33 of 33 kept** · 253 of 273 removed |
| PHT4;3 | 0/16 · 0 of 12 kept · 0 of 4 removed | 155/158 · **12 of 12 kept** · 143 of 146 removed |

Two things follow, and both belong in the manuscript:

1. **Dropping the elite-only sites costs almost nothing for this comparison.** 193 of 201,
   253 of 273 and 143 of 146 of them are monomorphic across the five lines, so they could
   not have separated the cultivars even if the wild set had called them.
2. **Dropping the wild-only sites does cost real wild variation** — 10 of 11, 8 of 8 and
   4 of 4 of them are polymorphic among the assigned accessions. This is the honest price
   of assuming nothing, and it is why the two filled versions are kept.
**Site by site → [`Table_removed_sites.tsv`](results/tables/Table_removed_sites.tsv)**
(`06_removed_sites_table.R`, 2026-09-23). The counts above say how many; that table says
which, with each site's allele tally over the rows that would have been drawn and a plain
`cost` column. Its totals reconcile with `Table_site_overlap.tsv` by construction, and the
script asserts nothing else — it only reads. The sites that cost real wild variation are
**22 in total**: 10 at GPAT6, 8 at GH17 and 4 at PHT4;3, each polymorphic among that gene's
assigned accessions (minor-allele counts 2–21). Both triallelic sites are among them.

3. **The five elite lines are effectively one haplotype at all three genes** — monomorphic
   at 42 of 43, 33 of 33 and 12 of 12 of the drawn columns. That is what
   `config/elite_lines.tsv` anticipated ("if all five carry the same haplotype, that says
   the modern malting germplasm is fixed for it"), but it means **the five cultivars are
   not five independent observations** and must not be written as if they were.

## The figure

One figure per gene per version, `results/figures/<version>/`.

- **Top** — violin plot of the wild haplotype groups (phenotype BLUP ~ group),
  boxplot inside, x-axis label `A` over `(n=85)`, and **pairwise significance
  brackets** (added 2026-09-18). Group **n** is shown; no mean marker.
- **Bottom** — aligned genotype barcodes: one row per wild haplotype group, a gap,
  then one row per elite line. Columns are SNPs in **genomic order**, shared with
  the panel above so the two read as one figure.

Group rows are a **per-SNP majority consensus** over the group's members, ignoring
missing calls. An exact 50/50 tie resolves to the **reference** allele when REF is one of
the tied states (changed 2026-09-24, user decision; until then ties were drawn as missing).
This affects only two GH17 cells: group A at `5H:462,728,332` (18 REF / 18 ALT) and group E
at `5H:462,729,123` (5 / 5); GPAT6 and PHT4;3 have no ties. A tie without REF (ALT vs HET)
would still be missing; none occurs. A consensus is used rather
than one representative accession because wild missingness is high (19% of calls
at PHT4;3), so a single plant's row would show grey tiles that are artefacts of
its sequencing rather than features of the haplotype. The closest real accession
to each consensus is still recorded, in `consensus_representatives.tsv`, as a
verification aid only — it is not published.

Colours are inherited from step 04's `plot_heatmaps.R` so these read the same as
the per-gene heatmaps already in the paper, plus three new states:

| state | colour | note |
|---|---|---|
| Reference | `#FFFACD` | |
| Alternate | `#2F4F4F` | |
| Heterozygous | `#C46210` | **elite only** — the wild VCFs contain no het call at all (verified) |
| Missing | `grey70` | |
| No elite record | white | `filled_marked` only |
| Triallelic (elite-only allele) | `#7B3FA0` purple | filled versions only — a line carrying an allele that is neither the wild REF nor the wild ALT |

Step 04's `gt_to_bin01()` recognised only `0/0`, `1/1` and `./.`, so a het call
would have been folded into "missing". Here het is its own state.

**Each legend lists only the states actually drawn** in that figure, so the
triallelic colour appears only if a shown line really carries the third allele.
**Barcode rows are separated by a white gap** (tiles drawn at 72% of the row
height, added 2026-09-14), so each row reads as its own barcode rather than all
rows fusing into one heatmap.

**Elite lines are not assigned to a haplotype group anywhere in this step** —
crosshap never saw them, so no assignment would be honest. The barcodes are
aligned precisely so the comparison can be left to the reader's eye.

### The brackets

The test is taken unchanged from step 04's combined PDFs
(`04_.../01_scripts/R/plot_combined_pdf.R`), so the two figure sets agree:

- **Wilcoxon rank-sum of every group against the largest group** — a deterministic
  reference, not all pairs, which keeps the number of comparisons at k − 1.
- **Holm correction** across those k − 1 tests, within the gene.
- Symbols `****` ≤ 1e-4, `***` ≤ 1e-3, `**` ≤ 0.01, `*` ≤ 0.05, else `ns`; `ns` is
  shown rather than hidden, as in step 04.
- The Kruskal–Wallis p of step 04 remains the overall test; these are the follow-up.
- Values are written to `Table_pairwise_group_tests.tsv` and into the prose file.

Drawn with plain `geom_segment`/`geom_text` rather than
`ggpubr::stat_pvalue_manual`: ggpubr 0.6.3 errors against ggplot2 4.0.1 here
(a NULL fontsize reaches `grid::gpar`). The statistics are unaffected.

Deliberately not drawn, decided 2026-09-10: gene-model strip, crosshap
marker-group annotation bar, group means, per-elite "% identity to group" column.

## Pipeline

```bash
bash scripts/run_all.sh                    # warm: seconds. cold: ~10 minutes
FORCE_REFETCH=1 bash scripts/run_all.sh    # ignore cached downloads
```

| # | script | does |
|---|---|---|
| 0 | `00_panel_metadata.sh` | DivBrowse sample map + EBI BioSamples → `inputs/panel_metadata.tsv`; **verifies the five configured lines** |
| 1 | `01_fetch_elite_vcfs.sh` | padded DivBrowse export of the **136-line spring elite pool** → integrity + bracket assertion → local trim → `intermediates/elite_vcfs_trimmed/<gene>.pool.vcf.gz` |
| 2 | `02_screen_elite_lines.R` | call rate of every pool line at the shared sites → `Table_elite_line_screen.tsv`; **fails if a configured line is below 85%** |
| 3 | `03_build_matrices.R` | subsets the configured lines from the pool; wild consensus, allele-concordance check, triallelic genotypes, the three versions |
| 4 | `04_figures.R` | violin + aligned barcodes, one figure per gene per version |
| 5 | `05_paper_tables.R` | every paper-facing table + `results_chapter_numbers.txt` |
| 6 | `06_removed_sites_table.R` | **which** sites `shared_sites` drops, one row each, and whether each still carried variation |

Steps 0 and 1 are resumable and cache to `intermediates/`. All parameters live in
**[`config/params.sh`](config/params.sh)** and nowhere else; bash sources it
directly and R reads the same file via `scripts/_load_params.R`, so the two can
never drift.

## Directory map

| path | contents |
|---|---|
| `config/params.sh` | **every parameter**, single source of truth |
| `config/elite_lines.tsv` | the five lines, with the selection rule written out |
| `scripts/` | the six scripts + the shared R param loader |
| `inputs/` | `panel_metadata.tsv` (all 1,315 genotypes), `divbrowse_sample_ids.tsv`, `elite_pool_samples.tsv` (the 136-line pool) |
| `intermediates/` | BioSamples JSON cache, padded + trimmed pool VCFs, matrices |
| `results/figures/<version>/` | one PDF + PNG per gene per version |
| `results/tables/` | all `Table_*.tsv` plus the prose file |
| `logs/` | one log per script, plus `run_all.log` |

## Output tables

| file | one row per | notes |
|---|---|---|
| `results_chapter_numbers.txt` | — | **start here**: everything as prose |
| `Table_elite_lines.tsv` | elite line | name, SAMEA, breeder, year, habit, rationale, call rate |
| `Table_elite_line_screen.tsv` | pool line (136) | call rate per gene and pooled — the evidence behind the line choice |
| `Table_triallelic_sites.tsv` | triallelic site × line | what each shown line carries at the two triallelic sites |
| `Table_pairwise_group_tests.tsv` | group comparison | Wilcoxon vs the largest group + Holm p — the bracket statistics |
| `Table_panel_composition.tsv` | panel subset | what the 1,315-genotype source panel is |
| `Table_gene_windows.tsv` | gene | coords, window, locus, lead SNP, q, η² from step 04 |
| `Table_site_overlap.tsv` | gene | wild / elite / shared / triallelic / wild-only / elite-only counts (wild-only and elite-only = no record at that position in the other file), **+ monomorphic-site counts per file** (added 2026-09-23, see below) |
| `Table_allele_concordance.tsv` | shared site | REF/ALT in both files + verdict |
| `Table_haplotype_groups.tsv` | gene × group | n, mean, median, sd of the phenotype |
| `Table_elite_genotypes_wide__<gene>.tsv` | elite line | genotype at every shared site |
| `Table_removed_sites.tsv` | removed site × file | **every site `shared_sites` drops**: position, alleles, indel flag, genotype tally over the rows that would have been drawn, `state` (polymorphic / monomorphic / undetermined) and a plain `cost` column |
| `Table_removed_sites_summary.tsv` | gene × file × reason | the same, counted: sites, polymorphic, monomorphic, undetermined, indels |
| `elite_vcf_provenance.tsv` | gene | requested vs returned span, pad, counts, timestamp |
| `consensus_representatives.tsv` | gene × group | closest real accession — verification only |

## What the figures show

Percent agreement between each elite line and each wild group consensus, over the
shared sites where the elite line has a call (from `Table_elite_genotypes_wide__*.tsv`):

| gene | result |
|---|---|
| **GPAT6** (fiber) | All five elites match the two **low-fibre** wild groups — **A at 95–98%** (mean −0.120, n=85) and **B at 88–91%** (mean −0.081, n=20) — and the **high-fibre** group **C at only 39–45%** (mean +0.281, n=18), the group carrying the entire GPAT6 signal. **The high-fibre haplotype is essentially absent from these five modern malting cultivars.** *(Corrected 2026-09-23: this row previously named only group A and omitted B. A and B are both low-fibre and are not distinguished by the elite lines, so the finding is about the absence of the high-fibre haplotype, not about group A specifically. Recomputed from `Table_elite_genotypes_wide__*.tsv`; no analysis changed.)* |
| **GH17** (fiber) | No clean match: 46–53% to groups A and D, 29–34% to B, C, E. The elites do not correspond to any one wild haplotype. *(Updated 2026-09-24 after the tie rule changed: B/C/E was 26–34%. A and D unchanged in range. Recomputed from the `shared_sites` matrices.)* |
| **PHT4;3** (starch) | All five match group **C at 100%** — but see the warning below before reading anything into that. |

### ⚠ Reference-identity is not haplotype sharing

**Morex is itself an elite six-row malting cultivar, so REF is an elite genome.**
Elite lines therefore tend to look reference-like by construction, and whichever
wild haplotype happens to sit closest to Morex will appear "most elite" whether or
not anything was shared. Measured here:

| gene | elite %REF (called sites) | reading |
|---|---|---|
| **GH17** | **100% REF**, all five lines | the "match" is reference identity, nothing more |
| **PHT4;3** | **100% REF**, all five lines | group C is the only 100%-REF wild consensus, so the 100% match to C is an artefact of that, **not** evidence that elites carry the low-starch haplotype |
| **GPAT6** | **12–13% REF** | **the informative case** — the elites differ from Morex at most sites and *still* land on the low-fibre groups A (14% REF) and B (12% REF) rather than group C (65% REF), so the match is real signal, not reference attraction |

**Only GPAT6 supports a claim about shared haplotypes.** For GH17 and PHT4;3 the
honest statement is that the elite lines carry the Morex reference haplotype in
this window. PHT4;3 now rests on 10–11 called sites of 12 per line (the first selection had
lines with as few as 3).

## Caveat for the manuscript

> The figure shows which haplotypes elite lines carry; it does not by itself show
> that breeding selected them.

## History

| date | change |
|---|---|
| 2026-09-10 | Step created. Five lines chosen on fame: RGT Planet, Laureate, LG Diablo, KWS Irina, Propino. Three SNP-matching versions built. Allele-concordance check added (0 swapped of 90 shared positions). |
| 2026-09-14 | **Call-rate screen added** (`02_screen_elite_lines.R`): the whole 136-line spring elite pool is downloaded once and scored; configured lines must reach ≥ 85% called. RGT Planet (53.4%) and Propino (62.5%) failed and were **replaced by Avalon (96.6%) and Odyssey (94.3%)**; Explorer (88.6% overall, 66.7% at PHT4;3) was passed over. Scripts renumbered 02→03, 03→04, 04→05. |
| 2026-09-14 | **Triallelic columns** (GPAT6 `3H:545,003,306`, GH17 `5H:462,728,384`) now drawn with each line's real genotype in the filled versions instead of "No elite record" / assumed REF. Checked across all 136 pool lines: **none carries the elite-only third allele** (all 0/0 or missing), so these columns draw as REF and the triallelic colour never appears. |
| 2026-09-14 | **White gap between barcode rows** (tile height 0.72); legends list only states present. |
| 2026-09-18 | **Pairwise significance brackets added to the violins**, using step 04's test (Wilcoxon vs the largest group, Holm-corrected); statistics also written to `Table_pairwise_group_tests.tsv`. `config/params.sh` `STEP_ROOT` updated after the directory was renamed to `07_USED_...`, without which no script could run. |
| 2026-09-24 | **Consensus tie rule: ties → REF** (user decision; was: ties → missing). `consensus_row()` in `03_build_matrices.R`; steps 03–06 re-run (steps 00–02 not re-run, cached downloads unchanged). Effect: 2 GH17 consensus cells (group A `5H:462,728,332`, group E `5H:462,729,123`) turn from missing to REF; GH17 agreement to B/C/E 26–34% → 29–34%; `consensus_representatives.tsv` GH17 rows A and E change (verification only). GPAT6 and PHT4;3 matrices, all elite genotypes, allele concordance, site overlap and every number used in the manuscript are unchanged. Fig. 4 (`08_USED_creating_figures/make_figure_4.R`) re-rendered from the new matrices. |
| 2026-09-27 | **Stale tie-rule text corrected** (user approval). Section [G] of `results_chapter_numbers.txt` still said "exact 50/50 ties are drawn as missing", because the text in `05_paper_tables.R` had not been updated on 2026-09-24. It now states the rule the code applies: a tie resolves to REF when REF is one of the tied states, otherwise missing. Only `05_paper_tables.R` was re-run: the other 20 tables are byte-identical, and only that line and the "Generated" timestamp changed in `results_chapter_numbers.txt`. The same rule was aligned in the 7H branch (`03_01_.../scripts/02_elite_comparison.R`, `03_violin_plus_barcode.R`) with no change to any output. |

## Conventions

- All paths absolute. `bcftools` = `/usr/bin/bcftools`. R 4.1.2, vcfR 1.16.0,
  pheatmap 1.0.13.
- `TMPDIR=/mnt/data/shahar/.tmp`; nothing written under `/tmp` or `$HOME`.
- Scripts write only named files under `intermediates/` and `results/`.
- Chromosome naming: step 04 uses `3H`, DivBrowse uses `chr3H`; renamed on import.
