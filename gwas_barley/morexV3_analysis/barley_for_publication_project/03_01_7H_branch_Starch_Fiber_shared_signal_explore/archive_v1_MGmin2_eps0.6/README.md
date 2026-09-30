# 03_01_7H_branch_Starch_Fiber_shared_signal_explore — ARCHIVED first version (MGmin = 2, ε = 0.6)

> **ARCHIVED 2026-09-24 (user decision).** This was the branch root until then. The branch root
> (`../`) now holds the **MGmin = 3, ε = 0.9** haplotype analysis that Results Chapter 4 and
> Fig. 5 are built on; read `../README.md` first.
>
> - **Still valid and still cited by Ch. 4 (MGmin-independent):** steps 00–04 — gene position and
>   distances, the InterPro / Swiss-Prot annotation, the local association profile (promoter vs
>   gene body), signal-SNP LD, and the direct genotype split. Fig. 5 reads `results/tables/signal_snps.tsv`.
> - **Superseded:** steps 05–06 — the MGmin = 2 crosshap run and the ε sweep below, and their
>   figures. Under MGmin = 2 the gene is null; do not quote that as the result.
> - Moved, not changed: only `scripts/config.sh` (`BRANCH` now points here). Steps 01, 02 and 04
>   were re-run from this location on 2026-09-24 and reproduced their outputs byte-for-byte.
> - **Do not run `scripts/run_all.sh` casually:** step `03b` starts a remote NCBI nr BLAST (not
>   cached; it previously ran > 1 h 40 min and was abandoned).

An **exploratory branch** off step 03. Asks what gene lies under the one place in the
study where **fiber and starch share a signal** — 7H ~573.6 Mb, where the two lead SNPs
are **154 bp apart** — a locus step 03 returned as *gene-empty* because both LD spans are
near-zero and `FLANK_BP = 0`.

Created 2026-09-17. 290 wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*,
Southern Levant), Morex V3.

> ⚠ **This branch feeds nothing.** It does not modify step 03, step 04 or any pipeline
> output, and no result here has entered the main tables. It exists to evaluate whether
> the gene is worth adopting; that decision has not been made.

**Start here → [`results/tables/genotype_group_tests.tsv`](results/tables/genotype_group_tests.tsv)** — the main result.

---

## Why this branch exists

| | |
|---|---|
| fiber lead | `7H:573606306`, −log10p **7.2451** (`fiber_L17`) — 2nd strongest fiber signal in the study |
| starch lead | `7H:573606460`, −log10p **6.4367** (`starch_L05`) |
| distance between them | **154 bp** |
| `fiber_L17` LD span | **185 bp** (3 SNPs) |
| `starch_L05` LD span | 1,563 bp (6 SNPs) |
| genes returned by step 03 | **0** — both spans are too small to contain one |

Two traits with a well-documented compositional trade-off, one shared signal, and no
candidate gene. That was worth a second look.

## The gene

**`HORVU.MOREX.r3.7HG0729030`** — 7H:573,604,051–573,605,728, **minus strand**, 1,678 bp,
369 aa.

Because it is on the minus strand its **TSS is the high coordinate (573,605,728)** and its
promoter runs upward — straight into the signal:

```
  573,604,051 ──────── gene body ──────── 573,605,728 │ ← TSS
       3' end       6 SNPs, ALL NULL          5' end  │   ┌──── promoter ────►
                    fiber p 0.08–1.39                 │   ▼
                    starch p 0.06–1.59                │  573,606,282 … 573,606,491
  ◄──── transcription ────                            │  fiber 7.25 / starch 6.44
                                                      │  next gene: 112 kb away
```

| | |
|---|---|
| gene → fiber lead | **578 bp** |
| gene → starch lead | **732 bp** |
| gene → `starch_L05` span edge | 554 bp |
| **next annotated gene in either direction** | **112.4 kb** (`7HG0729040`) |
| flank that would have caught it in step 03 | **≥ 1 kb** |

For comparison, the three genes the retired `04_.../07_fiber_starch_tradeoff_direction/`
was built on sit **156.4, 185.2 and 187.1 kb** from the same lead. This gene is ~270×
closer.

## Annotation

Checked with the **same three-source chain as step 05**, not Swiss-Prot alone.

### Check 1 — UniProt Swiss-Prot (DIAMOND BLASTP)

All five hits are GDSL esterase/lipases — 43–45% identity over ~88% coverage,
E ≈ 1e-89 to 1e-95 (`results/tables/gene_annotation_swissprot.tsv`). `--max-target-seqs 5`
matches step 05 check 1 exactly.

### Check 2 — EBI InterProScan 5 (all member databases + GO)

`results/tables/gene_annotation_interpro.tsv`

| database | signature | InterPro | E-value |
|---|---|---|---|
| **CDD** | `cd01837` **SGNH_plant_lipase_like** | **IPR035669** GDSL lipase/esterase-like, **plant** | **2.6e-120** |
| PANTHER | PTHR22835 | — | 3.3e-118 |
| Gene3D | G3DSA:3.40.50.1110 SGNH hydrolase | IPR036514 SGNH hydrolase superfamily | 1.9e-71 |
| **Pfam** | **PF00657** GDSL-like Lipase/Acylhydrolase | IPR001087 GDSL lipase/esterase | 8.0e-31 |
| SUPERFAMILY | SSF52266 SGNH hydrolase | — | 5.7e-11 |

**GO:** `GO:0016788` hydrolase activity, acting on ester bonds.

**Four independent databases converge on GDSL/SGNH**, and CDD's top signature is
**plant-specific** at E = 2.6e-120 — much stronger support than the 43–45% Swiss-Prot
identity alone implied.

> **The protein is secreted.** **Phobius** and **SignalP** independently predict a
> **signal peptide**, and Phobius calls the mature chain non-cytoplasmic. The enzyme is
> targeted to the **apoplast / cell wall** — precisely where an enzyme acting on
> arabinoxylan acetylation, cutin, or hull adhesion must be. This is not visible from
> Swiss-Prot alone and materially strengthens the fibre candidacy.

**One dissenting signature:** PANTHER labels it "ZINC FINGER FYVE DOMAIN CONTAINING
PROTEIN". PANTHER family names are inherited and are frequently wrong for plant proteins;
five other lines of evidence disagree. Noted, not weighted. The MetaCyc/Reactome pathway
list attached to the GO term is generic ester-hydrolase inheritance and carries no
information here.

### Check 3 — NCBI nr BLASTP: **not applicable, by protocol**

Step 05 defines check 3 as **rescue only**, run on a derived residual set:

```
residual = {Swiss-Prot WEAK_OR_NO_HIT}  INTERSECT  {no InterPro match}
```

This gene fails **both** criteria — a confident Swiss-Prot hit (45.3% id, E = 7.6e-95)
*and* strong InterPro matches from four databases. Under the project's own rules it is
**not a residual gene**, so check 3 is **not run** — exactly as it is not run for the 16
check-1-annotated genes in the main table.

A trial nr run was attempted anyway and abandoned after ~1 h 40 min with no return (NCBI's
remote queue was the bottleneck; the service itself was reachable, HTTP 200). It is not
needed and its absence is not a gap. `scripts/03b_interpro_and_nr.sh` keeps the nr block so
the step can be run on demand, but `run_all.sh` does not depend on its output.

> Still a **family-level** assignment, **not** a named ortholog. GDSL is large — 100+
> members in Arabidopsis, 114 in rice.

### Why the family is a strong fiber candidate

| organism | gene | what it does | relevance |
|---|---|---|---|
| rice | **DARX1** | deacetylates **arabinoxylan** side chains | arabinoxylan is the other major dietary-fibre polysaccharide in barley grain besides β-glucan |
| rice | **BS1** (Brittle Leaf Sheath 1) | xylan deacetylase, secondary wall patterning | direct cell-wall remodelling |
| tomato | **CD1 / GDSL1** (cutin synthase) | polymerises **cutin** | cutin reports into the insoluble-fibre fraction — same logic as GPAT6 in the main gene set |
| **barley** | GDSL-motif esterase/lipase | wax + cutin deposition, **hull–caryopsis attachment** | hull adherence directly changes measured grain fibre |
| Arabidopsis | GDSL esterase/lipase family | pectin/homogalacturonan acetylesterases | cell-wall ester modification |

Three independent routes to fibre — hemicellulose acetylation, cutin, and hull retention —
and one of them is documented **in barley itself**.

**For starch the link is indirect.** GDSL esterases are not starch enzymes. The most
defensible reading is **compositional**: if the cell-wall/hull fraction shifts, starch as a
percentage of grain shifts inversely. That is consistent with the direction observed below,
and with the two leads being 154 bp apart — one causal variant, two correlated readouts.

## Main result — direct genotype split

Splitting the 290 accessions on the two lead SNPs (no haplotype clustering):

| group | n | fiber mean | fiber SD | starch mean | starch SD |
|---|---|---|---|---|---|
| **REF/REF (major)** | **167** | −0.0679 | 0.265 | +1.0698 | 5.431 |
| **ALT/ALT (minor)** | **34** | **+0.0981** | 0.304 | **−2.3861** | 5.500 |
| discordant | 4 | | | | |
| missing call | 85 | | | | |

| trait | ALT − REF | in SD | KW H (df=1) | **KW p (raw)** | η² | `delta_top_bottom_sd` |
|---|---|---|---|---|---|---|
| **fiber** | +0.166 | **+0.60** | 10.588 | **0.001138** | 0.048 | 0.596 |
| **starch** | −3.456 | **−0.62** | 11.733 | **0.000614** | 0.054 | 0.619 |

**The same minor allele class raises fiber by 0.60 SD and lowers starch by 0.62 SD** — an
almost symmetric inverse, from one 4-SNP block, at one gene.

> **Statistics deliberately match step 04**, so this branch adds no new method: the same
> **Kruskal–Wallis** test, the same **η² = (H − k + 1)/(n − k)**, and the same
> **`delta_top_bottom_sd` = (mean_top − mean_bottom)/SD**. The only difference is how
> groups are *defined* — observed genotype at the two lead SNPs instead of DBSCAN
> haplotype clustering. (`delta_top_bottom_sd` is unsigned; the signed ALT−REF direction
> is in the "in SD" column.)
>
> p-values are **raw**: one pre-specified gene × 2 traits, so no BH correction applies
> (step 04 corrects within trait *across genes*). They are **not kinship-corrected** and
> are not comparable to the GWAS p-values (5.7e-08 / 3.7e-07) — they quantify direction
> and effect size only.

## The signal SNPs

Four SNPs, all in the promoter, in strong mutual LD
(`results/tables/signal_snp_ld_r2.tsv`):

| | 573606282 | 573606306 | 573606460 | 573606491 |
|---|---|---|---|---|
| **573606282** | 1.000 | 0.440 | 0.544 | 0.499 |
| **573606306** | 0.440 | 1.000 | 0.779 | 0.726 |
| **573606460** | 0.544 | 0.779 | 1.000 | **0.936** |
| **573606491** | 0.499 | 0.726 | **0.936** | 1.000 |

| SNP | offset from TSS | fiber −log10p | starch −log10p | MAF |
|---|---|---|---|---|
| 573606282 | +554 | 2.35 | 3.13 | 0.188 |
| **573606306** | **+578** | **7.25** ← fiber lead | 5.34 | 0.170 |
| **573606460** | **+732** | 5.34 | **6.44** ← starch lead | 0.147 |
| 573606491 | +763 | 5.71 | 6.26 | 0.140 |

MAF 0.14–0.19 — healthy, and above the median of the study's 52 significant SNPs (0.067).

## Haplotype analysis — negative, for a methodological reason

At the pipeline's fixed **ε = 0.6**, crosshap finds **nothing**: fiber KW p = 0.16,
starch p = 0.36.

The reason is not biological. **DBSCAN assigns all four GWAS signal SNPs to marker group 0
(noise).** The single marker group that forms is built from two *null* SNPs at the far end
of the window, so the KW test evaluates a partition unrelated to the signal
(`results/tables/crosshap_marker_groups.tsv`).

This is the documented ε trap: **ε is a Euclidean radius over each SNP's full r² *profile*,
not an r² cutoff**, so its stringency depends on SNP density. With only 27 SNPs in the
window, ε = 0.6 is too strict — even though the four SNPs are in real LD (r² 0.44–0.94).

| ε | marker groups | **signal SNPs kept** | assigned | fiber KW p | starch KW p |
|---|---|---|---|---|---|
| 0.4 | 3 | 3 of 4 | 151 | 0.106 | 0.058 |
| **0.6 (fixed)** | 1 | **0 of 4** | 109 | 0.160 | 0.361 |
| 0.8 | 2 | 0 of 4 | 96 | 0.069 | 0.660 |
| 1.0 | 3 | 0 of 4 | 76 | 0.107 | 0.685 |
| 1.2 | — | crosshap internal error | | | |
| 1.5 | 2 | **4 of 4** | 147 | 0.051 | 0.046 |
| 2.0 | — | crosshap internal error | | | |

> ⚠ **The sweep is diagnostic, not a result.** Selecting ε by outcome is exactly the
> forking path the main pipeline closed by fixing ε = 0.6. The ε = 1.5 row must not be
> reported as a finding. It is here only to show *why* the gene is missed.

**Even at ε = 1.5 the haplotype test reaches only raw p ≈ 0.05** — far weaker than the
direct genotype split (p = 0.0006–0.0011) on the same data. Two reasons: crosshap discards
143 of 290 accessions as unassigned, and its haplotype "C" (n = 11) is a strict **subset**
of the true minor-allele class — **34** accessions are ALT at both leads, but 20+ were
dumped into the unassigned bin over missing data elsewhere in the window.

**This is why the direct split is the primary analysis here, not crosshap.**

## Figures

`results/figures/crosshap/` and `results/figures/heatmaps/` — rendered with step 04's own
unmodified renderers at the **exploratory ε = 1.5**, headers stamped accordingly.
`results/figures/previews/` holds PNG previews of the informative panel (page 5).

Page 5 shows the structure clearly: marker group **MG2 = the four GWAS SNPs**
(meanR2 0.74), defining haplotype **C**, which has the **highest fiber and the lowest
starch** of the three groups — the trade-off in one picture.

*Page 1 of each PDF (the CrossHap tree) fails to render at this ε with a crosshap-internal
`subscript out of bounds`. A library issue, not a data problem.*

## Caveats, stated plainly

1. **Not causal.** These are tag SNPs. The promoter interpretation depends on the annotated
   TSS being correct, and the causal variant could be any correlated ungenotyped variant.
   What the data *do* exclude is a coding change in this gene — the six in-gene SNPs are flat.
2. **Family-level annotation**, not a named ortholog (43–45% identity).
3. **The starch link is compositional**, not mechanistic.
4. **High missingness** — 85 of 290 accessions lack a call at one or both leads.
5. **No kinship correction** in the group tests.
6. **Haplotype analysis does not support the gene**, and the ε value that comes closest was
   chosen after the fact.

## Layout

| path | contents |
|---|---|
| `scripts/config.sh` | every parameter, single source of truth |
| `scripts/00…06`, `run_all.sh` | the pipeline, ~2 min end to end |
| `inputs/` | the gene's protein sequence |
| `intermediates/` | per-gene raw + imputed VCFs, LD matrix, crosshap working dir |
| `results/tables/` | all output tables |
| `results/figures/` | crosshap + heatmap PDFs and PNG previews |
| `logs/` | one log per script |

### Tables

| file | contents |
|---|---|
| **`genotype_group_tests.tsv`** | **the main result** — effect sizes and p-values |
| `genotype_group_summary.tsv` | group means/SDs/medians |
| `genotype_phenotype_per_accession.tsv` | per-accession genotype + both phenotypes |
| `gene_position.tsv` | coordinates, strand, TSS |
| `distances_to_signals.tsv` | gene vs every 7H locus |
| `neighbouring_genes.tsv` | how isolated the gene is |
| `gene_annotation_swissprot.tsv` | the ten Swiss-Prot hits |
| `local_association_profile.tsv` | per-SNP p (both traits), region label, MAF |
| `signal_snps.tsv` | the SNPs carrying signal |
| `signal_snp_ld_r2.tsv`, `window_ld_r2_matrix.tsv` | LD |
| `crosshap_epsilon_sweep.tsv` | the ε diagnostic |
| `crosshap_marker_groups.tsv` | which SNPs DBSCAN kept vs discarded |
| `crosshap_haplotype_assignment_eps1.5.tsv` | per-accession haplotype at ε = 1.5 |

## Methods — nothing new for the M&M

Every method here is one the project already uses elsewhere, so this branch needs **no new
methods section**:

| step | method | already used in |
|---|---|---|
| per-gene VCFs | `bcftools view` (gene ±1 kb, `-m2 -M2 -v snps`), 290-sample keep-list, doubled→single reheader | step 04 `01_`/`02_` |
| LD | `plink --r2 square --keep-allele-order` | inside step 04's `run_crosshap.R` |
| association / MAF | EMMAX `.assoc`, `morexV3_290_freq.frq` | step 01 |
| gene extraction | Ensembl Plants r62 GFF3, col3 `gene`, 1H–7H, `bedtools intersect` | step 03 `02_extract_genes.sh` |
| annotation check 1 | DIAMOND BLASTP vs Swiss-Prot, `--very-sensitive --max-target-seqs 5` | step 05 check 1 |
| annotation checks 2–3 | EBI InterProScan 5 (vendored `iprscan5.py`), NCBI nr remote BLASTP | step 05 checks 2, 3 |
| haplotyping | crosshap, MGmin = 2, ε = 0.6, step 04's `run_crosshap.R` **unmodified** | step 04 |
| ε diagnostic | genotype-only parameter scan | step 04 `diagnostics/param_grid_scan.R` |
| figures | step 04's `plot_combined_pdf.R`, `plot_heatmaps.R` **unmodified** | step 04 |
| **statistics** | **Kruskal–Wallis**, **η²**, **`delta_top_bottom_sd`** | **step 04** |

**The single sentence the M&M needs:** *"the same Kruskal–Wallis test and effect-size
measures used for the haplotype analysis, applied to groups defined by observed genotype
at the two lead SNPs rather than by DBSCAN haplotype clustering."*

Two footnotes: p-values here are **raw** (one pre-specified gene, no BH), and the
InterProScan submission uses a **stop-codon-stripped** copy of the protein because the
PGSB proteome carries a trailing `*` that the EBI service rejects.

## Reproduce

```bash
bash scripts/run_all.sh          # ~2 min
```

Reads step 01's association files, step 04's staged imputed VCF and crosshap helper code,
and step 05's Swiss-Prot DIAMOND database — **all read-only**. Writes only under this
directory.

## Conventions

- All paths absolute; `TMPDIR=/mnt/data/shahar/.tmp`; nothing under `/tmp` or `$HOME`.
- `plink` must be `/usr/local/bin/plink` — the bare name resolves to a broken v0.76.
