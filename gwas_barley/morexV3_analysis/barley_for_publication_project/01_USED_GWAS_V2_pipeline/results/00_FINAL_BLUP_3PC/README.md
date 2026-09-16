# GWAS results — chosen configuration (v3)

Wild barley (*Hordeum vulgare* ssp. *spontaneum*), 290 accessions, Southern Levant.
Reference MorexV3. 7,110,996 SNPs. Traits: beta-glucan, fiber, protein, starch.

**This directory is the live analysis. Everything superseded is under
`_ARCHIVE_superseded_2026-09-08/`, never deleted.**

Reorganised 2026-09-08. Next step: candidate-gene search inside the loci in `02_loci_FINAL/`.

---

## Locked configuration

| parameter | value |
|---|---|
| Phenotype value | BLUP |
| PCs as EMMAX covariates | 3 |
| Kinship | EMMAX aIBS, 290 × 290 |
| SNPs tested | 7,110,996 |
| LD pruning (for PCA/kinship) | `--indep-pairwise 1000kb 1 0.2` → 111,017 SNPs |
| alpha | 0.10 |
| **Significance threshold** | **−log10p = 6.0454** (p < 9.008e-07), Bonferroni over 111,017 independent tests |
| Significant SNPs | 52 |

Only significant SNPs are used. There is no "marginal"/"suggestive" class — that
concept was dropped by decision, and every script and table referring to it has
been removed or archived.

---

## Locus definition (locked 2026-09-08)

| parameter | value | why |
|---|---|---|
| lead SNPs | `--clump-p1 9.008e-07` | only genome-wide-significant SNPs may lead a locus |
| membership | `--clump-p2 1` | **LD only, no p-value condition** — a causal variant need not itself be significant |
| LD | `--clump-r2 0.5` | |
| reach | `--clump-kb 2000` | ±2 Mb, so 4 Mb maximum span |
| **contiguity** | **severed at the first gap > 50 kb** | a locus must be *connected*; a correlated SNP across an empty stretch is not part of it |
| **iteration** | **repeat until no significant SNP is orphaned** | `--clump` is winner-take-all, so a severed significant SNP would otherwise belong to nothing and vanish from the tables and figures |

**Result: 36 loci, median span 93.5 kb (mean 219.8), 3,423 member SNPs, and all
52 significant SNPs inside a locus** — asserted by the script, converged in 3 passes.
Total gene-search space **7.9 Mb** (see `02_loci_FINAL/README.md`).

The gap rule is what makes these usable. Plain `--clump` reports min→max of every
member it absorbed, which gave 100–1,977 kb spans bridging uncorrelated gaps.
Severing at 50 kb brings the median locus to 93.5 kb. The gap tolerance sets the
total search space (7.9 Mb here), so state it explicitly in the methods.

---

## Directory map

| directory | contents |
|---|---|
| `00_snp_level/` | SNP-level results: the 52 significant SNPs, per-trait summary, analysis parameters, top-15-per-chromosome |
| `01_assoc/` | per-trait `SNP` + `P` files (7,110,996 rows each) — the shared GWAS output that feeds clumping |
| **`02_loci_FINAL/`** | **the live locus definition — start here.** Tables, figures, raw PLINK output, and the 2 extra peaks |
| `03_pc_selection/` | PC variance table used to justify 3 PCs |
| `04_diagnostics/qq/` | QQ plots and genomic inflation |
| `_ARCHIVE_superseded_2026-09-08/` | every superseded run, each with the scripts that produced it |

---

## Scripts

All in `../../scripts/`. Live flow, in order:

| script | does |
|---|---|
| `01`–`05` | VCF → PLINK → LD pruning → PCA/kinship → EMMAX 24-cell grid |
| `06`–`09` | QQ/lambda, Manhattans, summary tables, v2↔v3 comparison |
| `11_chosen_config_tables_v3.R` | SNP-level package → `00_snp_level/` |
| **`31_loci_clump_iterative.R`** | LD clumping + 50 kb gap rule + iteration to zero orphans → `02_loci_FINAL/` |

| **`33_manhattan_loci_painted.R`** | painted Manhattans → `02_loci_FINAL/figures/` |
| **`36_paper_tables_loci.R`** | publication tables → `02_loci_FINAL/tables/Table_*.tsv` |

Reproduce the locus set from scratch:

```bash
MAX_GAP=50000 KB=2000 Rscript scripts/06_loci_FINAL/31_loci_clump_iterative.R
Rscript scripts/06_loci_FINAL/36_paper_tables_loci.R
Rscript scripts/06_loci_FINAL/33_manhattan_loci_painted.R
```

Scripts for archived analyses live **inside** their archived result directory, so
each superseded run is self-contained and still reproducible.

`plink` must be `/usr/local/bin/plink` — the bare name resolves to a broken v0.76.
