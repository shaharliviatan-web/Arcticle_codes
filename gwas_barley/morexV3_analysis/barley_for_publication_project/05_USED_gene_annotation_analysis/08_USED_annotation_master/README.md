# 08_USED_annotation_master

Consolidated functional annotation for the **23 haplotype-significant genes** of the
wild-barley grain-quality GWAS. Merges the three evidence sources into one auditable
table, assigns a single best functional call, a source, and an evidence tier.

**Rebuilt 2026-09-10** against step-04 run `loci_LDspan_eps06_V4`.
The v1 (45-gene) version is in `_ARCHIVE_v1_45genes_2026-09-09/`.

**Start here → [`results/tables/Table_significant_genes_paper.tsv`](results/tables/Table_significant_genes_paper.tsv)**

## Inputs

| source | file | covers |
|---|---|---|
| Swiss-Prot BLASTP | `../05_USED_gene_annotation/results/tables/fdr_gene_annotation_swissprot.tsv` | 23 genes |
| InterProScan | `../06_USED_interpro_domains/results/tables/fdr_gene_interpro.tsv` | 23 genes |
| NCBI nr rescue | `../07_USED_BLASTP_genes_with_no_annotation_left/results/tables/fdr_residual_nr.tsv` | the 2 residual |
| gene list + statistics | `../05_USED_gene_annotation/inputs/fdr_genes_table.tsv` | locus, q, effect sizes |

Nothing is filtered or dropped here — all 23 rows are kept.

## How the call is chosen

`final_call` follows a fixed priority:
**confident Swiss-Prot → InterPro domain family → confident nr → none.**

## Where the annotation came from

There is **no confidence tier**. It was removed on 2026-09-10 because it conflated two
independent things — how sure we are *what the protein is*, and how plausible the gene is
*as a candidate for the trait* — and in collapsing them to HIGH/MEDIUM/LOW it hid the
evidence. What is reported instead is plain provenance.

The three sources are checked in a fixed order; that order is the **check number**:

| check | source | what it gives |
|---|---|---|
| **1** | UniProt Swiss-Prot BLASTP | curated protein names |
| **2** | InterProScan | domain families + GO |
| **3** | NCBI nr BLASTP | last-resort rescue, residual genes only |

| column | meaning |
|---|---|
| `final_call` | the annotation used, taken from the first check that produced one |
| `annotation_source` | `check1_swissprot` / `check2_interpro` / `check3_nr` / `none` |
| `annotation_from_check` | 1, 2, 3, or NA |
| `n_sources_with_call` | how many of the three produced a call |
| `call_check1_swissprot` | Swiss-Prot's own call — **always reported** |
| `call_check2_interpro` | InterPro's own call — **always reported** |
| `call_check3_nr` | nr's own call — **always reported** |
| `needs_review_two_calls` | TRUE when ≥2 sources gave a call and the texts differ |

Current set: **16 from check 1, 5 from check 2, 2 with no call.** 16 genes have calls
from two sources; **4 are flagged for joint review.**

### The 4 flagged for review

Both texts are kept side by side so they can be judged together — often they are the
same protein under different vocabulary:

| gene | check 1 (Swiss-Prot) | check 2 (InterPro) |
|---|---|---|
| 5HG0487060 | Glucan endo-1,3-β-glucosidase 5 | X8 domain; Glycoside hydrolase family |
| 3HG0301260 | Retrovirus-related Pol polyprotein | Reverse transcriptase, RNA-dep. DNA pol |
| 4HG0341560 | Protein VERNALIZATION 3 | PEBP-like superfamily |
| 5HG0487090 | Tuliposide A-converting enzyme 1 (46% id) | Carboxylesterase; α/β hydrolase fold |

`needs_review_two_calls` **decides nothing and never changes `final_call`.** The
comparison is textual and errs toward flagging.

`sp_pident`, `sp_qcovhsp`, `sp_evalue` are reported so identification quality can be
judged directly. `NEAR_FLOOR_PIDENT = 40` no longer alters any label.

## Identification is not biological candidacy

A confident identification means *"we are sure what the protein is"*, **not** *"this gene
plausibly affects the trait"*. Kinesin KIN-7D is a solid Swiss-Prot call and an unlikely
starch candidate; the two judgements are independent.

Trait candidacy is therefore **manual**, made on the reported evidence plus literature.
The columns `trait_candidate_strength` and `candidate_rationale` are written **empty**
by the pipeline. To record decisions, create
`results/tables/trait_candidate_calls.tsv` with columns
`gene_id, trait, trait_candidate_strength, candidate_rationale`; `01_build_paper_table.R`
picks it up automatically on the next run and will not overwrite it.

## Scripts

| script | does |
|---|---|
| `00_build_master.R` | merges the three sources → `fdr_annotation_master.tsv` |
| `01_build_paper_table.R` | assigns `serial_no` by ascending q, attaches haplotype statistics and figure stems → the two paper tables |

```bash
bash scripts/run_all.sh                    # 00 only
Rscript scripts/01_build_paper_table.R     # then the paper tables
```

## Outputs (`results/tables/`)

| file | one row per | contents |
|---|---|---|
| `fdr_annotation_master.tsv` | (trait, gene) | every field from all three sources |
| **`Table_significant_genes_annotated.tsv`** | (trait, gene) | the full paper table — 55 columns, all evidence |
| **`Table_significant_genes_paper.tsv`** | (trait, gene) | trimmed to the 22 columns a manuscript needs |

Both paper tables are **ordered by ascending BH q** and carry `serial_no` (1 = most
significant). **`serial_no` is the stable public identifier** for a gene in the manuscript
and in the `Significant_genes/` figure filenames.

Key columns: `kw_p_raw` (nominal Kruskal–Wallis), `fdr_q` (BH within trait),
`eta_squared` (share of phenotype rank variance explained), `delta_top_bottom_sd`
(highest- minus lowest-mean haplotype group, in phenotype SD), `group_sizes`,
`dist_to_lead_bp`, and the full Swiss-Prot / InterPro / nr evidence.
