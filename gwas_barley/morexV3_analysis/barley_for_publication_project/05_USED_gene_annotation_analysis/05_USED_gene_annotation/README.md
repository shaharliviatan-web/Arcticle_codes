# 05_USED_gene_annotation

Reproducible functional annotation of the haplotype-significant genes from the
wild-barley grain-quality GWAS (Morex V3 reference, PGSB r3 annotation). This is a
NEW top-level pipeline step; it does not reuse code from steps 01-04.

**Re-run 2026-09-10** on step-04 run `loci_LDspan_eps06_V4`: **23 genes**
(beta-glucan 6, fiber 12, starch 5; protein has no testable gene in that run).
The v1 45-gene results are in `_ARCHIVE_v1_45genes_2026-09-09/`.

Primary evidence here = **DIAMOND BLASTP vs UniProt Swiss-Prot** (curated protein
names). GO / InterPro / NCBI nr are the NEXT, online phase and are NOT done here.

## Input set

Genes that passed BH-FDR in step 04, locked by `scripts/00_lock_fdr_genes.R` from:
`../04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Significant_genes/significant_genes.tsv`

- **23 (trait, gene) rows = 23 unique genes**: beta-glucan 6, fiber 12, starch 5.
  No gene appears under two traits in this set.
- That source is already filtered to `significant_fdr == TRUE`, so no filtering
  happens here; the columns are asserted instead.
- The step-04 statistics (`fdr_q`, `eta_squared`, `delta_top_bottom_sd`, group counts)
  travel with the gene list and are carried through to the output table - a functional
  call is only interpretable next to the effect it explains.

> **Changed 2026-09-10.** v1 read a retired step-04 review CSV and carried a
> `legacy_annotation` column from it (manual Ensembl lookups, no recorded evidence,
> never used as the call). Both are gone.

## Sources and versions

- Proteome (LOCAL, not downloaded, not translated from GFF):
  `/mnt/data/Barley_2021/morexV3/gene_annotation/Hv_Morex.pgsb.Jul2020.HC.aa.fa`
  - Upstream PGSB Morex V3 r3 HC proteome that Ensembl Plants r62 redistributes;
    same gene models and `HORVU.MOREX.r3` IDs as the step-03 GFF3.
  - Headers `>HORVU.MOREX.r3.<gene>.<isoform>`. A copy is kept in `inputs/` for
    provenance. A single trailing stop-codon `*` is stripped from each peptide
    (PGSB carries it; Ensembl reports length without it; DIAMOND should not see it).
    Check: `HORVU.MOREX.r3.1HG0079280.1` = 206 aa (matches Ensembl).
  - One representative peptide per gene: the LONGEST isoform, ties broken by lowest
    isoform index (= canonical `.1` for single-isoform genes). One gene in the current
    23 is multi-isoform: `5HG0487060` -> `.1` chosen (equal lengths, lowest index).
- Swiss-Prot: UniProtKB/Swiss-Prot Release **2026_03** of 02-Sep-2026,
  `uniprot_sprot.fasta.gz` from
  `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/`
  (575,748 sequences), downloaded 2026-09-09 for the 23-gene run. See `intermediates/versions.txt`.
  (Corrected 2026-09-30, user-approved: this line said Release 2026_01 of 28-Jan-2026 with
  574,627 sequences, from the earlier 45-gene run; `intermediates/versions.txt` and the file
  itself show 2026_03.)
- DIAMOND **2.0.14** (`/usr/bin/diamond`). makeblastdb/blastp also present at
  `/usr/bin` but unused here.
- Ensembl Plants FTP is unreachable from this server, so the peptide file was NOT
  downloaded from Ensembl; the equivalent local PGSB file is used instead.

## Method

1. `00_lock_fdr_genes.R` - select `significant_fdr == TRUE`; write the unique gene
   list and the (trait, gene) table with `legacy_annotation`.
2. `01_make_protein_fasta.R` - copy the HC proteome, subset to the locked genes, one
   representative peptide per gene -> `inputs/fdr_proteins.faa` (45 sequences).
3. `02_swissprot_diamond.sh` - download Swiss-Prot, `diamond makedb`, then
   `diamond blastp --very-sensitive -e 1e-5 --max-target-seqs 5` with outfmt 6
   `qseqid sseqid stitle pident qcovhsp length evalue bitscore`. `stitle` carries
   the protein name, `OS=` organism and `OX=` taxid.
4. `03_build_annotation_table.R` - best hit per gene (lowest e-value, ties by
   highest bitscore), parse name/accession/organism, classify, build deliverables.

Reproduce all: `bash scripts/run_all.sh` (idempotent; logs to `logs/`).
Temp is pinned to `/mnt/data/shahar/.tmp` (TMPDIR) in every script.

## Classification thresholds (explicit, adjustable)

A gene is **CONFIDENT_SWISSPROT_HIT** when its best hit satisfies ALL of:
- `evalue <= 1e-5`
- `qcovhsp >= 50` (percent of the query covered)
- `pident >= 30`
- `characterized_flag == TRUE`

otherwise **WEAK_OR_NO_HIT** (includes genes with no Swiss-Prot hit at all).

`characterized_flag` = best-hit title does NOT match (case-insensitive)
`predicted|uncharacterized|hypothetical|DUF[0-9]|unknown function|putative uncharacterized`.

To change thresholds, edit the `TH_*` / `UNCHAR_RE` constants at the top of
`03_build_annotation_table.R` and re-run step 03 only.

## Results (this run)

- **23 unique genes; 23 (trait, gene) rows.**
- **17 of 23** genes had >=1 Swiss-Prot hit; 6 had none.
- **16 CONFIDENT_SWISSPROT_HIT**, 7 WEAK_OR_NO_HIT (the residual list handed to step 06/07).
- Unique-gene status: **29 CONFIDENT_SWISSPROT_HIT, 16 WEAK_OR_NO_HIT**.
- Per (trait, gene) rows CONFIDENT / WEAK: beta-glucan 18/12, fiber 6/3, starch 8/1.

## Deliverables

- `results/tables/fdr_gene_annotation_swissprot.tsv` - one row per (trait, gene):
  `trait, gene_id, lead_SNP, sp_name, sp_accession, sp_organism, pident, qcovhsp,
  evalue, bitscore, characterized_flag, status, legacy_annotation`.
- `results/tables/residual_weak_or_no_hit.txt` - the 16 WEAK_OR_NO_HIT genes; these
  feed the next (online) phase.
- `inputs/fdr_proteins.faa` - the 45 query proteins.

## What is NOT here (the online phase - do not run on this server)

- NCBI nr BLASTP fallback (Viridiplantae) for the residual WEAK_OR_NO_HIT genes.
- InterProScan (Pfam/InterPro domains + GO) on `fdr_proteins.faa`.
- Per-gene reconciliation (Swiss-Prot + nr + InterPro -> one call with confidence),
  comparison vs `legacy_annotation` and the grain-quality candidate-gene literature.
- Biological-relevance filter and GO enrichment (deferred until the set is final).

## Directory layout

```
05_USED_gene_annotation/
  scripts/        00_lock_fdr_genes.R  01_make_protein_fasta.R  02_swissprot_diamond.sh
                  03_build_annotation_table.R  run_all.sh
  inputs/         fdr_genes.txt  fdr_genes_table.tsv
                  Hv_Morex.pgsb.Jul2020.HC.aa.fa (provenance copy)  fdr_proteins.faa
  intermediates/  uniprot_sprot.fasta.gz  swissprot.dmnd
                  fdr_vs_swissprot.diamond.tsv  versions.txt
  logs/           per-script logs
  results/tables/ fdr_gene_annotation_swissprot.tsv  residual_weak_or_no_hit.txt
  README.md
```


## Set-size independence (2026-09-10)

Two v1 hard-codings were removed so this step survives a change of gene set:

- `01_make_protein_fasta.R` asserted a fixed count of 45 and a specific gene's length
  (`1HG0079280 == 206 aa`). It now checks the invariants that hold for **any** set -
  one representative per requested gene, no duplicates, all non-empty - and only runs
  the length cross-check when that reference gene happens to be present.
- `02_swissprot_diamond.sh` reported hits "of 45"; it now counts the query FASTA.
