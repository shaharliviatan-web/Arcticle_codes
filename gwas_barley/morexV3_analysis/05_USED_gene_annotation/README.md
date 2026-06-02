# 05_USED_gene_annotation

Reproducible functional annotation of the FDR-passing candidate genes from the
wild-barley grain-quality GWAS (Morex V3 reference, PGSB r3 annotation). This is a
NEW top-level pipeline step; it does not reuse code from steps 01-04.

Primary evidence here = **DIAMOND BLASTP vs UniProt Swiss-Prot** (curated protein
names). GO / InterPro / NCBI nr are the NEXT, online phase and are NOT done here.

## Input set

FDR-passing genes = rows with `significant_fdr == TRUE` in the step-04 review table:
`../04_USED_haplotype_analysis_crosshap/05_results/tables/gene_annotation_review.csv`

- 48 (trait, gene) rows = 45 unique genes: beta-glucan 30, fiber 9, starch 9
  (protein has 0 FDR genes).
- 3 genes are shared between fiber and starch (overlapping +/-200 kb loci):
  `HORVU.MOREX.r3.7HG0729020`, `...7HG0729090`, `...7HG0729100`. The protein is
  annotated once per unique gene and the annotation is replicated to each trait
  row; only the per-trait `lead_SNP` differs.

The legacy `annotation` column from that CSV is carried as `legacy_annotation` for
comparison ONLY. It is not reproducible (manual Ensembl-website lookups + BLASTX on
CDS, no recorded evidence) and is never used as the functional call.

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
    isoform index (= canonical `.1` for single-isoform genes). Only one of the 45
    genes is multi-isoform: `3HG0299460` -> `.2` chosen (440 aa vs 331 aa).
- Swiss-Prot: UniProtKB/Swiss-Prot Release **2026_01** of 28-Jan-2026,
  `uniprot_sprot.fasta.gz` from
  `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/`
  (574,627 sequences). See `intermediates/versions.txt`.
- DIAMOND **2.0.14** (`/usr/bin/diamond`). makeblastdb/blastp also present at
  `/usr/bin` but unused here.
- Ensembl Plants FTP is unreachable from this server, so the peptide file was NOT
  downloaded from Ensembl; the equivalent local PGSB file is used instead.

## Method

1. `00_lock_fdr_genes.R` - select `significant_fdr == TRUE`; write the unique gene
   list and the (trait, gene) table with `legacy_annotation`.
2. `01_make_protein_fasta.R` - copy the HC proteome, subset to the 45 genes, one
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

- 45 unique genes; 48 (trait, gene) rows.
- 36 of 45 genes had >=1 Swiss-Prot hit; 9 had none.
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
