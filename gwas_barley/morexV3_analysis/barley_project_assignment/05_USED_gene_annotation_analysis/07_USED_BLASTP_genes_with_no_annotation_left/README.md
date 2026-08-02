# 07_USED_BLASTP_genes_with_no_annotation_left

Final rescue step in the functional-annotation pipeline for the wild-barley
grain-quality GWAS FDR genes. After Swiss-Prot BLASTP (step 05) and InterProScan
(step 06), 5 genes still had no functional call. This step runs a last-resort
BLASTP of just those 5 proteins against NCBI nr restricted to Viridiplantae, to
find a named plant homolog that Swiss-Prot lacked, or to confirm they are
genuinely uncharacterized.

## The 5 residual genes (input)

Pulled from `../05_USED_gene_annotation/inputs/fdr_proteins.faa` into
`inputs/residual_5.faa`:

| gene_id | step-06 status | aa |
|---|---|---|
| HORVU.MOREX.r3.1HG0079220 | SIGNATURE_ONLY (disorder/TM/coiled-coil only) | 204 |
| HORVU.MOREX.r3.1HG0079290 | SIGNATURE_ONLY | 132 |
| HORVU.MOREX.r3.2HG0112840 | SIGNATURE_ONLY | 248 |
| HORVU.MOREX.r3.4HG0339140 | SIGNATURE_ONLY | 211 |
| HORVU.MOREX.r3.2HG0112890 | NO_MATCH | 78 |

## Method

NCBI nr is not installed locally, so BLASTP is run REMOTELY on NCBI servers.

1. `00_extract_residual_proteins.sh` - extract the 5 proteins into
   `inputs/residual_5.faa` (verifies 5 sequences, no stop codons).
2. `01_blastp_nr_remote.sh` - one batched submission of the 5 proteins:
   `blastp -remote -db nr -entrez_query "Viridiplantae[ORGN]" -evalue 1e-5
   -max_target_seqs 5 -outfmt 6 qseqid sseqid stitle pident qcovs evalue bitscore
   staxids sscinames`. Idempotent (a `.done` marker skips re-running); retries up
   to 3x on transient NCBI errors. `sscinames`/`staxids` come back as `N/A`
   because no local taxdb is installed - this is expected and handled in step 02.
3. `02_build_table.R` - best hit per gene (lowest e-value, ties by bitscore),
   organism parsed from the nr title's trailing `[Organism]`, characterized flag
   via the SAME regex as step 05.

Reproduce: `bash scripts/run_all.sh` (launch step 01 under GNU screen for a long
search - see `run_all.sh`). TMPDIR pinned to `/mnt/data/shahar/.tmp`; absolute
tool paths (`/usr/bin/blastp`); no `/tmp` or `$HOME` writes.

## Versions

See `logs/versions.txt`. This run: `blastp 2.12.0+`, NCBI nr (remote),
`Viridiplantae[ORGN]`, e-value <= 1e-5, top 5 hits, accessed 2026-06-02 (UTC).
The remote search uses NCBI's then-current nr; record the access date for Methods.

## Classification

`characterized_flag` = best-hit title does NOT match (case-insensitive)
`predicted|uncharacterized|hypothetical|DUF[0-9]|unknown function|putative uncharacterized`
(identical to step 05).

`status`:
- `CONFIDENT_NR_HIT` - e-value <= 1e-5 AND qcovs >= 50 AND pident >= 30 AND
  characterized_flag == TRUE.
- `WEAK_OR_UNCHARACTERIZED_NR_HIT` - has an nr hit but it is uncharacterized or
  fails a threshold.
- `NO_NR_HIT` - no Viridiplantae nr hit at e-value <= 1e-5.

## Results (this run)

| gene_id | best nr hit | organism | pident | qcovs | e-value | status |
|---|---|---|---|---|---|---|
| 1HG0079220 | antifreeze protein Maxi-like | H. vulgare subsp. vulgare | 100 | 100 | 1.97e-132 | CONFIDENT_NR_HIT |
| 1HG0079290 | uncharacterized protein LOC123416960 | H. vulgare subsp. vulgare | 100 | 100 | 5.02e-90 | WEAK_OR_UNCHARACTERIZED_NR_HIT |
| 2HG0112840 | uncharacterized protein LOC123424732 | H. vulgare subsp. vulgare | 100 | 100 | 9.68e-174 | WEAK_OR_UNCHARACTERIZED_NR_HIT |
| 4HG0339140 | uncharacterized protein LOC123451082 | H. vulgare subsp. vulgare | 100 | 100 | 3.19e-153 | WEAK_OR_UNCHARACTERIZED_NR_HIT |
| 2HG0112890 | (no hit) | - | - | - | - | NO_NR_HIT |

**1 of 5 got a real plant name** (`1HG0079220` -> "antifreeze protein Maxi-like");
3 stay uncharacterized; 1 (`2HG0112890`, 78 aa) has no Viridiplantae nr hit.

Caveat for interpretation: every top hit is 100% identical to *Hordeum vulgare*
subsp. *vulgare* (cultivated barley) - i.e. nr mostly just re-finds the same gene
model in cultivated barley, where it is also annotated "uncharacterized". So nr
adds little functional information here beyond confirming these are conserved
barley proteins; only `1HG0079220` carries a descriptive RefSeq name.

## Deliverables

- `results/tables/fdr_residual_nr.tsv` - one row per residual gene (5 rows):
  `gene_id, nr_title, nr_accession, nr_organism, pident, qcovs, evalue, bitscore,
  characterized_flag, status`.
- `results/tables/residual_nr_note.txt` - one-line named-vs-uncharacterized summary.

## Scope

Just these 5 genes vs nr (Viridiplantae). NO cross-evidence reconciliation and NO
biological-relevance filter yet - those are the following steps.

## Directory layout

```
07_USED_BLASTP_genes_with_no_annotation_left/
  scripts/        00_extract_residual_proteins.sh  01_blastp_nr_remote.sh
                  02_build_table.R  run_all.sh
  inputs/         residual_5.faa  residual_5_ids.txt
  intermediates/  residual_5_vs_nr.tsv  residual_5_vs_nr.done
  logs/           versions.txt  01_progress.log  01_blastp.log  *.log
  results/tables/ fdr_residual_nr.tsv  residual_nr_note.txt
  README.md
```
