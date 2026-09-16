# 07_USED_BLASTP_genes_with_no_annotation_left

Final rescue step in the functional-annotation pipeline. After Swiss-Prot BLASTP
(step 05) and InterProScan (step 06), any gene still without a functional call is
BLASTed against NCBI nr restricted to Viridiplantae, to find a named plant homolog
the earlier sources lacked - or to confirm it is genuinely uncharacterized.

**Re-run 2026-09-10** for step-04 run `loci_LDspan_eps06_V4`: **2 residual genes**.

## The residual set is DERIVED, not hard-coded (changed 2026-09-10)

v1 listed 5 gene IDs literally in the script, so it aborted the moment the gene set
changed. `00_extract_residual_proteins.sh` now computes:

```
residual = {Swiss-Prot WEAK_OR_NO_HIT}  INTERSECT  {no InterPro match}
```

A gene rescued by *either* source is not residual - taking the intersection rather
than the union is the whole point of running two sources. If the intersection is
empty the step exits cleanly with nothing to do, which is a valid outcome.

Files were renamed `residual_5*` -> `residual*` to match.

### Current run: 2 residual genes

Swiss-Prot left 7 without a confident call, InterPro left 2 without a domain match;
the intersection is 2:

| gene_id | trait | rank | step-06 status |
|---|---|---|---|
| HORVU.MOREX.r3.6HG0619810 | starch | **#1** | SIGNATURE_ONLY |
| HORVU.MOREX.r3.4HG0339140 | fiber | #10 | SIGNATURE_ONLY |

## Method

NCBI nr is not installed locally, so BLASTP is run REMOTELY on NCBI servers.

1. `00_extract_residual_proteins.sh` - derive the residual set and extract those
   proteins into `inputs/residual.faa` (verifies the count matches, no stop codons).
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
  inputs/         residual.faa  residual_ids.txt
  intermediates/  residual_vs_nr.tsv  residual_vs_nr.done
  logs/           versions.txt  01_progress.log  01_blastp.log  *.log
  results/tables/ fdr_residual_nr.tsv  residual_nr_note.txt
  README.md
```


## Result of the current run

Both genes returned strong nr hits, but **neither is a functional name**:

| gene | best nr hit | identity | verdict |
|---|---|---|---|
| 6HG0619810 | *hypothetical protein* ZWY2020_045447 (*H. vulgare*) | 85% | uncharacterized |
| 4HG0339140 | *uncharacterized protein* LOC123451082 (*H. vulgare*) | **100%** | uncharacterized |

The 100% hit is the same barley gene in another assembly: the gene is real and
conserved, but nobody has characterized it. All three evidence lines therefore agree
that these two are genuinely uncharacterized (`annotation_source = none` in step 08).

This matters for interpretation, because `6HG0619810` is the **most significant gene
in the whole study** (BH q = 5.5e-11, eta-squared = 0.261).

Timing note: the remote nr search took ~2 h of NCBI queue for 2 proteins. It is
resumable - `intermediates/residual_vs_nr.done` makes a re-run skip a completed search.
