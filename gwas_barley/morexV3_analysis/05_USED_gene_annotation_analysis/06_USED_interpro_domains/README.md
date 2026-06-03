# 06_USED_interpro_domains

Protein-domain and Gene Ontology (GO) annotation of the 45 FDR-passing candidate
genes from the wild-barley grain-quality GWAS, via EBI InterProScan 5. This is the
second line of functional evidence after the Swiss-Prot BLASTP in
`../05_USED_gene_annotation/`: domains validate the confident Swiss-Prot hits,
can rescue the no-hit genes, and supply GO terms the Morex V3 GFF3 lacks.

Run on ALL 45 FDR genes (not only the residual 16), so every gene has both lines
of evidence for the later reconciliation step.

## Input

`../05_USED_gene_annotation/inputs/fdr_proteins.faa` (45 representative proteins,
one per FDR gene; copied here to `inputs/fdr_proteins.faa` for provenance). Headers
are `HORVU.MOREX.r3.<gene>.<isoform>`; the per-gene id is the header minus the
trailing `.<isoform>`.

## Method

InterProScan is run remotely on the EBI servers through the official EBI Web
Services REST client (`iprscan5`); InterProScan itself is NOT installed locally.

1. `00_setup_and_versions.sh` - ensure the official `iprscan5.py` client and its
   `xmltramp2` dependency are present (both vendored under `scripts/`, nothing
   under `$HOME`), apply a small bytes->str patch to the client (see below), and
   record tool/service versions to `logs/versions.txt`.
2. `01_run_interproscan.sh` - split-free per-sequence loop: each of the 45 proteins
   is submitted as its own InterProScan job with the required `--email`, polite
   status polling (`--pollFreq 5`, 3 s between jobs), `--goterms` and `--pathways`,
   the default application set (includes Pfam; no `--appl` restriction), and
   `--stype p`. Results retrieved per sequence as TSV + JSON + GFF3 into
   `intermediates/raw/<seqid>.{tsv.tsv,json.json,gff.txt}`. Idempotent / resumable:
   a sequence whose valid JSON already exists is skipped, so re-running only fills
   gaps or retries failures.
3. `02_parse_interpro.R` - parse the per-sequence TSVs into one row per gene.

Reproduce: `bash scripts/run_all.sh`. For the ~30 min remote run, launch step 01
under GNU screen (see `run_all.sh` header). TMPDIR is pinned to
`/mnt/data/shahar/.tmp` and absolute paths are used throughout.

### Client patch (recorded for transparency)

The upstream `iprscan5.py` (rev 2024-07-11) crashes intermittently during status
polling: its `restRequest()` HTTPError fallback returns `requests`' `.content`
(bytes) and the client then does `'status: ' + status` (str + bytes -> TypeError).
`00_setup_and_versions.sh` applies an idempotent patch decoding that fallback to
str (binary bodies are left as bytes), matching the non-error path. No behavioural
change otherwise.

## Versions

See `logs/versions.txt`. This run: InterProScan **5.77-108.0** (EBI),
iprscan5 client rev **2024-07-11**, EBI REST
`https://www.ebi.ac.uk/Tools/services/rest/iprscan5`. `--goterms` and `--pathways`
enabled; default member-database set (Pfam, PANTHER, Gene3D, SUPERFAMILY, CDD,
SMART, PRINTS, ProSiteProfiles, NCBIfam, FunFam, etc.).

## Deliverable

`results/tables/fdr_gene_interpro.tsv` - one row per gene (45 rows). Columns:

| column | meaning |
|---|---|
| `gene_id` | HORVU.MOREX.r3 gene id |
| `protein_length` | aa length of the representative peptide |
| `n_signatures` | total member-database matches |
| `member_dbs` | distinct member databases that matched (`;`-separated) |
| `n_interpro` | number of distinct InterPro entries |
| `interpro_ids` | distinct InterPro accessions (`;`) |
| `interpro_descs` | `IPRxxxxxx:description` (`; `) |
| `n_pfam` | number of distinct Pfam entries |
| `pfam_ids` | distinct Pfam accessions (`;`) |
| `pfam_descs` | `PFxxxxx:description` (`; `) |
| `n_go` | number of distinct GO terms |
| `go_terms` | distinct GO accessions, e.g. `GO:0003924;GO:0005525` (`;`) |
| `status` | `INTERPRO_MATCH` / `SIGNATURE_ONLY` / `NO_MATCH` |

GO term *names* and pathway annotations are not in the TSV columns but are present
in the per-gene JSON (`intermediates/raw/<seqid>.json.json`) if needed later.

## Results (this run)

45/45 genes annotated successfully (0 job failures after the client patch).

- Per-gene status: **40 INTERPRO_MATCH, 4 SIGNATURE_ONLY, 1 NO_MATCH**.
- 40 genes have >=1 InterPro entry; 39 have >=1 Pfam; 35 have >=1 GO term.
- `SIGNATURE_ONLY` (member-db hit but no InterPro entry, i.e. only disorder /
  transmembrane / signal-peptide / coiled-coil / unintegrated PANTHER):
  `1HG0079220`, `1HG0079290`, `2HG0112840`, `4HG0339140`.
- `NO_MATCH`: `2HG0112890` (78 aa - too short for any domain signature).

Rescue of the Swiss-Prot residual set (the 16 WEAK_OR_NO_HIT genes from step 05):
**11 of 16 rescued** with an InterPro domain. Notable agreements/rescues:
- `3HG0299440` -> AP2/ERF domain (PF00847 / IPR001471) - confirms the low-coverage
  Swiss-Prot ERF110 hit that had failed the coverage threshold.
- `7HG0642350` -> Transferase family PF02458 - consistent with a BAHD acyltransferase.
- `2HG0191510/520/540` -> PF03140 (DUF247) - the three former "UPF0481" weak hits
  resolve to the DUF247 plant family.
- `6HG0627120` -> PF03109 (ABC1 family); `2HG0112900` -> PF00578 (AhpC/TSA,
  peroxiredoxin-like); `1HG0079230` -> PF01984.
Not rescued: the 4 SIGNATURE_ONLY genes above plus `2HG0112890` (NO_MATCH).

## Scope

InterProScan only. NCBI nr BLAST, cross-evidence reconciliation (Swiss-Prot + nr +
InterPro -> one functional call), the biological-relevance filter, and GO
enrichment are later steps and are NOT done here.

## Directory layout

```
06_USED_interpro_domains/
  scripts/        00_setup_and_versions.sh  01_run_interproscan.sh
                  02_parse_interpro.R  run_all.sh  iprscan5.py (vendored client)
                  pylib/ (vendored xmltramp2 + six)
  inputs/         fdr_proteins.faa (45, copy of the step-05 input)
  intermediates/  seqs/<seqid>.fasta (45)   raw/<seqid>.{tsv.tsv,json.json,gff.txt}
  logs/           versions.txt  01_progress.log  <seqid>.iprscan.log
  results/tables/ fdr_gene_interpro.tsv
  README.md
```
