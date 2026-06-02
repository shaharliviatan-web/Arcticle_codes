# 08_USED_annotation_master

Consolidated functional-annotation master table for the 45 FDR-passing candidate
genes of the wild-barley grain-quality GWAS. Merges the three evidence sources
gathered earlier into ONE auditable table, one row per (trait, gene), assigns a
single best functional call + source, an evidence tier, and a provisional
Swiss-Prot-vs-InterPro agreement flag. This is the paper's supplementary
annotation table and the input to the (manual) biological-relevance filter.

NOTHING is filtered or dropped here - all 48 (trait, gene) rows are kept,
including both trait rows for the 3 fiber/starch shared genes.

## Inputs

- `../05_USED_gene_annotation/results/tables/fdr_gene_annotation_swissprot.tsv`
  (Swiss-Prot BLASTP; 48 rows, per trait+gene)
- `../06_USED_interpro_domains/results/tables/fdr_gene_interpro.tsv`
  (InterProScan domains + GO; 45 rows, per gene)
- `../07_USED_BLASTP_genes_with_no_annotation_left/results/tables/fdr_residual_nr.tsv`
  (nr BLASTP for the 5 residual genes)
- `../04_USED_haplotype_analysis_crosshap/05_results/tables/gene_annotation_review.csv`
  - for the lead-SNP `class` (significant / marginal), joined on (trait, gene_id).

06 and 07 are per gene; their evidence is replicated onto each of a gene's trait
rows. The 3 shared genes therefore carry identical Swiss-Prot/InterPro/nr evidence
and the same call/tier on both trait rows, differing only in `trait`, `lead_SNP`
and `lead_SNP_class`.

## Method (`scripts/00_build_master.R`, `bash scripts/run_all.sh`)

Pure local merge, no network. TMPDIR pinned to `/mnt/data/shahar/.tmp`; absolute
`/usr/bin/Rscript`; no `/tmp` or `$HOME` writes; output pure ASCII.

### final_call + call_source - reconciliation priority

1. confident Swiss-Prot hit  -> `swissprot`, name = Swiss-Prot protein name
2. else InterPro domain match -> `interpro`, name = InterPro entry description(s)
3. else confident nr hit      -> `nr`, name = nr title
4. else                       -> `none`, name = "uncharacterized (conserved,
   unnamed)" if any homolog/signature exists, else "no homolog found"

"confident Swiss-Prot" = step-05 `CONFIDENT_SWISSPROT_HIT`; "InterPro domain
match" = step-06 `INTERPRO_MATCH` (an integrated InterPro entry, not merely a
disorder/TM/coiled-coil signature); "confident nr" = step-07 `CONFIDENT_NR_HIT`.

### evidence_tier

- `HIGH`   - confident Swiss-Prot AND a concordant InterPro domain AND Swiss-Prot
  identity >= 40% (i.e. not near the 30% floor).
- `MEDIUM` - a single confident source, OR a would-be HIGH whose Swiss-Prot
  identity is near-floor (`< 40%`), OR confident Swiss-Prot + InterPro that are
  discordant.
- `LOW`    - only uncharacterized hits, or nothing.

`NEAR_FLOOR_PIDENT = 40` is set at the top of the script and is adjustable.

### sp_vs_interpro_agreement (provisional - for manual review, not auto-resolved)

Best-effort keyword overlap between the Swiss-Prot name and the InterPro/Pfam
descriptions:
- `concordant`    - a meaningful keyword is shared.
- `discordant`    - both sides have a descriptor but no keyword overlaps.
- `single_source` - only one of the two sides has a usable descriptor.
- `none`          - neither side has a usable descriptor.

The comparison uses ANY characterized Swiss-Prot name (including names that failed
the step-05 coverage/identity cutoff, e.g. ERF110), so weak-but-named hits are
still evaluated against the domain.

The matcher is deliberately CONSERVATIVE (exact token match, or substring with
length >= 4; generic words like protein/domain/family/putative/unknown are
ignored). It is biased to NOT over-call concordance: several `discordant` rows are
in fact biologically concordant once a human reads them, e.g.
- `3HG0299460` Polypyrimidine tract-binding protein vs "RNA recognition motif"
  (PTB proteins contain RRMs),
- `3HG0299450` "Protein LIKE COV 2" vs "CONTINUOUS VASCULAR..." (COV = that name),
- the three `2HG0191510/520/540` "UPF0481" vs DUF247 "unknown function" families.
Treat every `discordant` as a flag to inspect, not as a contradiction.

## Results (this run)

48 rows, 45 unique genes (none dropped).

- call_source: swissprot 32 rows (29 genes), interpro 11, nr 1, none 4.
- evidence_tier: HIGH 24, MEDIUM 20, LOW 4.
- sp_vs_interpro_agreement: concordant 31, discordant 6, single_source 6, none 5.
- lead_SNP_class: significant 30, marginal 18.
- The 4 LOW/none genes: `1HG0079290`, `2HG0112840`, `4HG0339140`
  (conserved but uncharacterized) and `2HG0112890` (no homolog found).

## Deliverable

`results/tables/fdr_annotation_master.tsv` - 48 rows x 31 columns:

`trait, gene_id, lead_SNP, lead_SNP_class, final_call, call_source,
evidence_tier, sp_vs_interpro_agreement, sp_name, sp_accession, sp_organism,
sp_pident, sp_qcovhsp, sp_evalue, sp_status, interpro_status, n_interpro,
interpro_ids, interpro_descs, n_pfam, pfam_ids, pfam_descs, n_go, go_terms,
nr_title, nr_accession, nr_pident, nr_qcovs, nr_evalue, nr_status,
legacy_annotation`

All raw per-source evidence is carried so every call is auditable. `legacy_annotation`
is the old manual note, retained for comparison only - never used in the call.

## Scope

Build the table only. The biological-relevance filter (keeping genes whose function
plausibly connects to beta-glucan / fiber / starch pathways) and any manual
resolution of `discordant` rows are the NEXT step and are NOT done here.

## Directory layout

```
08_USED_annotation_master/
  scripts/        00_build_master.R  run_all.sh
  logs/           00_build_master.log
  results/tables/ fdr_annotation_master.tsv
  README.md
```
