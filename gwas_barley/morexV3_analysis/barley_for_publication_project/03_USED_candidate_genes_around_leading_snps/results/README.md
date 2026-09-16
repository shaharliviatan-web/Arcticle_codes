# results/

All output of step 03. Tables only — no plots. Everything here is regenerated from
scratch by `bash ../scripts/run_all.sh`; nothing is hand-edited.

**For writing the paper, start with [`tables/results_chapter_numbers.txt`](tables/results_chapter_numbers.txt)** —
every number in this step as prose, including a paste-ready sentence, a per-trait
block, a per-locus listing with each gene and its distance to the lead SNP, and the
list of loci that came back empty.

## tables/

| file | one row per | what it is for |
|---|---|---|
| **`candidate_genes.tsv`** | gene × locus | **the main result.** Also the step-04 input — column names are a contract, see the top-level README |
| **`lead_loci.tsv`** | locus | every locus with its coordinates, lead SNP and gene counts. Also step-04 input |
| `loci_without_genes.tsv` | locus | the 21 loci (of 36) whose interval contains no annotated gene; 4 are lead-only (span 0 -> 1 bp interval) |
| `genes_with_annotation.tsv` | gene | the 3 candidates carrying a functional description — sorted by distance to lead. **This is the shortlist you actually write about** |
| `per_trait_summary.tsv` | trait | counts per trait — the summary table for the results chapter |
| `analysis_parameters.tsv` | parameter | **every constant, step 01 + step 03 in one file** — the methods section |
| `flank_sensitivity.tsv` | flank size | what 0/25/50/100/200 kb would each return; `is_chosen` flags the live setting |
| `flank_sensitivity_by_trait.tsv` | flank × trait | the same, split by trait |
| `publication_gene_recovery.tsv` | curated gene | whether step 04's four publication genes are recovered, and the flank each would need |
| `results_chapter_numbers.txt` | — | all of the above as prose |

## Column notes

**`dist_to_lead_bp`** — signed nearest-edge distance from the lead SNP to the gene:
`0` = the lead SNP falls inside the gene body; `<0` = gene upstream (lower
coordinate); `>0` = gene downstream. Sort by `abs()` for "closest gene to the peak".

**`class`** — `significant_locus` for all 36. The 2 sub-threshold protein peaks that
earlier versions carried were **dropped by decision on 2026-09-09**; protein now contributes
only its genome-wide-significant loci. `INCLUDE_SUBTHRESHOLD_PEAKS` in `config/params.sh`
therefore has no effect on the current input. The two
sub-threshold protein peaks sit just under the Bonferroni lead threshold; step 01
kept them separate rather than lowering `--clump-p1`. They are included here because
the handoff table flags `include_in_gene_search=TRUE`. To exclude them, set
`INCLUDE_SUBTHRESHOLD_PEAKS=no` in `config/params.sh`.

**`in_locus_span`** — TRUE for every row while `FLANK_BP=0`, because the interval
*is* the span. The column exists so it stays meaningful if a flank is switched on.

**`has_annotation` / `description`** — only ~6.8% of Morex V3 genes carry a
`description=` attribute (projected from Arabidopsis/rice/UniProt). There are **no
GO terms in this GFF**. Most candidates therefore have a gene ID and nothing else;
deeper functional annotation is step 05's job, not this one.

**`locus_id`** — assigned by step 01 and **not stable across its re-runs**; join on
`lead_SNP` when comparing to an older result set.

**`n_member_SNPs`** — carried through from step 01: how many SNPs the LD clump held
after gap-severing. A locus with 1 member is a lone significant SNP.

## Caveat on the current numbers

`FLANK_BP=0` leaves 21 of 36 loci empty and misses 3 of the 4 curated publication genes. This is a deliberate setting, not a bug — see the warning section in the
top-level [`../README.md`](../README.md) and `flank_sensitivity.tsv` for the
alternatives.
