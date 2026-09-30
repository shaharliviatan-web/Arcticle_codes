# Loci table for the supplement — note (2026-09-16)

**Decision:** the 36-locus table goes to the supplement as an Online Resource, cited from Results
chapter 2. **Not built yet** — this note only records that it is needed and where the source is.

## The source table already exists

`01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv`
— 36 rows × 20 columns, one row per locus:

`locus_id, trait, chr, lead_SNP, lead_bp, locus_start, locus_end, span_kb, n_member_SNPs,`
`lead_A1, lead_A2, lead_MAF, lead_beta, lead_SE, lead_p, lead_neg_log10_p,`
`max_internal_gap_kb, clump_span_before_severing_kb, n_clump_members_before_severing, n_members_severed`

Companion tables in the same directory:

| file | one row per | use |
|---|---|---|
| `Table_loci_per_trait.tsv` | trait | the small in-text summary table |
| `Table_loci_members_full.tsv` | member SNP | too large for the body — supplement only |
| `Table_analysis_parameters.tsv` | parameter | Materials and methods |
| `Table_loci_for_gene_search.tsv` | locus | the step-03 handoff; not for the paper |

## What still has to be decided / done

- [ ] Which columns the published table keeps. The last four (`max_internal_gap_kb`,
      `clump_span_before_severing_kb`, `n_clump_members_before_severing`, `n_members_severed`)
      are clumping diagnostics — probably supplement-only, possibly dropped.
- [ ] Add a **carrier count** column next to `lead_MAF` (number of accessions carrying the minor
      allele, out of the accessions with a call at that SNP). Not a new analysis — the same MAF
      restated as a count, from `plink --freq counts` (run 2026-09-16); not yet in the table.
      More intuitive than MAF alone. **Keep-or-drop still open.**
- [ ] Decide whether the 55 candidate genes
      (`03_USED_candidate_genes_around_leading_snps/results/tables/candidate_genes.tsv`)
      are a second Online Resource or extra columns on the locus table.
- [ ] Export to `new_publishing_paper/supplementary/ESM_N.xlsx` once the ESM numbering is fixed.

## Why not in the body

36 rows × ~12 useful columns does not fit a printed TAG table. The body gets the per-trait summary
(`Table_loci_per_trait.tsv`), the supplement gets the full locus list.

**Nothing in the analysis directories was modified.** The carrier-count column above was derived in a
scratch directory for review only; if it is adopted it must be added by a script in step 01 and
documented there.

**Dropped 2026-09-17:** a lead-SNP-to-nearest-gene distance column was considered and rejected — not
a necessary analysis. Removed from this note and from the blueprint.
