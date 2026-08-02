# 03_USED_candidate_genes_around_leading_snps

Candidate genes around the GWAS lead SNPs for the wild-barley grain-quality TAG paper.
290 wild barley accessions (*Hordeum vulgare* ssp. *spontaneum*, Southern Levant), Morex **V3**
reference, 7,110,996 SNPs, 4 traits (beta-glucan, fiber, protein, starch). GWAS: EMMAX (aIBS
kinship), **BLUP phenotypes + 3 PCs**; significance = Bonferroni alpha=0.10 on the LD-pruned
590,462-SNP set (-log10p = 6.7712). Source GWAS:
`USED_GWAS_V2_pipeline/results/publication_BonfOnly_BLUP_3PC/`.

This step extracts candidate genes from the Morex V3 GFF3 within an LD-based window around the
lead SNPs, organised per trait / per locus. **Tables only, no plots.** The crosshap haplotype
analysis is the next step (not done here).

## Window and locus rule (locked)

- **Window = 200 kb** each side. The genome-wide LD decay reaches r^2=0.2 at ~188 kb in
  `USED_LD_decay_V2_wholegenome`, rounded to a clean 200 kb.
- **Loci by single-linkage clustering**, per trait, per chromosome: SNPs whose nearest neighbour
  is <= 200 kb away join one locus (chains of 3-4 SNPs generalise). Lead SNP = highest -log10p in
  the locus.
- **Candidate interval = cluster span +/- 200 kb** (leftmost member -200 kb .. rightmost member
  +200 kb). For a lone SNP this equals lead +/- 200 kb.

Methods phrasing: *"significant and marginal SNPs were grouped by single-linkage clustering at
200 kb (the genome-wide LD-decay distance to r^2=0.2); each locus interval spanned its clustered
SNPs extended by 200 kb on each side; the most significant SNP was the lead."*

## Inputs

- **Significant SNPs (15):** all `above_Bonf=TRUE` rows from
  `USED_GWAS_V2_pipeline/results/publication_BonfOnly_BLUP_3PC/tables/lead_snps.tsv`
  (beta-glucan 12, fiber 2, protein 1, starch 0). Pulled by `01_build_loci.R`.
- **Marginal SNPs (9, curated):** `inputs/marginal_snps.tsv` (BLUP x 3PC -log10p values, below
  6.7712 but near it). beta-glucan 3, fiber 2, protein 0, starch 4.
- **GFF3:** `../Hordeum_vulgare.MorexV3_pseudomolecules_assembly.62.gff3.gz`. Genes = col3 ==
  `gene`. SNP/GFF chromosome names both use `1H..7H` (no `chr` prefix). Unplaced `CAJHDD*`
  scaffolds are excluded. All 35,106 genes on 1H-7H are `biotype=protein_coding`; ~6.8% carry a
  `description=` (functional annotation projected from Arabidopsis/rice/UniProt; no GO). The
  `description=` value is percent-decoded (e.g. `%3B` -> `;`) and left NA when absent.

## Directory layout

```
inputs/         marginal_snps.tsv (curated); curated_snps.tsv (generated: sig + marginal, class-tagged)
scripts/        01_build_loci.R  02_extract_genes.sh  03_build_tables.R  run_all.sh
intermediates/  loci.tsv  loci_intervals.bed  genes_1to7H.gff  intersect_raw.tsv
logs/           per-step logs
results/tables/ candidate_genes.tsv  lead_loci.tsv  per_trait_summary.tsv  results_chapter_numbers.txt
```

## Pipeline (`bash scripts/run_all.sh`)

1. **01_build_loci.R** — read significant (`lead_snps.tsv`) + marginal (`marginal_snps.tsv`) SNPs,
   write the combined `inputs/curated_snps.tsv`, single-linkage cluster (200 kb) per trait/chr,
   emit `intermediates/loci.tsv` and `intermediates/loci_intervals.bed` (intervals = span +/- 200 kb).
2. **02_extract_genes.sh** — filter the GFF to `gene` features on 1H-7H
   (`intermediates/genes_1to7H.gff`), `bedtools intersect -wa -wb` against the interval BED
   (`intermediates/intersect_raw.tsv`).
3. **03_build_tables.R** — parse gene attributes (gene_id / biotype / description), compute the
   signed gene->lead distance, assemble the four output tables.

`loci_intervals.bed` is 0-based half-open (interval_start - 1 .. interval_end); the GFF is 1-based.
`bedtools` converts internally from the `.gff` extension, so the intersection is coordinate-correct.

## Output tables (`results/tables/`)

- **candidate_genes.tsv** (main, one row per gene x locus): trait, locus_id, lead_SNP, class,
  gene_id, biotype, chr, gene_start, gene_end, strand, dist_to_lead_bp (signed), in_cluster_span,
  description.
- **lead_loci.tsv**: trait, locus_id, chr, lead_SNP, lead_pos, lead_neg_log10p, class, n_SNPs,
  member_SNPs, span_start, span_end, interval_start, interval_end, n_genes, n_protein_coding,
  n_with_description.
- **per_trait_summary.tsv**: trait, n_loci (+ significant/marginal split), n_genes,
  n_protein_coding, n_with_description.
- **results_chapter_numbers.txt**: plain-sentence counts for the manuscript.

`dist_to_lead_bp` = signed nearest-edge distance from the lead SNP to the gene: 0 = lead inside the
gene; <0 = gene upstream (lower coordinate); >0 = gene downstream.

## Result summary (this run)

24 input SNPs (15 significant + 9 marginal) -> **20 loci** (11 significant, 9 marginal)
-> **108 candidate genes** (all protein-coding; 12 with a functional description).

| trait      | loci (sig/marg) | genes | with description |
|------------|-----------------|-------|------------------|
| beta-glucan| 11 (8/3)        | 72    | 6                |
| fiber      | 4 (2/2)         | 18    | 2                |
| protein    | 1 (1/0)         | 1     | 0                |
| starch     | 4 (0/4)         | 17    | 4                |

Notable clustered loci: `betaglucan_4H_1` = 4 significant SNPs spanning 146.5 kb;
`betaglucan_3H_1` = 2 SNPs spanning 53.7 kb. The two 1H beta-glucan SNPs (~3.9 Mb apart) form two
separate loci. `betaglucan_5H_1` (marginal) sits in a gene desert (0 genes in window). Two marginal
SNPs added later each form a new singleton locus: `betaglucan_6H_2` (lead 6H:545947347, 12 genes)
and `starch_7H_2` (lead 7H:573606460, 9 genes).

## What is NOT here

- No plots (tables only).
- Functional `description=` exists for only ~6.8% of genes; most candidates carry no annotation
  here (gene_id only). No GO terms in the GFF.
- Crosshap haplotype analysis is the next, separate step.

## Reproducibility / conventions

- All paths absolute. `bedtools` = `/usr/bin/bedtools` (v2.30.0). R = 4.1.2.
- `TMPDIR=/mnt/data/shahar/.tmp` exported in every script; nothing written to `/tmp` or `$HOME`.
  Scripts write only named files under `intermediates/`/`results/` (no `mktemp`), so no temp leak.
