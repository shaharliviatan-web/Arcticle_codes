# 09_USED_canonical_betaglucan_gene_check

Goal: test whether the canonical (1,3;1,4)-beta-glucan genes were missed by our
+/-200 kb candidate-gene window (the window is tight because LD decays fast in this
wild panel). We locate each canonical gene on Morex V3 by SEQUENCE (assembly-
independent) and measure its distance to every beta-glucan / fiber / starch GWAS
lead SNP.

Critical method constraint honoured: literature gene positions (genetic-map cM, or
Morex V1/V2 coordinates) are NOT used - they are not comparable to our V3 SNP
positions. Each gene is mapped from its reference protein sequence onto our local
Morex V3 files.

## Canonical genes + reference accessions (provenance)

10 genes: the (1,3;1,4)-beta-glucan synthases HvCslF6, HvCslH1, HvCslF9, HvCslF7,
HvCslF3, HvCslF4, HvCslF8, HvCslF10 and the endohydrolases HvGlb1 (EI) /
HvGlb2 (EII). Reference accessions (identified by gene-name + Hordeum vulgare
searches at UniProt REST / NCBI, then fetched by accession for reproducibility;
see `inputs/canonical_refs_accessions.tsv`, `inputs/canonical_refs_provenance.tsv`,
`logs/versions.txt`):

| gene | accession | source |
|---|---|---|
| HvCslF3/F4/F6/F7/F8/F9/F10 | B1P2T2 / B1P2T3 / B1P2T4 / B1P2T5 / B1P2T6 / B1P2T7 / B1P2T8 | UniProt (Burton et al.) |
| HvCslH1 | ACN67534.1 | NCBI (Doblin et al. 2009) |
| HvGlb1 (EI) | Q02345 | UniProt |
| HvGlb2 (EII) | P12257 | UniProt (reviewed; a 312-aa fragment - still locates the locus cleanly) |

Optional HvCslF11/F12/F13 have no clean H. vulgare reference protein and were not
included.

## Method

Local Morex V3 files: genome `Barley_MorexV3_pseudomolecules.fasta` (headers chr1H..),
HC proteome `Hv_Morex.pgsb.Jul2020.HC.aa.fa`, HC GFF3 `Hv_Morex.pgsb.Jul2020.HC.gff3`.
Chromosome names are normalised (chr1H -> 1H) to match the SNP tables.

1. `01_fetch_refs.sh` - fetch the 10 reference proteins -> `inputs/canonical_refs.faa`.
2. `02_makedb_and_blast.sh` - build BLAST dbs; locate by sequence two ways:
   - tblastn refs vs the genome -> genomic coordinates, annotation-independent;
   - blastp refs vs the HC proteome -> HORVU.MOREX.r3 id by best hit.
3. `03_reciprocal_blast.sh` - blastp the candidate V3 proteins back vs the refs
   (reciprocal best hit).
4. `04_tree.sh` - mafft MSA + ape neighbour-joining tree of the CslF/CslH refs +
   their V3 hits, to confirm each V3 gene clusters with the correct family member.
5. `05_gff3_crosscheck.sh` - grep the GFF3 description= fields for cellulose-
   synthase / glucan-synthase / glucanase terms (CAVEAT below).
6. `06_build_table.R` - assign orthologs, pull coordinates, cross-check tblastn vs
   blastp, and compute distances to the lead SNPs.

Reproduce: `bash scripts/run_all.sh`. TMPDIR pinned to `/mnt/data/shahar/.tmp`;
absolute tool paths; all writes under `/mnt/data/shahar`; ASCII only.

### Paralog resolution (essential - CslF is multi-copy)

References and Morex V3 are both Hordeum vulgare, so a true ortholog aligns at
~99-100% identity while paralogs are <85%. The ortholog is therefore taken as the
forward hit with the HIGHEST % identity, NOT the highest bitscore. This matters:
for **HvCslF4**, bitscore-best is the full-length 7H paralog `7HG0750850` (82.98%),
but the true ortholog is the partial-model `2HG0136610` (99.83%, in the 2H CslF
cluster). Reciprocal best hit (all 10 = Y) and the NJ tree both confirm the
identity-based assignment; the 2H cluster (CslF8/F4/F3/F10 at ~180.2-180.9 Mb)
matches known barley biology.

### GFF3 cross-check caveat

Only ~6-7% of HC genes carry a functional `description=`, so this grep can
CORROBORATE but its silence means nothing. It corroborated CslF6
(`7HG0698110` = "Cellulose synthase-like protein") and Glb2
(`7HG0750120` = "putative beta-1,3-glucanase") but was silent on the 2H CslF
cluster (those models have no description). BLAST is the primary, complete method.
(The companion-paper PDF cross-check was dropped - file not available.)

## Result (headline)

All 10 canonical genes locate cleanly and unambiguously in Morex V3 (reciprocal
best hit Y for all; >=97% identity for every ortholog). **Every one is FAR (>1 Mb)
from every beta-glucan / fiber / starch lead SNP.** The nearest any canonical gene
comes to a lead is ~51 Mb (HvGlb2/EII); the CslF synthases are 126-340 Mb away.

So the canonical (1,3;1,4)-beta-glucan biosynthesis/turnover genes were not merely
missed by the tight +/-200 kb window - they sit tens-to-hundreds of Mb away from
our association signals. No canonical gene falls within a window or just outside
it, so none is a candidate for direct haplotype testing on distance grounds.

## Deliverable

`results/tables/canonical_betaglucan_gene_distances.tsv` - one row per
(canonical gene x trait), 30 rows. Columns: canonical_gene, source_accession,
horvu_id, chr, gene_start, gene_end, strand, pct_identity, pct_coverage, evalue,
reciprocal_best_hit, tblastn_horvu, tblastn_vs_blastp_agree, trait,
nearest_lead_SNP, distance_bp, window_flag (IN_200kb_WINDOW / NEAR_200kb_to_1Mb /
FAR), tested_in_candidate_list, haplotype_significant_fdr, haplotype_trait_fdr_p.

Supporting: `results/tables/csl_family_tree.nwk`,
`results/tables/csl_tree_nearest_ref.tsv`,
`intermediates/gff3_description_hits.tsv`.

## Scope

Locate-and-measure only. No canonical gene is re-run through the crosshap
haplotype test here; since all are FAR, no "just outside the window" case needs the
follow-up decision.
