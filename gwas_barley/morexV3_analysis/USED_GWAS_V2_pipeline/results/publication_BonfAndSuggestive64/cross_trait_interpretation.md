# Cross-trait interpretation — BLUE × 10 PCs

**Status: SUGGESTIVE DRAFT PROSE.** Generated automatically from the analysis outputs; user must verify each claim and edit before publication. Use as a starting scaffold, not as final manuscript text.

## Per-trait headline summary

- **β-glucan**: 15 SNPs above Bonferroni (−log10p ≥ 6.77) (+12 above suggestive 6.4). Top SNP: `2H:41996626` on 2H at 41,996,626 bp (−log10p = 8.226). λ_GC = 0.9763.
- **Fiber**: 3 SNPs above Bonferroni (−log10p ≥ 6.77) (+1 above suggestive 6.4). Top SNP: `4H:24707376` on 4H at 24,707,376 bp (−log10p = 7.896). λ_GC = 0.9745.
- **Protein**: 1 SNPs above Bonferroni (−log10p ≥ 6.77) (+0 above suggestive 6.4). Top SNP: `3H:106623911` on 3H at 106,623,911 bp (−log10p = 6.899). λ_GC = 0.9820.
- **Starch**: 2 SNPs above Bonferroni (−log10p ≥ 6.77) (+1 above suggestive 6.4). Top SNP: `7H:573606460` on 7H at 573,606,460 bp (−log10p = 6.999). λ_GC = 0.9811.

## Key observations

1. **Headline β-glucan hit**: `2H:41996626` on chromosome 2H at 41,996,626 bp (−log10p = 8.226, MAF = 0.074, β = -0.4851, SE = 0.08078). This SNP exceeds the Bonferroni threshold (−log10p = 6.77) and is the strongest signal across all four traits. The minor allele (A1 = A) is associated with decreases β-glucan content.

2. **chr4H cluster (β-glucan)**: 4 SNPs co-located between positions 34,825,152 and 34,971,609 bp (~146 kb window), all near or above Bonferroni. The cluster fits inside a single ±188 kb LD-defined candidate region and likely tags one underlying causal locus rather than 4 independent associations.

3. **Cross-trait chromosomal overlap** (Bonferroni-passing hits on the same chromosome across multiple traits):

   - Chr 3H: traits = betaglucan, protein
   - Chr 4H: traits = betaglucan, fiber
   - Chr 7H: traits = betaglucan, fiber, starch

Worth checking whether the proximal lead SNPs in these regions are in LD (suggesting a shared causal locus) or independent (separate causal mechanisms).

## Genomic control summary

All four traits show slight deflation (λ_GC range 0.9745–0.9820, all above the 0.95 floor) under BLUE × 10 PCs. This is consistent with strong population-structure correction by the 10 PCs (capturing ~33% of total genetic variance) on a continuously structured wild barley sample. The deflation is mild enough that the strongest signals remain robust.

## Caveats for the paper

- Wild barley population structure is continuous (no sharp scree elbow); 10 PCs were chosen as the most defensible point on the λ-vs-N_PCs curve, supported by mean |λ − 1| minimisation across the four traits.
- No FDR-based analysis was performed (Bonferroni-only by design).
- The suggestive threshold (−log10p = 6.4) is provided for completeness; SNPs in the suggestive band are flagged via `above_suggestive_only = TRUE` in `lead_snps.tsv` and should be discussed only as supplementary candidates.

## Bridge to candidate-gene analysis

All 21 Bonferroni-passing lead SNPs define candidate regions of ±188 kb each (the empirical LD-decay window from the prior LD-decay analysis at r² = 0.2). The next analytical step is to count the number of annotated Morex V3 genes within each window using `bedtools intersect` against the V3 GFF, producing per-region gene counts and per-SNP candidate-gene lists.
