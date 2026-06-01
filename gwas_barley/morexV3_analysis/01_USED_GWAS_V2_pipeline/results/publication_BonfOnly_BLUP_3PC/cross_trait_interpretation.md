# Cross-trait interpretation — BLUP × 3 PCs

**Status: SUGGESTIVE DRAFT PROSE.** Generated automatically from the analysis outputs; verify each claim and edit before publication.

## Per-trait headline summary

- **β-glucan**: 12 SNP(s) above Bonferroni (−log10p ≥ 6.77). Top SNP: `2H:41996626` on 2H at 41,996,626 bp (−log10p = 8.906). λ_GC = 0.9735.
- **Fiber**: 2 SNP(s) above Bonferroni (−log10p ≥ 6.77). Top SNP: `7H:573606306` on 7H at 573,606,306 bp (−log10p = 7.146). λ_GC = 0.9696.
- **Protein**: 1 SNP(s) above Bonferroni (−log10p ≥ 6.77). Top SNP: `3H:106623911` on 3H at 106,623,911 bp (−log10p = 6.840). λ_GC = 0.9775.
- **Starch**: 0 SNP(s) above Bonferroni (−log10p ≥ 6.77). Top SNP: `6H:525776080` on 6H at 525,776,080 bp (−log10p = 6.444). λ_GC = 0.9846.

## Key observations

1. **Headline β-glucan hit**: `2H:41996626` on chromosome 2H at 41,996,626 bp (−log10p = 8.906, MAF = 0.074, β = -0.3913, SE = 0.06229) — the strongest signal across all four traits. The minor allele (A1 = A) decreases β-glucan content.

2. **chr4H cluster (β-glucan)**: 4 SNPs between 34,825,152 and 34,971,609 bp (~146 kb), all near or above Bonferroni; this fits inside a single ±188 kb LD-defined candidate region and likely tags one underlying locus rather than 4 independent associations.

3. **Cross-trait chromosomal overlap** (Bonferroni-passing hits shared across traits):

   - Chr 3H: traits = betaglucan, protein
   - Chr 4H: traits = betaglucan, fiber

Check whether the proximal lead SNPs are in LD (shared causal locus) or independent.

## Genomic control summary

All four traits show slight deflation (λ_GC 0.9696–0.9846, all above the 0.95 floor) under BLUP × 3 PCs, consistent with adequate population-structure correction on a continuously structured wild barley sample. Deflation is mild enough that the strongest signals remain robust.

## Caveats

- Wild barley structure is continuous (no sharp scree elbow); 3 PCs were chosen for parsimony with adequate genomic control, and the 24-run sensitivity grid (`supp/`) documents the alternatives.
- No FDR-based analysis was performed (Bonferroni-only by design).

## Bridge to candidate-gene analysis

The 15 Bonferroni-passing lead SNPs each define a ±188 kb candidate region (empirical LD-decay window at r² = 0.2). Next step: count annotated Morex V3 genes per window via `bedtools intersect` against the V3 GFF.
