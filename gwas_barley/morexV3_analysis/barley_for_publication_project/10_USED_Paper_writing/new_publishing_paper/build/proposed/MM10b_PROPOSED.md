# Review item MM10b: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM10b** = changes 1 and 4 of MM10, left over from your review.
>
> **Change 1 (Haplotype ¶3):** your wording is restored: "the BLUPs of the trait for which the gene was identified". Everything else in this change is as in MM10: the two effect sizes stated exactly, and the Wilcoxon sentence opening with its aim.
>
> **Change 2 (canonical genes):** unchanged from MM10 change 4. It had no answer in the chat.
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** **Change 1, the trait wording:** your "the trait for which the gene was identified" is in (default). If you want an exact alternative: "the trait whose locus contained the gene". It says how the gene was found (it lies inside that trait's locus), without implying the gene itself was identified by a test. Keep yours, or use this one?
>
> **Q2.** **Change 2:** approve, discard, or change? (It was MM10 change 4, which had no answer.)
>

---

## Change 1 of 2 · M&M · Haplotype analysis ¶3

*Why:* Two effect sizes stated exactly; aim-first Wilcoxon sentence. 'The trait for which the gene was identified' kept (user).

**With the changes marked:**

For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">haplotype </ins>group were excluded<del style="color:#c0392b;background:#fdecea">, and genes</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. Genes</ins> with fewer than two haplotype groups <del style="color:#c0392b;background:#fdecea">were not tested; these genes</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">could not be tested; they</ins> were left out of the multiple-testing correction and are listed with the reason (Online Resource N)<sup>(TODO 224)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. <del style="color:#c0392b;background:#fdecea">Effect size was measured as</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Two effect sizes were calculated:</ins> the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014)<sup>(TODO 225)</sup>, and <del style="color:#c0392b;background:#fdecea">as </del>the difference between the highest and <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">the </ins>lowest group <del style="color:#c0392b;background:#fdecea">means in standard deviations of the tested values</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">mean, divided by the standard deviation of the trait among the tested accessions</ins>. <del style="color:#c0392b;background:#fdecea">In Fig. 4, each</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To show which groups differed, each</ins> haplotype group was compared with the largest group <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">of its gene </ins>by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">; Fig. 4</ins>).<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource N)<sup>(TODO 224)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014)<sup>(TODO 225)</sup>, and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 4).

---

## Change 2 of 2 · M&M · Functional annotation ¶3 (canonical genes) — was MM10 change 4

*Why:* Gene list after a colon; reason in its own sentence; 'we'.

**With the changes marked:**

To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">we located </ins>ten canonical genes<del style="color:#c0392b;background:#fdecea">,</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> in MorexV3:</ins> the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">,</ins> and the endohydrolases *HvGlb1* and *HvGlb2*<del style="color:#c0392b;background:#fdecea">, were located in MorexV3 from their reference protein sequences, because their</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. Because their</ins> published positions refer to genetic maps or earlier genome assemblies<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, they were mapped from their reference protein sequences</ins>. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene (Online Resource N)<sup>(TODO 228)</sup>.

**Reads after approval:**

> To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, we located ten canonical genes in MorexV3: the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*, and the endohydrolases *HvGlb1* and *HvGlb2*. Because their published positions refer to genetic maps or earlier genome assemblies, they were mapped from their reference protein sequences. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene (Online Resource N)<sup>(TODO 228)</sup>.

---

**Your decision:** approve · discard · or tell me what to change.
