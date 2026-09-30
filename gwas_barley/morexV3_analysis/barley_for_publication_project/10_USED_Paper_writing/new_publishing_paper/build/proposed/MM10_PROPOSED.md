# Review item MM10: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM10 of the M&M voice review** (item 6, voice and clarity, part 3 of 3): Haplotype analysis ¶3 and the three paragraphs of Functional annotation and choice of candidate genes. No facts change.
>
> **Haplotype ¶3:** "the trait for which the gene was identified" becomes "the trait of its locus". Two effect sizes are now introduced as two, and the SD one is stated exactly: the difference between the highest and lowest group mean, divided by the SD of the trait among the tested accessions (this also covers the SD-wording item from the minor list). The Wilcoxon sentence now opens with its aim.
>
> **Annotation ¶1 (A9):** the "Because 3 of 55 … the 23 genes…" jump is split into two sentences, with "we". "It named a characterized protein" now spells out the filter the code applies. The nr sentence is split. "The call of the first source in this order that produced one" now names the order (Swiss-Prot, InterPro, nr). The Table 2 clause is removed, because the Table 2 footnote already says "Best UniProt Swiss-Prot match".
>
> **Annotation ¶2:** "chosen on protein function and the literature, as genes with…" now reads directly, with "we".
>
> **Annotation ¶3 (A7):** the gene list set off by commas mid-sentence moves after a colon, and the reason ("because their published positions…") gets its own sentence.
>
> Length: 441 → 454 words (a clarity pass).
>

## Change 1 of 4 · M&M · Haplotype analysis ¶3

*Why:* Trait of the locus; two effect sizes stated exactly; aim-first Wilcoxon sentence.

**With the changes marked:**

For each gene, the BLUPs of the trait <del style="color:#c0392b;background:#fdecea">for which the gene was identified</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">of its locus</ins> were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">haplotype </ins>group were excluded<del style="color:#c0392b;background:#fdecea">, and genes</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. Genes</ins> with fewer than two haplotype groups <del style="color:#c0392b;background:#fdecea">were not tested; these genes</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">could not be tested; they</ins> were left out of the multiple-testing correction and are listed with the reason (Online Resource N)<sup>(TODO 224)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. <del style="color:#c0392b;background:#fdecea">Effect size was measured as</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Two effect sizes were calculated:</ins> the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014)<sup>(TODO 225)</sup>, and <del style="color:#c0392b;background:#fdecea">as </del>the difference between the highest and <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">the </ins>lowest group <del style="color:#c0392b;background:#fdecea">means in standard deviations of the tested values</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">mean, divided by the standard deviation of the trait among the tested accessions</ins>. <del style="color:#c0392b;background:#fdecea">In Fig. 4, each</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To show which groups differed, each</ins> haplotype group was compared with the largest group <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">of its gene </ins>by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">; Fig. 4</ins>).<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> For each gene, the BLUPs of the trait of its locus were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource N)<sup>(TODO 224)</sup>. *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014)<sup>(TODO 225)</sup>, and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 4).

---

## Change 2 of 4 · M&M · Functional annotation ¶1

*Why:* Reason and action split; filter made explicit; nr sentence split; source order named; Table 2 clause removed (in its footnote).

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Because only</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Only</ins> 3 of the 55 candidate genes had a functional description in the MorexV3 annotation<del style="color:#c0392b;background:#fdecea">, the</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. We therefore annotated the</ins> 23 genes with significant haplotype associations <del style="color:#c0392b;background:#fdecea">were annotated </del>from <del style="color:#c0392b;background:#fdecea">their protein sequences, using the longest isoform of each gene</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">the protein sequence of their longest isoform</ins> in the MorexV3 high-confidence proteome (Mascher et al. 2021). Each protein was first searched against UniProtKB/Swiss-Prot release 2026_03 (UniProt Consortium 2025)<sup>(TODO 226)</sup> by BLASTP in DIAMOND v2.0.14 (Buchfink et al. 2021), and its best hit was accepted when the e-value was ≤ 10⁻⁵, the hit covered ≥ 50% of the query at ≥ 30% identity, and <del style="color:#c0392b;background:#fdecea">it named a characterized protein</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">the hit was not described as uncharacterized, hypothetical, predicted or of unknown function</ins>. Protein domains and Gene Ontology terms were then assigned with InterProScan v5.78-109.0 (Jones et al. 2014)<del style="color:#c0392b;background:#fdecea">, and proteins resolved by neither source were searched by BLASTP</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. Proteins without a result from either source were searched with BLASTP (BLAST+ v2.12.0; Camacho et al. 2009)</ins> against the NCBI non-redundant protein database, restricted to green plants (Viridiplantae; <del style="color:#c0392b;background:#fdecea">BLAST+ v2.12.0; Camacho et al. 2009; </del>accessed 9 September 2026), with the same criteria. Each gene <del style="color:#c0392b;background:#fdecea">received the call of the first source in this order that produced one</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">took its annotation from the first of these sources that gave one, in the order Swiss-Prot, InterPro and nr</ins> (Online Resource N)<sup>(TODO 227)</sup><del style="color:#c0392b;background:#fdecea">; the identity and query coverage in Table 2 are those of the best Swiss-Prot hit</del>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> Only 3 of the 55 candidate genes had a functional description in the MorexV3 annotation. We therefore annotated the 23 genes with significant haplotype associations from the protein sequence of their longest isoform in the MorexV3 high-confidence proteome (Mascher et al. 2021). Each protein was first searched against UniProtKB/Swiss-Prot release 2026_03 (UniProt Consortium 2025)<sup>(TODO 226)</sup> by BLASTP in DIAMOND v2.0.14 (Buchfink et al. 2021), and its best hit was accepted when the e-value was ≤ 10⁻⁵, the hit covered ≥ 50% of the query at ≥ 30% identity, and the hit was not described as uncharacterized, hypothetical, predicted or of unknown function. Protein domains and Gene Ontology terms were then assigned with InterProScan v5.78-109.0 (Jones et al. 2014). Proteins without a result from either source were searched with BLASTP (BLAST+ v2.12.0; Camacho et al. 2009) against the NCBI non-redundant protein database, restricted to green plants (Viridiplantae; accessed 9 September 2026), with the same criteria. Each gene took its annotation from the first of these sources that gave one, in the order Swiss-Prot, InterPro and nr (Online Resource N)<sup>(TODO 227)</sup>.

---

## Change 3 of 4 · M&M · Functional annotation ¶2 (choice of candidates)

*Why:* Direct wording, 'we'.

**With the changes marked:**

Among the 23 significant genes, <del style="color:#c0392b;background:#fdecea">candidates were chosen on protein function and the literature, as genes with a documented or plausible role</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">we selected as candidates the genes with a documented or plausible role, based on their protein function and the literature,</ins> in the biosynthesis, transport or regulation of the corresponding grain component<del style="color:#c0392b;background:#fdecea">; *q* values</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. The *q* values</ins> and effect sizes were not used <del style="color:#c0392b;background:#fdecea">to rank them</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">in this choice</ins>.

**Reads after approval:**

> Among the 23 significant genes, we selected as candidates the genes with a documented or plausible role, based on their protein function and the literature, in the biosynthesis, transport or regulation of the corresponding grain component. The *q* values and effect sizes were not used in this choice.

---

## Change 4 of 4 · M&M · Functional annotation ¶3 (canonical genes)

*Why:* Gene list after a colon; reason in its own sentence; 'we'.

**With the changes marked:**

To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">we located </ins>ten canonical genes<del style="color:#c0392b;background:#fdecea">,</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> in MorexV3:</ins> the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">,</ins> and the endohydrolases *HvGlb1* and *HvGlb2*<del style="color:#c0392b;background:#fdecea">, were located in MorexV3 from their reference protein sequences, because their</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. Because their</ins> published positions refer to genetic maps or earlier genome assemblies<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, they were mapped from their reference protein sequences</ins>. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene (Online Resource N)<sup>(TODO 228)</sup>.

**Reads after approval:**

> To relate the β-glucan loci to the known genes of (1,3;1,4)-β-glucan synthesis and hydrolysis, we located ten canonical genes in MorexV3: the synthases *HvCslF3*, *HvCslF4*, *HvCslF6*, *HvCslF7*, *HvCslF8*, *HvCslF9*, *HvCslF10* and *HvCslH1*, and the endohydrolases *HvGlb1* and *HvGlb2*. Because their published positions refer to genetic maps or earlier genome assemblies, they were mapped from their reference protein sequences. Each protein was searched against the MorexV3 genome with tBLASTn (e-value ≤ 10⁻¹⁰) and against the high-confidence proteome with BLASTP (BLAST+ v2.12.0), and its MorexV3 ortholog was taken as the protein hit with the highest identity. All ten orthologs (≥ 97% identity) were reciprocal best hits, and their genome and proteome positions agreed. The distance of each gene to the nearest β-glucan lead SNP on the same chromosome was measured from the nearest edge of the gene (Online Resource N)<sup>(TODO 228)</sup>.

---

**Your decision:** approve · discard · or tell me what to change.
