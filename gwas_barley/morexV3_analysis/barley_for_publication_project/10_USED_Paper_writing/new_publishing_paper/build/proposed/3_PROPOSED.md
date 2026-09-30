# Review item 3: Discussion reorganisation, proposed

Proposal only. `Results_Discussion.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where a paragraph moves.

> **❓ Questions for you, please answer together with your decision**
>
> **Q1.** Name of the new section: **"From associated loci to candidate genes"**. Another option: *"From loci to candidate genes: what the analysis can resolve"*. OK, or your own?
>
> **Q2.** Change B cuts *"Moreover, a single grouping setting was applied to all genes, and a setting suited to SNP-dense genes can fail to group the associated SNPs in sparse ones."* It alludes to the 7H parameter change, which is only in M&M, so the reader can't connect it. OK to cut?
>
> **Q3.** Change A drops *"it shared no genomic region with starch or fiber and contributed no candidate gene"* from the β-glucan opening. Both facts are reported in the Results (Ch. 2 ¶4, Ch. 3 ¶2), and the Discussion does not build on them. OK?

---

## 1. The new structure: 6 sections become 5

| now | proposed | text |
|---|---|---|
| §1 Ecological differentiation of grain composition | **§1**, same | unchanged |
| §2 A shared carbon-allocation axis… | **§2**, same | unchanged |
| §3 Heritability and GWAS… ¶1 (heritability) | **§3 ¶1** | unchanged |
| §3 ¶2 (rare alleles, direction, caveat) | **§3 ¶2** | unchanged |
| §5 β-glucan: strong associations… (own section) | **§3 ¶3**, moved | opening shortened (**change A**) |
| §3 ¶3 (statistical rank) **+** §4 What the design can and cannot resolve | **§4 ¶1**, new heading *From associated loci to candidate genes* | merged into one paragraph and shortened (**change B**) |
| §3 ¶4 (two genes + closing sentence) | **§4 ¶2** | one function sentence for each of the four genes, then Sariel's closing sentence (**change C**) |
| §6 Wild alleles for cultivated barley | **§5**, same | unchanged; leads into the Conclusions |

**Why this order:** §3 is now only about genetic architecture: what the GWAS found, including β-glucan's strong but gene-less signal, which is where Sariel discussed β-glucan's mapping. §4 is only about getting from loci to genes: how genes were chosen, what the two steps can miss, and what the four genes suggest. It ends on Sariel's closing sentence, which now closes the gene discussion rather than sitting in the middle of the Discussion. §5 then leads into the Conclusions.

**Length:** these paragraphs go from **929 to 906 visible words**.

---

## Change A · §3, last paragraph: β-glucan, moved here

<span style="color:#1f4e9c;font-weight:bold;font-style:italic">(→ moved from its own section, §5 "β-glucan: strong associations without candidate genes", whose heading is removed)</span>

*Why:* its first sentence restated the Results (Q3). The paragraph now starts from its own finding, the distance from the known pathway, and the 45-word sentence is split in two. Everything from *"The absence of a signal at HvCslF6…"* onward is unchanged, TODOs 408–410 included.

**With the changes marked:**

β-glucan produced some of the strongest associations of the study, yet <del style="color:#c0392b;background:#fdecea">it shared no genomic region with starch or fiber and contributed no candidate gene. Nor did its loci lie</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">none of its loci lay</ins> near the canonical mixed-linkage glucan pathway<del style="color:#c0392b;background:#fdecea">: the</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. The</ins> nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP <sup>(TODO 410)</sup>, and no β-glucan locus was found on chromosome 7H (Fig. 3a), which carries *HvCslF6* and *HvGlb2*. The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Burton et al. 2008) <sup>(TODO 409)</sup>, is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content, which the authors attributed to the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional, as other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021), whereas a panel of European spring barley detected its region (Shu and Rasmussen 2014) <sup>(TODO 408)</sup>, indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

**Reads after approval:**

> β-glucan produced some of the strongest associations of the study, yet none of its loci lay near the canonical mixed-linkage glucan pathway. The nearest of the ten canonical synthases and endohydrolases, *HvGlb1*, lay 83.5 Mb from any β-glucan lead SNP <sup>(TODO 410)</sup>, and no β-glucan locus was found on chromosome 7H (Fig. 3a), which carries *HvCslF6* and *HvGlb2*. The absence of a signal at *HvCslF6*, the major synthase of grain β-glucan (Burton et al. 2008) <sup>(TODO 409)</sup>, is consistent with the unusual conservation of this gene. A previous study of 1,336 barley accessions, 288 of them wild, found only three coding SNPs segregating in *HvCslF6*, none of them associated with grain β-glucan content, which the authors attributed to the indispensable role of the gene in β-glucan synthesis (Garcia-Gimenez et al. 2019). A gene can thus be essential for a trait and still show no association, because association requires segregating functional variation. This outcome is common rather than exceptional, as other association studies also failed to detect *HvCslF6* (Houston et al. 2014; Geng et al. 2021), whereas a panel of European spring barley detected its region (Shu and Rasmussen 2014) <sup>(TODO 408)</sup>, indicating that its detection depends on the population. The β-glucan loci identified here, all located away from these genes, therefore point to sources of variation in wild barley that remain to be characterized.

---

## Change B · new §4 "From associated loci to candidate genes", ¶1: the "rank" and "two filters" paragraphs merged

<span style="color:#1f4e9c;font-weight:bold;font-style:italic">(→ first part: the former §3 ¶3 "Statistical rank…"; second part: the former §4 "What the design can and cannot resolve", whose heading is removed)</span>

*Why:* both paragraphs answer one question: how the genes were chosen from the loci, and what that can miss. So they become one paragraph. The cuts: the clause on how loci are cut at large SNP gaps (method, in M&M); *"Accordingly, the absence of a candidate gene… is weak evidence…"*, which the last sentence already says; and the grouping-setting sentence (Q2). The long GDSL sentence is split in two.

**With the changes marked:**

Statistical rank was <del style="color:#c0392b;background:#fdecea">also </del>a poor guide to biological candidacy. The most significant genes were often the least interpretable, unannotated or derived from mobile elements, whereas the genes with a clear route to the trait ranked lower. Because genes within a locus share one linkage-disequilibrium block, their significance largely reflects linkage with the lead signal rather than their own function, so significance identifies the locus, and biological relevance has to choose the gene. <del style="color:#c0392b;background:#fdecea">The genes carried forward passed two successive filters, and each of them can miss genuine candidates.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Two steps of the analysis can also miss genuine candidates.</ins> <del style="color:#c0392b;background:#fdecea">The first, the</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The</ins> definition of loci<del style="color:#c0392b;background:#fdecea">,</del> delimits each association by linkage disequilibrium with its lead SNP<del style="color:#c0392b;background:#fdecea"> and ends it at large gaps between SNPs, so the size of a locus is influenced by both the underlying genetic signal and the thresholds used to define it. A</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, so a</ins> gene lying just outside a narrow locus is <del style="color:#c0392b;background:#fdecea">therefore </del>never considered<del style="color:#c0392b;background:#fdecea">, as illustrated by the GDSL esterase/lipase at 7H, which</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">. The GDSL esterase/lipase at 7H illustrates this: it</ins> lay beside two of the strongest signals of the study and was found only by examining the region directly. <del style="color:#c0392b;background:#fdecea">Accordingly, the absence of a candidate gene at a locus is weak evidence that none exists. </del><del style="color:#c0392b;background:#fdecea">The second filter, the</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The</ins> haplotype analysis<del style="color:#c0392b;background:#fdecea">,</del> requires sufficient variation within a gene to define haplotype groups, so a gene that could not be tested is not a gene without effect. <del style="color:#c0392b;background:#fdecea">Moreover, a single grouping setting was applied to all genes, and a setting suited to SNP-dense genes can fail to group the associated SNPs in sparse ones. </del>The three genes carried forward are therefore those that passed both filters, not a complete list of the genes underlying these traits.

**Reads after approval:**

> Statistical rank was a poor guide to biological candidacy. The most significant genes were often the least interpretable, unannotated or derived from mobile elements, whereas the genes with a clear route to the trait ranked lower. Because genes within a locus share one linkage-disequilibrium block, their significance largely reflects linkage with the lead signal rather than their own function, so significance identifies the locus, and biological relevance has to choose the gene. Two steps of the analysis can also miss genuine candidates. The definition of loci delimits each association by linkage disequilibrium with its lead SNP, so a gene lying just outside a narrow locus is never considered. The GDSL esterase/lipase at 7H illustrates this: it lay beside two of the strongest signals of the study and was found only by examining the region directly. The haplotype analysis requires sufficient variation within a gene to define haplotype groups, so a gene that could not be tested is not a gene without effect. The three genes carried forward are therefore those that passed both filters, not a complete list of the genes underlying these traits.

---

## Change C · new §4 ¶2: the four genes

<span style="color:#1f4e9c;font-weight:bold;font-style:italic">(→ the former §3 ¶4, now after the merged paragraph)</span>

*Why:* one function sentence for each of the four genes, in the same form ("X suggests that… consistent with… (citation)"), your rule from item 1.11. The closing *"These four genes therefore…"* now follows from what precedes it. The two new sentences cite the references already used for these genes in Results Ch. 3; TODO 304 and 305 cover their verification. The short names match the Results.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Two of the candidate genes connect to further evidence.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Each candidate gene connects its trait to a specific process in the developing grain. The acyltransferase *GPAT6* suggests that cutin deposition in the pericarp and testa contributes to the insoluble fiber fraction, consistent with the central role of GPAT6 in fruit cutin biosynthesis in tomato (Petit et al. 2016). The glucanase *GH17* suggests that the hydrolysis of cell-wall glucans contributes to grain fiber, consistent with the role of its family, which includes barley's own (1,3;1,4)-β-glucan endohydrolases, in degrading mixed-linkage glucan (Hrmova and Fincher 2024).</ins> The phosphate transporter <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">*PHT4;3* </ins>suggests that plastid phosphate homeostasis contributes to the metabolic conditions required for starch accumulation, consistent with the impaired grain filling reported in rice phosphate-transporter mutants (Ma et al. 2021). The GDSL esterase/lipase at 7H suggests that the modification of cell-wall esters contributes to grain fiber, consistent with the roles of this family in arabinoxylan deacetylation and cutin deposition (Zhang et al. 2019; Girard et al. 2012). Its associated SNPs were confined to the putative promoter, while variants within the gene body showed no association, pointing to regulatory rather than coding variation. These four genes therefore identify distinct biological pathways for future investigation: cutin deposition, glucan hydrolysis, plastid phosphate transport and cell-wall ester modification.

**Reads after approval:**

> Each candidate gene connects its trait to a specific process in the developing grain. The acyltransferase *GPAT6* suggests that cutin deposition in the pericarp and testa contributes to the insoluble fiber fraction, consistent with the central role of GPAT6 in fruit cutin biosynthesis in tomato (Petit et al. 2016). The glucanase *GH17* suggests that the hydrolysis of cell-wall glucans contributes to grain fiber, consistent with the role of its family, which includes barley's own (1,3;1,4)-β-glucan endohydrolases, in degrading mixed-linkage glucan (Hrmova and Fincher 2024). The phosphate transporter *PHT4;3* suggests that plastid phosphate homeostasis contributes to the metabolic conditions required for starch accumulation, consistent with the impaired grain filling reported in rice phosphate-transporter mutants (Ma et al. 2021). The GDSL esterase/lipase at 7H suggests that the modification of cell-wall esters contributes to grain fiber, consistent with the roles of this family in arabinoxylan deacetylation and cutin deposition (Zhang et al. 2019; Girard et al. 2012). Its associated SNPs were confined to the putative promoter, while variants within the gene body showed no association, pointing to regulatory rather than coding variation. These four genes therefore identify distinct biological pathways for future investigation: cutin deposition, glucan hydrolysis, plastid phosphate transport and cell-wall ester modification.

---

**Your decision:** approve · discard · or tell me what to change (by change letter A, B or C).
