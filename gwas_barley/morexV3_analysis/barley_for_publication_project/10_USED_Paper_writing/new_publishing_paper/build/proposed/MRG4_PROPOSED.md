# Review item MRG4: proposed changes

Proposal only. `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MRG4 — figure citations outside the merged chapter** (user decision 2026-10-17: the M&M, the Discussion and the open-items list follow the combined Fig. 3). Panel key of the combined figure: violins **a** *GPAT6*, **b** *GH17*, **c** *PHT4;3*, **d** GDSL fiber, **e** GDSL starch; barcodes **f** *GPAT6*, **g** *GH17*, **h** *PHT4;3*, **i** GDSL. Six changes in the manuscript (below); only figure numbers and panel letters change, plus "four figures" → "three" in the open-items list.
>
> **Also in this item: two Online Resource files cite the old Fig. 4** (found by searching the built files and their builders). On your approval I change the two builder scripts and rebuild the supplement (`bash supplementary/scripts/build_all.sh`, which also refreshes `supplementary/for_submission/`):<br>• **Online Resource 9** (PDF legend, `build_ESM_pdf.py`): "The haplotype groups are those of ~~Figs. 3 and 4; their tests are in Online Resource 6 and Fig. 4~~ **Fig. 3; their tests are in Fig. 3 and Online Resource 6**".<br>• **Online Resource 11** (xlsx, sheet "About", `build_ESM_xlsx.py`): "Groups as in Fig. 3 (GPAT6, GH17, PHT4;3~~;~~**:** crosshap ε = 0.6, MGmin = 2~~) and Fig. 4 (~~**; **GDSL gene~~;~~**:** ε = 0.9, MGmin = 3)."<br>No number and no other text in the supplement changes.
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** Change 2 makes the M&M haplotype-test citation more specific, "Fig. 3" → "Fig. 3a–c", because the combined figure now also shows the GDSL gene, which the next M&M paragraph describes separately (with its own citation, change 3). Keep "Fig. 3a–c", or leave it as "Fig. 3"?
>
> **Q2.** The review package sent to S. Hübner on 2026-10-04 (`supplementary/for_review_Hubner/`) also says "Fig. 4" in the same two places. I suggest leaving it as it is, a record of what he received. Agree, or should I rebuild it too?
>

---

## Change 1 of 6 · Open-items list, item 8 (Final Word file)

*Why:* The manuscript now has three figures.

**With the changes marked:**

1. **Title page (not written yet):** the author names and their affiliations (institution, department, city, country); the corresponding author's e-mail and ORCID (ORCID is mandatory); and the Acknowledgments, which TAG places on the title page, with the funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, affiliation and e-mail also go into the title block of every Online Resource file, where they are placeholders now; the files are then rebuilt.
2. **Author contribution statement (not written yet):** one short line per author, using initials. TAG prints it before the Acknowledgments, and it is also entered in the submission system.
3. **Statements and Declarations (not written yet; they go after the References):** TAG returns a manuscript that lacks any of them. Needed: Funding (agencies and grant numbers); Competing interests; Author contributions; Data availability; Ethics approval (not applicable to this study). The Data availability text is already drafted and names Online Resource 1 for the genotype BLUPs; it still needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub link (with its Zenodo archive) of the analysis code (comment 234).
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary): done.** The 13 files are built, numbered in the order they are first cited, cited in the text, and captioned in the last section ("Supplementary Information"). Left to do: put the author names, affiliation and e-mail into the files (item 1); at submission, upload the 13 files and enter each one's caption in the submission system; and decide whether the extra, uncited supplement of raw-phenotype checks should join the supplement (it would then need a number and a citation in the text).
6. **Other open comments** (each marked by a Word comment).
7. **Cover letter (not written yet):** give the accession number of the raw reads (ENA PRJEB79623); state any restriction on sharing the seeds of the accessions (the same statement belongs at the end of the Materials and methods); and suggest reviewers who are independent of the authors, from different countries, with institutional e-mail addresses.
8. **Final Word file:** 10-pt Times Roman; automatic page numbering; no field functions. The <del style="color:#c0392b;background:#fdecea">four</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">three</ins> figures are embedded in the text and are also uploaded as separate TIFF files, Fig1 to Fig<del style="color:#c0392b;background:#fdecea">4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3</ins> (600 dpi, 174 mm wide). Before submitting, delete this list and every remaining Word comment.
9. *(Optional)* **Figure accessibility:** TAG recommends patterns or symbols in addition to colour, but Figs. 2–4 tell loci and haplotype groups apart by colour alone.

**Reads after approval:**

> 1. **Title page (not written yet):** the author names and their affiliations (institution, department, city, country); the corresponding author's e-mail and ORCID (ORCID is mandatory); and the Acknowledgments, which TAG places on the title page, with the funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, affiliation and e-mail also go into the title block of every Online Resource file, where they are placeholders now; the files are then rebuilt.
2. **Author contribution statement (not written yet):** one short line per author, using initials. TAG prints it before the Acknowledgments, and it is also entered in the submission system.
3. **Statements and Declarations (not written yet; they go after the References):** TAG returns a manuscript that lacks any of them. Needed: Funding (agencies and grant numbers); Competing interests; Author contributions; Data availability; Ethics approval (not applicable to this study). The Data availability text is already drafted and names Online Resource 1 for the genotype BLUPs; it still needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub link (with its Zenodo archive) of the analysis code (comment 234).
4. **References that still need attention** (each marked by a Word comment): Friedman and Atsmon 1988 and Holland et al. 2003 could not be checked, because no PDF was available (comments 102, 401, 257, 501); the page numbers of Keeling et al. 1993 need confirming (502); two citations only partly support their sentences, which were kept as written (104, 428); and the source of the soil variables is not described in the cited paper (216).
5. **Online Resources (supplementary): done.** The 13 files are built, numbered in the order they are first cited, cited in the text, and captioned in the last section ("Supplementary Information"). Left to do: put the author names, affiliation and e-mail into the files (item 1); at submission, upload the 13 files and enter each one's caption in the submission system; and decide whether the extra, uncited supplement of raw-phenotype checks should join the supplement (it would then need a number and a citation in the text).
6. **Other open comments** (each marked by a Word comment).
7. **Cover letter (not written yet):** give the accession number of the raw reads (ENA PRJEB79623); state any restriction on sharing the seeds of the accessions (the same statement belongs at the end of the Materials and methods); and suggest reviewers who are independent of the authors, from different countries, with institutional e-mail addresses.
8. **Final Word file:** 10-pt Times Roman; automatic page numbering; no field functions. The three figures are embedded in the text and are also uploaded as separate TIFF files, Fig1 to Fig3 (600 dpi, 174 mm wide). Before submitting, delete this list and every remaining Word comment.
9. *(Optional)* **Figure accessibility:** TAG recommends patterns or symbols in addition to colour, but Figs. 2–4 tell loci and haplotype groups apart by colour alone.

---

## Change 2 of 6 · M&M, Haplotype analysis ¶3

*Why:* Panels a–c are the genes tested in this pipeline; the GDSL gene (d, e) is cited in the 7H paragraph (change 3).

**With the changes marked:**

For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource 6). *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014), and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 3<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">a–c</ins>).

**Reads after approval:**

> For each gene, the BLUPs of the trait for which the gene was identified were compared among haplotype groups with a Kruskal–Wallis test (Kruskal and Wallis 1952). Accessions not assigned to a haplotype group were excluded. Genes with fewer than two haplotype groups could not be tested; they were left out of the multiple-testing correction and are listed with the reason (Online Resource 6). *P* values were corrected by the Benjamini–Hochberg procedure within each trait (Benjamini and Hochberg 1995), and genes with *q* ≤ 0.05 were considered significant. Two effect sizes were calculated: the Kruskal–Wallis η², the proportion of the variance in trait ranks explained by haplotype group (Tomczak and Tomczak 2014), and the difference between the highest and the lowest group mean, divided by the standard deviation of the trait among the tested accessions. To show which groups differed, each haplotype group was compared with the largest group of its gene by a two-sided Wilcoxon rank-sum test (Wilcoxon 1945), with Holm correction within each gene (Holm 1979; Fig. 3a–c).

---

## Change 3 of 6 · M&M, Functional annotation and choice of candidate genes ¶4 (7H gene)

*Why:* The GDSL gene's violins are now Fig. 3d (fiber) and 3e (starch).

**With the changes marked:**

To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes (Online Resource 9). As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. <del style="color:#c0392b;background:#fdecea">4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3d, e</ins>).

**Reads after approval:**

> To search for a candidate gene at the shared fiber and starch loci on chromosome 7H, we examined the nearest annotated gene, *HORVU.MOREX.r3.7HG0729030*. It was annotated as described above, including the SignalP and Phobius predictions of InterProScan. Each SNP from 1 kb downstream to 3 kb upstream of the gene was then positioned relative to its annotated transcription start site. LD among the significant SNPs was taken from the r² matrix of imputed genotypes used for the haplotype analysis. At the crosshap settings used for all candidate genes, none of the significant SNPs of this gene was included in a marker group, and the gene was therefore analyzed with MGmin = 3 and ε = 0.9, at which the three significant SNPs were included in the marker group that defined the haplotypes (Online Resource 9). As a single pre-specified gene, it was tested without multiple-testing correction, and its two haplotype groups were also compared by a two-sided Wilcoxon rank-sum test (Fig. 3d, e).

---

## Change 4 of 6 · M&M, Comparison with elite cultivars ¶2

*Why:* Barcodes of the three genes are now Fig. 3f–h; the 7H gene's barcode is Fig. 3i.

**With the changes marked:**

For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. 3<del style="color:#c0392b;background:#fdecea">d–f</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">f–h</ins>). The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups. For the 7H gene (Fig. <del style="color:#c0392b;background:#fdecea">4c</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">3i</ins>), the same five cultivars and procedure were used without a new call-rate screen.

**Reads after approval:**

> For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. 3f–h). The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups. For the 7H gene (Fig. 3i), the same five cultivars and procedure were used without a new call-rate screen.

---

## Change 5 of 6 · M&M, Geographic origin of allele and haplotype carriers

*Why:* The four presented genes are all in Fig. 3.

**With the changes marked:**

To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (<del style="color:#c0392b;background:#fdecea">Figs. 3, 4</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Fig. 3</ins>) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure (Online Resources 11, 12). To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

**Reads after approval:**

> To test whether the carriers of associated alleles and haplotypes came from particular regions, we assigned each accession to its sampling site and to one of the six regions of the sampling design. The minor- and major-allele carriers of each of the 36 lead SNPs were compared, as was each haplotype group of the four presented genes (Fig. 3) with the other groups of its gene; accessions with a missing call or without a haplotype group were excluded. Enrichment was measured as the difference between the two sets in the proportion of desert accessions and tested by permuting the region labels among the 29 sites 10,000 times, keeping the accessions of each site together (two-sided); *P* values were corrected by the Benjamini–Hochberg procedure (Online Resources 11, 12). To separate the allelic effect from geographic origin, each minor-allele carrier was also compared with the major-allele carriers of its own site: we calculated the proportion of carriers whose BLUP differed from the mean of these site-mates in the direction of the allelic effect.

---

## Change 6 of 6 · Discussion, Wild alleles for cultivated barley

*Why:* Barcodes: PHT4;3 is now Fig. 3h (was 3f), GH17 is Fig. 3g (was 3e).

**With the changes marked:**

The five elite cultivars carried nearly identical genotypes at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every called site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 3<del style="color:#c0392b;background:#fdecea">f</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">h</ins>), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 3<del style="color:#c0392b;background:#fdecea">e</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">g</ins>). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resources 11, 12). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

**Reads after approval:**

> The five elite cultivars carried nearly identical genotypes at *GPAT6*, *GH17* and *PHT4;3*, as expected in modern germplasm. At *GH17* and *PHT4;3*, the cultivars carried the reference allele at every called site. Since the reference genome, Morex, is itself an elite cultivar, their match with the wild groups at these genes should be interpreted with caution: at *PHT4;3* the cultivars appeared identical to the low-starch group C at the shared sites (Fig. 3h), most likely because both carry the reference allele, whereas at *GH17* they matched none of the wild groups (Fig. 3g). In contrast, *GPAT6* and the GDSL esterase/lipase, two independent loci, provided a clear result. At both genes the cultivars carried non-reference alleles at many sites, so their genotypes do not simply reflect the reference genome, and the high-fiber wild haplotype was absent from all five cultivars. At *GPAT6*, this haplotype came from three desert and two coast–desert populations rather than from a single site (Online Resources 11, 12). These haplotypes may therefore represent wild variation for increasing fiber in cultivated barley.

---

**Your decision:** approve · discard · or tell me what to change.
