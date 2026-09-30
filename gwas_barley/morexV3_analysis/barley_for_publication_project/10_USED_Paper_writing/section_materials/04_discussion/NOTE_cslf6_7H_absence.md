# Why there is no β-glucan signal at *HvCslF6* on 7H — note for the Discussion

Written 2026-09-17. Background for a Discussion paragraph; **not yet written into any manuscript file.**
The Results blueprint reports only the observation and points here.

## The observation

β-glucan maps to 11 loci in this study, on 1H–6H. **There is no β-glucan locus on 7H**, where
*HvCslF6* — the main barley mixed-linkage-glucan synthase — is located.
<!-- src: 01_.../results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv -->

Related, and for the same paragraph: no canonical β-glucan gene sits near any β-glucan lead in this
study, and β-glucan contributes no gene to the final selection.
<!-- src: 06_USED_genes_selected_to_present/README.md § "Beta-glucan contributes no gene here" -->

## The explanation to develop

*HvCslF6* is essential for β-glucan synthesis but carries **very little sequence variation**, so in a
panel where the gene is effectively invariant it cannot be the source of β-glucan variation and will
not appear as a GWAS peak. Detection is therefore **population-dependent** — some panels find it,
most do not. Our absence is the common outcome, not an anomaly.

## The four papers

| paper | what it gives us |
|---|---|
| **Garcia-Gimenez G et al. (2019)** Barley grain (1,3;1,4)-β-glucan content: effects of transcript and sequence variation in genes encoding the corresponding synthase and endohydrolase enzymes. *Sci Rep* 9:17250. https://doi.org/10.1038/s41598-019-53798-8 | **The core citation.** Exome data on 1,336 accessions, 288 of them *spontaneum*. Only **3 SNPs in the coding region** (two synonymous, one non-synonymous A590T) and 12 SNPs across 3 kb of promoter in 35 genotypes. **"None of these SNPs appeared to associate with variation in (1,3;1,4)-β-glucan content."** The authors write of **"atypical conservation of the *HvCslF6* nucleotide sequence, possibly due to its indispensable role in (1,3;1,4)-β-glucan synthesis."** They also note that *cslf6* mutants lacking β-glucan show poor agronomic performance |
| **Geng L et al. (2021)** Identification of genetic loci and candidate genes related to β-glucan content in barley grain by GWAS in the International Barley Core Selected Collection. *Mol Breed* 41:6. https://doi.org/10.1007/s11032-020-01199-5 | Found **no association at *CslF6***, and say so explicitly: *"we were somewhat surprised at the failure of detecting an association of grain β-glucan content with CslF6."* Propose the same two explanations — the gene's essential role and high conservation, or effective alleles too rare for GWAS. Also cites four further studies that missed it (Houston et al. 2014; Oziel et al. 1996; Islamovic et al. 2013; Panozzo et al. 2007) |
| **Houston K et al. (2014)** A genome-wide association scan for (1,3;1,4)-β-glucan content in the grain of contemporary 2-row spring and winter barleys. *BMC Genomics* 15:907. https://doi.org/10.1186/1471-2164-15-907 | 603 cultivated accessions, **no significant association at *CslF6*** on 7H |
| **Shu X, Rasmussen SK (2014)** Quantification of amylose, amylopectin and β-glucan in search for genes controlling the three major quality traits in barley by genome-wide association studies. *Front Plant Sci* 5:197. https://doi.org/10.3389/fpls.2014.00197 | **The counter-example that makes the argument population-dependent rather than universal.** 254 European spring barleys; six significant β-glucan SNPs on 7H around 71 cM, 1–2 cM from *HvCslF6* |

## Still to check before writing

- [ ] **Do we have the variation ourselves?** The strongest version of this argument is not a
      literature claim but a measurement: how much variation segregates at *HvCslF6* in **our 290
      accessions**. The MorexV3 coordinates of *HvCslF6* and the SNP count and allele frequencies in
      that window would settle it directly. **Not yet done, and it would be a new analysis — needs
      approval.**
- [ ] Whether Garcia-Gimenez et al. (2019) report a separate Israeli wild subset. A secondary summary
      claimed ~80 Israeli wild accessions within the 1,336; **the paper text confirms 288
      *spontaneum* but the Israeli figure was not verified.** Do not write it without checking.
- [ ] Verify the Oziel et al. year — Geng et al. cite it as 1996 in one place and 2010 in another.
- [ ] Decide whether the parallel case belongs here too: *HvGlb1* / *HvGlb2*, the β-glucan
      endohydrolases, are also far from every β-glucan lead in this study.
      <!-- src: 05_.../09_USED_canonical_betaglucan_gene_check/ — note that step is stale, built on the retired ±200 kb window -->

## Caution

The absence of a 7H locus at n = 290 under a Bonferroni threshold is **weak evidence on its own** —
failing to detect is not demonstrating absence. Frame the paragraph as *consistent with* the
published pattern, never as proof that *HvCslF6* is invariant in wild barley.
