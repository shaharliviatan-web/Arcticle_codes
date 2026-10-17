# Review item REF3: proposed changes

Proposal only. `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item REF3: reference list** (References session, 2026-10-01; user request: "attach it also to the bottom of the main manuscript"). The 64 references of `05_references.md`, inserted after the Conclusions and before the Supplementary Information (TAG order: main text, References, Statements and Declarations; the Online Resource captions stay last). Word comments 501 (Holland publisher place) and 502 (Keeling pages) mark the two details that could not be checked. Submission checklist item 4 updated.
>

## Change 1 of 2 · Submission checklist (top of file), item 4

*Why:* Status after the citation audit and REF2/REF3.

**With the changes marked:**

1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References:** <del style="color:#c0392b;background:#fdecea">build the list in TAG style (author–year in the text with no comma; list alphabetical; journal names abbreviated (ISSN LTWA); full DOI links); verify every citation marked "verify" in the TODO comments, and add the β-amylase citation (TODO 411).</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">list built in TAG style and every citation checked (references/citation_audit.md; section "References" below, source 05_references.md). Still open: the Word comments 102 and 401 (Friedman and Atsmon 1988, no PDF), 257 and 501 (Holland et al. 2003, no PDF), 502 (Keeling et al. 1993 pages), 104 and 428 (partly supported claims kept as written), and 216 (source of the soil variables).</ins>
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

**Reads after approval:**

> 1. **Title page:** authors; affiliations (institution, department, city, country); corresponding author's e-mail and **ORCID (mandatory)**; **Acknowledgments** on the title page, with funding organizations written in full (including Equinom Ltd., which lent the NIR instrument). The same author names, and the corresponding author's affiliation and e-mail, go into the **title block of every Online Resource file** (placeholders `[AUTHOR NAMES]`, `[AFFILIATION]`, `[E-MAIL]`; set them in `supplementary/scripts/esm_common.py` and rebuild the files).
2. **Author contribution statement:** short, by initials, for all authors (placed before the Acknowledgments; also entered in the submission system).
3. **Statements and Declarations** (after the References; missing ones → the manuscript is returned as incomplete): Funding (agency + grant numbers); Competing interests; Author contributions; **Data availability** (draft in `build/DECLARATIONS_HANDOVER.md`: needs the Zenodo DOI of the MorexV3 SNP calls and the GitHub/Zenodo link of the code; BLUPs as an Online Resource); Ethics approval: not applicable.
4. **References:** list built in TAG style and every citation checked (references/citation_audit.md; section "References" below, source 05_references.md). Still open: the Word comments 102 and 401 (Friedman and Atsmon 1988, no PDF), 257 and 501 (Holland et al. 2003, no PDF), 502 (Keeling et al. 1993 pages), 104 and 428 (partly supported claims kept as written), and 216 (source of the soil variables).
5. **Online Resources (supplementary):** ✅ numbered 1–13 in citation order; files `supplementary/ESM_1.xlsx` … `ESM_13.pdf` built (`supplementary/scripts/build_all.sh`); ✅ every placeholder replaced; ✅ title block in each file (author names, affiliation and e-mail still placeholders: item 1); ✅ concise captions (last section, "Supplementary Information"; enter the same captions in the submission system); ✅ NIR-calibration sheets of Online Resource 2 (from N. Pintel's files); ✅ citations added in Results Ch. 2 ¶3, Ch. 3 ¶1 (Q13) and of Online Resource 9 in Ch. 3–4 (Q10). Map of all files: `supplementary/README.md`; **files to upload, with their captions: `supplementary/for_submission/`** (zip + `CAPTIONS.md`); list: `supplementary/supplementary.md`. Status of each citation: the "supmat" Word comments (❌ not resolved / ✅ resolved); delete them all before submitting. Not cited, for a later decision: `supplementary/unreferenced/`.
6. **Open TODO comments:** resolve and delete them (`grep -n 'TODO:'` on this file); the wet-lab facts (TODOs 201–213) need Noam Pintel / Evgenii Potapenko.
7. **Cover letter:** refer to the accession number (ENA PRJEB79623) and state any restriction on sharing the seeds; suggested reviewers (independent, mixed countries, institutional e-mails).
8. **Final Word formatting:** 10-pt Times Roman; automatic page numbering; no field functions; submit the figures as `Fig1.tif`–`Fig4.tif` (600 dpi, 174 mm), also embedded in the text.
9. *(Optional)* **Figure accessibility:** TAG asks for patterns in addition to colour; Figs. 2–4 rely on colour (loci, haplotype groups).

---

## Change 2 of 2 · End of the main text, before Supplementary Information

*Why:* New section "References" (user request).

**New text (all of it is added; shown once):**

> # References

AACC (2001) The definition of dietary fiber. Cereal Foods World 46:112–126. https://www.cerealsgrains.org/resources/definitions/Documents/DietaryFiber/DFDef.pdf

Alqudah AM, Sallam A, Baenziger PS, Börner A (2020) GWAS: fast-forwarding gene identification and characterization in temperate cereals: lessons from barley – a review. J Adv Res 22:119–135. https://doi.org/10.1016/j.jare.2019.10.013

Anker-Nilssen K, Sahlstrøm S, Knutsen SH, Holtekjølen AK, Uhlen AK (2008) Influence of growth temperature on content, viscosity and relative molecular weight of water-soluble β-glucans in barley (*Hordeum vulgare* L.). J Cereal Sci 48:670–677. https://doi.org/10.1016/j.jcs.2008.02.009

Badr A, Müller K, Schäfer-Pregl R, El Rabey H, Effgen S, Ibrahim HH, Pozzi C, Rohde W, Salamini F (2000) On the origin and domestication history of barley (*Hordeum vulgare*). Mol Biol Evol 17:499–510. https://doi.org/10.1093/oxfordjournals.molbev.a026330

Bates D, Mächler M, Bolker B, Walker S (2015) Fitting linear mixed-effects models using lme4. J Stat Softw 67:1–48. https://doi.org/10.18637/jss.v067.i01

Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Ser B Methodol 57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

Bogard M, Allard V, Brancourt-Hulmel M, Heumez E, Machet J-M, Jeuffroy M-H, Gate P, Martre P, Le Gouis J (2010) Deviation from the grain protein concentration–grain yield negative relationship is highly correlated to post-anthesis N uptake in winter wheat. J Exp Bot 61:4303–4312. https://doi.org/10.1093/jxb/erq238

Broad Institute (2019) Picard toolkit. Broad Institute, GitHub repository. https://broadinstitute.github.io/picard/. Accessed 1 October 2026

Buchfink B, Reuter K, Drost H-G (2021) Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat Methods 18:366–368. https://doi.org/10.1038/s41592-021-01101-x

Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL (2009) BLAST+: architecture and applications. BMC Bioinformatics 10:421. https://doi.org/10.1186/1471-2105-10-421

Corke H, Avivi N, Atsmon D (1989) Pre- and post-anthesis accumulation of dry matter and nitrogen in wild barley (*Hordeum spontaneum*) and in barley cultivars (*H. vulgare*) differing in final grain size and protein content. Euphytica 40:127–134. https://doi.org/10.1007/BF00023307

Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R, 1000 Genomes Project Analysis Group (2011) The variant call format and VCFtools. Bioinformatics 27:2156–2158. https://doi.org/10.1093/bioinformatics/btr330

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H (2021) Twelve years of SAMtools and BCFtools. GigaScience 10:giab008. https://doi.org/10.1093/gigascience/giab008

Dempewolf H, Baute G, Anderson J, Kilian B, Smith C, Guarino L (2017) Past and future use of wild relatives in crop breeding. Crop Sci 57:1070–1082. https://doi.org/10.2135/cropsci2016.10.0885

Farag MA, Xiao J, Abdallah HM (2022) Nutritional value of barley cereal and better opportunities for its processing as a value-added food: a comprehensive review. Crit Rev Food Sci Nutr 62:1092–1104. https://doi.org/10.1080/10408398.2020.1835817

Frichot E, François O (2015) LEA: an R package for landscape and ecological association studies. Methods Ecol Evol 6:925–929. https://doi.org/10.1111/2041-210X.12382

Friedman M, Atsmon D (1988) Comparison of grain composition and nutritional quality in wild barley (*Hordeum spontaneum*) and in a standard cultivar. J Agric Food Chem 36:1167–1172. https://doi.org/10.1021/jf00084a010

Garcia-Gimenez G, Russell J, Aubert MK, Fincher GB, Burton RA, Waugh R, Tucker MR, Houston K (2019) Barley grain (1,3;1,4)-β-glucan content: effects of transcript and sequence variation in genes encoding the corresponding synthase and endohydrolase enzymes. Sci Rep 9:17250. https://doi.org/10.1038/s41598-019-53798-8

GBD 2019 Risk Factors Collaborators (2020) Global burden of 87 risk factors in 204 countries and territories, 1990–2019: a systematic analysis for the Global Burden of Disease Study 2019. Lancet 396:1223–1249. https://doi.org/10.1016/S0140-6736(20)30752-2

Geng L, Li M, Xie S, Wu D, Ye L, Zhang G (2021) Identification of genetic loci and candidate genes related to β-glucan content in barley grain by genome-wide association study in International Barley Core Selected Collection. Mol Breed 41:6. https://doi.org/10.1007/s11032-020-01199-5

Girard A-L, Mounet F, Lemaire-Chamley M, Gaillard C, Elmorjani K, Vivancos J, Runavot J-L, Quemener B, Petit J, Germain V, Rothan C, Marion D, Bakan B (2012) Tomato GDSL1 is required for cutin deposition in the fruit cuticle. Plant Cell 24:3119–3134. https://doi.org/10.1105/tpc.112.101055

Godfray HCJ, Beddington JR, Crute IR, Haddad L, Lawrence D, Muir JF, Pretty J, Robinson S, Thomas SM, Toulmin C (2010) Food security: the challenge of feeding 9 billion people. Science 327:812–818. https://doi.org/10.1126/science.1185383

Harlan JR, Zohary D (1966) Distribution of wild wheats and barley. Science 153:1074–1080. https://doi.org/10.1126/science.153.3740.1074

Holland JB, Nyquist WE, Cervantes-Martínez CT (2003) Estimating and interpreting heritability for plant breeding: an update. In: Janick J (ed) Plant breeding reviews, vol 22. Wiley, New York<sup>(TODO 501)</sup>, pp 9–112. https://doi.org/10.1002/9780470650202.ch2

Holm S (1979) A simple sequentially rejective multiple test procedure. Scand J Stat 6:65–70

Houston K, Russell J, Schreiber M, Halpin C, Oakey H, Washington JM, Booth A, Shirley N, Burton RA, Fincher GB, Waugh R (2014) A genome wide association scan for (1,3;1,4)-β-glucan content in the grain of contemporary 2-row spring and winter barleys. BMC Genomics 15:907. https://doi.org/10.1186/1471-2164-15-907

Hrmova M, Zimmer J, Bulone V, Fincher GB (2024) Enzymes in 3D: synthesis, remodelling, and hydrolysis of cell wall (1,3;1,4)-β-glucans. Plant Physiol 194:33–50. https://doi.org/10.1093/plphys/kiad415

Huang K, Jahani M, Gouzy J et al (2023) The genomics of linkage drag in inbred lines of sunflower. Proc Natl Acad Sci USA 120:e2205783119. https://doi.org/10.1073/pnas.2205783119

Hübner S, Kantar MB (2021) Tapping diversity from the wild: from sampling to implementation. Front Plant Sci 12:626565. https://doi.org/10.3389/fpls.2021.626565

Hübner S, Höffken M, Oren E, Haseneyer G, Stein N, Graner A, Schmid K, Fridman E (2009) Strong correlation of wild barley (*Hordeum spontaneum*) population structure with temperature and precipitation variation. Mol Ecol 18:1523–1536. https://doi.org/10.1111/j.1365-294X.2009.04106.x

Hübner S, Bdolach E, Ein-Gedy S, Schmid KJ, Korol A, Fridman E (2013) Phenotypic landscapes: phenological patterns in wild and cultivated barley. J Evol Biol 26:163–174. https://doi.org/10.1111/jeb.12043

Irigoyen S, Karlsson PM, Kuruvilla J, Spetea C, Versaw WK (2011) The sink-specific plastidic phosphate transporter PHT4;2 influences starch accumulation and leaf size in *Arabidopsis*. Plant Physiol 157:1765–1777. https://doi.org/10.1104/pp.111.181925

Jayakodi M, Lu Q, Pidon H et al (2024) Structural variation in the pangenome of wild and domesticated barley. Nature 636:654–662. https://doi.org/10.1038/s41586-024-08187-1

Jones P, Binns D, Chang H-Y et al (2014) InterProScan 5: genome-scale protein function classification. Bioinformatics 30:1236–1240. https://doi.org/10.1093/bioinformatics/btu031

Kang HM, Sul JH, Service SK, Zaitlen NA, Kong S-Y, Freimer NB, Sabatti C, Eskin E (2010) Variance component model to account for sample structure in genome-wide association studies. Nat Genet 42:348–354. https://doi.org/10.1038/ng.548

Keeling PL, Bacon PJ, Holt DC (1993) Elevated temperature reduces starch deposition in wheat endosperm by reducing the activity of soluble starch synthase. Planta 191:342–348<sup>(TODO 502)</sup>. https://doi.org/10.1007/BF00195691

König P, Beier S, Mascher M, Stein N, Lange M, Scholz U (2023) DivBrowse—interactive visualization and exploratory data analysis of variant call matrices. GigaScience 12:giad025. https://doi.org/10.1093/gigascience/giad025

Kruskal WH, Wallis WA (1952) Use of ranks in one-criterion variance analysis. J Am Stat Assoc 47:583–621. https://doi.org/10.1080/01621459.1952.10483441

Lim WL, Collins HM, Byrt CS, Lahnstein J, Shirley NJ, Aubert MK, Tucker MR, Peukert M, Matros A, Burton RA (2020) Overexpression of *HvCslF6* in barley grain alters carbohydrate partitioning plus transfer tissue and endosperm development. J Exp Bot 71:138–153. https://doi.org/10.1093/jxb/erz407

Ma B, Zhang L, Gao Q, Wang J, Li X, Wang H, Liu Y, Lin H, Liu J, Wang X, Li Q, Deng Y, Tang W, Luan S, He Z (2021) A plasma membrane transporter coordinates phosphate reallocation and grain filling in cereals. Nat Genet 53:906–915. https://doi.org/10.1038/s41588-021-00855-6

Marsh JI, Petereit J, Johnston BA, Bayer PE, Tay Fernandez CG, Al-Mamun HA, Batley J, Edwards D (2023) crosshap: R package for local haplotype visualization for trait association analysis. Bioinformatics 39:btad518. https://doi.org/10.1093/bioinformatics/btad518

Mascher M, Wicker T, Jenkins J et al (2021) Long-read sequence assembly: a technical evaluation in barley. Plant Cell 33:1888–1906. https://doi.org/10.1093/plcell/koab077

McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res 20:1297–1303. https://doi.org/10.1101/gr.107524.110

Nevo E, Chen G (2010) Drought and salt tolerances in wild relatives for wheat and barley improvement. Plant Cell Environ 33:670–685. https://doi.org/10.1111/j.1365-3040.2009.02107.x

Pedersen TL (2025) patchwork: the composer of plots. R package version 1.3.2. https://doi.org/10.32614/CRAN.package.patchwork

Petit J, Bres C, Mauxion J-P, Wong Jun Tai F, Martin LBB, Fich EA, Joubès J, Rose JKC, Domergue F, Rothan C (2016) The glycerol-3-phosphate acyltransferase GPAT6 from tomato plays a central role in fruit cutin biosynthesis. Plant Physiol 171:894–913. https://doi.org/10.1104/pp.16.00409

Potapenko E, Shermeister B, Mandel T, Hübner S (2026a) Intraspecific genome size variation is attributed to adaptive silencing of transposable elements in *Hordeum* species. Mol Biol Evol 43:msag051. https://doi.org/10.1093/molbev/msag051

Potapenko E, Schwarts D, Mandel T, Ashkenazy N, Fuerst D, Atsmon G, Snir S, Korol AB, Kantar MB, Bar-Massada A, Hübner S (2026b) Haplotype blocks are associated with rapid local adaptation to environmental shifts in wild barley. Mol Ecol 35:e70540. https://doi.org/10.1111/mec.70540

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ, Sham PC (2007) PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet 81:559–575. https://doi.org/10.1086/519795

Quinlan AR, Hall IM (2010) BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26:841–842. https://doi.org/10.1093/bioinformatics/btq033

R Core Team (2021) R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna. https://www.R-project.org/

Savin R, Nicolas ME (1996) Effects of short periods of drought and high temperature on grain growth and starch accumulation of two malting barley cultivars. Aust J Plant Physiol 23:201–210. https://doi.org/10.1071/PP9960201

Shu X, Rasmussen SK (2014) Quantification of amylose, amylopectin, and β-glucan in search for genes controlling the three major quality traits in barley by genome-wide association studies. Front Plant Sci 5:197. https://doi.org/10.3389/fpls.2014.00197

Simmonds NW (1995) The relation between yield and protein in cereal grain. J Sci Food Agric 67:309–315. https://doi.org/10.1002/jsfa.2740670306

Tomczak M, Tomczak E (2014) The need to report effect size estimates revisited. An overview of some recommended measures of effect size. Trends Sport Sci 21:19–25

UniProt Consortium (2025) UniProt: the Universal Protein Knowledgebase in 2025. Nucleic Acids Res 53:D609–D617. https://doi.org/10.1093/nar/gkae1010

Van der Auwera GA, O'Connor BD (2020) Genomics in the cloud: using Docker, GATK, and WDL in Terra. O'Reilly Media, Sebastopol

Vasimuddin M, Misra S, Li H, Aluru S (2019) Efficient architecture-aware acceleration of BWA-MEM for multicore systems. In: 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS), Rio de Janeiro, Brazil, pp 314–324. https://doi.org/10.1109/IPDPS.2019.00041

Wickham H (2016) ggplot2: elegant graphics for data analysis. Springer, New York. https://doi.org/10.1007/978-3-319-24277-4

Wilcoxon F (1945) Individual comparisons by ranking methods. Biometrics Bull 1:80–83. https://doi.org/10.2307/3001968

Xu Y, Zhang X-Q, Harasymow S, Westcott S, Zhang W, Li C (2018) Molecular marker-assisted backcrossing breeding: an example to transfer a thermostable β-amylase gene from wild barley. Mol Breed 38:63. https://doi.org/10.1007/s11032-018-0828-8

Yang W, Pollard M, Li-Beisson Y, Beisson F, Feig M, Ohlrogge J (2010) A distinct type of glycerol-3-phosphate acyltransferase with *sn*-2 preference and phosphatase activity producing 2-monoacylglycerol. Proc Natl Acad Sci USA 107:12040–12045. https://doi.org/10.1073/pnas.0914149107

Yates AD, Austine-Orimoloye O, Azov AG et al (2026) Ensembl 2026. Nucleic Acids Res 54:D1053–D1060. https://doi.org/10.1093/nar/gkaf1239

Zhang L, Gao C, Mentink-Vigier F, Tang L, Zhang D, Wang S, Cao S, Xu Z, Liu X, Wang T, Zhou Y, Zhang B (2019) Arabinosyl deacetylase modulates the arabinoxylan acetylation profile and secondary wall formation. Plant Cell 31:1113–1126. https://doi.org/10.1105/tpc.18.00894



# Supplementary Information

---

**Your decision:** approve · discard · or tell me what to change.
