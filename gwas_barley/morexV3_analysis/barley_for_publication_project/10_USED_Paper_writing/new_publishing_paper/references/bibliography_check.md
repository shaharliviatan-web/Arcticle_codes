# Bibliographic check (step 2)

Every cited work checked against Crossref on 2026-10-01 (`scripts/01_crossref_check.py`; raw records in
`meta/crossref/`). The 7 works without a Crossref record were checked by hand (sources in the table). This step
checks only that each reference exists and has the right authors, year, journal, volume, pages and DOI. Whether
it supports the sentence is the audit (`citation_audit.md`).

**Result:** all 66 works exist (65 cited, plus the candidate for the missing β-amylase citation). The details the
TODO comments drafted from memory were right, except for the findings below.

## Findings for the manuscript

Findings 1–4 were approved by the user and **applied on 2026-10-01 as review item REF1** (log at the end of
`build/REVIEW_2026-09-27_Results_Discussion.md`).

1. ✅ **Hrmova and Fincher 2024 → Hrmova et al. 2024.** Four authors (Hrmova, Zimmer, Bulone, Fincher). 2 places.
2. ✅ **Yates et al. 2025 → 2026.** *Ensembl 2026*, *Nucleic Acids Res* 54(D1), January 2026 issue (online
   2025-11-26); the issue year, as for the other references. 1 place.
3. ✅ **Potapenko letters swapped to TAG's chronological order.** 2026a = *Mol Biol Evol* (10 March 2026), 2026b =
   *Mol Ecol* (September 2026). 11 places in the text, plus the TODO anchors and hidden notes. In this folder the two
   papers have letter-free keys: `PotapenkoMBE2026`, `PotapenkoMolEcol2026`.
4. ✅ **(Zhang et al. 2019; Girard et al. 2012) → (Girard et al. 2012; Zhang et al. 2019).** 2 places.
5. **Holland et al. 2003:** Crossref dates *Plant Breeding Reviews* vol. 22 to 30 October 2002; it is usually cited as
   2003. The chapter could not be obtained, so the year stays 2003 and the check is in Word comment 257 (REF2).
6. **Added in REF2:** AACC 2001 (no DOI; free PDF held) and Keeling et al. 1993 (Planta 191(3):342-348; pages from
   the publisher listing, Crossref has none).

## All works


<!-- AUTO:table start -->
| key | in the text | record found (first author, title, journal volume:pages) | DOI | note |
|---|---|---|---|---|
| Alqudah2020 | Alqudah et al. 2020 | Alqudah et al — GWAS: Fast-forwarding gene identification and characterization in temperate Cereals: lessons from Barley – A review — Journal of Advanced Research 22:119-135 | https://doi.org/10.1016/j.jare.2019.10.013 | ✅ |
| AnkerNilssen2008 | Anker-Nilssen et al. 2008 | Anker-Nilssen et al — Influence of growth temperature on content, viscosity and relative molecular weight of water-soluble β-glucans in barley (Hordeum vulgare L.) — Journal of Cereal Science 48:670-677 | https://doi.org/10.1016/j.jcs.2008.02.009 | ✅ |
| Badr2000 | Badr et al. 2000 | Badr et al — On the Origin and Domestication History of Barley (Hordeum vulgare) — Molecular Biology and Evolution 17:499-510 | https://doi.org/10.1093/oxfordjournals.molbev.a026330 | ✅ |
| Bates2015 | Bates et al. 2015 | Bates et al — Fitting Linear Mixed-Effects Models Using lme4 — J. Stat. Soft. 67 | https://doi.org/10.18637/jss.v067.i01 | ✅ |
| Beljkas2010 | Beljkaš et al. 2010 | Beljkaš et al — Rapid method for determination of protein content in cereals and oilseeds: validation, measurement uncertainty and comparison with the Kjeldahl method — Accred Qual Assur 15:555-561 | https://doi.org/10.1007/s00769-010-0677-6 | ✅ |
| Benjamini1995 | Benjamini and Hochberg 1995 | Benjamini and Hochberg — Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing — Journal of the Royal Statistical Society Series B: Statistical Methodology 57:289-300 | https://doi.org/10.1111/j.2517-6161.1995.tb02031.x | ✅ |
| Bogard2010 | Bogard et al. 2010 | Bogard et al — Deviation from the grain protein concentration–grain yield negative relationship is highly correlated to post-anthesis N uptake in winter wheat — Journal of Experimental Botany 61:4303-4312 | https://doi.org/10.1093/jxb/erq238 | ✅ |
| BroadInstitute2019 | Broad Institute 2019 | Picard toolkit (2019). Broad Institute, GitHub repository. https://broadinstitute.github.io/picard/ | — | the citation the Picard README asks for |
| Buchfink2021 | Buchfink et al. 2021 | Buchfink et al — Sensitive protein alignments at tree-of-life scale using DIAMOND — Nat Methods 18:366-368 | https://doi.org/10.1038/s41592-021-01101-x | ✅ |
| Burton2008 | Burton et al. 2008 | Burton et al — The Genetics and Transcriptional Profiles of the Cellulose Synthase-Like HvCslF Gene Family in Barley — Plant Physiol. 146:1821-1833 | https://doi.org/10.1104/pp.107.114694 | ✅ |
| Camacho2009 | Camacho et al. 2009 | Camacho et al — BLAST+: architecture and applications — BMC Bioinformatics 10:421 | https://doi.org/10.1186/1471-2105-10-421 | ✅ |
| Corke1989 | Corke et al. 1989 | Corke et al — Pre- and post-anthesis accumulation of dry matter and nitrogen in wild barley (Hordeum spontaneum) and in barley cultivars (H. vulgare) differing in final grain size and protein content — Euphytica 40:127-134 | https://doi.org/10.1007/bf00023307 | ✅ |
| Danecek2011 | Danecek et al. 2011 | Danecek et al — The variant call format and VCFtools — Bioinformatics 27:2156-2158 | https://doi.org/10.1093/bioinformatics/btr330 | ✅ |
| Danecek2021 | Danecek et al. 2021 | Danecek et al — Twelve years of SAMtools and BCFtools — GigaScience 10:giab008 | https://doi.org/10.1093/gigascience/giab008 | ✅ |
| Dempewolf2017 | Dempewolf et al. 2017 | Dempewolf et al — Past and Future Use of Wild Relatives in Crop Breeding — Crop Science 57:1070-1082 | https://doi.org/10.2135/cropsci2016.10.0885 | ✅ |
| FAO2025 | FAO 2025 | FAOSTAT statistical database. FAO, Rome. https://www.fao.org/faostat/ (access date needed) | — | online database; the 4th-cereal statement to be checked in FAOSTAT (TODO 101) |
| Farag2022 | Farag et al. 2022 | Farag et al — Nutritional value of barley cereal and better opportunities for its processing as a value-added food: a comprehensive review — Critical Reviews in Food Science and Nutrition 62:1092-1104 | https://doi.org/10.1080/10408398.2020.1835817 | online 2020, issue 62(4) 2022: 2022 correct |
| Frichot2015 | Frichot and François 2015 | Frichot and François — LEA: An R package for landscape and ecological association studies — Methods Ecol Evol 6:925-929 | https://doi.org/10.1111/2041-210x.12382 | ✅ |
| Friedman1988 | Friedman and Atsmon 1988 | Friedman and Atsmon — Comparison of grain composition and nutritional quality in wild barley (Hordeum spontaneum) and in a standard cultivar — J. Agric. Food Chem. 36:1167-1172 | https://doi.org/10.1021/jf00084a010 | ✅ |
| GBD2020 | GBD 2019 Risk Factors Collaborators 2020 | Murray et al — Global burden of 87 risk factors in 204 countries and territories, 1990–2019: a systematic analysis for the Global Burden of Disease Study 2019 — The Lancet 396:1223-1249 | https://doi.org/10.1016/s0140-6736(20)30752-2 | corporate author; 940 named authors |
| GarciaGimenez2019 | Garcia-Gimenez et al. 2019 | Garcia-Gimenez et al — Barley grain (1,3;1,4)-β-glucan content: effects of transcript and sequence variation in genes encoding the corresponding synthase and endohydrolase enzymes — Sci Rep 9:17250 | https://doi.org/10.1038/s41598-019-53798-8 | ✅ |
| Geng2021 | Geng et al. 2021 | Geng et al — Identification of genetic loci and candidate genes related to β-glucan content in barley grain by genome-wide association study in International Barley Core Selected Collection — Mol Breeding 41:6 | https://doi.org/10.1007/s11032-020-01199-5 | ✅ |
| Girard2012 | Girard et al. 2012 | Girard et al — Tomato GDSL1 Is Required for Cutin Deposition in the Fruit Cuticle — Plant Cell 24:3119-3134 | https://doi.org/10.1105/tpc.112.101055 | ✅ |
| Godfray2010 | Godfray et al. 2010 | Godfray et al — Food Security: The Challenge of Feeding 9 Billion People — Science 327:812-818 | https://doi.org/10.1126/science.1185383 | ✅ |
| Guo2008 | Guo et al. 2008 | Guo et al — Functional analysis of the Arabidopsis PHT4 family of intracellular phosphate transporters — New Phytologist 177:889-898 | https://doi.org/10.1111/j.1469-8137.2007.02331.x | online 2007, issue 177(4) 2008: 2008 correct |
| Harlan1966 | Harlan and Zohary 1966 | Harlan and Zohary — Distribution of Wild Wheats and Barley — Science 153:1074-1080 | https://doi.org/10.1126/science.153.3740.1074 | ✅ |
| Holland2003 | Holland et al. 2003 | Holland et al — Estimating and Interpreting Heritability for Plant Breeding: An Update — Plant Breeding Reviews 22:9-112 | https://doi.org/10.1002/9780470650202.ch2 | Crossref dates the volume 2002 (published 30 Oct 2002); usually cited 2003; no PDF obtainable, kept 2003 (Word comment 257) |
| Holm1979 | Holm 1979 | A simple sequentially rejective multiple test procedure. Scand J Stat 6:65-70 | — | no DOI (JSTOR 4615733; the 10.2307 DOI does not resolve) |
| Houston2014 | Houston et al. 2014 | Houston et al — A genome wide association scan for (1,3;1,4)-β-glucan content in the grain of contemporary 2-row Spring and Winter barleys — BMC Genomics 15:907 | https://doi.org/10.1186/1471-2164-15-907 | ✅ |
| Hrmova2024 | Hrmova et al. 2024 | Hrmova et al — Enzymes in 3D: Synthesis, remodelling, and hydrolysis of cell wall (1,3;1,4)-β-glucans — Plant Physiology 194:33-50 | https://doi.org/10.1093/plphys/kiad415 | 4 authors (Hrmova, Zimmer, Bulone, Fincher): was *Hrmova and Fincher 2024*, fixed in REF1 |
| Huang2023 | Huang et al. 2023 | Huang et al — The genomics of linkage drag in inbred lines of sunflower — Proc. Natl. Acad. Sci. U.S.A. 120:e2205783119 | https://doi.org/10.1073/pnas.2205783119 | PNAS 120(14):e2205783119 (the article number ends in 119, not 120) |
| HubnerKantar2021 | Hübner and Kantar 2021 | Hübner and Kantar — Tapping Diversity From the Wild: From Sampling to Implementation — Front. Plant Sci. 12:626565 | https://doi.org/10.3389/fpls.2021.626565 | ✅ |
| Hubner2009 | Hübner et al. 2009 | Hübner et al — Strong correlation of wild barley ( Hordeum spontaneum ) population structure with temperature and precipitation variation — Molecular Ecology 18:1523-1536 | https://doi.org/10.1111/j.1365-294x.2009.04106.x | ✅ |
| Hubner2013 | Hübner et al. 2013 | Hübner et al — Phenotypic landscapes: phenological patterns in wild and cultivated barley — J of Evolutionary Biology 26:163-174 | https://doi.org/10.1111/jeb.12043 | online 2012, issue 26(1) 2013: 2013 correct |
| Irigoyen2011 | Irigoyen et al. 2011 | Irigoyen et al — The Sink-Specific Plastidic Phosphate Transporter PHT4;2 Influences Starch Accumulation and Leaf Size in Arabidopsis  — Plant Physiology 157:1765-1777 | https://doi.org/10.1104/pp.111.181925 | ✅ |
| Jayakodi2024 | Jayakodi et al. 2024 | Jayakodi et al — Structural variation in the pangenome of wild and domesticated barley — Nature 636:654-662 | https://doi.org/10.1038/s41586-024-08187-1 | ✅ |
| Jones2014 | Jones et al. 2014 | Jones et al — InterProScan 5: genome-scale protein function classification — Bioinformatics 30:1236-1240 | https://doi.org/10.1093/bioinformatics/btu031 | ✅ |
| Kang2010 | Kang et al. 2010 | Kang et al — Variance component model to account for sample structure in genome-wide association studies — Nat Genet 42:348-354 | https://doi.org/10.1038/ng.548 | ✅ |
| Konig2023 | König et al. 2023 | König et al — DivBrowse—interactive visualization and exploratory data analysis of variant call matrices — GigaScience 12:giad025 | https://doi.org/10.1093/gigascience/giad025 | Crossref print date 2022-12-28, online 2023-04-21, volume 12 (2023): 2023 correct |
| Kruskal1952 | Kruskal and Wallis 1952 | Kruskal and Wallis — Use of Ranks in One-Criterion Variance Analysis — Journal of the American Statistical Association 47:583-621 | https://doi.org/10.1080/01621459.1952.10483441 | ✅ |
| Lim2020 | Lim et al. 2020 | Lim et al — Overexpression of HvCslF6 in barley grain alters carbohydrate partitioning plus transfer tissue and endosperm development — Journal of Experimental Botany 71:138-153 | https://doi.org/10.1093/jxb/erz407 | online 2019, issue 71(1) January 2020: 2020 correct |
| Ma2021 | Ma et al. 2021 | Ma et al — A plasma membrane transporter coordinates phosphate reallocation and grain filling in cereals — Nat Genet 53:906-915 | https://doi.org/10.1038/s41588-021-00855-6 | ✅ |
| Marsh2023 | Marsh et al. 2023 | Marsh et al — crosshap: R package for local haplotype visualization for trait association analysis — Bioinformatics 39:btad518 | https://doi.org/10.1093/bioinformatics/btad518 | ✅ |
| Mascher2021 | Mascher et al. 2021 | Mascher et al — Long-read sequence assembly: a technical evaluation in barley — The Plant Cell 33:1888-1906 | https://doi.org/10.1093/plcell/koab077 | ✅ |
| McKenna2010 | McKenna et al. 2010 | McKenna et al — The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data — Genome Res. 20:1297-1303 | https://doi.org/10.1101/gr.107524.110 | ✅ |
| Nevo2010 | Nevo and Chen 2010 | Nevo and Chen — Drought and salt tolerances in wild relatives for wheat and barley improvement — Plant Cell & Environment 33:670-685 | https://doi.org/10.1111/j.1365-3040.2009.02107.x | ✅ |
| Pedersen2025 | Pedersen 2025 | patchwork: the composer of plots. R package version 1.3.2. https://doi.org/10.32614/CRAN.package.patchwork | — | CRAN: v1.3.2 published 2025-08-25 |
| Petit2016 | Petit et al. 2016 | Petit et al — The glycerol-3-phosphate acyltransferase GPAT6 from tomato plays a central role in fruit cutin biosynthesis — Plant Physiol. 171:894-913 | https://doi.org/10.1104/pp.16.00409 | Crossref holds the early-online record; final pages 171(2):894-913 (Europe PMC) |
| PotapenkoMolEcol2026 | Potapenko et al. 2026b | Potapenko et al — Haplotype Blocks Are Associated With Rapid Local Adaptation to Environmental Shifts in Wild Barley — Molecular Ecology 35:e70540 | https://doi.org/10.1111/mec.70540 | Mol Ecol 35(17), September 2026 (online 2026-09-01): 2026b (was 2026a; letters swapped in REF1) |
| PotapenkoMBE2026 | Potapenko et al. 2026a | Potapenko et al — Intraspecific genome size variation is attributed to adaptive silencing of transposable elements in Hordeum species — Molecular Biology and Evolution 43:msag051 | https://doi.org/10.1093/molbev/msag051 | Mol Biol Evol 43(3), 10 March 2026: 2026a (was 2026b; letters swapped in REF1) |
| Purcell2007 | Purcell et al. 2007 | Purcell et al — PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses — The American Journal of Human Genetics 81:559-575 | https://doi.org/10.1086/519795 | ✅ |
| Quinlan2010 | Quinlan and Hall 2010 | Quinlan and Hall — BEDTools: a flexible suite of utilities for comparing genomic features — Bioinformatics 26:841-842 | https://doi.org/10.1093/bioinformatics/btq033 | ✅ |
| RCoreTeam2021 | R Core Team 2021 | R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna. https://www.R-project.org/ | — | R 4.1.2 was released 2021-11-01 |
| Savin1996 | Savin and Nicolas 1996 | Savin and Nicolas — Effects of Short Periods of Drought and High Temperature on Grain Growth and Starch Accumulation of Two Malting Barley Cultivars — Australian Journal of Plant Physiology 23:201-210 | https://doi.org/10.1071/pp9960201 | ✅ |
| Shu2014 | Shu and Rasmussen 2014 | Shu and Rasmussen — Quantification of amylose, amylopectin, and Î²-glucan in search for genes controlling the three major quality traits in barley by genome-wide association studies — Front. Plant Sci. 5 | https://doi.org/10.3389/fpls.2014.00197 | ✅ |
| Simmonds1995 | Simmonds 1995 | Simmonds — The relation between yield and protein in cereal grain — J Sci Food Agric 67:309-315 | https://doi.org/10.1002/jsfa.2740670306 | ✅ |
| Tomczak2014 | Tomczak and Tomczak 2014 | The need to report effect size estimates revisited. An overview of some recommended measures of effect size. Trends Sport Sci 21(1):19-25 | — | no DOI; PDF held (Wielkopolska Digital Library) |
| UniProt2025 | UniProt Consortium 2025 | The UniProt Consortium et al — UniProt: the Universal Protein Knowledgebase in 2025 — Nucleic Acids Research 53:D609-D617 | https://doi.org/10.1093/nar/gkae1010 | online 2024, issue 53(D1) January 2025: 2025 correct |
| VanderAuwera2020 | Van der Auwera and O'Connor 2020 | Genomics in the cloud: using Docker, GATK, and WDL in Terra. O'Reilly Media, Sebastopol. ISBN 978-1-491-97519-0 | — | Open Library record |
| Vasimuddin2019 | Vasimuddin et al. 2019 | Vasimuddin et al — Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems — 2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS) :314-324 | https://doi.org/10.1109/ipdps.2019.00041 | ✅ |
| Wickham2016 | Wickham 2016 | Wickham — ggplot2 — Use R!  | https://doi.org/10.1007/978-3-319-24277-4 | ✅ |
| Wilcoxon1945 | Wilcoxon 1945 | Wilcoxon — Individual Comparisons by Ranking Methods — Biometrics Bulletin 1:80 | https://doi.org/10.2307/3001968 | ✅ |
| Yang2010 | Yang et al. 2010 | Yang et al — A distinct type of glycerol-3-phosphate acyltransferase with sn -2 preference and phosphatase activity producing 2-monoacylglycerol — Proc. Natl. Acad. Sci. U.S.A. 107:12040-12045 | https://doi.org/10.1073/pnas.0914149107 | ✅ |
| Yates2026 | Yates et al. 2026 | Yates et al — Ensembl 2026 — Nucleic Acids Research 54:D1053-D1060 | https://doi.org/10.1093/nar/gkaf1239 | title *Ensembl 2026*; online 2025-11-26, issue 54(D1) January 2026: was *Yates et al. 2025*, fixed in REF1 |
| Zhang2019 | Zhang et al. 2019 | Zhang et al — Arabinosyl Deacetylase Modulates the Arabinoxylan Acetylation Profile and Secondary Wall Formation — Plant Cell 31:1113-1126 | https://doi.org/10.1105/tpc.18.00894 | ✅ |
| Xu2018 | [TODO 411 candidate] Xu et al. 2018 | Xu et al — Molecular marker-assisted backcrossing breeding: an example to transfer a thermostable β-amylase gene from wild barley — Mol Breeding 38:63 | https://doi.org/10.1007/s11032-018-0828-8 | candidate for the missing β-amylase citation (TODO 411), not yet in the text |
| AACC2001 | AACC 2001 | AACC (2001) The definition of dietary fiber. Report of the Dietary Fiber Definition Committee to the Board of Directors of the American Association of Cereal Chemists. Cereal Foods World 46:112-126. https://www.cerealsgrains.org/resources/definitions/Documents/DietaryFiber/DFDef.pdf | — | no DOI; free PDF held (pdfs/AACC2001.pdf); added in REF2 |
| Keeling1993 | Keeling et al. 1993 | Keeling et al — Elevated temperature reduces starch deposition in wheat endosperm by reducing the activity of soluble starch synthase — Planta 191:342-348 | https://doi.org/10.1007/bf00195691 | pages 342-348 from the publisher listing (Crossref has no pages); added in REF2 |
<!-- AUTO:table end -->
