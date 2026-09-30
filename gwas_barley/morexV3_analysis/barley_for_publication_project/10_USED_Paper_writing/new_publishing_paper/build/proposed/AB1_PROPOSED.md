# Review item AB1: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item AB1: title, Key message, Abstract, Keywords and Introduction**, drafted from the mini paper (Sariel's text) and updated to this study. They go at the top of the working file, before the M&M, in TAG order (Q1). Every number comes from the Results of this file.
>
> **TAG requirements checked:** the Key message is 27 words (limit 30) and adds to the title (the trade-off and the cultivars). The Abstract is 238 words (150–250), with no references and no undefined abbreviations. There are 6 Keywords (4–6). **Not here:** the title page (authors, affiliations, the corresponding author's ORCID, Acknowledgments) and the author-contribution statement. Those belong to the Declarations session.
>
> **Title:** Sariel's mini-paper title, still accurate. The species name follows the paper (*Hordeum spontaneum*, as in the M&M). See Q2.
>
> **Abstract:** Sariel's first and last sentences are kept in substance; the middle is updated. The loci count changes from 20 to 36, and the genes change from Pho, PHT, AP2/ERF and BAHD to *GPAT6*, *GH17*, *PHT4;3* and the GDSL gene. Rare alleles, loci vs heritability, the 7H shared signal and the cultivars are added.
>
> **Introduction:**
> - Sariel's five paragraphs, with his wording wherever the facts hold. His one-sentence opening paragraph merges into the next.
> - Small fixes: "genetic drag" → "linkage drag" (as in the Conclusions); "I analyzed" → "we"; grammar; TAG citation style.
> - **Added**, in Evgeny's style, because the mini paper lacks it: one sentence on why the haplotypes of the genes within a locus are analyzed; the sampling design and the earlier studies of this collection (Potapenko et al. 2026a local adaptation, 2026b genome size); objective (iv), the elite cultivars.
> - **QTL is defined here** ("quantitative trait loci (QTLs)", Sariel's wording), because the Ch. 2 heading keeps "QTLs".
> - **Hypothesis and objectives** (the mini paper's separate section) become the closing paragraph, in the past tense, with the objectives matched to what the paper did.
>
> **Word comments added:**
> - TODO 101: FAO citation format and the "fourth cereal" check.
> - TODO 102: verify Friedman and Atsmon's percentages in the PDF.
> - TODO 103: whether the crosshap paper is the right citation for the haplotype sentence.
>
> The other citations come from the mini paper's reference list; their DOIs go in in the References session (listed in a hidden note).
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** **Q1, placement:** put these sections into the working file, before the M&M (as proposed; one file for the whole manuscript, as you asked)? If yes, the file name `Methods_Results_Discussion_Conclusions.md` no longer describes it. Rename it (e.g. `Manuscript.md`, updating the review tool and CLAUDE.md), or keep the name?
>
> **Q2.** **Q2, title:** keep Sariel's title (default), or one that names the main result, e.g. "Genome-wide association and haplotype analysis link a starch–fiber trade-off to candidate genes for grain composition in wild barley (*Hordeum spontaneum*)"?
>
> **Q3.** **Q3, keywords:** OK with the six (Wild barley · Crop wild relatives · Grain nutritional quality · β-glucan · Genome-wide association · Haplotype analysis)?
>

---

## Change 1 of 1 · Top of the manuscript (before Materials and methods)

*Why:* New: title, Key message, Abstract, Keywords, Introduction (TAG order).

**With the changes marked:**

<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">

**Genome-wide association and haplotype analysis identify candidate genes for grain nutritional quality in wild barley (*Hordeum spontaneum*)**

# Key message

Genome-wide association and haplotype analysis in 290 wild barley accessions link a starch–fiber trade-off to shared loci and identify high-fiber wild haplotypes absent from elite malting cultivars.

# Abstract

Wild barley (*Hordeum spontaneum*), the progenitor of cultivated barley, retains genetic diversity that was reduced in cultivated germplasm through domestication and breeding, including variation in grain nutritional quality. Here, 290 wild barley accessions sampled across the environmental gradients of Israel were grown in a common garden over three seasons and phenotyped for grain protein, starch, β-glucan and fiber by near-infrared spectroscopy. Starch had the highest broad-sense heritability of the four traits, and β-glucan the largest genotype-by-season component. Starch was negatively correlated with fiber and β-glucan, and these traits followed opposite gradients of soil texture and precipitation at the sites of origin, suggesting a carbon-allocation trade-off between storage starch and cell-wall components, whereas protein varied largely independently. Genome-wide association mapping with 7.1 million single-nucleotide polymorphisms identified 36 loci, most of them carried by rare alleles, and the number of loci did not follow heritability. Haplotype analysis of the 55 genes within these loci, followed by functional annotation, pointed to a glycerol-3-phosphate acyltransferase (*GPAT6*) and a glycoside hydrolase family 17 glucanase (*GH17*) for fiber, and to a plastid phosphate transporter (*PHT4;3*) for starch. The tightest fiber–starch co-localization, on chromosome 7H, resolved to a GDSL esterase/lipase gene, whose haplotype, defined by variants in its putative promoter, was associated with higher fiber and lower starch. The high-fiber wild haplotypes of *GPAT6* and the GDSL gene were absent from five elite malting cultivars, providing targets for improving the nutritional value of cultivated barley. 

**Keywords** Wild barley · Crop wild relatives · Grain nutritional quality · β-glucan · Genome-wide association · Haplotype analysis 

# Introduction

Climate change and population growth are intensifying pressure on the global food supply, making food security an urgent priority (Godfray et al. 2010). For decades, crop breeding has targeted yield, often at the expense of grain nutritional quality and of the diversity available for further improvement (Dempewolf et al. 2017). Higher grain yield is also frequently associated with lower grain protein concentration, partly because carbohydrate and dry-matter accumulation can outpace nitrogen accumulation during grain filling (Simmonds 1995; Bogard et al. 2010). Dietary risks remain major contributors to the global burden of disease (GBD 2019 Risk Factors Collaborators 2020). Improving the nutritional quality of cereal grains is therefore an important breeding target.

Crop wild relatives (CWR) are a rich reservoir of beneficial alleles for improving disease resistance, abiotic-stress tolerance and yield in major crops, yet their potential for improving nutritional value remains largely untapped (Dempewolf et al. 2017). A principal obstacle is linkage drag, whereby an introgressed wild allele carries linked deleterious variants that degrade agronomic performance (Hübner and Kantar 2021; Huang et al. 2023). High-resolution mapping mitigates this by pinpointing quantitative trait loci (QTLs) at fine resolution, so that they can be introgressed while minimizing linkage drag (Hübner and Kantar 2021). The identification of beneficial alleles at high resolution is therefore essential for exploiting CWR in breeding.

Barley (*Hordeum vulgare*) was domesticated in the Fertile Crescent about 10,000 years ago and is now the fourth most important cereal worldwide (Badr et al. 2000; FAO 2025)<sup>(TODO 101)</sup>. Its wild progenitor, *H. vulgare* ssp. *spontaneum* (hereafter *H. spontaneum*), is distributed across a wide range of environments in the Fertile Crescent and adjacent regions (Harlan and Zohary 1966). Wild barley is characterized by broad genetic and phenotypic variation and has been used successfully in breeding for disease resistance and for drought and salt tolerance (Nevo and Chen 2010). The grain is dominated by starch, with protein, dietary fiber and the soluble fiber β-glucan as the other major components (Farag et al. 2022). Crucially, wild barley grain is nutritionally superior to that of cultivated barley, with on average about 50% more protein, 38% less starch and 65% more fiber<sup>(TODO 102)</sup> (Friedman and Atsmon 1988), making wild barley a potential genetic resource for improving grain quality in cultivated barley.

Advances in genomics now allow the genetic basis of target traits to be identified at high resolution. Genome-wide association studies (GWAS) are widely used to link genetic to phenotypic variation and have identified trait-associated loci and candidate genes in barley and many other crops (Alqudah et al. 2020). An associated single-nucleotide polymorphism (SNP), however, rarely identifies the causal gene, because an associated locus often spans several genes in linkage disequilibrium; analyzing the haplotypes of the genes within a locus can help resolve an association to candidate genes (Marsh et al. 2023)<sup>(TODO 103)</sup>. To dissect the genetic basis of grain nutritional traits in wild barley, we analyzed a collection of wild barley sampled across a wide range of environments in Israel, under a design that decouples environmental from geographic distance. This collection has been characterized genetically and phenotypically and used to study local adaptation (Potapenko et al. 2026a) and genome-size variation (Potapenko et al. 2026b), making it suitable for identifying candidate genes for key nutritional traits.

We hypothesized that the nutritional composition of wild barley grain is shaped by ecological constraints, so heritable variation is expected among populations from contrasting environments. To test this hypothesis, we (i) characterized the variation in grain nutritional traits among wild barley populations representing the different ecotypes of Israel; (ii) examined the correlations of these traits with the environmental gradients of the sites of origin and with plant morphology and phenology; (iii) identified loci and candidate genes associated with the nutritional traits by genome-wide association and haplotype analysis; and (iv) examined whether the wild haplotypes associated with these traits are present in elite cultivars. 

</ins># Materials and methods

**Reads after approval:**

> **Genome-wide association and haplotype analysis identify candidate genes for grain nutritional quality in wild barley (*Hordeum spontaneum*)**

# Key message

Genome-wide association and haplotype analysis in 290 wild barley accessions link a starch–fiber trade-off to shared loci and identify high-fiber wild haplotypes absent from elite malting cultivars.

# Abstract

Wild barley (*Hordeum spontaneum*), the progenitor of cultivated barley, retains genetic diversity that was reduced in cultivated germplasm through domestication and breeding, including variation in grain nutritional quality. Here, 290 wild barley accessions sampled across the environmental gradients of Israel were grown in a common garden over three seasons and phenotyped for grain protein, starch, β-glucan and fiber by near-infrared spectroscopy. Starch had the highest broad-sense heritability of the four traits, and β-glucan the largest genotype-by-season component. Starch was negatively correlated with fiber and β-glucan, and these traits followed opposite gradients of soil texture and precipitation at the sites of origin, suggesting a carbon-allocation trade-off between storage starch and cell-wall components, whereas protein varied largely independently. Genome-wide association mapping with 7.1 million single-nucleotide polymorphisms identified 36 loci, most of them carried by rare alleles, and the number of loci did not follow heritability. Haplotype analysis of the 55 genes within these loci, followed by functional annotation, pointed to a glycerol-3-phosphate acyltransferase (*GPAT6*) and a glycoside hydrolase family 17 glucanase (*GH17*) for fiber, and to a plastid phosphate transporter (*PHT4;3*) for starch. The tightest fiber–starch co-localization, on chromosome 7H, resolved to a GDSL esterase/lipase gene, whose haplotype, defined by variants in its putative promoter, was associated with higher fiber and lower starch. The high-fiber wild haplotypes of *GPAT6* and the GDSL gene were absent from five elite malting cultivars, providing targets for improving the nutritional value of cultivated barley. 

**Keywords** Wild barley · Crop wild relatives · Grain nutritional quality · β-glucan · Genome-wide association · Haplotype analysis 

# Introduction

Climate change and population growth are intensifying pressure on the global food supply, making food security an urgent priority (Godfray et al. 2010). For decades, crop breeding has targeted yield, often at the expense of grain nutritional quality and of the diversity available for further improvement (Dempewolf et al. 2017). Higher grain yield is also frequently associated with lower grain protein concentration, partly because carbohydrate and dry-matter accumulation can outpace nitrogen accumulation during grain filling (Simmonds 1995; Bogard et al. 2010). Dietary risks remain major contributors to the global burden of disease (GBD 2019 Risk Factors Collaborators 2020). Improving the nutritional quality of cereal grains is therefore an important breeding target.

Crop wild relatives (CWR) are a rich reservoir of beneficial alleles for improving disease resistance, abiotic-stress tolerance and yield in major crops, yet their potential for improving nutritional value remains largely untapped (Dempewolf et al. 2017). A principal obstacle is linkage drag, whereby an introgressed wild allele carries linked deleterious variants that degrade agronomic performance (Hübner and Kantar 2021; Huang et al. 2023). High-resolution mapping mitigates this by pinpointing quantitative trait loci (QTLs) at fine resolution, so that they can be introgressed while minimizing linkage drag (Hübner and Kantar 2021). The identification of beneficial alleles at high resolution is therefore essential for exploiting CWR in breeding.

Barley (*Hordeum vulgare*) was domesticated in the Fertile Crescent about 10,000 years ago and is now the fourth most important cereal worldwide (Badr et al. 2000; FAO 2025)<sup>(TODO 101)</sup>. Its wild progenitor, *H. vulgare* ssp. *spontaneum* (hereafter *H. spontaneum*), is distributed across a wide range of environments in the Fertile Crescent and adjacent regions (Harlan and Zohary 1966). Wild barley is characterized by broad genetic and phenotypic variation and has been used successfully in breeding for disease resistance and for drought and salt tolerance (Nevo and Chen 2010). The grain is dominated by starch, with protein, dietary fiber and the soluble fiber β-glucan as the other major components (Farag et al. 2022). Crucially, wild barley grain is nutritionally superior to that of cultivated barley, with on average about 50% more protein, 38% less starch and 65% more fiber<sup>(TODO 102)</sup> (Friedman and Atsmon 1988), making wild barley a potential genetic resource for improving grain quality in cultivated barley.

Advances in genomics now allow the genetic basis of target traits to be identified at high resolution. Genome-wide association studies (GWAS) are widely used to link genetic to phenotypic variation and have identified trait-associated loci and candidate genes in barley and many other crops (Alqudah et al. 2020). An associated single-nucleotide polymorphism (SNP), however, rarely identifies the causal gene, because an associated locus often spans several genes in linkage disequilibrium; analyzing the haplotypes of the genes within a locus can help resolve an association to candidate genes (Marsh et al. 2023)<sup>(TODO 103)</sup>. To dissect the genetic basis of grain nutritional traits in wild barley, we analyzed a collection of wild barley sampled across a wide range of environments in Israel, under a design that decouples environmental from geographic distance. This collection has been characterized genetically and phenotypically and used to study local adaptation (Potapenko et al. 2026a) and genome-size variation (Potapenko et al. 2026b), making it suitable for identifying candidate genes for key nutritional traits.

We hypothesized that the nutritional composition of wild barley grain is shaped by ecological constraints, so heritable variation is expected among populations from contrasting environments. To test this hypothesis, we (i) characterized the variation in grain nutritional traits among wild barley populations representing the different ecotypes of Israel; (ii) examined the correlations of these traits with the environmental gradients of the sites of origin and with plant morphology and phenology; (iii) identified loci and candidate genes associated with the nutritional traits by genome-wide association and haplotype analysis; and (iv) examined whether the wild haplotypes associated with these traits are present in elite cultivars. 

# Materials and methods

---

**Your decision:** approve · discard · or tell me what to change.
