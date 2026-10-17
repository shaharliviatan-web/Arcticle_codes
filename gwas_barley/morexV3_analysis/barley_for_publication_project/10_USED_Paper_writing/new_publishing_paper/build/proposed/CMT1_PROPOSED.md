# Review item CMT1: proposed changes

Proposal only. `Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item CMT1** (user request, 2026-10-01): TODO comments 201, 202, 205, 206, 207, 208, 209, 210 and 211 removed from the manuscript; their text is kept word for word, still unresolved, in `build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md`. The highlighted text they were attached to stays, as plain text. TODO 401 shortened to the user's wording and kept.
>

## Change 1 of 10 · Plant material and genotyping (TODO 201)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil (Online Resource 1). At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through <del style="color:#c0392b;background:#fdecea">three rounds of single-seed descent and selfing<sup>(TODO 201)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">three rounds of single-seed descent and selfing</ins>.

**Reads after approval:**

> We used a collection of 300 wild barley (*Hordeum spontaneum*) accessions sampled across Israel under a design that decouples environmental from geographic distance, as described in detail by Potapenko et al. (2026b). Briefly, the sampling area was divided into six regions, three corresponding to the northern, coastal and desert ecotypes of wild barley and three intermediate zones between them, and in each region five sites were selected to represent its climate and soil (Online Resource 1). At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through three rounds of single-seed descent and selfing.

---

## Change 2 of 10 · Plant material and genotyping (TODO 202)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026b). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed <del style="color:#c0392b;background:#fdecea">as described by Potapenko et al. (2026b)<sup>(TODO 202)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">as described by Potapenko et al. (2026b)</ins>. Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow (McKenna et al. 2010; Van der Auwera and O'Connor 2020). For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, with at most two alternate alleles per site. To reduce the computational load of the large genome, calling was run separately in non-overlapping 60-Mb intervals.

**Reads after approval:**

> Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026b). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed as described by Potapenko et al. (2026b). Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow (McKenna et al. 2010; Van der Auwera and O'Connor 2020). For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, with at most two alternate alleles per site. To reduce the computational load of the large genome, calling was run separately in non-overlapping 60-Mb intervals.

---

## Change 3 of 10 · Plant material and genotyping (TODO 205)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or <del style="color:#c0392b;background:#fdecea">more than 30% missing data were removed<sup>(TODO 205)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">more than 30% missing data were removed</ins>. Applied to all 300 accessions, these filters retained 7,110,996 SNPs. The ten accessions of site 04 (Hachola) were then removed because of an error in the sampling site<sup>(TODO 204)</sup>; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions.

**Reads after approval:**

> Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or more than 30% missing data were removed. Applied to all 300 accessions, these filters retained 7,110,996 SNPs. The ten accessions of site 04 (Hachola) were then removed because of an error in the sampling site<sup>(TODO 204)</sup>; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions.

---

## Change 4 of 10 · Common-garden experiment and phenotyping (TODO 206)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

The accessions were grown in a common-garden experiment over three growing seasons (2019–20, 2020–21 and 2021–22) in an insect-proof net-house in the Hula Valley, Israel (33°09′08.4″ N 35°37′15.8″ E), in 5-L pots under full irrigation, <del style="color:#c0392b;background:#fdecea">as described by Potapenko et al. (2026b)<sup>(TODO 206)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">as described by Potapenko et al. (2026b)</ins>. Each season was laid out as a randomized complete block design with one plant per accession in each block: four blocks in 2019–20 and 2020–21, and five in 2021–22.

**Reads after approval:**

> The accessions were grown in a common-garden experiment over three growing seasons (2019–20, 2020–21 and 2021–22) in an insect-proof net-house in the Hula Valley, Israel (33°09′08.4″ N 35°37′15.8″ E), in 5-L pots under full irrigation, as described by Potapenko et al. (2026b). Each season was laid out as a randomized complete block design with one plant per accession in each block: four blocks in 2019–20 and 2020–21, and five in 2021–22.

---

## Change 5 of 10 · Common-garden experiment and phenotyping (TODO 207)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): <del style="color:#c0392b;background:#fdecea">flowering time (days from sowing to heading)<sup>(TODO 207)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">flowering time (days from sowing to heading)</ins>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading), plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 6 of 10 · Common-garden experiment and phenotyping (TODO 208)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, <del style="color:#c0392b;background:#fdecea">plant height (length of the tallest tiller)<sup>(TODO 208)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">plant height (length of the tallest tiller)</ins>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller), spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 7 of 10 · Common-garden experiment and phenotyping (TODO 209)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, <del style="color:#c0392b;background:#fdecea">grain weight<sup>(TODO 209)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">grain weight</ins> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 8 of 10 · Common-garden experiment and phenotyping (TODO 210)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and <del style="color:#c0392b;background:#fdecea">number of grains per spike<sup>(TODO 210)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">number of grains per spike</ins>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four nutritional traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 9 of 10 · Common-garden experiment and phenotyping (TODO 211)

*Why:* User request (2026-10-01): comment removed from the manuscript; kept, unresolved, in build/UNRESOLVED_COMMENTS_OUTSIDE_MANUSCRIPT.md. The anchored text stays, without highlight.

**With the changes marked:**

At the end of each season, five spikes per plant were bagged and oven-dried at 35 °C for 30 days. The dried spikes were <del style="color:#c0392b;background:#fdecea">threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22<sup>(TODO 211)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22</ins>, and the grain was cleaned with a South Dakota seed blower (Seedburo Equipment Company) and stored for nutritional analysis.

**Reads after approval:**

> At the end of each season, five spikes per plant were bagged and oven-dried at 35 °C for 30 days. The dried spikes were threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22, and the grain was cleaned with a South Dakota seed blower (Seedburo Equipment Company) and stored for nutritional analysis.

---

## Change 10 of 10 · A shared carbon-allocation axis for starch, fiber and β-glucan (TODO 401, comment text only)

*Why:* User request: shortened to the user's wording; the comment stays in the manuscript.

**With the changes marked:**

Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; Online Resource 13), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

**Reads after approval:**

> Wild barley grain composition reveals two partly distinct biological dimensions: an independent protein axis, and a coordinated carbohydrate axis defining an inverse relationship between storage starch and the cell-wall polysaccharides, fiber and β-glucan. This trade-off was consistently observed across genotypes and along the north-to-desert gradient (Fig. 1c, d; Online Resource 13), indicating that grain composition reflects the competitive allocation of a common sucrose supply within the developing endosperm. This mechanism is supported by functional evidence showing that endosperm-specific overexpression of *HvCslF6* increases β-glucan accumulation while reducing starch content (Lim et al. 2020), demonstrating that these components are governed by interconnected physiological and genetic systems. The carbohydrate axis further integrates with a broader whole-plant resource-allocation strategy linked to plant architecture and phenology: larger, later-flowering genotypes with heavier grains prioritize storage starch, whereas smaller, faster-developing plants distribute resources toward higher grain numbers and structural cell-wall material. While domestication and breeding have historically shifted cultivated barley toward the starch-rich end of this axis, with higher starch and lower protein and fiber, wild barley remains a potential source of genetic variation for enriching the crop in fiber and β-glucan (Friedman and Atsmon 1988)<sup>(TODO 401)</sup>.

---

**Your decision:** approve · discard · or tell me what to change.
