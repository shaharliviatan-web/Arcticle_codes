# Review item MM8: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM8 of the M&M voice review** (item 6, voice and clarity, part 1 of 3): Plant material ¶2–3, Common garden ¶2 and Grain composition by NIR. No facts change; the TODO comments stay where they are.
>
> **Plant material ¶2:** "its GVCF was reblocked with ReblockGVCF" (the tool name repeats the verb) → "a GVCF file was produced … and condensed with ReblockGVCF". The 60-Mb intervals sentence now opens with its aim, as in Evgeny's "To reduce the load…".
>
> **Plant material ¶3 (B17):** "without re-filtering the SNPs, leaving 290 accessions and 7,110,996 SNPs" is stated as a plain fact. The wording of the Hachola reason still waits for TODO 204.
>
> **Common garden ¶2:** the list set off by commas mid-sentence ("The four grain-composition traits, protein, starch, β-glucan and fiber, were…") gets parentheses.
>
> **NIR:** the one-sentence ¶1 joins ¶2, and "measured in percentages on the cleaned grain … which scans the grain without destroying it" becomes "contents (%) … measured non-destructively". The calibration sentences now open with their aim and use "we". The stacked qualifiers on the kits are split (what was adapted, and how the mean was taken). The calibration text (TODO 213) is reworded; "calibration package" becomes "calibration".
>
> Unchanged on purpose: the sequencing and alignment sentences (already in Evgeny's form), the filter thresholds, the trait list, the post-harvest paragraph.
>
> Length: 313 → 312 words (the aim is clarity; the length stays about the same).
>

## Change 1 of 5 · M&M · Plant material and genotyping ¶2 (end)

*Why:* 'reblocked with ReblockGVCF'; aim-first interval sentence.

**With the changes marked:**

Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026a). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed as described by Potapenko et al. (2026a)<sup>(TODO 202)</sup>. Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow (McKenna et al. 2010; Van der Auwera and O'Connor 2020)<sup>(TODO 203)</sup>. <del style="color:#c0392b;background:#fdecea">Each accession was first called in GVCF mode with HaplotypeCaller (ploidy 2) and its GVCF was reblocked with ReblockGVCF</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF</ins>; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, <del style="color:#c0392b;background:#fdecea">allowing</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">with</ins> at most two alternate alleles per site. <del style="color:#c0392b;background:#fdecea">Because of the size of the genome, calling</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To reduce the computational load of the large genome, calling</ins> was run separately in non-overlapping 60-Mb intervals.

**Reads after approval:**

> Genomic DNA was extracted from a single seedling of each accession with a cetyltrimethylammonium bromide (CTAB) protocol, and NEBNext libraries were sequenced on an Illumina NovaSeq 6000 (2 × 150 bp), targeting a coverage of at least 5× (Potapenko et al. 2026a). The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed as described by Potapenko et al. (2026a)<sup>(TODO 202)</sup>. Cleaned reads were aligned to the barley cv. Morex reference genome assembly v3 (MorexV3; Mascher et al. 2021) with BWA-MEM2 v2.2.1 under default parameters (Vasimuddin et al. 2019), and duplicate reads were marked with Picard MarkDuplicates v3.4.0 (Broad Institute 2019). Variants were called with GATK v4.6.2.0 following the best-practices workflow (McKenna et al. 2010; Van der Auwera and O'Connor 2020)<sup>(TODO 203)</sup>. For each accession, a GVCF file was produced with HaplotypeCaller (ploidy 2) and condensed with ReblockGVCF; all accessions were then genotyped jointly with GenomicsDBImport and GenotypeGVCFs, with at most two alternate alleles per site. To reduce the computational load of the large genome, calling was run separately in non-overlapping 60-Mb intervals.

---

## Change 2 of 5 · M&M · Plant material and genotyping ¶3 (end)

*Why:* 'without re-filtering …, leaving …' stated as a plain fact.

**With the changes marked:**

Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or more than 30% missing data were removed<sup>(TODO 205)</sup>. <del style="color:#c0392b;background:#fdecea">These filters were applied to all 300 accessions and yielded 7,110,996 SNPs.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Applied to all 300 accessions, these filters retained 7,110,996 SNPs.</ins> The ten accessions of site 04 (Hachola) were then <del style="color:#c0392b;background:#fdecea">excluded</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">removed</ins> because of an error in the sampling site<sup>(TODO 204)</sup><del style="color:#c0392b;background:#fdecea">, without re-filtering the SNPs, leaving 290 accessions and 7,110,996 SNPs for all analyses</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions</ins>.

**Reads after approval:**

> Variants were filtered with bcftools v1.13 (Danecek et al. 2021) and VCFtools v0.1.15 (Danecek et al. 2011). Only biallelic single-nucleotide polymorphisms (SNPs) were kept, and SNPs were removed when QualByDepth < 5, MappingQuality < 45, FisherStrand > 60, StrandOddsRatio > 3, MappingQualityRankSum < −2.5, QUAL < 140, or total site read depth was below 900 or above 1,800. Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or more than 30% missing data were removed<sup>(TODO 205)</sup>. Applied to all 300 accessions, these filters retained 7,110,996 SNPs. The ten accessions of site 04 (Hachola) were then removed because of an error in the sampling site<sup>(TODO 204)</sup>; the SNPs were not re-filtered, so all analyses used the same 7,110,996 SNPs in 290 accessions.

---

## Change 3 of 5 · M&M · Common garden ¶2 (last sentence)

*Why:* Comma-appositive list → parentheses.

**With the changes marked:**

Six morphological and phenological traits were recorded as described by Potapenko et al. (2026a): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. <del style="color:#c0392b;background:#fdecea">The four grain-composition traits, protein, starch, β-glucan and fiber,</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Four grain-composition traits (protein, starch, β-glucan and fiber)</ins> were measured by near-infrared spectroscopy, as described below.

**Reads after approval:**

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026a): flowering time (days from sowing to heading)<sup>(TODO 207)</sup>, plant height (length of the tallest tiller)<sup>(TODO 208)</sup>, spike length, number of tillers, grain weight<sup>(TODO 209)</sup> and number of grains per spike<sup>(TODO 210)</sup>. Flowering time, plant height and spike length were recorded in all three seasons, number of tillers in 2019–20 and 2020–21, grain weight in 2019–20 and 2021–22, and number of grains per spike in 2019–20 only. Four grain-composition traits (protein, starch, β-glucan and fiber) were measured by near-infrared spectroscopy, as described below.

---

## Change 4 of 5 · M&M · NIR ¶1

*Why:* One-sentence paragraph merged into ¶2.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Protein, starch, β-glucan and fiber contents were measured in percentages on the cleaned grain of each plant with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments), which scans the grain without destroying it.</del><span style="color:#1f4e9c;font-weight:bold;font-style:italic"> (→ opens the next paragraph, shortened) </span>

**Reads after approval:**

> *(paragraph removed; nothing is left of it in Word)*

---

## Change 5 of 5 · M&M · NIR ¶2

*Why:* Opening sentence from ¶1, shortened; aim first + 'we'; kit sentence split; calibration sentences plain.

**With the changes marked:**

<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). </ins><del style="color:#c0392b;background:#fdecea">The instrument was calibrated for wild barley on 30 accessions, one randomly chosen from each sampling site.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites.</ins> Their grain was scanned and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), <del style="color:#c0392b;background:#fdecea">adapted to small samples, as the mean of three subsamples of about 200 mg per accession</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">with the protocols adapted to small samples; each accession was assayed in three subsamples of about 200 mg, and their mean was used</ins>. Total protein was determined by the Kjeldahl method on 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College (Beljkaš et al. 2010)<sup>(TODO 212)</sup>. <del style="color:#c0392b;background:#fdecea">For protein and starch, the manufacturer's calibration package for cultivated barley was adjusted to the laboratory values, whereas the β-glucan calibration was developed from the kit values and added to the package. Fiber was estimated with the cultivated-barley calibration, as no laboratory reference was available for this trait.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.</ins><sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (Online Resource N<sup>(TODO 214)</sup>).<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> The protein, starch, β-glucan and fiber contents (%) of the cleaned grain of each plant were measured non-destructively with a DA 7250 near-infrared (NIR) analyzer (Perten Instruments). To calibrate the instrument for wild barley, we chose one accession at random from each of the 30 sampling sites. Their grain was scanned and then milled with an IKA A11 analytical mill and sieved through a 50-mesh screen for laboratory analysis. Starch and β-glucan were determined with the Megazyme Total Starch and Mixed-Linkage β-Glucan assay kits (Megazyme), with the protocols adapted to small samples; each accession was assayed in three subsamples of about 200 mg, and their mean was used. Total protein was determined by the Kjeldahl method on 200 mg of ground grain at the Analytical Biochemistry Laboratory of Tel-Hai Academic College (Beljkaš et al. 2010)<sup>(TODO 212)</sup>. The manufacturer's calibration for cultivated barley was then adjusted to the laboratory values for protein and starch, and a new β-glucan calibration was built from the kit values and added to it. Fiber had no laboratory reference and was predicted with the cultivated-barley calibration.<sup>(TODO 213)</sup> The calibration was validated by Pearson correlation of the NIR predictions with the laboratory values (Online Resource N<sup>(TODO 214)</sup>).

---

**Your decision:** approve · discard · or tell me what to change.
