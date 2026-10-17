# Unresolved points kept outside the manuscript

Open questions that were Word comments (TODOs) in the manuscript and were **removed from it at the user's request**, so that
the copy sent for review shows fewer comments. **They are not resolved**: each still needs an answer, and the manuscript text
they concern may change when it comes. When one is answered, fix the text through the usual review item, then mark it resolved
here (date, what changed). Their comment IDs are retired: never reuse them.

## Removed 2026-10-01 (item CMT1)

| ID | section | the text the comment was attached to | status |
|---|---|---|---|
| 201 | Plant material and genotyping | "three rounds of single-seed descent and selfing" | unresolved |
| 202 | Plant material and genotyping | "as described by Potapenko et al. (2026b)" | unresolved |
| 205 | Plant material and genotyping | "more than 30% missing data were removed" | unresolved |
| 206 | Common-garden experiment and phenotyping | "as described by Potapenko et al. (2026b)" | unresolved |
| 207 | Common-garden experiment and phenotyping | "flowering time (days from sowing to heading)" | unresolved |
| 208 | Common-garden experiment and phenotyping | "plant height (length of the tallest tiller)" | unresolved |
| 209 | Common-garden experiment and phenotyping | "grain weight" | unresolved |
| 210 | Common-garden experiment and phenotyping | "number of grains per spike" | unresolved |
| 211 | Common-garden experiment and phenotyping | "threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22" | unresolved |

### TODO 201 — Plant material and genotyping

**Sentence in the manuscript** (the commented words in bold):

> At each of the 30 sites, spikes were collected from ten plants at least 10 m apart, and each accession was derived from one plant through **three rounds of single-seed descent and selfing**.

**Comment, as it was in the manuscript:**

> TODO - number of rounds of single-seed descent. The text says three rounds, as in the paper cited for the plant material, Potapenko et al. 2026b (Mol Ecol): "This procedure (single-seed-descent and selfing) was repeated three times". The genome-size paper on the same collection, Potapenko et al. 2026a (Mol Biol Evol), says two rounds, and so do our earlier Materials and methods drafts (May and June 2026). N. Pintel's thesis does not say. Please confirm the number of rounds.

**Status:** unresolved.

### TODO 202 — Plant material and genotyping

**Sentence in the manuscript** (the commented words in bold):

> The raw reads (European Nucleotide Archive, PRJEB79623) were quality-checked and trimmed **as described by Potapenko et al. (2026b)**.

**Comment, as it was in the manuscript:**

> TODO - how the reads were cleaned before the MorexV3 re-call. The text says the raw reads were quality-checked and trimmed as described by Potapenko et al. 2026b (Mol Ecol), that is, with FastQC v0.11.5 and fastp v0.20.0 as in that paper. This is right only if the MorexV3 re-call started from the reads already cleaned for that paper. If the reads were cleaned again for this study, please give the tools, versions and settings, and they will replace the citation.

**Status:** unresolved.

### TODO 205 — Plant material and genotyping

**Sentence in the manuscript** (the commented words in bold):

> Heterozygous genotypes and genotypes supported by fewer than three reads were set to missing, and SNPs with a minor allele frequency (MAF) below 0.05 or **more than 30% missing data were removed**.

**Comment, as it was in the manuscript:**

> TODO - the missing-data filter of the SNPs. The text says that SNPs with more than 30% missing data were removed, but the final SNP set contains SNPs with up to 31.7% missing data (40 SNPs above 30% over the 300 accessions; 155 over the 290 analysed). The likely reason is the order of the filters: the 30% filter was applied before genotypes supported by fewer than three reads were set to missing, which added missing calls afterwards. Please confirm the order and the threshold, so the sentence can state them correctly. The other filter, minor-allele frequency 0.05, holds exactly.

**Status:** unresolved.

### TODO 206 — Common-garden experiment and phenotyping

**Sentence in the manuscript** (the commented words in bold):

> The accessions were grown in a common-garden experiment over three growing seasons (2019–20, 2020–21 and 2021–22) in an insect-proof net-house in the Hula Valley, Israel (33°09′08.4″ N 35°37′15.8″ E), in 5-L pots under full irrigation, **as described by Potapenko et al. (2026b)**.

**Comment, as it was in the manuscript:**

> TODO - the pre-sowing treatment of the seeds. The text gives no details and cites Potapenko et al. 2026b (Mol Ecol), whose protocol is: seeds dried at 30 °C for 14 days, dipped in 4% sodium hypochlorite, sown in mid-October on germination paper, kept 10 days in the dark at 4 °C, then moved to 5-L pots with a commercial soil mix. The sources disagree: Potapenko et al. 2026a (Mol Biol Evol) gives a 3% chlorine solution, sowing trays, 10 days in the dark at 4 °C and 14 days in a growth chamber at 20 °C before transplanting; N. Pintel's thesis and our old drafts give 3% sodium hypochlorite for 5 minutes, a 5-minute rinse in distilled water and germination paper at 16 °C for 18 days. Please confirm which protocol was used in this experiment; if it differs from the cited paper, the details will be written here.

**Status:** unresolved.

### TODO 207 — Common-garden experiment and phenotyping

**Sentence in the manuscript** (the commented words in bold):

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): **flowering time (days from sowing to heading)**, plant height (length of the tallest tiller), spike length, number of tillers, grain weight and number of grains per spike.

**Comment, as it was in the manuscript:**

> TODO - the definition of flowering time. The text follows the cited paper, Potapenko et al. 2026b (Mol Ecol): days from sowing to heading. Other sources define it differently: Potapenko et al. 2026a (Mol Biol Evol) counts days from transplanting to the emergence of the awns (Zadoks stage 49), and N. Pintel's thesis and the old draft count days from transplanting to the emergence of the first spike from the flag-leaf sheath. Our data fit the cited definition: in 2020-21 they range from 73 to 147 days, as reported in Mol Ecol. Please confirm.

**Status:** unresolved.

### TODO 208 — Common-garden experiment and phenotyping

**Sentence in the manuscript** (the commented words in bold):

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading), **plant height (length of the tallest tiller)**, spike length, number of tillers, grain weight and number of grains per spike.

**Comment, as it was in the manuscript:**

> TODO - the definition of plant height. The text follows the cited paper, Potapenko et al. 2026b (Mol Ecol): the length of the tallest tiller, which is the measurement used in our analyses. N. Pintel's thesis defines plant height instead as the height from the ground to the base of the highest flag leaf, which Potapenko et al. 2026a (Mol Biol Evol) records as a separate trait (flag-leaf height). Please confirm that the length of the tallest tiller is right.

**Status:** unresolved.

### TODO 209 — Common-garden experiment and phenotyping

**Sentence in the manuscript** (the commented words in bold):

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading), plant height (length of the tallest tiller), spike length, number of tillers, **grain weight** and number of grains per spike.

**Comment, as it was in the manuscript:**

> TODO - the definition of grain weight. Neither cited paper defines it, so the text gives no definition. Our old draft (May 2026) says the cleaned grain of the five harvested spikes was weighed per plant, but in the data the grains were counted on three spikes per plant (two for 13 plants in 2019-20), so it is not clear how many spikes were weighed. Please give the definition: what was weighed, on how many spikes, and in which unit.

**Status:** unresolved.

### TODO 210 — Common-garden experiment and phenotyping

**Sentence in the manuscript** (the commented words in bold):

> Six morphological and phenological traits were recorded as described by Potapenko et al. (2026b): flowering time (days from sowing to heading), plant height (length of the tallest tiller), spike length, number of tillers, grain weight and **number of grains per spike**.

**Comment, as it was in the manuscript:**

> TODO - the definition of grain number. The text follows the cited paper, Potapenko et al. 2026b (Mol Ecol): number of grains per spike. The values used in our analyses, however, are the total number of grains on the spikes counted per plant (three spikes; two for 13 plants in 2019-20); grains per spike is a separate column that we did not use, and the Mol Ecol figures themselves label the trait "Number of grains". Please confirm which it is. If it is the total, the text should say "number of grains on three spikes per plant".

**Status:** unresolved.

### TODO 211 — Common-garden experiment and phenotyping

**Sentence in the manuscript** (the commented words in bold):

> The dried spikes were **threshed by hand in 2019–20 and with a Haldrup LT-21 thresher in 2020–21 and 2021–22**, and the grain was cleaned with a South Dakota seed blower (Seedburo Equipment Company) and stored for nutritional analysis.

**Comment, as it was in the manuscript:**

> TODO - how the grain was threshed. The text says by hand in 2019-20 and with a Haldrup LT-21 thresher in 2020-21 and 2021-22. N. Pintel's thesis, which covers only the first two seasons, says by hand in the first season and with the Haldrup LT-21 in the second. For 2021-22 the only source is our old draft (May 2026). Please confirm the threshing in 2021-22.

**Status:** unresolved.

## Shortened comment, for the record

TODO 401 (Discussion, "A shared carbon-allocation axis for starch, fiber and β-glucan") stays in the manuscript, shortened to the
user's wording on 2026-10-01. The part removed from it:

> Our experiment grew two cultivated checks, Morex and Clipper, but they are not in the analysed data. N. Pintel's NIR files contain one scan of each, from 17 November 2020: Morex 7.0% protein, 64.7% starch, 4.1% beta-glucan and 2.7% fiber; Clipper 12.9%, 59.7%, 5.4% and 4.2%; the wild accessions averaged about 16% protein and 33-43% starch in the first two seasons. One scan per cultivar is too little to replace the citation, but it points the same way. Decide whether to keep the citation alone, add these values, or find a second source.

