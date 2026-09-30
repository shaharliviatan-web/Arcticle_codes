# Review item MM4: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM4 of the M&M voice review** (item 4 in the chat): results and defences of method choices come out of the M&M. The 7H parts of item 4 (the *P* values 0.16 / 0.36) go into MM5, which rewrites the 7H section anyway. The same paragraphs also get the clarity fixes (B11, B16 from the chat) and the voice fixes.
>
> **LD decay (Change 1):** the text keeps only the method. The result (r² = 0.2 at about 188 kb) and two settings move to the caption of the LD-decay Online Resource; TODO 219 says so.
>
> **ε and MGmin (Change 2):** the ε comparison is now one clause, with its numbers kept in the Online Resource. The Results already give 30 of 55. ε gets a short gloss ("DBSCAN clustering radius"), which matches what crosshap does, as the new Haplotype ¶1 now says. "Without reference to the phenotypes" stays, because it is the guarantee the reader needs.
>
> **Elite cultivars (Changes 3–4):** the old sentence read as if all 136 lines passed the 85% call rate; in fact 91 did (checked in 07_…/results_chapter_numbers.txt). "Compared sites" was used in ¶1 but defined in ¶2, so the definition moves up. "Was not rerun … not assigned" is now stated positively. The 7H call-rate result ("only two of the cultivars met 85%") is removed; the per-gene call rates are already in the elite-line Online Resource (TODO 231).
>
> Length: 381 → 333 words (−48). The elite section stays about the same length (the fix there is clarity, not length).
>

> **❓ Questions for you — please answer these together with your decision**
>
> **Q1.** **LD decay:** keep the method sentence (as proposed), or delete the paragraph and its Online Resource completely? LD decay is not used anywhere in the Results or Discussion.
>
> **Q2.** **ε gloss:** OK to add "a DBSCAN clustering radius (ε)"? On 23.9 you decided not to explain ε. This is three words, and without them the reader cannot tell what 0.6 is.
>

---

## Change 1 of 4 · M&M · Genome-wide association ¶3 (LD decay)

*Why:* The result (188 kb) and two settings go to the Online Resource caption; the text keeps the method.

**With the changes marked:**

Genome-wide LD decay was estimated from the pairwise r² <del style="color:#c0392b;background:#fdecea">values computed in PLINK for</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">of</ins> all SNP pairs within 2.5 Mb<del style="color:#c0392b;background:#fdecea"> with r² ≥ 0.001, across the 290 accessions and the seven chromosomes. The mean r² was calculated</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, computed in PLINK, averaged</ins> in 1-kb distance bins and smoothed with a LOESS curve<del style="color:#c0392b;background:#fdecea"> (span 0.10)</del> weighted by the number of <del style="color:#c0392b;background:#fdecea">SNP </del>pairs per bin<del style="color:#c0392b;background:#fdecea">. The smoothed r² decreased to 0.2 at about 188 kb</del> <del style="color:#c0392b;background:#fdecea">(Online Resource N)<sup>(TODO 219)</sup></del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">(Online Resource N)<sup>(TODO 219)</sup></ins>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> Genome-wide LD decay was estimated from the pairwise r² of all SNP pairs within 2.5 Mb, computed in PLINK, averaged in 1-kb distance bins and smoothed with a LOESS curve weighted by the number of pairs per bin (Online Resource N)<sup>(TODO 219)</sup>.

---

## Change 2 of 4 · M&M · Haplotype analysis ¶2 (ε, MGmin)

*Why:* The ε comparison becomes one clause, and its numbers (Results numbers) are removed; ε gets a 3-word gloss; the sentence no longer starts with lowercase 'crosshap'; 'based on genotypes alone' → 'without reference to the phenotypes'.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">crosshap was run with the same settings for every gene: ε = 0.6</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">All genes were analyzed with the same crosshap settings: a DBSCAN clustering radius (ε) of 0.6</ins>, a minimum marker-group size (MGmin) of <del style="color:#c0392b;background:#fdecea">2</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">two SNPs</ins> and a minimum haplotype size (minHap) of nine accessions, with all other settings at their defaults. MGmin = 2<del style="color:#c0392b;background:#fdecea">, the smallest marker-group size,</del> allowed genes with few SNPs to be tested<del style="color:#c0392b;background:#fdecea">. In a comparison of fixed ε values from 0.2 to 1.0 based on genotypes alone, larger values made more genes testable but assigned fewer accessions to haplotype groups; ε = 0.6 balanced the two, with 30 of the 55 genes testable and a median of 65% of the accessions assigned</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">, and ε = 0.6 was chosen without reference to the phenotypes, as the value that balanced the number of testable genes against the proportion of accessions assigned to haplotype groups</ins> (Online Resource N)<sup>(TODO 223)</sup>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> All genes were analyzed with the same crosshap settings: a DBSCAN clustering radius (ε) of 0.6, a minimum marker-group size (MGmin) of two SNPs and a minimum haplotype size (minHap) of nine accessions, with all other settings at their defaults. MGmin = 2 allowed genes with few SNPs to be tested, and ε = 0.6 was chosen without reference to the phenotypes, as the value that balanced the number of testable genes against the proportion of accessions assigned to haplotype groups (Online Resource N)<sup>(TODO 223)</sup>.

---

## Change 3 of 4 · M&M · Comparison with elite cultivars ¶1

*Why:* Aim first, 'we'; 'compared sites' defined before use; the 136-line sentence read as if all 136 passed; names no longer set off by commas mid-sentence.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Genotypes</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To compare the wild haplotypes with modern cultivars, we retrieved the genotypes</ins> of elite cultivars in the same gene windows (gene ± 1 kb) <del style="color:#c0392b;background:#fdecea">were retrieved </del>from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024)<sup>(TODO 230)</sup> in DivBrowse (König et al. 2023). <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource N)<sup>(TODO 232)</sup>. </ins><del style="color:#c0392b;background:#fdecea">Five widely grown spring malting cultivars released between 2012 and 2018, Avalon, KWS Irina, Odyssey, Laureate and LG Diablo, were chosen among the 136 spring elite lines of the panel that were called at ≥ 85% of the compared sites, pooled over the three genes; their names were verified in their EBI BioSamples records</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records</ins> (Online Resource N)<sup>(TODO 231)</sup>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> To compare the wild haplotypes with modern cultivars, we retrieved the genotypes of elite cultivars in the same gene windows (gene ± 1 kb) from the unimputed MorexV3 SNP calls of the barley pangenome panel (Jayakodi et al. 2024)<sup>(TODO 230)</sup> in DivBrowse (König et al. 2023). Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource N)<sup>(TODO 232)</sup>. Of the 136 spring elite lines of the panel, 91 were called at ≥ 85% of the compared sites of the three genes together, and among them we chose five widely grown spring malting cultivars released between 2012 and 2018: Avalon, KWS Irina, Odyssey, Laureate and LG Diablo. Their names were verified in their EBI BioSamples records (Online Resource N)<sup>(TODO 231)</sup>.

---

## Change 4 of 4 · M&M · Comparison with elite cultivars ¶2

*Why:* Shared-site sentence moved to ¶1; 'not rerun' stated positively; the 7H call-rate result removed (→ Online Resource).

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">Only sites with the same position and the same reference and alternate alleles in both call sets were compared (Online Resource N)<sup>(TODO 232)</sup>. </del><span style="color:#1f4e9c;font-weight:bold;font-style:italic"> (→ ¶1) </span>For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. 4d–f). <del style="color:#c0392b;background:#fdecea">The haplotype analysis was not rerun with the cultivars, and the cultivars were not assigned to haplotype groups.</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups.</ins> For the 7H gene (Fig. 5c), the same five cultivars and procedure were used without a new call-rate screen<del style="color:#c0392b;background:#fdecea">, and only two of the cultivars met the 85% call-rate threshold at this gene</del>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> For each wild haplotype group, a consensus genotype was taken at each site as the majority allele among the accessions of the group, ignoring missing calls, with ties resolved to the reference allele (Fig. 4d–f). The cultivars were compared with these consensus genotypes only; they were not included in the haplotype analysis or assigned to haplotype groups. For the 7H gene (Fig. 5c), the same five cultivars and procedure were used without a new call-rate screen.

---

**Your decision:** approve · discard · or tell me what to change.
