# Review item MM9b: proposed changes

Proposal only. `Methods_Results_Discussion_Conclusions.md` is unchanged until you approve.
<del style="color:#c0392b;background:#fdecea">Red, struck through</del> = to be deleted · <ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">green</ins> = to be inserted · <span style="color:#1f4e9c;font-weight:bold;font-style:italic">(blue)</span> = where removed text goes; preview only, never in the manuscript. Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders marked).

> **ℹ️ Notes**
>
> **Item MM9b** = change 3 of MM9 (Loci ¶1), corrected after your check. You were right that a SNP belongs to one locus only. Checked in the outputs: no SNP is a member of two loci, and no two loci of the same trait overlap (`02_loci_FINAL/tables/loci_members.tsv`, `loci_summary.tsv`).
>
> **What the later rounds did:** only significant SNPs left outside every locus could lead. My clause "allowed to join but not to lead" described a code detail: already-assigned SNPs are kept as possible members, but in this run that never produced an overlap. That detail is dropped, and the rounds are described by what they did. "Non-overlapping" is added, because it is now verified and it answers the question a reader would ask.
>
> Rounds, for the record (hidden note): 31 loci in round 1, 4 in round 2, and 1 in round 3 (a lead-only locus).
>
> The rest of the paragraph is as in MM9 (aim first + "we", the gap rule in plain words).
>

## Change 1 of 1 · M&M · Loci and candidate genes ¶1

*Why:* Clumping and gap rule in plain words; aim first + 'we'; the repeat rounds described by what they did (leftover significant SNPs lead); loci verified non-overlapping.

**With the changes marked:**

<del style="color:#c0392b;background:#fdecea">For each trait, the significant SNPs were grouped into loci</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">To define loci, we grouped the significant SNPs of each trait</ins> by LD clumping in PLINK (Purcell et al. 2007)<del style="color:#c0392b;background:#fdecea">, using the genotypes of the 290 accessions</del>. Only <del style="color:#c0392b;background:#fdecea">SNPs above the genome-wide threshold could lead a locus</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">significant SNPs could serve as lead SNPs</ins>, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then <del style="color:#c0392b;background:#fdecea">limited to the contiguous run of member SNPs around its lead SNP, ending at the first gap longer than 50 kb on either side</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs</ins>. <del style="color:#c0392b;background:#fdecea">Clumping was repeated, with SNPs already assigned to a locus excluded as leads, until every significant SNP lay within a locus, giving 36 loci</del><ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">Significant SNPs left outside every locus then served as lead SNPs in further rounds of clumping, until every significant SNP lay within a locus, giving 36 non-overlapping loci</ins> (Online Resource N)<sup>(TODO 220)</sup>.<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none"> </ins>

**Reads after approval:**

> To define loci, we grouped the significant SNPs of each trait by LD clumping in PLINK (Purcell et al. 2007). Only significant SNPs could serve as lead SNPs, and any SNP within 2 Mb of a lead SNP and in LD with it (r² ≥ 0.5) joined its locus, regardless of its *P* value. Each locus was then trimmed to a continuous block around its lead SNP: moving outward on each side, it ended at the first gap of more than 50 kb between consecutive member SNPs. Significant SNPs left outside every locus then served as lead SNPs in further rounds of clumping, until every significant SNP lay within a locus, giving 36 non-overlapping loci (Online Resource N)<sup>(TODO 220)</sup>.

---

**Your decision:** approve · discard · or tell me what to change.
