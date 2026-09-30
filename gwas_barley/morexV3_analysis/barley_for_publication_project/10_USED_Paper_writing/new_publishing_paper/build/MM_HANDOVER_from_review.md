# Text moved out of Results_Discussion.md for Materials and methods

Written by `proposed/review_tool.py apply` when an approved review item moves text to M&M. Each entry is the
removed text word for word, with where it came from. The M&M session writes it in its own words, and
checks every number against the project files first.

## Item 2.1 · from Results Ch. 2 ¶1 ("To identify the genomic regions…"), sentences 2–3 (2.1a)

*→ M&M; the threshold value also stays in the Fig. 3 caption*

> Genome-wide significance was set by a Bonferroni correction at α = 0.10 over the 111,017 LD-pruned SNPs (−log~10~(*P*) = 6.0454).

## Item 2.1 · from Results Ch. 3, elite paragraph ("To ask which of these haplotypes…"; ¶8 in the review's numbering), sentences 1–2 (2.1c)

*→ M&M: the shared-site rule; the per-gene site counts stay in Online Resource 3, now cited in the sentence above*

> Because the two call sets were produced independently, only positions carrying a record with identical reference and alternate alleles in both were compared, leaving 43 of 54, 33 of 41 and 12 of 16 wild SNPs at the three genes (Online Resource 3 <sup>(TODO 306)</sup>).

## Item 3b · a Discussion sentence that M&M must back (added 2026-09-27)

*Not moved text: a requirement.* Discussion §4 ¶1 now says: "Likewise, one fixed grouping setting for all genes can leave the
associated SNPs of a sparse gene ungrouped, so a real association can be missed even in a gene that was tested. This was the case at
the GDSL esterase/lipase on 7H, whose associated SNPs were grouped only at a different setting."

> M&M must report both settings. At the Ch. 3 settings (MGmin 2, ε 0.6) the three significant SNPs of `7HG0729030` were not grouped,
> and the test was null (P = 0.16 fiber, 0.36 starch). The Ch. 4 result uses MGmin 3, ε 0.9, chosen by the GWAS-informed two-step
> procedure: the smallest MGmin that groups the signal SNPs, then the ε that maximises assignment, with no haplotype-test P used.
> Sources: `03_01_7H_branch_Starch_Fiber_shared_signal_explore/README.md` §2; `archive_v1_MGmin2_eps0.6/results/tables/crosshap_epsilon_sweep.tsv`.

## Item 5 · from Results Ch. 2 ¶5 ("To search for potential candidate genes…"), last sentence (5.4)

*→ M&M, as the reason the genes were annotated from their protein sequences*

> , and only three of the 55 genes carried a functional description in the MorexV3 annotation
