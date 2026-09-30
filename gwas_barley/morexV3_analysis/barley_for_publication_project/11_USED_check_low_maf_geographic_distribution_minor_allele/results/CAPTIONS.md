# Captions — geographic origin of the lead-SNP alleles and the haplotype groups

<!-- Ready-to-paste captions for the figures and tables of step 11, written 2026-09-30 in the manuscript's TAG
style (build/UPDATED_Results_Discussion_Conclusions.md): bold "Fig." / "Table" + number, bold lower-case panel
letters, no punctuation after the number and none at the end. Numbering is left open ("S1", "S2", "[N]"): the
user fixes the Online Resource numbers when the supplement is arranged. Figures: 174 mm wide, Liberation Sans
(Arial metric) 9 pt as Figs. 4-5, 600 dpi, RGB; PNG for the docx build and LZW TIFF for submission, pixel-identical.
Each caption has a hidden src note; the "Columns" lists under the tables are column definitions for the ESM
sheet header, not part of the caption. US spelling, as in the manuscript ("gray", "colored"). -->

---

## Figures

![](figures/Fig_carrier_origin.png){width="6.85in"}

**Fig. S1** Geographic origin of the accessions carrying the lead-SNP alleles and the haplotypes of the four candidate genes. **a** Genotype at the lead SNP of each of the 36 genome-wide significant loci, grouped by trait: minor-allele carrier (black), major-allele carrier (light tan) and missing call (gray). Each row is labeled with the chromosome and position (Mb) of the lead SNP, its minor-allele frequency, the direction of the minor allele's effect on the trait (+, increase; −, decrease) and, where one was carried forward, the candidate gene. **b** Haplotype groups of the four candidate genes: member of the group (black), member of another group of the same gene (light tan) and accession not assigned to a group (gray). Each row is labeled with the group, its size and its mean genotype BLUP (% of grain; fiber / starch for the GDSL esterase/lipase). Columns are the 290 *H. spontaneum* accessions, ordered by ecological region (HZ1, North–Coast; HZ2, North–Desert; HZ3, Coast–Desert transitional zones) and by sampling site
<!-- src: scripts/08_figure_carrier_origin.R (2026-09-30), results/figures/Fig_carrier_origin.{png,tif} (174 x 232 mm).
Data: intermediates/lead_genotype_long.tsv (script 02; lead SNPs and minor allele from 01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv, genotypes from 01_.../intermediates/morexV3_290.bed; no heterozygous calls) and intermediates/haplotype_group_long.tsv (script 03; GPAT6, GH17, PHT4;3 from the step-04 crosshap cache, run loci_LDspan_eps06_V4, MGmin 2, epsilon 0.6 = the groups of Fig. 4 and Table 2; GDSL from 03_01_.../results/tables/mgmin3_haplotype_assignment.tsv, MGmin 3, epsilon 0.9 = the groups of Fig. 5). Regions and sites: 00_THIN_.../outputs/subsection_1/tables/A4_site_boxplot_data.csv (the Fig. 2a assignment). Group means: the BLUPs of the GWAS (data/inputs/<trait>_corrected_V3.pheno). Region colors are those of Fig. 2a; every region is also named in the band (the palette alone is not color-vision safe). -->

![](figures/Fig_desert_enrichment.png){width="6.85in"}

**Fig. S2** Desert enrichment of the lead-SNP alleles and the haplotypes of the four candidate genes. **a** For each lead SNP, the percentage of minor-allele carriers (filled) and of major-allele carriers (open) that originated from desert sites. **b** For each haplotype group, the percentage of its members (filled) and of the members of the other groups of the same gene (open) that originated from desert sites. The vertical line marks the desert share of the panel (50 of 290 accessions). Right, Benjamini–Hochberg *q* of two site-permutation tests: Desert, difference in the desert share between the two sets (two-sided); regions, difference in their distribution across the six ecological regions (in **b**, one test per gene comparing all its groups, shown on the gene row). *q* was computed within trait in **a**, within gene for the Desert test in **b**, and across the four genes for the regions test in **b**; *q* ≤ 0.05 in bold. Row labels as in Fig. S1
<!-- src: scripts/09_figure_desert_enrichment.R (2026-09-30), results/figures/Fig_desert_enrichment.{png,tif} (174 x 207 mm); values from results/tables/T02_lead_enrichment_tests.tsv (desert_pct_minor, desert_pct_major, q_desert, q_regions) and T04_haplotype_group_enrichment_tests.tsv (desert_pct_group, desert_pct_other_groups, q_desert, q_regions_gene). Site permutation: accessions kept in their site, region labels shuffled among the 29 sites (North 4 sites, other regions 5), 10,000 permutations, seed 20260925 (scripts/00_config.R, perm_tests). M&M must state the test and that it is conservative for carriers concentrated in one site (one site cannot distinguish region from site). -->

---

## Tables

**Table S[N]** Regional origin of the minor- and major-allele carriers at the lead SNPs of the 36 genome-wide significant loci, and the enrichment tests. For each lead SNP: minor-allele frequency, direction of the minor allele's effect, numbers of minor- and major-allele carriers and of missing calls, minor-allele carriers per ecological region, the numbers of regions, sites and desert sites they came from, the site contributing most carriers, the percentage of minor- and major-allele carriers from desert sites, the odds ratio of desert origin (minor vs major carriers, accession level), and *P* and Benjamini–Hochberg *q* (within trait) of the two site-permutation tests
<!-- src: results/tables/T02_lead_enrichment_tests.tsv (script 02). The per-region detail, including the major-allele carriers and missing calls per region, is T01_lead_carriers_by_region.tsv; it can be a second sheet of the same Online Resource. -->

Columns (T02): `trait`; `locus_id`; `lead_SNP` (chromosome:position, MorexV3); `gene` (candidate gene carried forward, if any); `lead_neg_log10_p` (GWAS); `MAF` (minor-allele frequency among called accessions); `minor_effect` (raises / lowers the trait); `n_minor`, `n_major`, `n_missing`; `minor_North` … `minor_HZ3` (minor-allele carriers per region); `n_regions_minor`, `n_sites_minor`, `n_desert_sites_minor`; `top_site`, `top_site_n`; `desert_pct_minor`, `desert_pct_major` (% of carriers from desert sites); `desert_OR` (accession-level odds ratio, effect size only: accessions of one site are not independent); `desert_direction`; `p_desert`, `q_desert` (site permutation, two-sided; BH within trait); `p_regions`, `q_regions` (site permutation, chi-square statistic over the six regions; BH within trait)

**Table S[N]** Regional origin of the haplotype groups of the four candidate genes, and the enrichment tests. For each group: size, mean genotype BLUP, members per ecological region, the numbers of regions, sites and desert sites represented, the site contributing most members, the percentage of its members and of the other assigned accessions of the same gene from desert sites, the odds ratio of desert origin (accession level), and *P* and Benjamini–Hochberg *q* of the site-permutation tests: Desert (group vs the other groups, two-sided) and regions (group vs the other groups), both corrected within gene, and the gene-level regions test comparing all groups of a gene, corrected across the four genes. Accessions not assigned to a group are excluded from the tests
<!-- src: results/tables/T04_haplotype_group_enrichment_tests.tsv (script 03); per-region counts including the unassigned accessions: T03_haplotype_group_members_by_region.tsv. Groups: GPAT6, GH17, PHT4;3 = Fig. 4 / Table 2 (step 04, MGmin 2, epsilon 0.6); GDSL = Fig. 5 (03_01 branch, MGmin 3, epsilon 0.9; one grouping for fiber and starch). -->

Columns (T04): `gene`, `gene_id`, `trait`; `group`; `n`; `n_assigned_gene`, `n_groups_gene`; `mean_BLUP_fiber`, `mean_BLUP_starch` (% of grain); `North` … `HZ3` (members per region); `n_regions`, `n_sites`, `n_desert_sites`; `top_site`, `top_site_n`; `desert_pct_group`, `desert_pct_other_groups`; `desert_OR`; `desert_direction`; `p_desert`, `q_desert` (two-sided, BH within gene); `p_regions_group`, `q_regions_group` (group vs other groups, BH within gene); `p_regions_gene`, `q_regions_gene` (all groups of the gene x six regions, BH across the four genes)

**Table S[N]** Minor-allele carriers at the lead SNPs compared with a matched genome-wide background. For each lead SNP, 500 SNPs matched on minor-allele frequency (± 0.005) and on the number of called accessions (± 10), and more than 2 Mb from any lead SNP: the percentage of minor-allele carriers from desert sites at the lead and its matched distribution, the ecological region contributing most carriers, and the number of sites the carriers came from
<!-- src: results/tables/T05_lead_matched_background.tsv (script 04; 17,956 distinct background SNPs from 01_.../intermediates/morexV3_290_freq.frq and morexV3_290.bed). p_bg_desert = share of the matched SNPs with a desert share at least as high as the lead's; p_bg_fewer_sites = share with carriers in as few sites or fewer. -->

**Table S[N]** Number of lead SNPs whose minor-allele carriers came mainly from each ecological region, observed and expected from the matched genome-wide background, per trait. *P* by simulation of the sum of the per-lead probabilities (10⁵ draws); lead SNPs of one trait share carriers, so *P* is optimistic
<!-- src: results/tables/T06_top_region_vs_background.tsv (script 04). Rows "fiber+betaglucan" and "all_but_protein" pool traits. Background composition by MAF class (all background SNPs pooled): T07_background_top_region_by_MAF.tsv - a supporting table, possibly not for the ESM. -->

**Table S[N]** Phenotypic difference between minor- and major-allele carriers at each lead SNP, and the part of it predicted by the carriers' sites of origin. Observed difference in genotype BLUPs (trait standard deviations); difference expected if each carrier had the mean of the major-allele carriers of its own site (or region); percentage of carriers whose value exceeds that of their own site's major-allele carriers in the direction of the allelic effect; number of carriers with at least one major-allele carrier at their site; and the observed difference without the two most extreme accessions of the trait. Descriptive; no association was re-tested
<!-- src: results/tables/T08_lead_within_site.tsv (script 05). raw_sign_agrees_with_beta = FALSE only at starch_L04 (6H:525,776,080): the raw carrier difference is opposite to the GWAS beta (3 PCs + kinship). -->

**Table S[N]** Minor-allele carriers shared between loci of the same trait. For each pair of loci: whether they lie on the same chromosome and their distance, the number of minor-allele carriers of each among the accessions called at both, the number shared, the number expected by chance, and a hypergeometric *P*
<!-- src: results/tables/T09_lead_carrier_sharing.tsv (script 06). -->

**Table S[N]** Number of lead-SNP minor alleles carried by each accession, per trait, with its genotype BLUPs, trait ranks and genome-wide rare-allele load (share of 12,460 background SNPs with minor-allele frequency < 0.10 at which it carries the minor allele)
<!-- src: results/tables/T10_accession_minor_allele_load.tsv (script 06; rare-allele load from script 04). -->

**Table S[N]** Genetic relatedness (EMMAX aIBS kinship) between accessions from the same site, from different sites of the same ecological region and from different regions
<!-- src: results/tables/T11_kinship_by_region.tsv (script 07; 01_.../intermediates/morexV3_kinship.aIBS.kinf, the kinship of the GWAS). -->

**Table S[N]** Carriers of the minor and major alleles and missing calls at each lead SNP, per ecological region
<!-- src: results/tables/T01_lead_carriers_by_region.tsv (script 02); may be merged with the T02 Online Resource as a second sheet. -->

**Table S[N]** Members of each haplotype group of the four candidate genes, and accessions not assigned to a group, per ecological region
<!-- src: results/tables/T03_haplotype_group_members_by_region.tsv (script 03); may be merged with the T04 Online Resource as a second sheet. -->
