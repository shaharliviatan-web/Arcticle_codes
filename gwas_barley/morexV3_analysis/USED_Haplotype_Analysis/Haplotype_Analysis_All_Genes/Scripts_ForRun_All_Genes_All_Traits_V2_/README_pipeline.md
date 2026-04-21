# Barley CrossHap All-Genes Pipeline

## Purpose

This tool is a new all-genes version of the wild-only barley CrossHap haplotype-analysis pipeline.

It is designed to:

- auto-discover all gene-level harmonized VCF targets for three traits
- run CrossHap for both MGmin 2 and MGmin 3 in one pipeline
- generate cached analysis objects, combined PDFs, and heatmap PDFs
- produce publication-oriented omnibus statistics tables
- preserve a full audit trail with both human-readable logging and timestamped structured event logs

Traits handled by this pipeline:

- betaglucan
- fiber
- starch

All outputs are written only inside this tool directory.

## Folder Structure

The tool is organized as:

- Scripts
- Scripts/R
- Runs/<run_id>
- Runs/<run_id>/Cache
- Runs/<run_id>/CombinedPDF
- Runs/<run_id>/Heatmaps
- Runs/<run_id>/Logs
- Runs/<run_id>/Stats

### Scripts

Main files:

- config.yaml
- targets.yaml
- run_pipeline.R
- run_pipeline_logged.sh
- README_pipeline.md

R helpers:

- R/utils.R
- R/run_crosshap_harmonized_wild.R
- R/plot_combined_pdf.R
- R/plot_heatmaps.R

### Runs/<run_id>

This is the run output directory for the current pipeline configuration.

The exact run folder name is controlled by:

- run_id in config.yaml

#### Cache

Stores cached HapObject results:

- Cache/trait/gene/MGmin_2/HapObject.rds
- Cache/trait/gene/MGmin_3/HapObject.rds

#### CombinedPDF

Stores one combined PDF per gene per MGmin.

Each combined PDF contains:

- CrossHap tree page
- raw omnibus Kruskal summary page across epsilon values for that MGmin
- corrected test summary page for that MGmin
- gene-level corrected summary page
- one CrossHap visualization page per available epsilon
- one violin plot page per available epsilon

The per-epsilon pages still show the epsilon-specific raw omnibus Kruskal p-value.

The corrected values are shown on the dedicated summary pages.

#### Heatmaps

Stores one heatmap PDF per gene per MGmin per epsilon.

Each heatmap PDF contains one page per haplotype group for that epsilon.

#### Logs

Stores:

- text run logs from the shell wrapper
- one timestamped event log per run

Event log filename pattern:

- event_log_<run_id>_<YYYYmmdd_HHMMSS>.csv

#### Stats

Stores the main output tables:

- per_test_stats.csv
- gene_summary.csv
- failures_notes.csv

## How Targets Are Generated

Targets are not manually curated in this tool.

At the start of each run, the pipeline scans the raw harmonized VCF directories for:

- betaglucan
- fiber
- starch

Only files ending in .vcf.gz are used as targets.

Files such as .csi indexes are ignored.

For each discovered VCF, the pipeline writes one target entry into targets.yaml with:

- trait
- gene_file
- title

At this stage the title is simply the gene filename itself. No functional annotation is guessed.

## Input Data Layout

Trait base naming convention:

- betaglucan_corrected_V3
- fiber_corrected_V3
- starch_corrected_V3

For each trait:

- raw VCF directory = trait_corrected_V3_targets_gene_220000bp_HARMONIZED
- imputed VCF directory = trait_corrected_V3_imputed_VCFs_for_LD_targets_gene_220000bp_HARMONIZED
- phenotype file = trait_corrected_V3.pheno

## What MGmin Means Here

MGmin is passed into CrossHap and represents the minimum marker-group size threshold used during haplotype construction.

This pipeline runs both:

- MGmin = 2
- MGmin = 3

The goal is to evaluate whether the same gene produces meaningful and stable haplotype grouping under both marker-group thresholds.

## What Epsilon Means Here

The epsilon vector is:

- 0.05
- 0.2
- 0.4
- 0.5
- 0.6
- 0.8
- 0.85

CrossHap is run once per gene and MGmin, but it produces haplotype objects across all listed epsilon values.

Different epsilon values can sometimes produce effectively identical downstream groupings. This is why duplicate collapsing is required before within-gene multiple-testing correction.

## How CrossHap Is Run

The analysis logic follows the earlier wild-only tool.

For each trait, gene, and MGmin:

1. phenotype is loaded from the trait-specific .pheno file
2. the raw harmonized VCF is loaded
3. the imputed harmonized VCF is loaded
4. variants are intersected between raw and imputed sets using CHROM:POS
5. PLINK LD is computed from the imputed VCF
6. CrossHap haplotyping is run on the raw VCF using:
   - the LD matrix
   - the phenotype data
   - the epsilon vector
   - MGmin
   - minHap
   - hetmiss_as
   - keep_outliers

Important design choice:

- raw VCF is used for haplotyping and genotype display
- imputed VCF is used for LD estimation

This matches the earlier tool logic.

## What Counts As A Valid Test

The primary discovery unit is gene within trait.

Traits are handled separately:

- betaglucan
- fiber
- starch

A valid omnibus test is defined at the level of one gene, one MGmin, and one epsilon.

A test is valid only if all of the following are true:

- the HapObject entry for that epsilon exists
- the grouping is non-empty after filtering
- haplotype 0 is excluded
- samples with missing phenotype are excluded
- more than one haplotype group remains
- a Kruskal-Wallis p-value can be computed

If any of those conditions fail:

- the test is excluded from multiple-testing correction
- no fake p-value is assigned
- the event is written to failures_notes.csv
- the event is also logged

If a gene has zero valid tests across all MGmin and epsilon values:

- it is excluded from gene_summary.csv
- it appears only in failures_notes.csv

## Duplicate Test Collapsing Within Gene

Different epsilon values, and sometimes different MGmin values, can produce effectively identical omnibus results.

To avoid over-counting the same effective test, duplicate tests are collapsed within each gene before within-gene correction.

Duplicate identity is defined using:

- raw Kruskal p-value
- number of groups
- sorted group sizes

MGmin is intentionally not part of the duplicate identity.

This means two tests are treated as duplicates if they have the same:

- raw omnibus p-value
- number of groups
- group-size signature

Example:

- 10,7,5
- 5,10,7

These are treated as the same group-size signature after sorting.

The retained representative test is chosen using an explicit rule:

1. sort valid tests by epsilon ascending
2. if epsilon ties, sort by MGmin descending
3. collapse duplicate signatures
4. keep the first row after that explicit ordering

This means the retained representative always uses:

- smallest epsilon
- and, if needed, larger MGmin as the secondary tie-breaker

The detailed audit file keeps one canonical row per unique signature and records duplicate-collapse notes when relevant.

## Within-Gene Holm Correction

After duplicate collapsing, each gene has a set of unique valid tests.

Holm correction is applied within the gene across those unique tests only.

Examples:

- 1 unique test -> Holm-adjusted p equals raw p
- 5 unique tests -> Holm correction across 5
- 9 unique tests -> Holm correction across 9

The representative result for a gene is the unique valid test with the smallest within-gene Holm-adjusted omnibus p-value.

The representative result keeps:

- MGmin
- epsilon
- raw p-value
- Holm-adjusted p-value
- number of individuals
- number of groups
- group-size signature

If more than one unique valid test shares the same smallest within-gene Holm-adjusted p-value, the tie is resolved deterministically using:

1. smallest raw p-value
2. if still tied, smallest epsilon
3. if still tied, larger MGmin

This tie resolution is logged explicitly in the event log as:

- AMBIGUOUS_BEST_TEST_RESOLVED

## Trait-Level Across-Gene Correction

After within-gene correction, each valid gene contributes one number:

- its best within-gene Holm-adjusted omnibus p-value

Then, separately within each trait, the pipeline applies:

- BH/FDR
- Bonferroni

There is no pooling across traits.

Each valid gene therefore ends with:

- best raw Kruskal p-value
- best Holm-adjusted p-value within gene
- trait-level FDR-adjusted p-value
- trait-level Bonferroni-adjusted p-value
- significant_fdr flag
- significant_bonferroni flag

The significance threshold is 0.05.

## Combined PDF Statistics

The combined PDFs now show both raw and corrected values.

They include:

- epsilon-specific raw omnibus Kruskal p-values on the per-epsilon pages
- raw omnibus Kruskal summary across epsilon values
- within-gene Holm-adjusted omnibus values where the epsilon and MGmin combination was retained in the final unique-test set
- final gene-level values:
  - selected best MGmin
  - selected best epsilon
  - best raw p
  - best within-gene Holm p
  - trait-level FDR p
  - trait-level Bonferroni p
  - significant_fdr
  - significant_bonferroni

The corrected values are added after the full statistics workflow is complete, by regenerating the combined PDFs using the final statistics tables.

The corrected test summary pages are matched MGmin-specifically, so each MGmin PDF shows only the corrected rows relevant to that MGmin.

## Violin Plot Statistics

The violin plots are visual and post hoc.

They are not the main discovery criterion.

The main discovery criterion is always the omnibus Kruskal workflow described above.

Violin plots use:

- reference-vs-others comparisons only
- reference group chosen as the haplotype group farthest from the global phenotype mean
- Wilcoxon pairwise tests
- Holm correction across the displayed pairwise comparisons

The plot displays Holm-adjusted significance symbols, not raw pairwise significance.

Interpretation:

- the violin annotations are useful for explaining group differences visually
- they should not be used as the primary gene-discovery statistic

## Heatmaps

Heatmaps keep the same behavior as the earlier tool unless a consistency requirement forced otherwise.

For each valid HapObject epsilon:

- genotype values are reconstructed from the raw VCF
- variants are matched back to the CrossHap variant IDs
- rows are clustered within haplotype group using Hamming distance that ignores missing values
- columns are annotated by marker group

There are no heterozygotes expected in this dataset, so the existing genotype binarization approach remains appropriate.

## Output Tables

## per_test_stats.csv

One row per unique valid test after duplicate collapsing.

This is the detailed audit table.

Typical columns include:

- trait
- gene_file
- gene_name
- title
- MGmin
- epsilon
- test_signature
- n_ind
- n_groups
- group_sizes
- kw_p_raw
- kw_p_holm_within_gene
- is_best_within_gene
- combined_pdf_path
- heatmap_pdf_path
- notes

Use this file when you want to inspect the full within-gene testing structure.

## gene_summary.csv

One row per gene with at least one valid omnibus test and a representative result.

This is the main file for gene-level discovery.

Typical columns include:

- trait
- gene_file
- gene_name
- title
- n_total_valid_tests_before_collapse
- n_unique_valid_tests
- best_MGmin
- best_epsilon
- best_n_ind
- best_n_groups
- best_group_sizes
- best_kw_p_raw
- best_kw_p_holm_within_gene
- trait_fdr_p
- trait_bonferroni_p
- significant_fdr
- significant_bonferroni

## How To Use gene_summary.csv To Fish Significant Genes

This is the main downstream filtering file.

A practical workflow is:

1. sort within each trait by trait_fdr_p
2. inspect significant_fdr
3. inspect significant_bonferroni for stricter hits
4. use best_MGmin and best_epsilon to find the corresponding plots
5. inspect the matching combined PDF and heatmap PDF for biological interpretation

Recommended interpretation order:

- start with significant_fdr
- then inspect Bonferroni hits
- then visually inspect the representative gene plots

## failures_notes.csv

This file contains all failures, invalid tests, and no-valid-result cases.

Examples include:

- crosshap_failed
- cache_read_failed
- no_valid_haplotype_grouping
- no_valid_kw_pvalue
- combined_pdf_failed
- combined_pdf_postprocess_failed
- heatmap_failed
- no_valid_tests_for_gene

This file is important for auditability.

A gene that never produced a valid omnibus test appears only here and does not appear in gene_summary.csv.

## Event Logs

Each run creates a new structured event log file.

Filename pattern:

- event_log_<run_id>_<YYYYmmdd_HHMMSS>.csv

This complements the human-readable text log.

Columns include:

- timestamp
- level
- trait
- gene_file
- MGmin
- epsilon
- stage
- event
- message

Typical events include:

- run start
- target discovered
- MGmin start
- valid KW test found
- duplicate collapsed
- Holm correction complete
- ambiguous_best_test_resolved
- trait-level correction complete
- PDF generated
- heatmap generated
- failure encountered
- run end

Use this file when you need programmatic auditability or want to trace failures systematically.

## How To Rerun Later

From the Scripts directory, run:

```bash
bash run_pipeline_logged.sh