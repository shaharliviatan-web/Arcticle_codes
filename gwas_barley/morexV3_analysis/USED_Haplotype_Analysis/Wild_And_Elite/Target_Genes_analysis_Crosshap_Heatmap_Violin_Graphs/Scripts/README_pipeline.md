# Wild + Elite CrossHap Target-Gene Pipeline

## Purpose

This tool is the simplified wild+elite barley haplotype-analysis pipeline for the fixed 8 target genes.

It is designed to:

- use the merged intersection wild+elite VCFs for haplotyping
- use the harmonized imputed wild VCFs for LD estimation
- run one fixed epsilon and one fixed MGmin per gene
- generate CrossHap pages, violin plots, heatmap PDFs, and elite haplotype assignment tables

This version does not do an all-genes scan, epsilon sweeps, or across-gene correction.

## Inputs

- merged intersection VCFs from Wild_And_Elite_VCFs__IntersectionMethod
- harmonized imputed wild VCFs from the trait-specific LD folders
- wild+elite phenotype files from phenotype_files_wild&elite

Phenotype files already contain only the 6 selected elite lines:

- Corvette
- Triumph
- Iris
- Gunnar
- Manchuria_A
- Rapid

## Folder Structure

- Scripts
- Scripts/R
- Runs/<run_id>
- Runs/<run_id>/Cache
- Runs/<run_id>/CombinedPDF
- Runs/<run_id>/Heatmaps
- Runs/<run_id>/Logs
- Runs/<run_id>/Stats

## Main Files

- config.yaml
- targets.yaml
- run_pipeline.R
- run_pipeline_logged.sh
- README_pipeline.md

R helpers:

- R/utils.R
- R/run_crosshap_intersection_wild_elite.R
- R/plot_combined_pdf.R
- R/plot_heatmaps.R

## Methods

For each target gene:

1. load the merged intersection wild+elite VCF
2. load the matching harmonized imputed wild VCF
3. intersect variants between merged and LD VCFs by CHROM:POS
4. compute LD with PLINK from the imputed wild VCF
5. run CrossHap on the merged intersection VCF using the agreed epsilon and MGmin
6. compute phenotype statistics on wild lines only after excluding hap 0
7. draw violin plots with elite lines shown as labeled points below the plot
8. apply Holm correction to the pairwise violin significance symbols
9. export elite haplotype assignments to CSV
10. draw haplotype-group heatmaps from the merged intersection VCF

Important design choice:

- merged intersection VCF is used for haplotyping and genotype display
- imputed harmonized wild VCF is used for LD estimation

## Running

From the Scripts directory:

```bash
bash run_pipeline_logged.sh
```

Or directly:

```bash
Rscript run_pipeline.R config.yaml targets.yaml
```

## Main Outputs

For each run:

- one cached HapObject per gene
- one combined PDF per gene
- one heatmap PDF per gene
- one elite assignment CSV per gene
- one summary CSV with Kruskal statistics and output paths
- one failures CSV and one event log
