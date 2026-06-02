#!/usr/bin/env Rscript

# -----------------------------------------------------------------
# SCRIPT START: imputation_core.R
# -----------------------------------------------------------------
# This script performs the core imputation logic using the LEA package.
# It expects 'morexV3_samples_removed.vcf' to exist in the same directory.
# -----------------------------------------------------------------

# --- Step 1: Load Libraries ---
# (We assume packages are already installed)
library(LEA)
library(vcfR)
library(data.table)

print("[R] Libraries loaded.")

# --- Step 2: Data Conversion ---
print("[R] Starting VCF to LFMM conversion...")
vcf2lfmm("morexV3_samples_removed.vcf", "morexV3_samples_removed.lfmm")
print("[R] LFMM conversion finished.")

# --- Step 3: Run sNMF for K=3 ---
print("[R] Starting sNMF (K=3, 10 reps)... This will take a long time.")
project_k3 = snmf("morexV3_samples_removed.lfmm",
                  K = 3,
                  entropy = TRUE,
                  repetitions = 10,
                  project = "new")
print("[R] sNMF run finished.")

# --- Step 4: Select the Best Run ---
best_run = which.min(cross.entropy(project_k3, K = 3))
print(paste("[R] The best run selected for K=3 is:", best_run))

# --- Step 5: Perform Imputation ---
print("[R] Starting the imputation process... This will also take a long time.")
impute(project_k3,
       "morexV3_samples_removed.lfmm",
       method = 'mode',
       K = 3,
       run = best_run)

print("[R] Imputation process finished!")
print("[R] Output file: morexV3_samples_removed.lfmm_imputed.lfmm")
print("[R] R script complete.")