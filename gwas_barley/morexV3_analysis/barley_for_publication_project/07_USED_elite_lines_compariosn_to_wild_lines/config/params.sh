# =============================================================================
# 07_elite_lines_comparison_to_wild_lines -- EVERY parameter, single source of truth
# =============================================================================
# Bash sources this file directly; R reads the same file via scripts/_load_params.R,
# so the two can never drift. Nothing in scripts/ may hard-code a value that lives
# here. Convention copied from step 03 (config/params.sh).
#
# Created 2026-09-10.
# =============================================================================

# ---- Roots ------------------------------------------------------------------
PROJECT_ROOT="/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project"
STEP_ROOT="${PROJECT_ROOT}/07_elite_lines_compariosn_to_wild_lines"

DIR_CONFIG="${STEP_ROOT}/config"
DIR_SCRIPTS="${STEP_ROOT}/scripts"
DIR_INPUTS="${STEP_ROOT}/inputs"
DIR_INTERMEDIATES="${STEP_ROOT}/intermediates"
DIR_RESULTS="${STEP_ROOT}/results"
DIR_TABLES="${DIR_RESULTS}/tables"
DIR_FIGURES="${DIR_RESULTS}/figures"
DIR_LOGS="${STEP_ROOT}/logs"

DIR_ELITE_RAW="${DIR_INTERMEDIATES}/elite_vcfs_raw"
DIR_ELITE_TRIMMED="${DIR_INTERMEDIATES}/elite_vcfs_trimmed"
DIR_MATRICES="${DIR_INTERMEDIATES}/matrices"

# ---- Upstream inputs (step 04 is the ONLY source of wild haplotypes) --------
STEP04_ROOT="${PROJECT_ROOT}/04_USED_haplotype_analysis_crosshap"
STEP04_RUN_ID="loci_LDspan_eps06_V4"
STEP04_RUN_ROOT="${STEP04_ROOT}/04_runs/${STEP04_RUN_ID}"
STEP04_CACHE="${STEP04_RUN_ROOT}/Cache"
STEP04_GENE_WINDOWS="${STEP04_ROOT}/00_config/gene_windows.tsv"
STEP04_GENE_RESULTS="${STEP04_RUN_ROOT}/Stats/gene_results.tsv"
STEP04_RAW_VCF_DIR="${STEP04_ROOT}/03_per_gene_vcfs/raw_1000bp"

# crosshap parameters of that run -- used only to locate the cached HapObject.
# They are NOT re-run here: the wild haplotype groups are consumed as given.
CROSSHAP_MGMIN=2
CROSSHAP_EPSILON="0.6"

# Step-05 annotation, for the gene's functional name in figure titles/tables.
STEP05_PAPER_TABLE="${PROJECT_ROOT}/05_USED_gene_annotation_analysis/08_USED_annotation_master/results/tables/Table_significant_genes_paper.tsv"

# ---- The three genes carried forward by step 06 -----------------------------
# "gene_id:short_name" -- short_name is used in filenames and figure titles.
# These are the three tier-1 candidates of 06_USED_genes_selected_to_present.
TARGET_GENES=(
  "HORVU.MOREX.r3.3HG0301300:GPAT6"
  "HORVU.MOREX.r3.5HG0487060:GH17"
  "HORVU.MOREX.r3.3HG0301710:PHT4-3"
)

# ---- Window -----------------------------------------------------------------
# Must match step 04's config.yaml window_bp, or the wild and elite windows would
# not describe the same sequence. Asserted at run time by 01_fetch_elite_vcfs.sh.
WINDOW_BP=1000

# ---- DivBrowse (IPK barley pangenome v2) ------------------------------------
DIVBROWSE_BASE="https://divbrowse.ipk-gatersleben.de/barley_pangenome_v2"
DIVBROWSE_INDEX="${DIVBROWSE_BASE}/"

# The /vcf_export endpoint returns the CORRECT NUMBER of variants but reads them
# from a SHIFTED position range (verified 2026-09-10: a request for
# 546,472,709-546,478,564 returned 158 records spanning 546,476,514-546,481,290;
# the drift is not constant). The endpoint is therefore always called with a
# padded window and the result is trimmed locally, and 01_fetch_elite_vcfs.sh
# ASSERTS that the returned span brackets the target window on both sides.
# Without that assertion the failure is silent data loss.
DIVBROWSE_PAD_BP=50000
DIVBROWSE_PAD_MAX_BP=400000     # widen-and-retry ceiling
DIVBROWSE_RETRY_FACTOR=4        # pad multiplier per retry
DIVBROWSE_TIMEOUT=900        # a padded window on 5 samples takes ~3 min server-side

# EBI BioSamples -- the authoritative source for what each SAMEA accession IS.
BIOSAMPLES_BASE="https://www.ebi.ac.uk/biosamples/samples"
BIOSAMPLES_PARALLEL=12

# ---- Elite line selection ---------------------------------------------------
# The five lines themselves live in config/elite_lines.tsv (name + SAMEA + why).
ELITE_LINES_TSV="${DIR_CONFIG}/elite_lines.tsv"

# Candidate POOL the lines are chosen from, and the call-rate screen applied to it.
# Added 2026-09-14. The first selection (2026-09-10) turned out to include two lines
# with many no-calls at our sites (RGT Planet 53%, Propino 63%), which show up as
# grey tiles. The DivBrowse panel is unimputed low-coverage WGS, so missingness is
# a property of each sample's sequencing depth and can be screened for.
#
# The pool is downloaded ONCE (01_fetch_elite_vcfs.sh); 02_screen_elite_lines.R
# scores every pool line; the configured lines are then SUBSET from that same
# download, so the screen and the figures are guaranteed to use identical data.
#
# CALL RATE = called genotypes (REF, ALT or HET) / shared sites, where shared sites
# are the wild SNPs that also carry an elite record with identical REF and ALT --
# i.e. exactly the columns the elite rows can be compared on. It is judged pooled
# over the three genes; per-gene rates are reported alongside.
#
# Selecting on call rate is a genotype-QUALITY criterion only. It never looks at
# which allele or haplotype a line carries, so it cannot bias the comparison.
ELITE_POOL_ACCESSION_TYPE="elite lines"
ELITE_POOL_ANNUALITY="spring type"   # spring-only, see config/elite_lines.tsv
CALL_RATE_MIN="0.85"                 # every configured line must pass; asserted

# ---- The three SNP-matching versions ----------------------------------------
# The two VCFs do not carry the same SNPs, in BOTH directions. All three
# treatments are built; the choice between them is deliberately left open.
#
#   shared_sites  (a) plot only SNPs present in wild AND elite (CHROM:POS:REF:ALT).
#                     No assumption. Costs columns.
#   filled_marked (b) plot ALL wild SNPs; elite tiles with no elite record get
#                     their own colour "no record". Honest, keeps the haplotype whole.
#   filled_silent (c) the earlier behaviour: bcftools merge -0, missing elite
#                     records silently become REF. Kept for comparison only.
VERSIONS=("shared_sites" "filled_marked" "filled_silent")

# ---- Allele concordance between the two VCFs (HARD CHECK) -------------------
# Both call sets are called against Morex V3 and the wild set was produced without
# any PLINK allele-flipping, so at a shared position REF and ALT must be IDENTICAL,
# never swapped. If they were swapped, a 0/0 in one file and a 0/0 in the other
# would mean OPPOSITE alleles and every barcode comparison in this step would be
# silently inverted. 02_build_matrices.R classifies every shared site as
# identical / swapped / different and writes results/tables/allele_concordance.tsv.
#
#   stop  -- abort the run if any shared site is swapped or carries different
#            alleles (the default; a swap is a data problem, not something to
#            paper over by flipping genotypes automatically)
#   warn  -- log and continue, marking the affected sites in the output tables
ALLELE_MISMATCH_ACTION="stop"

# ---- Genotype encoding / colours -------------------------------------------
# REF/ALT/missing match 04_.../01_scripts/R/plot_heatmaps.R exactly, so the new
# figures read the same as the per-gene heatmaps already in the paper.
# HET is NEW here: the wild VCFs contain no heterozygous call at all (verified
# 2026-09-10 -- only ./. 0/0 1/1 1|1), but the elite panel does, and step 04's
# gt_to_bin01() would have silently folded het into "missing".
COL_REF="#FFFACD"        # light yellow  -- reference allele (Morex V3)
COL_ALT="#2F4F4F"        # dark teal     -- alternate allele
COL_MISS="grey70"        # grey          -- ./. no call
COL_HET="#C46210"        # orange        -- 0/1, elite only
COL_NORECORD="#FFFFFF"   # white         -- version (b) only: no elite record at this site
COL_TRIALLELIC="#7B3FA0" # purple        -- elite line carries a THIRD allele (see below)
#
# Triallelic sites (added 2026-09-14). At two wild SNPs the elite file has a record
# with the SAME reference base but a DIFFERENT alternate allele (GPAT6
# 3H:545,003,306 wild G/A vs elite G/C; GH17 5H:462,728,384 wild G/A vs elite G/T).
# Those sites are not "shared" (key CHROM:POS:REF:ALT), so version shared_sites does
# not draw them. Versions filled_marked and filled_silent DO draw every wild SNP, and
# there the elite line's real genotype is used, never a fill:
#   0/0                 -> Reference   (same base as the wild REF, so comparable)
#   1/1 or 0/1 of the   -> Triallelic  (the line carries the elite-only third allele,
#   elite-only ALT                      which is neither the wild REF nor the wild ALT)
#   ./.                 -> Missing
# Legends list only the states actually drawn in each figure, so the triallelic
# colour appears only if some line really carries the third allele.
#
# Row gap: barcode tiles are drawn at this fraction of the row height, so each
# row reads as its own barcode rather than all rows fusing into one heatmap.
TILE_HEIGHT="0.72"

# ---- Figure geometry --------------------------------------------------------
FIG_WIDTH_IN=11
FIG_HEIGHT_IN=8.5
FIG_DPI=400

# ---- Tools ------------------------------------------------------------------
# bare 'plink'/'bcftools' can resolve to the wrong build on this host; absolute.
BCFTOOLS_BIN="/usr/bin/bcftools"
BGZIP_BIN="/usr/bin/bgzip"
TABIX_BIN="/usr/bin/tabix"
RSCRIPT_BIN="/usr/bin/Rscript"

# ---- Temp -------------------------------------------------------------------
# Never /tmp, never $HOME.
TMPDIR="/mnt/data/shahar/.tmp"
export TMPDIR
