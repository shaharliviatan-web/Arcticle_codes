#!/usr/bin/env bash
# 05_publication_plots.sh -- render publication PNG+PDF for a curated final-gene list.
# Usage: bash 05_publication_plots.sh [config.yaml] [final_genes.tsv]
# Create the curated list first (columns: trait gene_id MGmin epsilon [title_override]),
# e.g. by copying the chosen rows from gene_shortlist.csv.

set -euo pipefail
export TMPDIR=/mnt/data/shahar/.tmp
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="${1:-$script_dir/../00_config/config.yaml}"
FINAL="${2:-$script_dir/../00_config/final_genes.tsv}"
/usr/bin/Rscript "$script_dir/05_publication_plots.R" "$CFG" "$FINAL"
