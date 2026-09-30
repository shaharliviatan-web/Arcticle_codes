#!/usr/bin/env bash
# run_all.sh - reproduce the whole 7H shared-signal branch (~2 min).
set -euo pipefail
S="$(dirname "$(readlink -f "$0")")"; source "$S/config.sh"; mkdir -p "$BRANCH/logs"
bash     "$S/00_make_gene_window_vcfs.sh"        2>&1 | tee "$BRANCH/logs/00.log"
Rscript  "$S/01_gene_position_and_distances.R"   2>&1 | tee "$BRANCH/logs/01.log"
Rscript  "$S/02_local_association_profile.R"     2>&1 | tee "$BRANCH/logs/02.log"
bash     "$S/03_annotate_gene.sh"                2>&1 | tee "$BRANCH/logs/03.log"
bash     "$S/03b_interpro_and_nr.sh"          2>&1 | tee "$BRANCH/logs/03b.log"
Rscript  "$S/04_genotype_phenotype_split.R"      2>&1 | tee "$BRANCH/logs/04.log"
Rscript  "$S/05_crosshap_and_epsilon_sweep.R"    2>&1 | tee "$BRANCH/logs/05.log"
Rscript  "$S/06_make_crosshap_figures.R"         2>&1 | tee "$BRANCH/logs/06.log"
echo "run_all complete -> $BRANCH/results"
