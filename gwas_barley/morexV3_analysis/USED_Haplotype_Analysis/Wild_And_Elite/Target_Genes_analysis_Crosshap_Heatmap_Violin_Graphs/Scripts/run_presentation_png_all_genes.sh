#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="$SCRIPT_DIR/make_one_gene_presentation_png_NO_EliteInViolinPlot.R"

# Format: "TRAIT|GENE_FILE|TITLE|EPSILON|MGMIN"
GENES=(
    "betaglucan|HORVU.MOREX.r3.1HG0077860.vcf.gz|NAC domain-containing transcription factor|0.8|3"
    "betaglucan|HORVU.MOREX.r3.2HG0112770.vcf.gz|bHLH (ACT-like) transcription factor|0.5|3"
    "betaglucan|HORVU.MOREX.r3.2HG0112910.vcf.gz|Pectinesterase-like (PME)|0.8|2"
    "betaglucan|HORVU.MOREX.r3.3HG0299440.vcf.gz|APETALA2-like protein 5 (AP2/ERF TF)|0.6|3"
    "fiber|HORVU.MOREX.r3.7HG0729010.vcf.gz|Peroxidase 45-like (Class III peroxidase)|0.5|2"
    "fiber|HORVU.MOREX.r3.7HG0642350.vcf.gz|Putative BAHD Acyltransferase|0.5|3"
    "starch|HORVU.MOREX.r3.3HG0301750.vcf.gz|Alpha-glucan phosphorylase|0.4|2"
    "starch|HORVU.MOREX.r3.3HG0301770.vcf.gz|Alpha-glucan phosphorylase alpha 1-4|0.8|3"
    "starch|HORVU.MOREX.r3.3HG0301710.vcf.gz|Phosphate transporter 43 / PHT / Anion transporter|0.85|2"
)

TOTAL=${#GENES[@]}
SUCCESS=0
FAIL=0

for entry in "${GENES[@]}"; do
    TRAIT=$(echo "$entry" | cut -d'|' -f1)
    GENE_FILE=$(echo "$entry" | cut -d'|' -f2)
    GENE_TITLE=$(echo "$entry" | cut -d'|' -f3)
    EPSILON=$(echo "$entry" | cut -d'|' -f4)
    MGMIN=$(echo "$entry" | cut -d'|' -f5)

    echo "------------------------------------------------------------"
    echo "[$((SUCCESS + FAIL + 1))/$TOTAL] $GENE_FILE ($TRAIT)"

    TMP_SCRIPT=$(mktemp /tmp/png_script_XXXXXX.R)
    sed \
        -e "s|^TRAIT <- .*|TRAIT <- \"${TRAIT}\"|" \
        -e "s|^GENE_FILE <- .*|GENE_FILE <- \"${GENE_FILE}\"|" \
        -e "s|^GENE_TITLE <- .*|GENE_TITLE <- \"${GENE_TITLE}\"|" \
        -e "s|^EPSILON <- .*|EPSILON <- ${EPSILON}|" \
        -e "s|^MGMIN <- .*|MGMIN <- ${MGMIN}|" \
        "$R_SCRIPT" > "$TMP_SCRIPT"

    if Rscript "$TMP_SCRIPT"; then
        echo "[SUCCESS] $GENE_FILE"
        SUCCESS=$((SUCCESS + 1))
    else
        echo "[FAILED]  $GENE_FILE"
        FAIL=$((FAIL + 1))
    fi

    rm -f "$TMP_SCRIPT"
done

echo "============================================================"
echo "DONE: $SUCCESS/$TOTAL succeeded, $FAIL failed."
