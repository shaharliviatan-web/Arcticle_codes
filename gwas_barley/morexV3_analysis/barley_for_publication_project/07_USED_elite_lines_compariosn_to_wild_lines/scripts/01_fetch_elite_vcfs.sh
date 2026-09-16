#!/bin/bash
# =============================================================================
# 01_fetch_elite_vcfs.sh -- pull the elite CANDIDATE POOL for the three genes
# =============================================================================
# Replaces the manual "custom export" downloads of the earlier analysis with a
# reproducible API call, so the window and the line set are recorded rather than
# remembered.
#
# WHAT IS DOWNLOADED (changed 2026-09-14)
#   Not only the five configured lines, but the whole candidate POOL defined in
#   config/params.sh (ELITE_POOL_ACCESSION_TYPE + ELITE_POOL_ANNUALITY, i.e. the
#   136 spring elite lines). 02_screen_elite_lines.R scores every pool line for
#   call rate, and 03_build_matrices.R subsets the configured lines from this same
#   file -- so the screen and the figures are guaranteed to use identical data,
#   and changing the line selection never needs a new download.
#   The configured lines are asserted to be members of the pool.
#
# Endpoint (verified 2026-09-10):
#     POST <base>/vcf_export        Content-Type: application/x-www-form-urlencoded
#     chrom=chr3H&startpos=..&endpos=..&samples=["SAMEA..",..]
#   A JSON body returns HTTP 500 -- it must be form-encoded.
#
# THE ONE TRAP, and why this script is longer than it looks:
#   /vcf_export returns the CORRECT NUMBER of variants but reads them from a
#   SHIFTED position range, and the drift is not constant. Measured:
#     request 546,472,709-546,478,564 -> 158 records spanning 546,476,514-546,481,290
#     request 546,473,000-546,474,000 ->  52 records spanning 546,478,016-546,480,350
#   A naive call therefore returns real data for the WRONG PART OF THE GENOME,
#   silently. The fix, verified to recover 158/158 of the sites in the earlier
#   manual export with zero missing:
#     1. request a PADDED window (DIVBROWSE_PAD_BP either side),
#     2. ASSERT the returned span brackets the target window on BOTH sides,
#        widening the pad and retrying while it does not,
#     3. trim locally with bcftools to the exact gene +/- WINDOW_BP window.
#   Two further guards, both added after the failure actually occurred:
#     - a curl timeout can leave a TRUNCATED file that still has a valid header;
#       vcf_ok() checks every data line carries 9 + n_samples fields.
#     - `grep | head` under `set -o pipefail` dies of SIGPIPE; first/last/count
#       are read with awk directly from the file instead.
#
# The window is taken from step 04's gene_windows.tsv, and WINDOW_BP is asserted
# to equal the window step 04 actually used -- otherwise the wild and elite
# barcodes would not describe the same stretch of sequence.
#
# Output -> intermediates/elite_vcfs_raw/<gene>.pool.padded.vcf          (as returned)
#           intermediates/elite_vcfs_trimmed/<gene>.pool.vcf.gz (+ .csi)  (exact window)
#           inputs/elite_pool_samples.tsv                                (the pool)
#           results/tables/elite_vcf_provenance.tsv                      (audit trail)
#
# RESUMABLE: a cached padded download is reused if it passes vcf_ok(), carries the
# full pool, and brackets its window. FORCE_REFETCH=1 overrides.
#
# Author : Shahar Liviatan
# Created: 2026-09-10   Pool download: 2026-09-14
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/params.sh"

mkdir -p "$DIR_ELITE_RAW" "$DIR_ELITE_TRIMMED" "$DIR_TABLES" "$DIR_LOGS" "$DIR_INPUTS"

PANEL_TSV="${DIR_INPUTS}/panel_metadata.tsv"
POOL_TSV="${DIR_INPUTS}/elite_pool_samples.tsv"
PROV="${DIR_TABLES}/elite_vcf_provenance.tsv"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# vcf_ok <file> <n_samples>
# A curl that times out mid-transfer leaves a file that still has a valid #CHROM
# header and plausible first/last positions, so a first/last check alone will
# happily accept a TRUNCATED download (this bit us once). Verify instead that the
# header declares exactly n_samples and that EVERY data line carries 9+n fields.
vcf_ok () {
  awk -F'\t' -v ns="$2" '
    /^#CHROM/ { hdr=1; if (NF-9 != ns) bad=1; next }
    /^#/      { next }
                { n++; if (NF != ns+9) bad=1 }
    END       { if (!hdr || bad || n==0) exit 1; exit 0 }
  ' "$1"
}

# ---------------------------------------------------------------------------
# The candidate pool, from the harvested BioSamples metadata (step 00)
# ---------------------------------------------------------------------------
[ -s "$PANEL_TSV" ] || { echo "ERROR: ${PANEL_TSV} missing -- run 00_panel_metadata.sh" >&2; exit 1; }

awk -F'\t' -v at="$ELITE_POOL_ACCESSION_TYPE" -v an="$ELITE_POOL_ANNUALITY" '
  NR==1 {for(i=1;i<=NF;i++) h[$i]=i; print "sample_id\taccession_name\tbreeder\tyear_of_release\tannuality"; next}
  $h["accession_type"]==at && $h["annuality"]==an {
    print $h["sample_id"]"\t"$h["accession_name"]"\t"$h["breeder"]"\t"$h["year_of_release"]"\t"$h["annuality"]
  }' "$PANEL_TSV" > "$POOL_TSV"

N_SAMPLES=$(( $(wc -l < "$POOL_TSV") - 1 ))
log "Candidate pool: ${N_SAMPLES} lines (accession type '${ELITE_POOL_ACCESSION_TYPE}', '${ELITE_POOL_ANNUALITY}') -> ${POOL_TSV}"
[ "$N_SAMPLES" -ge 1 ] || { echo "ERROR: empty candidate pool" >&2; exit 1; }

# every configured line must be in the pool
MISSING=$(awk -F'\t' 'NR==FNR{if(FNR>1) pool[$1]=1; next}
                      !/^#/ && NF>1 && $1!="line_name" && !($2 in pool){print $1" ("$2")"}' \
          "$POOL_TSV" "$ELITE_LINES_TSV")
if [ -n "$MISSING" ]; then
  echo "ERROR: configured elite lines not in the candidate pool:" >&2; echo "$MISSING" >&2; exit 1
fi
log "All configured lines in config/elite_lines.tsv are members of the pool"

SAMPLES_JSON=$(awk -F'\t' 'NR>1 {printf "%s\"%s\"", (n++?",":""), $1} END{print ""}' "$POOL_TSV")
SAMPLES_JSON="[${SAMPLES_JSON}]"

printf 'gene_id\tshort_name\tchr\tgene_start\tgene_end\twin_start\twin_end\tpad_bp\treturned_first\treturned_last\tn_returned\tn_in_window\tn_indel\tn_multiallelic\tn_samples\tfetched_utc\n' > "$PROV"

# ---------------------------------------------------------------------------
# Assert our WINDOW_BP matches the window step 04 actually used
# ---------------------------------------------------------------------------
for entry in "${TARGET_GENES[@]}"; do
  GENE_ID="${entry%%:*}"
  awk -F'\t' -v g="$GENE_ID" -v w="$WINDOW_BP" '
    NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
    $h["gene_id"]==g {
      lo=$h["gene_start"]-$h["win_start"]; hi=$h["win_end"]-$h["gene_end"];
      if (lo!=w || hi!=w) {
        printf "ERROR: %s step-04 window is -%d/+%d bp but params.sh WINDOW_BP=%d\n", g, lo, hi, w > "/dev/stderr";
        exit 1
      }
      found=1
    }
    END { if (!found) { printf "ERROR: %s not found in gene_windows.tsv\n", g > "/dev/stderr"; exit 1 } }
  ' "$STEP04_GENE_WINDOWS"
done
log "WINDOW_BP=${WINDOW_BP} matches step 04 for all target genes"

# ---------------------------------------------------------------------------
# Per gene: padded fetch -> integrity + bracket assertion -> trim
# ---------------------------------------------------------------------------
for entry in "${TARGET_GENES[@]}"; do
  GENE_ID="${entry%%:*}"
  SHORT="${entry##*:}"

  read -r CHR GSTART GEND WSTART WEND < <(awk -F'\t' -v g="$GENE_ID" '
    NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
    $h["gene_id"]==g {print $h["chr"], $h["gene_start"], $h["gene_end"], $h["win_start"], $h["win_end"]; exit}
  ' "$STEP04_GENE_WINDOWS")

  DVB_CHR="chr${CHR}"     # step 04 uses "3H"; DivBrowse uses "chr3H"

  log "-------------------------------------------------------------"
  log "${SHORT} (${GENE_ID})  ${CHR}:${WSTART}-${WEND}  [gene ${GSTART}-${GEND}]"

  PAD="$DIVBROWSE_PAD_BP"
  RAW="${DIR_ELITE_RAW}/${GENE_ID}.pool.padded.vcf"
  OK=0

  # --- resume: reuse a cached padded download if it is complete and brackets ---
  if [ "${FORCE_REFETCH:-0}" != "1" ] && [ -s "$RAW" ] && vcf_ok "$RAW" "$N_SAMPLES"; then
    CF=$(awk '!/^#/{print $2; exit}' "$RAW")
    CL=$(awk '!/^#/{p=$2} END{print p}' "$RAW")
    if [ -n "${CF:-}" ] && [ -n "${CL:-}" ] && [ "$CF" -le "$WSTART" ] && [ "$CL" -ge "$WEND" ]; then
      log "  cached download reused: ${CF}-${CL}, ${N_SAMPLES} samples (FORCE_REFETCH=1 to re-download)"
      FIRST="$CF"; LAST="$CL"; NRET=$(awk '!/^#/{n++} END{print n+0}' "$RAW"); OK=1
    fi
  fi

  while [ "$OK" -ne 1 ] && [ "$PAD" -le "$DIVBROWSE_PAD_MAX_BP" ]; do
    REQ_START=$(( WSTART - PAD )); [ "$REQ_START" -lt 1 ] && REQ_START=1
    REQ_END=$(( WEND + PAD ))

    log "  requesting pad=${PAD} bp  (${REQ_START}-${REQ_END}), ${N_SAMPLES} samples ..."
    curl -s -m "$DIVBROWSE_TIMEOUT" -X POST "${DIVBROWSE_BASE}/vcf_export" \
      --data-urlencode "chrom=${DVB_CHR}" \
      --data-urlencode "startpos=${REQ_START}" \
      --data-urlencode "endpos=${REQ_END}" \
      --data-urlencode "samples=${SAMPLES_JSON}" \
      -o "$RAW" || true
    # curl failure (timeout, reset) must NOT abort under `set -e`: fall through to
    # the integrity/bracket checks below, which drive the widen-and-retry loop.

    if ! vcf_ok "$RAW" "$N_SAMPLES"; then
      log "  no usable VCF returned (missing header, wrong sample count, or truncated); retrying"
      rm -f "$RAW"
      PAD=$(( PAD * DIVBROWSE_RETRY_FACTOR )); continue
    fi

    FIRST=$(awk '!/^#/{print $2; exit}' "$RAW")
    LAST=$(awk '!/^#/{p=$2} END{print p}' "$RAW")
    NRET=$(awk '!/^#/{n++} END{print n+0}' "$RAW")

    # THE ASSERTION: returned span must bracket the target window on both sides.
    if [ "$FIRST" -le "$WSTART" ] && [ "$LAST" -ge "$WEND" ]; then
      log "  returned ${NRET} records spanning ${FIRST}-${LAST}  -> brackets target, OK"
      OK=1; break
    fi
    log "  returned span ${FIRST}-${LAST} does NOT bracket ${WSTART}-${WEND}; widening"
    PAD=$(( PAD * DIVBROWSE_RETRY_FACTOR ))
  done

  if [ "$OK" -ne 1 ]; then
    echo "ERROR: ${GENE_ID}: could not obtain a complete window bracketing ${WSTART}-${WEND} up to pad ${DIVBROWSE_PAD_MAX_BP} bp." >&2
    echo "       Refusing to continue -- a non-bracketing export means silent data loss." >&2
    exit 1
  fi

  # --- normalise: DivBrowse writes FILTER=NA and no contig line -------------
  TRIMMED="${DIR_ELITE_TRIMMED}/${GENE_ID}.pool.vcf.gz"
  WORK="${DIR_INTERMEDIATES}/.work_${GENE_ID}"
  mkdir -p "$WORK"

  grep '^##' "$RAW" > "${WORK}/hdr.vcf" || true
  echo "##contig=<ID=${DVB_CHR}>" >> "${WORK}/hdr.vcf"
  grep '^#CHROM' "$RAW" >> "${WORK}/hdr.vcf"
  grep -v '^#' "$RAW" | awk 'BEGIN{FS=OFS="\t"} {$7="."; print}' > "${WORK}/body.vcf"

  echo -e "${DVB_CHR}\t${CHR}" > "${WORK}/chr_map.txt"

  # -t (targets), not -r (regions): -r needs a random-access index, which a
  # stream on stdin does not have. -t streams and filters, which is what we want.
  cat "${WORK}/hdr.vcf" "${WORK}/body.vcf" \
    | "$BCFTOOLS_BIN" view -Ou \
    | "$BCFTOOLS_BIN" annotate --rename-chrs "${WORK}/chr_map.txt" -Ou \
    | "$BCFTOOLS_BIN" view -t "${CHR}:${WSTART}-${WEND}" -Oz -o "$TRIMMED"
  "$BCFTOOLS_BIN" index -f "$TRIMMED"

  N_IN=$("$BCFTOOLS_BIN" query -f '%POS\n' "$TRIMMED" | wc -l)
  N_INDEL=$("$BCFTOOLS_BIN" query -f '%REF\t%ALT\n' "$TRIMMED" \
            | awk -F'\t' '{n=split($2,a,","); ind=(length($1)>1); for(i=1;i<=n;i++) if(length(a[i])>1) ind=1; if(ind) c++} END{print c+0}')
  N_MULTI=$("$BCFTOOLS_BIN" query -f '%ALT\n' "$TRIMMED" | awk '/,/{n++} END{print n+0}')

  log "  trimmed to ${CHR}:${WSTART}-${WEND} -> ${N_IN} records (${N_INDEL} indel, ${N_MULTI} multiallelic)"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$GENE_ID" "$SHORT" "$CHR" "$GSTART" "$GEND" "$WSTART" "$WEND" "$PAD" \
    "$FIRST" "$LAST" "$NRET" "$N_IN" "$N_INDEL" "$N_MULTI" "$N_SAMPLES" \
    "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$PROV"

  rm -rf "$WORK"
done

log "============================================================="
log "DONE. Provenance -> ${PROV}"
column -t -s$'\t' "$PROV"
