#!/bin/bash
# =============================================================================
# 00_panel_metadata.sh -- what every sample in the DivBrowse panel actually IS
# =============================================================================
# The DivBrowse barley pangenome v2 instance labels most samples by an IPK
# genebank number (HOR #####), not by cultivar name. That is why the earlier,
# pre-publication version of this analysis could not verify four of its six
# "elite" lines. This script removes the guesswork:
#
#   1. scrape the SAMEA <-> panel-label mapping that the DivBrowse landing page
#      embeds as a `sampleIdMapping` JavaScript array (all 1315 samples), and
#   2. pull the EBI BioSamples record for each SAMEA accession, which carries the
#      real `accession name`, `accession type` (elite lines / precision
#      collection), organism, breeder, year of release and growth habit.
#
# Output -> inputs/panel_metadata.tsv        one row per panel sample
#           inputs/divbrowse_sample_ids.tsv  raw SAMEA <-> panel label mapping
#
# RESUMABLE: one cached JSON per accession under intermediates/biosamples_cache/.
# Re-running only fills gaps. A full cold run is ~1315 requests, a few minutes.
#
# This step is READ-ONLY against public APIs and is not needed to reproduce the
# figures if inputs/panel_metadata.tsv is already present -- it exists so the
# elite-line selection in config/elite_lines.tsv is auditable rather than asserted.
#
# Author : Shahar Liviatan
# Created: 2026-09-10
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/params.sh"

CACHE_DIR="${DIR_INTERMEDIATES}/biosamples_cache"
mkdir -p "$CACHE_DIR" "$DIR_INPUTS" "$DIR_LOGS"

INDEX_HTML="${DIR_INTERMEDIATES}/divbrowse_index.html"
ID_MAP="${DIR_INPUTS}/divbrowse_sample_ids.tsv"
OUT_TSV="${DIR_INPUTS}/panel_metadata.tsv"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# 1. DivBrowse landing page -> SAMEA <-> panel label
# ---------------------------------------------------------------------------
log "Fetching DivBrowse landing page ..."
curl -s -m "$DIVBROWSE_TIMEOUT" "$DIVBROWSE_INDEX" -o "$INDEX_HTML"
[ -s "$INDEX_HTML" ] || { echo "ERROR: empty DivBrowse index page" >&2; exit 1; }

python3 - "$INDEX_HTML" "$ID_MAP" <<'PY'
import re, sys
html, out = sys.argv[1], sys.argv[2]
s = open(html, encoding="utf-8", errors="replace").read()
pairs = re.findall(r'\{"NAME":"([^"]+)","ID":"([^"]+)"\}', s)
if not pairs:
    sys.exit("ERROR: could not find sampleIdMapping in the DivBrowse page. "
             "The site layout may have changed -- check before trusting any label.")
with open(out, "w") as fh:
    fh.write("sample_id\tpanel_label\tpanel_group\n")
    for name, sid in pairs:
        grp = name.split(":")[0] if ":" in name else ""
        fh.write(f"{sid}\t{name}\t{grp}\n")
print(f"  mapped {len(pairs)} samples")
PY

N_IDS=$(( $(wc -l < "$ID_MAP") - 1 ))
log "DivBrowse panel size: ${N_IDS} samples -> ${ID_MAP}"

# ---------------------------------------------------------------------------
# 2. EBI BioSamples record per accession (cached, parallel)
# ---------------------------------------------------------------------------
FETCH_ONE="${DIR_INTERMEDIATES}/.fetch_one_biosample.sh"
cat > "$FETCH_ONE" <<EOS
#!/bin/bash
id="\$1"
out="${CACHE_DIR}/\${id}.json"
[ -s "\$out" ] && exit 0
curl -s -m 45 -H "Accept: application/json" "${BIOSAMPLES_BASE}/\${id}" -o "\$out"
EOS
chmod +x "$FETCH_ONE"

log "Fetching BioSamples records (parallel ${BIOSAMPLES_PARALLEL}, cached) ..."
tail -n +2 "$ID_MAP" | cut -f1 | xargs -P "$BIOSAMPLES_PARALLEL" -n 1 "$FETCH_ONE"
log "Cache now holds $(ls "$CACHE_DIR" | wc -l) records"

# ---------------------------------------------------------------------------
# 3. Flatten to one table
# ---------------------------------------------------------------------------
python3 - "$ID_MAP" "$CACHE_DIR" "$OUT_TSV" <<'PY'
import json, os, sys, csv, re

id_map, cache, out = sys.argv[1], sys.argv[2], sys.argv[3]

rows = []
with open(id_map) as fh:
    next(fh)
    for line in fh:
        sid, label, grp = line.rstrip("\n").split("\t")
        rec = {"sample_id": sid, "panel_label": label, "panel_group": grp}
        p = os.path.join(cache, sid + ".json")
        ch = {}
        name = ""
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                d = json.load(open(p))
                ch = d.get("characteristics", {}) or {}
                name = d.get("name", "") or ""
            except Exception:
                pass
        g = lambda k: (ch.get(k, [{}])[0].get("text", "") if ch.get(k) else "")
        biomat = g("bio material") or g("biological material id")
        m = re.search(r"HOR[: ]\s*(\d+)", biomat or "")
        rec.update(
            biosample_name       = name,
            accession_name       = g("accession name"),
            accession_type       = g("accession type"),
            organism             = g("organism"),
            infraspecific_name   = g("infraspecific name"),
            shape_panel          = g("shape panel"),
            project_name         = g("project name"),
            row_type             = g("row type"),
            annuality            = g("annuality"),
            breeder              = g("breeder"),
            year_of_release      = g("year of release"),
            country              = g("material source geographic location")
                                   or g("geographic location (country and/or sea)"),
            biological_material  = biomat,
            hor_number           = m.group(1) if m else "",
        )
        rows.append(rec)

cols = ["sample_id", "panel_label", "panel_group", "biosample_name", "accession_name",
        "accession_type", "organism", "infraspecific_name", "shape_panel",
        "project_name", "row_type", "annuality", "breeder", "year_of_release",
        "country", "biological_material", "hor_number"]

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

from collections import Counter
print(f"  wrote {len(rows)} rows -> {out}")
print("  accession_type:", dict(Counter(r["accession_type"] or "(none)" for r in rows)))
wild = [r for r in rows if "spontaneum" in (r["organism"] + r["infraspecific_name"]).lower()]
print(f"  wild (ssp. spontaneum) accessions in panel: {len(wild)}")
PY

# ---------------------------------------------------------------------------
# 4. Verify every configured elite line resolves, and is what config claims
# ---------------------------------------------------------------------------
log "Verifying config/elite_lines.tsv against the harvested metadata ..."
python3 - "$ELITE_LINES_TSV" "$OUT_TSV" <<'PY'
import sys, csv
cfg_path, meta_path = sys.argv[1], sys.argv[2]

cfg = [l for l in open(cfg_path) if not l.startswith("#") and l.strip()]
cfg = list(csv.DictReader(cfg, delimiter="\t"))
meta = {r["sample_id"]: r for r in csv.DictReader(open(meta_path), delimiter="\t")}

bad = []
for c in cfg:
    m = meta.get(c["sample_id"])
    if m is None:
        bad.append(f"{c['line_name']}: {c['sample_id']} is not in the panel"); continue
    if m["accession_name"] != c["line_name"]:
        bad.append(f"{c['line_name']}: BioSamples accession name is '{m['accession_name']}'")
    if m["accession_type"] != "elite lines":
        bad.append(f"{c['line_name']}: accession type is '{m['accession_type']}', not 'elite lines'")
    if m["year_of_release"] != c["year_of_release"]:
        bad.append(f"{c['line_name']}: year {m['year_of_release']} != config {c['year_of_release']}")
    print(f"  OK  {c['line_name']:12s} {c['sample_id']:16s} "
          f"{m['accession_type']:14s} {m['year_of_release']:5s} {m['annuality']:12s} {m['breeder']}")

if bad:
    print("\nELITE LINE VERIFICATION FAILED:")
    for b in bad: print("   -", b)
    sys.exit(1)
print(f"\n  all {len(cfg)} configured elite lines verified against EBI BioSamples")
PY

log "DONE -> ${OUT_TSV}"
