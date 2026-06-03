#!/usr/bin/env bash
# 00_setup_and_versions.sh
# Step 00 of 06_USED_interpro_domains.
# Idempotent setup: make sure the official EBI iprscan5 client and its xmltramp2
# dependency are present (vendored project-local, nothing under $HOME), and record
# tool/service versions for the Methods.
set -euo pipefail

export TMPDIR=/mnt/data/shahar/.tmp
mkdir -p "$TMPDIR"

BASE=/mnt/data/shahar/gwas_barley/morexV3_analysis/06_USED_interpro_domains
SCRIPTS="$BASE/scripts"
PYLIB="$SCRIPTS/pylib"
LOGS="$BASE/logs"
mkdir -p "$PYLIB" "$LOGS"

CLIENT="$SCRIPTS/iprscan5.py"
CLIENT_URL="https://raw.githubusercontent.com/ebi-wp/webservice-clients/master/python/iprscan5.py"
BASEURL="https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
VERS="$LOGS/versions.txt"

echo "=== Step 00: setup + versions ==="

# Official EBI iprscan5 REST client (download if missing).
if [ ! -s "$CLIENT" ]; then
  echo "Fetching official iprscan5 client ..."
  curl -sS --connect-timeout 30 --max-time 120 -o "$CLIENT" "$CLIENT_URL"
fi
echo "Client present: $CLIENT ($(stat -c%s "$CLIENT") bytes)"

# Local patch (idempotent): the upstream client's restRequest() returns bytes on
# the HTTPError fallback (requests .content) and then crashes on 'str + bytes'
# during status polling. Decode that fallback to str. Applied only if missing.
if ! grep -q "Local patch:" "$CLIENT"; then
  echo "Applying restRequest bytes->str patch ..."
  python3 - "$CLIENT" <<'PY'
import sys, io
p = sys.argv[1]
src = io.open(p, encoding="utf-8").read()
old = "    except HTTPError as ex:\n        result = requests.get(url).content\n"
new = ("    except HTTPError as ex:\n"
       "        # Local patch: decode requests .content (bytes) to str to match\n"
       "        # the success path; keep bytes only for genuinely binary bodies.\n"
       "        resp = requests.get(url).content\n"
       "        try:\n"
       "            result = resp.decode(u'utf-8')\n"
       "        except (UnicodeDecodeError, AttributeError):\n"
       "            result = resp\n")
if old in src:
    io.open(p, "w", encoding="utf-8").write(src.replace(old, new, 1))
    print("  patched")
else:
    print("  WARNING: expected buggy block not found (client may have changed)")
PY
else
  echo "restRequest patch already present."
fi

# xmltramp2 dependency, vendored under scripts/pylib (NOT under $HOME).
if ! PYTHONPATH="$PYLIB" python3 -c "import xmltramp2" 2>/dev/null; then
  echo "Installing xmltramp2 into project-local pylib ..."
  python3 -m pip install --quiet --target "$PYLIB" xmltramp2
fi
echo "xmltramp2 import: $(PYTHONPATH="$PYLIB" python3 -c 'import xmltramp2,os;print("OK", os.path.dirname(xmltramp2.__file__))')"

# Record versions / provenance.
{
  echo "06_USED_interpro_domains - tool/service versions"
  echo "date_run:               $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  echo "EBI iprscan5 client URL: $CLIENT_URL"
  echo "EBI iprscan5 client rev: $(grep -E "^version = " "$CLIENT" | head -1 | sed "s/version = u\?'//; s/'.*//")"
  echo "EBI REST baseUrl:        $BASEURL"
  echo "python3:                 $(python3 --version 2>&1)"
  echo "requests:                $(python3 -c 'import requests;print(requests.__version__)' 2>/dev/null)"
  echo "xmltramp2:               $(PYTHONPATH="$PYLIB" python3 -c 'import xmltramp2,importlib.metadata as m;print(m.version("xmltramp2"))' 2>/dev/null)"
  # InterProScan engine version, taken from any retrieved JSON if present.
  J=$(ls "$BASE"/intermediates/raw/*.json.json 2>/dev/null | head -1 || true)
  if [ -n "${J:-}" ]; then
    echo "interproscan-version:    $(python3 -c "import json,sys;print(json.load(open('$J')).get('interproscan-version'))" 2>/dev/null)"
  else
    echo "interproscan-version:    (recorded after first result is retrieved in step 01)"
  fi
} > "$VERS"

echo "Recorded -> $VERS"
cat "$VERS"
