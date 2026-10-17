#!/usr/bin/env python3
"""Convert the full texts in pdfs/ to plain text in text/<key>.txt, for the audit (step 5).

PDF -> one block per page, headed '=== page N (printed p. X) ===' (N = PDF page; the printed page number is taken
from the journal's page range when the first page is known, else left out). XML (Europe PMC JATS) -> text with the
section titles as '=== section: ... ===' (no page numbers). Run from references/: python3 scripts/04_to_text.py
"""
import csv, os, re, subprocess, html
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
cr = {r["key"]: r for r in csv.DictReader(open(os.path.join(ROOT, "meta", "crossref_summary.tsv")), delimiter="\t")}
for f in sorted(os.listdir(os.path.join(ROOT, "pdfs"))):
    key, ext = os.path.splitext(f)
    src, dst = os.path.join(ROOT, "pdfs", f), os.path.join(ROOT, "text", key + ".txt")
    if ext == ".pdf":
        n = int(re.search(r"Pages:\s+(\d+)", subprocess.run(["pdfinfo", src], capture_output=True, text=True).stdout).group(1))
        first = re.match(r"(\d+)-", cr.get(key, {}).get("page", "") or "")
        out = []
        for p in range(1, n + 1):
            t = subprocess.run(["pdftotext", "-f", str(p), "-l", str(p), src, "-"], capture_output=True, text=True).stdout
            out.append(f"=== page {p} ===\n{t}")
        open(dst, "w").write("\n".join(out))
    elif ext == ".xml":
        x = open(src, encoding="utf8").read()
        x = re.sub(r"<ref-list.*?</ref-list>", "\n=== reference list omitted ===\n", x, flags=re.S)
        x = re.sub(r"<title>(.*?)</title>", lambda m: "\n=== section: " + re.sub(r"<[^>]+>", "", m.group(1)) + " ===\n", x, flags=re.S)
        x = re.sub(r"</p>|<p[ >]", "\n\n", x)
        x = re.sub(r"<xref[^>]*>(.*?)</xref>", r"\1", x, flags=re.S)
        x = html.unescape(re.sub(r"<[^>]+>", "", x))
        x = re.sub(r"\n{3,}", "\n\n", x)
        open(dst, "w").write(x)
    elif ext == ".md":  # full text saved as markdown (e.g. from the publisher's HTML); no page numbers
        open(dst, "w").write(open(src, encoding="utf8").read())
    else:
        continue
    print(key, ext, os.path.getsize(dst))
