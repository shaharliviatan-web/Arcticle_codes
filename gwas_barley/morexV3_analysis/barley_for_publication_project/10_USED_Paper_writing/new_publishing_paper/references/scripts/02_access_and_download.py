#!/usr/bin/env python3
"""Step 3 of the citation audit: find a free full text for each tier-A/B work and download it.

For each work in meta/works.tsv with tier A or B (plus a few tier-C works that carry a specific
claim; EXTRA below): look up Europe PMC (PMCID, open-access flag) and OpenAlex (OA locations),
then try, in order: the Europe PMC rendered PDF (publisher PDF of PMC articles), OpenAlex PDF
URLs, and the Europe PMC full-text XML (PMC open-access subset; no page numbers).
A download counts only if it is a real PDF/XML whose text contains words of the title.
Files: pdfs/<key>.pdf or pdfs/<key>.xml; log: meta/access.tsv. Works already held locally
(LOCAL) are copied, not downloaded. Run from references/:  python3 scripts/02_access_and_download.py
"""
import csv, json, os, re, shutil, subprocess, sys, time, unicodedata
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PDFS = os.path.join(ROOT, "pdfs")
TMP = os.environ.get("TMPDIR", "/mnt/data/shahar/.tmp")
WRITING = os.path.dirname(os.path.dirname(ROOT))  # 10_USED_Paper_writing
LOCAL = {
    "PotapenkoMolEcol2026": os.path.join(WRITING, "reference_sources/Potapenko_2026_MolEcol/Potapenko_2026_MolEcol.pdf"),
    "PotapenkoMBE2026": os.path.join(WRITING, "reference_sources/Potapenko_2026_MBE/Potapenko_2026_MBE.pdf"),
}
EXTRA = {"Jayakodi2024", "Konig2023"}  # tier C, but the text makes a specific claim about them (TODO 230)
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"
S = requests.Session()
S.headers["User-Agent"] = UA


def norm(s):
    s = unicodedata.normalize("NFKD", s or "").encode("ascii", "ignore").decode().lower()
    return re.sub(r"[^a-z0-9 ]", " ", s)


def title_words(title):
    return [w for w in norm(title).split() if len(w) > 4][:8]


def text_of(path):
    if path.endswith(".pdf"):
        out = subprocess.run(["pdftotext", "-l", "3", path, "-"], capture_output=True, text=True)
        return out.stdout
    return re.sub(r"<[^>]+>", " ", open(path, errors="ignore").read())


def matches(path, title):
    words = title_words(title)
    txt = norm(text_of(path))
    hit = sum(w in txt for w in words)
    return hit >= max(2, int(0.6 * len(words))), f"{hit}/{len(words)}"


def fetch(url, dest, kind):
    try:
        r = S.get(url, timeout=60, allow_redirects=True)
    except Exception as e:
        return False, f"error {type(e).__name__}"
    if r.status_code != 200:
        return False, f"http {r.status_code}"
    head = r.content[:5]
    if kind == "pdf" and head != b"%PDF-":
        return False, f"not a pdf ({r.headers.get('content-type', '')[:30]})"
    if kind == "xml" and b"<article" not in r.content[:3000]:
        return False, "not article xml"
    open(dest, "wb").write(r.content)
    return True, f"{len(r.content) // 1024} kB"


def epmc(doi):
    try:
        r = S.get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                  params={"query": f'DOI:"{doi}"', "format": "json", "resultType": "core"}, timeout=30).json()
        res = r.get("resultList", {}).get("result", [])
        return res[0] if res else {}
    except Exception:
        return {}


def openalex(doi):
    try:
        return S.get(f"https://api.openalex.org/works/doi:{doi}", timeout=30).json()
    except Exception:
        return {}


works = [r for r in csv.DictReader(open(os.path.join(ROOT, "meta", "works.tsv")), delimiter="\t")]
cr = {r["key"]: r for r in csv.DictReader(open(os.path.join(ROOT, "meta", "crossref_summary.tsv")), delimiter="\t")}
log = []
for w in works:
    key, doi = w["key"], w["doi"]
    if w["tier"] not in ("A", "B") and key not in EXTRA:
        continue
    title = cr.get(key, {}).get("title", "")
    rec = {"key": key, "tier": w["tier"], "doi": doi, "pmcid": "", "epmc_oa": "", "oa_status": "",
           "file": "", "how": "", "title_check": "", "tried": ""}
    if key in LOCAL:
        dest = os.path.join(PDFS, key + ".pdf")
        shutil.copy(LOCAL[key], dest)
        rec.update(file=os.path.basename(dest), how="local copy (reference_sources/)",
                   title_check=matches(dest, title)[1])
        log.append(rec); print(key, "local"); continue
    if not doi:
        rec["how"] = "no DOI: manual"
        log.append(rec); print(key, "no doi"); continue
    e = epmc(doi)
    oa = openalex(doi)
    rec["pmcid"] = e.get("pmcid", "")
    rec["epmc_oa"] = e.get("isOpenAccess", "")
    rec["oa_status"] = (oa.get("open_access") or {}).get("oa_status", "")
    cands = []
    if rec["pmcid"]:
        cands.append(("pdf", f"https://europepmc.org/backend/ptpmcrender.fcgi?accid={rec['pmcid']}&blobtype=pdf", "Europe PMC PDF"))
    for loc in [oa.get("best_oa_location") or {}] + (oa.get("locations") or []):
        if loc.get("pdf_url"):
            cands.append(("pdf", loc["pdf_url"], "OA PDF (" + ((loc.get("source") or {}).get("display_name") or "repository")[:40] + ")"))
    if rec["pmcid"]:
        cands.append(("xml", f"https://www.ebi.ac.uk/europepmc/webservices/rest/{rec['pmcid']}/fullTextXML", "Europe PMC full-text XML"))
    seen, tried = set(), []
    for kind, url, how in cands:
        if url in seen:
            continue
        seen.add(url)
        dest = os.path.join(PDFS, f"{key}.{kind}")
        ok, msg = fetch(url, dest, kind)
        if ok:
            good, score = matches(dest, title)
            if good:
                rec.update(file=os.path.basename(dest), how=how, title_check=score)
                break
            os.remove(dest)
            msg = f"title mismatch {score}"
        tried.append(f"{how}: {msg}")
        time.sleep(0.5)
    rec["tried"] = " | ".join(tried)
    log.append(rec)
    print(key, rec["oa_status"], rec["pmcid"], "->", rec["file"] or "NONE", rec["how"], rec["title_check"])
    time.sleep(0.5)

cols = ["key", "tier", "doi", "pmcid", "epmc_oa", "oa_status", "file", "how", "title_check", "tried"]
with open(os.path.join(ROOT, "meta", "access.tsv"), "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
    wr.writeheader()
    wr.writerows(log)
