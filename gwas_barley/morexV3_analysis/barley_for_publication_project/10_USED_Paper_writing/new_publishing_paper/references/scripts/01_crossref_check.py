#!/usr/bin/env python3
"""Step 2 of the citation audit: bibliographic check of every cited work against Crossref.

For each row of meta/works.tsv: fetch the Crossref record of the DOI (if given) and run a
bibliographic search (top 3 hits). Raw JSON goes to meta/crossref/<key>.json; a summary table to
meta/crossref_summary.tsv. Nothing in the manuscript is read or changed.
Run from references/:  python3 scripts/01_crossref_check.py
"""
import csv, json, os, sys, time, unicodedata
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "meta", "crossref")
os.makedirs(OUT, exist_ok=True)
S = requests.Session()
S.headers["User-Agent"] = "citation-audit-script/1.0"


def ascii_lower(s):
    return unicodedata.normalize("NFKD", s or "").encode("ascii", "ignore").decode().lower()


def get(url, params=None):
    for _ in range(3):
        try:
            r = S.get(url, params=params, timeout=30)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            return r.json()["message"]
        except Exception as e:  # retry on network hiccups
            err = e
            time.sleep(2)
    print("FAILED", url, err, file=sys.stderr)
    return None


def summarize(m):
    if not m:
        return {}
    def yr(k):
        d = m.get(k, {}).get("date-parts", [[None]])
        return d[0][0] if d and d[0] else None
    authors = m.get("author", [])
    return {
        "doi": m.get("DOI", ""),
        "type": m.get("type", ""),
        "title": " ".join(m.get("title", [""])),
        "first_author": authors[0].get("family", authors[0].get("name", "")) if authors else "",
        "n_authors": len(authors),
        "authors": "; ".join(f'{a.get("family", a.get("name", ""))} {a.get("given", "")}'.strip() for a in authors),
        "journal": " ".join(m.get("container-title", [""])),
        "journal_short": " ".join(m.get("short-container-title", [""])),
        "volume": m.get("volume", ""),
        "issue": m.get("issue", ""),
        "page": m.get("page", ""),
        "article_number": m.get("article-number", ""),
        "year_issued": yr("issued"),
        "year_print": yr("published-print"),
        "year_online": yr("published-online"),
        "publisher": m.get("publisher", ""),
        "isbn": "; ".join(m.get("ISBN", [])),
    }


rows = list(csv.DictReader(open(os.path.join(ROOT, "meta", "works.tsv")), delimiter="\t"))
summary = []
for r in rows:
    key = r["key"]
    rec = {"key": key, "cite": r["cite"], "tier": r["tier"], "doi_given": r["doi"]}
    by_doi = get(f"https://api.crossref.org/works/{r['doi']}") if r["doi"] else None
    hits = []
    if r["query"]:
        q = get("https://api.crossref.org/works", {"query.bibliographic": r["query"], "rows": 3})
        hits = q.get("items", []) if q else []
    json.dump({"by_doi": by_doi, "search": hits}, open(os.path.join(OUT, key + ".json"), "w"), indent=1)
    best = by_doi or (hits[0] if hits else None)
    s = summarize(best)
    rec.update({"source": "doi" if by_doi else ("search" if hits else "none")})
    rec.update(s)
    # does the search's top hit agree with the DOI record?
    if by_doi and hits:
        rec["search_top_doi"] = hits[0].get("DOI", "")
    # expected first author / year from the in-text citation
    exp_author = ascii_lower(r["cite"].replace("[TODO 411 candidate] ", "").split(" ")[0]).rstrip(",")
    exp_year = "".join(ch for ch in r["cite"] if ch.isdigit())[-4:]
    rec["author_match"] = exp_author[:5] in ascii_lower(s.get("first_author", "")) if s else ""
    years = {str(s.get(k)) for k in ("year_issued", "year_print", "year_online")} if s else set()
    rec["year_match"] = exp_year in years if s else ""
    summary.append(rec)
    print(key, rec["source"], s.get("doi", ""), rec["author_match"], rec["year_match"], s.get("title", "")[:70])
    time.sleep(0.3)

cols = ["key", "cite", "tier", "source", "doi_given", "doi", "search_top_doi", "author_match", "year_match",
        "first_author", "n_authors", "year_issued", "year_print", "year_online", "title", "journal",
        "journal_short", "volume", "issue", "page", "article_number", "type", "publisher", "isbn", "authors"]
with open(os.path.join(ROOT, "meta", "crossref_summary.tsv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(summary)
