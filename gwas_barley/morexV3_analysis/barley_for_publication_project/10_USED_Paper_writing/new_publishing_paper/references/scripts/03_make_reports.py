#!/usr/bin/env python3
"""Regenerate the status tables of the citation audit from meta/*.tsv and the files in pdfs/.

Writes the auto-generated blocks of citation_audit.md (claim inventory with source status) and
to_download.md (what is still missing). Run again after adding PDFs to pdfs/:
    python3 scripts/03_make_reports.py        (from references/)
Hand-written text (findings, verdicts) lives outside the AUTO blocks and is kept.
"""
import csv, html, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
works = {r["key"]: r for r in csv.DictReader(open(os.path.join(ROOT, "meta", "works.tsv")), delimiter="\t")}
cr = {r["key"]: r for r in csv.DictReader(open(os.path.join(ROOT, "meta", "crossref_summary.tsv")), delimiter="\t")}
claims = list(csv.DictReader(open(os.path.join(ROOT, "meta", "claims.tsv")), delimiter="\t"))
VPATH = os.path.join(ROOT, "meta", "verdicts.tsv")  # step 5: one row per claim x work (id, key, verdict, where, quote, note)
verdicts = list(csv.DictReader(open(VPATH), delimiter="\t")) if os.path.exists(VPATH) else []
VD = {(v["id"], v["key"]): v for v in verdicts}

# free copies the scripts could not fetch (publisher / PMC bot checks): readable in a browser without login
FREE = {
    "Badr2000": "https://doi.org/10.1093/oxfordjournals.molbev.a026330 (free at OUP)",
    "Bogard2010": "https://doi.org/10.1093/jxb/erq238 (free at OUP)",
    "Burton2008": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2287353/",
    "Dempewolf2017": "https://doi.org/10.2135/cropsci2016.10.0885 (open access, Wiley)",
    "Geng2021": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10236047/ (author manuscript; the Springer version via your institution is better)",
    "Girard2012": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3426136/",
    "Irigoyen2011": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3327177/",
    "Nevo2010": "https://doi.org/10.1111/j.1365-3040.2009.02107.x (free, Wiley)",
    "Petit2016": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4902622/",
    "Yang2010": "https://pmc.ncbi.nlm.nih.gov/articles/PMC2900678/",
    "Zhang2019": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6533017/",
}
# Crossref holds an early-online or volume-less record for these; values checked in Europe PMC / OpenAlex
PAGES = {"Petit2016": "171:894-913", "Holland2003": "22:9-112", "Keeling1993": "191:342-348"}
OPTIONAL = {"Farag2022", "Godfray2010", "Harlan1966", "Hubner2009", "Hubner2013"}  # tier B, abstract saved


def held(key):
    if works[key]["tier"] == "C" and key not in ("Marsh2023", "Jayakodi2024", "Konig2023", "Tomczak2014"):
        return "tier C: bibliographic check only"
    for ext, label in (("pdf", "PDF"), ("xml", "XML, no page numbers"), ("md", "full text as markdown, no page numbers")):
        if os.path.exists(os.path.join(ROOT, "pdfs", f"{key}.{ext}")):
            return f"`pdfs/{key}.{ext}` ({label})"
    if os.path.exists(os.path.join(ROOT, "abstracts", f"{key}.txt")):
        return "abstract only"
    return "none"


def short_ref(key):
    c = cr.get(key, {})
    if not c.get("title"):
        return works[key]["cite"]
    j = c.get("journal_short") or c.get("journal")
    vol = c.get("volume", "")
    pg = c.get("page") or c.get("article_number") or ""
    if key in PAGES:
        vol, pg = PAGES[key].split(":")
    j = html.unescape(j)
    title = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", c["title"]))
    return f"{works[key]['cite']}. *{title[:90]}{'…' if len(title) > 90 else ''}*. {j} {vol}{':' + pg if pg else ''}"


def replace_block(path, name, body):
    txt = open(path).read() if os.path.exists(path) else ""
    start, end = f"<!-- AUTO:{name} start -->", f"<!-- AUTO:{name} end -->"
    block = f"{start}\n{body}\n{end}"
    if start in txt:
        txt = re.sub(re.escape(start) + r".*?" + re.escape(end), lambda m: block, txt, flags=re.S)
    else:
        txt += "\n" + block + "\n"
    open(path, "w").write(txt)


# ---- claim inventory -------------------------------------------------------------------------
sections = []
for c in claims:
    sec = c["section"].split(" ¶")[0].split(",")[0].split(";")[0]
    sec = {"M&M Plant material": "Materials and methods", "M&M Common garden": "Materials and methods"}.get(sec, sec)
    if sec.startswith("M&M"):
        sec = "Materials and methods"
    if sec.startswith("Results"):
        sec = "Results"
    if sec.startswith("Discussion"):
        sec = "Discussion"
    if sec.startswith("Introduction"):
        sec = "Introduction"
    if sec.startswith("Conclusions"):
        sec = "Conclusions"
    if sec not in sections:
        sections.append(sec)
    c["_sec"] = sec
lines = []
for sec in sections:
    lines += [f"### {sec}", "", "| ID | where | claim in the manuscript | cited | tier | TODO | source held | verdict |",
              "|---|---|---|---|---|---|---|---|"]
    for c in claims:
        if c["_sec"] != sec:
            continue
        keys = c["keys"].split("; ")
        cites = "; ".join(works[k]["cite"] for k in keys)
        src = "<br>".join(f"{k}: {held(k)}" for k in keys) if len(keys) > 1 else held(keys[0])
        vs = [f"{VD[(c['id'], k)]['verdict']}" + (f" {k}" if len(keys) > 1 else "") for k in keys if (c['id'], k) in VD]
        lines.append(f"| {c['id']} | {c['section']} | {c['claim']} | {cites} | {c['tier']} | {c['todo'] or ''} | {src} | {'<br>'.join(vs)} |")
    lines.append("")
n_full = sum(1 for k in works if held(k).startswith("`"))
ab = [k for k in works if works[k]["tier"] in "AB" and not held(k).startswith("`")]
lines.insert(0, f"*Generated by `scripts/03_make_reports.py`. Works with a full text held: {n_full} of {len(works)}; "
                f"tier-A/B works still without one: {len(ab)}.*\n")
replace_block(os.path.join(ROOT, "citation_audit.md"), "inventory", "\n".join(lines))

# ---- verdict details ---------------------------------------------------------------------------
order = {c["id"]: i for i, c in enumerate(claims)}
det = [f"*Generated by `scripts/03_make_reports.py` from `meta/verdicts.tsv`: {len(verdicts)} claim–work pairs checked "
       f"({sum(v['verdict'] == '✅' for v in verdicts)} ✅, {sum(v['verdict'] == '⚠️' for v in verdicts)} ⚠️, "
       f"{sum(v['verdict'] == '❌' for v in verdicts)} ❌, {sum(v['verdict'] == '❓' for v in verdicts)} ❓ cannot verify). Quotes are copied from the files in `pdfs/` "
       "(… = text left out); page = printed page of the journal, or the section for XML full texts.*", ""]
for v in sorted(verdicts, key=lambda v: (order.get(v["id"], 999), v["key"])):
    c = next(c for c in claims if c["id"] == v["id"])
    det += [f"**{v['id']} · {works[v['key']]['cite']} · {v['verdict']}**  ", f"*Manuscript:* {c['claim']}  ",
            f"*Source ({v['where']}):* “{v['quote']}”  "]
    if v["note"]:
        det.append(f"*Note:* {v['note']}")
    det.append("")
replace_block(os.path.join(ROOT, "citation_audit.md"), "details", "\n".join(det))

# ---- download list ---------------------------------------------------------------------------
def claims_of(key):
    return ", ".join(c["id"] for c in claims if key in c["keys"].split("; "))

free, pay, opt = [], [], []
for k, w in works.items():
    if w["tier"] not in "AB" or held(k).startswith("`"):
        continue
    doi = w["doi"]
    row = f"| `{k}.pdf` | {short_ref(k)} | {'https://doi.org/' + doi if doi else '—'} | {claims_of(k)} |"
    if k in FREE:
        free.append(row + f" {FREE[k]} |")
    elif k in OPTIONAL:
        opt.append(row)
    else:
        pay.append(row)
hdr = "| save as | reference | DOI | supports |"
out = [f"*Generated by `scripts/03_make_reports.py`; rows disappear as the files arrive in `pdfs/`.*", "",
       f"### 1. Free: open in your browser, no login needed ({len(free)})", "",
       "The publisher or PubMed Central blocks scripted downloads, but the article is free to read.", "",
       hdr + " free copy |", "|---|---|---|---|---|", *free, "",
       f"### 2. Paywalled: through your institution ({len(pay)})", "", hdr, "|---|---|---|---|", *pay, "",
       f"### 3. Optional: background (tier B), abstract already saved ({len(opt)})", "",
       "Download only if convenient; the abstract may be enough for these claims.", "", hdr, "|---|---|---|---|", *opt]
replace_block(os.path.join(ROOT, "to_download.md"), "list", "\n".join(out))
print("reports written:", len(free), "free,", len(pay), "paywalled,", len(opt), "optional")

# ---- bibliographic check table ---------------------------------------------------------------
# works without a Crossref record, checked by hand (2026-10-01; source of each check given)
MANUAL = {
    "AACC2001": "AACC (2001) The definition of dietary fiber. Report of the Dietary Fiber Definition Committee to the Board of Directors of the American Association of Cereal Chemists. Cereal Foods World 46:112-126. https://www.cerealsgrains.org/resources/definitions/Documents/DietaryFiber/DFDef.pdf | no DOI; free PDF held (pdfs/AACC2001.pdf); added in REF2",
    "BroadInstitute2019": "Picard toolkit (2019). Broad Institute, GitHub repository. https://broadinstitute.github.io/picard/ | the citation the Picard README asks for",
    "FAO2025": "FAOSTAT statistical database. FAO, Rome. https://www.fao.org/faostat/ (access date needed) | online database; the 4th-cereal statement to be checked in FAOSTAT (TODO 101)",
    "Holm1979": "A simple sequentially rejective multiple test procedure. Scand J Stat 6:65-70 | no DOI (JSTOR 4615733; the 10.2307 DOI does not resolve)",
    "Pedersen2025": "patchwork: the composer of plots. R package version 1.3.2. https://doi.org/10.32614/CRAN.package.patchwork | CRAN: v1.3.2 published 2025-08-25",
    "RCoreTeam2021": "R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna. https://www.R-project.org/ | R 4.1.2 was released 2021-11-01",
    "Tomczak2014": "The need to report effect size estimates revisited. An overview of some recommended measures of effect size. Trends Sport Sci 21(1):19-25 | no DOI; PDF held (Wielkopolska Digital Library)",
    "VanderAuwera2020": "Genomics in the cloud: using Docker, GATK, and WDL in Terra. O'Reilly Media, Sebastopol. ISBN 978-1-491-97519-0 | Open Library record",
}
NOTES = {
    "Keeling1993": "pages 342-348 from the publisher listing (Crossref has no pages); added in REF2",
    "Hrmova2024": "4 authors (Hrmova, Zimmer, Bulone, Fincher): was *Hrmova and Fincher 2024*, fixed in REF1",
    "Yates2026": "title *Ensembl 2026*; online 2025-11-26, issue 54(D1) January 2026: was *Yates et al. 2025*, fixed in REF1",
    "PotapenkoMolEcol2026": "Mol Ecol 35(17), September 2026 (online 2026-09-01): 2026b (was 2026a; letters swapped in REF1)",
    "PotapenkoMBE2026": "Mol Biol Evol 43(3), 10 March 2026: 2026a (was 2026b; letters swapped in REF1)",
    "Holland2003": "Crossref dates the volume 2002 (published 30 Oct 2002); usually cited 2003; no PDF obtainable, kept 2003 (Word comment 257)",
    "Petit2016": "Crossref holds the early-online record; final pages 171(2):894-913 (Europe PMC)",
    "Huang2023": "PNAS 120(14):e2205783119 (the article number ends in 119, not 120)",
    "Lim2020": "online 2019, issue 71(1) January 2020: 2020 correct",
    "Farag2022": "online 2020, issue 62(4) 2022: 2022 correct",
    "Hubner2013": "online 2012, issue 26(1) 2013: 2013 correct",
    "Guo2008": "online 2007, issue 177(4) 2008: 2008 correct",
    "Konig2023": "Crossref print date 2022-12-28, online 2023-04-21, volume 12 (2023): 2023 correct",
    "GBD2020": "corporate author; 940 named authors",
    "Xu2018": "candidate for the missing β-amylase citation (TODO 411), not yet in the text",
    "UniProt2025": "online 2024, issue 53(D1) January 2025: 2025 correct",
}
rows = ["| key | in the text | record found (first author, title, journal volume:pages) | DOI | note |", "|---|---|---|---|---|"]
for k, w in works.items():
    c = cr.get(k, {})
    if k in MANUAL:
        rec, note = MANUAL[k].split(" | ")
        rows.append(f"| {k} | {w['cite']} | {rec} | — | {note} |")
        continue
    j = html.unescape(c.get("journal_short") or c.get("journal") or "")
    vol, pg = c.get("volume", ""), c.get("page") or c.get("article_number") or ""
    if k in PAGES:
        vol, pg = PAGES[k].split(":")
    title = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", c.get("title", "")))
    fa = c.get("first_author", "").title() if c.get("first_author", "").isupper() else c.get("first_author", "")
    n = int(c.get("n_authors") or 0)
    second = c.get("authors", "").split("; ")[1].split(" ")[0] if n == 2 else ""
    second = second.title() if second.isupper() else second
    rec = f"{fa}{' et al' if n > 2 else (' and ' + second if n == 2 else '')} — {title} — {j} {vol}{':' + pg if pg else ''}"
    rows.append(f"| {k} | {w['cite']} | {rec} | https://doi.org/{c.get('doi', '')} | {NOTES.get(k, '✅')} |")
replace_block(os.path.join(ROOT, "bibliography_check.md"), "table", "\n".join(rows))
print("bibliography table written:", len(rows) - 2, "works")
