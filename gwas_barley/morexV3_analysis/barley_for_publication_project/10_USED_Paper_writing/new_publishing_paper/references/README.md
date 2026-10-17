# references/: citation audit and bibliography

Started 2026-10-01 (References session). Goal: check that every citation in the manuscript supports its sentence,
then build the TAG reference list in `../05_references.md`. That file stays separate until the user approves it,
and is then merged into the working file. Nothing here changes the manuscript. Fixes go through `../build/proposed/`.

## Workflow

| step | who | status |
|---|---|---|
| 1. Inventory: every claim–citation pair, with a tier (A full text / B abstract / C tool citation) | Claude | ✅ 66 claims, 66 works (`citation_audit.md`); +2 works and 1 claim added in REF2 |
| 2. Bibliographic check: authors, year, journal, volume, pages, DOI (Crossref; by hand for 7) | Claude | ✅ `bibliography_check.md`: findings 1–4 fixed in review item REF1 (2026-10-01); Holland year open (PDF) |
| 3. Access: free full texts downloaded to `pdfs/`; abstracts to `abstracts/` | Claude | ✅ 15 full texts fetched by Claude |
| 4. Download the remaining PDFs into `pdfs/` | user | ✅ 20 full texts supplied (Zotero export in `incoming/`, copied to `pdfs/` under their keys); 6 not obtainable (Friedman 1988, Holland 2003, Savin 1996, Simmonds 1995, Farag 2022, Harlan and Zohary 1966): abstracts only |
| 5. Audit each tier-A/B pair against the local file: verdict ✅/⚠️/❌, with the exact passage and page | Claude | ✅ 2026-10-01: 63 pairs, 45 ✅ 13 ⚠️ 1 ❌ 4 ❓ (incl. Keeling 1993 and AACC 2001, added in REF2); findings 1–13 in `citation_audit.md` |
| 6. Fixes to the manuscript (`build/proposed/` cycle), then `05_references.md` in TAG style | user decides | ✅ REF2 applied; reference list `../05_references.md` built (`scripts/05_build_reference_list.py`) and copied into the working file (REF3), 2026-10-01 |

Rule: verdicts come only from files held here (`pdfs/`, `abstracts/`), never from reading a web page.

## Layout

```
README.md               this file
citation_audit.md       claim inventory + verdicts (the main audit file)
bibliography_check.md   step 2: findings + the verified record of every work
to_download.md          what is still missing (regenerated)
pdfs/<key>.pdf|.xml     full texts (key = FirstAuthorYear, e.g. Burton2008); .xml = Europe PMC full text, no pages
abstracts/<key>.txt     abstracts (Europe PMC / OpenAlex) of works without a full text
meta/works.tsv          master list: key, in-text citation, tier, DOI, search string
meta/claims.tsv         claim–citation pairs (source of the inventory)
meta/crossref/, meta/crossref_summary.tsv, meta/access.tsv   raw outputs of steps 2–3
scripts/01_crossref_check.py        step 2
scripts/02_access_and_download.py   step 3 (re-running overwrites the downloads)
scripts/03_make_reports.py          regenerates the AUTO blocks of the three .md files (verdicts from meta/verdicts.tsv); run after adding PDFs
scripts/04_to_text.py               pdfs/ -> text/<key>.txt (page-marked text used for quoting)
scripts/05_build_reference_list.py  builds ../05_references.md (TAG style) and checks it one to one against the in-text citations
meta/verdicts.tsv                   step 5: one row per claim x work (verdict, page or section, quote, note)
text/                               plain-text versions of the full texts
incoming/                           the user's Zotero export, as received (source of the 20 supplied PDFs)
```

Run the scripts from `references/`. The scripts get no e-mail address and do not get around publisher or PubMed
Central bot checks. Papers blocked that way are listed as "free, open in a browser" for the user.
