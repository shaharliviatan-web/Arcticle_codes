# TAG (Theoretical and Applied Genetics) — manuscript requirements

Journal requirements for the manuscript. Writing rules and the project map are in [`CLAUDE.md`](CLAUDE.md).

**Verified 2026-09-15** against the journal pages saved in `journal_TAG/`:
- `journal_TAG/TAG_submission_guidelines_2026-09-15.pdf` (86 pp): **the source**
- `journal_TAG/TAG_how_to_publish_fees_2026-09-15.pdf`: publishing model and fees

Not in the saved PDF: the body of "Mistakes to avoid during manuscript preparation" (generic) and the
linked Springer Nature "list of mandated data types". Fetch them if needed.

---

## ⚠ Implications for THIS manuscript (action items)

| # | issue | action |
|---|---|---|
| 1 | **LLM use must be documented in Methods.** Only "AI-assisted copy editing" of human-written text is exempt; "generative editorial work and autonomous content creation" is not | user decides the wording; add an LLM-use statement to Materials and methods |
| 2 | Old paper citations use a comma and `&`: `(Godfray et al., 2010)`, `Hübner & Kantar, 2021` | TAG: `(Godfray et al. 2010)`, `(Hübner and Kantar 2021)`, `(Abbott 1991; Barakat et al. 1995a, b)` |
| 3 | Old reference list has full journal names and no DOIs | abbreviate journal names (ISSN LTWA; full title if unsure) and add full DOI links |
| 4 | Old figures use uppercase panel letters (A, B) and "Figure 1." | panels **a, b, c**; captions **`Fig. 1`** in bold; no punctuation after the number or at the end of the caption; no titles inside figure images |
| 5 | **Plant genetic resources:** wild barley accessions = PGR, so genotype/sequencing data deposition in a public repository is **mandatory** | ✅ **Resolved 2026-09-15.** Raw WGS reads for all 300 *H. spontaneum* accessions are public at **ENA PRJEB79623** (ERP163757; 300 NovaSeq WGS runs, public since 2026-03-27). All 290 GWAS accessions are verified present: VCF `HS0103` = ENA alias `01_03`; the 10 excluded are site 04. Published with Potapenko et al. (2026) *Mol Ecol*, https://doi.org/10.1111/mec.70540 (also used in Potapenko et al. 2026 *MBE*, https://doi.org/10.1093/molbev/msag051). Cite PRJEB79623 at the end of M&M and in Data availability. |
| 5b | The published SNP calls are against **MorexV2**; ours are a GATK4 **re-call against MorexV3**, which is not published | Methods must describe the MorexV3 re-call (commands and filters are in the header of `/mnt/data/shahar/gwas_barley/data/inputs/morexV3_with_ids.vcf.gz`). Share the MorexV3 VCF (6.9 GB, too large for ESM) through Zenodo/Dryad, or state "available on request" |
| 6 | QTL/association studies: genotype, phenotype and mapping data are expected as supplementary data | open: BLUPs + locus/SNP tables as Online Resources; MorexV3 genotypes via 5b |
| 7 | Code essential to the results → download links + formatted example data | open: GitHub/Zenodo release of this project's scripts (the lab publishes under https://github.com/hubner-lab) |
| 8 | Minimum study requirements | N = 290 > 100 ✓; 3 growing seasons = 3 environments ✓ ("locations and/or years") |
| 9 | Declarations need real details | funding agency + grant numbers, competing interests, author contributions (initials): ask the supervisor |
| 10 | **Guidelines contradict each other on Declarations position**: "Statements and Declarations … placed **after** the References" vs "Declarations section **before** the reference list" (Competing-interests summary) | follow the primary instruction (after References) unless a recent TAG article shows otherwise |
| 11 | Supplementary files must be cited as "Online Resource N" | name them `ESM_1.xlsx`, `ESM_2.pdf` … in citation order |

---

## Article structure and order

1. **Title page**
   - Concise, informative title
   - Authors; affiliations (institution, department, city, country); corresponding author e-mail
   - ORCID: **mandatory for the corresponding author**
   - **Acknowledgments** go in a separate section on the title page; funding organizations written in full
2. **Author contribution statement**: short, per author, **initials**, all authors. Published in front of the Acknowledgments. Also entered in the submission interface.
3. **Key message**: **≤ 30 words**. States the main achievement *beyond the meaning of the title*. Required for original research.
4. **Abstract**: **150–250 words**. No undefined abbreviations, no unspecified references.
5. **Keywords**: **4–6**.
6. **Introduction**: purpose of the investigation plus a short review of the pertinent literature.
7. **Materials and methods**: enough detail to repeat the work. Sequences essential to repeat it must be disclosed. **Accession numbers go at the end of M&M.** State any material-sharing restrictions here too.
8. **Results**: the outcome, as concisely as possible.
9. **Discussion**: interpretation and significance, with reference to other authors' work.
10. **References**
11. **Statements and Declarations** (see item 10 above on position)

Headings: **no more than 3 levels.**

## Text formatting

- Submit in **Word (.docx)**; editable source files are required at every submission and revision
- Plain font, e.g. **10-pt Times Roman**; italics for emphasis
- **Automatic page numbering**; **no field functions**
- Indents with tab stops, not spaces
- Tables with the Word table function, not spreadsheets
- Equations with the equation editor or MathType
- Abbreviations defined at first mention, then used consistently
- Footnotes, not endnotes. Footnotes are never only a citation and never contain bibliographic details.
- **Genus and species names in italics**
- No word limit for the main text is stated
- No line-numbering or double-spacing requirement is stated

## References

- **In text:** name and year in parentheses, with no comma: `(Thompson 1990)`, `Becker and Seligman (1996)`, `(Abbott 1991; Barakat et al. 1995a, b; Kelso and Smith 1998)`
- **List:** only cited, published or accepted works. Personal communications and unpublished work are mentioned in the text only.
- **Order:** alphabetical by first author's last name.
  - one author: by name, then chronologically
  - two authors: by first author, then coauthor, then chronologically
  - more than two: by first author, then chronologically
- **Always include DOIs as full links** (`https://doi.org/...`)
- **Journal abbreviations** per the ISSN LTWA; use the full title if unsure
- Formats:
  - Journal: `Gamelin FX, Baquet G, Berthoin S (2009) Title. Eur J Appl Physiol 105:731-738. https://doi.org/10.1007/s00421-008-0955-8`
  - Long author lists: `Smith J, Jones M Jr, Houghton L et al (1999) ...`
  - Book: `South J, Blass B (2001) The future of modern genomics. Blackwell, London`
  - Chapter: `Brown B, Aaron M (2001) Title. In: Smith J (ed) Book, 3rd edn. Wiley, New York, pp 230-257`
  - Online: `Author (2007) Title. Publisher. URL. Accessed 26 June 2007`
  - Dataset (DataCite): creator, title, publisher [repository], year, identifier (DOI/accession)

## Statements and Declarations (required; missing ones → returned as incomplete)

| heading | content for this paper |
|---|---|
| Funding | agency + grant numbers, or a "no funding" statement |
| Competing interests | e.g. "The authors have no relevant financial or non-financial interests to disclose." Also entered in the submission interface. |
| Author contributions | free text or CRediT (Conceptualization, Methodology, Formal analysis, Writing…) |
| Data availability | **required**: repository + persistent links, or "available on reasonable request" |
| Ethics approval / consent | human/animal research only. Not applicable here; can be omitted or stated as not applicable |

Material availability: materials (genetic stocks, software…) should be freely available for
non-commercial use within 60 days of request. Restrictions go in the cover letter **and** M&M.

## Tables

- Arabic numerals, cited in text in consecutive order
- Every table has a caption explaining its components
- Previously published material: source reference at the end of the caption
- Footnotes: superscript lower-case letters (or asterisks for significance), placed below the table body

## Figures

| rule | requirement |
|---|---|
| Placement | **inside the body of the text**; separate files only if upload size is a problem |
| Formats | vector → **EPS** (fonts embedded); halftone → **TIFF**; MS Office also accepted |
| File names | `Fig1.eps`, `Fig2.tif` … |
| Resolution | line art ≥ **1200 dpi**; halftone ≥ **300 dpi**; combination (plots with lettering, colour diagrams) ≥ **600 dpi** |
| Colour | **RGB, 8 bits/channel**; free online |
| Size | double-column layout → **84 mm** (1 column) or **174 mm** (full width), height ≤ **234 mm** (TAG assumed to be a large-format journal; small-format would be 119 × 195 mm) |
| Lettering | **Helvetica/Arial**, **8–12 pt** at final size; minimal size variation; no shading or outline effects |
| Lines | ≥ 0.1 mm (0.3 pt) |
| Numbering | Arabic numerals, cited consecutively; **panels a, b, c**; appendix figures continue the main numbering; SI figures numbered separately |
| Captions | in the manuscript text, not in the image. **`Fig. 1`** bold; **no punctuation after the number or at the end**; identify every element |
| Content | **no titles or captions inside the image**; scale bars if magnified |
| Accessibility | descriptive captions; **patterns in addition to colour**; lettering contrast ≥ 4.5:1 |
| Other | state the graphics program used; permission needed for previously published figures |

## Supplementary Information (ESM)

- Cited in text as **"Online Resource 1"** etc.; files named **`ESM_1.xlsx`, `ESM_2.pdf`** consecutively
- Each file contains: article title, journal name, author names, and the corresponding author's affiliation and e-mail
- A concise caption for each file (in the manuscript)
- Formats: spreadsheets **.csv/.xlsx**; text/figure collections **PDF**; multiple files may be zipped
- Published as received, with no editing

## Submission logistics

- Preprints on non-commercial servers (bioRxiv) are **allowed**
- Suggested reviewers: independent, mixed countries and institutions, institutional e-mails
- Revisions: **track changes** (or coloured text) plus an itemized response letter
- Author list cannot change after acceptance
- Hybrid journal: subscription route has **no fee**; open access APC £3590 / $5390 / €4290 (+VAT); OA licences CC BY or CC BY-NC-ND
- Colour in print costs extra; online colour is free
