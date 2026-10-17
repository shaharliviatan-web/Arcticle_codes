#!/usr/bin/env python3
"""package_for_review.py -- the supplement reorganized for S. Hübner's review (his request of 2026-10-04: all the
material in one Excel file, as sheets, and one Word file, each table and figure with its legend, as in E. Potapenko's
papers; the PDFs could not be edited). Writes ../for_review_Hubner/:

  Online_Resources_tables.xlsx                              Online Resources 1-3, 5-8, 10, 11 (the cited tables)
  Online_Resources_figures.docx                             Online Resources 4, 9, 12, 13 (the cited figures)
  Raw_phenotype_checks_uncited_in_manuscript_tables.xlsx    the unnumbered supplement (../unreferenced/), tables
  Raw_phenotype_checks_uncited_in_manuscript_figures.docx   the unnumbered supplement, figures
  Supplementary_for_review.zip                              the four files

Layout modelled on the supplementary data of Potapenko et al. 2026a, Mol Biol Evol
(reference_sources/Potapenko_2026_MBE/msag051_supplementary_data (1).zip): one workbook with one sheet per table, the
bold label and the legend in the first rows, then the table; one document with each figure followed by its legend.
Here: sheet 1 "Contents" (title block, the manuscript's captions, links to the sheets); every sheet opens with the
caption of its Online Resource, the sheet's own description and its column definitions, then the table (autofilter on
the header). Sheet names start with the Online Resource number (OR1 ...), so the numbering is the manuscript's.
The unnumbered supplement carries descriptive sheet and figure names instead of numbers (user, 2026-10-04: no
numbers, and no statement inside the files that they are not cited; "uncited_in_manuscript" is in the file names).

Inputs (nothing is recomputed or retyped):
  tables  - the built files ../ESM_N.xlsx and ../unreferenced/ESM_extra_phenotype_checks.xlsx, copied cell by cell
            (values, number formats, header styles, column widths; checked after copying); captions, sheet
            descriptions and column definitions parsed from their "About" sheets.
  figures - the same source figures, legends and identification lines as build_ESM_pdf.py and build_ESM_extra.py,
            captured by running their page functions with the LaTeX steps replaced. Vector PDF pages are rendered
            at 300 dpi (pdftoppm); the 600-dpi PNG figures are downscaled to 300 dpi.
  captions on the Contents sheet / title page - the "Supplementary Information" section of the manuscript.
The TAG submission package (../for_submission/) is not touched.

Requires python-docx (1.2.0, installed with pip --target in /mnt/data/shahar/gwas_barley/tools/pylib, 2026-10-04).
Run: python3 package_for_review.py   (after build_all.sh, so the ESM files are current). Created 2026-10-04.
"""
import math, os, re, shutil, subprocess, sys, tempfile, zipfile, pathlib
from copy import copy

sys.path.insert(0, "/mnt/data/shahar/gwas_barley/tools/pylib")
os.environ["TMPDIR"] = "/mnt/data/shahar/.tmp"

import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from PIL import Image
from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

import esm_common as C
import build_ESM_pdf as BP
import build_ESM_extra as BX

OUT = C.SUPP / "for_review_Hubner"
POPPLER = BP.POPPLER
XLSX_CITED = "Online_Resources_tables.xlsx"
DOCX_CITED = "Online_Resources_figures.docx"
XLSX_EXTRA = "Raw_phenotype_checks_uncited_in_manuscript_tables.xlsx"
DOCX_EXTRA = "Raw_phenotype_checks_uncited_in_manuscript_figures.docx"
ZIP = "Supplementary_for_review.zip"
TABLES, FIGURES = [1, 2, 3, 5, 6, 7, 8, 10, 11], [4, 9, 12, 13]

# sheet names longer than Excel's 31 characters once "OR11 " is added
SHEET_RENAME = {(11, "Haplotype groups – tests"): "Haplotypes – tests",
                (11, "Haplotype groups – by region"): "Haplotypes – by region"}
# the unnumbered supplement: descriptive names instead of Fig. E1-E7 and of the short sheet names
FIG_NAME = {"E1": "Nutritional traits by season", "E2": "Nutritional traits by block",
            "E3": "Correlations between seasons", "E4": "Correlations between traits",
            "E5": "Correlations between blocks", "E6": "Distributions of raw values",
            "E7": "Distributions of genotype BLUPs"}
EXTRA_SHEET = {"Summary by season": "Raw values by season", "Summary by block": "Raw values by block",
               "Correlations seasons": "Correlations between seasons", "Correlations traits": "Correlations between traits",
               "Correlations blocks": "Correlations between blocks", "BLUP summary": "Genotype BLUPs"}


# ---- text helpers -------------------------------------------------------------------------------------------------
def md_parts(s):
    """'**bold**' and '*italic*' markup -> [(text, bold, italic)]"""
    out = []
    for tok in re.split(r"(\*\*.+?\*\*|\*.+?\*)", s):
        if not tok:
            continue
        if tok.startswith("**"):
            out.append((tok[2:-2], True, False))
        elif tok.startswith("*") and len(tok) > 1:
            out.append((tok[1:-1], False, True))
        else:
            out.append((tok, False, False))
    return out


def plain(s):
    return "".join(t for t, _, _ in md_parts(s))


def rich(parts):
    """[(text, bold, italic)] -> openpyxl rich text"""
    return CellRichText([TextBlock(InlineFont(b=b or None, i=i or None), t) if (b or i) else t for t, b, i in parts])


def unnumber(s):
    """Fig. E1 / Figs. E1, E6 / Figs. E3–E5 -> the figure names of the unnumbered supplement"""
    def rep(m):
        ids = []
        for part in re.split(r", | and ", m.group(1)):
            if "–" in part:
                a, b = (int(x[1:]) for x in part.split("–"))
                ids += [f"E{k}" for k in range(a, b + 1)]
            else:
                ids.append(part)
        q = [f"“{FIG_NAME[i]}”" for i in ids]
        return f"the figure {q[0]}" if len(q) == 1 else "the figures " + ", ".join(q[:-1]) + " and " + q[-1]
    s = re.sub(r"Figs?\. (E\d(?:(?:, |–| and )E\d)*)", rep, s).replace("(the figure", "(figure")
    assert not re.search(r"\bE\d\b", s), s
    return s


# Online Resource 12 holds two figures, each with its own panels a and b. The PDF builder (and the manuscript's
# caption) call them "Fig. S1" and "Fig. S2", a numbering used nowhere else (user, 2026-10-04: not good). In the
# review files they are named by content, under the number of their Online Resource. No effect once the labels are
# removed at the source (build_ESM_pdf.py esm12 and the manuscript caption).
OR12_FIG = {"S1": "carrier origin", "S2": "desert enrichment"}


def no_fig_s(s):
    s = re.sub(r"^\*\*Fig\. (S\d)\*\* ", lambda m: f"**Online Resource 12, {OR12_FIG[m.group(1)]}.** ", s)
    s = re.sub(r" \(Fig\. S\d\)", "", s)
    s = re.sub(r"Fig\. (S\d)", lambda m: f"the {OR12_FIG[m.group(1)].replace(' ', '-')} figure", s)
    assert "Fig. S" not in s, s
    return s


def manuscript_captions():
    t = C.WORKING_FILE.read_text(encoding="utf-8")
    caps = dict((int(n), c.strip()) for n, c in
                re.findall(r"\*\*Online Resource (\d+)\*\* ([^\n]+)", t[t.index("# Supplementary Information"):]))
    assert sorted(caps) == list(range(1, 14)), sorted(caps)
    for n, c in caps.items():
        if "Fig. S" in c:
            caps[n] = no_fig_s(c)
            print(f"NOTE: Online Resource {n}: the manuscript caption still names Fig. S1/S2; shown here as: {plain(caps[n])}")
    return caps


# ---- the About sheet of a built ESM file ----------------------------------------------------------------------------
def parse_about(path):
    wb = openpyxl.load_workbook(path)
    rows = [list(r) for r in wb["About"].iter_rows(values_only=True)]
    label, caption = rows[5][0], rows[6][0]
    i = next(k for k, r in enumerate(rows) if r[0] == "Sheets") + 1
    sheets = []
    while i < len(rows) and rows[i][0] is not None:
        sheets.append((rows[i][0], rows[i][1])); i += 1
    defs, notes, cur = {}, [], None
    for r in rows[i:]:
        a = r[0]
        if a is None and (len(r) < 2 or r[1] is None):
            cur = None
        elif isinstance(a, str) and a.startswith("Columns of sheet “"):
            cur = a[len("Columns of sheet “"):-1]; defs[cur] = []
        elif a == "Notes":
            cur = "__notes__"
        elif cur == "__notes__":
            notes.append(r[1])
        elif cur:
            defs[cur].append((a, r[1]))
    assert [s for s, _ in sheets] == wb.sheetnames[1:], (path, sheets, wb.sheetnames)
    assert set(defs) <= set(wb.sheetnames), (path, defs.keys())
    return wb, label, caption, sheets, defs, notes


# ---- Excel ----------------------------------------------------------------------------------------------------------
def title_rows(ws):
    ws.append([CellRichText([TextBlock(InlineFont(b=True, i=it or None), p) for p, it in C.TITLE_PARTS])])
    ws.append([C.JOURNAL]); ws.append([C.AUTHORS]); ws.append([C.CORR_LINE]); ws.append([])


def legend_row(ws, r, value, span, n_chars, bold=False):
    ws.cell(r, 1).value = value
    ws.cell(r, 1).alignment = Alignment(wrap_text=True, vertical="top")
    if bold:
        ws.cell(r, 1).font = Font(bold=True)
    if span > 1:
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=span)
    return n_chars


def copy_table_sheet(wb_out, src, name, legend, rename=None):
    """legend = [(parts, bold_row)]: rows above the table; the source sheet (header in row 1) follows a blank row"""
    ws = wb_out.create_sheet(name)
    ncol = src.max_column
    wmap = {}
    for dim in src.column_dimensions.values():
        for j in range(dim.min or 0, (dim.max or 0) + 1):
            wmap[j] = dim.width
    widths = [wmap.get(j) or 8.43 for j in range(1, ncol + 1)]
    # legend block: merged over the first columns up to about 125 character widths (more columns if the table is narrow)
    span, tot = 0, 0.0
    while tot < 125 and span < 60:
        span += 1
        tot += widths[span - 1] if span <= ncol else 8.43
    for r, (parts, _) in enumerate(legend, 1):
        if parts is None:
            continue
        text = "".join(t for t, _, _ in parts)
        legend_row(ws, r, rich(parts) if any(b or i for _, b, i in parts) else text, span, len(text))
        lines = sum(max(1, math.ceil(len(p) / (tot * 0.95))) for p in text.split("\n"))
        ws.row_dimensions[r].height = max(15, 14.5 * lines + 2)
    h = len(legend) + 2                                     # one blank row, then the header
    for i, row in enumerate(src.iter_rows(), start=h):
        for c in row:
            d = ws.cell(i, c.column, c.value)
            if c.has_style:
                d.font, d.alignment, d.border, d.number_format = copy(c.font), copy(c.alignment), copy(c.border), c.number_format
    ws.row_dimensions[h].height = src.row_dimensions[1].height or 45
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    for j in range(ncol + 1, span + 1):
        ws.column_dimensions[get_column_letter(j)].width = 8.43
    last = h + src.max_row - 1
    ws.auto_filter.ref = f"A{h}:{get_column_letter(ncol)}{last}"
    # check: the table is copied unchanged
    for r in range(1, src.max_row + 1):
        for j in range(1, ncol + 1):
            assert src.cell(r, j).value == ws.cell(h + r - 1, j).value, (name, r, j)
    return ws


def contents_sheet(wb, heading, intro, rows, col_b="Caption (as in the manuscript)", widths=(18, 95, 36)):
    """rows: [(label, caption parts, [(sheet name, description or None)])]"""
    ws = wb.active
    ws.title = "Contents"
    title_rows(ws)
    ws.append([heading]); ws.cell(ws.max_row, 1).font = Font(bold=True, size=12)
    for p in intro:
        ws.append([p]); r = ws.max_row
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=3)
        ws.cell(r, 1).alignment = Alignment(wrap_text=True, vertical="top")
        ws.row_dimensions[r].height = 15 * max(1, math.ceil(len(p) / (sum(widths) * 0.95))) + 2
    ws.append([])
    ws.append(["", col_b, "Sheet"] if col_b else ["Sheet", "Content"])
    for c in ws[ws.max_row]:
        c.font = Font(bold=True)
    for label, cap, sheets in rows:
        first = True
        for s, desc in sheets or [(None, None)]:
            r = ws.max_row + 1
            if first:
                ws.cell(r, 1, label).font = Font(bold=True)
                ws.cell(r, 2, rich(cap) if any(b or i for _, b, i in cap) else "".join(t for t, _, _ in cap))
                ws.cell(r, 2).alignment = Alignment(wrap_text=True, vertical="top")
                ws.row_dimensions[r].height = 15 * max(1, math.ceil(len(plain("".join(t for t, _, _ in cap))) / (widths[1] * 0.95))) + 2
            if s and desc is None:                        # a sheet of this workbook: link to it
                ws.cell(r, 3, s).hyperlink = f"#'{s}'!A1"
                ws.cell(r, 3).font = Font(color="0563C1", underline="single")
            elif s:
                ws.cell(r, 3, s)
            ws.cell(r, 1).alignment = ws.cell(r, 3).alignment = Alignment(vertical="top", wrap_text=True)
            first = False
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    for r in range(1, 5):
        ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=3)
    ws.row_dimensions[1].height = 32
    ws.cell(1, 1).alignment = Alignment(wrap_text=True, vertical="top")
    return ws


def check_sheet_refs(wb):
    names = set(wb.sheetnames)
    for ws in wb.worksheets:
        for row in ws.iter_rows(values_only=True):
            for v in row:
                for ref in re.findall(r"sheets? “([^”]+)”", str(v or "")):
                    assert ref in names, (ws.title, ref)


def build_cited_xlsx(caps):
    wb = openpyxl.Workbook()
    contents = []
    for n in range(1, 14):
        if n in FIGURES:
            contents.append((f"Online Resource {n}", md_parts(caps[n]), [(f"Word file {DOCX_CITED}", "")]))
            continue
        src, label, caption, sheets, defs, notes = parse_about(C.SUPP / f"ESM_{n}.xlsx")
        assert label == f"Online Resource {n}", (n, label)
        new = {s: f"OR{n} {SHEET_RENAME.get((n, s), s)}" for s, _ in sheets}
        assert all(len(v) <= 31 for v in new.values()), new
        fix = lambda t: re.sub(r"(sheets? )“([^”]+)”", lambda m: f"{m.group(1)}“{new.get(m.group(2), m.group(2))}”", t)
        for k, (s, desc) in enumerate(sheets):
            legend = [([(f"Online Resource {n}. ", True, False), (fix(caption), False, False)], False),
                      ([(f"{SHEET_RENAME.get((n, s), s)}. ", True, False), (fix(desc), False, False)], False)]
            if defs.get(s):
                legend += [([("Columns", True, False)], True)]
                legend += [([(f"{col}: ", True, False), (fix(d), False, False)], False) for col, d in defs[s]]
            if notes and k == 0:
                legend += [([("Notes", True, False)], True)] + [([(fix(x), False, False)], False) for x in notes]
            copy_table_sheet(wb, src[s], new[s], legend)
        contents.append((f"Online Resource {n}", md_parts(caps[n]), [(new[s], None) for s, _ in sheets]))
    contents_sheet(wb, "Online Resources: tables", [
        f"The tables of the Online Resources, one sheet per table, in the numbering of the manuscript (sheet names start "
        f"with the number: OR1 = Online Resource 1). Each sheet opens with the caption of its Online Resource, a "
        f"description of the sheet and the definitions of its columns; the table follows. The figures (Online Resources "
        f"4, 9, 12 and 13) are in the Word file {DOCX_CITED}."], contents)
    check_sheet_refs(wb)
    C.check_no_screen_categories(wb)
    wb.save(OUT / XLSX_CITED)
    print("wrote", OUT / XLSX_CITED, f"({len(wb.sheetnames) - 1} table sheets)")


def build_extra_xlsx():
    src, _label, _caption, sheets, defs, notes = parse_about(C.SUPP / "unreferenced/ESM_extra_phenotype_checks.xlsx")
    assert [s for s, _ in sheets] == list(EXTRA_SHEET), sheets
    wb = openpyxl.Workbook()
    panel = unnumber(BX.PANEL)
    for s, desc in sheets:
        legend = [([(f"{EXTRA_SHEET[s]}. ", True, False), (unnumber(desc), False, False)], False),
                  ([(panel, False, False)], False)]
        if defs.get(s):
            legend += [([("Columns", True, False)], True)]
            legend += [([(f"{col}: ", True, False), (unnumber(d), False, False)], False) for col, d in defs[s]]
        assert not notes
        copy_table_sheet(wb, src[s], EXTRA_SHEET[s], legend)
    contents_sheet(wb, "Raw phenotype checks: tables", [
        "Descriptive checks of the raw phenotypic data of the common-garden experiment (three seasons, 2019–20 to "
        "2021–22), one sheet per table. The figures are in the accompanying Word file.", panel],
        [("", [(unnumber(d), False, False)], [(EXTRA_SHEET[s], None)]) for s, d in sheets],
        col_b="Content", widths=(4, 95, 30))
    check_sheet_refs(wb)
    C.check_no_screen_categories(wb)
    wb.save(OUT / XLSX_EXTRA)
    print("wrote", OUT / XLSX_EXTRA, f"({len(wb.sheetnames) - 1} table sheets)")


# ---- figures: capture what the PDF builders put on their pages -----------------------------------------------------
def capture(module, call):
    rec = {"legend": None, "panels": []}
    saved = {k: getattr(module, k) for k in ("title_page", "figure_page", "compile_pages", "finish")}
    module.title_page = lambda n, legend, pages: rec.update(legend=legend, n=n) or ""
    module.figure_page = lambda n, label, graphic, ident, paper, opts: rec["panels"].append(
        dict(label=label, graphic=graphic, ident=ident, landscape="landscape" in paper)) or ""
    module.compile_pages = lambda *a, **k: []
    module.finish = lambda *a, **k: None
    try:
        call()
    finally:
        for k, v in saved.items():
            setattr(module, k, v)
    return rec


def render(graphic, work, tag):
    """figure -> PNG at about 300 dpi; returns (path, width px, height px)"""
    path, page = graphic
    path = pathlib.Path(path)
    out = work / f"{tag}.png"
    if path.suffix.lower() == ".pdf":
        subprocess.run([f"{POPPLER}/pdftoppm", "-png", "-r", "300", "-singlefile", "-f", str(page or 1),
                        "-l", str(page or 1), str(path), str(out.with_suffix(""))], check=True)
    else:
        im = Image.open(path)
        dpi = (im.info.get("dpi") or (300, 300))[0]
        if dpi > 310:
            im = im.resize((round(im.width * 300 / dpi), round(im.height * 300 / dpi)), Image.LANCZOS)
        im.save(out, dpi=(300, 300))
    w, h = Image.open(out).size
    return out, w, h


# ---- Word ----------------------------------------------------------------------------------------------------------
PORTRAIT = (21.0, 29.7)
MARGIN = 2.0


def new_doc(title):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Arial"; st.font.size = Pt(10)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Arial")
    st.paragraph_format.space_after = Pt(4)
    for name, size in (("Heading 1", 12),):
        h = doc.styles[name]
        h.font.name = "Arial"; h.font.size = Pt(size); h.font.bold = True; h.font.color.rgb = RGBColor(0, 0, 0)
        rf = h.element.rPr.rFonts
        for a in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
            if rf.get(qn(a)) is not None:
                del rf.attrib[qn(a)]
        rf.set(qn("w:ascii"), "Arial"); rf.set(qn("w:hAnsi"), "Arial")
        h.paragraph_format.space_before = Pt(0); h.paragraph_format.space_after = Pt(6)
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(PORTRAIT[0]), Cm(PORTRAIT[1])
    for a in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, a, Cm(MARGIN))
    page_number(sec)
    doc.core_properties.title = title
    doc.core_properties.author = C.AUTHORS
    doc.core_properties.comments = ""
    return doc


def page_number(sec):
    p = sec.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    fld = OxmlElement("w:fldSimple"); fld.set(qn("w:instr"), "PAGE")
    r = OxmlElement("w:r"); t = OxmlElement("w:t"); t.text = "1"; r.append(t); fld.append(r)
    p._p.append(fld)


def add_par(doc, parts, size=None, align=WD_ALIGN_PARAGRAPH.JUSTIFY, space_after=None):
    p = doc.add_paragraph()
    p.alignment = align
    for t, b, i in parts:
        r = p.add_run(t); r.bold = b or None; r.italic = i or None
        if size:
            r.font.size = Pt(size)
    if space_after is not None:
        p.paragraph_format.space_after = Pt(space_after)
    return p


def title_page(doc, heading, intro, items):
    p = add_par(doc, [(t, True, it) for t, it in C.TITLE_PARTS], size=13, align=WD_ALIGN_PARAGRAPH.LEFT, space_after=8)
    for line in (C.JOURNAL, C.AUTHORS, C.CORR_LINE):
        add_par(doc, [(line, False, False)], align=WD_ALIGN_PARAGRAPH.LEFT)
    doc.add_paragraph()
    doc.add_paragraph(heading, style="Heading 1")
    for x in intro:
        add_par(doc, md_parts(x) if isinstance(x, str) else x)
    for parts in items:
        q = add_par(doc, parts, align=WD_ALIGN_PARAGRAPH.LEFT, space_after=3)
        q.paragraph_format.left_indent = Cm(0.5)


def est_lines(text, width_cm, size=10):
    return max(1, math.ceil(len(text) / (width_cm * 4.9 * 10 / size)))


def new_section(doc, want):
    """new page with the given orientation. python-docx puts the closing sectPr of the previous section in a new
    empty paragraph, which can spill onto a blank page: move it into the paragraph before"""
    doc.add_section(WD_SECTION.NEW_PAGE)
    body = doc.element.body
    holder = body.findall(qn("w:p"))[-1]
    prev = holder.getprevious()
    sp = holder.find(qn("w:pPr") + "/" + qn("w:sectPr"))
    if sp is not None and prev is not None and prev.tag == qn("w:p") and not holder.findall(".//" + qn("w:t")):
        prev.get_or_add_pPr().append(sp)
        body.remove(holder)
    sec = doc.sections[-1]
    sec.orientation = want
    w, h = (PORTRAIT[1], PORTRAIT[0]) if want == WD_ORIENT.LANDSCAPE else PORTRAIT
    sec.page_width, sec.page_height = Cm(w), Cm(h)


def add_figures(doc, heading, panels, caption, work, tag):
    """one page per panel (landscape where the PDF builder used landscape); heading on the first page; each panel's
    identification line under it; the caption (list of paragraphs, each a list of parts) after the last panel"""
    for k, pnl in enumerate(panels):
        want = WD_ORIENT.LANDSCAPE if pnl["landscape"] else WD_ORIENT.PORTRAIT
        brk = doc.sections[-1].orientation == want          # same orientation: page break; else a new section
        if not brk:
            new_section(doc, want)
        sec = doc.sections[-1]
        first = None
        tw = sec.page_width.cm - 2 * MARGIN
        th = sec.page_height.cm - 2 * MARGIN - 0.8
        if k == 0 and heading:
            first = doc.add_paragraph(heading, style="Heading 1")
            th -= 1.0
        ident = pnl["ident"] if (len(panels) > 1 or pnl["label"]) else ""
        if ident:
            th -= 0.45 * est_lines(plain(ident), tw, 9 if len(ident) <= 200 else 10) + 0.3
        if k == len(panels) - 1:
            need = sum(0.42 * est_lines(plain("".join(t for t, _, _ in c)), tw) + 0.2 for c in caption)
            if th - need >= 0.8 * th:                    # the legend fits under the figure without shrinking it
                th -= need                                  # by more than 20%; otherwise it flows to the next page
        png, wpx, hpx = render(pnl["graphic"], work, f"{tag}_{k}")
        w = min(tw, th * wpx / hpx)
        p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_after = Pt(2)
        p.add_run().add_picture(str(png), width=Cm(w))
        (first or p).paragraph_format.page_break_before = brk
        if ident and len(ident) > 200:
            add_par(doc, md_parts(ident))
        elif ident:
            add_par(doc, md_parts(ident), size=9, align=WD_ALIGN_PARAGRAPH.LEFT)
    for c in caption:
        add_par(doc, c)


def build_cited_docx(caps, work):
    doc = new_doc(f"{C.TITLE_PLAIN} - Online Resources (figures)")
    title_page(doc, "Online Resources: figures",
               [f"The figures of the Online Resources, in the numbering of the manuscript, each followed by its legend. "
                f"The tables (Online Resources 1–3, 5–8, 10 and 11) are in the Excel file {XLSX_CITED}."],
               [[(f"Online Resource {n} ", True, False)] + md_parts(caps[n]) for n in FIGURES])
    for n in FIGURES:
        rec = capture(BP, lambda: getattr(BP, f"esm{n}")(None))
        assert rec["n"] == n, (n, rec["n"])
        rec["legend"] = [no_fig_s(x) for x in rec["legend"]]
        rec["panels"] = [dict(q, ident=no_fig_s(q["ident"])) for q in rec["panels"]]
        caption = [[(f"Online Resource {n}. ", True, False)] + md_parts(rec["legend"][0])] + [md_parts(x) for x in rec["legend"][1:]]
        add_figures(doc, f"Online Resource {n}", rec["panels"], caption, work, f"or{n}")
    doc.save(OUT / DOCX_CITED)
    print("wrote", OUT / DOCX_CITED)


def build_extra_docx(work):
    rec = capture(BX, BX.build_pdf)
    assert len(rec["panels"]) == len(FIG_NAME)
    doc = new_doc(f"{C.TITLE_PLAIN} - raw phenotype checks (figures)")
    overview = re.sub(r" \(Figs?\. E[^)]*\)", "", rec["legend"][0])
    overview = re.sub(r" Tables: ESM_extra_phenotype_checks\.xlsx$", " The tables are in the accompanying Excel file.", overview)
    title_page(doc, "Raw phenotype checks: figures", [unnumber(overview)],
               [[(name, False, False)] for name in FIG_NAME.values()])
    for (lab, *_), pnl in zip(BX.FIGS, rec["panels"]):
        e = lab.split()[-1]
        cap = pnl["ident"]
        assert cap.startswith(f"**Fig. {e}** "), cap
        body = unnumber(cap[len(f"**Fig. {e}** "):])
        pnl = dict(pnl, ident="", label="")
        add_figures(doc, FIG_NAME[e], [pnl], [[(f"{FIG_NAME[e]}. ", True, False)] + md_parts(body)], work, f"extra_{e}")
    doc.save(OUT / DOCX_EXTRA)
    print("wrote", OUT / DOCX_EXTRA)


# ---- final checks ---------------------------------------------------------------------------------------------------
def docx_text(path):
    d = Document(str(path))
    return "\n".join(p.text for p in d.paragraphs)


def xlsx_text(path):
    wb = openpyxl.load_workbook(path)
    return "\n".join(str(v) for ws in wb.worksheets for row in ws.iter_rows(values_only=True) for v in row if v is not None) \
        + "\n" + "\n".join(wb.sheetnames)


def final_checks(caps):
    for f in (XLSX_EXTRA, DOCX_EXTRA):                      # no numbers, no "not cited" statement in the unnumbered files
        t = (xlsx_text if f.endswith("xlsx") else docx_text)(OUT / f)
        assert not re.search(r"\bE\d\b|ESM_extra|not cited|unnumbered|uncited", t), f
    t = docx_text(OUT / DOCX_CITED)
    assert "Fig. S" not in t and "Fig. S" not in xlsx_text(OUT / XLSX_CITED)
    for n in FIGURES:
        assert f"Online Resource {n}." in t, n
    t = xlsx_text(OUT / XLSX_CITED)
    for n in TABLES:
        assert f"Online Resource {n}." in t, n
    for n in range(1, 14):                                  # the manuscript's captions, word for word, on the Contents sheet
        assert plain(caps[n]) in t, n


if __name__ == "__main__":
    C.check_title()
    caps = manuscript_captions()
    if OUT.exists():
        shutil.rmtree(OUT)                                  # only ever holds files written by this script
    OUT.mkdir()
    work = pathlib.Path(tempfile.mkdtemp(prefix="esm_review_", dir="/mnt/data/shahar/.tmp"))
    build_cited_xlsx(caps)
    build_extra_xlsx()
    build_cited_docx(caps, work)
    build_extra_docx(work)
    shutil.rmtree(work)
    final_checks(caps)
    with zipfile.ZipFile(OUT / ZIP, "w", zipfile.ZIP_DEFLATED) as z:
        for f in (XLSX_CITED, DOCX_CITED, XLSX_EXTRA, DOCX_EXTRA):
            z.write(OUT / f, f)
    print(f"zip {OUT / ZIP} ({round((OUT / ZIP).stat().st_size / 1e6, 1)} MB)")
