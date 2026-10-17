#!/usr/bin/env python3
"""esm_common.py -- shared pieces of the Online Resource (ESM) builders: the TAG title block, the "About" sheet of
every .xlsx, a table-sheet writer, and the guard against the relevance-screen categories.

TAG (journal_TAG/TAG_submission_guidelines_2026-09-15.pdf, "Supplementary Information", pp. 37-39): every ESM file
carries the article title, the journal name, the author names and the corresponding author's affiliation and e-mail;
the manuscript carries a concise caption for each file; files are published as received.

Layout of every ESM_N.xlsx (user decision, ESM0 Q8/Q9, 2026-10-01): sheet 1 "About" = title block, "Online Resource N",
the caption, the list of sheets and the column definitions; the data sheets follow, each starting with its header row
(machine-readable). The project paths of the sources are NOT written into the published files; they are in the
builders' docstrings and in README.md.

AUTHOR NAMES, AFFILIATION AND E-MAIL ARE PLACEHOLDERS until the title page is settled (manuscript checklist item 1,
item CHK1). Set them below and rebuild all files (bash build_all.sh).

Created 2026-10-01 (supplementary-material session; build/proposed/ESM0_QUESTIONS.md).
"""
import pathlib, re
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont

ROOT = pathlib.Path("/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project")
NP = ROOT / "10_USED_Paper_writing/new_publishing_paper"
SUPP = NP / "supplementary"
WORKING_FILE = NP / "Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md"

# ---- TAG title block -------------------------------------------------------------------------------------------
TITLE_PLAIN = ("Genome-wide association and haplotype analysis identify candidate genes for grain nutritional quality "
               "in wild barley (Hordeum spontaneum)")
TITLE_PARTS = [("Genome-wide association and haplotype analysis identify candidate genes for grain nutritional "
                "quality in wild barley (", False), ("Hordeum spontaneum", True), (")", False)]
JOURNAL = "Theoretical and Applied Genetics"
AUTHORS = "[AUTHOR NAMES]"                                   # placeholder (checklist item 1)
CORR_NAME = "[CORRESPONDING AUTHOR]"                         # placeholder
CORR_AFFILIATION = "[AFFILIATION: institution, department, city, country]"   # placeholder
CORR_EMAIL = "[E-MAIL]"                                      # placeholder
CORR_LINE = f"Corresponding author: {CORR_NAME}, {CORR_AFFILIATION}; e-mail: {CORR_EMAIL}"


def check_title():
    """the title in the files must be the manuscript's title, word for word"""
    t = WORKING_FILE.read_text(encoding="utf-8")
    md_title = "**Genome-wide association and haplotype analysis identify candidate genes for grain nutritional quality in wild barley (*Hordeum spontaneum*)**"
    assert md_title in t, "the manuscript title changed: update TITLE_PARTS in esm_common.py"
    assert "".join(p for p, _ in TITLE_PARTS) == TITLE_PLAIN


# ---- the relevance-screen categories must not appear in any supplementary file (user, 2026-09-30) ------------------
BANNED_CELLS = {"strong", "plausible", "unlikely", "no_annotation", "no annotation", "trait_candidate_strength",
                "candidate_rationale", "tier"}
BANNED_WORDS = ("candidacy", "candidate strength", "rationale", "relevance screen", "relevance-screen")


def check_no_screen_categories(wb):
    for ws in wb.worksheets:
        for row in ws.iter_rows(values_only=True):
            for v in row:
                if v is None:
                    continue
                s = str(v).strip().lower()
                assert s not in BANNED_CELLS, (ws.title, v)
                assert not any(w in s for w in BANNED_WORDS), (ws.title, v)


# ---- styles ------------------------------------------------------------------------------------------------------
BOLD = Font(bold=True)
THIN = Side(style="thin")
WRAP_TOP = Alignment(wrap_text=True, vertical="top")


def about_sheet(wb, n, caption, sheets, definitions, notes=()):
    """sheet 1: title block, 'Online Resource n', caption, sheet list, column definitions, notes.
    sheets: [(sheet name, one-line description)]; definitions: {sheet name: [(column, definition)]}"""
    check_title()
    ws = wb.active
    ws.title = "About"
    title = CellRichText([TextBlock(InlineFont(b=True, i=it), p) for p, it in TITLE_PARTS])
    ws.append([title]); ws.append([JOURNAL]); ws.append([AUTHORS]); ws.append([CORR_LINE]); ws.append([])
    ws.append([f"Online Resource {n}" if isinstance(n, int) else str(n)]); ws.cell(ws.max_row, 1).font = Font(bold=True, size=12)
    ws.append([caption]); ws.cell(ws.max_row, 1).alignment = WRAP_TOP
    ws.row_dimensions[ws.max_row].height = max(30, 15 * (len(caption) // 120 + 1))
    ws.append([])
    ws.append(["Sheets"]); ws.cell(ws.max_row, 1).font = BOLD
    for name, desc in sheets:
        ws.append([name, desc]); ws.cell(ws.max_row, 2).alignment = WRAP_TOP
    for name, defs in definitions.items():
        ws.append([]); ws.append([f"Columns of sheet “{name}”"]); ws.cell(ws.max_row, 1).font = BOLD
        for col, d in defs:
            ws.append([col, d])
            ws.cell(ws.max_row, 1).alignment = WRAP_TOP
            ws.cell(ws.max_row, 2).alignment = WRAP_TOP
    if notes:
        ws.append([]); ws.append(["Notes"]); ws.cell(ws.max_row, 1).font = BOLD
        for k in notes:
            ws.append(["", k]); ws.cell(ws.max_row, 2).alignment = WRAP_TOP
    for r in range(1, 8):          # the title block and caption span both columns
        if ws.cell(r, 1).value is not None:
            ws.merge_cells(start_row=r, start_column=1, end_row=r, end_column=2)
            ws.cell(r, 1).alignment = WRAP_TOP
    ws.row_dimensions[1].height = 32
    ws.column_dimensions["A"].width = 34
    ws.column_dimensions["B"].width = 120
    return ws


def table_sheet(wb, name, header, rows, fmt=None, widths=None, freeze="B2", default_width=12):
    """a data sheet: header row (bold, wrapped), then the rows; fmt = {column: number format}"""
    ws = wb.create_sheet(name)
    ws.append(header)
    for c in ws[1]:
        c.font = BOLD
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        c.border = Border(bottom=THIN)
    for r in rows:
        ws.append([None if (isinstance(v, float) and v != v) else v for v in r])   # NaN -> empty cell
    fmt, widths = fmt or {}, widths or {}
    for j, h in enumerate(header, 1):
        if h in fmt:
            for i in range(2, len(rows) + 2):
                ws.cell(i, j).number_format = fmt[h]
        ws.column_dimensions[get_column_letter(j)].width = widths.get(h, default_width)
    ws.row_dimensions[1].height = 45
    if freeze:
        ws.freeze_panes = freeze
    return ws


def new_workbook():
    return Workbook()


def save(wb, n):
    check_no_screen_categories(wb)
    out = SUPP / f"ESM_{n}.xlsx"
    wb.save(out)
    return out


TRAIT_LABEL = {"betaglucan": "β-glucan", "fiber": "Fiber", "starch": "Starch", "protein": "Protein"}


def trait_label(t):
    return TRAIT_LABEL.get(str(t), str(t))


def strip_comment_header(path):
    """read a TSV whose first lines are '#' comments (e.g. epsilon_choice_supplementary.tsv)"""
    lines = [l for l in pathlib.Path(path).read_text(encoding="utf-8").splitlines() if not l.startswith("#")]
    return lines
