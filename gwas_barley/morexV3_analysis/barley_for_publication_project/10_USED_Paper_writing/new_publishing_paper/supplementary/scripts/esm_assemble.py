#!/usr/bin/env python3
"""esm_assemble.py -- helpers that assemble published workbooks from the intermediate ones: markup parsing
(**bold**, *italic*), reading an "About" sheet, copying a table sheet cell by cell under a legend block, and a Contents
sheet. Moved unchanged from package_for_review.py on 2026-10-17, when build_ESM_tables.py (ESM_1.xlsx, Supplementary
Tables S1-S22) started to use them too.
"""
import math, re
from copy import copy
import openpyxl
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
import esm_common as C


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
