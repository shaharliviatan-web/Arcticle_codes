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

RENUMBERED 2026-10-17 (S. Hübner's comment; user decision): the 13 "Online Resources" became two running series,
Supplementary Tables S1-S22 and Supplementary Figures S1-S5, each numbered by first citation in the manuscript, one
number per table or figure (as in Potapenko et al. 2026a, Mol Biol Evol). Two files: ESM_1.xlsx (the tables, one sheet
each; build_ESM_tables.py) and ESM_2.pdf (the figures, each with its legend; build_ESM_pdf.py). The content builders
(build_ESM_xlsx.py) still group the tables as before; those workbooks are intermediate (../parts/part_N.xlsx, N = the
former Online Resource number) and are never published. SUPP_TABLES / SUPP_FIGURES below are the one map from the
former grouping to the new numbers.
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
    ws.append([f"Part {n} (intermediate workbook; published as Supplementary Tables in ESM_1.xlsx)" if isinstance(n, int) else str(n)])
    ws.cell(ws.max_row, 1).font = Font(bold=True, size=12)
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
    """an intermediate workbook (one former Online Resource): ../parts/part_N.xlsx, read by build_ESM_tables.py"""
    check_no_screen_categories(wb)
    PARTS.mkdir(exist_ok=True)
    out = PARTS / f"part_{n}.xlsx"
    wb.save(out)
    return out


# ---- the supplementary numbering (2026-10-17) -------------------------------------------------------------------
PARTS = SUPP / "parts"                       # intermediate workbooks, never published
TABLES_FILE = SUPP / "ESM_1.xlsx"            # Supplementary Tables S1-S22
FIGURES_FILE = SUPP / "ESM_2.pdf"            # Supplementary Figures S1-S5
# Table S1, S2, ... in the order of first citation: (part = former Online Resource, sheet of that part, caption).
# Captions are taken from the former Online Resource captions and sheet descriptions (*italic* marks gene names).
SUPP_TABLES = [
    (1, "Sites", "Sampling sites: the 30 sites with their ecological region, coordinates and ten environmental variables"),
    (2, "NIR calibration", "Near-infrared (NIR) calibration: for the 30 calibration accessions (one per sampling site), the "
        "laboratory values, the NIR predictions of the same grain with the manufacturer's calibration for cultivated barley, "
        "and the NIR values of the same accessions in the 2019–20 season"),
    (2, "Calibration correlations", "Validation of the NIR calibration: Pearson correlations among the laboratory values, "
        "the NIR predictions and the 2019–20 NIR values of the 30 calibration accessions (Table S2)"),
    (2, "Blocks 2019–20", "Pearson correlations between the blocks of the 2019–20 season for the four nutritional traits, "
        "from the raw NIR values of the accessions grown in both blocks of each pair; block 4 correlated poorly with the "
        "other blocks for β-glucan and was excluded from the analyses"),
    (1, "Accessions", "The 290 analyzed accessions: identifiers in the European Nucleotide Archive (project PRJEB79623), "
        "site, ecological region and the genotype BLUPs of the ten traits"),
    (3, "Variance components", "Variance components and broad-sense heritability (H²) of the four nutritional and four "
        "morphological and phenological traits, from a linear mixed model with genotype, season-by-block and "
        "genotype-by-season random effects fitted to the values centered within season (290 accessions). For each "
        "component, its variance and its percentage of the total phenotypic variance"),
    (1, "Environment correlations", "Pearson correlations among the ten environmental variables across the 29 sites of "
        "the analyzed accessions"),
    (5, "Loci", "The 36 loci of the genome-wide association analysis of the four nutritional traits (EMMAX; 290 "
        "accessions; MorexV3 coordinates): position, span, lead SNP and its statistics, and the number of candidate genes"),
    (5, "Significant SNPs", "The 52 SNPs that exceeded the genome-wide threshold (Bonferroni, α = 0.10; "
        "−log10(P) = 6.0454), with their locus, allele frequencies and effects"),
    (6, "Candidate genes", "The 55 candidate genes of the 36 loci, one row per gene: locus and position, whether the "
        "gene's haplotypes could be tested and why not, the Kruskal–Wallis test of the trait BLUPs among haplotype groups, "
        "the functional annotation of the 23 genes with significant haplotype associations, and the three genes carried "
        "forward. Empty cells: step not applied to the gene"),
    (7, "Summary", "Choice of the crosshap clustering radius (ε), summary: for five values of ε (MGmin = 2), the number "
        "of the 55 candidate genes whose haplotypes could be tested, the proportion of accessions assigned to haplotype "
        "groups and the number of groups. The analysis used ε = 0.6. Phenotypes were not used"),
    (7, "Per gene", "Choice of the crosshap clustering radius (ε), per gene: the same quantities for each of the 55 "
        "candidate genes at each of the five values of ε (MGmin = 2; Table S11)"),
    (8, "Canonical genes", "The ten canonical (1,3;1,4)-β-glucan synthase and endohydrolase genes located on the MorexV3 "
        "genome from their reference protein sequences, and the distance of each to the nearest β-glucan lead SNP on the "
        "same chromosome"),
    (10, "Site accounting", "Comparison with elite cultivars: per gene, the SNP sites of the wild and elite MorexV3 call "
        "sets in the gene window (gene ± 1 kb), the sites compared and the sites excluded"),
    (10, "Excluded sites", "Comparison with elite cultivars: every SNP site excluded from the comparison, with the reason"),
    (10, "Screen", "Call-rate screen of the 136 spring elite lines of the barley pangenome panel, by pooled call rate "
        "over the compared sites of *GPAT6*, *GH17* and *PHT4;3*"),
    (10, "Cultivars", "The five spring malting cultivars compared with the wild haplotype groups, and their call rates at "
        "the four presented genes"),
    (11, "Lead SNPs – tests", "Geographic origin of the carriers of the lead-SNP alleles of the 36 loci: carriers per "
        "ecological region, and the site-permutation tests of desert enrichment and of the distribution over the six regions"),
    (11, "Lead SNPs – by region", "Minor- and major-allele carriers and missing calls of each lead SNP, per ecological region"),
    (11, "Haplotype groups – tests", "Geographic origin of the haplotype groups of the four presented genes: members per "
        "ecological region, and the site-permutation tests of desert enrichment and of the distribution over the six regions"),
    (11, "Haplotype groups – by region", "Members of each haplotype group of the four presented genes per ecological "
        "region, and the unassigned accessions"),
    (11, "Carriers vs site-mates", "Comparison of the minor-allele carriers of each lead SNP with the major-allele "
        "carriers of their own site"),
]
# Fig. S1, S2, ... in the order of first citation: (key in build_ESM_pdf.FIGS, short title for the contents lists)
SUPP_FIGURES = [
    ("ld", "Decay of linkage disequilibrium with physical distance across the genome (**a**) and per chromosome (**b**)"),
    ("crosshap", "crosshap haplotype plots of *GPAT6*, *GH17*, *PHT4;3* and the GDSL esterase/lipase gene"),
    ("origin", "Geographic origin of the accessions carrying the lead-SNP alleles and the haplotypes of the four "
               "candidate genes"),
    ("enrichment", "Desert enrichment of the lead-SNP alleles and the haplotypes of the four candidate genes"),
    ("sites", "Genotype BLUPs of the four nutritional traits across the 29 sampling sites, by ecological region"),
]
N_TABLES, N_FIGURES = len(SUPP_TABLES), len(SUPP_FIGURES)


def table_no(part, sheet):
    """the Table S number of a sheet of an intermediate workbook"""
    k = [i for i, (p, s, _) in enumerate(SUPP_TABLES, 1) if (p, s) == (part, sheet)]
    assert len(k) == 1, (part, sheet)
    return k[0]


def fig_no(key):
    return [k for k, _ in SUPP_FIGURES].index(key) + 1


def sheet_refs_to_tables(text, part):
    """'sheet “Loci”' in a description of part 5 -> 'Table S8'"""
    return re.sub(r"sheet “([^”]+)”", lambda m: f"Table S{table_no(part, m.group(1))}", text)


def check_supp_refs(text, where):
    """no 'Online Resource' left; every Table S / Fig. S reference exists"""
    assert "Online Resource" not in text, (where, "Online Resource")
    for k in re.findall(r"Tables? S(\d+)", text) + re.findall(r"Tables? S\d+–S(\d+)", text):
        assert 1 <= int(k) <= N_TABLES, (where, f"Table S{k}")
    for k in re.findall(r"Figs?\. S(\d+)", text):
        assert 1 <= int(k) <= N_FIGURES, (where, f"Fig. S{k}")


TRAIT_LABEL = {"betaglucan": "β-glucan", "fiber": "Fiber", "starch": "Starch", "protein": "Protein"}


def trait_label(t):
    return TRAIT_LABEL.get(str(t), str(t))


def strip_comment_header(path):
    """read a TSV whose first lines are '#' comments (e.g. epsilon_choice_supplementary.tsv)"""
    lines = [l for l in pathlib.Path(path).read_text(encoding="utf-8").splitlines() if not l.startswith("#")]
    return lines
