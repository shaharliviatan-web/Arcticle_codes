#!/usr/bin/env python3
r"""build_ESM_pdf.py -- ESM_2.pdf, the Supplementary Figures S1-S5 of the TAG manuscript (one file; numbered by first
citation; S. Hübner's comment, user decision 2026-10-17). Until 2026-10-17 it built four separate figure Online Resources
(ESM_4, ESM_9, ESM_12, ESM_13); the figures, their sources and their legends are unchanged, only the numbering and the
assembly into one file.

Page 1 = TAG title block (article title, journal, authors, corresponding author's affiliation and e-mail),
"Supplementary Figures" and the list of figures; then each figure on its own page(s): header "Fig. Sk" (and the panel
letter), the figure, a short identification line for multi-panel figures, and on the figure's first page its full legend.
THE FIGURES ARE USED AS PRODUCED by their analysis steps (vector PDF pages are embedded unchanged; the two PNG figures of
step 11 at their native 600 dpi, size 174 mm); nothing is redrawn. The map key -> Fig. S number is SUPP_FIGURES in
esm_common.py; figs() below holds the legends and panels and is also read by package_for_review.py (Word copy).

  Fig. S1 (key "ld", former ESM 4)  02_USED_LD_decay_V2_wholegenome/results/figure_ld_decay_genomewide.pdf (a),
          figure_ld_decay_perchromosome.pdf (b) (both from the same PLINK run as the GWAS SNP set: 7,110,996 SNPs,
          290 accessions; checked 2026-10-01)
  Fig. S2 (key "crosshap", former ESM 9)  page 5 (the full crosshap plot; page 6 is the violin) of the combined
          crosshap PDFs at the settings used: 04_.../loci_LDspan_eps06_V4/CombinedPDF/<trait>/<gene>/
          CrosshapTree+Violin__<gene>__MGmin2__Eps0p6.pdf (GPAT6, GH17 fiber; PHT4;3 starch) and 03_01_7H_.../results/
          figures/crosshap/CrosshapTree+Violin__7HG0729030__{fiber,starch}__MGmin3_eps0.9.pdf (GDSL: both traits, user
          ESM0 Q10a, 2026-10-01; the haplotypes are the same, the phenotype panels differ)
  Figs. S3, S4 (keys "origin", "enrichment", former ESM 12)  11_USED_.../results/figures/Fig_carrier_origin.png,
          Fig_desert_enrichment.png (600 dpi; legends from results/CAPTIONS.md, step 11's draft "Figs. S1-S2")
  Fig. S5 (key "sites", former ESM 13)  08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf (former Fig. 2a)

The page helpers title_page / figure_page / compile_pages / finish are also used by build_ESM_extra.py (the unnumbered
supplement), unchanged.

Method: each page is a small LaTeX document (pdflatex, Helvetica = URW Nimbus Sans, Arial metrics; Greek letters
and math symbols mapped with \DeclareUnicodeCharacter; lualatex cannot be used: luaotfload is not installed), the pages are joined with
pdfunite, and a final Ghostscript pass embeds every font (the crosshap PDFs reference non-embedded Helvetica/Courier),
keeps images lossless (Flate, no downsampling) and writes the document title. Intermediate files go to a temporary
folder under /mnt/data/shahar/.tmp; the .tex sources are kept in scripts/tex/ for review.

Run: python3 build_ESM_pdf.py   (pdflatex, pdfunite, gs). Created 2026-10-01; one file of Supplementary Figures since 2026-10-17.
"""
import os, re, shutil, subprocess, sys, tempfile, pathlib
from esm_common import (ROOT, SUPP, TITLE_PARTS, TITLE_PLAIN, JOURNAL, AUTHORS, CORR_LINE, check_title, FIGURES_FILE,
                        SUPP_FIGURES, N_FIGURES, N_TABLES, check_supp_refs, fig_no)

POPPLER = "/mnt/data/shahar/gwas_barley/tools/miniconda3/bin"
TEXDIR = pathlib.Path(__file__).resolve().parent / "tex"
TEXDIR.mkdir(exist_ok=True)
os.environ["TMPDIR"] = "/mnt/data/shahar/.tmp"
# luaotfload's font cache defaults to ~/.texlive*/texmf-var: keep it off the home directory (user rule)
os.environ["TEXMFVAR"] = os.environ["TEXMFCACHE"] = "/mnt/data/shahar/.tmp/texmf-var"

LD = ROOT / "02_USED_LD_decay_V2_wholegenome/results"
CH = ROOT / "04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/CombinedPDF"
CH7 = ROOT / "03_01_7H_branch_Starch_Fiber_shared_signal_explore/results/figures/crosshap"
OR = ROOT / "11_USED_check_low_maf_geographic_distribution_minor_allele/results/figures"
SITE = ROOT / "08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf"


def tex(s):
    """escape LaTeX specials, then *italic* and **bold** markup"""
    for a, b in [("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"), ("_", r"\_"),
                 ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")]:
        s = s.replace(a, b)
    s = s.replace(r"\textbackslash\{\}", r"\textbackslash{}").replace(r"\textasciitilde\{\}", r"\textasciitilde{}") \
         .replace(r"\textasciicircum\{\}", r"\textasciicircum{}")
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    return re.sub(r"\*(.+?)\*", r"\\textit{\1}", s)


PRE = r"""\documentclass[10pt]{article}
\usepackage[T1]{fontenc}\usepackage[utf8]{inputenc}\usepackage{textcomp}
\usepackage[scaled]{helvet}\renewcommand{\familydefault}{\sfdefault}
\DeclareUnicodeCharacter{03B5}{\ensuremath{\varepsilon}}\DeclareUnicodeCharacter{03B2}{\ensuremath{\beta}}
\DeclareUnicodeCharacter{03B7}{\ensuremath{\eta}}\DeclareUnicodeCharacter{03B1}{\ensuremath{\alpha}}
\DeclareUnicodeCharacter{03C7}{\ensuremath{\chi}}\DeclareUnicodeCharacter{2212}{\ensuremath{-}}
\DeclareUnicodeCharacter{2265}{\ensuremath{\geq}}\DeclareUnicodeCharacter{2264}{\ensuremath{\leq}}
\DeclareUnicodeCharacter{2192}{\ensuremath{\rightarrow}}
\usepackage[%s]{geometry}
\usepackage{graphicx}
\pagestyle{empty}\setlength{\parindent}{0pt}\setlength{\parskip}{3pt}
\begin{document}
"""


def head(n):
    """page heading: "Online Resource N" for a numbered file, else the given text (unnumbered supplement)"""
    return f"Online Resource {n}" if isinstance(n, int) else str(n)


def title_page(n, legend, pages):
    title = "".join(r"\textit{%s}" % tex(p) if it else tex(p) for p, it in TITLE_PARTS)
    body = (r"{\large\bfseries %s\par}\vspace{3mm}" % title + "\n" + tex(JOURNAL) + r"\par" + "\n" + tex(AUTHORS)
            + r"\par" + "\n" + tex(CORR_LINE) + r"\par\vspace{10mm}" + "\n"
            + r"{\large\bfseries %s\par}\vspace{3mm}" % tex(head(n)) + "\n"
            + "\n\n".join(tex(p) for p in legend) + r"\par\vspace{6mm}" + "\n"
            + r"{\small " + r"\\ ".join(tex(p) for p in pages) + "}\n")
    return PRE % "a4paper,margin=20mm" + body + r"\end{document}" + "\n"


def figure_page(n, label, graphic, ident, paper, fig_opts):
    """paper = 'a4paper' or 'a4paper,landscape'; graphic = (path, page or None)"""
    path, page = graphic
    opt = (f"page={page}," if page else "") + fig_opts
    body = (r"{\small %s\hfill %s}\par\vspace{2mm}" % (tex(head(n)), tex(label)) + "\n"
            + r"\begin{center}\includegraphics[%s]{%s}\end{center}" % (opt, str(path)) + "\n"
            + r"\vspace{-2mm}{\small %s\par}" % tex(ident) + "\n")
    return PRE % (paper + ",margin=10mm") + body + r"\end{document}" + "\n"


def compile_pages(n, docs, work):
    pdfs = []
    for i, src in enumerate(docs, 1):
        name = f"ESM_{n}_p{i}"
        (TEXDIR / f"{name}.tex").write_text(src, encoding="utf-8")
        shutil.copy(TEXDIR / f"{name}.tex", work / f"{name}.tex")
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{name}.tex"], cwd=work,
                           capture_output=True, text=True, errors="replace")
        log = (work / f"{name}.log").read_text(encoding="utf-8", errors="replace")
        assert r.returncode == 0, f"pdflatex failed for {name}:\n{r.stdout[-2000:]}"
        assert "Missing character" not in log, f"{name}: a character is missing from the font"
        pages = int(re.search(r"Pages:\s+(\d+)", subprocess.run([f"{POPPLER}/pdfinfo", str(work / f"{name}.pdf")],
                    capture_output=True, text=True).stdout).group(1))
        assert pages == 1, f"{name}: {pages} pages (content overflows one page)"
        pdfs.append(work / f"{name}.pdf")
    return pdfs


def finish(n, pdfs, work, doc_title, out=None):
    joined = work / f"ESM_{n}_joined.pdf"
    subprocess.run([f"{POPPLER}/pdfunite", *map(str, pdfs), str(joined)], check=True)
    marks = work / "docinfo.ps"
    esc = lambda s: s.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
    # PostScript strings are Latin-1: keep the document title ASCII
    marks.write_text("[ /Title (%s) /Author (%s) /Subject (%s) /DOCINFO pdfmark\n" %
                     (esc(doc_title.encode("ascii", "replace").decode()), esc(AUTHORS), esc(JOURNAL)), encoding="latin-1")
    out = out or SUPP / f"ESM_{n}.pdf"
    subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-sDEVICE=pdfwrite", "-dCompatibilityLevel=1.6",
                    "-dEmbedAllFonts=true", "-dSubsetFonts=true", "-dAutoRotatePages=/None",
                    "-dDownsampleColorImages=false", "-dDownsampleGrayImages=false", "-dDownsampleMonoImages=false",
                    "-dAutoFilterColorImages=false", "-dColorImageFilter=/FlateEncode",
                    "-dAutoFilterGrayImages=false", "-dGrayImageFilter=/FlateEncode",
                    f"-sOutputFile={out}", str(joined), str(marks)], check=True)
    fonts = subprocess.run([f"{POPPLER}/pdffonts", str(out)], capture_output=True, text=True).stdout.splitlines()[2:]
    assert fonts and all(l.split()[-5] == "yes" for l in fonts), f"{out.name}: a font is not embedded:\n" + "\n".join(fonts)
    info = subprocess.run([f"{POPPLER}/pdfinfo", str(out)], capture_output=True, text=True).stdout
    print(f"wrote {out}  ({re.search(r'Pages:\s+(\d+)', info).group(1)} pages, all fonts embedded)")
    return out


# ---- the Supplementary Figures (one file, ESM_2.pdf; 2026-10-17) ----------------------------------------------
def _cap(n):
    cap = (ROOT / "11_USED_check_low_maf_geographic_distribution_minor_allele/results/CAPTIONS.md").read_text(encoding="utf-8")
    text = re.search(r"\*\*Fig\. S%d\*\* (.+?)\n" % n, cap).group(1)       # step 11's own draft numbering (S1, S2)
    # its cross-references use that draft numbering too ("labels as in Fig. S1"): map them to the manuscript's numbers
    draft = {"1": fig_no("origin"), "2": fig_no("enrichment")}
    return re.sub(r"Fig\. S([12])\b", lambda m: f"Fig. S{draft[m.group(1)]}", text)


def _crosshap_genes():
    return [("a", "fiber", "HORVU.MOREX.r3.3HG0301300", "GPAT6", CH / "fiber/HORVU.MOREX.r3.3HG0301300/CrosshapTree+Violin__HORVU.MOREX.r3.3HG0301300__MGmin2__Eps0p6.pdf", "ε = 0.6, MGmin = 2"),
            ("b", "fiber", "HORVU.MOREX.r3.5HG0487060", "GH17", CH / "fiber/HORVU.MOREX.r3.5HG0487060/CrosshapTree+Violin__HORVU.MOREX.r3.5HG0487060__MGmin2__Eps0p6.pdf", "ε = 0.6, MGmin = 2"),
            ("c", "starch", "HORVU.MOREX.r3.3HG0301710", "PHT4;3", CH / "starch/HORVU.MOREX.r3.3HG0301710/CrosshapTree+Violin__HORVU.MOREX.r3.3HG0301710__MGmin2__Eps0p6.pdf", "ε = 0.6, MGmin = 2"),
            ("d", "fiber", "HORVU.MOREX.r3.7HG0729030", "GDSL esterase/lipase", CH7 / "CrosshapTree+Violin__7HG0729030__fiber__MGmin3_eps0.9.pdf", "ε = 0.9, MGmin = 3"),
            ("e", "starch", "HORVU.MOREX.r3.7HG0729030", "GDSL esterase/lipase", CH7 / "CrosshapTree+Violin__7HG0729030__starch__MGmin3_eps0.9.pdf", "ε = 0.9, MGmin = 3")]


def figs():
    """the five figures: key -> legend paragraphs and panels (letter, graphic, paper, LaTeX options, identification line;
    'opts_legend' = options for the first page, which also carries the legend). Read by the Word builder of the review
    package (package_for_review.py) as well, so both show the same figures and legends."""
    s = (LD / "ld_decay_summary_for_results_chapter.txt").read_text()
    assert "r^2 = 0.20  ->  188.0 kb" in s and "Samples used                          : 290" in s, "LD summary changed"
    genes = _crosshap_genes()
    for *_, f, _s in genes:
        info = subprocess.run([f"{POPPLER}/pdfinfo", str(f)], capture_output=True, text=True).stdout
        assert "Pages:           6" in info, f"{f}: expected 6 pages (page 5 = crosshap plot, page 6 = violin)"
    name = lambda g: g[3] if g[3] == "GDSL esterase/lipase" else f"*{g[3]}*"
    tables_behind = " The tables behind the figure, and the details of the tests, are in Tables S18–S22"
    return {
        "ld": dict(legend=[
            "Decay of linkage disequilibrium (LD) with physical distance in the 290 wild barley accessions: **a** across "
            "the genome; **b** per chromosome. The r² of every pair of the 7,110,996 SNPs within 2.5 Mb (pairs with "
            "r² ≥ 0.001) was computed in PLINK v1.90, averaged in 1-kb distance bins (gray points) and smoothed with a "
            "LOESS curve (span 0.10) weighted by the number of pairs per bin (red line). Dashed line, r² = 0.2; dotted "
            "line in **a**, the distance at which the smoothed genome-wide r² falls to 0.2 (about 188 kb)"],
            panels=[dict(letter="a", graphic=(LD / "figure_ld_decay_genomewide.pdf", None), paper="a4paper",
                         opts="width=174mm", opts_legend="width=174mm,height=200mm,keepaspectratio", ident="a  Genome-wide LD decay"),
                    dict(letter="b", graphic=(LD / "figure_ld_decay_perchromosome.pdf", None), paper="a4paper,landscape",
                         opts="width=260mm", ident="b  LD decay per chromosome")]),
        "crosshap": dict(legend=[
            "Haplotypes of the four presented genes, as plotted by crosshap v1.4.0 at the settings used in the analysis "
            "(SNPs within the gene ± 1 kb): **a** *GPAT6* (HORVU.MOREX.r3.3HG0301300), fiber; **b** *GH17* "
            "(HORVU.MOREX.r3.5HG0487060), fiber; **c** *PHT4;3* (HORVU.MOREX.r3.3HG0301710), starch (ε = 0.6, "
            "MGmin = 2); the GDSL esterase/lipase gene HORVU.MOREX.r3.7HG0729030, **d** fiber and **e** starch (ε = 0.9, "
            "MGmin = 3; the same haplotypes for both traits)",
            "In each plot: top, the number of accessions in each haplotype group (A, B, …); center, the marker-group "
            "(MG) alleles that define each haplotype (open circle, reference allele; filled circle, alternate allele); "
            "left, the genotype counts of the SNPs of each marker group (reference, missing, heterozygous, alternate); "
            "right, the absolute difference in the trait between alternate- and reference-allele carriers at each SNP "
            "of the marker group (Pheno association), shaded by alternate-allele frequency; bottom, the genotype BLUPs "
            "of the trait (% of grain) in each haplotype group. Tables: for each marker group, its number of SNPs "
            "(nSNP), phenotypic difference (phenodiff), mean r² (meanR2) and alternate-allele frequency (AltAF); for "
            "each haplotype group, its mean BLUP (Pheno) and number of accessions (nTotal). Accessions not assigned to "
            "a haplotype group are not shown. The haplotype groups are those of Fig. 3; their tests are in Fig. 3 and "
            "Table S10"],
            panels=[dict(letter=g[0], graphic=(g[4], 5), paper="a4paper,landscape",
                         opts="height=168mm,width=\\linewidth,keepaspectratio",
                         opts_legend="height=112mm,width=\\linewidth,keepaspectratio",
                         ident=f"{g[0]}  {name(g)} ({g[2]}), {g[1]}; {g[5]}") for g in genes]),
        "origin": dict(legend=[_cap(1) + "." + tables_behind],
                       panels=[dict(letter="", graphic=(OR / "Fig_carrier_origin.png", None), paper="a4paper",
                                    opts="width=174mm,height=226mm,keepaspectratio",
                                    opts_legend="width=174mm,height=205mm,keepaspectratio", ident="")]),
        "enrichment": dict(legend=[_cap(2) + "." + tables_behind],
                           panels=[dict(letter="", graphic=(OR / "Fig_desert_enrichment.png", None), paper="a4paper",
                                        opts="width=174mm", opts_legend="width=174mm,height=200mm,keepaspectratio", ident="")]),
        "sites": dict(legend=[
            "Genotype BLUPs of the four nutritional traits (centered; % of grain) of the 290 accessions across the 29 "
            "sampling sites, ordered and colored by ecological region (HZ1, HZ2 and HZ3, the North–Coast, "
            "North–Desert and Coast–Desert intermediate zones). Boxes, median and interquartile range; whiskers, "
            "1.5 × the interquartile range; points, accessions (ten per site). Sites are numbered as in Table S1 "
            "(every fifth site labeled; site 4 was not analyzed)"],
            panels=[dict(letter="", graphic=(SITE, None), paper="a4paper", opts="width=174mm",
                         opts_legend="width=174mm,height=200mm,keepaspectratio", ident="")]),
    }


def fig_page(k, pnl, legend=None):
    """one page of Fig. Sk: header (Fig. Sk, panel letter), the graphic, its identification line, and on the figure's
    first page the legend"""
    path, page = pnl["graphic"]
    opt = (f"page={page}," if page else "") + (pnl.get("opts_legend") or pnl["opts"] if legend else pnl["opts"])
    body = (r"{\small\textbf{Fig. S%d}\hfill %s}\par\vspace{2mm}" % (k, tex(pnl["letter"])) + "\n"
            + r"\begin{center}\includegraphics[%s]{%s}\end{center}" % (opt, str(path)) + "\n")
    if pnl["ident"]:
        body += r"\vspace{-2mm}{\small %s\par}" % tex(pnl["ident"]) + "\n"
    if legend:
        paras = [r"\textbf{Fig. S%d} " % k + tex(legend[0])] + [tex(x) for x in legend[1:]]
        body += r"\vspace{1mm}{\small " + "\n\n".join(paras) + r"\par}" + "\n"
    return PRE % (pnl["paper"] + ",margin=10mm") + body + r"\end{document}" + "\n"


def figures_title_page(entries):
    title = "".join(r"\textit{%s}" % tex(p) if it else tex(p) for p, it in TITLE_PARTS)
    body = (r"{\large\bfseries %s\par}\vspace{3mm}" % title + "\n" + tex(JOURNAL) + r"\par" + "\n" + tex(AUTHORS)
            + r"\par" + "\n" + tex(CORR_LINE) + r"\par\vspace{10mm}" + "\n"
            + r"{\large\bfseries Supplementary Figures\par}\vspace{3mm}" + "\n"
            + tex(f"Supplementary Figures S1–S{N_FIGURES} of the manuscript, numbered in the order in which they are first "
                  f"cited, each with its legend. The tables (Tables S1–S{N_TABLES}) are in the accompanying file of "
                  f"Supplementary Tables.") + r"\par\vspace{4mm}" + "\n"
            + r"{\small " + "\n\n".join(r"\textbf{Fig. S%d} %s" % (k, tex(t) + r"\hfill " + tex(pg)) for k, t, pg in entries)
            + r"\par}" + "\n")
    return PRE % "a4paper,margin=20mm" + body + r"\end{document}" + "\n"


def build_figures(work):
    spec = figs()
    assert set(spec) == {k for k, _ in SUPP_FIGURES}, "FIGS and esm_common.SUPP_FIGURES disagree"
    pages, entries, page = [], [], 2
    for k, (key, title) in enumerate(SUPP_FIGURES, 1):
        f = spec[key]
        n = len(f["panels"])
        entries.append((k, title, f"page {page}" if n == 1 else f"pages {page}–{page + n - 1}"))
        pages += [fig_page(k, pnl, f["legend"] if j == 0 else None) for j, pnl in enumerate(f["panels"])]
        page += n
    for old in TEXDIR.glob("ESM_[0-9]*_p*.tex"):             # the former per-Online-Resource page sources
        old.unlink()
    docs = [figures_title_page(entries)] + pages
    out = finish(2, compile_pages(2, docs, work), work, f"{TITLE_PLAIN} - Supplementary Figures", out=FIGURES_FILE)
    text = subprocess.run([f"{POPPLER}/pdftotext", str(out), "-"], capture_output=True, text=True).stdout
    check_supp_refs(text, out.name)
    for k in range(1, N_FIGURES + 1):
        assert f"Fig. S{k}" in text, f"Fig. S{k} missing from {out.name}"
    return out


if __name__ == "__main__":
    check_title()
    work = pathlib.Path(tempfile.mkdtemp(prefix="esm_figs_", dir="/mnt/data/shahar/.tmp"))
    build_figures(work)
    shutil.rmtree(work)
