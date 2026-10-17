#!/usr/bin/env python3
r"""build_ESM_pdf.py -- builds the figure Online Resources of the TAG manuscript as PDF: ESM_4.pdf (LD decay),
ESM_9.pdf (crosshap plots of the four presented genes), ESM_12.pdf (carrier-origin figures), ESM_13.pdf (site figure).

Each PDF: page 1 = TAG title block (article title, journal, authors, corresponding author's affiliation and e-mail),
"Online Resource N" and the full legend; then one page per figure with a header line, the figure and a short
identification line. THE FIGURES ARE USED AS PRODUCED by their analysis steps (vector PDF pages are embedded unchanged;
the two PNG figures of step 11 at their native 600 dpi, size 174 mm); nothing is redrawn.

  ESM 4   02_USED_LD_decay_V2_wholegenome/results/figure_ld_decay_genomewide.pdf (a), figure_ld_decay_perchromosome.pdf (b)
          (both from the same PLINK run as the GWAS SNP set: 7,110,996 SNPs, 290 accessions; checked 2026-10-01)
  ESM 9   page 5 (the full crosshap plot; page 6 is the violin) of the combined crosshap PDFs at the settings used:
          04_.../loci_LDspan_eps06_V4/CombinedPDF/<trait>/<gene>/CrosshapTree+Violin__<gene>__MGmin2__Eps0p6.pdf
          (GPAT6, GH17 fiber; PHT4;3 starch) and 03_01_7H_.../results/figures/crosshap/
          CrosshapTree+Violin__7HG0729030__{fiber,starch}__MGmin3_eps0.9.pdf (GDSL: both pages, user ESM0 Q10a,
          2026-10-01; the haplotypes are the same, the phenotype panels differ)
  ESM 12  11_USED_.../results/figures/Fig_carrier_origin.png, Fig_desert_enrichment.png (600 dpi; captions from
          results/CAPTIONS.md, Figs. S1-S2)
  ESM 13  08_USED_creating_figures/Figure_ESM_site_BLUPs/ESM_site_BLUPs.pdf (former Fig. 2a)

Method: each page is a small LaTeX document (pdflatex, Helvetica = URW Nimbus Sans, Arial metrics; Greek letters
and math symbols mapped with \DeclareUnicodeCharacter; lualatex cannot be used: luaotfload is not installed), the pages are joined with
pdfunite, and a final Ghostscript pass embeds every font (the crosshap PDFs reference non-embedded Helvetica/Courier),
keeps images lossless (Flate, no downsampling) and writes the document title. Intermediate files go to a temporary
folder under /mnt/data/shahar/.tmp; the .tex sources are kept in scripts/tex/ for review.

Run: python3 build_ESM_pdf.py [4 9 12 13]   (pdflatex, pdfunite, gs). Created 2026-10-01.
"""
import os, re, shutil, subprocess, sys, tempfile, pathlib
from esm_common import ROOT, SUPP, TITLE_PARTS, TITLE_PLAIN, JOURNAL, AUTHORS, CORR_LINE, check_title

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


# ---------------------------------------------------------------------------------------------------------------
def esm4(work):
    s = (LD / "ld_decay_summary_for_results_chapter.txt").read_text()
    assert "r^2 = 0.20  ->  188.0 kb" in s and "Samples used                          : 290" in s, "LD summary changed"
    legend = ["Decay of linkage disequilibrium (LD) with physical distance in the 290 wild barley accessions: **a** across "
              "the genome; **b** per chromosome. The r² of every pair of the 7,110,996 SNPs within 2.5 Mb (pairs with "
              "r² ≥ 0.001) was computed in PLINK v1.90, averaged in 1-kb distance bins (gray points) and smoothed with a "
              "LOESS curve (span 0.10) weighted by the number of pairs per bin (red line). Dashed line, r² = 0.2; dotted "
              "line in **a**, the distance at which the smoothed genome-wide r² falls to 0.2 (about 188 kb)"]
    docs = [title_page(4, legend, ["Page 2: a, genome-wide", "Page 3: b, per chromosome"]),
            figure_page(4, "a", (LD / "figure_ld_decay_genomewide.pdf", None), "a  Genome-wide LD decay",
                        "a4paper", "width=174mm"),
            figure_page(4, "b", (LD / "figure_ld_decay_perchromosome.pdf", None), "b  LD decay per chromosome",
                        "a4paper,landscape", "width=260mm")]
    return finish(4, compile_pages(4, docs, work), work, f"{TITLE_PLAIN} - Online Resource 4")


def esm9(work):
    genes = [("a", "fiber", "HORVU.MOREX.r3.3HG0301300", "GPAT6", CH / "fiber/HORVU.MOREX.r3.3HG0301300/CrosshapTree+Violin__HORVU.MOREX.r3.3HG0301300__MGmin2__Eps0p6.pdf", "ε = 0.6, MGmin = 2"),
             ("b", "fiber", "HORVU.MOREX.r3.5HG0487060", "GH17", CH / "fiber/HORVU.MOREX.r3.5HG0487060/CrosshapTree+Violin__HORVU.MOREX.r3.5HG0487060__MGmin2__Eps0p6.pdf", "ε = 0.6, MGmin = 2"),
             ("c", "starch", "HORVU.MOREX.r3.3HG0301710", "PHT4;3", CH / "starch/HORVU.MOREX.r3.3HG0301710/CrosshapTree+Violin__HORVU.MOREX.r3.3HG0301710__MGmin2__Eps0p6.pdf", "ε = 0.6, MGmin = 2"),
             ("d", "fiber", "HORVU.MOREX.r3.7HG0729030", "GDSL esterase/lipase", CH7 / "CrosshapTree+Violin__7HG0729030__fiber__MGmin3_eps0.9.pdf", "ε = 0.9, MGmin = 3"),
             ("e", "starch", "HORVU.MOREX.r3.7HG0729030", "GDSL esterase/lipase", CH7 / "CrosshapTree+Violin__7HG0729030__starch__MGmin3_eps0.9.pdf", "ε = 0.9, MGmin = 3")]
    for *_, f, _s in genes:
        info = subprocess.run([f"{POPPLER}/pdfinfo", str(f)], capture_output=True, text=True).stdout
        assert "Pages:           6" in info, f"{f}: expected 6 pages (page 5 = crosshap plot, page 6 = violin)"
    legend = ["Haplotypes of the four presented genes, as plotted by crosshap v1.4.0 at the settings used in the analysis "
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
              "a haplotype group are not shown. The haplotype groups are those of Figs. 3 and 4; their tests are in "
              "Online Resource 6 and Fig. 4"]
    pages = [f"Page {i + 2}: {g[0]}, {g[3]} ({g[1]})" for i, g in enumerate(genes)]
    name = lambda g: g[3] if g[3] == "GDSL esterase/lipase" else f"*{g[3]}*"
    docs = [title_page(9, legend, pages)] + [
        figure_page(9, g[0], (g[4], 5), f"{g[0]}  {name(g)} ({g[2]}), {g[1]}; {g[5]}", "a4paper,landscape",
                    "height=168mm,width=\\linewidth,keepaspectratio") for g in genes]
    return finish(9, compile_pages(9, docs, work), work, f"{TITLE_PLAIN} - Online Resource 9")


def esm12(work):
    cap = (ROOT / "11_USED_check_low_maf_geographic_distribution_minor_allele/results/CAPTIONS.md").read_text(encoding="utf-8")
    s1 = re.search(r"\*\*Fig\. S1\*\* (.+?)\n", cap).group(1)
    s2 = re.search(r"\*\*Fig\. S2\*\* (.+?)\n", cap).group(1)
    md2plain = lambda s: re.sub(r"\*\*(.+?)\*\*", r"**\1**", s)
    legend = ["Geographic origin of the carriers of the lead-SNP alleles of the 36 loci and of the haplotype groups of "
              "the four presented genes (Fig. S1), and their enrichment in desert sites (Fig. S2). The tables behind both "
              "figures, and the details of the tests, are in Online Resource 11"]
    docs = [title_page(12, legend, ["Page 2: Fig. S1, carrier origin", "Page 3: Fig. S2, desert enrichment"]),
            figure_page(12, "Fig. S1", (OR / "Fig_carrier_origin.png", None), "**Fig. S1** " + md2plain(s1),
                        "a4paper", "width=174mm,height=226mm,keepaspectratio"),
            figure_page(12, "Fig. S2", (OR / "Fig_desert_enrichment.png", None), "**Fig. S2** " + md2plain(s2),
                        "a4paper", "width=174mm")]
    return finish(12, compile_pages(12, docs, work), work, f"{TITLE_PLAIN} - Online Resource 12")


def esm13(work):
    legend = ["Genotype BLUPs of the four nutritional traits (centered; % of grain) of the 290 accessions across the 29 "
              "sampling sites, ordered and colored by ecological region (HZ1, HZ2 and HZ3, the North–Coast, "
              "North–Desert and Coast–Desert intermediate zones). Boxes, median and interquartile range; whiskers, "
              "1.5 × the interquartile range; points, accessions (ten per site). Sites are numbered as in Online "
              "Resource 1 (every fifth site labeled; site 4 was not analyzed)"]
    docs = [title_page(13, legend, ["Page 2: the figure"]),
            figure_page(13, "", (SITE, None), "Genotype BLUPs by sampling site and ecological region", "a4paper", "width=174mm")]
    return finish(13, compile_pages(13, docs, work), work, f"{TITLE_PLAIN} - Online Resource 13")


if __name__ == "__main__":
    check_title()
    only = set(sys.argv[1:])
    for n, f in [(4, esm4), (9, esm9), (12, esm12), (13, esm13)]:
        if not only or str(n) in only:
            work = pathlib.Path(tempfile.mkdtemp(prefix=f"esm{n}_", dir="/mnt/data/shahar/.tmp"))
            f(work)
            shutil.rmtree(work)
