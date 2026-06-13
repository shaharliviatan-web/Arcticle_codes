#!/usr/bin/env python3
"""
Publication figure assembly for barley GWAS paper.

Converts source PDFs → high-resolution images, assembles multi-panel
figures, adds bold panel labels, and saves both PNG and PDF.

Layout rules
  - All multi-panel figures are horizontal (panels side by side).
  - Maximum 2 panels per row.
  - 4-panel figures use a 2 × 2 grid (2 panels per row).
  - Final width = A4 text-column width (160 mm = 1890 px at 300 DPI),
    suitable for direct insertion into a standard Word A4 document.

Usage:
    python3 assemble_figures.py 1        # Figure 1 only
    python3 assemble_figures.py 2        # Figure 2 only
    python3 assemble_figures.py all      # All figures
"""

import os
import sys
from pathlib import Path
from pdf2image import convert_from_path
from PIL import Image, ImageDraw, ImageFont

# ── Paths ─────────────────────────────────────────────────────────────────────

REPO = Path("/mnt/data/shahar/gwas_barley/morexV3_analysis")
OUT  = REPO / "06_USED_figures"

S1   = REPO / "00_THIN_Generate_Plots_For_Publication/outputs/subsection_1"
C2   = S1   / "C2_env_LMM_r08/figures"
GWAS = REPO / "01_USED_GWAS_V2_pipeline/results/publication_BonfOnly_BLUP_3PC"
LD   = REPO / "02_USED_LD_decay_V2_wholegenome/results"
HAP  = REPO / "04_USED_haplotype_analysis_crosshap/06_publication_figures/results/figures"

FONT_BOLD = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"

# ── Layout constants ──────────────────────────────────────────────────────────

DPI          = 300
# A4 page (210 mm) minus 2 × 25 mm Word margins = 160 mm text width
A4_TEXT_W    = int(160 / 25.4 * DPI)   # 1890 px
PANEL_GAP    = 35                       # px between panels (≈ 3 mm)
ROW_GAP      = 40                       # px between rows (≈ 3.4 mm)
PANEL_LABEL_PX = 62                     # fixed panel-letter size (uniform A/B/C/D)

# ── Helpers ───────────────────────────────────────────────────────────────────

def pdf_to_img(path: Path, border: int = 18) -> Image.Image:
    """Convert PDF page to image with a white border to prevent edge clipping."""
    imgs = convert_from_path(str(path), dpi=DPI, first_page=1, last_page=1)
    img  = imgs[0].convert("RGB")
    canvas = Image.new("RGB", (img.width + 2 * border, img.height + 2 * border), "white")
    canvas.paste(img, (border, border))
    return canvas

def get_font(size: int) -> ImageFont.FreeTypeFont:
    if os.path.exists(FONT_BOLD):
        return ImageFont.truetype(FONT_BOLD, size)
    return ImageFont.load_default()

def resize_to_width(img: Image.Image, width: int) -> Image.Image:
    if img.width == width:
        return img
    return img.resize((width, int(img.height * width / img.width)), Image.LANCZOS)

def resize_to_height(img: Image.Image, height: int) -> Image.Image:
    if img.height == height:
        return img
    return img.resize((int(img.width * height / img.height), height), Image.LANCZOS)

def add_label(img: Image.Image, label: str) -> Image.Image:
    """Add a white strip above the image and stamp the bold panel label there."""
    font_size  = max(48, int(img.width * 0.045))
    font       = get_font(font_size)
    strip_h    = int(font_size * 1.6)   # white strip height above the plot
    pad_left   = int(font_size * 0.25)
    pad_top    = (strip_h - font_size) // 2

    canvas = Image.new("RGB", (img.width, img.height + strip_h), "white")
    canvas.paste(img, (0, strip_h))

    draw = ImageDraw.Draw(canvas)
    draw.text((pad_left, pad_top), label, font=font, fill="black")
    return canvas

def stamp_label(canvas: Image.Image, x: int, y: int, letter: str,
                size: int = PANEL_LABEL_PX) -> None:
    """Draw a fixed-size bold panel letter at (x, y) on an already-composed
    canvas, on a small white box so it reads over any panel content.
    Using a constant size keeps A/B/C/D identical across panels and figures."""
    draw = ImageDraw.Draw(canvas)
    font = get_font(size)
    bbox = draw.textbbox((0, 0), letter, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    pad = 8
    ix, iy = x + 4, y + 4
    draw.rectangle([ix, iy, ix + tw + 2 * pad, iy + th + 2 * pad], fill="white")
    draw.text((ix + pad - bbox[0], iy + pad - bbox[1]), letter, font=font, fill="black")

def hrow(panels: list, total_w: int = A4_TEXT_W, gap: int = PANEL_GAP) -> Image.Image:
    """Place panels side by side, each scaled to an equal share of total_w."""
    n      = len(panels)
    pw     = (total_w - gap * (n - 1)) // n
    panels = [resize_to_width(p, pw) for p in panels]
    h      = max(p.height for p in panels)
    # Compute exact canvas width from actual panel widths (avoids rounding gaps)
    actual_w = sum(p.width for p in panels) + gap * (n - 1)
    canvas   = Image.new("RGB", (actual_w, h), "white")
    x = 0
    for p in panels:
        canvas.paste(p, (x, 0))
        x += p.width + gap
    return canvas

def hrow_fill(panels: list, total_w: int = A4_TEXT_W, gap: int = PANEL_GAP) -> Image.Image:
    """Place panels side by side at a COMMON height so widths fill total_w exactly.
    Unlike hrow (equal widths), this matches heights, so panels of different
    aspect ratios sit flush with no vertical whitespace."""
    ars = [p.width / p.height for p in panels]            # width/height per panel
    h   = int(round((total_w - gap * (len(panels) - 1)) / sum(ars)))
    panels = [resize_to_height(p, h) for p in panels]
    actual_w = sum(p.width for p in panels) + gap * (len(panels) - 1)
    canvas = Image.new("RGB", (actual_w, h), "white")
    x = 0
    for p in panels:
        canvas.paste(p, (x, 0))
        x += p.width + gap
    return canvas

def vstack_rows(rows: list, gap: int = ROW_GAP) -> Image.Image:
    """Stack pre-built rows vertically."""
    w = max(r.width for r in rows)
    h = sum(r.height for r in rows) + gap * (len(rows) - 1)
    canvas = Image.new("RGB", (w, h), "white")
    y = 0
    for r in rows:
        canvas.paste(r, (0, y))
        y += r.height + gap
    return canvas

def save(img: Image.Image, out_dir: Path, name: str):
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{name}.png"
    pdf = out_dir / f"{name}.pdf"
    img.save(str(png), "PNG",  dpi=(DPI, DPI))
    img.save(str(pdf), "PDF",  resolution=DPI)
    print(f"  -> {png}  ({img.width}×{img.height} px)")
    print(f"  -> {pdf}")

# ── Figure builders ───────────────────────────────────────────────────────────

def build_figure_1():
    """
    Figure 1 — Genetic architecture
      A  B1_variance_partitioning_8_traits
      B  B2_reaction_norm_per_trait_highlights
    Layout: A | B  (side by side, equal width). Panel letters stamped at a
    fixed size after composition (uniform with Figure 2).
    """
    print("Building Figure 1...")
    imgA = pdf_to_img(S1 / "pdf/B1_variance_partitioning_8_traits.pdf")
    imgB = pdf_to_img(S1 / "pdf/B2_reaction_norm_per_trait_highlights.pdf")

    pw = (A4_TEXT_W - PANEL_GAP) // 2
    a, b = resize_to_width(imgA, pw), resize_to_width(imgB, pw)
    h = max(a.height, b.height)
    fig = Image.new("RGB", (A4_TEXT_W, h), "white")
    fig.paste(a, (0, 0))
    fig.paste(b, (pw + PANEL_GAP, 0))
    stamp_label(fig, 0,              0, "A")
    stamp_label(fig, pw + PANEL_GAP, 0, "B")
    save(fig, OUT / "Figure_1", "Figure_1")
    print("Figure 1 done.\n")


def build_figure_2():
    """
    Figure 2 — Ecological architecture
      A  A4_site_boxplots_centered          (wide, landscape)
      B  A12_panelA_nutri_heatmap_rawP      (square)
      C  A12_panelB_nutri_morpho_dotplot    (square)
      D  C2r08_effects_standardized         (tall, portrait)
    Layout: 2 x 2 grid
            A (boxplots) | D (forest)     -- top row
            B (heatmap)  | C (dotplot)    -- bottom row
            Each row uses common-height sizing so the differently-shaped
            panels (wide A, tall D) sit flush with no whitespace.
    """
    print("Building Figure 2...")
    imgA = pdf_to_img(S1 / "pdf/A4_site_boxplots_centered.pdf")
    imgB = pdf_to_img(S1 / "pdf/A12_panelA_nutri_heatmap_rawP.pdf")
    imgC = pdf_to_img(S1 / "pdf/A12_panelB_nutri_morpho_dotplot_localFDR.pdf")
    imgD = pdf_to_img(C2  / "C2r08_effects_standardized.pdf")

    def fit_row(panels):
        """Resize panels to a common height so widths fill A4_TEXT_W; return
        (resized panels, their x-offsets, row height)."""
        ars = [p.width / p.height for p in panels]
        h   = int(round((A4_TEXT_W - PANEL_GAP * (len(panels) - 1)) / sum(ars)))
        rs  = [resize_to_height(p, h) for p in panels]
        xs, x = [], 0
        for p in rs:
            xs.append(x)
            x += p.width + PANEL_GAP
        return rs, xs, h

    (a, d), x1, h1 = fit_row([imgA, imgD])     # top:    A | D
    (b, c), x2, h2 = fit_row([imgB, imgC])     # bottom: B | C
    fig = Image.new("RGB", (A4_TEXT_W, h1 + ROW_GAP + h2), "white")
    fig.paste(a, (x1[0], 0));            fig.paste(d, (x1[1], 0))
    fig.paste(b, (x2[0], h1 + ROW_GAP)); fig.paste(c, (x2[1], h1 + ROW_GAP))
    stamp_label(fig, x1[0], 0,            "A")
    stamp_label(fig, x1[1], 0,            "D")
    stamp_label(fig, x2[0], h1 + ROW_GAP, "B")
    stamp_label(fig, x2[1], h1 + ROW_GAP, "C")
    save(fig, OUT / "Figure_2", "Figure_2")
    print("Figure 2 done.\n")


def build_figure_3():
    """
    Figure 3 — GWAS + population/LD diagnostics
      A  Figure_main_manhattan_qq_plain   (left, full height)
      B  Figure_S1b_PC_scree_withoutPC1   (right top)
      C  figure_ld_decay_genomewide       (right bottom)
    Layout: A occupies the full left column; B stacked above C on the right.
    Column widths are solved so both sides reach the same total height.
    """
    print("Building Figure 3...")
    imgA = add_label(pdf_to_img(GWAS / "main_figure/Figure_main_manhattan_qq_plain.pdf"),       "A")
    imgB = add_label(pdf_to_img(GWAS / "pc_selection/Figure_S1b_PC_scree_withoutPC1.pdf"),      "B")
    imgC = add_label(pdf_to_img(LD   / "figure_ld_decay_genomewide.pdf"),                       "C")

    # Aspect ratios (height/width) after label strips are added
    ar_a = imgA.height / imgA.width
    ar_b = imgB.height / imgB.width
    ar_c = imgC.height / imgC.width

    # Solve for right-column width so that:
    #   left_height  == right_height
    #   pw_a * ar_a  == pw_bc * (ar_b + ar_c) + ROW_GAP
    #   pw_a + PANEL_GAP + pw_bc == A4_TEXT_W
    # => pw_bc = ((A4_TEXT_W - PANEL_GAP) * ar_a - ROW_GAP) / (ar_a + ar_b + ar_c)
    pw_bc = int(((A4_TEXT_W - PANEL_GAP) * ar_a - ROW_GAP) / (ar_a + ar_b + ar_c))
    pw_bc = max(200, min(pw_bc, A4_TEXT_W - PANEL_GAP - 200))   # sanity clamp
    pw_a  = A4_TEXT_W - PANEL_GAP - pw_bc

    imgA = resize_to_width(imgA, pw_a)
    imgB = resize_to_width(imgB, pw_bc)
    imgC = resize_to_width(imgC, pw_bc)

    total_h = max(imgA.height, imgB.height + ROW_GAP + imgC.height)
    canvas  = Image.new("RGB", (A4_TEXT_W, total_h), "white")
    canvas.paste(imgA, (0, 0))
    canvas.paste(imgB, (pw_a + PANEL_GAP, 0))
    canvas.paste(imgC, (pw_a + PANEL_GAP, imgB.height + ROW_GAP))

    save(canvas, OUT / "Figure_3", "Figure_3")
    print("Figure 3 done.\n")


def build_figure_4():
    """
    Figure 4 — Haplotype structure of four candidate genes  (was Figure 5)
    Single R-generated 2×2 composite (panels a–d already labelled).
    Layout: resize to A4 text width, no re-assembly.
    """
    print("Building Figure 4...")
    img = pdf_to_img(HAP / "composite_4panel.pdf")
    img = resize_to_width(img, A4_TEXT_W)
    save(img, OUT / "Figure_4", "Figure_4")
    print("Figure 4 done.\n")


# ── Entry point ───────────────────────────────────────────────────────────────

BUILDERS = {
    "1": build_figure_1,
    "2": build_figure_2,
    "3": build_figure_3,
    "4": build_figure_4,
}

if __name__ == "__main__":
    target = sys.argv[1].lower() if len(sys.argv) > 1 else "1"
    if target == "all":
        for fn in BUILDERS.values():
            fn()
    elif target in BUILDERS:
        BUILDERS[target]()
    else:
        print(f"Unknown target '{target}'. Use 1–5 or 'all'.")
        sys.exit(1)
