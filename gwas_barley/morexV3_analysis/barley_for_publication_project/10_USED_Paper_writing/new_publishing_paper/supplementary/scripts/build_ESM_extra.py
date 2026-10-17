#!/usr/bin/env python3
"""build_ESM_extra.py -- the UNNUMBERED supplement, not cited in the manuscript (user request, 2026-10-01): descriptive
checks of the raw phenotypic data. Written to ../unreferenced/ESM_extra_phenotype_checks.{pdf,xlsx}, with the same TAG
title block as the numbered files. Kept apart so the numbered Online Resources stay TAG-compliant (every ESM must be
cited in the text); the user decides later whether any of it joins the numbered set.

Sources (read-only):
  12_USED_adittional_checks_for_supplementary/02_raw_phenotype_checks/results/ (figures E1-E4, E6, E7; tables T1-T5)
  12_USED_adittional_checks_for_supplementary/01_block4_2019_interblock_correlations/results/ (figure E5 =
      Fig_interblock_correlation_heatmaps.png; table interblock_correlations_all.tsv, raw rows)
Figure numbering E1-E7 follows the user's list: 1 boxplots (by season, by block); 2 correlations (between seasons,
between traits, between blocks); 3 distributions (raw values, BLUPs).

The starch values of 2021-22 carry the correction used in the analyses; it is not mentioned in the files (user, ESM0 Q12a).
Run: python3 build_ESM_extra.py   (after 12_.../02_raw_phenotype_checks/scripts/run_all.sh). Created 2026-10-01.
"""
import pathlib, shutil, tempfile
import pandas as pd
from esm_common import ROOT, SUPP, new_workbook, about_sheet, table_sheet, check_no_screen_categories
from build_ESM_pdf import title_page, figure_page, compile_pages, finish, TITLE_PLAIN

C2 = ROOT / "12_USED_adittional_checks_for_supplementary/02_raw_phenotype_checks/results"
C1 = ROOT / "12_USED_adittional_checks_for_supplementary/01_block4_2019_interblock_correlations/results"
OUT = SUPP / "unreferenced"
OUT.mkdir(exist_ok=True)
HEAD = "Supplementary material not cited in the text (draft, unnumbered): raw phenotype checks"

PANEL = ("Analysis panel: 290 accessions (site 04 and the cultivated checks removed); block 4 of 2019–20 removed except "
         "in Fig. E2. Values before centering within season and before the ±3 SD outlier rule of the analyses "
         "(values centered within season, pooled over seasons, kept within the mean ± 3 SD).")
FIGS = [  # label, file, paper, caption
    ("Fig. E1", C2 / "figures/Fig_E1_nutritional_traits_by_season.pdf", "a4paper",
     "**Fig. E1** Raw values of the four nutritional traits (% of grain, near-infrared spectroscopy) by season: **a** protein, "
     "**b** starch, **c** β-glucan, **d** fiber. Boxes, median and interquartile range; whiskers, 1.5 × the interquartile "
     "range; points, values beyond the whiskers; n, plants. Dashed lines, the limits of the ±3 SD outlier rule"),
    ("Fig. E2", C2 / "figures/Fig_E2_nutritional_traits_by_block.pdf", "a4paper",
     "**Fig. E2** Raw values of the four nutritional traits by block within each season (rows, traits; columns, seasons). "
     "Block 4 of 2019–20 (shaded) was excluded from the analyses because its β-glucan values correlated poorly with "
     "those of the other blocks (Fig. E5). Boxes as in Fig. E1"),
    ("Fig. E3", C2 / "figures/Fig_E3_correlations_between_seasons.pdf", "a4paper",
     "**Fig. E3** Pearson correlations of accession means between seasons, per nutritional trait: **a** protein, **b** "
     "starch, **c** β-glucan, **d** fiber. For each season, the accession means of the raw values in that season; all "
     "seasons, the accession means of the values centered within season. Accessions with values in both members of a "
     "pair (285–290)"),
    ("Fig. E4", C2 / "figures/Fig_E4_correlations_between_traits.pdf", "a4paper",
     "**Fig. E4** Pearson correlations between the four nutritional traits, from accession means: **a** over all seasons "
     "(values centered within season; 290 accessions), **b**–**d** within each season (raw values; 287–290 accessions)"),
    ("Fig. E5", C1 / "figures/Fig_interblock_correlation_heatmaps.png", "a4paper",
     "**Fig. E5** Pearson correlations between blocks within each season for the four nutritional traits (raw values; "
     "block 4 of 2019–20 included): **a** 2019–20, **b** 2020–21, **c** 2021–22. Each pair of blocks over the accessions "
     "with values in both (171–288)"),
    ("Fig. E6", C2 / "figures/Fig_E6_raw_distributions_by_season.pdf", "a4paper",
     "**Fig. E6** Distributions of the raw values of the four nutritional and the six morphological and phenological "
     "traits, by season (frequency polygons, 40 bins). The x axis spans the widest ±3 SD limits of the three seasons; "
     "each panel gives the number of values outside it. Number of tillers was recorded in 2019–20 and 2020–21, grain "
     "weight in 2019–20 and 2021–22, and number of grains in 2019–20 only"),
    ("Fig. E7", C2 / "figures/Fig_E7_BLUP_distributions.pdf", "a4paper",
     "**Fig. E7** Distributions of the genotype BLUPs (centered; units of the trait) of the ten traits, the phenotypes of "
     "all analyses (290 accessions)"),
]


def build_pdf():
    work = pathlib.Path(tempfile.mkdtemp(prefix="esm_extra_", dir="/mnt/data/shahar/.tmp"))
    legend = ["Descriptive checks of the raw phenotypic data of the common-garden experiment (three seasons, 2019–20 to "
              "2021–22): the four nutritional traits by season and by block (Figs. E1, E2), correlations between seasons, "
              "between traits and between blocks (Figs. E3–E5), and the distributions of the raw values and of the "
              "genotype BLUPs (Figs. E6, E7). " + PANEL + " Tables: ESM_extra_phenotype_checks.xlsx"]
    docs = [title_page(HEAD, legend, [f"Page {i + 2}: {f[0]}" for i, f in enumerate(FIGS)])] + \
           [figure_page(HEAD, lab, (path, None), cap, paper, "width=174mm,height=222mm,keepaspectratio") for lab, path, paper, cap in FIGS]
    out = finish("extra", compile_pages("extra", docs, work), work, f"{TITLE_PLAIN} - raw phenotype checks (unnumbered)",
                 out=OUT / "ESM_extra_phenotype_checks.pdf")
    shutil.rmtree(work)
    return out


def build_xlsx():
    T = {k: pd.read_csv(C2 / f"tables/{k}.tsv", sep="\t") for k in
         ["T1_summary_by_season", "T2_summary_by_block", "T3_correlations_between_seasons", "T4_correlations_between_traits", "T5_BLUP_summary"]}
    B = pd.read_csv(C1 / "tables/interblock_correlations_all.tsv", sep="\t")
    B = B[B.analysis == "raw"].drop(columns=["analysis"])
    B["season"] = B.season.map({2019: "2019–20", 2020: "2020–21", 2021: "2021–22"})
    B["trait"] = B.trait.replace({"Beta-glucan": "β-glucan"})
    wb = new_workbook()
    sheets = [("Summary by season", "Raw values of the ten traits per season (Figs. E1, E6)."),
              ("Summary by block", "Raw values of the four nutritional traits per season and block (Fig. E2)."),
              ("Correlations seasons", "Correlations of accession means between seasons, per nutritional trait (Fig. E3)."),
              ("Correlations traits", "Correlations between the nutritional traits, all seasons and per season (Fig. E4)."),
              ("Correlations blocks", "Correlations between blocks within each season (Fig. E5)."),
              ("BLUP summary", "The genotype BLUPs of the ten traits (Fig. E7).")]
    defs = {"Summary by season": [("n_plants, n_accessions", "Plants with a value, and their accessions."),
                                  ("mean, sd, cv_pct, min, q1, median, q3, max", "Of the raw values (units in column unit)."),
                                  ("sd3_lower, sd3_upper", "Limits of the ±3 SD outlier rule of the analyses, on the raw scale of the season."),
                                  ("n_below_sd3, n_above_sd3", "Values outside those limits (removed in the analyses).")],
            "Correlations seasons": [("scope_1, scope_2", "Season (accession means of the raw values) or all seasons (accession means of the values centered within season)."),
                                     ("r, n, p", "Pearson r, accessions, two-sided P (no correction for multiple testing).")],
            "Correlations traits": [("scope", "All seasons (season-centered values) or one season (raw values)."),
                                    ("r, n, p", "As above.")],
            "Correlations blocks": [("block_a, block_b", "The two blocks (one plant per accession in each)."),
                                    ("n, r, ci_low, ci_high, p", "Accessions with values in both blocks, Pearson r with its 95% confidence interval (Fisher z), two-sided P.")]}
    cap = ("Descriptive checks of the raw phenotypic data (tables of Figs. E1–E7 of the PDF of the same name). " + PANEL)
    about_sheet(wb, HEAD, cap, sheets, defs)
    for (name, _), key in zip(sheets, ["T1_summary_by_season", "T2_summary_by_block", "T3_correlations_between_seasons",
                                       "T4_correlations_between_traits", None, "T5_BLUP_summary"]):
        df = B if key is None else T[key]
        num = [c for c in df.columns if pd.api.types.is_float_dtype(df[c])]
        table_sheet(wb, name, list(df.columns), df.values.tolist(),
                    fmt={c: ("0.00E+00" if c == "p" else "0.000") for c in num}, widths={"trait": 16, "season": 10}, freeze="B2")
    check_no_screen_categories(wb)
    out = OUT / "ESM_extra_phenotype_checks.xlsx"
    wb.save(out)
    print("wrote", out)


if __name__ == "__main__":
    build_xlsx()
    build_pdf()
