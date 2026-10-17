#!/usr/bin/env python3
"""04_results_numbers.py -- check 02: key numbers in plain prose (results/tables/results_numbers.txt), from the tables
of scripts 01-02 and the raw phenotype table (read-only). Run after 01-03 (run_all.sh does)."""
import pathlib
import pandas as pd
HERE = pathlib.Path(__file__).resolve().parent.parent
TAB = HERE / "results/tables"
RAW = HERE.parent.parent / "00_THIN_Generate_Plots_For_Publication/all years barley.csv"
T1 = pd.read_csv(TAB / "T1_summary_by_season.tsv", sep="\t")
T3 = pd.read_csv(TAB / "T3_correlations_between_seasons.tsv", sep="\t")
T4 = pd.read_csv(TAB / "T4_correlations_between_traits.tsv", sep="\t")
raw = pd.read_csv(RAW, low_memory=False)
raw = raw[~raw.short_Tag.str.match(r"^[cCM]") & ~raw.short_Tag.str.startswith("HS04") & ~((raw.season == 2019) & (raw.Block_no == 4))]
neg = raw[raw.season == 2021]["StarchAsis."].dropna()
L = ["Check 02 -- raw phenotype checks: key numbers (scripts/04_results_numbers.py; tables in results/tables/).",
     "Panel: 290 accessions; 870 / 1,160 / 1,450 plants in 2019-20 / 2020-21 / 2021-22 (block 4 of 2019-20 removed); "
     "862 / 1,139 / 1,163 with NIR values (= M&M; asserted in 01_boxplots.R)."]
for tr in ["Protein", "Starch", "β-glucan", "Fiber"]:
    s = T1[T1.trait == tr]
    L.append(f"{tr}: season means " + " / ".join(f"{m:.1f}" for m in s["mean"]) + " %; values outside the ±3 SD limits per season "
             + " / ".join(str(int(a + b)) for a, b in zip(s.n_below_sd3, s.n_above_sd3)))
L.append(f"Starch 2021-22 (internal note, not for the supplement): as recorded, {int((neg < 0).sum())} of {len(neg)} plants below 0 % "
         f"(min {neg.min():.1f} %); with the +8.85 used here and in the analyses, {int(((neg + 8.85) < 0).sum())} remain below 0 "
         f"(min {neg.min() + 8.85:.1f} %).")
for tr in ["Protein", "Starch", "β-glucan", "Fiber"]:
    s = T3[(T3.trait == tr) & (T3.scope_2 != "All seasons")]; a = T3[(T3.trait == tr) & (T3.scope_2 == "All seasons")]
    L.append(f"{tr}: between-season r {s.r.min():.2f}-{s.r.max():.2f}; season vs all-season mean r {a.r.min():.2f}-{a.r.max():.2f}")
a = T4[T4.scope == "All seasons"].set_index(["trait_1", "trait_2"]).r
L.append("Between traits, all seasons (accession means of season-centered raw values): starch-fiber %.2f, starch-β-glucan %.2f, "
         "β-glucan-fiber %.2f, protein-starch %.2f (BLUP-based Fig. 1c of the manuscript: -0.78, -0.53, 0.35, -0.12)"
         % (a[("Starch", "Fiber")], a[("Starch", "β-glucan")], a[("β-glucan", "Fiber")], a[("Protein", "Starch")]))
(TAB / "results_numbers.txt").write_text("\n".join(L) + "\n", encoding="utf-8")
print("\n".join(L))
