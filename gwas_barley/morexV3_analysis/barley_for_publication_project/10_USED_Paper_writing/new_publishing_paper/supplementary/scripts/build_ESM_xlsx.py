#!/usr/bin/env python3
"""build_ESM_xlsx.py -- builds the spreadsheet Online Resources of the TAG manuscript, ESM_1.xlsx ... ESM_11.xlsx
(numbering = order of first citation in the manuscript; final list approved by the user 2026-10-01, ESM0 Q1).

NOTHING IS RE-COMPUTED. Every value is copied from a read-only project output; joins use gene_id / lead_SNP, never
locus_id alone (10_USED_Paper_writing/CLAUDE.md, hard rule 3). Column renaming, unit conversion of call rates from
fractions to %, and dropping columns are the only transformations. Each builder asserts the numbers the manuscript
cites from it (the script stops if one differs).

  ESM 1  sampling sites, accessions (ENA, genotype BLUPs), site environments, environmental correlations
  ESM 2  NIR data quality: validation of the NIR calibration (30 calibration accessions; ESM0 Q17) and inter-block
         correlations 2019-20
  ESM 3  variance components and H² (table from make_ESM_variance_components.py)
  ESM 5  the 36 loci and the 52 significant SNPs
  ESM 6  the 55 candidate genes (table from make_ESM_candidate_genes.py)
  ESM 7  choice of the crosshap epsilon
  ESM 8  canonical (1,3;1,4)-beta-glucan genes on MorexV3
  ESM 10 elite cultivars: lines, 136-line screen, site accounting, removed sites
  ESM 11 geographic origin of allele and haplotype carriers (tables)
(ESM 4, 9, 12, 13 are PDFs: build_ESM_pdf.py.)

Sources (read-only; project root = esm_common.ROOT) are listed in each builder below and in README.md.
User decisions applied (build/proposed/ESM0_QUESTIONS.md, 2026-10-01): Q7 desert_OR and the extreme-accession check
dropped from ESM 11; Q11(d) electrical conductivity in mS/cm (as the source table and Fig. 1d); no project paths in
the published files; relevance-screen categories never written (esm_common.check_no_screen_categories).
ESM0 Q11 (user, 2026-10-01): ENA sample/run accessions kept in ESM 1; no unit for plant height and grain weight (as in
Potapenko et al. 2026b, Mol Ecol); BLUPs on the season-centered scale, as used in all analyses.

Run: python3 build_ESM_xlsx.py   (after make_ESM_variance_components.py and make_ESM_candidate_genes.py; see
build_all.sh). Created 2026-10-01.
"""
import csv, statistics, sys
import pandas as pd
from esm_common import (ROOT, SUPP, new_workbook, about_sheet, table_sheet, save, trait_label, strip_comment_header)

T00 = ROOT / "00_THIN_Generate_Plots_For_Publication"
TAB00 = T00 / "outputs/subsection_1/tables"
LOCI = ROOT / "01_USED_GWAS_V2_pipeline/results/00_FINAL_BLUP_3PC"
EPS = ROOT / "04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Diagnostics"
CANON = ROOT / "05_USED_gene_annotation_analysis/09_USED_canonical_betaglucan_gene_check/results/tables"
ELITE = ROOT / "07_USED_elite_lines_compariosn_to_wild_lines/results/tables"
B7H = ROOT / "03_01_7H_branch_Starch_Fiber_shared_signal_explore/results/tables"
ORIG = ROOT / "11_USED_check_low_maf_geographic_distribution_minor_allele/results/tables"
BLOCK = ROOT / "12_USED_adittional_checks_for_supplementary/01_block4_2019_interblock_correlations/results/tables"
ENA = SUPP / "sources/ENA_PRJEB79623_read_run_2026-10-01.tsv"

rd = lambda p, **k: pd.read_csv(p, sep="\t", **k)
CHECKS = []


def rh(x, nd=2):
    """round half up, as the manuscript text does (0.605 -> 0.61); float round() would give 0.6"""
    from decimal import Decimal, ROUND_HALF_UP
    return float(Decimal(str(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_HALF_UP))


def ok(cond, what):
    assert cond, "CHECK FAILED: " + what
    CHECKS.append(what)


# trait display names and units, as in the Materials and methods. Units only where Potapenko et al. (2026b, Mol Ecol)
# give them (user, ESM0 Q11b, 2026-10-01: "keep it as in Mol Ecol"): no unit for plant height and grain weight.
TRAITS10 = [  # BLUP_10_traits.csv column, label, unit
    ("ProteinCenterd", "Protein", "%"), ("StarchCenterd", "Starch", "%"), ("BetaglucansCenterd", "β-glucan", "%"),
    ("FiberCenterd", "Fiber", "%"), ("flowering_timeCenterd", "Flowering time", "days"),
    ("Tiller.LengthCenterd", "Plant height", ""), ("Spike.LengthCenterd", "Spike length", "cm"),
    ("Number.Of.TillersCenterd", "Number of tillers", "count"), ("Grain_weightCenterd", "Grain weight", ""),
    ("number_of_grainsCenterd", "Number of grains", "count")]
REGION_LABEL = {"North": "North", "Coast": "Coast", "Desert": "Desert", "HZ1": "HZ1 (North–Coast)",
                "HZ2": "HZ2 (North–Desert)", "HZ3": "HZ3 (Coast–Desert)"}


# =================================================================================================================
def esm1():
    """sites, accessions, environment. Sources: 00_THIN_.../01_phenotypic_analysis_no_GxE_v3_streamlined.R
    (SITE_NAME_MAP, SITE_REGION_MAP: copied below), all years barley.csv (longitude, latitude),
    merged_Evgeny_env_fitness_environmental_30_sites.csv (environment, 30 sites), outputs/subsection_1/tables/
    BLUP_10_traits.csv, C2_trait_environment_correlations/tables/C2corr_site_level_dataset.csv (check),
    C2_environmental_predictor_collinearity/tables/C2col_environment_correlation_matrix.csv, ENA report."""
    SITE_NAME = {"01": "Meizar", "02": "Beit Zeida", "03": "Mount Hermon", "04": "Hachola", "05": "Misgav River",
                 "06": "Zichron Yaakov", "07": "Mount Tayasim", "08": "Rohama", "09": "Nahala", "10": "Zikim",
                 "11": "Masciot", "12": "Mount Ramon", "13": "Masada", "14": "Maale Akrabim", "15": "Arad",
                 "16": "Cabri", "17": "Meron", "18": "Alon Hagalil", "19": "Afek", "20": "Rechasim", "21": "Meirav",
                 "22": "Sde Ilan", "23": "Gesher", "24": "Givat Hamoreh", "25": "Tirat Zvi", "26": "Megido Forest",
                 "27": "Revivim", "28": "Orim", "29": "Hagosh Junction", "30": "Amona"}   # 04 from the 30-site file
    env = pd.read_csv(T00 / "merged_Evgeny_env_fitness_environmental_30_sites.csv")
    env["site"] = env.site_id.str[2:]
    ph = pd.read_csv(T00 / "all years barley.csv", low_memory=False)
    ph = ph[ph.site.notna()]
    ph["site2"] = ph.site.astype(int).map("{:02d}".format)
    coords = ph.groupby("site2")[["latitude", "longitude"]].agg(lambda x: sorted(set(x.dropna())))
    ok(all(len(v) == 1 for v in coords.latitude) and all(len(v) == 1 for v in coords.longitude) and len(coords) == 30,
       "ESM1: one latitude/longitude pair per site, 30 sites")
    blup = pd.read_csv(TAB00 / "BLUP_10_traits.csv")
    blup["site"] = blup.Genotype.str[:2]
    ok(len(blup) == 290 and blup.Genotype.is_unique and (blup.site.value_counts() == 10).all()
       and blup.site.nunique() == 29 and "04" not in set(blup.site), "ESM1: 290 accessions, 29 sites x 10, no site 04")
    # BLUPs = the GWAS phenotypes (data/inputs/<trait>_corrected_V3.pheno)
    for col, tr in [("ProteinCenterd", "protein"), ("StarchCenterd", "starch"), ("BetaglucansCenterd", "betaglucan"),
                    ("FiberCenterd", "fiber")]:
        ph_g = pd.read_csv(f"/mnt/data/shahar/gwas_barley/data/inputs/{tr}_corrected_V3.pheno", sep=r"\s+",
                           header=None, names=["fid", "iid", "y"])
        ph_g["acc"] = ph_g.iid.str[2:4] + "_" + ph_g.iid.str[4:6]
        m = blup.merge(ph_g, left_on="Genotype", right_on="acc")
        ok(len(m) == 290 and (m[col] - m.y).abs().max() < 1e-9, f"ESM1: {tr} BLUPs = GWAS phenotype file (290)")
    # site means = the Fig. 1d input
    c2 = pd.read_csv(T00 / "outputs/subsection_1/C2_trait_environment_correlations/tables/C2corr_site_level_dataset.csv")
    sm = blup.groupby("site")[["ProteinCenterd", "StarchCenterd", "BetaglucansCenterd", "FiberCenterd"]].mean()
    sm.index = "HS" + sm.index
    c2i = c2.set_index("site_id")
    ok(max((sm[t] - c2i.loc[sm.index, t]).abs().max() for t in sm.columns) < 1e-8, "ESM1: site-mean BLUPs = Fig. 1d input")
    ENVCOLS = [("Elevation", "Elevation (m)"), ("March..temp", "March maximum temperature (°C)"),
               ("Precipitation", "Annual precipitation (mm)"), ("pH", "pH"),
               ("Electrical.Conductivity..mS.cm.", "Electrical conductivity (mS/cm)"), ("Clay....", "Clay (%)"),
               ("Silt....", "Silt (%)"), ("Sand....", "Sand (%)"), ("Organic.Carbon..mg.Kg.", "Organic carbon (mg/kg)"),
               ("Total.N..mg.Kg.", "Total N (mg/kg)")]
    C2MAP = {"pH": "pH", "Electrical.Conductivity..mS.cm.": "Electrical_Cond", "Clay....": "Clay_pct",
             "Sand....": "Sand_pct", "March..temp": "March_Temp", "Precipitation": "Precipitation",
             "Organic.Carbon..mg.Kg.": "Organic_Carbon", "Total.N..mg.Kg.": "Total_N"}
    e29 = env.set_index("site_id").loc[c2.site_id]
    ok(all((e29[s] - c2i[c]).abs().max() < 1e-9 for s, c in C2MAP.items()), "ESM1: 8 Fig. 1d variables identical in both sources")
    region = env.set_index("site").Polygon.to_dict()
    n_acc = blup.site.value_counts().to_dict()
    site_rows = []
    for _, r in env.sort_values("site").iterrows():
        s = r.site
        site_rows.append([s, SITE_NAME[s], REGION_LABEL[region[s]], coords.loc[s, "latitude"][0],
                          coords.loc[s, "longitude"][0], n_acc.get(s, 0)] + [r[c] for c, _ in ENVCOLS])
    ok(sum(1 for r in site_rows if r[5] == 10) == 29 and [r for r in site_rows if r[5] == 0][0][1] == "Hachola",
       "ESM1: 29 sites with 10 accessions; site 04 (Hachola) none")
    # accessions + ENA
    ena = rd(ENA)
    ok(ena.sample_alias.is_unique and len(ena) == 300, "ESM1: ENA PRJEB79623, 300 runs, unique aliases")
    ena = ena.set_index("sample_alias")
    ok(set(blup.Genotype) <= set(ena.index) and sorted(set(ena.index) - set(blup.Genotype))
       == ["04_01", "04_02", "04_04", "04_05", "04_06", "04_07", "04_08", "04_09", "04_10", "04_11"],
       "ESM1: all 290 accessions are ENA sample aliases; the 10 others are site 04")
    acc_rows = []
    for _, r in blup.sort_values("Genotype").iterrows():
        a = r.Genotype
        acc_rows.append([a, "HS" + a.replace("_", ""), ena.loc[a, "sample_accession"], ena.loc[a, "run_accession"],
                         r.site, REGION_LABEL[region[r.site]]] + [r[c] for c, _, _ in TRAITS10])
    # environmental correlation matrix (29 sites); recomputed from the 30-site file as a check only
    cm = pd.read_csv(T00 / "outputs/subsection_1/C2_environmental_predictor_collinearity/tables/C2col_environment_correlation_matrix.csv")
    CMNAME = {"pH": "pH", "Electrical_Cond": "Electrical conductivity", "Clay_pct": "Clay", "Silt_pct": "Silt",
              "Sand_pct": "Sand", "March_Temp": "March maximum temperature", "Precipitation": "Annual precipitation",
              "Organic_Carbon": "Organic carbon", "Total_N": "Total N", "Elevation": "Elevation"}
    SRC = {"pH": "pH", "Electrical_Cond": "Electrical.Conductivity..mS.cm.", "Clay_pct": "Clay....", "Silt_pct": "Silt....",
           "Sand_pct": "Sand....", "March_Temp": "March..temp", "Precipitation": "Precipitation",
           "Organic_Carbon": "Organic.Carbon..mg.Kg.", "Total_N": "Total.N..mg.Kg.", "Elevation": "Elevation"}
    rec = env[env.site != "04"][[SRC[v] for v in cm.predictor]].corr()
    ok(all(abs(rec.loc[SRC[a], SRC[b]] - cm.set_index("predictor").loc[a, b]) < 1e-5 for a in cm.predictor for b in cm.predictor),
       "ESM1: correlation matrix = Pearson r over the 29 analyzed sites")
    cmi = cm.set_index("predictor")
    ok(round(cmi.loc["March_Temp", "Elevation"], 3) == -0.951 and round(cmi.loc["Sand_pct", "Silt_pct"], 3) == -0.934,
       "ESM1: |r| > 0.8 pairs: March temperature-elevation -0.951, sand-silt -0.934 (M&M)")
    big = [(a, b) for i, a in enumerate(cm.predictor) for b in list(cm.predictor)[i + 1:] if abs(cmi.loc[a, b]) > 0.8]
    ok(sorted(map(sorted, big)) == sorted(map(sorted, [("March_Temp", "Elevation"), ("Silt_pct", "Sand_pct")])),
       "ESM1: exactly two pairs with |r| > 0.8")
    cor_rows = [[CMNAME[a]] + [float(cmi.loc[a, b]) for b in cm.predictor] for a in cm.predictor]

    wb = new_workbook()
    sheets = [("Sites", "The 30 sampling sites: ecological region, coordinates, accessions analyzed and ten environmental variables."),
              ("Accessions", "The 290 accessions: ENA identifiers, site, region and the genotype BLUPs of the ten traits."),
              ("Environment correlations", "Pearson correlations among the ten environmental variables across the 29 sites of the analyzed accessions.")]
    defs = {
        "Sites": [("Site", "Sampling-site number (first two digits of the accession ID)."),
                  ("Site name", "Locality of the sampling site, Israel."),
                  ("Region", "Ecological region of the sampling design: North, Coast and Desert (the three ecotype regions) "
                             "and the intermediate zones HZ1 (North–Coast), HZ2 (North–Desert) and HZ3 (Coast–Desert)."),
                  ("Latitude (°N), Longitude (°E)", "Coordinates of the sampling site (decimal degrees, WGS84)."),
                  ("Accessions analyzed", "Accessions of the site in the analysis panel. Site 04 (Hachola): its accessions were "
                                          "removed before the analyses (Materials and methods); its environment is listed for completeness."),
                  ("Elevation; March maximum temperature; Annual precipitation", "Elevation (m above sea level), maximum temperature "
                   "in March (°C) and average annual precipitation (mm), from the Israeli Meteorological Service."),
                  ("pH; Electrical conductivity; Clay; Silt; Sand; Organic carbon; Total N", "Soil properties measured in a "
                   "soil sample from each site. Elevation and silt were not used in the trait–environment correlations "
                   "(|r| > 0.8 with March maximum temperature and sand; sheet “Environment correlations”).")],
        "Accessions": [("Accession", "Accession ID used throughout the study; identical to the sample alias in the "
                                     "European Nucleotide Archive (ENA), project PRJEB79623."),
                       ("Sample name in the MorexV3 VCF", "Name of the accession in the MorexV3 SNP calls."),
                       ("ENA sample accession; ENA run accession", "ENA identifiers of the accession's whole-genome reads "
                                                                    "(project PRJEB79623; retrieved 1 October 2026)."),
                       ("Site; Region", "Sampling site and ecological region (sheet “Sites”)."),
                       ("<trait> BLUP (unit)", "Genotype best linear unbiased prediction of the trait, from a linear mixed model "
                        "with genotype and season-by-block as random effects fitted to the values centered within season "
                        "(Materials and methods); a deviation from the season mean, in the units of the trait. These are the "
                        "phenotypes of all analyses (genome-wide association, haplotype tests, correlations). Protein, starch, "
                        "β-glucan and fiber: % of grain, by near-infrared spectroscopy. Number of grains: the grains counted on "
                        "the spikes scored per plant (2019–20 only). Number of tillers: 2019–20 and 2020–21; grain weight: 2019–20 "
                        "and 2021–22; the other traits: all three seasons.")],
        "Environment correlations": [("Matrix", "Pearson r between each pair of environmental variables across the 29 "
                                                "sampling sites of the analyzed accessions (site 04 excluded).")],
    }
    cap = ("Sampling sites, accessions and site environments. The 30 sampling sites with their ecological region, "
           "coordinates and ten environmental variables; the 290 analyzed accessions with their identifiers in the "
           "European Nucleotide Archive (project PRJEB79623), site and genotype BLUPs of the ten traits; and the Pearson "
           "correlations among the environmental variables across the 29 sites")
    about_sheet(wb, 1, cap, sheets, defs)
    H1 = ["Site", "Site name", "Region", "Latitude (°N)", "Longitude (°E)", "Accessions analyzed"] + [h for _, h in ENVCOLS]
    table_sheet(wb, "Sites", H1, site_rows, fmt={"Latitude (°N)": "0.0000", "Longitude (°E)": "0.0000",
                "Elevation (m)": "0", "March maximum temperature (°C)": "0.00", "Annual precipitation (mm)": "0.0", "pH": "0.00",
                "Electrical conductivity (mS/cm)": "0", "Clay (%)": "0.0", "Silt (%)": "0.0", "Sand (%)": "0.0",
                "Organic carbon (mg/kg)": "0.0", "Total N (mg/kg)": "0.0"}, widths={"Site name": 17, "Region": 19}, freeze="C2")
    H2 = ["Accession", "Sample name in the MorexV3 VCF", "ENA sample accession", "ENA run accession", "Site", "Region"] + \
         [f"{lab} BLUP ({u})" if u not in ("count", "") else f"{lab} BLUP" for _, lab, u in TRAITS10]
    table_sheet(wb, "Accessions", H2, acc_rows, fmt={h: "0.0000" for h in H2[6:]}, widths={"Region": 19,
                "ENA sample accession": 16, "ENA run accession": 14}, freeze="B2")
    H3 = ["Variable"] + [CMNAME[v] for v in cm.predictor]
    table_sheet(wb, "Environment correlations", H3, cor_rows, fmt={h: "0.000" for h in H3[1:]}, widths={"Variable": 26}, freeze="B2")
    return save(wb, 1)


# =================================================================================================================
def pearson(x, y):
    """Pearson r, n and two-sided P (t test, df = n - 2), as R's cor.test; complete pairs only. Written out because
    scipy is not installed; checked against cor.test (thesis Table 1 values)."""
    import math
    pairs = [(a, b) for a, b in zip(x, y) if a == a and b == b]
    n = len(pairs); mx = sum(a for a, _ in pairs) / n; my = sum(b for _, b in pairs) / n
    sxy = sum((a - mx) * (b - my) for a, b in pairs); sxx = sum((a - mx) ** 2 for a, _ in pairs); syy = sum((b - my) ** 2 for _, b in pairs)
    r = sxy / math.sqrt(sxx * syy); df = n - 2
    if abs(r) >= 1:
        return r, n, 0.0
    t2 = r * r * df / (1 - r * r)
    def betacf(a, b, x, itmax=300, eps=3e-16):
        qab, qap, qam = a + b, a + 1, a - 1; c, d = 1.0, 1 - qab * x / qap
        d = 1 / (d if abs(d) > 1e-300 else 1e-300); h = d
        for m in range(1, itmax + 1):
            m2 = 2 * m; aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / c if abs(1 + aa / c) > 1e-300 else 1e-300; h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1 + aa * d; d = 1 / (d if abs(d) > 1e-300 else 1e-300); c = 1 + aa / c if abs(1 + aa / c) > 1e-300 else 1e-300
            de = d * c; h *= de
            if abs(de - 1) < eps:
                break
        return h
    def ibeta(a, b, x):
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
        return bt * betacf(a, b, x) / a if x < (a + 1) / (a + b + 2) else 1 - bt * betacf(b, a, 1 - x) / b
    return r, n, ibeta(df / 2, 0.5, df / (df + t2))


def esm2():
    """NIR data quality. Sheets:
    (1) "NIR calibration" + (2) "Calibration correlations" (user, ESM0 Q17 option A, 2026-10-01): the 30 calibration
        accessions with their laboratory values, the NIR predictions of the same grain with the factory calibration
        (mean of 5 scans) and their 2019-20 values, and the Pearson correlations. Sources (copied unchanged from N. Pintel's
        files, ../sources/NIR_calibration/, README there): kits.results.xlsx, calibration.scans.xlsx, season2020.xlsx,
        DA7250_Extensive Report_..._avg.xlsx. Left out (Q17): beta-glucan lab vs calibration scans (r = 1 by construction:
        that column is a copy of the kit values) and fiber (no laboratory reference; the thesis r = 0.526 does not
        reproduce from the source data, r = 0.25).
    (3) "Blocks 2019-20": 12_USED_adittional_checks_for_supplementary/01_block4_2019_interblock_correlations/results/tables/
        Table_OnlineResource_2019_block_correlations.tsv (raw values; run 2026-09-29)."""
    NIRSRC = SUPP / "sources/NIR_calibration"
    K = pd.read_excel(NIRSRC / "kits.results.xlsx"); C = pd.read_excel(NIRSRC / "calibration.scans.xlsx")
    S = pd.read_excel(NIRSRC / "season2020.xlsx")
    ok(len(K) == 30 and list(K.GI) == list(C.GI) == list(S.GI) and list(K["Sample ID"]) == list(C["Sample ID"]) == list(S["Sample ID"]),
       "ESM2: 30 calibration accessions, same order in the three files")
    ids = K["Sample ID"].str.extract(r"(\d+)_(\d+)$").astype(int)
    acc = [f"{a:02d}_{b:02d}" for a, b in zip(ids[0], ids[1])]
    ok(sorted(ids[0]) == list(range(1, 31)), "ESM2: one calibration accession per site (sites 1-30)")
    # calibration scans = mean of the 5 device scans (factory calibration, bias 0, slope 1)
    X = pd.read_excel(NIRSRC / "DA7250_Extensive Report_Barley_20200107_0000_To_20200107_1324_avg.xlsx")
    X["site"] = X["Sample ID"].astype(str).str.extract(r"(\d+)_\d+")[0].astype(int)
    g = X.groupby("site").agg(n=("site", "size"), prot=("Protein As is %, Predicted", "mean"), starch=("Starch As is %, Predicted", "mean"))
    ok((g.n == 5).all() and all(abs(g.loc[s, "prot"] - v) < 0.006 for s, v in zip(ids[0], C["Protein %"]))
       and all(abs(g.loc[s, "starch"] - v) < 0.006 for s, v in zip(ids[0], C["Starch %"])),
       "ESM2: calibration-scan values = mean of the 5 DA 7250 scans per accession")
    ok(set(X["Calib File Name Protein"]) == {"Barley_Protein_Asis_7July2011_Micro.41M"} and set(X["Bias Protein"]) == {0}
       and set(X["Slope Starch"]) == {1}, "ESM2: calibration scans used the factory cultivated-barley calibrations, unadjusted")
    ok((K["Beta Glucan %"] == C["Beta Glucan %"]).all(), "ESM2: beta-glucan 'calibration scan' column = kit values (left out)")
    # 2019-20 values = one plant (block 1) of the study's data
    ph = pd.read_csv(T00 / "all years barley.csv", low_memory=False)
    ph = ph[(ph.season == 2019) & (ph.Block_no == 1)].set_index("short_Tag")
    tags = ["HS" + a.replace("_", "") for a in acc]
    ok(all(abs(ph.loc[tg, "ProteinAsis."] - pv) < 0.006 and abs(ph.loc[tg, "StarchAsis."] - sv) < 0.051 and abs(ph.loc[tg, "BetaglucansAsis."] - bv) < 0.006
           for tg, pv, sv, bv in zip(tags, S["Protein %"], S["Starch %"], S["Beta Glucan %"])),
       "ESM2: 2019-20 values = block-1 plant of the same accession in the study's data (30 of 30)")
    rows = [[a, int(s), K.loc[i, "Protein %"], C.loc[i, "Protein %"], S.loc[i, "Protein %"], K.loc[i, "Starch %"], C.loc[i, "Starch %"],
             S.loc[i, "Starch %"], K.loc[i, "Beta Glucan %"], S.loc[i, "Beta Glucan %"]] for i, (a, s) in enumerate(zip(acc, ids[0]))]
    rows.sort(key=lambda r: r[0])
    HC = ["Accession", "Site", "Protein, laboratory (%)", "Protein, NIR calibration scan (%)", "Protein, NIR 2019–20 (%)",
          "Starch, laboratory (%)", "Starch, NIR calibration scan (%)", "Starch, NIR 2019–20 (%)", "β-glucan, laboratory (%)",
          "β-glucan, NIR 2019–20 (%)"]
    THESIS = {("Protein", "lab", "cal"): 0.915, ("Protein", "lab", "s20"): 0.556, ("Protein", "s20", "cal"): 0.649,
              ("Starch", "lab", "cal"): 0.785, ("Starch", "lab", "s20"): 0.466, ("Starch", "s20", "cal"): 0.548,
              ("β-glucan", "lab", "s20"): 0.514}
    col = {"Protein": "Protein %", "Starch": "Starch %", "β-glucan": "Beta Glucan %"}
    src = {"lab": K, "cal": C, "s20": S}
    LAB = {"lab": "laboratory", "cal": "NIR calibration scan", "s20": "NIR 2019–20"}
    crow = []
    for (tr, a, b), rt in THESIS.items():
        r, n, pv = pearson(list(src[a][col[tr]]), list(src[b][col[tr]]))
        ok(rh(r, 3) == rt and n == 30, f"ESM2: {tr} {a} vs {b} r = {rt} (thesis Table 1), n = 30")
        crow.append([tr, LAB[a], LAB[b], r, n, pv])
    r, n, pv = pearson(list(K["Protein %"]), list(C["Protein %"]))
    ok(abs(pv - 1.49e-12) / 1.49e-12 < 0.01, "ESM2: P of the t test = R cor.test (protein, 1.49e-12)")
    HR = ["Trait", "Values 1", "Values 2", "Pearson r", "n", "P"]
    t = rd(BLOCK / "Table_OnlineResource_2019_block_correlations.tsv", dtype={"block_pair": str})
    b4 = t[t.block_pair.str.contains("4")]; b13 = t[~t.block_pair.str.contains("4")]
    ok(len(t) == 6 and rh(b4["Beta-glucan_r"].min()) == 0.27 and rh(b4["Beta-glucan_r"].max()) == 0.34
       and rh(b13["Beta-glucan_r"].min()) == 0.61 and rh(b13["Beta-glucan_r"].max()) == 0.64,
       "ESM2: beta-glucan r = 0.27-0.34 with block 4 vs 0.61-0.64 among blocks 1-3 (M&M)")
    brows = [[f"Block {p.split('-')[0]} – block {p.split('-')[1]}"] + [r_[f"{tr}_{k}"] for tr in ["Protein", "Starch", "Beta-glucan", "Fiber"] for k in ("r", "n")]
             for p, (_, r_) in zip(t.block_pair, t.iterrows())]
    HB = ["Block pair"] + [f"{lab} {k}" for lab in ["Protein", "Starch", "β-glucan", "Fiber"] for k in ("r", "n")]
    wb = new_workbook()
    cap = ("Quality of the near-infrared (NIR) data. Validation of the NIR calibration: for the 30 calibration accessions "
           "(one per sampling site), the laboratory values, the NIR predictions of the same grain with the manufacturer's "
           "calibration for cultivated barley, and the NIR values of the same accessions in the 2019–20 season, with their "
           "Pearson correlations. Pearson correlations between the blocks of the 2019–20 season for the four nutritional "
           "traits, from the raw NIR values of the accessions grown in both blocks of each pair; block 4 correlated poorly "
           "with the other blocks for β-glucan and was excluded from the analyses")
    about_sheet(wb, 2, cap, [
        ("NIR calibration", "The 30 calibration accessions: laboratory values and NIR predictions (%)."),
        ("Calibration correlations", "Pearson correlations between the laboratory values and the NIR predictions."),
        ("Blocks 2019–20", "Pearson r and number of accessions for each pair of blocks and each trait.")], {
        "NIR calibration": [("Accession, Site", "Calibration accession (one chosen at random per sampling site; ID as in Online Resource 1) and its site. "
                             "Site 04 is included: the calibration used all 30 sites."),
                            ("laboratory", "Protein: Kjeldahl, nitrogen × 6.25, mean of three subsamples of about 200 mg. Starch and β-glucan: "
                             "Megazyme Total Starch and Mixed-Linkage β-Glucan kits, mean of four subsamples (two for the starch of one accession)."),
                            ("NIR calibration scan", "Prediction for the same grain by the DA 7250 with the manufacturer's calibration for "
                             "cultivated barley, mean of five scans. There was no β-glucan calibration for this step: the β-glucan "
                             "calibration was built from the laboratory values."),
                            ("NIR 2019–20", "NIR value of the same accession in the 2019–20 season (one plant, block 1).")],
        "Calibration correlations": [("Values 1, Values 2", "The two sets of values compared (sheet “NIR calibration”)."),
                                     ("Pearson r, n, P", "Pearson correlation over the 30 accessions and its two-sided P.")],
        "Blocks 2019–20": [("Block pair", "The two blocks compared (one plant per accession in each block)."),
                           ("<trait> r", "Pearson correlation coefficient between the values of the same accessions in the two blocks (raw values, % of grain)."),
                           ("<trait> n", "Number of accessions with a value in both blocks.")]})
    table_sheet(wb, "NIR calibration", HC, rows, fmt={h: "0.00" for h in HC[2:]}, freeze="C2")
    table_sheet(wb, "Calibration correlations", HR, crow, fmt={"Pearson r": "0.000", "P": "0.0E+00"},
                widths={"Values 1": 22, "Values 2": 24}, freeze=None)
    table_sheet(wb, "Blocks 2019–20", HB, brows, fmt={h: "0.000" for h in HB if h.endswith(" r")}, widths={"Block pair": 18})
    return save(wb, 2)


# =================================================================================================================
def esm3():
    """variance components. Source: ../tables/ESM_variance_components.tsv, written and checked by
    make_ESM_variance_components.py from 00_THIN_.../tables/Variance_components_GxE.csv, _WIDE.csv, H2_GxE.csv."""
    t = list(csv.reader(open(SUPP / "tables/ESM_variance_components.tsv", encoding="utf-8"), delimiter="\t"))
    H, rows = t[0], t[1:]
    conv = lambda v: v if v.startswith("<") else float(v)
    rows = [r[:2] + [conv(v) for v in r[2:]] for r in rows]
    h2 = {r[0]: r[-1] for r in rows}
    ok(h2["Starch"] == 0.474 and h2["Fiber"] == 0.248 and h2["Flowering time"] == 0.747 and h2["Tillers"] == "< 0.001",
       "ESM3: H² starch 0.474, fiber 0.248, flowering time 0.747, tillers < 0.001 (Results Ch. 1-2)")
    gxe = {r[0]: r[7] for r in rows if r[1] == "Nutritional"}
    ok(round(min(gxe.values()), 1) == 4.5 and round(max(gxe.values()), 1) == 15.2 and max(gxe, key=gxe.get) == "β-glucan",
       "ESM3: G×E 4.5-15.2% of the nutritional traits, largest for beta-glucan (Results Ch. 1)")
    rows = [[{"Tillers": "Number of tillers"}.get(r[0], r[0])] + r[1:] for r in rows]
    wb = new_workbook()
    cap = ("Variance components and broad-sense heritability (H²) of the four nutritional and four morphological and "
           "phenological traits, from a linear mixed model with genotype, season-by-block and genotype-by-season random "
           "effects fitted to the values centered within season (290 accessions). For each component, its variance and "
           "its percentage of the total phenotypic variance")
    about_sheet(wb, 3, cap, [("Variance components", "One row per trait; traits within each group ordered by decreasing H².")],
                {"Variance components": [
                    ("V_G", "Genetic variance (genotype)."), ("V_E (season & block)", "Environmental variance: season-by-block."),
                    ("V_G×E", "Genotype-by-season interaction variance."), ("V_R", "Residual variance."),
                    ("(%)", "Percentage of the total phenotypic variance."),
                    ("Total variance", "V_G + V_E + V_G×E + V_R, in the squared units of the trait: % of grain for the "
                     "nutritional traits, days² for flowering time and cm² for spike length."),
                    ("H²", "V_G / (V_G + V_E + V_G×E + V_R) (Holland et al. 2003). Number of tillers: V_G = 2.0 × 10⁻⁸ of a total of 76.26, so H² < 0.001.")]})
    table_sheet(wb, "Variance components", H, rows, fmt={h: ("0.00" if h.endswith("(%)") else "0.0000") for h in H[2:-1]} | {"H²": "0.000"},
                widths={"Trait": 18, "Trait group": 14})
    return save(wb, 3)


# =================================================================================================================
def esm5():
    """loci and significant SNPs. Sources: 01_.../00_FINAL_BLUP_3PC/02_loci_FINAL/tables/Table_loci_master.tsv,
    Table_loci_members_full.tsv; 00_snp_level/significant_snps.tsv (beta = minor allele A1, sign fixed 2026-09-22);
    candidate genes per locus from ../tables/ESM_candidate_genes.tsv (joined on trait + lead SNP)."""
    L = rd(LOCI / "02_loci_FINAL/tables/Table_loci_master.tsv")
    S = rd(LOCI / "00_snp_level/significant_snps.tsv")
    M = rd(LOCI / "02_loci_FINAL/tables/Table_loci_members_full.tsv")
    G = rd(SUPP / "tables/ESM_candidate_genes.tsv")
    ok(len(L) == 36 and L.trait.value_counts().to_dict() == {"fiber": 18, "betaglucan": 11, "starch": 5, "protein": 2},
       "ESM5: 36 loci: fiber 18, beta-glucan 11, starch 5, protein 2 (Results Ch. 2)")
    ok(len(S) == 52 and S.trait.value_counts().to_dict() == {"fiber": 24, "betaglucan": 20, "starch": 6, "protein": 2},
       "ESM5: 52 significant SNPs: fiber 24, beta-glucan 20, starch 6, protein 2")
    sig = S.merge(M[["trait", "SNP", "locus_id"]], left_on=["trait", "SNP_id"], right_on=["trait", "SNP"], how="left")
    ok(sig.locus_id.notna().all() and len(sig) == 52, "ESM5: every significant SNP lies in exactly one locus")
    lead = L.merge(S, left_on=["trait", "lead_SNP"], right_on=["trait", "SNP_id"])
    ok(len(lead) == 36 and (lead.lead_beta - lead.beta).abs().max() == 0 and (lead.lead_A1 == lead.A1).all(),
       "ESM5: lead SNP effects identical in the locus and SNP tables (minor allele)")
    ok(rh(statistics.median(L.span_kb), 1) == 93.5 and (L.n_member_SNPs == 1).sum() == 4 and round(L.span_kb.max(), 1) == 1727.8,
       "ESM5: median span 93.5 kb, 4 lead-only loci, largest 1.73 Mb (Results Ch. 2)")
    ok((L.lead_MAF < 0.10).sum() == 25 and round(statistics.median(L.lead_MAF), 3) == 0.068,
       "ESM5: 25 of 36 lead MAF < 0.10, median 0.068")
    eff = {t: (L[L.trait == t].lead_beta.min(), L[L.trait == t].lead_beta.max()) for t in L.trait.unique()}
    ok([round(x, 2) for x in eff["betaglucan"]] == [0.15, 0.36] and [round(x, 2) for x in eff["fiber"]] == [0.12, 0.25]
       and [round(x, 2) for x in eff["starch"]] == [-2.69, -1.31] and [round(x, 2) for x in eff["protein"]] == [-0.40, 0.68],
       "ESM5: allelic effects beta-glucan 0.15-0.36, fiber 0.12-0.25, starch -1.31 to -2.69, protein -0.40 / 0.68")
    top = L.loc[L.groupby("trait").lead_neg_log10_p.idxmax()].set_index("trait")
    ok(top.loc["betaglucan", "lead_SNP"] == "2H:41996626" and round(top.loc["betaglucan", "lead_neg_log10_p"], 2) == 7.60
       and top.loc["fiber", "lead_SNP"] == "7H:573606306" and top.loc["starch", "lead_SNP"] == "7H:573606460"
       and top.loc["protein", "lead_SNP"] == "3H:106623911", "ESM5: strongest peaks 2H:41,996,626 (7.60), 7H fiber/starch, 3H protein")
    ok(set(L[L.locus_id.isin(["betaglucan_L05", "betaglucan_L06"])].lead_SNP) == {"4H:34858354", "4H:34971609"},
       "ESM5: both 4H beta-glucan loci listed (TODO 309)")
    G["trait_raw"] = G.Trait.map({v: k for k, v in {"betaglucan": "β-glucan", "fiber": "Fiber", "starch": "Starch", "protein": "Protein"}.items()})
    ng = G.groupby(["trait_raw", "Lead SNP"]).size().to_dict()
    nsig = sig.groupby("locus_id").size().to_dict()
    ok(sum(ng.values()) == 55 and len(ng) == 15, "ESM5: 55 candidate genes at 15 loci")
    # minor-allele carrier counts (user, ESM0 Q16a, 2026-10-01; note of 2026-09-16): from step 11, T02 (genotype counts of
    # the 36 lead SNPs in the 290 accessions; inbred lines, no heterozygous calls), joined on trait + lead SNP
    T02 = rd(ORIG / "T02_lead_enrichment_tests.tsv")
    car = {(r.trait, r.lead_SNP): (int(r.n_minor), int(r.n_minor + r.n_major)) for _, r in T02.iterrows()}
    ok(len(car) == 36 and all(k in car for k in zip(L.trait, L.lead_SNP)), "ESM5: carrier counts for all 36 lead SNPs (step 11, T02)")
    ok(all(abs(car[(r.trait, r.lead_SNP)][0] / car[(r.trait, r.lead_SNP)][1] - r.lead_MAF) < 0.00051 for _, r in L.iterrows()),
       "ESM5: carriers / accessions called = lead MAF (4 decimals)")
    locus_rows = []
    for _, r in L.iterrows():
        locus_rows.append([r.locus_id, trait_label(r.trait), r.chr, r.lead_SNP, r.lead_bp, r.locus_start, r.locus_end,
                           r.span_kb, r.n_member_SNPs, nsig.get(r.locus_id, 0), r.lead_A1, r.lead_A2, r.lead_MAF,
                           car[(r.trait, r.lead_SNP)][0], car[(r.trait, r.lead_SNP)][1],
                           r.lead_beta, r.lead_SE, r.lead_p, r.lead_neg_log10_p, r.max_internal_gap_kb,
                           r.clump_span_before_severing_kb, r.n_clump_members_before_severing, r.n_members_severed,
                           ng.get((r.trait, r.lead_SNP), 0)])
    HL = ["Locus", "Trait", "Chromosome", "Lead SNP", "Lead SNP position (bp)", "Locus start (bp)", "Locus end (bp)",
          "Span (kb)", "Member SNPs", "Significant SNPs", "Lead minor allele (A1)", "Lead major allele (A2)",
          "Lead MAF", "Minor-allele carriers", "Accessions called", "Lead β (%)", "Lead SE", "Lead P", "Lead −log10(P)", "Largest gap between member SNPs (kb)",
          "Clump span before the gap rule (kb)", "Clump members before the gap rule", "Members removed by the gap rule",
          "Candidate genes"]
    sig = sig.sort_values(["trait", "chr", "position_bp"])
    lead_set = set(zip(L.trait, L.lead_SNP))
    snp_rows = [[trait_label(r.trait), r.locus_id, r.SNP_id, r.chr, r.position_bp, r.A1, r.A2, r.MAF, r.beta, r.SE,
                 r.p_value, r.neg_log10_p, "Y" if (r.trait, r.SNP_id) in lead_set else "N"] for _, r in sig.iterrows()]
    HS = ["Trait", "Locus", "SNP", "Chromosome", "Position (bp)", "Minor allele (A1)", "Major allele (A2)", "MAF",
          "β (%)", "SE", "P", "−log10(P)", "Lead SNP of its locus"]
    wb = new_workbook()
    cap = ("The 36 loci of the genome-wide association analysis of the four nutritional traits, and the 52 SNPs that "
           "exceeded the genome-wide threshold (Bonferroni, α = 0.10; −log10(P) = 6.0454), with their allele frequencies "
           "and effects (EMMAX; 290 accessions; MorexV3 coordinates)")
    about_sheet(wb, 5, cap, [("Loci", "One row per locus: position, span, lead SNP and its statistics, and the number of candidate genes."),
                             ("Significant SNPs", "One row per genome-wide significant SNP, with its locus.")], {
        "Loci": [("Locus", "Locus ID (trait and number); the same IDs are used in Online Resources 6 and 11."),
                 ("Lead SNP", "Chromosome:position (MorexV3) of the most significant SNP of the locus."),
                 ("Locus start, Locus end", "First and last member SNP: the interval searched for candidate genes."),
                 ("Span (kb)", "Locus end − locus start; 0 for the four loci that consist of the lead SNP alone."),
                 ("Member SNPs", "SNPs in linkage disequilibrium with the lead SNP (r² ≥ 0.5, within 2 Mb) in the contiguous block "
                                 "around it, whatever their P value."),
                 ("Significant SNPs", "Member SNPs above the genome-wide threshold."),
                 ("Lead minor / major allele", "A1 = minor allele, A2 = major allele (PLINK)."),
                 ("Lead MAF", "Minor-allele frequency of the lead SNP in the 290 accessions."),
                 ("Minor-allele carriers, Accessions called", "Accessions carrying the minor allele of the lead SNP, and accessions "
                  "with a genotype call at it (the accessions are inbred: no heterozygous calls; MAF = carriers / called)."),
                 ("Lead β (%)", "Effect of the minor allele on the trait BLUP, in % of grain (EMMAX)."),
                 ("Lead SE, Lead P, Lead −log10(P)", "Standard error of β and the association P value (EMMAX)."),
                 ("Largest gap between member SNPs (kb)", "Largest distance between consecutive member SNPs of the locus (at most 50 kb)."),
                 ("Clump span / members before the gap rule", "Extent and size of the LD clump before it was cut at the first gap of more "
                  "than 50 kb on each side of the lead SNP; members removed by the gap rule = the difference."),
                 ("Candidate genes", "High-confidence MorexV3 genes overlapping the locus interval (Online Resource 6).")],
        "Significant SNPs": [("Locus", "Locus that contains the SNP (sheet “Loci”)."),
                             ("Minor / major allele, MAF", "As in sheet “Loci”."),
                             ("β (%)", "Effect of the minor allele, in % of grain."),
                             ("SE, P, −log10(P)", "EMMAX statistics."),
                             ("Lead SNP of its locus", "Y for the 36 lead SNPs.")]})
    table_sheet(wb, "Loci", HL, locus_rows, fmt={"Lead SNP position (bp)": "#,##0", "Locus start (bp)": "#,##0",
                "Locus end (bp)": "#,##0", "Span (kb)": "0.0", "Lead MAF": "0.000", "Lead β (%)": "0.000", "Lead SE": "0.000",
                "Lead P": "0.00E+00", "Lead −log10(P)": "0.00", "Largest gap between member SNPs (kb)": "0.0",
                "Clump span before the gap rule (kb)": "0.0"}, widths={"Locus": 15, "Lead SNP": 14}, freeze="B2")
    table_sheet(wb, "Significant SNPs", HS, snp_rows, fmt={"Position (bp)": "#,##0", "MAF": "0.000", "β (%)": "0.000",
                "SE": "0.000", "P": "0.00E+00", "−log10(P)": "0.00"}, widths={"Locus": 15, "SNP": 14}, freeze="C2")
    return save(wb, 5)


# =================================================================================================================
def esm6():
    """55 candidate genes. Source: ../tables/ESM_candidate_genes.tsv, written and checked by make_ESM_candidate_genes.py
    (sources and checks there; the 3 carried-forward genes reproduce the former Table 2 cell by cell)."""
    G = rd(SUPP / "tables/ESM_candidate_genes.tsv", keep_default_na=False)
    H = list(G.columns)
    ok(len(G) == 55 and (G.Testable == "Y").sum() == 30 and (G["Significant (q ≤ 0.05)"] == "Y").sum() == 23
       and (G["Carried forward"] == "Y").sum() == 3, "ESM6: 55 genes, 30 tested, 23 significant, 3 carried forward (Results Ch. 3)")
    q = G.set_index("Gene ID")
    for gid, qv, e2 in [("HORVU.MOREX.r3.3HG0301300", "5.3e-06", 0.208), ("HORVU.MOREX.r3.5HG0487060", "7.5e-05", 0.222),
                        ("HORVU.MOREX.r3.3HG0301710", "3.8e-04", 0.078), ("HORVU.MOREX.r3.6HG0619810", "5.5e-11", None),
                        ("HORVU.MOREX.r3.3HG0301210", "1.1e-09", None)]:
        ok(f"{float(q.loc[gid, 'BH q']):.1e}" == qv and (e2 is None or round(float(q.loc[gid, 'η²']), 3) == e2),
           f"ESM6: {gid} q = {qv}" + (f", η² = {e2}" if e2 else "") + " (Results Ch. 3)")
    num = {"Gene start (bp)", "Gene end (bp)", "Distance to lead SNP (bp)", "SNPs in gene window", "Accessions grouped",
           "Haplotype groups", "Kruskal–Wallis P", "BH q", "η²", "Difference (SD)", "Identity (%)", "Coverage (%)",
           "nr identity (%)", "nr coverage (%)"}
    rows = [[(float(v) if (h in num and v != "") else (v if v != "" else None)) for h, v in zip(H, r)] for r in G.values.tolist()]
    defs = [("Trait, Locus, Lead SNP", "Trait of the locus, locus ID and its lead SNP (Online Resource 5)."),
            ("Gene ID", "MorexV3 high-confidence gene (Ensembl Plants release 62)."),
            ("Distance to lead SNP (bp)", "From the lead SNP to the nearest edge of the gene; negative when the gene lies at "
                                          "lower coordinates than the lead SNP. No lead SNP lies inside a gene."),
            ("MorexV3 functional description", "Description in the Ensembl Plants release 62 annotation (3 of 55 genes)."),
            ("SNPs in gene window", "SNPs within the gene ± 1 kb, the window of the haplotype analysis."),
            ("Testable, Reason not testable", "Y if crosshap (ε = 0.6, MGmin = 2, minHap = 9) defined at least two haplotype "
                                              "groups; otherwise the reason."),
            ("Accessions grouped, Haplotype groups, Group sizes", "Accessions assigned to a haplotype group (of 290), number of "
                                                                 "groups, and accessions per group, largest first."),
            ("Kruskal–Wallis P", "Kruskal–Wallis test of the trait BLUPs among haplotype groups."),
            ("BH q", "Benjamini–Hochberg correction within each trait (untestable genes excluded)."),
            ("η²", "Kruskal–Wallis η², (H − k + 1) / (n − k)."),
            ("Difference (SD)", "Highest minus lowest haplotype-group mean, in standard deviations of the trait among the tested accessions."),
            ("Significant (q ≤ 0.05)", "Y for the 23 genes with a significant haplotype association."),
            ("Annotation call, Annotation source", "Call taken from the first source that gave one, in the order Swiss-Prot, "
                                                   "InterPro, NCBI nr (significant genes only)."),
            ("Swiss-Prot / InterPro / NCBI nr call", "The call of each source; nr was searched only for proteins without a call "
                                                     "from Swiss-Prot or InterPro."),
            ("Best Swiss-Prot hit, accession, Identity (%), Coverage (%), accepted", "Best DIAMOND BLASTP hit in UniProtKB/Swiss-Prot "
             "2026_03, with its percent identity and query coverage; accepted when e-value ≤ 10⁻⁵, coverage ≥ 50%, identity ≥ 30% "
             "and the hit is a characterized protein."),
            ("Pfam domains", "Pfam signatures among the InterProScan 5.78-109.0 matches."),
            ("Best nr hit, nr identity (%), nr coverage (%)", "Best BLASTP hit in NCBI nr (Viridiplantae)."),
            ("Carried forward", "Y for the three genes presented in the Results (GPAT6, GH17, PHT4;3).")]
    wb = new_workbook()
    cap = ("The 55 candidate genes of the 36 loci, one row per gene: locus and position, whether the gene's haplotypes "
           "could be tested and why not, the Kruskal–Wallis test of the trait BLUPs among haplotype groups, the functional "
           "annotation of the 23 genes with significant haplotype associations, and the three genes carried forward. "
           "Empty cells: step not applied to the gene")
    about_sheet(wb, 6, cap, [("Candidate genes", "One row per gene, in the order of the loci.")], {"Candidate genes": defs})
    table_sheet(wb, "Candidate genes", H, rows, fmt={"Kruskal–Wallis P": "0.00E+00", "BH q": "0.00E+00", "η²": "0.000",
                "Difference (SD)": "0.00", "Identity (%)": "0.0", "Coverage (%)": "0.0", "nr identity (%)": "0.0",
                "nr coverage (%)": "0.0", "Gene start (bp)": "#,##0", "Gene end (bp)": "#,##0", "Distance to lead SNP (bp)": "#,##0"},
                widths={"Gene ID": 27, "MorexV3 functional description": 40, "Reason not testable": 34, "Annotation call": 36,
                        "Swiss-Prot call": 30, "InterPro call": 40, "Best Swiss-Prot hit": 36, "Best nr hit (nr accession)": 40,
                        "Pfam domains": 22, "Lead SNP": 14, "Locus": 14, "Group sizes": 14}, default_width=11, freeze="E2")
    return save(wb, 6)


# =================================================================================================================
def esm7():
    """epsilon choice. Sources: 04_.../loci_LDspan_eps06_V4/Diagnostics/epsilon_choice_supplementary.tsv (summary) and
    epsilon_choice_per_gene.tsv (55 genes x 5 epsilon), written by 01_scripts/diagnostics/parameter_choice_tables.R
    (2026-09-30; genotype-only, no P values)."""
    lines = strip_comment_header(EPS / "epsilon_choice_supplementary.tsv")
    S = pd.DataFrame([l.split("\t") for l in lines[1:]], columns=lines[0].split("\t"))
    P = rd(EPS / "epsilon_choice_per_gene.tsv")
    G = rd(SUPP / "tables/ESM_candidate_genes.tsv", keep_default_na=False)
    s6 = S[S.epsilon == "0.6"].iloc[0]
    ok(s6.genes_testable == "30" and s6.median_assignment_pct == "65" and s6.setting_used == "TRUE",
       "ESM7: at epsilon 0.6, 30 of 55 genes testable, median 65% of accessions assigned (M&M, Results)")
    ok(len(P) == 275 and P.gene.nunique() == 55 and set(P.MGmin) == {2}, "ESM7: 55 genes x 5 epsilon, MGmin 2")
    p6 = P[P.epsilon == 0.6].set_index("gene")
    g = G.set_index("Gene ID")
    ok(set(p6[p6.testable].index) == set(g[g.Testable == "Y"].index)
       and all(int(p6.loc[i, "n_assigned"]) == int(float(g.loc[i, "Accessions grouped"])) and int(p6.loc[i, "n_groups"]) == int(float(g.loc[i, "Haplotype groups"]))
               for i in g[g.Testable == "Y"].index), "ESM7: the 30 rows at epsilon 0.6 reproduce the published run gene by gene")
    srows = [[float(r.epsilon), int(r.genes_testable), int(r.pct_of_candidates), int(r.total_accessions_assigned),
              float(r.median_assignment_pct), float(r.mean_assignment_pct), float(r.median_n_haplotype_groups),
              "Y" if r.setting_used == "TRUE" else "N"] for _, r in S.iterrows()]
    HS_ = ["ε", "Genes testable (of 55)", "Genes testable (%)", "Accessions assigned, summed over genes",
           "Median accessions assigned per testable gene (%)", "Mean accessions assigned per testable gene (%)",
           "Median haplotype groups per testable gene", "Setting used"]
    NOTE = {"ok": "", "no_marker_groups": "no marker group formed", "no_or_too_few_variants": "no or too few SNPs"}
    prow = [[r.gene, trait_label(r.trait), r.n_snps, r.epsilon, "Y" if r.testable else "N",
             (int(r.n_groups) if r.testable else None), (int(r.n_assigned) if r.testable else None),
             (round(100 * r.assign_rate, 1) if r.testable else None), NOTE[r.note]] for _, r in P.sort_values(["gene", "epsilon"]).iterrows()]
    HP = ["Gene ID", "Trait", "SNPs in gene window", "ε", "Testable", "Haplotype groups", "Accessions assigned",
          "Accessions assigned (%)", "Reason not testable"]
    wb = new_workbook()
    cap = ("Choice of the crosshap clustering radius (ε): for five values of ε (MGmin = 2), the number of the 55 candidate "
           "genes whose haplotypes could be tested, the proportion of accessions assigned to haplotype groups and the "
           "number of groups, summarized and per gene. The analysis used ε = 0.6. Phenotypes were not used")
    about_sheet(wb, 7, cap, [("Summary", "One row per ε value."), ("Per gene", "One row per gene and ε value (55 × 5).")], {
        "Summary": [("ε", "DBSCAN clustering radius of crosshap; MGmin = 2 and minHap = 9 throughout."),
                    ("Genes testable", "Genes with at least two haplotype groups."),
                    ("Accessions assigned, summed over genes", "Accessions assigned to a haplotype group, summed over all 55 genes (an untestable gene counts 0)."),
                    ("Median / mean accessions assigned per testable gene (%)", "Over the testable genes, the percentage of the 290 accessions assigned to a group."),
                    ("Median haplotype groups per testable gene", "Over the testable genes."),
                    ("Setting used", "Y for the value used in the analysis.")],
        "Per gene": [("SNPs in gene window", "SNPs within the gene ± 1 kb."),
                     ("Testable", "Y if crosshap defined at least two haplotype groups."),
                     ("Accessions assigned (%)", "Of 290 accessions.")]})
    table_sheet(wb, "Summary", HS_, srows, fmt={"ε": "0.0", "Median accessions assigned per testable gene (%)": "0.0",
                "Mean accessions assigned per testable gene (%)": "0.0", "Median haplotype groups per testable gene": "0.0"}, freeze=None)
    table_sheet(wb, "Per gene", HP, prow, fmt={"ε": "0.0", "Accessions assigned (%)": "0.0"},
                widths={"Gene ID": 27, "Reason not testable": 24}, freeze="B2")
    return save(wb, 7)


# =================================================================================================================
def esm8():
    """canonical beta-glucan genes. Sources: 05_.../09_USED_canonical_betaglucan_gene_check/results/tables/
    canonical_betaglucan_gene_distances.tsv (mapping columns only; its distance columns are stale) and
    canonical_gene_distances_current_loci.tsv (distances to the current 36 loci, 2026-09-23; beta-glucan rows)."""
    A = rd(CANON / "canonical_betaglucan_gene_distances.tsv")
    D = rd(CANON / "canonical_gene_distances_current_loci.tsv")
    mapc = ["canonical_gene", "source_accession", "horvu_id", "chr", "gene_start", "gene_end", "strand", "pct_identity",
            "pct_coverage", "evalue", "reciprocal_best_hit", "tblastn_horvu", "tblastn_vs_blastp_agree"]
    A = A[mapc].drop_duplicates()
    ok(len(A) == 10 and A.canonical_gene.is_unique, "ESM8: ten canonical genes, one mapping each")
    ok((A.pct_identity >= 97).all() and (A.reciprocal_best_hit == "Y").all() and (A.tblastn_vs_blastp_agree == "Y").all(),
       "ESM8: all orthologs >= 97% identity, reciprocal best hits, genome and proteome positions agree (M&M)")
    Db = D[D.trait == "betaglucan"].set_index("horvu_id")
    ok(len(Db) == 10, "ESM8: one beta-glucan distance row per gene")
    near = Db.dropna(subset=["distance_bp"]).sort_values("distance_bp").iloc[0]
    ok(near.canonical_gene == "HvGlb1_EI" and near.distance_Mb == 83.5, "ESM8: nearest = HvGlb1, 83.5 Mb (Discussion)")
    ok(set(Db[Db.distance_bp.isna()].canonical_gene) == {"HvCslF6", "HvGlb2_EII"} and set(Db[Db.distance_bp.isna()].chr) == {"7H"},
       "ESM8: HvCslF6 and HvGlb2 on 7H, no beta-glucan locus there (Discussion)")
    NAME = {"HvGlb1_EI": "HvGlb1", "HvGlb2_EII": "HvGlb2"}
    ENZ = lambda g: ("(1,3;1,4)-β-glucan endohydrolase " + ("EI" if "Glb1" in g else "EII")) if "Glb" in g else \
        "Cellulose synthase-like (CslF/CslH), (1,3;1,4)-β-glucan synthase"
    order = ["HvCslF3", "HvCslF4", "HvCslF6", "HvCslF7", "HvCslF8", "HvCslF9", "HvCslF10", "HvCslH1", "HvGlb1_EI", "HvGlb2_EII"]
    rows = []
    for g in order:
        r = A[A.canonical_gene == g].iloc[0]; d = Db.loc[r.horvu_id]
        db = "NCBI Protein" if "." in r.source_accession else "UniProtKB"
        has = pd.notna(d.distance_bp)
        rows.append([NAME.get(g, g), ENZ(g), f"{r.source_accession} ({db})", r.horvu_id, r.chr, int(r.gene_start),
                     int(r.gene_end), r.strand, float(r.pct_identity), float(r.pct_coverage), float(r.evalue),
                     r.reciprocal_best_hit, r.tblastn_vs_blastp_agree,
                     d.nearest_lead_SNP if has else None, d.nearest_locus_id if has else None,
                     int(d.distance_bp) if has else None, float(d.distance_Mb) if has else None,
                     None if has else "No β-glucan locus on this chromosome"])
    H = ["Gene", "Enzyme", "Reference protein", "MorexV3 gene", "Chromosome", "Start (bp)", "End (bp)", "Strand",
         "Identity (%)", "Query coverage (%)", "E-value", "Reciprocal best hit", "tBLASTn and BLASTP positions agree",
         "Nearest β-glucan lead SNP", "β-glucan locus", "Distance (bp)", "Distance (Mb)", "Note"]
    wb = new_workbook()
    cap = ("The ten canonical (1,3;1,4)-β-glucan synthase and endohydrolase genes located on the MorexV3 genome from their "
           "reference protein sequences, and the distance of each to the nearest β-glucan lead SNP on the same chromosome")
    about_sheet(wb, 8, cap, [("Canonical genes", "One row per gene.")], {"Canonical genes": [
        ("Reference protein", "Accession of the published protein sequence used as the query."),
        ("MorexV3 gene", "MorexV3 high-confidence gene of the best BLASTP hit (highest identity) in the MorexV3 proteome."),
        ("Identity (%), Query coverage (%), E-value", "Of that BLASTP hit (BLAST+ 2.12.0)."),
        ("Reciprocal best hit", "Y if the MorexV3 protein's best hit among the reference proteins is the query."),
        ("tBLASTn and BLASTP positions agree", "Y if the best tBLASTn hit in the MorexV3 genome (e-value ≤ 10⁻¹⁰) lies at the same gene."),
        ("Nearest β-glucan lead SNP, β-glucan locus, Distance", "Distance from the nearest edge of the gene to the nearest β-glucan "
         "lead SNP on the same chromosome (Online Resource 5).")]})
    table_sheet(wb, "Canonical genes", H, rows, fmt={"Start (bp)": "#,##0", "End (bp)": "#,##0", "Identity (%)": "0.00",
                "Query coverage (%)": "0.0", "E-value": "0.0E+00", "Distance (bp)": "#,##0", "Distance (Mb)": "0.0"},
                widths={"Enzyme": 34, "Reference protein": 26, "MorexV3 gene": 27, "Note": 30}, freeze="B2")
    return save(wb, 8)


# =================================================================================================================
def esm10():
    """elite cultivars. Sources: 07_USED_elite_lines_compariosn_to_wild_lines/results/tables/Table_elite_lines.tsv,
    Table_elite_line_screen.tsv, Table_site_overlap.tsv, Table_removed_sites.tsv; 03_01_7H_branch_.../results/tables/
    Table_elite_line_screen.tsv, Table_site_overlap.tsv, Table_removed_sites.tsv."""
    E = rd(ELITE / "Table_elite_lines.tsv"); SC = rd(ELITE / "Table_elite_line_screen.tsv")
    E7 = rd(B7H / "Table_elite_line_screen.tsv").set_index("line_name")
    O = pd.concat([rd(ELITE / "Table_site_overlap.tsv"), rd(B7H / "Table_site_overlap.tsv")], ignore_index=True)
    R3 = rd(ELITE / "Table_removed_sites.tsv"); R7 = rd(B7H / "Table_removed_sites.tsv")
    ok(len(SC) == 136 and SC.passes_threshold.sum() == 91, "ESM10: 136 spring elite lines scored, 91 pass the pooled 85% call rate (M&M)")
    ok(set(SC[SC.configured].accession_name) == set(E.line_name) and len(E) == 5 and SC[SC.configured].passes_threshold.all(),
       "ESM10: the five cultivars are the configured lines and all pass")
    ok(sorted(E.year_of_release) == [2012, 2012, 2012, 2016, 2018], "ESM10: released 2012-2018 (M&M)")
    pooled = dict(zip(E.line_name, (100 * E.call_rate_all).round(1)))
    ok(pooled == {"Avalon": 96.6, "KWS Irina": 94.3, "Odyssey": 94.3, "Laureate": 86.4, "LG Diablo": 98.9},
       "ESM10: pooled call rates 86.4-98.9% (hidden note, M&M)")
    ok(E7.call_rate_pct.to_dict() == {"Avalon": 95.5, "KWS Irina": 95.5, "Odyssey": 81.8, "LG Diablo": 77.3, "Laureate": 54.5},
       "ESM10: GDSL-gene call rates over 22 shared sites (M&M hidden note)")
    o = O.set_index("short_name")
    ok([int(o.loc[g, "n_shared"]) for g in ["GPAT6", "GH17", "PHT4-3", "GDSL"]] == [43, 33, 12, 22]
       and [int(o.loc[g, "n_wild_snps"]) for g in ["GPAT6", "GH17", "PHT4-3", "GDSL"]] == [54, 41, 16, 27]
       and (O.n_alleles_swapped == 0).all(), "ESM10: shared sites 43/33/12/22 of 54/41/16/27 wild SNPs, no swapped alleles")
    GN = {"GPAT6": "GPAT6", "GH17": "GH17", "PHT4-3": "PHT4;3", "GDSL": "GDSL esterase/lipase"}
    cult = [[r.line_name, r.sample_id, int(r.year_of_release), r.breeder, r.country_first_registration, r.annuality,
             round(100 * r.call_rate_GPAT6, 1), round(100 * r["call_rate_GH17"], 1), round(100 * r["call_rate_PHT4-3"], 1),
             round(100 * r.call_rate_all, 1), float(E7.loc[r.line_name, "call_rate_pct"])] for _, r in E.iterrows()]
    HC = ["Cultivar", "BioSamples accession", "Year of release", "Breeder", "Country of first registration", "Growth habit",
          "Call rate GPAT6 (%)", "Call rate GH17 (%)", "Call rate PHT4;3 (%)", "Call rate, three genes pooled (%)",
          "Call rate GDSL esterase/lipase gene (%)"]
    scr = [[r.sample_id, r.accession_name, r.breeder, r.year_of_release, r.annuality, r.called_GPAT6, round(100 * r.call_rate_GPAT6, 1),
            r.called_GH17, round(100 * r.call_rate_GH17, 1), r["called_PHT4-3"], round(100 * r["call_rate_PHT4-3"], 1), r.called_all,
            r.shared_all, round(100 * r.call_rate_all, 1), "Y" if r.passes_threshold else "N", "Y" if r.configured else "N"]
           for _, r in SC.sort_values("call_rate_all", ascending=False).iterrows()]
    HSC = ["BioSamples accession", "Line", "Breeder", "Year of release", "Growth habit", "Sites called GPAT6", "Call rate GPAT6 (%)",
           "Sites called GH17", "Call rate GH17 (%)", "Sites called PHT4;3", "Call rate PHT4;3 (%)", "Sites called, three genes",
           "Compared sites, three genes", "Call rate, pooled (%)", "Pooled call rate ≥ 85%", "Selected"]
    OC = [("gene_id", "Gene ID"), ("trait", "Trait"), ("chr", "Chromosome"), ("win_start", "Window start (bp)"),
          ("win_end", "Window end (bp)"), ("n_wild_snps", "Wild SNPs in window"), ("n_elite_records", "Elite records in window"),
          ("n_shared_positions", "Positions in both call sets"), ("n_alleles_identical", "Same REF and ALT"),
          ("n_alleles_swapped", "REF and ALT swapped"), ("n_alleles_triallelic", "Third allele (excluded)"),
          ("n_shared", "Compared sites"), ("n_wild_only", "Wild-only sites (excluded)"), ("n_elite_only", "Elite-only records (excluded)"),
          ("n_elite_only_indel", "of which indels"), ("n_wild_assigned_accessions", "Wild accessions in haplotype groups"),
          ("n_wild_mono_all", "Wild: monomorphic, all sites"), ("n_wild_mono_kept", "Wild: monomorphic, compared sites"),
          ("n_wild_mono_removed", "Wild: monomorphic, excluded sites"), ("n_elite_lines_shown", "Elite cultivars"),
          ("n_elite_mono_all", "Elite: monomorphic, all records"), ("n_elite_mono_kept", "Elite: monomorphic, compared sites"),
          ("n_elite_mono_removed", "Elite: monomorphic, excluded records")]
    O["trait"] = O.trait.map(lambda x: "Fiber and starch" if x == "fiber+starch" else trait_label(x))
    acc = [[GN[r.short_name]] + [r[c] for c, _ in OC] for _, r in O.iterrows()]
    HO = ["Gene"] + [h for _, h in OC]
    REASON = {"elite_only": "only in the elite call set", "wild_only": "only in the wild call set", "triallelic": "third allele"}
    rem = []
    for _, r in R3.iterrows():
        rem.append([GN[r.short_name], r.gene_id, "Wild" if r.file == "wild" else "Elite", REASON[r.removed_from], r.chr, int(r.pos),
                    r.ref, r.alt, "Y" if r.is_indel else "N", r.state])
    for _, r in R7.iterrows():
        _, _, ref, alt = r.site_key.split(":")
        rem.append([GN[r.short_name], r.gene_id, "Wild" if r.file == "wild" else "Elite", REASON[r.removed_from], r.chr, int(r.pos),
                    ref, alt, "Y" if r.is_indel else "N", r.polymorphic])
    ok(len(rem) == len(R3) + len(R7) == 643 + 114, "ESM10: 757 excluded site records (643 + 114)")
    HR = ["Gene", "Gene ID", "Call set", "Excluded because", "Chromosome", "Position (bp)", "REF", "ALT", "Indel",
          "Variation in the compared genotypes"]
    wb = new_workbook()
    cap = ("Comparison with elite cultivars: the five spring malting cultivars and their call rates at the four presented "
           "genes, the call-rate screen of the 136 spring elite lines of the barley pangenome panel, and the accounting of "
           "the SNP sites compared between the wild and elite MorexV3 call sets (gene ± 1 kb), with every excluded site")
    about_sheet(wb, 10, cap, [
        ("Cultivars", "The five cultivars compared with the wild haplotype groups."),
        ("Screen", "All 136 spring elite lines, by pooled call rate over the compared sites of GPAT6, GH17 and PHT4;3."),
        ("Site accounting", "Per gene: SNP sites in each call set, compared, and excluded."),
        ("Excluded sites", "Every site excluded from the comparison, with the reason.")], {
        "Cultivars": [("BioSamples accession", "EBI BioSamples record of the line (barley pangenome panel; names verified there)."),
                      ("Call rate (%)", "Called genotypes (reference, alternate or heterozygous) / compared sites of the gene. "
                       "Three genes pooled: over the 88 compared sites of GPAT6, GH17 and PHT4;3 (selection threshold 85%). "
                       "GDSL esterase/lipase gene: over its 22 compared sites; no new screen was applied at this gene.")],
        "Screen": [("Sites called, Compared sites", "Called genotypes and compared sites, per gene and pooled."),
                   ("Pooled call rate ≥ 85%", "Y for the 91 lines that pass."),
                   ("Selected", "Y for the five cultivars chosen among the passing lines (widely grown spring malting cultivars).")],
        "Site accounting": [("Window", "Gene ± 1 kb, as in the haplotype analysis."),
                            ("Compared sites", "Sites with the same position and the same REF and ALT alleles in both call sets."),
                            ("Wild-only / elite-only", "Positions with a record in one call set only; excluded."),
                            ("Third allele", "Same position and REF, different ALT; excluded."),
                            ("Monomorphic", "No variation among the compared genotypes: wild = the accessions assigned to haplotype "
                             "groups of the gene; elite = the five cultivars.")],
        "Excluded sites": [("Variation in the compared genotypes", "polymorphic / monomorphic among the genotypes of its call set "
                            "used here (wild: accessions assigned to haplotype groups; elite: the five cultivars); undetermined = "
                            "no genotype called.")]})
    table_sheet(wb, "Cultivars", HC, cult, widths={"Breeder": 18, "Country of first registration": 16}, freeze="B2")
    table_sheet(wb, "Screen", HSC, scr, widths={"Line": 18, "Breeder": 22, "BioSamples accession": 16}, freeze="C2")
    table_sheet(wb, "Site accounting", HO, acc, fmt={"Window start (bp)": "#,##0", "Window end (bp)": "#,##0"},
                widths={"Gene": 20, "Gene ID": 27}, freeze="B2")
    table_sheet(wb, "Excluded sites", HR, rem, fmt={"Position (bp)": "#,##0"}, widths={"Gene": 20, "Gene ID": 27,
                "Excluded because": 24}, freeze="C2")
    return save(wb, 10)


# =================================================================================================================
def esm11():
    """carrier origin. Source: 11_USED_check_low_maf_geographic_distribution_minor_allele/results/tables/T01-T04, T08
    (scripts 02, 03, 05; rebuilt 2026-09-30). User decision ESM0 Q7: desert_OR and the extreme-accession check
    (top2_extreme, n_top2_extreme_carriers, n_top10_extreme_carriers, obs_diff_SD_without_top2,
    pct_diff_kept_without_top2) are not included."""
    T01, T02, T03, T04, T08 = (rd(ORIG / f"{n}.tsv", keep_default_na=False, na_values=["NA"]) for n in
                               ["T01_lead_carriers_by_region", "T02_lead_enrichment_tests", "T03_haplotype_group_members_by_region",
                                "T04_haplotype_group_enrichment_tests", "T08_lead_within_site"])
    t2 = T02.set_index("locus_id"); t4 = T04.set_index(["gene", "group"])
    ok(len(T02) == 36 and len(T08) == 36, "ESM11: 36 lead SNPs")
    ok((t2.loc["fiber_L17", ["n_minor", "n_sites_minor"]].tolist() == [41, 11]) and t2.loc["fiber_L17", "desert_pct_minor"] == 17.1
       and t2.loc["fiber_L17", "desert_pct_major"] == 15.0 and round(t2.loc["fiber_L17", "q_desert"], 2) == 0.81
       and t2.loc["starch_L05", "n_minor"] == 35 and t2.loc["starch_L05", "n_sites_minor"] == 10
       and round(t2.loc["starch_L05", "q_desert"], 2) == 0.82, "ESM11: 7H leads: 41 / 35 carriers, 11 / 10 sites, not desert-enriched (Discussion §2)")
    ok(t4.loc[("GPAT6", "C"), "desert_pct_group"] == 66.7 and t4.loc[("GPAT6", "C"), "desert_pct_other_groups"] == 6.7
       and round(t4.loc[("GPAT6", "C"), "q_desert"], 3) == 0.034 and t4.loc[("GPAT6", "C"), "n_sites"] == 5
       and t4.loc[("GPAT6", "C"), "n_desert_sites"] == 3, "ESM11: GPAT6 group C 66.7% vs 6.7% desert, q 0.034, 5 sites, 3 desert (Discussion §3, §5)")
    ok(t4.loc[("GH17", "C"), "desert_pct_group"] == 93.8 and round(t4.loc[("GH17", "C"), "q_desert"], 3) == 0.043
       and t4.loc[("PHT4;3", "A"), "desert_pct_group"] == 5.1 and round(t4.loc[("PHT4;3", "A"), "q_desert"], 4) == 0.0024
       and t4.loc[("GDSL", "B"), "desert_pct_group"] == 20.6 and round(t4.loc[("GDSL", "B"), "q_desert"], 2) == 0.67,
       "ESM11: GH17 C 93.8% (q 0.043), PHT4;3 A 5.1% (q 0.0024), GDSL B 20.6% (q 0.67) (Discussion)")
    drop2 = ["desert_OR"]; drop8 = ["top2_extreme", "n_top2_extreme_carriers", "n_top10_extreme_carriers",
                                     "obs_diff_SD_without_top2", "pct_diff_kept_without_top2"]
    T02 = T02.drop(columns=drop2); T04 = T04.drop(columns=drop2); T08 = T08.drop(columns=drop8)
    for T in (T01, T02, T08):
        T["trait"] = T.trait.map(trait_label)
    for T in (T03, T04):
        T["trait"] = T.trait.map(lambda x: "Fiber and starch" if x in ("fiber+starch", "fiber; starch") else trait_label(x))
    for T in (T01, T03):
        T["region"] = T.region.map(lambda x: REGION_LABEL.get(x, x))
    RN = {"trait": "Trait", "locus_id": "Locus", "lead_SNP": "Lead SNP", "gene": "Candidate gene", "lead_neg_log10_p": "Lead −log10(P)",
          "MAF": "MAF", "minor_effect": "Minor allele effect", "n_minor": "Minor-allele carriers", "n_major": "Major-allele carriers",
          "n_missing": "Missing calls", "minor_North": "Carriers North", "minor_Coast": "Carriers Coast", "minor_Desert": "Carriers Desert",
          "minor_HZ1": "Carriers HZ1", "minor_HZ2": "Carriers HZ2", "minor_HZ3": "Carriers HZ3", "n_regions_minor": "Regions (carriers)",
          "n_sites_minor": "Sites (carriers)", "n_desert_sites_minor": "Desert sites (carriers)", "top_site": "Site with most carriers",
          "top_site_n": "Carriers at that site", "desert_pct_minor": "Carriers from desert (%)", "desert_pct_major": "Major-allele carriers from desert (%)",
          "desert_direction": "Desert direction", "p_desert": "P desert", "q_desert": "q desert", "p_regions": "P regions", "q_regions": "q regions",
          "region": "Region", "pct_minor_within_region": "Minor-allele carriers within region (%)", "pct_of_minor_carriers": "Share of all minor-allele carriers (%)",
          "pct_of_major_carriers": "Share of all major-allele carriers (%)", "gene_id": "Gene ID", "group": "Haplotype group", "n": "Accessions",
          "n_assigned_gene": "Accessions assigned (gene)", "n_groups_gene": "Haplotype groups (gene)", "mean_BLUP_fiber": "Mean fiber BLUP (%)",
          "mean_BLUP_starch": "Mean starch BLUP (%)", "North": "North", "Coast": "Coast", "Desert": "Desert", "HZ1": "HZ1", "HZ2": "HZ2", "HZ3": "HZ3",
          "n_regions": "Regions", "n_sites": "Sites", "n_desert_sites": "Desert sites", "desert_pct_group": "Group members from desert (%)",
          "desert_pct_other_groups": "Other groups from desert (%)", "p_regions_group": "P regions (group)", "q_regions_group": "q regions (group)",
          "p_regions_gene": "P regions (gene)", "q_regions_gene": "q regions (gene)", "pct_of_group": "Share of the group (%)",
          "obs_diff_SD": "Observed difference (SD)", "raw_sign_agrees_with_beta": "Sign agrees with GWAS β",
          "exp_diff_region_SD": "Expected from region (SD)", "exp_diff_site_SD": "Expected from site (SD)",
          "frac_expected_region": "Expected / observed, region", "frac_expected_site": "Expected / observed, site",
          "n_minor_with_sitemates": "Carriers with site-mates", "pct_carriers_beyond_sitemates": "Carriers beyond site-mates (%)"}
    def sheet_rows(T):
        H = [RN[c] for c in T.columns]
        rows = [[(None if (isinstance(v, float) and v != v) else ("Y" if v is True else "N" if v is False else v)) for v in r] for r in T.values.tolist()]
        return H, rows
    wb = new_workbook()
    cap = ("Geographic origin of the carriers of the lead-SNP alleles of the 36 loci and of the haplotype groups of the four "
           "presented genes: carriers per ecological region and site, the site-permutation tests of desert enrichment and "
           "of the distribution over the six regions, and the comparison of the minor-allele carriers with the major-allele "
           "carriers of their own site")
    about_sheet(wb, 11, cap, [
        ("Lead SNPs – tests", "One row per lead SNP: carriers, their origin, and the two site-permutation tests."),
        ("Lead SNPs – by region", "Minor- and major-allele carriers and missing calls per lead SNP and region."),
        ("Haplotype groups – tests", "One row per haplotype group of GPAT6, GH17, PHT4;3 and the GDSL esterase/lipase gene."),
        ("Haplotype groups – by region", "Members of each haplotype group per region (and the unassigned accessions)."),
        ("Carriers vs site-mates", "One row per lead SNP: minor-allele carriers compared with the major-allele carriers of their own site.")], {
        "Lead SNPs – tests": [("Minor allele effect", "Direction of the minor allele's effect on the trait (raises / lowers)."),
                              ("Carriers <region>", "Minor-allele carriers per region (North, Coast, Desert, HZ1 North–Coast, HZ2 North–Desert, HZ3 Coast–Desert)."),
                              ("Carriers from desert (%)", "Percentage of the minor-allele carriers (and of the major-allele carriers) from the five desert sites; the panel's share is 50 of 290 (17%)."),
                              ("P desert, q desert", "Two-sided site-permutation test of the difference in desert share between minor- and major-allele carriers: "
                               "region labels permuted among the 29 sites (10,000 permutations), the accessions of each site kept together; "
                               "Benjamini–Hochberg q within trait."),
                              ("P regions, q regions", "Site-permutation test of the distribution of the minor-allele carriers over the six regions "
                               "(chi-square statistic, upper tail); Benjamini–Hochberg q within trait.")],
        "Lead SNPs – by region": [("Share of all minor- / major-allele carriers (%)", "Of the carriers of that allele, the percentage from the region.")],
        "Haplotype groups – tests": [("Haplotype groups", "Groups as in Fig. 3 (GPAT6, GH17, PHT4;3; crosshap ε = 0.6, MGmin = 2) and Fig. 4 "
                                      "(GDSL gene; ε = 0.9, MGmin = 3). Unassigned accessions excluded."),
                                     ("Mean fiber / starch BLUP (%)", "Mean genotype BLUP of the group members."),
                                     ("Group members / other groups from desert (%)", "Desert share of the group and of the other assigned accessions of the same gene."),
                                     ("P desert, q desert", "As for the lead SNPs, group vs the other groups of its gene; Benjamini–Hochberg q within gene."),
                                     ("P regions (group), q regions (group)", "Distribution over the six regions, group vs the other groups; q within gene."),
                                     ("P regions (gene), q regions (gene)", "All groups of the gene together (chi-square); q across the four genes.")],
        "Carriers vs site-mates": [("Observed difference (SD)", "Mean BLUP of the minor-allele carriers minus that of the major-allele carriers, in trait standard deviations."),
                                   ("Sign agrees with GWAS β", "N where the raw difference has the opposite sign to the GWAS effect (estimated with three principal components and kinship)."),
                                   ("Expected from site / region (SD)", "The difference expected if each carrier had the mean of the major-allele carriers of its own site "
                                    "(or region); a site without major-allele carriers is replaced by its region."),
                                   ("Expected / observed", "Share of the observed difference predicted by the carriers' origin."),
                                   ("Carriers with site-mates", "Carriers whose site has at least one major-allele carrier."),
                                   ("Carriers beyond site-mates (%)", "Carriers whose BLUP differs from the mean of the major-allele carriers of their own site in the direction of the allelic effect.")]})
    for name, T in [("Lead SNPs – tests", T02), ("Lead SNPs – by region", T01), ("Haplotype groups – tests", T04),
                    ("Haplotype groups – by region", T03), ("Carriers vs site-mates", T08)]:
        H, rows = sheet_rows(T)
        table_sheet(wb, name, H, rows, fmt={h: "0.0000" for h in H if h.startswith(("P ", "q "))},
                    widths={"Locus": 15, "Lead SNP": 14, "Gene ID": 27, "Site with most carriers": 16, "Region": 19}, freeze="C2")
    return save(wb, 11)


if __name__ == "__main__":
    only = set(sys.argv[1:])
    for n, f in [(1, esm1), (2, esm2), (3, esm3), (5, esm5), (6, esm6), (7, esm7), (8, esm8), (10, esm10), (11, esm11)]:
        if not only or str(n) in only:
            print("wrote", f())
    print(f"{len(CHECKS)} checks passed:")
    for c in CHECKS:
        print("  ✓", c)
