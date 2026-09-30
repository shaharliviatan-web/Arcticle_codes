#!/usr/bin/env python3
"""make_ESM_candidate_genes.py -- Online Resource (ESM number fixed later): all 55 candidate genes, one row per
gene, through every step of the pipeline: locus and position, testability, haplotype test, annotation, and whether
the gene was carried forward. Replaces Table 2 of the main text.

WHY (S. Hübner's comment on Table 2, 2026-09-30: "put all genes in one supmat"; user decision: no tables remain in
the main text). This one Online Resource replaces the three planned separately: the 55-gene table (M&M TODOs 222,
224), the 23 annotated significant genes (M&M TODO 227, Results TODO 303, "Online Resource 2") and Table 2.
Everything is put in and nothing is filtered; a cell is empty where a step does not apply to the gene.

NOTHING IS RE-COMPUTED. Values are copied from the project outputs (read-only), joined on gene_id (and lead_SNP),
never on locus_id (10_USED_Paper_writing/CLAUDE.md, hard rule 3):
  03_USED_candidate_genes_around_leading_snps/results/tables/candidate_genes.tsv       55 genes: trait, locus, lead SNP,
                                                                                         position, strand, distance, MorexV3 description
  04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv     30 tested genes: SNPs in the
                                                                                         window, accessions grouped, groups, KW P, BH q, eta2, SD difference
  04_.../Stats/genes_not_tested.tsv                                                      25 untestable genes: SNPs in the window, reason
  05_USED_gene_annotation_analysis/08_USED_annotation_master/results/tables/Table_significant_genes_paper.tsv
                                                                                         23 significant genes: final call, source, call of each source
  05_.../results/tables/Table_significant_genes_annotated.tsv                            best Swiss-Prot hit (name, accession,
                                                                                         identity, coverage, accepted), Pfam IDs, nr best hit
  06_USED_genes_selected_to_present/results/tables/Table_2_genes_carried_forward.tsv     the 3 genes carried forward
                                                                                         (written by 06_.../scripts/01_make_table2_paper.R)

Checks (the script stops if any fails): 55 unique genes; 30 tested + 25 untestable = 55, no overlap; every tested and
untestable gene has the same trait, lead SNP and locus as in candidate_genes.tsv; 23 significant (BH q <= 0.05) in 12
loci; per trait as in the Results (genes / loci with genes: fiber 32 / 6, starch 10 / 3, beta-glucan 7 / 4, protein
6 / 2; tested / significant: fiber 19 / 12, beta-glucan 6 / 6, starch 5 / 5, protein 0 / 0); 3 genes with a MorexV3
functional description; the 23 annotated genes are exactly the 23 significant ones; and the three genes carried
forward reproduce every cell of Table 2 (Table_2_genes_carried_forward.tsv, the source of the manuscript's Table 2).

NOT IN THIS TABLE (user decision 2026-09-30, S0_QUESTIONS.md Q6): the categories of the biological-relevance screen
(trait_candidate_strength: strong / plausible / unlikely / no annotation; candidate_rationale) must not appear in any
supplementary file. The script never reads those columns, and it stops if any of those words or column names appears
in a cell or header of the output.

Distance to the lead SNP: from the lead SNP to the nearest edge of the gene, as in candidate_genes.tsv; negative when
the gene lies at lower coordinates than the lead SNP (no lead SNP lies inside a gene).

Outputs (this folder's ../tables/): ESM_candidate_genes.xlsx (sheet 1 the table, sheet 2 notes and column
definitions), ESM_candidate_genes.tsv (same values, for checking and diffs).
The ESM title block TAG asks for is added when the Online Resources are numbered and assembled, not here.

Created 2026-09-30 (S. Hübner's structural comments; hand-over build/STRUCTURE_CHANGES_2026-10.md).
Run: python3 make_ESM_candidate_genes.py   (seconds; openpyxl)
"""
import csv, math, pathlib
from collections import Counter
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, Border, Side
from openpyxl.utils import get_column_letter

ROOT = pathlib.Path("/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project")
CAND = ROOT / "03_USED_candidate_genes_around_leading_snps/results/tables/candidate_genes.tsv"
STATS = ROOT / "04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Stats"
ANN = ROOT / "05_USED_gene_annotation_analysis/08_USED_annotation_master/results/tables"
T2 = ROOT / "06_USED_genes_selected_to_present/results/tables/Table_2_genes_carried_forward.tsv"
OUT = pathlib.Path(__file__).resolve().parent.parent / "tables"
OUT.mkdir(exist_ok=True)


def rd(path):
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def na(x):
    return x is None or x in ("", "NA")


cand = rd(CAND)
tested = {r["gene_id"]: r for r in rd(STATS / "gene_results.tsv")}
untested = {r["gene_id"]: r for r in rd(STATS / "genes_not_tested.tsv")}
paper = {r["gene_id"]: r for r in rd(ANN / "Table_significant_genes_paper.tsv")}
annot = {r["gene_id"]: r for r in rd(ANN / "Table_significant_genes_annotated.tsv")}
t2 = {r["gene_id"]: r for r in rd(T2)}

ids = [r["gene_id"] for r in cand]
assert len(ids) == 55 == len(set(ids)), "expected 55 unique candidate genes"
assert len(tested) == 30 and len(untested) == 25 and not set(tested) & set(untested)
assert set(tested) | set(untested) == set(ids)
for r in cand:                                   # same trait, lead SNP and locus in every source (join on gene_id)
    s = tested.get(r["gene_id"]) or untested[r["gene_id"]]
    assert (s["trait"], s["lead_SNP"], s["locus_id"]) == (r["trait"], r["lead_SNP"], r["locus_id"]), r["gene_id"]
sig = {g for g, r in tested.items() if r["significant_fdr"] == "TRUE"}
assert all((float(tested[g]["fdr_q"]) <= 0.05) == (g in sig) for g in tested)
assert len(sig) == 23 and len({tested[g]["locus_id"] for g in sig}) == 12
assert set(paper) == set(annot) == sig, "the annotated genes must be exactly the 23 significant ones"
for g in sig:
    assert paper[g]["lead_SNP"] == tested[g]["lead_SNP"] == annot[g]["lead_SNP"], g
assert set(t2) <= sig and len(t2) == 3

# per-trait counts reported in Results Ch. 2 and Ch. 3
n_genes = Counter(r["trait"] for r in cand)
n_loci = {t: len({r["locus_id"] for r in cand if r["trait"] == t}) for t in n_genes}
assert n_genes == {"fiber": 32, "starch": 10, "betaglucan": 7, "protein": 6}, n_genes
assert n_loci == {"fiber": 6, "starch": 3, "betaglucan": 4, "protein": 2}, n_loci
n_t = Counter(tested[g]["trait"] for g in tested); n_s = Counter(tested[g]["trait"] for g in sig)
assert n_t == {"fiber": 19, "betaglucan": 6, "starch": 5} and n_s == {"fiber": 12, "betaglucan": 6, "starch": 5}
assert sum(r["has_annotation"] == "TRUE" for r in cand) == 3

TRAIT = {"betaglucan": "β-glucan", "fiber": "Fiber", "protein": "Protein", "starch": "Starch"}
REASON = {"no_variants_in_window": "No SNP in the gene window",
          "fewer_than_MGmin_variants_after_raw_imputed_intersection":
              "Fewer than two SNPs shared by the observed and imputed genotypes",
          "no_marker_groups_at_epsilon_0.6": "No marker group formed (ε = 0.6, MGmin = 2)"}
SOURCE = {"check1_swissprot": "Swiss-Prot", "check2_interpro": "InterPro", "none": "None"}


def reason(s):
    if s.startswith("crosshap_error"):
        return "crosshap stopped with an internal error"
    return REASON[s]


HEAD = ["Trait", "Locus", "Lead SNP", "Gene ID", "Chromosome", "Gene start (bp)", "Gene end (bp)", "Strand",
        "Distance to lead SNP (bp)", "MorexV3 functional description", "SNPs in gene window", "Testable",
        "Reason not testable", "Accessions grouped", "Haplotype groups", "Group sizes", "Kruskal–Wallis P", "BH q",
        "η²", "Difference (SD)", "Significant (q ≤ 0.05)", "Annotation call", "Annotation source",
        "Swiss-Prot call", "InterPro call", "NCBI nr call", "Best Swiss-Prot hit", "Swiss-Prot accession",
        "Identity (%)", "Coverage (%)", "Swiss-Prot hit accepted", "Pfam domains", "Best nr hit (nr accession)",
        "nr identity (%)", "nr coverage (%)", "Carried forward"]
TRAIT_ORDER = ["betaglucan", "fiber", "protein", "starch"]

rows = []
for r in sorted(cand, key=lambda r: (TRAIT_ORDER.index(r["trait"]), r["chr"], int(r["gene_start"]))):
    g = r["gene_id"]
    t, u, p, a = tested.get(g), untested.get(g), paper.get(g), annot.get(g)
    row = {h: None for h in HEAD}
    row.update({"Trait": TRAIT[r["trait"]], "Locus": r["locus_id"], "Lead SNP": r["lead_SNP"], "Gene ID": g,
                "Chromosome": r["chr"], "Gene start (bp)": int(r["gene_start"]), "Gene end (bp)": int(r["gene_end"]),
                "Strand": r["strand"], "Distance to lead SNP (bp)": int(r["dist_to_lead_bp"]),
                "MorexV3 functional description": None if na(r["description"]) else r["description"],
                "Testable": "Y" if t else "N",
                "Carried forward": "Y" if g in t2 else "N"})
    if u:
        row["SNPs in gene window"] = int(u["n_snps_window"])
        row["Reason not testable"] = reason(u["reason"])
    if t:
        assert int(t["dist_to_lead_bp"]) == int(r["dist_to_lead_bp"]), g
        sizes = [int(x) for x in t["group_sizes"].split("|")]
        assert int(t["n_groups"]) == len(sizes) >= 2 and sizes == sorted(sizes, reverse=True) \
            and sum(sizes) == int(t["n_ind"]), g
        row.update({"SNPs in gene window": int(t["n_snps_window"]), "Accessions grouped": int(t["n_ind"]),
                    "Haplotype groups": int(t["n_groups"]), "Group sizes": t["group_sizes"].replace("|", ", "),
                    "Kruskal–Wallis P": float(t["kw_p_raw"]), "BH q": float(t["fdr_q"]),
                    "η²": float(t["eta_squared"]), "Difference (SD)": float(t["delta_top_bottom_sd"]),
                    "Significant (q ≤ 0.05)": "Y" if g in sig else "N"})
    if p:
        assert p["final_call"] == a["final_call"], g
        row.update({"Annotation call": p["final_call"], "Annotation source": SOURCE[p["annotation_source"]],
                    "Swiss-Prot call": None if na(p["call_check1_swissprot"]) else p["call_check1_swissprot"],
                    "InterPro call": None if na(p["call_check2_interpro"]) else p["call_check2_interpro"],
                    "NCBI nr call": None if na(p["call_check3_nr"]) else p["call_check3_nr"]})
        if not na(a["sp_name"]):
            row.update({"Best Swiss-Prot hit": a["sp_name"], "Swiss-Prot accession": a["sp_accession"],
                        "Identity (%)": float(a["sp_pident"]), "Coverage (%)": float(a["sp_qcovhsp"]),
                        "Swiss-Prot hit accepted": "Y" if a["sp_status"] == "CONFIDENT_SWISSPROT_HIT" else "N"})
        if not na(a["pfam_ids"]):
            row["Pfam domains"] = a["pfam_ids"].replace(";", ", ")
        if not na(a["nr_title"]):
            row.update({"Best nr hit (nr accession)": f'{a["nr_title"]} ({a["nr_accession"]})',
                        "nr identity (%)": float(a["nr_pident"]), "nr coverage (%)": float(a["nr_qcovs"])})
    rows.append(row)

# ---- the three genes carried forward must reproduce Table 2 cell by cell --------------------------------
SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def sci(x):                                      # as 01_make_table2_paper.R: 5.26e-06 -> "5.3 × 10⁻⁶"
    e = math.floor(math.log10(x)); m = round(x / 10 ** e, 1)
    if m >= 10:
        m, e = m / 10, e + 1
    return f"{m:.1f} × 10{str(e).translate(SUP)}"


for row in rows:
    if row["Carried forward"] != "Y":
        continue
    c = t2[row["Gene ID"]]
    got = {"position_cell": f'{row["Gene start (bp)"] / 1e6:.3f}–{row["Gene end (bp)"] / 1e6:.3f}',
           "lead_cell": f'{int(row["Lead SNP"].split(":")[1]) / 1e6:.3f}',
           "dist_cell": f'{round(abs(row["Distance to lead SNP (bp)"]) / 1000)} kb',
           "grouped_cell": f'{row["Accessions grouped"]} ({row["Haplotype groups"]})',
           "q_cell": sci(row["BH q"]), "eta_cell": f'{row["η²"]:.3f}', "delta_cell": f'{row["Difference (SD)"]:.2f}',
           "ident_cell": f'{row["Identity (%)"]:.1f} / {row["Coverage (%)"]:.1f}'}
    for k, v in got.items():
        assert v == c[k], f'{row["Gene ID"]} {k}: {v} != Table 2 {c[k]}'
    assert (row["Chromosome"], row["Trait"].lower()) == (c["chr"], c["trait"])

# ---- the screen categories must not leak into the Online Resource (user, Q6) ------------------------------
BANNED = {"strong", "plausible", "unlikely", "no_annotation", "no annotation", "trait_candidate_strength",
          "candidate_rationale", "tier"}
for row in rows:
    for h in HEAD:
        v = row[h]
        assert h.lower() not in BANNED, h
        assert not (isinstance(v, str) and v.strip().lower() in BANNED), (row["Gene ID"], h, v)
assert not any(w in h.lower() for h in HEAD for w in ("candidacy", "strength", "rationale", "screen"))

# ---- TSV -----------------------------------------------------------------------------------------------
with open(OUT / "ESM_candidate_genes.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(HEAD)
    w.writerows([["" if row[h] is None else row[h] for h in HEAD] for row in rows])

# ---- XLSX ----------------------------------------------------------------------------------------------
CAPTION = ("The 55 candidate genes of the 36 loci, one row per gene: the locus and lead SNP, the position of the gene, "
           "whether its haplotypes could be tested and why not, the Kruskal–Wallis test of the trait BLUPs among "
           "haplotype groups (Benjamini–Hochberg q within each trait), the functional annotation of the 23 genes with "
           "significant haplotype associations, and the three genes carried forward. Empty cells: step not applied "
           "to the gene.")
DEFS = [
    ("Locus", "Locus ID as in the 36-locus table (Online Resource N)."),
    ("Lead SNP", "Chromosome:position (MorexV3) of the lead SNP of the locus."),
    ("Gene ID", "MorexV3 high-confidence gene (Ensembl Plants release 62)."),
    ("Distance to lead SNP (bp)", "From the lead SNP to the nearest edge of the gene; negative when the gene lies at "
                                  "lower coordinates than the lead SNP. No lead SNP lies inside a gene."),
    ("MorexV3 functional description", "Description in the Ensembl Plants release 62 annotation (3 of 55 genes)."),
    ("SNPs in gene window", "SNPs within the gene ± 1 kb, the window of the haplotype analysis."),
    ("Testable", "Y if crosshap (ε = 0.6, MGmin = 2, minHap = 9) defined at least two haplotype groups."),
    ("Accessions grouped", "Accessions assigned to a haplotype group (of 290)."),
    ("Group sizes", "Accessions per haplotype group, largest first."),
    ("Kruskal–Wallis P", "Kruskal–Wallis test of the trait BLUPs among haplotype groups."),
    ("BH q", "Benjamini–Hochberg correction of the Kruskal–Wallis P values within each trait (untestable genes excluded)."),
    ("η²", "Kruskal–Wallis η², (H − k + 1) / (n − k)."),
    ("Difference (SD)", "Difference between the highest and the lowest haplotype-group mean, in standard deviations "
                        "of the trait among the tested accessions."),
    ("Annotation call / source", "Call taken from the first source that gave one, in the order Swiss-Prot, InterPro, "
                                 "NCBI nr (significant genes only)."),
    ("Swiss-Prot / InterPro / NCBI nr call", "The call of each source. nr was searched only for proteins without a call "
                                             "from Swiss-Prot or InterPro."),
    ("Best Swiss-Prot hit, identity, coverage", "Best DIAMOND BLASTP hit in UniProtKB/Swiss-Prot 2026_03, with its "
                                                "percent identity and query coverage. Accepted: e-value ≤ 10⁻⁵, "
                                                "coverage ≥ 50%, identity ≥ 30% and a characterized protein."),
    ("Pfam domains", "Pfam signatures among the InterProScan 5.78-109.0 matches."),
    ("Best nr hit", "Best BLASTP hit in NCBI nr (Viridiplantae), with its percent identity and query coverage."),
    ("Carried forward", "Y for the three genes selected for their biological relevance and presented in the Results."),
]
SOURCES = ("03_USED_candidate_genes_around_leading_snps/results/tables/candidate_genes.tsv; "
           "04_USED_haplotype_analysis_crosshap/04_runs/loci_LDspan_eps06_V4/Stats/gene_results.tsv and genes_not_tested.tsv; "
           "05_USED_gene_annotation_analysis/08_USED_annotation_master/results/tables/Table_significant_genes_paper.tsv and "
           "Table_significant_genes_annotated.tsv; 06_USED_genes_selected_to_present/results/tables/"
           "Table_2_genes_carried_forward.tsv")

wb = Workbook()
ws = wb.active
ws.title = "Candidate genes"
bold, thin = Font(bold=True), Side(style="thin")
ws.append(HEAD)
for c in ws[1]:
    c.font = bold
    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    c.border = Border(bottom=thin)
for row in rows:
    ws.append([row[h] for h in HEAD])
FMT = {"Kruskal–Wallis P": "0.00E+00", "BH q": "0.00E+00", "η²": "0.000", "Difference (SD)": "0.00",
       "Identity (%)": "0.0", "Coverage (%)": "0.0", "nr identity (%)": "0.0", "nr coverage (%)": "0.0",
       "Gene start (bp)": "#,##0", "Gene end (bp)": "#,##0", "Distance to lead SNP (bp)": "#,##0"}
for j, h in enumerate(HEAD, 1):
    if h in FMT:
        for i in range(2, len(rows) + 2):
            ws.cell(i, j).number_format = FMT[h]
WIDTH = {"Gene ID": 27, "MorexV3 functional description": 40, "Reason not testable": 34, "Annotation call": 36,
         "Swiss-Prot call": 30, "InterPro call": 40, "NCBI nr call": 14, "Best Swiss-Prot hit": 36,
         "Best nr hit (nr accession)": 40, "Pfam domains": 22, "Lead SNP": 14, "Locus": 14, "Group sizes": 14}
for j, h in enumerate(HEAD, 1):
    ws.column_dimensions[get_column_letter(j)].width = WIDTH.get(h, 11)
ws.row_dimensions[1].height = 45
ws.freeze_panes = "E2"

wn = wb.create_sheet("Notes")
for k, v in [("Caption (draft)", CAPTION)] + DEFS + [
        ("Sources", SOURCES),
        ("Built by", "10_USED_Paper_writing/new_publishing_paper/supplementary/scripts/make_ESM_candidate_genes.py"),
        ("To add", "TAG ESM title block (article title, journal, authors, corresponding author's affiliation and "
                   "e-mail) when the Online Resources are numbered; the locus-table Online Resource number.")]:
    wn.append([k, v])
    wn.cell(wn.max_row, 1).font = bold
    wn.cell(wn.max_row, 1).alignment = Alignment(vertical="top", wrap_text=True)
    wn.cell(wn.max_row, 2).alignment = Alignment(vertical="top", wrap_text=True)
wn.column_dimensions["A"].width = 30
wn.column_dimensions["B"].width = 110
wb.save(OUT / "ESM_candidate_genes.xlsx")

n_test = sum(r["Testable"] == "Y" for r in rows)
print(f"wrote {OUT / 'ESM_candidate_genes.xlsx'} and .tsv: {len(rows)} genes, {n_test} tested, "
      f"{sum(r['Significant (q ≤ 0.05)'] == 'Y' for r in rows)} significant, "
      f"{sum(r['Carried forward'] == 'Y' for r in rows)} carried forward (match Table 2); all checks passed")
