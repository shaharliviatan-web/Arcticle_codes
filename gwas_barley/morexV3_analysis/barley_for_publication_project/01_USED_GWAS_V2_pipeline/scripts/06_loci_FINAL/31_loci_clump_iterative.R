#!/usr/bin/env Rscript
# 31_loci_clump_iterative.R
# Locus definition, end to end: LD clumping, the gap rule, and iteration until no
# Bonferroni-significant SNP is left outside a locus.
#
# RULES (locked 2026-09-09):
#   leads        --clump-p1 9.008e-07  (Bonferroni alpha=0.10 over 111,017 pruned SNPs)
#   members      --clump-p2 1          membership on LD alone, no p-value condition
#   LD           --clump-r2 0.5
#   reach        --clump-kb 2000       +/- 2 Mb, 4 Mb maximum span
#   contiguity   severed at the first gap > 50 kb walking outward from the lead
#   ITERATION    any Bonferroni-significant SNP left outside every locus after
#                severing becomes a candidate lead in the next pass; repeat until
#                none are left.
#
# WHY THE ITERATION EXISTS. --clump is winner-take-all: each SNP joins exactly one
# clump. A significant SNP that is not a lead, and that sits beyond a gap from the
# lead which owns it, is severed and then belongs to nothing -- it would be absent
# from the tables and from the painted Manhattans despite passing the genome-wide
# threshold. The iteration promotes any such SNP to a lead of its own.
#
# HOW A SNP IS BARRED FROM LEADING IN LATER PASSES: PLINK has no "only these SNPs
# may lead" option, so each pass writes an assoc file in which every ALREADY
# ASSIGNED significant SNP is set to P = 1. It then cannot pass --clump-p1, while
# --clump-p2 1 still lets it join a clump as a member. Orphans keep their real
# p-value and so become the only eligible leads.
#
# Termination is guaranteed: an orphan promoted to lead always anchors its own
# walk and is therefore always assigned (worst case as a lead-only locus).
#
# The script asserts at the end that EVERY Bonferroni-significant SNP sits inside
# a locus, so this class of silent loss cannot regress.
#
# Outputs -> 02_loci_FINAL/tables/{loci_summary.tsv, loci_members.tsv}
#            02_loci_FINAL/plink/pass<k>/
# Created 2026-09-09.
Sys.setenv(TMPDIR="/mnt/data/shahar/.tmp")
suppressPackageStartupMessages(library(data.table))
PIPE  <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/barley_for_publication_project/01_USED_GWAS_V2_pipeline"
SET   <- file.path(PIPE,"results","00_FINAL_BLUP_3PC")
ASSOC <- file.path(SET,"01_assoc")
OUT   <- file.path(SET,"02_loci_FINAL"); TAB <- file.path(OUT,"tables"); PLK <- file.path(OUT,"plink")
BFILE <- file.path(PIPE,"intermediates","morexV3_290")
PLINK <- "/usr/local/bin/plink"
dir.create(TAB, showWarnings=FALSE, recursive=TRUE)

P1 <- as.numeric(Sys.getenv("P1","9.008e-07")); P2 <- as.numeric(Sys.getenv("P2","1"))
R2 <- as.numeric(Sys.getenv("R2","0.5"));       KB <- as.numeric(Sys.getenv("KB","2000"))
MAX_GAP  <- as.numeric(Sys.getenv("MAX_GAP","50000"))
MAX_PASS <- as.integer(Sys.getenv("MAX_PASS","8"))
TRAITS <- c("betaglucan","fiber","protein","starch")
NEED <- 7110997L
cat(sprintf("[31] p1=%s p2=%s r2=%.2f kb=%.0f (+/-%.0f kb) | gap<=%.0f kb | iterate to zero orphans\n",
            format(P1,scientific=TRUE), P2, R2, KB, KB, MAX_GAP/1000))

sig <- fread(file.path(SET,"00_snp_level","significant_snps.tsv"))[, .(trait, SNP=SNP_id, bp=position_bp)]
cat(sprintf("[31] %d Bonferroni-significant SNPs must all end up inside a locus\n", nrow(sig)))

# gap-rule walk: the connected run of clump members around the lead
walk_block <- function(lead_bp, member_bp) {
  bp <- sort(unique(c(lead_bp, member_bp)))
  ai <- which(bp == lead_bp)[1]
  hi <- ai; while (hi < length(bp) && (bp[hi+1]-bp[hi]) <= MAX_GAP) hi <- hi + 1L
  lo <- ai; while (lo > 1        && (bp[lo]-bp[lo-1]) <= MAX_GAP) lo <- lo - 1L
  list(bp = bp[lo:hi], n_all = length(bp), raw_span = max(bp)-min(bp))
}

all_loci <- list(); all_mem <- list(); assigned <- character(0)
for (pass in seq_len(MAX_PASS)) {
  pdir <- file.path(PLK, sprintf("pass%d", pass)); dir.create(pdir, showWarnings=FALSE, recursive=TRUE)
  made <- 0L
  for (tr in TRAITS) {
    a0 <- file.path(ASSOC, paste0(tr, ".assoc"))
    stopifnot(file.exists(a0), length(readLines(a0, n=1)) == 1)
    A <- fread(a0)
    if (nrow(A) + 1L != NEED) stop(sprintf("[31] FAIL: %s has %d rows, expected %d", a0, nrow(A)+1L, NEED))
    if (length(assigned)) A[SNP %chin% assigned, P := 1]        # bar from leading, keep as member
    af <- file.path(pdir, paste0(tr, ".assoc")); fwrite(A, af, sep="\t")
    system2(PLINK, c("--bfile",BFILE,"--allow-extra-chr","--clump",af,
                     "--clump-p1",P1,"--clump-p2",P2,"--clump-r2",R2,"--clump-kb",KB,
                     "--clump-field","P","--clump-snp-field","SNP",
                     "--out",file.path(pdir,tr),"--silent"), stdout=FALSE, stderr=FALSE)
    cf <- file.path(pdir, paste0(tr,".clumped")); if (!file.exists(cf)) next
    cl <- fread(cf, fill=TRUE); cl <- cl[!is.na(BP) & SNP != ""]
    for (i in seq_len(nrow(cl))) {
      sp <- cl$SP2[i]
      ms <- if (is.na(sp) || sp %chin% c("NONE","")) character(0) else trimws(gsub("\\(1\\)","",strsplit(sp,",")[[1]]))
      ms <- ms[grepl(":", ms, fixed=TRUE)]
      mbp <- suppressWarnings(as.numeric(sub(".*:","",ms))); mbp <- mbp[is.finite(mbp)]
      w <- walk_block(as.numeric(cl$BP[i]), mbp)
      g <- if (length(w$bp) > 1) diff(w$bp) else 0
      id <- sprintf("%s_P%d_%02d", tr, pass, i)
      all_loci[[length(all_loci)+1L]] <- data.table(
        trait=tr, locus_id=id, pass=pass, lead_SNP=cl$SNP[i], chr=as.character(cl$CHR[i]),
        lead_bp=as.numeric(cl$BP[i]), lead_p=cl$P[i],
        start=min(w$bp), end=max(w$bp), span_kb=(max(w$bp)-min(w$bp))/1000,
        n_members=length(w$bp), max_internal_gap_kb=max(g)/1000,
        raw_clump_span_kb=w$raw_span/1000, n_clump_members=w$n_all,
        n_severed=w$n_all-length(w$bp))
      all_mem[[length(all_mem)+1L]] <- data.table(
        trait=tr, locus_id=id, pass=pass,
        SNP=sprintf("%s:%d", as.character(cl$CHR[i]), w$bp),   # chr already carries the H
        BP=w$bp, is_index=(w$bp == as.numeric(cl$BP[i])))
      made <- made + 1L
    }
  }
  L <- rbindlist(all_loci)
  sig[, inloc := FALSE]
  for (i in seq_len(nrow(L)))
    sig[trait==L$trait[i] & bp>=L$start[i] & bp<=L$end[i] &
        sub(":.*","",SNP)==L$chr[i], inloc := TRUE]
  orph <- sig[inloc == FALSE]
  cat(sprintf("[31] pass %d: +%d loci (%d total) | significant SNPs covered %d/%d | orphans %d\n",
              pass, made, nrow(L), sum(sig$inloc), nrow(sig), nrow(orph)))
  if (!nrow(orph)) break
  assigned <- unique(c(assigned, sig[inloc == TRUE, SNP]))
  if (pass == MAX_PASS) warning("[31] hit MAX_PASS with orphans remaining")
}

S <- rbindlist(all_loci); M <- rbindlist(all_mem)
setorder(S, trait, chr, lead_bp)
S[, new_id := sprintf("%s_L%02d", trait, rowid(trait))]
map <- S[, .(trait, locus_id, new_id)]
M <- merge(M, map, by=c("trait","locus_id"), all.x=TRUE); stopifnot(!any(is.na(M$new_id)))
M[, locus_id := new_id][, new_id := NULL]; S[, locus_id := new_id][, new_id := NULL]
setorder(S, trait, chr, lead_bp); setorder(M, trait, locus_id, BP)
stopifnot(setequal(S$locus_id, unique(M$locus_id)))

# HARD INVARIANT: every Bonferroni-significant SNP is inside a locus
sig[, inloc := FALSE]
for (i in seq_len(nrow(S)))
  sig[trait==S$trait[i] & bp>=S$start[i] & bp<=S$end[i] &
      sub(":.*","",SNP)==S$chr[i], inloc := TRUE]
if (any(!sig$inloc)) { print(sig[inloc==FALSE]); stop("[31] FAIL: significant SNPs outside every locus") }

# single source of truth: every downstream script reads its parameter strings from
# here, so figure titles / footnotes / tables can never state a stale value again
fwrite(data.table(parameter=c("clump_p1","clump_p2","clump_r2","clump_kb",
                              "max_span_kb","gap_rule_kb","n_passes","n_loci",
                              "n_significant_snps_covered"),
                  value=c(format(P1,scientific=TRUE), P2, R2, KB, 2*KB,
                          MAX_GAP/1000, max(S$pass), nrow(S), nrow(sig))),
       file.path(TAB,"locus_definition_params.tsv"), sep="\t")
fwrite(S, file.path(TAB,"loci_summary.tsv"), sep="\t")
fwrite(M, file.path(TAB,"loci_members.tsv"), sep="\t")
cat(sprintf("\n[31] %d loci over %d pass(es) | median span %.1f kb | mean %.1f kb | ALL %d significant SNPs covered\n",
            nrow(S), max(S$pass), median(S$span_kb), mean(S$span_kb), nrow(sig)))
print(S[, .(n_loci=.N, median_span_kb=round(median(span_kb),1), max_span_kb=round(max(span_kb),1),
            from_later_passes=sum(pass>1)), by=trait])
cat(sprintf("[31] OK -> %s\n", TAB))
