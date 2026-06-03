#!/usr/bin/env Rscript
# 02_build_table.R
# Step 02 of 07_USED_BLASTP_genes_with_no_annotation_left.
# Build the deliverable from the remote nr BLASTP output: one row per residual
# gene (5 rows), best hit = lowest e-value (ties -> highest bitscore). Organism is
# parsed from the trailing "[Organism]" of the nr title (sscinames is often blank
# on -remote). characterized_flag uses the SAME regex as step 05.
#
# Outputs:
#   results/tables/fdr_residual_nr.tsv   (5 rows)
#   results/tables/residual_nr_note.txt  (one-line named-vs-uncharacterized summary)

Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

base <- "/mnt/data/shahar/gwas_barley/morexV3_analysis/07_USED_BLASTP_genes_with_no_annotation_left"
faa  <- file.path(base, "inputs/residual_5.faa")
raw  <- file.path(base, "intermediates/residual_5_vs_nr.tsv")
out  <- file.path(base, "results/tables/fdr_residual_nr.tsv")
note <- file.path(base, "results/tables/residual_nr_note.txt")

UNCHAR_RE <- "predicted|uncharacterized|hypothetical|DUF[0-9]|unknown function|putative uncharacterized"

strip_iso <- function(x) sub("\\.[0-9]+$", "", x)

# residual gene list (from the 5-seq FASTA), preserve as gene_id
ids <- sub("^>", "", sub("\\s.*$", "", readLines(faa)[grep("^>", readLines(faa))]))
genes <- data.frame(seqid = ids, gene_id = strip_iso(ids), stringsAsFactors = FALSE)

bcols <- c("qseqid","sseqid","stitle","pident","qcovs","evalue","bitscore","staxids","sscinames")
if (file.exists(raw) && file.info(raw)$size > 0) {
  bl <- read.delim(raw, header = FALSE, quote = "", comment.char = "",
                   stringsAsFactors = FALSE, fill = TRUE)
  # pad/trim to expected width
  while (ncol(bl) < length(bcols)) bl[[paste0("V", ncol(bl)+1)]] <- ""
  bl <- bl[, seq_along(bcols)]; names(bl) <- bcols
  bl$pident   <- as.numeric(bl$pident)
  bl$qcovs    <- as.numeric(bl$qcovs)
  bl$evalue   <- as.numeric(bl$evalue)
  bl$bitscore <- as.numeric(bl$bitscore)
  bl$gene_id  <- strip_iso(bl$qseqid)
} else {
  bl <- data.frame(qseqid=character(), gene_id=character())
}

best_for <- function(g) {
  sub <- bl[bl$gene_id == g, , drop = FALSE]
  if (nrow(sub) == 0) {
    return(data.frame(nr_title=NA, nr_accession=NA, nr_organism=NA,
                      pident=NA, qcovs=NA, evalue=NA, bitscore=NA,
                      characterized_flag=NA, status="NO_NR_HIT",
                      stringsAsFactors = FALSE))
  }
  sub <- sub[order(sub$evalue, -sub$bitscore), , drop = FALSE]
  b <- sub[1, ]
  # In NCBI BLAST outfmt 6 the stitle column is ALREADY just the description
  # (the subject id is the separate sseqid column), so use it as-is. Do NOT strip
  # the first token or the leading descriptive word ("uncharacterized",
  # "antifreeze", ...) would be lost and the characterized regex defeated.
  title <- b$stitle
  # accession from sseqid (e.g. "ref|XP_...|" or "gi|..|ref|XP_..|"); take last non-empty
  parts <- strsplit(b$sseqid, "\\|")[[1]]; parts <- parts[nzchar(parts)]
  acc <- if (length(parts)) parts[length(parts)] else b$sseqid
  if (grepl("^(ref|gb|emb|dbj|sp|pdb|pir|tpg|tpe|prf)$", acc) && length(parts) >= 2)
    acc <- parts[length(parts)-1]   # guard if last token was a db tag
  # organism: prefer sscinames, but blastp writes literal "N/A" when no local
  # taxdb is installed (as here on -remote), so treat that as missing and fall
  # back to the last bracketed [Organism] in the nr title.
  sci <- b$sscinames
  valid_sci <- nzchar(sci) && !grepl("^\\s*(N/A)(\\s*;\\s*N/A)*\\s*$", sci)
  org <- if (valid_sci) sci else {
    m <- regmatches(title, gregexpr("\\[[^][]+\\]", title))[[1]]
    if (length(m)) gsub("^\\[|\\]$", "", m[length(m)]) else NA
  }
  # title without the trailing organism bracket, for the name
  name <- trimws(sub("\\s*\\[[^][]+\\]\\s*$", "", title))
  charf <- !grepl(UNCHAR_RE, name, ignore.case = TRUE)
  status <- if (b$evalue <= 1e-5 && b$qcovs >= 50 && b$pident >= 30 && charf)
              "CONFIDENT_NR_HIT" else "WEAK_OR_UNCHARACTERIZED_NR_HIT"
  data.frame(nr_title=name, nr_accession=acc, nr_organism=org,
             pident=b$pident, qcovs=b$qcovs, evalue=b$evalue, bitscore=b$bitscore,
             characterized_flag=charf, status=status, stringsAsFactors = FALSE)
}

res <- cbind(gene_id = genes$gene_id,
             do.call(rbind, lapply(genes$gene_id, best_for)))
res <- res[order(res$gene_id), ]
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# one-line note: which got a real name vs stayed uncharacterized / no hit
named <- res$gene_id[res$status == "CONFIDENT_NR_HIT"]
unc   <- res$gene_id[res$status == "WEAK_OR_UNCHARACTERIZED_NR_HIT"]
none  <- res$gene_id[res$status == "NO_NR_HIT"]
line <- sprintf("nr Viridiplantae rescue of 5 residual genes: %d got a named plant homolog (%s); %d uncharacterized/weak (%s); %d no hit (%s).",
                length(named), ifelse(length(named), paste(sub("HORVU.MOREX.r3.","",named),collapse=","), "-"),
                length(unc),   ifelse(length(unc),   paste(sub("HORVU.MOREX.r3.","",unc),collapse=","), "-"),
                length(none),  ifelse(length(none),  paste(sub("HORVU.MOREX.r3.","",none),collapse=","), "-"))
writeLines(line, note)

cat("=== Step 02: build nr residual table ===\n")
cat("rows:", nrow(res), " (expect 5)\n\n")
print(res[, c("gene_id","nr_title","pident","qcovs","evalue","characterized_flag","status")])
cat("\nstatus counts:\n"); print(table(res$status))
cat("\nNOTE: ", line, "\n", sep = "")
cat("\nWrote:\n  ", out, "\n  ", note, "\n")
