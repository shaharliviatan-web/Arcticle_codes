#!/usr/bin/env python3
# make_reading_copy.py - writes a reading copy of a manuscript md: Word-comment markup (TODO comments)
# and hidden <!-- --> notes removed; highlighted anchor text kept as plain text. Added 2026-09-30
# (M&M session). Never edit a reading copy; regenerate it. Run from new_publishing_paper/:
#   python3 build/make_reading_copy.py Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley_reading_copy.md
import re, sys

src, dst = sys.argv[1], sys.argv[2]
t = open(src, encoding="utf-8").read()
t = re.sub(r"<!--.*?-->", "", t, flags=re.S)                                  # hidden notes
t = re.sub(r"\[[^\[\]]*\]\{\.comment-start[^}]*\}", "", t)                    # comment bodies
t = re.sub(r"\[\]\{\.comment-end[^}]*\}", "", t)                              # comment ends
t = re.sub(r"\[\[([^\[\]]*)\]\]\{\.mark\}", r"[\1]", t)                       # placeholders keep brackets
t = re.sub(r"\[([^\[\]]*)\]\{\.mark\}", r"\1", t)                             # anchored text
t = re.sub(r"[ \t]+\n", "\n", t); t = re.sub(r"\n{3,}", "\n\n", t).strip() + "\n"
hdr = f"*Reading copy of `{src}`, comments and hidden notes removed. Do not edit; regenerated from the original.*\n\n"
open(dst, "w", encoding="utf-8").write(hdr + t)
