#!/usr/bin/env python3
"""Proposal / approval tool for the review items of build/REVIEW_2026-09-27_Results_Discussion.md.

Each review item has an edit file, <ID>_edits.json: a list of edits, each with
  where : where the paragraph is, e.g. "Results Ch. 2 ¶3" (as in the review's paragraph key)
  why   : one line on the reason
  old   : the exact text in the working file (must occur exactly once)
  new   : the replacement, with {-deleted text-} and {+inserted text+} marking the change

  python3 review_tool.py preview <ID>   -> writes <ID>_PROPOSED.md (open it with the markdown preview:
                                           deletions red and struck through, insertions green).
                                           The working file is NOT touched.
  python3 review_tool.py apply <ID>     -> applies the approved edits to the working file and
                                           checks that nothing else in the file changed.

Written 2026-09-27 (user request: review changes in the md preview, approve or discard each item).
2026-09-28 (user decision): the working file is build/UPDATED_Results_Discussion_Conclusions.md, converted from the
user's Word-edited docx. It carries the user's own Word comments and reply threads besides the TODO comments, so
readable() now parses every comment span instead of assuming one TODO comment per anchor.
2026-09-30 (user decision): the working file is new_publishing_paper/Methods_Results_Discussion_Conclusions.md, the
Materials and methods merged with the Results, Discussion and Conclusions; build/UPDATED_Results_Discussion_Conclusions.md
is superseded. `apply` also regenerates the reading copy (Methods_Results_Discussion_Conclusions_reading_copy.md).
2026-10-01 (user decision): the title, Key message, Abstract, Keywords and Introduction were added to the working file (item AB1),
and the file was renamed after the title (WORK / READING_COPY below).
"""
import json, re, sys, pathlib

HERE = pathlib.Path(__file__).resolve().parent
WORK = HERE.parent.parent / "Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley.md"  # renamed 2026-10-01 (user: the file is named after the title); was "Methods_Results_Discussion_Conclusions.md"   # since 2026-09-30 (M&M + R + D + C); was build/UPDATED_Results_Discussion_Conclusions.md (2026-09-28 to 2026-09-30), now superseded
READING_COPY = HERE.parent.parent / "Genome-wide_association_and_haplotype_analysis_identify_candidate_genes_for_grain_nutritional_quality_in_wild_barley_reading_copy.md"
MAKE_COPY = HERE.parent / "make_reading_copy.py"
DEL = '<del style="color:#c0392b;background:#fdecea">{}</del>'
INS = '<ins style="color:#1e8449;background:#e9f7ef;text-decoration:none">{}</ins>'
MARK = re.compile(r'\{-(.*?)-\}|\{\+(.*?)\+\}', re.S)
NOTE = re.compile(r'\{@(.*?)@\}', re.S)   # preview-only label, e.g. {@→ M&M@}; never written to the manuscript
LABEL = '<span style="color:#1f4e9c;font-weight:bold;font-style:italic"> ({}) </span>'


def clean(new):            # the text that goes into the manuscript on approval
    new = NOTE.sub('', new)
    return MARK.sub(lambda m: m.group(2) if m.group(2) is not None else '', new)


def render(new):           # the text shown in the preview
    new = NOTE.sub(lambda m: LABEL.format(m.group(1)), new)
    return MARK.sub(lambda m: INS.format(m.group(2)) if m.group(2) is not None else DEL.format(m.group(1)), new)


def _drop_comment_starts(par):
    """remove every [..]{.comment-start ..} span (bracket-matched, escapes respected); return text + {id: author}"""
    authors, out, last = {}, [], 0
    for m in re.finditer(r'\]\{\.comment-start id="(\d+)" author="([^"]*)"[^}]*\}', par):
        if m.start() < last:
            continue
        i, depth = m.start(), 0
        while i >= 0:
            if par[i] in '[]' and not (i > 0 and par[i - 1] == '\\'):
                depth += 1 if par[i] == ']' else -1
                if depth == 0:
                    break
            i -= 1
        out.append(par[last:i]); last = m.end()
        authors[m.group(1)] = m.group(2)
    out.append(par[last:])
    return ''.join(out), authors


def readable(par):         # hide what the Word file hides; mark where each Word comment sits
    par = re.sub(r'<!--.*?-->', '', par, flags=re.S)
    par, authors = _drop_comment_starts(par)
    par = re.sub(r'\[\[(.*?)\]\]\{\.mark\}', r'\1', par, flags=re.S)     # TODO placeholder [[x]]{.mark}
    par = re.sub(r'\[([^\[\]]*)\]\{\.mark\}', r'\1', par, flags=re.S)     # user-edited highlight [x]{.mark}
    par = re.sub(r'\[\]\{\.comment-end id="(\d+)"\}',
                 lambda m: f'<sup>({"TODO" if authors.get(m.group(1)) == "TODO" else "comment"} {m.group(1)})</sup>', par)
    return re.sub(r'[ \t]{2,}', ' ', par).strip()


def load(item, with_questions=False):
    # the edit file is either a list of edits, or {"questions": [...], "edits": [...]}.
    # Questions are shown at the top of the preview: the user reads only the preview file,
    # so every decision Claude needs goes there, never only in chat (user, 2026-09-27).
    spec = json.loads((HERE / f"{item}_edits.json").read_text(encoding="utf-8"))
    edits, questions = (spec, []) if isinstance(spec, list) else (spec["edits"], spec.get("questions", []))
    notes = [] if isinstance(spec, list) else spec.get("notes", [])
    text = WORK.read_text(encoding="utf-8")
    for e in edits:
        n = text.count(e["old"])
        if n != 1:
            sys.exit(f'ERROR [{e["where"]}]: the old text occurs {n} times in {WORK.name} '
                     f'(expected 1). The file may have changed; the proposal must be rebuilt.')
    return (edits, text, questions, notes) if with_questions else (edits, text)


def preview(item):
    edits, text, questions, notes = load(item, with_questions=True)
    paras = re.split(r'\n\s*\n', text)
    out = [f"# Review item {item}: proposed changes",
           "",
           f"Proposal only. `{WORK.name}` is unchanged until you approve.",
           "<del style=\"color:#c0392b;background:#fdecea\">Red, struck through</del> = to be deleted · "
           "<ins style=\"color:#1e8449;background:#e9f7ef;text-decoration:none\">green</ins> = to be inserted · "
           "<span style=\"color:#1f4e9c;font-weight:bold;font-style:italic\">(blue)</span> = where removed text goes; preview only, never in the manuscript. "
           "Each paragraph is shown in full, as it will read in Word (source comments hidden, TODO placeholders "
           "marked).",
           ""]
    if notes:
        out += ["> **ℹ️ Notes**", ">"] + sum(([f"> {n}", ">"] for n in notes), []) + [""]
    if questions:
        out += ["> **❓ Questions for you — please answer these together with your decision**", ">"]
        for qi, q in enumerate(questions, 1):
            out += [f"> **Q{qi}.** {q}", ">"]
        out += ["", "---", ""]
    for i, e in enumerate(edits, 1):
        par = next(p for p in paras if e["old"] in p)
        shown = readable(par.replace(e["old"], render(e["new"])))
        after = readable(par.replace(e["old"], clean(e["new"]))) or "*(paragraph removed; nothing is left of it in Word)*"
        # a change that only adds new text (no deletions) is shown once: the marked and the clean
        # version would be identical apart from the green colour (user, 2026-10-01, item AB1)
        pure_insert = "{-" not in e["new"] and e["new"].startswith("{+") and e["new"].endswith(e["old"])
        if pure_insert:
            out += [f"## Change {i} of {len(edits)} · {e['where']}", "", f"*Why:* {e['why']}", "",
                    "**New text (all of it is added; shown once):**", "", "> " + after, "", "---", ""]
            continue
        out += [f"## Change {i} of {len(edits)} · {e['where']}", "", f"*Why:* {e['why']}", "",
                "**With the changes marked:**", "", shown, "",
                "**Reads after approval:**", "", "> " + after, "", "---", ""]
    out += ["**Your decision:** approve · discard · or tell me what to change."]
    dst = HERE / f"{item}_PROPOSED.md"
    dst.write_text("\n".join(out) + "\n", encoding="utf-8")
    print("written", dst)


def apply(item):
    edits, text = load(item)
    new_text = text
    for e in edits:
        new_text = new_text.replace(e["old"], clean(e["new"]))
    # exact check: the result must equal the untouched text between the edits + the approved
    # replacements, with no two edits overlapping
    spans = sorted((text.index(e["old"]), text.index(e["old"]) + len(e["old"]), clean(e["new"])) for e in edits)
    if any(spans[i][1] > spans[i + 1][0] for i in range(len(spans) - 1)):
        sys.exit("ERROR: two edits overlap; nothing written.")
    rebuilt, pos = [], 0
    for s, t, repl in spans:
        rebuilt += [text[pos:s], repl]
        pos = t
    rebuilt.append(text[pos:])
    if "".join(rebuilt) != new_text:
        sys.exit("ERROR: the result differs from the approved edits alone; nothing written.")
    WORK.write_text(new_text, encoding="utf-8")
    # keep the comment-free reading copy in step with the working file (2026-09-30)
    import subprocess
    subprocess.run([sys.executable, str(MAKE_COPY), str(WORK), str(READING_COPY)], check=True)
    print(f"reading copy regenerated: {READING_COPY.name}")
    # text removed "→ M&M" is saved word for word for the Materials and methods session
    moved = [(e["where"], readable(m.group(1)), m.group(2)) for e in edits
             for m in re.finditer(r'\{-((?:(?!-\}).)*)-\}\{@(→ M&M[^@]*)@\}', e["new"], re.S)]
    if moved:
        ho = HERE.parent / "MM_HANDOVER_from_review.md"
        head = "" if ho.exists() else (
            f"# Text moved out of {WORK.name} for Materials and methods\n\n"
            "Written by `proposed/review_tool.py apply` when an approved review item moves text to M&M. Each entry is the\n"
            "removed text word for word, with where it came from. The M&M session writes it in its own words, and\n"
            "checks every number against the project files first.\n")
        with ho.open("a", encoding="utf-8") as f:
            f.write(head)
            for where, txt, label in moved:
                f.write(f"\n## Item {item} · from {where}\n\n*{label}*\n\n> {txt}\n")
        print(f"saved {len(moved)} passage(s) for M&M to {ho.name}")
    print(f"applied {len(edits)} edit(s) of item {item} to {WORK.name}; "
          f"{len(text.split())} -> {len(new_text.split())} words (incl. hidden comments)")


if __name__ == "__main__":
    {"preview": preview, "apply": apply}[sys.argv[1]](sys.argv[2])
