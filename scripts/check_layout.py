#!/usr/bin/env python3
"""Text report on float placement in the compiled manuscript (run after `make paper`).

Reads manuscript/main.aux and main.log written by the real RevTeX build, and the source files in
the order main.tex inputs them.  For every figure and table it prints the page on which it was set
and the pages of the section that first cites it, and flags

  LATE     the float appears after the page on which the next section begins and more than one
           page after its own section begins (its section is the later of the one that first
           cites it and the one in which it is defined),
  EARLY    the float appears more than one page before its section begins,
  STUCK    RevTeX/LaTeX reported floats that could not be placed.  RevTeX prints "A float is
           stuck" whenever a single-column float is still waiting when a new page starts; the float
           is then placed normally in the next column.  Without the class option floatfix these
           notices are information only (the table above shows where each float went); with
           floatfix they trigger emergency float processing and count as problems,
  OVERFULL / STRETCHED COLUMN  a column that did not come out at the page height,
  FLUSHED  the \\FloatBarrier before the bibliography had to flush pending floats onto a page of
           their own (the references then start on a later page than the end of the appendices).

With pdftotext (poppler) available, the report also lists
  SPARSE   a column with few words away from the first page and the references: possibly a column
           holding only floats, which is worth a look in the PDF (not counted as a problem).

The report is plain text so that it can be pasted into a conversation instead of a PDF.

  python3 scripts/check_layout.py              # manuscript/main.{aux,log}
  python3 scripts/check_layout.py --job NAME   # manuscript/NAME.{aux,log}
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
MAN = ROOT / "manuscript"


def expand(path):
    """Source text with \\input files expanded in place and comments removed."""
    text = re.sub(r"(?<!\\)%.*", "", path.read_text())

    def repl(m):
        name = m.group(1)
        target = MAN / (name if name.endswith(".tex") else name + ".tex")
        return expand(target) if target.exists() else ""
    return re.sub(r"\\input\{([^}]*)\}", repl, text)


def main():
    job = sys.argv[sys.argv.index("--job") + 1] if "--job" in sys.argv else "main"
    aux, log = MAN / f"{job}.aux", MAN / f"{job}.log"
    version = re.search(r'__version__ = "([^"]+)"', (ROOT / "src/python/ggaa/__init__.py").read_text()).group(1)
    print(f"package version {version}; files {aux.name}, {log.name}")
    if "inlineblock" in (MAN / "main.tex").read_text():
        print("WARNING: main.tex still inputs inlineblock.tex -- this is an old copy of the package.")
    if not aux.exists():
        sys.exit("no aux file: run `make paper` first")

    pages = {k: int(p) for k, p in re.findall(r"\\newlabel\{([^}]*)\}\{\{[^}]*\}\{(\d+)\}", aux.read_text())}
    src = expand(MAN / "main.tex")

    # sections in order (top-level \section with a label, and the appendices), with their start pages
    sections = []
    for m in re.finditer(r"\\section\*?\{((?:[^{}]|\{[^{}]*\})*)\}\s*\\label\{([^}]*)\}", src):
        sections.append((m.start(), m.group(2), m.group(1)[:40]))
    order = [s[1] for s in sections] + ["start-of-references"]
    start = {lab: pages.get(lab) for lab in order}

    # the section that first cites each float
    floats = [(m.group(1), m.group(2), m.start())
              for m in re.finditer(r"\\begin\{(figure|table)\*?\}.*?\\label\{([^}]*)\}", src, re.S)]
    first_cite = {}
    for m in re.finditer(r"\\ref\{([^}]*)\}", src):
        for key in m.group(1).split(","):
            if key not in first_cite:
                home = [s for s in sections if s[0] < m.start()]
                first_cite[key] = home[-1][1] if home else order[0]

    problems = 0
    print(f"\n{'float':34s} {'no.':>4s} {'page':>4s}  first cited in (starts p.)   next section (p.)  verdict")
    for kind, key, pos in floats:
        cited = first_cite.get(key)
        defined = ([sec for sec in sections if sec[0] < pos] or [sections[0]])[-1][1]
        # the float belongs to the later of the section that first cites it and the one that defines it
        home = None if cited is None else max(cited, defined, key=order.index)
        page = pages.get(key)
        if home is None or page is None:
            print(f"{key:34s}  -- not cited or not in aux")
            continue
        nxt = order[order.index(home) + 1]
        h0, n0 = start.get(home), start.get(nxt)
        verdict = "ok"
        if n0 is not None and h0 is not None and page > max(n0, h0 + 1):
            verdict = "LATE"
        elif h0 is not None and page < h0 - 1:
            verdict = "EARLY"
        problems += verdict != "ok"
        num = re.search(r"\\newlabel\{%s\}\{\{([^}]*)\}" % re.escape(key), aux.read_text()).group(1)
        print(f"{key:34s} {num:>4s} {page:4d}  {home:24s} ({h0})   {nxt:18s} ({n0})  {verdict}")

    text = log.read_text(errors="ignore") if log.exists() else ""
    stuck = re.findall(r"(A float is stuck[^\n]*\n?[^\n]*|Float too large[^\n]*|Deferred float stuck[^\n]*)", text)
    options = re.search(r"\\documentclass\[([^]]*)\]", (MAN / "main.tex").read_text())
    floatfix = options is not None and "floatfix" in options.group(1).split(",")
    too_large = [s for s in stuck if "too large" in s]
    for s in too_large:
        print("STUCK:", " ".join(s.split()))
    problems += len(too_large)
    notices = len(stuck) - len(too_large)
    if notices and floatfix:
        print(f"STUCK: {notices} stuck-float notice(s) with floatfix on: remove floatfix from \\documentclass")
        problems += notices
    elif notices:
        print(f"note: {notices} stuck-float notice(s) from RevTeX; without floatfix these are information "
              "only (a float waiting at a page break), see the pages above")
    # columns that did not come out at the page height: overfull (text runs into the margin) or
    # stretched (\flushbottom had to open up gaps); pages are the ones being built at the time
    for m in re.finditer(r"(Overfull|Underfull) \\vbox \((?:badness (\d+)|([\d.]+)pt too high)\) has occurred while \\output is active", text):
        shipped = re.findall(r"\[(\d+)(?=[\s\]{<])", text[:m.start()])
        page = int(shipped[-1]) + 1 if shipped else 1
        if m.group(1) == "Overfull":
            print(f"OVERFULL COLUMN on p. {page}: {m.group(3)}pt too high")
            problems += 1
        elif int(m.group(2)) >= 10000:
            print(f"STRETCHED COLUMN on p. {page} (gaps opened to fill the column)")
            problems += 1
    end_app, refs = pages.get("end-of-appendices"), pages.get("start-of-references")
    if end_app is not None and refs is not None:
        flushed = refs > end_app
        print(f"\nappendices end on p. {end_app}; references start on p. {refs}"
              + ("  FLUSHED: floats were put on a page of their own before the references" if flushed else ""))
        problems += flushed
    sparse = sparse_columns(MAN / f"{job}.pdf")
    for page, side, words in sparse:
        print(f"SPARSE: p. {page} {side} column has {words} words -- check it is not a column of floats only")
    print(f"\n{problems} problem(s)" if problems else "\nno problems found")


def sparse_columns(pdf, threshold=150):
    """(page, side, words) for columns with fewer than `threshold` words, from pdftotext -bbox.

    The first page (title and abstract) and the pages from the start of the references on are skipped.
    Returns an empty list if pdftotext is not installed or the PDF is missing."""
    import shutil
    import subprocess
    if not pdf.exists() or shutil.which("pdftotext") is None:
        return []
    html = subprocess.run(["pdftotext", "-bbox", str(pdf), "-"], capture_output=True, text=True).stdout
    aux = pdf.with_suffix(".aux")
    labels = {k: int(p) for k, p in re.findall(r"\\newlabel\{([^}]*)\}\{\{[^}]*\}\{(\d+)\}", aux.read_text())} \
        if aux.exists() else {}
    last = labels.get("start-of-references", 10**6)
    out = []
    for n, page in enumerate(re.findall(r"<page (.*?)</page>", html, re.S), start=1):
        if n == 1 or n >= last:
            continue
        width = float(re.search(r'width="([\d.]+)"', page).group(1))
        xs = [float(x) for x in re.findall(r'<word xMin="([\d.]+)"', page)]
        for side, count in (("left", sum(x < width / 2 for x in xs)), ("right", sum(x >= width / 2 for x in xs))):
            if count < threshold:
                out.append((n, side, count))
    return out


if __name__ == "__main__":
    main()
