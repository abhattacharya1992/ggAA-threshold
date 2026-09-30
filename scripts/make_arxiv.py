#!/usr/bin/env python3
"""Build the arXiv submission tarball (run after `make paper`, which writes manuscript/main.bbl).

arXiv does not run BibTeX, so the tarball carries main.bbl.  The code and data of this package go
into anc/ as ancillary files; following arXiv's rules, no TeX files are placed there, so anc/ holds
everything except manuscript/, and the downloaded files of arXiv:1909.08993 (inputs/arXiv-1909.08993,
fetched by inputs/fetch_inputs.sh and checked against inputs/SHA256SUMS) are left out.
Writes build/ggAA-threshold-arxiv.tar.gz and lists what it contains.
"""

import pathlib
import re
import sys
import tarfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
MAN = ROOT / "manuscript"
OUT = ROOT / "build" / "ggAA-threshold-arxiv.tar.gz"
SKIP_DIRS = {"manuscript", "build", "__pycache__", ".git", "arXiv-1909.08993"}


def manuscript_files():
    """main.tex, main.bbl and every file main.tex reads, found by following \\input and \\includegraphics."""
    files, todo = {"main.tex", "main.bbl"}, ["main.tex"]
    while todo:
        text = re.sub(r"(?<!\\)%.*", "", (MAN / todo.pop()).read_text())
        for name in re.findall(r"\\input\{([^}]*)\}", text):
            name = name if name.endswith(".tex") else name + ".tex"
            if name not in files:
                files.add(name)
                todo.append(name)
        for name in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]*)\}", text):
            files.add(name if pathlib.Path(name).suffix else name + ".pdf")
    return sorted(files)


def ancillary_files():
    out = []
    for path in sorted(ROOT.rglob("*")):
        rel = path.relative_to(ROOT)
        if path.is_dir() or SKIP_DIRS & set(rel.parts) or path.suffix in {".pyc", ".tex", ".aux", ".log"}:
            continue
        out.append(rel)
    return out


def main():
    if not (MAN / "main.bbl").exists():
        sys.exit("manuscript/main.bbl is missing: run `make paper` first (arXiv does not run BibTeX)")
    paper = manuscript_files()
    missing = [f for f in paper if not (MAN / f).exists()]
    if missing:
        sys.exit(f"files read by main.tex are missing: {missing}")
    anc = ancillary_files()
    OUT.parent.mkdir(exist_ok=True)
    with tarfile.open(OUT, "w:gz") as tar:
        for f in paper:
            tar.add(MAN / f, arcname=f)
        for f in anc:
            tar.add(ROOT / f, arcname=str(pathlib.Path("anc") / f))
    tex_in_anc = [f for f in anc if f.suffix == ".tex"]
    assert not tex_in_anc, tex_in_anc
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size / 1e6:.1f} MB)")
    print(f"  manuscript: {len(paper)} files ({', '.join(f for f in paper if '/' not in f)}, "
          f"{sum('/' in f for f in paper)} in sections/, tables/, figures/)")
    print(f"  anc/: {len(anc)} files, no TeX files")
    print("Upload the tarball at arxiv.org; the anc/ directory becomes the ancillary files.")


if __name__ == "__main__":
    main()
