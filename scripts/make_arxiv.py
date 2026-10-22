#!/usr/bin/env python3
"""Build the arXiv submission tarball (run after `make paper`, which writes manuscript/main.bbl).

The tarball carries main.bbl and the JHEP style files.  The code and data of this package go
into anc/ as ancillary files; following arXiv's rules, no TeX files are placed there, so anc/ holds
the maintained code and data, excluding manuscript/, historical manuscript snapshots, build
products, and the downloaded files of arXiv:1909.08993 (inputs/arXiv-1909.08993, fetched by
inputs/fetch_inputs.sh and checked against inputs/SHA256SUMS).
Writes build/ggAA-threshold-arxiv.tar.gz and lists what it contains.
"""

import pathlib
import re
import sys
import tarfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
MAN = ROOT / "manuscript"
OUT = ROOT / "build" / "ggAA-threshold-arxiv.tar.gz"
SKIP_DIRS = {"manuscript", "history", "manuscript_arXiv_v1", "build", "__pycache__", ".git", "arXiv-1909.08993"}


def manuscript_files():
    """main.tex, main.bbl and every file main.tex reads, found by following \\input and \\includegraphics."""
    files, todo = {'JHEP.bst', 'jheppub.sty', 'main.bbl', 'main.tex', 'references.bib'}, ["main.tex"]
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


# local caches that .gitignore also leaves out
SKIP_FILES = {"results/exact_amplitude_cache.json", "results/li22_store.json"}
# local files at the top level of the working copy (also in .gitignore)
LOCAL_SUFFIXES = {".log", ".xz", ".gz", ".tgz"}
LOCAL_FILES = {"status.sh"}


def ancillary_files():
    """Package files for anc/: no current or historical manuscripts, hidden files or caches.

    Data files are kept whatever their extension (data/born/original/BornSqr.log is FORM output)."""
    out = []
    for path in sorted(ROOT.rglob("*")):
        rel = path.relative_to(ROOT)
        if (path.is_dir() or SKIP_DIRS & set(rel.parts) or any(part.startswith(".") for part in rel.parts)
                or path.suffix in {".pyc", ".tex", ".zip"} or rel.as_posix() in SKIP_FILES):
            continue
        if len(rel.parts) == 1 and (path.suffix in LOCAL_SUFFIXES or rel.name in LOCAL_FILES):
            continue                                   # run logs, monitoring scripts and archives at the top level
        out.append(rel)
    return out


def main():
    if not (MAN / "main.bbl").exists():
        sys.exit("manuscript/main.bbl is missing: run `make paper` first")
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
