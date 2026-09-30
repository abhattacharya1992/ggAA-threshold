r"""Look-alike layout check (NOT RevTeX).

The sandbox used for development has no RevTeX and no newtx fonts, so this script builds the manuscript
with a two-column article-class stand-in (same page geometry and float parameters).  Page breaks will
not match a real RevTeX build; the script is a stress test of float placement, not a preview.

  python3 scripts/check_layout_standin.py            # shifts 0..1320 pt in 120 pt steps
  python3 scripts/check_layout_standin.py 0 240 480  # chosen shifts

For each shift K the extra space is inserted after Sec. II in 100-pt pieces (a single large \vspace
would be discarded at a column break), which moves every later page and column break.  Checks:
Tables VII, VIII and Fig. 6 lie within the appendices, no float-only page before the references
(the \FloatBarrier did not flush pending floats), all main-text floats precede the appendices.
"""
import re, subprocess, sys, shutil, pathlib
MAN = pathlib.Path(__file__).resolve().parents[1] / 'manuscript'
PRE = r'''\documentclass[10pt,twocolumn]{article}
\usepackage[paperwidth=8.5in,paperheight=11in,textwidth=510pt,textheight=682pt,columnsep=18pt,top=0.8in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{mathptmx}
\usepackage{amsmath,amssymb,bm}
\usepackage{graphicx}
\input{diagram_setup.tex}
\usepackage{booktabs}
\usepackage{microtype}
\usepackage{placeins}
\usepackage{xcolor}
\usepackage[hidelinks]{hyperref}
\newenvironment{ruledtabular}{\par\hrule\smallskip}{\smallskip\hrule}
\newenvironment{acknowledgments}{\section*{Acknowledgments}}{}
\newcommand{\affiliation}[1]{\gdef\theaffil{#1}}
\hyphenpenalty=900
\exhyphenpenalty=900
\emergencystretch=1.25em
\flushbottom
\setcounter{topnumber}{4}
\setcounter{dbltopnumber}{2}
\renewcommand{\topfraction}{0.95}
\renewcommand{\dbltopfraction}{0.95}
\renewcommand{\textfraction}{0.05}
\renewcommand{\floatpagefraction}{0.82}
\renewcommand{\dblfloatpagefraction}{0.82}
'''
def build(K, name):
    s = (MAN / 'main.tex').read_text()
    body = s[s.index(r'\input{metadata.tex}'):]
    body = body.replace(r'\thanks{Contact author: \href{mailto:\PaperEmail}{\PaperEmail}}', '')
    body = body.replace(r'\begin{abstract}', r'\twocolumn[\maketitle\begin{abstract}').replace(r'\end{abstract}' + '\n\n' + r'\maketitle', r'\end{abstract}\vspace{1em}]')
    chunks = ''.join('\\par\\vspace{%dpt}\\mbox{}\\par\n' % min(100, K - 100 * i) for i in range((K + 99) // 100))
    body = body.replace('\\input{sections/02_setup.tex}\n', '\\input{sections/02_setup.tex}\n' + chunks)
    body = body.replace('\\bibliography{references}', '\\bibliographystyle{unsrt}\\bibliography{references}')
    (MAN / f'{name}.tex').write_text(PRE + body)
    for cmd in (['pdflatex'], ['bibtex'], ['pdflatex'], ['pdflatex']):
        subprocess.run(cmd + (['-interaction=nonstopmode', f'{name}.tex'] if cmd[0] == 'pdflatex' else [name]), cwd=MAN, capture_output=True)
    aux = (MAN / f'{name}.aux').read_text()
    pages = {k: int(p) for k, p in re.findall(r'\\newlabel\{([^}]*)\}\{\{[^}]*\}\{(\d+)\}', aux)}
    npages = int(re.search(r'Pages:\s+(\d+)', subprocess.run(['pdfinfo', f'{name}.pdf'], cwd=MAN, capture_output=True, text=True).stdout).group(1))
    log = (MAN / f'{name}.log').read_text(errors='ignore')
    return pages, npages, log

if __name__ == '__main__':
    shifts = [int(a) for a in sys.argv[1:]] or list(range(0, 1321, 120))
    report = pathlib.Path(__file__).resolve().parent / 'check_layout.py'
    allok = True
    for K in shifts:
        P, n, log = build(K, 'shift')
        out = subprocess.run([sys.executable, str(report), '--job', 'shift'], capture_output=True, text=True).stdout
        bad = [l for l in out.splitlines() if any(w in l for w in ('LATE', 'EARLY', 'STUCK', 'FLUSHED'))]
        errs = len(re.findall(r'^!', log, re.M)) + log.count('Overfull')
        ok = not bad and errs == 0
        allok &= ok
        print(f" shift {K:5d} pt: {n} pages, appendices p.{P['app:constants']}-{P['end-of-appendices']}, "
              f"Fig. 5 p.{P['fig:cross-section']}, Fig. 6 p.{P['fig:flavour-components']}, Table VIII p.{P['tab:timing']}: "
              + ("ok" if ok else "; ".join(bad) + (f" [{errs} errors/overfull]" if errs else "")))
    for f in MAN.glob('shift.*'):
        f.unlink()
    print("ALL OK" if allok else "PROBLEMS FOUND")
