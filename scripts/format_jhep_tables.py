#!/usr/bin/env python3
"""Format complete tables and legacy table fragments for the JHEP manuscript.

Both balanced and partly removed REVTeX wrappers are handled. Numerical
entries are retained. Run by the manuscript Makefile before compilation.
"""

from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "manuscript"
if not MANUSCRIPT.is_dir() and (ROOT.parent / "main.tex").is_file():
    MANUSCRIPT = ROOT.parent

BEGIN_RULED = r"\begin{ruledtabular}"
END_RULED = r"\end{ruledtabular}"


def group_end(text, start):
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{" and (i == 0 or text[i - 1] != "\\"):
            depth += 1
        elif text[i] == "}" and (i == 0 or text[i - 1] != "\\"):
            depth -= 1
            if depth == 0:
                return i + 1
    raise ValueError("Unbalanced TeX braces")


def tabular_rules(text):
    """Add booktabs rules inside tabular environments in a table block."""
    pattern = re.compile(r"\\begin\{tabular\}(?:\[[^]]*\])?\s*\{")
    end_pattern = re.compile(r"\\end\{tabular\}")
    cursor = 0
    while True:
        begin = pattern.search(text, cursor)
        if begin is None:
            return text
        columns_end = group_end(text, begin.end() - 1)
        end = end_pattern.search(text, columns_end)
        if end is None:
            raise ValueError("Unclosed tabular environment")
        body = text[columns_end:end.start()]
        body = re.sub(r"\\hline\b", lambda _: r"\midrule", body).strip()
        if not re.match(r"\\toprule\b", body):
            body = "\\toprule\n" + body
        bottom = re.search(r"\\bottomrule(?:\[[^]]*\])?\s*$", body)
        if bottom:
            rows = body[:bottom.start()].rstrip()
            final_rule = body[bottom.start():].strip()
        else:
            rows = re.sub(r"\\midrule(?:\[[^]]*\])?\s*$", "", body).rstrip()
            final_rule = r"\bottomrule"
        lines = rows.splitlines()
        if lines and "&" in lines[-1] and not lines[-1].rstrip().endswith("\\\\"):
            lines[-1] = lines[-1].rstrip() + " \\\\"
        body = "\n" + "\n".join(lines) + "\n" + final_rule + "\n"
        text = text[:columns_end] + body + text[end.start():]
        cursor = columns_end + len(body) + len(end.group(0))


def format_table(match):
    text = match.group(0).replace(BEGIN_RULED, "").replace(END_RULED, "")
    text = tabular_rules(text)
    cap = text.find(r"\caption{")
    tab = text.find(r"\begin{tabular}")
    if cap >= 0 and tab >= 0 and cap < tab:
        end = group_end(text, cap + len(r"\caption"))
        label = re.match(r"\s*\\label\{[^}]+\}", text[end:])
        if label:
            end += label.end()
        caption = text[cap:end].strip()
        text = text[:cap] + text[end:]
        end_table = text.rfind(r"\end{table")
        text = text[:end_table].rstrip() + "\n" + caption + "\n" + text[end_table:]
    return text


def format_text(text):
    # Complete ruledtabular fragments also occur outside table floats.
    text = re.sub(
        r"\\begin\{ruledtabular\}(.*?)\\end\{ruledtabular\}",
        lambda m: tabular_rules(m.group(1)), text, flags=re.S,
    )
    # Partly converted generators may leave only the opening or closing token.
    text = text.replace(BEGIN_RULED, "").replace(END_RULED, "")
    text = re.sub(
        r"\\begin\{table\*?\}.*?\\end\{table\*?\}",
        format_table, text, flags=re.S,
    )
    if BEGIN_RULED in text or END_RULED in text:
        raise ValueError("REVTeX table wrappers remain")
    return text


def main():
    if not MANUSCRIPT.is_dir():
        raise SystemExit("No manuscript directory found")
    changed = 0
    for path in sorted(MANUSCRIPT.rglob("*.tex")):
        before = path.read_text()
        after = format_text(before)
        if after != before:
            path.write_text(after)
            changed += 1
    print(f"JHEP table formatting: {changed} file(s) updated")


if __name__ == "__main__":
    main()
