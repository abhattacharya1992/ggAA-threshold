.PHONY: all clean format

all: main.pdf

# latexmk tracks all included sources, figures and bibliography dependencies.
# The formatting step runs after any table-generation targets in the root Makefile.
main.pdf: format
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

format:
	python3 ../scripts/format_jhep_tables.py

clean:
	latexmk -C main.tex
