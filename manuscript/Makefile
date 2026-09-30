.PHONY: all clean

all: main.pdf

main.pdf: main.tex numbers.tex macros.tex metadata.tex references.bib sections/*.tex tables/*.tex figures/*.pdf
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

clean:
	latexmk -C main.tex
	rm -f mainNotes.bib
