# Reproduce the results of the paper.  `make help` lists the targets.
# pipefail: a script that fails inside `... | tee file` stops make instead of passing silently.
SHELL := /bin/bash
.SHELLFLAGS := -o pipefail -c
PYTHON ?= python3
RUN = $(PYTHON) scripts

.PHONY: jz2 nlp-dps m1-forward provenance nf1-inputs reproduce arxiv crossover dps-hh renorm-2019 layout archive log-share diagrams help all inputs series constants leading-logs one-loop exact-check validation tables figures physics-figures all-figures numbers manuscript-data dps fig4-2019 rate-logs top-mass lo virtual normalization test paper clean

help:
	@echo "make inputs        download the arXiv:1909.08993 ancillary files (needed for one-loop)"
	@echo "make series        reduce data/construction to data/series (exact cancellations)"
	@echo "make constants     recompute the 13 endpoint constants and their reductions"
	@echo "make leading-logs  leading-power structure in terms of calL"
	@echo "make one-loop      one-loop threshold and forward limits (needs inputs)"
	@echo "make exact-check   independent evaluation of the exact two-loop remainders (needs inputs, ~45 min)"
	@echo "make validation    series against the 1080 validation records; Table 3"
	@echo "make lo            LO partonic cross section of the Type-II amplitude"
	@echo "make virtual       finite virtual corrections V1, V2 relative to LO (needs inputs)"
	@echo "make normalization checks of the tree, one-loop and two-loop normalization (needs inputs)"
	@echo "make tables        forward-angle and endpoint-constant tables"
	@echo "make figures       Figs. 1 (diagrams), 2, 4 and 6"
	@echo "make physics-figures  Figs. 3 and 5 (needs inputs; Fig. 5 needs make virtual)"
	@echo "make all-figures   regenerate all manuscript figures; runs virtual first for Fig. 5"
	@echo "make test          quick tests (under a minute)"
	@echo "make dps           one-loop ln p_T term against the double-parton-scattering singularity (needs inputs)"
	@echo "make rate-logs     ln p_T terms of the finite virtual corrections V1, V2 (needs inputs)"
	@echo "make top-mass      top quark mass factor of the Type II tree amplitude at threshold"
	@echo "make fig4-2019     ratio curves for comparison with Fig. 4 of arXiv:1909.08993"
	@echo "make numbers       manuscript/numbers.tex: every number quoted in the text, from data/ and results/ (needs inputs)"
	@echo "make manuscript-data  tables and numbers (fast; uses the stored results)"
	@echo "make paper         regenerate tables and numbers, build manuscript/main.pdf, print the layout report (needs inputs)"
	@echo "make arxiv         after make paper: build/ggAA-threshold-arxiv.tar.gz with main.bbl and anc/"
	@echo "make crossover     velocity below which the two-loop term of M1 exceeds the one-loop term"
	@echo "make dps-hh        one-loop gg -> HH: forward logarithm against F_S (data/dps_hh; Fig. 4, Sec. 4.3)"
	@echo "make renorm-2019   earlier 2019 run vs published renormalisation; June 2019 pole record (App. B)"
	@echo "make jz2           forward limit of the two-loop |J_z|=2 amplitude (Sec. 4.4), from the stored points; needs inputs"
	@echo "make nlp-dps       two-loop DPS logarithms beyond leading power in the threshold series (Sec. 4.4)"
	@echo "make m1-forward    forward limit of the two-loop M_1 (Sec. 4.4, Fig. 4), from the stored points; needs inputs"
	@echo "make provenance    data/construction and data/validation against the outputs of the construction runs (data/runs)"
	@echo "make nf1-inputs    the two inputs of the P1 n_f construction against the published remainder"
	@echo "make reproduce     rerun the construction from the published files (needs inputs, Mathematica 13.2, g++; several hours)"
	@echo "make all           series, constants, leading-logs, validation, tables, figures, numbers (needs inputs)"

all: series constants leading-logs validation tables figures numbers

inputs:
	bash inputs/fetch_inputs.sh

series:
	@mkdir -p results
	$(RUN)/reduce_series.py | tee results/reduce_series.txt

constants:
	@mkdir -p results
	$(RUN)/check_mb_constants.py | tee results/mb_constants.txt

leading-logs:
	@mkdir -p results
	$(RUN)/check_leading_logs.py | tee results/leading_logs.txt

one-loop:
	@mkdir -p results
	$(RUN)/check_one_loop.py | tee results/one_loop.txt

exact-check:
	@mkdir -p results
	$(RUN)/check_exact_amplitude.py | tee results/exact_amplitude.txt

validation:
	@mkdir -p results
	$(RUN)/compare_validation.py | tee results/compare_validation.txt

lo:
	@mkdir -p results
	$(RUN)/lo_cross_section.py | tee results/lo_cross_section.txt

virtual:
	@mkdir -p results
	$(RUN)/virtual_correction.py | tee results/virtual_correction.txt

normalization:
	@mkdir -p results
	$(RUN)/check_normalization.py | tee results/normalization.txt

tables:
	$(RUN)/make_tables.py

figures:
	$(RUN)/make_figures.py
	$(RUN)/make_diagrams.py

diagrams:
	$(RUN)/make_diagrams.py

physics-figures:
	$(RUN)/make_physics_figures.py

all-figures: figures virtual
	$(RUN)/make_physics_figures.py
	
test:
	$(PYTHON) -m unittest discover -s tests -v

dps:
	$(RUN)/check_dps.py

rate-logs:
	$(RUN)/check_rate_logs.py

crossover:
	$(RUN)/check_amplitude_crossover.py

dps-hh:
	$(RUN)/check_dps_hh.py

jz2:
	$(RUN)/check_jz2_forward.py --dps 140 --exponents 5 6 7 --workers 1
	cp results/jz2_forward.txt results/jz2_forward_nc3_nf5_e567.txt
	$(RUN)/check_jz2_forward.py --dps 140 --exponents 5 6 7 --Nc 4 --nf 2 --workers 1
	cp results/jz2_forward.txt results/jz2_forward_nc4_nf2_e567.txt

nlp-dps:
	$(RUN)/check_nlp_dps.py | tee results/nlp_dps.txt

m1-forward:
	$(RUN)/check_m1_forward.py --dps 180 --check-extra 40 --exponents 7 8 9 10 --workers 1
	$(RUN)/check_m1_forward.py --dps 180 --check-extra 40 --exponents 7 8 9 10 --workers 1 --Nc 4 --nf 2

provenance:
	$(RUN)/check_provenance.py

nf1-inputs:
	$(RUN)/check_nf1_inputs.py

reproduce:
	bash src/construction/reproduce_series.sh --inputs inputs/arXiv-1909.08993/anc \
		--nf1-inputs src/construction/inputs --work build/reproduce --package .

renorm-2019:
	$(RUN)/check_renorm_2019.py

log-share:
	$(RUN)/check_log_share.py

archive:
	$(RUN)/check_archive_2019.py

top-mass:
	$(RUN)/check_top_mass.py

fig4-2019:
	$(RUN)/check_fig4_2019.py

numbers:
	$(RUN)/make_numbers.py

manuscript-data: tables numbers

paper: manuscript-data
	@if [ -d manuscript ]; then $(MAKE) -C manuscript && $(PYTHON) scripts/check_layout.py; \
	else echo "arXiv layout: the manuscript is in the directory above anc/; run pdflatex main.tex there three times."; fi

layout:
	$(PYTHON) scripts/check_layout.py

arxiv:
	$(RUN)/make_arxiv.py

clean:
	$(MAKE) -C manuscript clean
