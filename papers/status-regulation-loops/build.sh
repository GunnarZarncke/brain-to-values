#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
paper="status_regulation_as_free_energy_loops"

pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex" >/dev/null
bibtex "$paper" >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex" >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex" >/dev/null

echo "Built $paper.pdf"
