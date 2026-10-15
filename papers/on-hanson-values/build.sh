#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
paper="on-hanson-values"

log=$(mktemp)
trap 'rm -f "$log"' EXIT

run() { "$@" >>"$log" 2>&1 || { cat "$log"; exit 1; }; }

py=python3
[[ -x ../../.venv/bin/python3.14 ]] && py=../../.venv/bin/python3.14
MPLBACKEND=Agg run "$py" make_figures.py

run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"
run bibtex "$paper"
run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"
run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"

echo "Built $paper.pdf"
