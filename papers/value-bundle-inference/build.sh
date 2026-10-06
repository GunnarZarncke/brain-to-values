#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
paper="value-bundle-inference"

log=$(mktemp)
trap 'rm -f "$log"' EXIT

run() { "$@" >>"$log" 2>&1 || { cat "$log"; exit 1; }; }

run dot -Tpdf -o figures/capability-observability.pdf figures/capability-observability.dot
run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"
run bibtex "$paper"
run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"
run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"

echo "Built $paper.pdf"
