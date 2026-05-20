#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
paper="value-bundle-drift"

pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex" >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex" >/dev/null

echo "Built $paper.pdf"
