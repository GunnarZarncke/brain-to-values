#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"
paper="loop-hub-value-model"
repo_root="$(cd ../.. && pwd)"
figure="$repo_root/papers/loop-hub-value-model/figures/brain-values-schematic.pdf"
brain_script="$repo_root/viz/brain_values/brain.py"

mkdir -p figures

if [[ -x "$repo_root/.venv/bin/python" ]]; then
  python="$repo_root/.venv/bin/python"
else
  python="python3"
fi

MPLBACKEND=Agg "$python" "$brain_script" --paper
cp "$repo_root/viz/brain_values/output/brain_values.pdf" "$figure"

log=$(mktemp)
trap 'rm -f "$log"' EXIT

run() { "$@" >>"$log" 2>&1 || { cat "$log"; exit 1; }; }

run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"
run pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex"

echo "Built $paper.pdf"
