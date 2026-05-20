# Papers

LaTeX sources and built PDFs for each paper. Rebuild with the corresponding `build.sh` script.

| Paper | Source | PDF | Build |
|-------|--------|-----|-------|
| Stratification of free-energy loops | [free_energy_loops.tex](free-energy-loops/free_energy_loops.tex) | [PDF](free-energy-loops/free_energy_loops.pdf) | `free-energy-loops/build.sh` |
| Loop–hub–control–value (LHCV) | [lhcv-model-v2.tex](loop-hub-control-value/lhcv-model-v2.tex) | [PDF](loop-hub-control-value/lhcv-model-v2.pdf) | `loop-hub-control-value/build.sh` |
| Loop–hub–value model (LHV) | [loop-hub-value-model.tex](loop-hub-value-model/loop-hub-value-model.tex) | [PDF](loop-hub-value-model/loop-hub-value-model.pdf) | `loop-hub-value-model/build.sh` |
| Status regulation loops | [status_regulation_as_free_energy_loops.tex](status-regulation-loops/status_regulation_as_free_energy_loops.tex) | [PDF](status-regulation-loops/status_regulation_as_free_energy_loops.pdf) | `status-regulation-loops/build.sh` |
| Value bundle drift | [value-bundle-drift.tex](value-bundle-drift/value-bundle-drift.tex) | [PDF](value-bundle-drift/value-bundle-drift.pdf) | `value-bundle-drift/build.sh` |
| Unit of caring | [unit-of-caring.tex](unit-of-caring/unit-of-caring.tex) | [PDF](unit-of-caring/unit-of-caring.pdf) | `unit-of-caring/build.sh` |

Build from the repo root:

```bash
papers/free-energy-loops/build.sh
papers/loop-hub-control-value/build.sh
papers/loop-hub-value-model/build.sh   # also regenerates the brain–value schematic
papers/status-regulation-loops/build.sh         # pdflatex + bibtex + pdflatex ×2
papers/value-bundle-drift/build.sh
papers/unit-of-caring/build.sh         # pdflatex + bibtex + pdflatex x2
```

LaTeX build artifacts (`*.aux`, `*.log`, etc.) are gitignored. Paper PDFs are tracked in git; regenerate after source changes.
