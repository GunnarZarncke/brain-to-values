# AGENTS.md

Guidance for AI agents working in **brain-to-values**: visualizations and LaTeX papers on mapping brain regions, processing stages, and agent architectures to intrinsic values and caring functionals.

## Project layout

```
brain-to-values/
├── README.md              # User-facing overview and quickstart
├── requirements.txt       # Python deps (matplotlib, numpy)
├── viz/                   # Matplotlib figure generators
│   ├── brain_values/      # Hub-centric brain–value diagram
│   │   ├── brain.py
│   │   ├── assets/        # Input images (e.g. brain.png)
│   │   └── output/        # Generated PNG/SVG/PDF
│   └── senses/            # Sensory bandwidth funnel diagram
│       ├── senses.py
│       └── output/
└── papers/                # One subdirectory per paper
    ├── loop-hub-value-model/
    ├── value-bundle-drift/
    └── unit-of-caring/    # Draft; no build.sh yet
```

Each paper directory should contain:

- `<name>.tex` — source
- `<name>.pdf` — built output (commit when publishing)
- `figures/` — static or generated figures referenced by TeX
- `build.sh` — reproducible build script (when ready)

## Commands

### Visualizations

```bash
pip3 install -r requirements.txt
python3 viz/brain_values/brain.py          # interactive; writes to output/
python3 viz/brain_values/brain.py --paper  # print-friendly; used by LHV build
python3 viz/senses/senses.py
```

Use `MPLBACKEND=Agg` for headless runs. Prefer `.venv/bin/python` if the venv exists.

### Papers

```bash
papers/loop-hub-value-model/build.sh   # generates schematic, then pdflatex ×2
papers/value-bundle-drift/build.sh     # pdflatex ×2
```

**LHV paper dependency:** `loop-hub-value-model/build.sh` runs `brain.py --paper` and copies `viz/brain_values/output/brain_values.pdf` to `papers/loop-hub-value-model/figures/brain-values-schematic.pdf` before compiling.

**Unit of caring:** draft only. Has `unit-of-caring.tex`, `unit-of-caring.bib`, and reference PDFs; no `build.sh` yet. Will need `bibtex` when bibliography is wired up.

## Conventions

### Python (`viz/`)

- Resolve paths relative to the script with `Path(__file__).resolve().parent`.
- Write outputs to `output/` under the script directory, not the repo root.
- Keep CONFIG/data near the top of scripts; match existing matplotlib style.
- `brain.py --paper` omits click handlers and `plt.show()` for LaTeX embedding.

### LaTeX (`papers/`)

- Place figures in `figures/` and reference as `figures/<file>`.
- Inline bibliographies (LHV) need only `pdflatex` twice.
- External `.bib` files (unit-of-caring) will need `pdflatex → bibtex → pdflatex ×2`.
- Add `\usepackage{graphicx}` when including images.
- LaTeX build artifacts (`*.aux`, `*.log`, `*.out`, `*.bak`) are gitignored.

### Build scripts

Follow the existing pattern:

```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
paper="<name>"
# optional: generate figures from viz/
pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex" >/dev/null
pdflatex -interaction=nonstopmode -halt-on-error "$paper.tex" >/dev/null
echo "Built $paper.pdf"
```

### README

Update `README.md` when adding a new viz or paper (table row + build command). Keep paths relative to repo root.

## Adding a new paper

1. Create `papers/<slug>/` with `<slug>.tex`, `figures/` as needed.
2. Add `build.sh` (executable).
3. Build and verify the PDF.
4. Add a row to the README papers table.

## Adding a new visualization

1. Create `viz/<name>/` with script, optional `assets/`, and `output/`.
2. Use path-relative I/O (see `viz/brain_values/brain.py`).
3. Document the run command in README.

## Git

- **Do not commit unless the user asks.**
- Do not commit venv directories, LaTeX aux/log files, or secrets.
- Commit source `.tex`, `.py`, `build.sh`, and publishable PDFs/figures when requested.

## Known issues

- **`unit-of-caring.tex`** contains `:contentReference[oaicite:…]` placeholders from an export; remove or replace when editing.
- **`unit-of-caring.tex`** cites `\cite{…}` but has no `\bibliography{…}` block yet; bibliography is incomplete.
- **README** does not yet list the unit-of-caring paper; add when it has a working build.
- **`brain.py` docstring** and interactive title differ from `--paper` export; intentional.

## Thematic links between artifacts

| Artifact | Role |
|----------|------|
| `viz/brain_values/` | Empirical hub→value map; feeds LHV paper Figure 1 |
| `viz/senses/` | Illustrative bandwidth compression across processing stages |
| `papers/loop-hub-value-model/` | Loop–hub–value (LHV) formal model |
| `papers/value-bundle-drift/` | Cultural selection vs value alignment |
| `papers/unit-of-caring/` | Integrity pressure, suffering, aggregation scaffold |

When changing hub–value mappings in `brain.py`, rebuild the LHV paper so the schematic stays in sync.

## Scope discipline

- Minimize diffs; do not refactor unrelated code.
- Do not add markdown files the user did not ask for (except this file).
- Prefer extending existing scripts and build patterns over new tooling.
