# From Orienting Signals to Self-Pointers

LaTeX paper and companion PyTorch experiment on emergent self-binding from privileged organism-centered cues.

## Paper

| File | Role |
|------|------|
| [`self-pointers.tex`](self-pointers.tex) | Source |
| [`self-pointers.bib`](self-pointers.bib) | Bibliography |
| [`build.sh`](build.sh) | Build PDF (regenerates figures, then `pdflatex` + `bibtex` ×2) |
| [`figures/`](figures/) | Paper figures |

### Figures

| Figure | Source |
|--------|--------|
| Matched counterfactual geometry | [`figures/self_pointer_diagram.py`](figures/self_pointer_diagram.py) → `self_pointer_counterfactual.pdf` |
| Cue-ablation learning curves | [`self_pointer_experiment/plot_learning_curves.py`](self_pointer_experiment/plot_learning_curves.py) → copy or regenerate into `figures/` |

Build from repo root:

```bash
papers/self-pointers/build.sh
```

## Experiment

Code lives in [`self_pointer_experiment/`](self_pointer_experiment/). See [`self_pointer_experiment/README.md`](self_pointer_experiment/README.md) for environment, ablations, probes, and run commands.

Quickstart:

```bash
cd papers/self-pointers/self_pointer_experiment
pip install -r requirements.txt
PYTHONPATH=. pytest -q
PYTHONPATH=. python train.py --experiment all --steps 2000
PYTHONPATH=. python plot_learning_curves.py --steps 300 --eval-every 10 --n-seeds 10
```

## Current results (paper Table/Figure)

Ten-seed learning curves to 300 updates (`eval-every 10`, held-out accuracy):

- **Controls:** vision-only and proprioception-only ≈ 50%
- **Orienting-related:** delayed phase transition, then >93% by round 300
- **Combinations:** modest lead over orient-only at convergence; clearer advantage early in training
- **Explicit self:** ≈99% (oracle upper bound)

Mechanistic readouts (probes, state-swap) are specified in the paper but not yet reported at scale in the public results JSON.
