# Brain to Values

Visualizations and papers on mapping brain regions and processing stages to intrinsic human values.

## Brain–value diagram

![Brain to values](viz/brain_values/output/brain_values.png)

[Open brain to values as clickable PDF](viz/brain_values/output/brain_values.pdf)

```bash
pip3 install -r requirements.txt
python3 viz/brain_values/brain.py
```

## Sensory bandwidth diagram

![Sensory processing](viz/senses/output/senses.png)

```bash
python3 viz/senses/senses.py
```

## Papers

| Paper | Source | PDF |
|-------|--------|-----|
| Stratification of free-energy loops | [papers/free-energy-loops/free_energy_loops.tex](papers/free-energy-loops/free_energy_loops.tex) | [PDF](papers/free-energy-loops/free_energy_loops.pdf) |
| Loop–hub–control–value (LHCV) | [papers/loop-hub-control-value/lhcv-model-v2.tex](papers/loop-hub-control-value/lhcv-model-v2.tex) | [PDF](papers/loop-hub-control-value/lhcv-model-v2.pdf) |
| Loop–hub–value model | [papers/loop-hub-value-model/loop-hub-value-model.tex](papers/loop-hub-value-model/loop-hub-value-model.tex) | [PDF](papers/loop-hub-value-model/loop-hub-value-model.pdf) |
| Status regulation loops | [papers/status-regulation-loops/status_regulation_as_free_energy_loops.tex](papers/status-regulation-loops/status_regulation_as_free_energy_loops.tex) | [PDF](papers/status-regulation-loops/status_regulation_as_free_energy_loops.pdf) |
| Value bundle drift | [papers/value-bundle-drift/value-bundle-drift.tex](papers/value-bundle-drift/value-bundle-drift.tex) | [PDF](papers/value-bundle-drift/value-bundle-drift.pdf) |
| Unit of caring | [papers/unit-of-caring/unit-of-caring.tex](papers/unit-of-caring/unit-of-caring.tex) | [PDF](papers/unit-of-caring/unit-of-caring.pdf) |

Build a paper:

```bash
papers/free-energy-loops/build.sh
papers/loop-hub-control-value/build.sh
papers/loop-hub-value-model/build.sh   # also generates the brain–value schematic
papers/status-regulation-loops/build.sh
papers/value-bundle-drift/build.sh
papers/unit-of-caring/build.sh         # pdflatex + bibtex + pdflatex ×2
```
