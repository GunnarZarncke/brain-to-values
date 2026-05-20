# Visualizations

Matplotlib figure generators for the brain-to-values project. Outputs are written to each script's `output/` directory.

## Brain–value diagram

![Brain to values](brain_values/output/brain_values.png)

Hub-centric schematic mapping neuromodulatory hubs and processing stages to intrinsic value readouts. Used as Figure 1 in the [loop–hub–value model](../papers/loop-hub-value-model/loop-hub-value-model.pdf).

```bash
pip3 install -r ../requirements.txt
python3 brain_values/brain.py              # interactive; writes PNG/SVG/PDF
python3 brain_values/brain.py --paper      # print-friendly PDF for LaTeX embed
```

## Sensory bandwidth diagram

![Sensory processing](senses/output/senses.png)

Illustrative funnel of bandwidth compression across sensory processing stages.

```bash
python3 senses/senses.py
```

Use `MPLBACKEND=Agg` for headless runs. Prefer `.venv/bin/python` from the repo root if a venv exists.
