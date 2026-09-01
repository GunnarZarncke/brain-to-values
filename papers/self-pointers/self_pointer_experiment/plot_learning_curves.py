"""Train cue ablations and plot held-out accuracy over training rounds."""

import argparse
import json
import os
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import torch
from torch.nn import functional as F

from self_pointer.analysis import task_accuracy
from self_pointer.env import CueConfig, generate_paired_batch
from self_pointer.model import ModelConfig, SelfPointerNet

EXPERIMENTS = {
    "vision_only": CueConfig(),
    "orient_only": CueConfig(orient=True),
    "efference_only": CueConfig(efference=True),
    "proprio_only": CueConfig(proprio=True),
    "orient_efference": CueConfig(orient=True, efference=True),
    "orient_proprio": CueConfig(orient=True, proprio=True),
    "all_internal": CueConfig(orient=True, efference=True, proprio=True),
    "explicit_self": CueConfig(explicit_self=True),
}

LABELS = {
    "vision_only": "Vision only",
    "orient_only": "Orienting only",
    "efference_only": "Efference only",
    "proprio_only": "Proprioception only",
    "orient_efference": "Orienting + efference",
    "orient_proprio": "Orienting + proprioception",
    "all_internal": "All internal cues",
    "explicit_self": "Explicit self (oracle)",
}


def train_curve(
    name,
    steps,
    batch_worlds,
    horizon,
    lr,
    device,
    seed,
    eval_worlds,
    eval_every,
):
    torch.manual_seed(seed)
    cues = EXPERIMENTS[name]
    model = SelfPointerNet(ModelConfig(recirculate=True)).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr)

    rounds = []
    accuracies = []
    for step in range(steps):
        batch = generate_paired_batch(
            batch_worlds, horizon, cues, device=device, seed=seed + step
        )
        logits = model(batch.visual, batch.cues)
        loss = F.cross_entropy(logits, batch.label)
        opt.zero_grad()
        loss.backward()
        opt.step()

        if (step + 1) % eval_every == 0 or step + 1 == steps:
            with torch.no_grad():
                acc = task_accuracy(
                    model,
                    cues,
                    n_worlds=eval_worlds,
                    horizon=horizon,
                    device=device,
                    seed=seed + 100_000,
                )
            rounds.append(step + 1)
            accuracies.append(acc)

    del model, opt
    if device.startswith("cuda"):
        torch.cuda.empty_cache()

    return {"seed": seed, "rounds": rounds, "held_out_accuracy": accuracies}


def aggregate_runs(name, seed_runs):
    rounds = seed_runs[0]["rounds"]
    matrix = np.array([run["held_out_accuracy"] for run in seed_runs], dtype=float)
    return {
        "name": name,
        "rounds": rounds,
        "n_seeds": len(seed_runs),
        "mean_accuracy": matrix.mean(axis=0).tolist(),
        "std_accuracy": matrix.std(axis=0, ddof=0).tolist(),
        "seed_runs": seed_runs,
    }


def plot_curves(curves, outpath, title, show_std=True):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.grid(True, alpha=0.3)

    order = list(EXPERIMENTS)
    cmap = plt.get_cmap("tab10")

    for i, name in enumerate(order):
        series = next(c for c in curves if c["name"] == name)
        rounds = np.array(series["rounds"])
        if "mean_accuracy" in series:
            mean = np.array(series["mean_accuracy"])
            std = np.array(series.get("std_accuracy", []))
        else:
            mean = np.array(series["held_out_accuracy"])
            std = np.zeros_like(mean)

        color = cmap(i % 10)
        lw = 2.2 if name in {"all_internal", "explicit_self"} else 1.7
        ax.plot(rounds, mean, label=LABELS[name], linewidth=lw, color=color)
        if show_std and len(std) and np.any(std > 0):
            ax.fill_between(
                rounds,
                np.clip(mean - std, 0.0, 1.0),
                np.clip(mean + std, 0.0, 1.0),
                color=color,
                alpha=0.18,
                linewidth=0,
            )

    ax.axhline(0.5, color="0.45", linestyle="--", linewidth=1.2, label="Chance (50%)")
    ax.set_xlim(1, max(curves[0]["rounds"]))
    ax.set_ylim(0.45, 1.02)
    ax.set_xlabel("Training round")
    ax.set_ylabel("Held-out accuracy")
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=8, frameon=True)
    fig.tight_layout()
    for suffix in (".png", ".svg", ".pdf"):
        fig.savefig(outpath.with_suffix(suffix), dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--batch-worlds", type=int, default=128)
    parser.add_argument("--eval-worlds", type=int, default=256)
    parser.add_argument("--eval-every", type=int, default=10)
    parser.add_argument("--horizon", type=int, default=8)
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--n-seeds", type=int, default=1)
    parser.add_argument(
        "--device", default="cuda" if torch.cuda.is_available() else "cpu"
    )
    parser.add_argument("--outdir", default="results/figures")
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Plot from existing ablation_learning_curves.json without retraining",
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data_path = outdir / "ablation_learning_curves.json"

    if args.plot_only:
        curves = json.loads(data_path.read_text())
    else:
        seeds = [args.seed + k for k in range(args.n_seeds)]
        curves = []
        for name in EXPERIMENTS:
            seed_runs = []
            for seed in seeds:
                print(f"Training {name} seed={seed} ...")
                seed_runs.append(
                    train_curve(
                        name,
                        args.steps,
                        args.batch_worlds,
                        args.horizon,
                        args.lr,
                        args.device,
                        seed,
                        args.eval_worlds,
                        args.eval_every,
                    )
                )
            curves.append(aggregate_runs(name, seed_runs))
            data_path.write_text(json.dumps(curves, indent=2))

    n_seeds = curves[0].get("n_seeds", 1)
    title = "Cue ablations: held-out accuracy over training"
    if n_seeds > 1:
        title += f" (mean ± std, {n_seeds} seeds)"

    plot_curves(curves, outdir / "ablation_learning_curves", title)
    print(f"Wrote {data_path}")
    print(f"Wrote {outdir / 'ablation_learning_curves.png'}")
    print(f"Wrote {outdir / 'ablation_learning_curves.svg'}")


if __name__ == "__main__":
    main()
