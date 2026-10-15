#!/usr/bin/env python3
"""Figures for on-hanson-values.tex: Study 1 from core_both, Study 2 from abc_extended_r3."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "core_both"
DATA2 = HERE / "abc_extended_r3"
OUT = HERE / "figures"

CONDITIONS = ["A_1930", "B_transplant", "C_elites", "D_full"]
COND_LABELS = ["A\n1930", "B\nfacts", "C\n+elites", "D\n+social"]
CONTRASTS = ["ecology_B_minus_A", "elite_C_minus_B", "social_D_minus_C", "total_D_minus_A"]
CONTRAST_LABELS = ["B$-$A ecology", "C$-$B elites", "D$-$C social", "D$-$A total"]
METRICS = [("recommendation", "Recommends B (turn 1)"),
           ("equalized_preference", "Prefers B at equal payoff (turn 3)"),
           ("elicited_B_moralization", "Moralizes B when asked (turn 2)")]
ARMS = [("talkie", "Talkie-1930", "#8c2d04", "o"), ("qwen1930", "Qwen3-14B as 1930", "#2b6c9e", "s")]
N_BOOT = 5000


def case_bootstrap(case_means, rng):
    values = case_means.dropna().to_numpy()
    draws = rng.choice(values, size=(N_BOOT, len(values)), replace=True).mean(axis=1)
    return values.mean(), *np.percentile(draws, [2.5, 97.5])


def trajectories(coded):
    rng = np.random.default_rng(1930)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.1), sharey=True)
    for ax, (metric, title) in zip(axes, METRICS):
        for k, (arm, label, color, marker) in enumerate(ARMS):
            rows = []
            for cond in CONDITIONS:
                sub = coded[(coded.arm == arm) & (coded.condition == cond)]
                rows.append(case_bootstrap(sub.groupby("case_id")[metric].mean(), rng))
            mean, lo, hi = map(np.array, zip(*rows))
            x = np.arange(len(CONDITIONS)) + (k - 0.5) * 0.12
            ax.errorbar(x, mean, yerr=[mean - lo, hi - mean], color=color, marker=marker,
                        capsize=2, lw=1.4, ms=5, label=label)
        ax.set_xticks(range(len(CONDITIONS)), COND_LABELS, fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_ylim(-0.02, 1.02)
        ax.axhline(0.5, color="0.8", lw=0.6, zorder=0)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("share coded B (A=0, mixed=.5, B=1)", fontsize=8)
    axes[0].legend(fontsize=7, loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "trajectories.pdf")
    plt.close(fig)


def contrasts(effects):
    fig, axes = plt.subplots(1, 3, figsize=(10, 2.9), sharey=True)
    for ax, (metric, title) in zip(axes, METRICS):
        for k, (arm, label, color, marker) in enumerate(ARMS):
            sub = effects[(effects.arm == arm) & (effects.metric == metric)].set_index("contrast")
            sub = sub.reindex(CONTRASTS)
            y = np.arange(len(CONTRASTS))[::-1] + (0.5 - k) * 0.22
            ax.errorbar(sub.mean_delta, y,
                        xerr=[sub.mean_delta - sub.ci_low, sub.ci_high - sub.mean_delta],
                        fmt=marker, color=color, capsize=2, ms=5, label=label)
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set_yticks(np.arange(len(CONTRASTS))[::-1], CONTRAST_LABELS, fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_xlim(-0.75, 1.0)
        ax.grid(axis="x", alpha=0.25)
        ax.set_xlabel("paired difference (8-domain bootstrap 95% CI)", fontsize=7)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=7, loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "contrasts.pdf")
    plt.close(fig)


EXT_CONTRASTS = CONTRASTS + ["hanson_C_opaque_minus_explained", "elite_direction_C_minus_opposite",
                             "explicit_moral_message", "peer_exposure", "media_exposure",
                             "institutional_exposure"]
EXT_LABELS = CONTRAST_LABELS + ["C$-$C$_{\\mathrm{expl}}$ opacity", "C$-$C$_{\\mathrm{opp}}$ direction",
                                "C$_{\\mathrm{moral}}-$C moral claim", "D$_{\\mathrm{peers}}-$C peers",
                                "D$_{\\mathrm{media}}-$D$_{\\mathrm{peers}}$ media",
                                "D$-$D$_{\\mathrm{media}}$ institutions"]
BATTERIES = [("A", "Battery A: original 8\n(new facts favor B)"),
             ("B", "Battery B: 8 new\n(new facts favor B)"),
             ("C", "Battery C: 8 drift tests\n(new facts favor A)")]


def battery_map():
    ids = [c["id"] for c in json.loads((DATA2 / "manifest.json").read_text())["cases"]]
    return {cid: "ABC"[k // 8] for k, cid in enumerate(ids)}


def study2_effects(case_effects):
    rng = np.random.default_rng(1930)
    ce = case_effects.copy()
    ce["battery"] = ce.case_id.map(battery_map())
    rows = []
    for (contrast, metric), g in ce.groupby(["contrast", "metric"]):
        wide = g.pivot(index="case_id", columns="arm", values="delta")
        series = {arm: wide[arm] for arm, *_ in ARMS if arm in wide}
        if len(series) == 2:
            series["talkie_minus_qwen"] = wide["talkie"] - wide["qwen1930"]
        for arm, values in series.items():
            values = values.dropna()
            if len(values) < 3:
                continue
            mean, lo, hi = case_bootstrap(values, rng)
            rows.append(dict(arm=arm, contrast=contrast, metric=metric, mean_delta=mean,
                             ci_low=lo, ci_high=hi, n_cases=len(values)))
    return pd.DataFrame(rows)


def study2_trajectories(coded):
    rng = np.random.default_rng(1930)
    coded = coded.assign(battery=coded.case_id.map(battery_map()))
    panels = [("recommendation", "Recommends B"), ("equalized_preference", "Prefers B at equal payoff")]
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.4), sharey=True)
    for row, (metric, ylabel) in enumerate(panels):
        for col, (bat, title) in enumerate(BATTERIES):
            ax = axes[row, col]
            for k, (arm, label, color, marker) in enumerate(ARMS):
                stats = []
                for cond in CONDITIONS:
                    sub = coded[(coded.arm == arm) & (coded.battery == bat) & (coded.condition == cond)]
                    stats.append(case_bootstrap(sub.groupby("case_id")[metric].mean(), rng))
                mean, lo, hi = map(np.array, zip(*stats))
                x = np.arange(len(CONDITIONS)) + (k - 0.5) * 0.12
                ax.errorbar(x, mean, yerr=[mean - lo, hi - mean], color=color, marker=marker,
                            capsize=2, lw=1.4, ms=5, label=label)
            ax.set_xticks(range(len(CONDITIONS)), COND_LABELS, fontsize=8)
            if row == 0:
                ax.set_title(title, fontsize=9)
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=8)
            ax.set_ylim(-0.02, 1.02)
            ax.axhline(0.5, color="0.8", lw=0.6, zorder=0)
            ax.grid(axis="y", alpha=0.25)
    axes[0, 0].legend(fontsize=7, loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "trajectories_abc.pdf")
    plt.close(fig)


def study2_contrasts(effects):
    fig, axes = plt.subplots(1, 3, figsize=(10, 4.4), sharey=True)
    for ax, (metric, title) in zip(axes, METRICS):
        for k, (arm, label, color, marker) in enumerate(ARMS):
            sub = effects[(effects.arm == arm) & (effects.metric == metric)].set_index("contrast")
            sub = sub.reindex(EXT_CONTRASTS)
            y = np.arange(len(EXT_CONTRASTS))[::-1] + (0.5 - k) * 0.25
            ax.errorbar(sub.mean_delta, y,
                        xerr=[sub.mean_delta - sub.ci_low, sub.ci_high - sub.mean_delta],
                        fmt=marker, color=color, capsize=2, ms=4, label=label)
        ax.axvline(0, color="0.4", lw=0.8)
        ax.axhline(len(EXT_CONTRASTS) - len(CONTRASTS) - 0.5, color="0.7", lw=0.6, ls=":")
        ax.set_yticks(np.arange(len(EXT_CONTRASTS))[::-1], EXT_LABELS, fontsize=7)
        ax.set_title(title, fontsize=9)
        ax.set_xlim(-0.5, 0.85)
        ax.grid(axis="x", alpha=0.25)
        ax.set_xlabel("paired difference (24-domain bootstrap 95% CI)", fontsize=7)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=7, loc="upper center", ncol=2, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "contrasts_abc.pdf")
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    trajectories(pd.read_csv(DATA / "coded_observations.csv"))
    contrasts(pd.read_csv(DATA / "effects.csv"))
    effects2 = study2_effects(pd.read_csv(DATA2 / "case_effects.csv"))
    effects2.to_csv(OUT / "study2_effects.csv", index=False)
    study2_trajectories(pd.read_csv(DATA2 / "coded_observations.csv"))
    study2_contrasts(effects2)


if __name__ == "__main__":
    main()
