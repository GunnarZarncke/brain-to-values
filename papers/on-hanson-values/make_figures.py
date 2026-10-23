#!/usr/bin/env python3
"""Figures, tables and in-text numbers for on-hanson-values.tex.

Study 1: 72 domains in nine nested worlds (abc_extended_r3, d_extended_r3, e_all_r3).
Study 2: the pre-registered moral/practical x real/fictional test on 40 domains (e_five_r3).
Writes figures/*.pdf, tables/*.tex and tables/numbers.tex; the prose cites numbers as \\n{key}.
Until the battery E extension is merged into e_all_r3, its five worlds are read from e_five_r3.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from historical_value_experiment import (CASES_BATTERY_A, CASES_BATTERY_B, CASES_BATTERY_C,  # noqa: E402
                                         CASES_BATTERY_D, CASES_BATTERY_E, load_jsonl)

FIG, TAB = HERE / "figures", HERE / "tables"
E_DIR = HERE / "e_all_r3" if (HERE / "e_all_r3" / "coded_observations.csv").exists() else HERE / "e_five_r3"
RUNS = [HERE / "abc_extended_r3", HERE / "d_extended_r3", E_DIR]
N_BOOT = 5000
ARMS = [("talkie", "Talkie", "#8c2d04", "o"), ("qwen1930", "Qwen", "#2b6c9e", "s")]
SETS = [("favorB", "Facts favor the modern practice"), ("favorA", "Facts favor the older practice"),
        ("hanson", "Changes Hanson names, real evidence"), ("balanced", "Balanced facts")]
CELLS = [("moral_real", "moral, real"), ("practical_real", "practical, real"),
         ("moral_fictional", "moral, fictional"), ("practical_fictional", "practical, fictional")]
WORLDS = ["A_1930", "B_transplant", "C_elites", "C_explained", "C_opposite", "C_moral",
          "D_peers", "D_media", "D_full"]
SPECS = {"eco": ("B_transplant", "A_1930"), "eli": ("C_elites", "B_transplant"),
         "soc": ("D_full", "C_elites"), "tot": ("D_full", "A_1930"),
         "opacity": ("C_elites", "C_explained"), "direction": ("C_elites", "C_opposite"),
         "pullback": ("B_transplant", "C_opposite"),
         "moralclaim": ("C_moral", "C_elites"), "peers": ("D_peers", "C_elites"),
         "media": ("D_media", "D_peers"), "inst": ("D_full", "D_media")}
NUM = {}
rng = np.random.default_rng(1930)


def domains():
    val = pd.read_csv(HERE / "domain_validation.csv").set_index("domain").verdict
    rows = []
    for bat, cases, dset in [("A", CASES_BATTERY_A, "favorB"), ("B", CASES_BATTERY_B, "favorB"),
                             ("C", CASES_BATTERY_C, "favorA"), ("D", CASES_BATTERY_D, "hanson"),
                             ("E", CASES_BATTERY_E, "balanced")]:
        for c in cases:
            cell = c.get("cell", "")
            if bat == "D":
                valid = "real"
            elif bat == "E":
                valid = "real" if cell.endswith("real") else "fictional"
            else:
                v = val[c["id"]]
                valid = ("real" if v.startswith("real") else
                         "partial" if v.startswith(("partial", "mostly", "close")) else "notreal")
            rows.append(dict(case_id=c["id"], battery=bat, set=dset, cell=cell, valid=valid))
    return pd.DataFrame(rows).set_index("case_id")


def observations(meta):
    obs = pd.concat([pd.read_csv(d / "coded_observations.csv") for d in RUNS], ignore_index=True)
    obs["elicited_unstated"] = (obs.elicited_moral_label == "not_stated").astype(float)
    return obs.join(meta, on="case_id")


def case_deltas(obs, metric="recommendation"):
    t = obs.pivot_table(index=["arm", "case_id", "rep"], columns="condition", values=metric, aggfunc="mean")
    return pd.DataFrame({k: (t[h] - t[lo]).groupby(level=[0, 1]).mean()
                         for k, (h, lo) in SPECS.items() if h in t and lo in t})


def boot(values):
    v = pd.Series(values).dropna().to_numpy()
    if len(v) < 2:
        return (v.mean() if len(v) else np.nan), np.nan, np.nan, len(v)
    draws = rng.choice(v, size=(N_BOOT, len(v)), replace=True).mean(axis=1)
    return v.mean(), *np.percentile(draws, [2.5, 97.5]), len(v)


def sgn(x):
    return "---" if not np.isfinite(x) else f"{x:+.2f}".replace("+0.00", "0.00").replace("-0.00", "0.00")


def ci(m, lo, hi):
    return f"${sgn(m)}$ $[{sgn(lo)},{sgn(hi)}]$" if np.isfinite(lo) else f"${sgn(m)}$"


def rho(d):
    both = d[["eco", "tot"]].dropna()
    return 1 - both.eco.mean() / both.tot.mean() if len(both) and abs(both.tot.mean()) >= 0.10 else np.nan


def num(key, value, fmt="{:.2f}"):
    NUM[key] = value if isinstance(value, str) else ("---" if not np.isfinite(value) else fmt.format(value))


def contrast_cell(d, k, key=None):
    m, lo, hi, n = boot(d[k]) if k in d else (np.nan,) * 4
    if key:
        num(key, sgn(m))
        num(key + "-ci", f"[{sgn(lo)},{sgn(hi)}]" if np.isfinite(lo) else "")
        num(key + "-n", str(n))
    return ci(m, lo, hi)


# ---------------------------------------------------------------- Study 1
def core_table(obs, cd):
    meta = domains()
    groups = [(label, s, cd.index.get_level_values(1).map(meta.set) == s) for s, label in SETS]
    groups += [("All domains", "all", np.ones(len(cd), bool)),
               ("Real changes since 1930", "real", cd.index.get_level_values(1).map(meta.valid) == "real")]
    lines = []
    for arm, name, *_ in ARMS:
        lines.append(f"\\multicolumn{{6}}{{@{{}}l}}{{\\emph{{{name}}}}}\\\\")
        for label, key, mask in groups:
            d = cd[mask & (cd.index.get_level_values(0) == arm)]
            n = d.index.get_level_values(1).nunique()
            cells = [contrast_cell(d, k, f"{arm}-{k}-{key}") for k in ("eco", "eli", "soc", "tot")]
            r = rho(d)
            num(f"{arm}-rho-{key}", r)
            lines.append(f"{label} ({n}) & " + " & ".join(cells) + f" & {'---' if not np.isfinite(r) else f'{r:.2f}'}\\\\")
        lines.append("\\midrule")
    (TAB / "core.tex").write_text("\n".join(lines[:-1]) + "\n")


def bundle_table(cd):
    rows = [("eli", "$C-B$, elites adopt the modern practice"), ("pullback", "$B-C_{\\text{opp}}$, elites keep the older practice"),
            ("direction", "$C-C_{\\text{opp}}$, direction of elite example"), ("opacity", "$C-C_{\\text{expl}}$, reasons withheld"),
            ("moralclaim", "$C_{\\text{moral}}-C$, elites call it good character"), ("peers", "$D_{\\text{peers}}-C$, peers"),
            ("media", "$D_{\\text{media}}-D_{\\text{peers}}$, favorable media"), ("inst", "$D-D_{\\text{media}}$, institutions")]
    t, q = cd.loc["talkie"], cd.loc["qwen1930"]
    diff = (t - q.reindex(t.index))
    lines = []
    for k, label in rows:
        cells = [contrast_cell(t, k, f"b-talkie-{k}"), contrast_cell(q, k, f"b-qwen1930-{k}"),
                 contrast_cell(diff, k, f"b-diff-{k}")]
        lines.append(f"{label} & " + " & ".join(cells) + f" & {t[k].notna().sum()}\\\\")
    (TAB / "bundle.tex").write_text("\n".join(lines) + "\n")
    for s, _ in SETS:
        ids = domains().index[domains().set == s]
        for k in ("direction", "pullback", "peers", "moralclaim"):
            contrast_cell(t[t.index.isin(ids)], k, f"talkie-{k}-{s}")
            contrast_cell(q[q.index.isin(ids)], k, f"qwen1930-{k}-{s}")


def adapt_table(obs, cd):
    meta = domains()
    groups = [("Facts favor the older practice", meta.set == "favorA"),
              ("Changes that did not happen as framed", meta.valid == "notreal"),
              ("Fictional moral practices", meta.cell == "moral_fictional"),
              ("Fictional practical practices", meta.cell == "practical_fictional"),
              ("Real practical changes", meta.cell == "practical_real"),
              ("Real changes since 1930", meta.valid == "real")]
    keys = ["favorA", "notreal", "moralfict", "practfict", "practreal", "real"]
    lvl = obs.groupby(["arm", "case_id", "condition"]).recommendation.mean()
    lines = []
    for (label, mask), key in zip(groups, keys):
        ids = meta.index[mask]
        cells = []
        for arm, *_ in ARMS:
            d = cd.loc[arm]
            d = d[d.index.isin(ids)]
            a = lvl.loc[arm].unstack().reindex(ids)
            lv = {w: a[w].mean() if w in a else np.nan for w in ("A_1930", "D_peers", "D_full")}
            for w, v in lv.items():
                num(f"ad-{arm}-{w.replace('_', '')}-{key}", v)
            cells += [f"{lv['A_1930']:.2f}", f"{lv['D_full']:.2f}".replace("nan", "---"),
                      contrast_cell(d, "eco", f"ad-{arm}-eco-{key}"), contrast_cell(d, "tot", f"ad-{arm}-tot-{key}")]
            contrast_cell(d, "peers", f"ad-{arm}-peers-{key}")
            pt = (lvl.loc[arm].unstack().reindex(ids))
            if "D_peers" in pt:
                totp = obs[(obs.arm == arm) & obs.case_id.isin(ids)].pivot_table(
                    index=["case_id", "rep"], columns="condition", values="recommendation")
                contrast_cell((totp.D_peers - totp.A_1930).groupby(level=0).mean().to_frame("x"), "x",
                              f"ad-{arm}-totpeers-{key}")
        lines.append(f"{label} ({len(ids)}) & " + " & ".join(cells) + "\\\\")
    (TAB / "adapt.tex").write_text("\n".join(lines) + "\n")


def levels_table(obs):
    lvl = obs.groupby(["arm", "set", "condition", "case_id"]).recommendation.mean().groupby(level=[0, 1, 2]).mean()
    lines = []
    for arm, name, *_ in ARMS:
        lines.append(f"\\multicolumn{{10}}{{@{{}}l}}{{\\emph{{{name}}}}}\\\\")
        for s, label in SETS + [("all", "All domains")]:
            if s == "all":
                row = obs[obs.arm == arm].groupby(["condition", "case_id"]).recommendation.mean().groupby(level=0).mean()
            else:
                row = lvl.loc[(arm, s)]
            lines.append(f"{label} & " + " & ".join(f"{row[w]:.2f}" if w in row else "---" for w in WORLDS) + "\\\\")
            for w in WORLDS:
                if w in row:
                    num(f"{arm}-level-{s}-{w.replace('_', '')}", row[w])
        lines.append("\\midrule")
    (TAB / "levels.tex").write_text("\n".join(lines[:-1]) + "\n")


def battery_table(cd):
    meta = domains()
    lines = []
    for bat, label in [("A", "A"), ("B", "B"), ("C", "C"), ("D", "D"), ("E", "E")]:
        ids = meta.index[meta.battery == bat]
        for arm, name, *_ in ARMS:
            d = cd.loc[arm]
            d = d[d.index.isin(ids)]
            r = rho(d)
            cells = [contrast_cell(d, k) for k in ("eco", "eli", "soc", "tot", "direction", "peers")]
            lines.append(f"{label if arm == 'talkie' else ''} & {name} & " + " & ".join(cells)
                         + f" & {'---' if not np.isfinite(r) else f'{r:.2f}'}\\\\")
    (TAB / "batteries.tex").write_text("\n".join(lines) + "\n")


def moral_table(obs):
    worlds = [("A_1930", "A"), ("B_transplant", "B"), ("C_elites", "C"), ("D_peers", "$D_{\\text{peers}}$"), ("D_full", "D")]
    metrics = [("spontaneous_any_moral", "any"), ("elicited_unstated", "unstated"), ("elicited_moral", "mdir")]
    lines = []
    for w, wl in worlds:
        cells = []
        for arm, *_ in ARMS:
            sub = obs[(obs.arm == arm) & (obs.condition == w)]
            for m, key in metrics:
                v = sub.groupby("case_id")[m].mean().mean()
                num(f"{arm}-{key}-{w.replace('_', '')}", v)
                cells.append(f"{v:.2f}")
        lines.append(f"{wl} & " + " & ".join(cells) + "\\\\")
    (TAB / "moral.tex").write_text("\n".join(lines) + "\n")


def reasons_table():
    df = pd.read_csv(HERE / "cue_reasons.csv")
    df["reason"], df["absent"], df["rejected"] = df.use == "reason", df.use == "absent", df.use == "rejected"
    per = df.groupby(["run", "trial_id", "arm", "coder", "follows"]).agg(
        reason=("reason", "max"), absent=("absent", "min"), rejected=("rejected", "max"),
        judge=("judge_social", "first")).reset_index()
    w = per.pivot_table(index=["run", "trial_id", "follows"], columns="coder", values="reason").dropna()
    a, b = w.columns
    num("cue-agree", 100 * ((w[a] > 0) == (w[b] > 0)).mean(), "{:.0f}")
    use = df.pivot_table(index=["run", "trial_id", "cue"], columns="coder", values="use", aggfunc="first").dropna()
    num("cue-agree4", 100 * (use.iloc[:, 0] == use.iloc[:, 1]).mean(), "{:.0f}")
    num("cue-n", str(df.drop_duplicates(["run", "trial_id"]).shape[0]))
    lines = []
    for arm, name, *_ in ARMS:
        withc = per[(per.arm == arm) & (per.follows == "with")]
        against = per[(per.arm == arm) & (per.follows == "against")]
        n = withc.trial_id.nunique()
        r, ab, j = (withc.groupby("coder")[x].mean() for x in ("reason", "absent", "judge"))
        rej, ab2 = (against.groupby("coder")[x].mean() for x in ("rejected", "absent"))
        rng_ = lambda s: f"{100 * s.min():.0f}" if round(100 * s.min()) == round(100 * s.max()) else f"{100 * s.min():.0f}--{100 * s.max():.0f}"
        num(f"{arm}-cue-reason", rng_(r)); num(f"{arm}-cue-absent", rng_(ab)); num(f"{arm}-cue-judge", f"{100 * j.mean():.0f}")
        num(f"{arm}-cue-n", str(n)); num(f"{arm}-cue-rejected", rng_(rej)); num(f"{arm}-cue-absent-against", rng_(ab2))
        lines.append(f"{name} & {n} & {rng_(r)}\\% & {rng_(ab)}\\% & {100 * j.mean():.0f}\\% & "
                     f"{against.trial_id.nunique()} & {rng_(rej)}\\% & {rng_(ab2)}\\%\\\\")
    (TAB / "reasons.tex").write_text("\n".join(lines) + "\n")
    # Domain-level difference in giving the cue as a reason, averaged over coders
    obs_ids = pd.concat([pd.read_csv(d / "coded_observations.csv")[["trial_id", "case_id"]] for d in RUNS])
    x = per[per.follows == "with"].merge(obs_ids, on="trial_id").groupby(["arm", "case_id"]).reason.mean().unstack(0).dropna()
    m, lo, hi, n = boot(x.talkie - x.qwen1930)
    num("cue-diff", sgn(m)); num("cue-diff-ci", f"[{sgn(lo)},{sgn(hi)}]"); num("cue-diff-n", str(n))


# ---------------------------------------------------------------- Study 2
def study2_tables(cd):
    c = pd.read_csv(HERE / "e_five_r3" / "battery_e_contrasts.csv").set_index(["arm", "test"])
    rows = [("X_primary", "Primary: $(E-P)_{\\text{moral}}-(E-P)_{\\text{practical}}$"),
            ("X_without_telephone2", "Primary without the replaced domain"),
            ("E_moral_minus_practical", "Elite step $E$, moral $-$ practical"),
            ("P_moral_minus_practical", "Peer step $P$, moral $-$ practical"),
            ("R_moral_minus_practical", "Elite direction $R$, moral $-$ practical"),
            ("E_real_minus_fictional", "Elite step $E$, real $-$ fictional"),
            ("E_known_failed_minus_not", "$E$, fictional: known failed by 1930 $-$ not"),
            ("E_pre1931_minus_not", "$E$, fictional: in print before 1931 $-$ not")]
    lines = []
    for test, label in rows:
        cells = []
        for arm in ("talkie", "qwen1930"):
            r = c.loc[(arm, test)]
            cells.append(ci(r.estimate, r.ci_low, r.ci_high))
            num(f"s2-{arm}-{test.replace('_', '')}", sgn(r.estimate))
            num(f"s2-{arm}-{test.replace('_', '')}-ci", f"[{sgn(r.ci_low)},{sgn(r.ci_high)}]")
        lines.append(f"{label} & " + " & ".join(cells) + f" & {int(c.loc[('talkie', test), 'n_high'])}/{int(c.loc[('talkie', test), 'n_low'])}\\\\")
    r = c.loc[("talkie_minus_qwen", "X_primary")]
    num("s2-diff-X", sgn(r.estimate)); num("s2-diff-X-ci", f"[{sgn(r.ci_low)},{sgn(r.ci_high)}]")
    lines.append(f"Talkie $-$ Qwen, primary & \\multicolumn{{2}}{{c}}{{{ci(r.estimate, r.ci_low, r.ci_high)}}} & \\\\")
    (TAB / "study2_tests.tex").write_text("\n".join(lines) + "\n")
    cells = pd.read_csv(HERE / "e_five_r3" / "battery_e_cells.csv").set_index(["arm", "step", "cell"])
    lines = []
    for arm, name, *_ in ARMS:
        lines.append(f"\\multicolumn{{5}}{{@{{}}l}}{{\\emph{{{name}}}}}\\\\")
        for step, label in [("E", "Elite step $E=C-B$"), ("P", "Peer step $P=D_{\\text{peers}}-C$"),
                            ("R", "Elite direction $R=C-C_{\\text{opp}}$"), ("EP", "$E-P$")]:
            vals = []
            for cell, _ in CELLS:
                r = cells.loc[(arm, step, cell)]
                vals.append(ci(r["mean"], r.ci_low, r.ci_high))
                num(f"s2-{arm}-{step}-{cell.replace('_', '')}", sgn(r["mean"]))
                num(f"s2-{arm}-{step}-{cell.replace('_', '')}-ci", f"[{sgn(r.ci_low)},{sgn(r.ci_high)}]")
            lines.append(f"{label} & " + " & ".join(vals) + "\\\\")
        meta = domains()
        for k, label in [("pullback", "Pullback $B-C_{\\text{opp}}$ (exploratory)"),
                         ("moralclaim", "Moral claim $C_{\\text{moral}}-C$ (exploratory)")]:
            vals = []
            for cell, _ in CELLS:
                d = cd.loc[arm]
                vals.append(contrast_cell(d[d.index.isin(meta.index[meta.cell == cell])], k,
                                          f"s2-{arm}-{k}-{cell.replace('_', '')}"))
            lines.append(f"{label} & " + " & ".join(vals) + "\\\\")
        lines.append("\\midrule")
    (TAB / "study2_cells.tex").write_text("\n".join(lines[:-1]) + "\n")
    return cells


def study2_figure(cells):
    steps = [("E", "Elite step $E=C-B$"), ("P", "Peer step $P=D_{peers}-C$"), ("R", "Elite direction $R=C-C_{opp}$")]
    fig, axes = plt.subplots(1, 3, figsize=(10, 2.8), sharey=True)
    for ax, (step, title) in zip(axes, steps):
        for k, (arm, label, color, marker) in enumerate(ARMS):
            sub = cells.sort_index().loc[(arm, step)].reindex([c for c, _ in CELLS])
            y = np.arange(len(CELLS))[::-1] + (0.5 - k) * 0.25
            ax.errorbar(sub["mean"], y, xerr=[sub["mean"] - sub.ci_low, sub.ci_high - sub["mean"]],
                        fmt=marker, color=color, capsize=2, ms=5, label=label)
        ax.axvline(0, color="0.4", lw=0.8)
        ax.set_yticks(np.arange(len(CELLS))[::-1], [lbl for _, lbl in CELLS], fontsize=8)
        ax.set_title(title, fontsize=9)
        ax.set_xlim(-0.45, 1.1)
        ax.grid(axis="x", alpha=0.25)
        ax.set_xlabel("paired difference (10-domain bootstrap 95% CI)", fontsize=7)
    axes[0].legend(fontsize=7, loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "study2_cells.pdf")
    plt.close(fig)


def trajectories(obs):
    worlds = ["A_1930", "B_transplant", "C_elites", "D_peers", "D_full"]
    labels = ["A\n1930", "B\nfacts", "C\n+elites", "D$_{peers}$\n+peers", "D\n+media,\ninst."]
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.2), sharey=True)
    for ax, (s, title) in zip(axes, SETS):
        for k, (arm, label, color, marker) in enumerate(ARMS):
            stats = []
            for w in worlds:
                sub = obs[(obs.arm == arm) & (obs.set == s) & (obs.condition == w)]
                stats.append(boot(sub.groupby("case_id").recommendation.mean())[:3] if len(sub) else (np.nan,) * 3)
            mean, lo, hi = map(np.array, zip(*stats))
            x = np.arange(len(worlds)) + (k - 0.5) * 0.12
            ax.errorbar(x, mean, yerr=[mean - lo, hi - mean], color=color, marker=marker, capsize=2, lw=1.4, ms=5, label=label)
        n = obs[obs.set == s].case_id.nunique()
        ax.set_title(f"{title}\n({n} domains)", fontsize=8)
        ax.set_xticks(range(len(worlds)), labels, fontsize=7)
        ax.set_ylim(-0.02, 1.02)
        ax.axhline(0.5, color="0.8", lw=0.6, zorder=0)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("recommends the modern practice\n(older 0, mixed .5, modern 1)", fontsize=8)
    axes[0].legend(fontsize=7, loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "trajectories.pdf")
    plt.close(fig)


def domains_e_table():
    for part in ("real", "fictional"):
        lines = []
        for cell, label in CELLS:
            if cell.endswith(part):
                lines.append(f"\\multicolumn{{3}}{{@{{}}l}}{{\\emph{{{label}}}}}\\\\")
                lines += [f"{c['id']} & {c['a']} & {c['b']}\\\\" for c in CASES_BATTERY_E if c["cell"] == cell]
        (TAB / f"domains_e_{part}.tex").write_text("\n".join(lines) + "\n")


def answer_lengths(obs):
    for arm, *_ in ARMS:
        rows = [r for d in RUNS for r in load_jsonl(d / "interviews.jsonl")
                if r.get("status") == "ok" and r["arm"] == arm]
        num(f"cap-{arm}-n", str(len(rows)))
        num(f"cap-{arm}-hit", str(sum(r["answers"][0].get("finish_reason") == "length" for r in rows)))
        num(f"words-{arm}", f"{np.median([len(r['answers'][0]['text'].split()) for r in rows]):.0f}")
        num(f"unclear-{arm}", 100 * (obs[obs.arm == arm].recommendation_label == "unclear").mean(), "{:.0f}")


def counts(obs):
    meta = domains()
    num("n-domains", str(len(meta)))
    num("n-interviews", str(len(obs)))
    for s, _ in SETS:
        num(f"n-{s}", str((meta.set == s).sum()))
    for v in ("real", "partial", "notreal", "fictional"):
        num(f"n-{v}", str((meta.valid == v).sum()))
    planned = len(meta) * len(WORLDS) * 3 * 2
    num("n-planned", str(planned))
    num("e-source", E_DIR.name)


def main():
    FIG.mkdir(exist_ok=True)
    TAB.mkdir(exist_ok=True)
    meta = domains()
    obs = observations(meta)
    cd = case_deltas(obs)
    counts(obs)
    answer_lengths(obs)
    domains_e_table()
    core_table(obs, cd)
    bundle_table(cd)
    adapt_table(obs, cd)
    levels_table(obs)
    battery_table(cd)
    moral_table(obs)
    reasons_table()
    study2_figure(study2_tables(cd))
    trajectories(obs)
    with open(TAB / "numbers.tex", "w") as f:
        f.write("% Generated by make_figures.py; do not edit.\n")
        for k, v in sorted(NUM.items()):
            f.write(f"\\expandafter\\def\\csname n@{k}\\endcsname{{{v}}}\n")
    print(f"{len(NUM)} numbers; battery E from {E_DIR.name}; {len(obs)} coded interviews")


if __name__ == "__main__":
    main()
