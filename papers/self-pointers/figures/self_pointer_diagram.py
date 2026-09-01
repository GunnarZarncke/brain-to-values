"""Matched counterfactual decision schematic for the self-pointers paper."""

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parent
OUT_STEM = ROOT / "self_pointer_counterfactual"


def draw_panel(ax, self_body, panel_label):
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.35, 2.05)
    ax.axis("off")

    frame = FancyBboxPatch(
        (0.35, 0.15),
        9.3,
        1.55,
        boxstyle="round,pad=0.04,rounding_size=0.16",
        linewidth=1.0,
        edgecolor="0.35",
        facecolor="0.97",
        transform=ax.transData,
        clip_on=False,
    )
    ax.add_patch(frame)

    xA, xT, xB, y = 2.5, 5.0, 7.5, 0.72

    for x, label in [(xA, "A"), (xB, "B")]:
        ax.add_patch(Circle((x, y), 0.22, fill=False, linewidth=1.4, edgecolor="0.15"))
        ax.text(x, 0.34, label, ha="center", va="top", fontsize=10)

    ax.plot(
        xT,
        y,
        marker="*",
        markersize=15,
        markerfacecolor="white",
        markeredgewidth=1.4,
        color="0.15",
    )
    ax.text(xT, 1.18, "target", ha="center", va="bottom", fontsize=9.5)

    sx = xA if self_body == "A" else xB
    ax.add_patch(Circle((sx, y), 0.36, fill=False, linewidth=1.8, edgecolor="0.15"))
    ax.text(sx, 1.52, r"hidden $S$", ha="center", va="bottom", fontsize=9.5, color="0.25")
    ax.plot([sx, sx], [1.48, y + 0.40], linestyle="--", linewidth=0.9, color="0.35")

    if self_body == "A":
        start = (xA + 0.38, y)
        end = (xT - 0.38, y)
        action = "RIGHT"
    else:
        start = (xB - 0.38, y)
        end = (xT + 0.38, y)
        action = "LEFT"

    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=1.7,
            color="0.15",
        )
    )

    ax.text(0.45, 1.88, panel_label, ha="left", va="center", fontsize=11, fontweight="bold")
    ax.text(
        5.0,
        -0.12,
        rf"identical final vision  $\Rightarrow$  label = {action}",
        ha="center",
        va="center",
        fontsize=9.5,
    )


def main():
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.35), constrained_layout=True)
    draw_panel(axes[0], "A", r"(a) $S=0$")
    draw_panel(axes[1], "B", r"(b) $S=1$")

    fig.savefig(f"{OUT_STEM}.pdf", bbox_inches="tight")
    fig.savefig(f"{OUT_STEM}.png", dpi=240, bbox_inches="tight")
    fig.savefig(f"{OUT_STEM}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_STEM}.pdf")


if __name__ == "__main__":
    main()
