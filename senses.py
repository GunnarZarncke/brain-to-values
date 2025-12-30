import numpy as np
import matplotlib.pyplot as plt
import math

# --- Data (illustrative)
stages = ["Peripheral\nTransduction", "Primary Maps\n(V1, A1, S1...)", "Higher Processing\n(V2-4, A2, S2...)", "Cortical Integration\n(inputs to PPC, STS, IT)"]
x_stage = np.array([0.0, 1.0, 2.0, 3.0])

def format_scientific(value: float) -> str:
    exp = int(math.floor(math.log10(value)))
    mantissa = value / (10 ** exp)
    if abs(mantissa - round(mantissa)) < 1e-6:
        mantissa = round(mantissa)
    return f"{mantissa}×10^{exp}"


sense_models = {
    "Vision":        [1e7, 1e6, 1e5, 1e4],
    "Hearing":       [1e5, 3e4, 1e4, 3e3],
    "Somatosensory": [1e5, 3e4, 1e4, 3e3],
    "Proprioception":[3e4, 1e4, 3e3, 1e3],
    "Olfaction":     [1e4, 3e3, 1e3, 3e2],
    "Gustation":     [1e3, 3e2, 1e2, 3e1],
    "Vestibular":    [3e3, 1e3, 3e2, 1e2],
}

# Shared downstream after merge
shared_stages = ["", "Cortical Integration\n(PPC, STS, IT)", "Feature selection\n(IPS, AI,...)", "Global Workspace\n(PFC, ACC, Thal)", "Reportable\n(mPFC, INS, PPC)"]
x_shared = np.array([3.0, 4, 5, 6, 7])
shared_bits = np.array([
    1.0,
    sum(v[-1] for v in sense_models.values()),  # integrated pool (sum of "Central" per-sense)
    1e3,                                        # feature selection output (illustrative)
    3e2,                                        # workspace ignition output (illustrative)
    5e1,                                        # serial stabilization (illustrative)
])

# --- Width mapping: half-width ∝ normalized log10(bits/s) ---
def log10_bits(bits: float) -> float:
    return math.log10(bits)

# Per-sense widths at the four sensory stages
sense_logw = {s: np.array([log10_bits(b) for b in vals]) for s, vals in sense_models.items()}
all_logw = np.concatenate([v for v in sense_logw.values()] + [np.array([log10_bits(b) for b in shared_bits])])
w_min, w_max = float(all_logw.min()), float(all_logw.max())

#def scale_halfwidth(wval: float, lo: float = 0.045, hi: float = 0.35) -> float:
#    t = 0.0 if w_max == w_min else (wval - w_min) / (w_max - w_min)
#    return lo + t * (hi - lo)

def scale_halfwidth(v): return v/20.0

sense_half = {s: np.array([scale_halfwidth(v) for v in ws]) for s, ws in sense_logw.items()}

# Shared funnel half-widths (use a slightly larger range so it reads as a "main channel")
shared_half = np.array([scale_halfwidth(log10_bits(b)) for b in shared_bits])

# --- Vertical layout: separate ribbons on the left ---
"""Mapping from senses to requested colors (case-insensitive keys)."""

sense_colors = {
    "Vision": "blue",
    "Hearing": "violet",
    "Somatosensory": "pink",
    "Proprioception": "yellow",
    "Olfaction": "brown",
    "Gustation": "red",
    "Vestibular": "green",
}
senses = sorted(
    sense_models.keys(),
    key=lambda s: sense_models[s][0],
)
gap = 0.5   # increase to keep ribbons distinct (adjust to taste)
base_centers = np.linspace(0, (len(senses) - 1) * gap, len(senses))
base_centers = base_centers - base_centers.mean()  # center around y=0

# --- Merge geometry (Central -> Integrated) ---
# Merge happens over x in [3.0, 4.2]. Centerlines converge smoothly toward y=0.
x_merge0, x_merge1 = x_shared[0], x_shared[1]
x_merge = np.linspace(x_merge0, x_merge1, 80)
merge_t = (x_merge - x_merge0) / (x_merge1 - x_merge0)
smoothstep = merge_t * merge_t * (3 - 2 * merge_t)                  # smooth 0->1
centers_merge = np.outer(1 - smoothstep, base_centers)  # shape: (len(x_merge), n_senses)

# --- Plot ---
plt.figure(figsize=(12, 6))
ax = plt.gca()

# Per-sense tapered ribbons from Peripheral->Central, plus colored merge segment
for i, sense in enumerate(senses):
    c0 = base_centers[i]
    hw = sense_half[sense]

    color = sense_colors[sense]
    # Sensory stages ribbon (tapered between stage points)
    top = c0 + hw
    bot = c0 - hw
    ax.fill(
        np.concatenate([x_stage, x_stage[::-1]]),
        np.concatenate([top, bot[::-1]]),
        alpha=0.65,
        label=sense,
        color=color,
    )
    for stage_idx, stage_x in enumerate(x_stage):
        stage_bits = sense_models[sense][stage_idx]
        ax.text(
            stage_x,
            c0,
            format_scientific(stage_bits),
            ha="center",
            va="center",
            fontsize="small",
            color="black",
        )

    # Colored merge segment (keep color identity through the merge)
    hw_central = hw[-1]  # hold constant in the merge region (simple + visually stable)
    top_m = centers_merge[:, i] + hw_central
    bot_m = centers_merge[:, i] - hw_central
    ax.fill(
        np.concatenate([x_merge, x_merge[::-1]]),
        np.concatenate([top_m, bot_m[::-1]]),
        alpha=0.65,
        color=color,
    )

# Shared downstream funnel (gray)
x_post = np.linspace(x_shared[1], x_shared[-1], 300)
hw_post = np.interp(x_post, x_shared[1:], shared_half[1:])
ax.fill(
    np.concatenate([x_post, x_post[::-1]]),
    np.concatenate([hw_post, -hw_post[::-1]]),
    alpha=0.30,
    color="gray",
    label="Merged stream"
)

for stage_idx, stage_x in enumerate(x_shared):
    stage_bits = shared_bits[stage_idx]
    ax.text(
        stage_x,
        0.0,
        format_scientific(stage_bits),
        ha="center",
        va="center",
        fontsize="small",
        color="black",
    )

# Stage markers + labels
stage_tick_padding = 0.08 * (ax.get_ylim()[1] - ax.get_ylim()[0])
stage_label_y = ax.get_ylim()[0] - stage_tick_padding
ax.set_xticks(list(x_stage))
ax.set_xticklabels([""] * len(x_stage))
for x, name in zip(x_stage, stages):
    ax.axvline(x, linewidth=1, alpha=0.25)
    ax.text(x, stage_label_y, name, ha="center", va="top")
ax.set_ylim(ax.get_ylim()[0] - stage_tick_padding, ax.get_ylim()[1])

for x, name in zip(x_shared, [""] + shared_stages[1:]):
    ax.axvline(x, linewidth=1, alpha=0.25)
    ax.text(x, stage_label_y, name, ha="center", va="top")

ax.set_title("Data reduction during sensory processing from raw sense to reportable sensation (width ∝ log10(bits/s) illustrative)")
#ax.set_xlabel("Processing stages (left → right)")
ax.set_yticks([])
ax.set_xlim(-0.2, x_shared[-1] + 0.2)
ax.legend(loc="upper right", frameon=False, ncol=2)

plt.tight_layout()
plt.savefig("senses.png", dpi=120, bbox_inches="tight")
plt.show()

