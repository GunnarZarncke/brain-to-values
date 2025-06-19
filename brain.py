"""
Dynamic brain-value diagram (matplotlib version).

Usage
-----
1. Drop a sagittal brain image somewhere (e.g. brain.png).
2. Adjust CONFIG block if you want different radii, fonts, or the mapping itself.
3. run  python brain_values_plot.py
"""

import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import textwrap
from matplotlib.axes import Axes
import matplotlib as mpl   
import webbrowser
# Use standard matplotlib font configuration
mpl.rcParams['font.family'] = 'sans-serif'

# Global variable for click detection
value_positions = {}

# ── CONFIG ────────────────────────────────────────────────────────────────────
brain_img_path = "brain.png"            # local file or URL
figsize         = (10, 10)              # inches
brain_img_scale = 0.5                   # fraction of figure width
region_radius   = 0.18                  # fraction of figure width
value_radius    = 0.47                  # fraction of figure width
region_fontsz   = 9
value_fontsz    = 12                    # 50 % larger than before
summary_fontsz  = 12                    # ditto
angle_offset    = -90                   # start at 12 o’clock

regions = [
    "Superior Colliculus", "Periaqueductal Gray", "Ventral Tegmental Area",
    "Locus Coeruleus", "Nucleus Accumbens", "Dorsal Raphe Nucleus",
    "Hypothalamus", "Amygdala", "Anterior Cingulate Cortex",
    "Insula", "Septal Nuclei", "Primary Visual Cortex"
]

# intrinsic values mapped to brain region (abbreviations are unnecessary here)
values = {
    "Superior Colliculus": ["Nature", "Truth"],
    "Periaqueductal Gray": ["Spirituality", "Freedom", "Non-suffering", "Protection"],
    "Ventral Tegmental Area": ["Learning", "Legacy", "Diversity"],
    "Locus Coeruleus": ["Achievement"],
    "Nucleus Accumbens": ["Happiness"],
    "Dorsal Raphe Nucleus": ["Pleasure"],
    "Hypothalamus": ["Longevity", "Reputation", "Caring"],
    "Amygdala": ["Virtue", "Respect"],
    "Anterior Cingulate Cortex": ["Justice"],
    "Insula": ["Fairness", "Purity"],
    "Septal Nuclei": ["Loyalty"],
    "Primary Visual Cortex": ["Beauty"]
}

# References for clickable links
references = {
  "Nature": {
    "title": "Bringing nature indoors: characterizing the unique contribution of fractal structure and the effects of Euclidean context on perception of fractal patterns",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10488714/"
  },
  "Beauty": {
    "title": "Unit 4 Overview: Symmetry Research in Neuroaesthetics",
    "url": "https://link.springer.com/chapter/10.1007/978-3-031-42323-9_11"
  },
  "Purity": {
    "title": "Core, social and moral disgust are bounded: A review on behavioral and neural bases of repugnance in clinical disorders",
    "url": "https://pubmed.ncbi.nlm.nih.gov/28506923/"
  },
  "Spirituality": {
    "title": "Brain mechanisms in religion and spirituality: An integrative predictive processing framework",
    "url": "https://pubmed.ncbi.nlm.nih.gov/28041787/"
  },
  "Truth": {
    "title": "Neural substrates of appetitive and aversive prediction error",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7933120/"
  },
  "Learning": {
    "title": "Dopamine reward prediction error coding",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4826767/"
  },
  "Achievement": {
    "title": "Locus coeruleus: a new look at the blue spot",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8991985/"
  },
  "Freedom": {
    "title": "Serotonin, the periaqueductal gray and panic",
    "url": "https://pubmed.ncbi.nlm.nih.gov/15225969/"
  },
  "Happiness": {
    "title": "Feeding releases endogenous opioids in humans",
    "url": "https://www.jneurosci.org/content/37/34/8284"
  },
  "Pleasure": {
    "title": "Neuronal Reward and Decision Signals: From Theories to Data",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4491543/"
  },
  "Non-suffering": {
    "title": "The contribution of periaqueductal gray in the regulation of pain, defensive and aggressive behaviors, anxiety and depression",
    "url": "https://www.frontiersin.org/articles/10.3389/fnins.2024.1380171/full"
  },
  "Longevity": {
    "title": "Understanding the aging hypothalamus, one cell at a time",
    "url": "https://www.cell.com/trends/neurosciences/fulltext/S0166-2236(22)00192-8"
  },
  "Virtue": {
    "title": "Integrative Moral Judgment: Dissociating the Roles of the Amygdala and Ventromedial Prefrontal Cortex",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6608126/"
  },
  "Reputation": {
    "title": "Oxytocin and social motivation",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3185363/"
  },
  "Legacy": {
    "title": "Pharmacological Modulation of Temporal Discounting: A Systematic Review",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10093895/"
  },
  "Loyalty": {
    "title": "Oxytocin and social relationships: From attachment to bond disruption",
    "url": "https://www.researchgate.net/publication/319133187_Oxytocin_and_Social_Relationships_From_Attachment_to_Bond_Disruption"
  },
  "Justice": {
    "title": "Neural Correlates of Advantageous and Disadvantageous Inequity in Sharing Decisions",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4169616/"
  },
  "Fairness": {
    "title": "Fairness and Inequity Aversion",
    "url": "https://www.researchgate.net/publication/283677766_Fairness_and_Inequity_Aversion"
  },
  "Diversity": {
    "title": "Dopamine Modulates Novelty Seeking Behavior During Decision Making",
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5861725/"
  },
  "Respect": {
    "title": "Convergence of oxytocin and dopamine signalling in neuronal circuits",
    "url": "https://www.sciencedirect.com/science/article/pii/S0149763424001441"
  },
  "Caring": {
    "title": "Oxytocin, Dopamine, and Opioid Interactions Underlying Pair Bonding",
    "url": "https://academic.oup.com/endo/article/162/2/bqaa223/6046188"
  },
  "Protection": {
    "title": "Periaqueductal gray activates antipredatory neural responses in the presence of predator threat",
    "url": "https://elifesciences.org/articles/88733"
  }
}

# five-word explanations (two lines are handled later)
summaries = {
    ("Superior Colliculus","Nature"):   "fractal visuals evoke nature affinity",
    ("Superior Colliculus","Truth"):    "errors may spur truth-seeking",
    ("Periaqueductal Gray","Spirituality"): "low error states may create unity",
    ("Periaqueductal Gray","Freedom"):  "escape circuits motivate freedom",
    ("Periaqueductal Gray","Non-suffering"): "homeostasis underpins comfort",
    ("Periaqueductal Gray","Protection"):    "threat response drives defense",
    ("Ventral Tegmental Area","Learning"):   "dopamine tags surprise events",
    ("Ventral Tegmental Area","Legacy"):     "low discount may foster identity",
    ("Ventral Tegmental Area","Diversity"):  "novelty-seeking encourages variety",
    ("Locus Coeruleus","Achievement"):       "norepinephrine may signal mastery",
    ("Nucleus Accumbens","Happiness"):       "opioids encode basic need fulfillment",
    ("Dorsal Raphe Nucleus","Pleasure"):     "serotonin rewards basic needs",
    ("Hypothalamus","Longevity"):            "vital sensing drives survival value",
    ("Hypothalamus","Reputation"):           "social bonding signals status",
    ("Hypothalamus","Caring"):               "oxytocin circuits foster nurturing",
    ("Amygdala","Virtue"):                   "acceptance learning may shape morality",
    ("Amygdala","Respect"):                  "hierarchy detection may inform respect",
    ("Anterior Cingulate Cortex","Justice"): "inequity distress flags violations",
    ("Insula","Fairness"):                   "inequity signals underpin fairness",
    ("Insula","Purity"):                     "disgust reflex grounds purity",
    ("Septal Nuclei","Loyalty"):             "attachment bonds build loyalty",
    ("Primary Visual Cortex","Beauty"):      "visual symmetry triggers beauty"
}


# ──────────────────────────────────────────────────────────────────────────────


def polar_to_cart(r, theta_deg):
    """Convenience: polar-to-Cartesian with center at (0,0)."""
    th = math.radians(theta_deg)
    return r * math.cos(th), r * math.sin(th)


def compute_positions(n, radius):
    """Evenly space n labels on a circle of given radius; return list of (x,y)."""
    return [polar_to_cart(radius, angle_offset + i * 360 / n) for i in range(n)]


def split_summary(text):
    """Half-split summary into ≤ 2 roughly equal lines."""
    words = text.split()
    half = len(words) // 2
    return " ".join(words[:half]), " ".join(words[half:])

def label_on_line(ax: Axes,
                  x1: float, y1: float,
                  x2: float, y2: float,
                  text: str,
                  *,
                  fontsize: int = 12,
                  max_words_per_line: int = 3,
                  offset_frac: float = 0.01,
                  color: str = "#555"):
    """
    Draw `text` along the connector (x1,y1)–(x2,y2), rotated with the line.
    For single line text, place it on one side of the line.
    For two lines of text, place one line on each side of the line.
    """
    # --- wrap text (≤ max_words_per_line per row) ---------------------------
    words = text.split()
    wrapped_lines = [" ".join(words[i:i + max_words_per_line])
                     for i in range(0, len(words), max_words_per_line)]

    # --- mid-point and angle ------------------------------------------------
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1
    angle = math.degrees(math.atan2(dy, dx))

    # keep upright
    if angle < -90 or angle > 90:
        angle += 180

    # --- perpendicular offset ----------------------------------------------
    # unit normal vector (−dy, dx)
    norm = np.array([-dy, dx])
    if np.linalg.norm(norm) > 0:
        norm = norm / np.linalg.norm(norm)
    else:
        norm = np.array([0.0, 0.0])

    # push label outward (away from centre (0,0))
    centre_vec = np.array([mx, my])
    if np.dot(norm, centre_vec) < 0:
        norm = -norm

    # figure width in data coords ≈ 1 (because we use radii 0-1), adjust if needed
    ox, oy = norm * offset_frac

    # --- draw based on number of lines -------------------------------------
    if len(wrapped_lines) == 1:
        # Single line: place on one side
        mx_off, my_off = mx + ox, my + oy
        ax.text(mx_off, my_off, wrapped_lines[0],
                ha="center", va="center",
                rotation=angle, rotation_mode="anchor",
                fontsize=fontsize, color=color,
                transform_rotates_text=True, zorder=3)
    elif len(wrapped_lines) == 2:
        # Two lines: place one on each side
        mx_off1, my_off1 = mx + ox, my + oy          # First line on one side
        mx_off2, my_off2 = mx - ox, my - oy          # Second line on other side
        
        ax.text(mx_off1, my_off1, wrapped_lines[0],
                ha="center", va="center",
                rotation=angle, rotation_mode="anchor",
                fontsize=fontsize, color=color,
                transform_rotates_text=True, zorder=3)
        ax.text(mx_off2, my_off2, wrapped_lines[1],
                ha="center", va="center",
                rotation=angle, rotation_mode="anchor",
                fontsize=fontsize, color=color,
                transform_rotates_text=True, zorder=3)
    else:
        # More than two lines: fall back to single-side placement with line breaks
        wrapped = "\n".join(wrapped_lines)
        mx_off, my_off = mx + ox, my + oy
        ax.text(mx_off, my_off, wrapped,
                ha="center", va="center",
                rotation=angle, rotation_mode="anchor",
                fontsize=fontsize, color=color,
                transform_rotates_text=True, zorder=3)

def on_click(event):
    """Handle mouse clicks to open URLs for value labels."""
    if event.inaxes is None or event.button != 1:  # Only left clicks on axes
        return
    
    # Check if click is near any value position
    click_x, click_y = event.xdata, event.ydata
    min_distance = float('inf')
    clicked_value = None
    
    # Find the closest value label
    for (reg, val), (x, y) in value_positions.items():
        distance = np.sqrt((click_x - x)**2 + (click_y - y)**2)
        if distance < 0.08 and distance < min_distance:  # Within reasonable click range
            min_distance = distance
            clicked_value = val
    
    # Open URL if value has a reference
    if clicked_value and clicked_value in references:
        print(f"Opening: {references[clicked_value]['title']}")
        webbrowser.open(references[clicked_value]['url'])

def draw_brain_values():
    global value_positions  # Store positions for click detection
    value_positions = {}
    
    # --- figure & axes -------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect("equal")
    ax.axis("off")

    # --- brain image ---------------------------------------------------------
    try:
        brain = mpimg.imread(brain_img_path)
        h, w = brain.shape[:2]
        # scale brain to brain_img_scale of figure width
        extent_w = brain_img_scale
        extent_h = brain_img_scale * h / w
        ax.imshow(brain,
                  extent=[-extent_w/2, +extent_w/2,
                          -extent_h/2, +extent_h/2],
                  zorder=0)
    except FileNotFoundError:
        print(f"Brain image '{brain_img_path}' not found; skipping.")

    # --- positions -----------------------------------------------------------
    region_pos = dict(zip(regions, compute_positions(len(regions), region_radius)))
    # flatten region-value pairs for outer circle
    value_items = [(reg, val) for reg in regions for val in values[reg]]
    value_angles = []
    k_per_region = {r: len(values[r]) for r in regions}
    total_vals = len(value_items)

    # To avoid overlapping, distribute each region’s values in ±15°
    for reg_i, reg in enumerate(regions):
        base_angle = angle_offset + reg_i * 360 / len(regions)
        span = 30 if k_per_region[reg] > 1 else 0
        for j in range(k_per_region[reg]):
            if k_per_region[reg] == 1:
                value_angles.append(base_angle)
            else:
                offset = -span/2 + j * (span / (k_per_region[reg] - 1))
                value_angles.append(base_angle + offset)

    value_pos = {item: polar_to_cart(value_radius, ang)
                 for item, ang in zip(value_items, value_angles)}

    # redistribute space

    value_items = [(reg, val) for reg in regions for val in values[reg]]
    total_vals = len(value_items)
    angle_step  = 360 / total_vals

    value_pos = {}
    for idx, (reg, val) in enumerate(value_items):
        ang = angle_offset + idx * angle_step
        value_pos[(reg, val)] = polar_to_cart(value_radius, ang - 30)


    # --- plot regions --------------------------------------------------------
    for reg, (x, y) in region_pos.items():
        ax.text(x, y, reg, ha="center", va="center",
                fontsize=region_fontsz, bbox=dict(boxstyle="round,pad=0.2",
                                                  fc="#eef", ec="#336"))

    # --- plot values & connectors -------------------------------------------
    for (reg, val), (xv, yv) in value_pos.items():
        xr, yr = region_pos[reg]

        # connector line
        ax.plot([xv, xr], [yv, yr], lw=0.8, color="gray", zorder=1)

        # label on connector (two lines)
        label_on_line(ax, xr, yr, xv, yv,
                  summaries[(reg, val)],
                  fontsize=summary_fontsz)
        #if (reg, val) in summaries:
        #    midx, midy = (xv + xr) / 2, (yv + yr) / 2
        #    t1, t2 = split_summary(summaries[(reg, val)])
        #    ax.text(midx, midy + 0.01, t1, ha="center", va="bottom",
        #            fontsize=summary_fontsz, color="#555")
        #    ax.text(midx, midy - 0.01, t2, ha="center", va="top",
        #            fontsize=summary_fontsz, color="#555")

        # Store position for click detection
        value_positions[(reg, val)] = (xv, yv)
        
        # value label with clickable link (works in SVG output)
        text_kwargs = {
            "ha": "center", "va": "center",
            "fontsize": value_fontsz,
            "bbox": dict(boxstyle="round,pad=0.25", fc="#fee", ec="#933", lw=1.2),
            "zorder": 3
        }
        
        # Add hyperlink if reference exists (only works in SVG)
        if val in references:
            text_kwargs["url"] = references[val]["url"]
            # Make clickable values more visually distinct
            text_kwargs["bbox"]["fc"] = "#ffe"  # Slightly different background
            text_kwargs["bbox"]["ec"] = "#b55"  # Slightly different border
            
        ax.text(xv, yv, val, **text_kwargs)

    # Connect click handler for interactive display
    fig.canvas.mpl_connect('button_press_event', on_click)
    
    # Add title with instruction
    fig.suptitle("Brain-Value Mapping\nClick on value labels to open research papers", 
                 fontsize=14, y=0.95)
    
    plt.tight_layout()
    plt.savefig("brain_values.svg", format='svg', bbox_inches='tight', dpi=300)
    plt.savefig("brain_values.png", format='png', bbox_inches='tight', dpi=300)
    print("Saved brain_values.svg (with clickable links) and brain_values.png")
    print("Click on value labels in the interactive display to open research papers!")
    plt.show()


if __name__ == "__main__":
    draw_brain_values()

