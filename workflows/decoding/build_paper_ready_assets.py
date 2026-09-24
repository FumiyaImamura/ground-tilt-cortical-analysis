"""Build the neuron-count and shuffled-target decoding panels.

The neuron-count panels use the paper's session-level aggregation. The ten
neuron-sampling iterations are aggregated within hemisphere/session. Each area
is compared with the within-session mean of the other seven areas using a
one-sided Wilcoxon signed-rank test, Bonferroni corrected across eight areas
within each motion and neuron count. MATLAB's default signrank p-value behavior
is reproduced (exact for n <= 15; normal approximation with continuity
correction otherwise). Error bars follow the paper's plotting convention:
sample SD across available sessions divided by sqrt(20).

The empirical-null panels use the same mouse aggregation and a paired one-sided
Wilcoxon signed-rank test of real versus shuffled-target decoding, with the same
eight-area Bonferroni family. Summary tables are written to `results/`, and
panels are written as separate SVG and PNG files to `panels/`.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
import warnings

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.stats import wilcoxon


REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / "source_data" / "decoding"
RESULTS = REPO / "outputs" / "decoding" / "results"
FIGURES = REPO / "outputs" / "decoding" / "panels"
RESULTS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

PERREP_FILES = {
    100: SOURCE / "perrep_scores_nc100.mat",
    200: SOURCE / "perrep_scores.mat",
    300: SOURCE / "perrep_scores_nc300.mat",
}
SHUFFLE_FILE = SOURCE / "perrep_scores_shuffle.mat"
COMPARABLE_FILE = SOURCE / "comparable_session_scores.mat"

MOTIONS = [
    ("p8", "8DT", True, 3 / 8),
    ("rot", "T-Rot", True, 3 / 8),
    ("tilt", "Roll", False, 0.0),
    ("yaw", "Yaw", False, 0.0),
    ("zvert", "Vert", False, 0.0),
]
REGIONS = ["M2", "M1", "Whisker", "UL", "LL", "Trunk", "PPC", "RSC"]
LABELS = ["M2", "M1", "Whisker", "FL", "HL", "Trunk", "PPC", "RSC"]
DISPLAY = dict(zip(REGIONS, LABELS))
COUNTS = [100, 200, 300]
COUNT_COLORS = {100: "#A8DADC", 200: "#457B9D", 300: "#1D3557"}
REAL_COLOR = "#5E5E5E"
SHUFFLE_COLOR = "#C7C7C7"
REFERENCE_COLOR = "#B2182B"


def load_perrep(path: Path):
    return sio.loadmat(path, struct_as_record=False, squeeze_me=True)["perrep"]


def session_values(perrep, motion: str, region: str, reference: float = 0.0):
    matrix = np.asarray(getattr(getattr(perrep, motion), region), dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        values = np.nanmean(matrix, axis=1)
    return values - reference


def mouse_values(session_vector):
    vector = np.asarray(session_vector, dtype=float)
    if vector.size != 20:
        raise ValueError(f"Expected 20 hemisphere/session rows, found {vector.size}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(vector.reshape(10, 2), axis=1)


def mean_sem(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    mean = float(np.mean(values))
    sem = float(np.std(values, ddof=1) / math.sqrt(values.size)) if values.size > 1 else math.nan
    return mean, sem, int(values.size)


def paper_session_sem(values, total_sessions=20):
    """Paper convention: sample SD across sessions divided by sqrt(20)."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    mean = float(np.mean(values))
    sem = float(np.std(values, ddof=1) / math.sqrt(total_sessions)) if values.size > 1 else math.nan
    return mean, sem, int(values.size)


def matlab_default_signrank_greater(differences):
    """Match MATLAB signrank(...,'tail','right') default method selection."""
    differences = np.asarray(differences, dtype=float)
    differences = differences[np.isfinite(differences) & (differences != 0)]
    if differences.size < 3:
        return math.nan
    if differences.size <= 15:
        return float(wilcoxon(differences, alternative="greater", method="exact").pvalue)
    return float(wilcoxon(
        differences, alternative="greater", method="approx", correction=True,
    ).pvalue)


def star(p):
    if not np.isfinite(p):
        return "n.s."
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def write_csv(name, header, rows):
    path = RESULTS / name
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def style_axis(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=3)


def add_s1_bracket(ax):
    transform = ax.get_xaxis_transform()
    ax.plot([2 - 0.38, 5 + 0.38], [-0.28, -0.28], color="black", lw=0.9,
            transform=transform, clip_on=False)
    ax.text(3.5, -0.37, "S1", ha="center", va="top", fontsize=9,
            transform=transform, clip_on=False)


perrep = {count: load_perrep(path) for count, path in PERREP_FILES.items()}
shuffle = load_perrep(SHUFFLE_FILE)
comparable = sio.loadmat(
    COMPARABLE_FILE, struct_as_record=False, squeeze_me=True,
)["comparable"]


# ---------- Neuron-count statistics ----------
neuron_rows = []
neuron_summary_rows = []
neuron_cache = {}
for motion, title, classification, reference in MOTIONS:
    for count in COUNTS:
        session_by_region = {
            region: np.asarray(
                getattr(getattr(getattr(comparable, f"nc{count}"), motion), region),
                dtype=float,
            ).ravel() - reference
            for region in REGIONS
        }
        mouse_by_region = {
            region: mouse_values(values)
            for region, values in session_by_region.items()
        }
        means = {region: paper_session_sem(values)[0]
                 for region, values in session_by_region.items()}
        order = sorted(REGIONS, key=lambda region: -means[region])
        ranks = {region: index + 1 for index, region in enumerate(order)}
        significant = []
        for region in REGIONS:
            area = session_by_region[region]
            contrasts = []
            for session_index in range(20):
                if not np.isfinite(area[session_index]):
                    continue
                other = [session_by_region[o][session_index] for o in REGIONS if o != region
                         and np.isfinite(session_by_region[o][session_index])]
                if other:
                    contrasts.append(area[session_index] - np.mean(other))
            contrasts = np.asarray(contrasts, dtype=float)
            raw_p = matlab_default_signrank_greater(contrasts)
            corrected_p = min(1.0, raw_p * 8) if np.isfinite(raw_p) else math.nan
            symbol = star(corrected_p)
            if symbol != "n.s.":
                significant.append(DISPLAY[region])
            session_mean, session_sem, n_sessions = paper_session_sem(
                session_by_region[region]
            )
            mouse_mean, mouse_sem, n_mice = mean_sem(mouse_by_region[region])
            neuron_rows.append([
                motion, title, count, DISPLAY[region], n_sessions, n_mice,
                f"{session_mean:.8f}", f"{session_sem:.8f}",
                f"{mouse_mean:.8f}", f"{mouse_sem:.8f}", ranks[region],
                contrasts.size, f"{raw_p:.10f}", f"{corrected_p:.10f}", symbol,
            ])
            neuron_cache[(motion, count, region)] = (
                session_mean, session_sem, symbol, corrected_p,
            )
        neuron_summary_rows.append([
            motion, title, count,
            DISPLAY[order[0]], DISPLAY[order[1]], DISPLAY[order[2]],
            ranks["LL"], ranks["Trunk"],
            str(order[0] in ("LL", "Trunk")),
            str(ranks["LL"] <= 2 and ranks["Trunk"] <= 2),
            ";".join(significant) if significant else "none",
        ])

write_csv(
    "neuron_count_stats.csv",
    ["motion_key", "motion", "neurons", "region", "n_sessions", "n_mice",
     "session_mean_delta", "session_sem", "mouse_mean_delta", "mouse_sem",
     "rank", "contrast_n_sessions", "p_raw_one_sided_signed_rank",
     "p_bonferroni_8", "star_bonferroni_8"],
    neuron_rows,
)
write_csv(
    "neuron_count_summary.csv",
    ["motion_key", "motion", "neurons", "top1", "top2", "top3", "hl_rank",
     "trunk_rank", "hl_or_trunk_top1", "both_hl_trunk_top2",
     "significant_regions_session_level"],
    neuron_summary_rows,
)

# ---------- Empirical-shuffle statistics ----------
shuffle_rows = []
shuffle_summary_rows = []
shuffle_cache = {}
for motion, title, classification, analytic_reference in MOTIONS:
    pooled_shuffle = []
    for region in REGIONS:
        real_mouse = mouse_values(session_values(perrep[200], motion, region, 0.0))
        shuffle_mouse = mouse_values(session_values(shuffle, motion, region, 0.0))
        valid = np.isfinite(real_mouse) & np.isfinite(shuffle_mouse)
        real_valid = real_mouse[valid]
        shuffle_valid = shuffle_mouse[valid]
        pooled_shuffle.extend(shuffle_valid.tolist())
        real_mean, real_sem, n = mean_sem(real_valid)
        shuffle_mean, shuffle_sem, _ = mean_sem(shuffle_valid)
        raw_p = float(wilcoxon(real_valid, shuffle_valid, alternative="greater").pvalue)
        p8 = min(1.0, raw_p * 8)
        p40 = min(1.0, raw_p * 40)
        symbol = star(p8)
        positive_mice = int(np.sum(real_valid > shuffle_valid))
        shuffle_rows.append([
            motion, title, DISPLAY[region], n, positive_mice,
            f"{real_mean:.8f}", f"{real_sem:.8f}",
            f"{shuffle_mean:.8f}", f"{shuffle_sem:.8f}",
            f"{analytic_reference:.8f}", f"{raw_p:.10f}", f"{p8:.10f}",
            f"{p40:.10f}", symbol,
        ])
        shuffle_cache[(motion, region)] = (
            real_mean, real_sem, shuffle_mean, shuffle_sem, symbol, p8,
        )
    shuffle_summary_rows.append([
        motion, title, len(pooled_shuffle), f"{np.mean(pooled_shuffle):.8f}",
        f"{analytic_reference:.8f}",
    ])

write_csv(
    "empirical_shuffle_stats.csv",
    ["motion_key", "motion", "region", "n_mice", "positive_mice",
     "real_mean", "real_mouse_sem", "shuffle_mean", "shuffle_mouse_sem",
     "analytic_reference", "p_raw_paired_one_sided", "p_bonferroni_8",
     "p_bonferroni_40", "star_bonferroni_8"],
    shuffle_rows,
)
write_csv(
    "empirical_shuffle_summary.csv",
    ["motion_key", "motion", "pooled_region_mouse_values", "pooled_shuffle_mean",
     "analytic_reference"],
    shuffle_summary_rows,
)


# ---------- One-panel neuron-count figures ----------
for motion, title, classification, reference in MOTIONS:
    fig, ax = plt.subplots(figsize=(4.45, 3.75))
    x = np.arange(len(REGIONS))
    width = 0.24
    all_tops = []
    for count_index, count in enumerate(COUNTS):
        means = np.asarray([neuron_cache[(motion, count, region)][0] for region in REGIONS])
        sems = np.asarray([neuron_cache[(motion, count, region)][1] for region in REGIONS])
        symbols = [neuron_cache[(motion, count, region)][2] for region in REGIONS]
        x_count = x + (count_index - 1) * width
        bars = ax.bar(
            x_count, means, width, yerr=sems, color=COUNT_COLORS[count],
            edgecolor="0.25", linewidth=0.55, label=f"{count} neurons",
            error_kw=dict(ecolor="0.25", lw=0.65, capsize=1.5, capthick=0.65),
        )
        all_tops.extend((means + sems).tolist())
        for bar, mean, sem, symbol in zip(bars, means, sems, symbols):
            if symbol != "n.s.":
                ax.annotate(
                    symbol,
                    xy=(bar.get_x() + bar.get_width() / 2, mean + sem),
                    xytext=(0, 2), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7,
                    fontweight="bold", color="black",
                )
    ax.set_title(title, fontsize=12, pad=8)
    ax.set_ylabel(r"$\Delta$ Decoding accuracy", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylim(0, max(all_tops) * 1.25)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    add_s1_bracket(ax)
    style_axis(ax)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.90, bottom=0.31)
    stem = f"neuron_count_{motion}"
    fig.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


# ---------- One-panel empirical-null figures ----------
for motion, title, classification, analytic_reference in MOTIONS:
    fig, ax = plt.subplots(figsize=(4.45, 3.75))
    x = np.arange(len(REGIONS))
    width = 0.36
    real_means = np.asarray([shuffle_cache[(motion, region)][0] for region in REGIONS])
    real_sems = np.asarray([shuffle_cache[(motion, region)][1] for region in REGIONS])
    shuffle_means = np.asarray([shuffle_cache[(motion, region)][2] for region in REGIONS])
    shuffle_sems = np.asarray([shuffle_cache[(motion, region)][3] for region in REGIONS])
    symbols = [shuffle_cache[(motion, region)][4] for region in REGIONS]
    ax.bar(
        x - width / 2, real_means, width, yerr=real_sems, color=REAL_COLOR,
        edgecolor="0.2", linewidth=0.55, label="Real",
        error_kw=dict(ecolor="0.2", lw=0.65, capsize=1.5, capthick=0.65),
    )
    ax.bar(
        x + width / 2, shuffle_means, width, yerr=shuffle_sems, color=SHUFFLE_COLOR,
        edgecolor="0.2", linewidth=0.55, label="Shuffle",
        error_kw=dict(ecolor="0.2", lw=0.65, capsize=1.5, capthick=0.65),
    )
    ax.axhline(analytic_reference, color=REFERENCE_COLOR, ls="--", lw=0.9,
               label="Analytic reference")
    span = max(real_means + real_sems) - min(shuffle_means - shuffle_sems)
    for xi, mean, sem, symbol in zip(x, real_means, real_sems, symbols):
        if symbol != "n.s.":
            ax.text(xi, mean + sem + 0.015 * span, symbol, ha="center",
                    va="bottom", fontsize=8, fontweight="bold")
    ax.set_title(title, fontsize=12, pad=8)
    ax.set_ylabel("P(within scoring tolerance)" if classification else r"$R^2$", fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, rotation=45, ha="right", rotation_mode="anchor")
    if classification:
        lower = max(0.0, min(shuffle_means - shuffle_sems) - 0.04)
        upper = max(real_means + real_sems) + 0.10
    else:
        lower = min(shuffle_means - shuffle_sems) - 0.025
        upper = max(real_means + real_sems) + 0.045
    ax.set_ylim(lower, upper)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    style_axis(ax)
    fig.subplots_adjust(left=0.17, right=0.98, top=0.90, bottom=0.25)
    stem = f"empirical_null_{motion}"
    fig.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

print(f"Wrote results to {RESULTS}")
print(f"Wrote manuscript panels to {FIGURES}")
