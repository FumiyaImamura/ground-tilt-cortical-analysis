"""Render Figure 3h and Extended Data Figure 10k."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from fig3h_common import NBIN, PLOT_ORDER


REPOSITORY = Path(__file__).resolve().parents[2]
SOURCE = REPOSITORY / "source_data" / "frontoparietal_coordination"
OUTPUT = REPOSITORY / "outputs" / "frontoparietal_coordination" / "panels"
OUTPUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.65,
    "xtick.major.width": 0.65,
    "ytick.major.width": 0.65,
    "xtick.major.size": 2.6,
    "ytick.major.size": 2.6,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

CHANCE = 1.0 / NBIN
X9 = np.arange(9)
XTICK_POS = [0, 2, 4, 6, 8]
XTICK_LABELS = ["-180", "-90", "0", "90", "180"]
XLABEL = "Offset (Parietal) - Offset (Frontal)"
BAR_GREY = "0.65"


def read_csv(name: str) -> list[dict[str, str]]:
    with (SOURCE / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def to9(values: np.ndarray) -> np.ndarray:
    ordered = values[PLOT_ORDER]
    return np.concatenate([ordered, ordered[:1]])


def duplicate_endpoint(values: np.ndarray) -> np.ndarray:
    return np.concatenate([values, values[:1]])


def style_histogram_axis(ax: plt.Axes, ylabel: str = "Probability") -> None:
    ax.set_xticks(XTICK_POS)
    ax.set_xticklabels(XTICK_LABELS)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(-0.7, 8.7)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.20,
        1.05,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
    )


def significance_text(p_value: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUTPUT / f"{stem}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(OUTPUT / f"{stem}.png", bbox_inches="tight", dpi=450)
    plt.close(fig)
    print(f"Wrote {stem}")


def render_figure_3h(animal_probabilities: np.ndarray) -> None:
    rows = read_csv("animal_level_stats.csv")
    means = np.array([float(row["mean"]) for row in rows])
    sems = np.array([float(row["sem"]) for row in rows])
    corrected = np.array([float(row["wilcoxon_p_bonf"]) for row in rows])
    means9 = duplicate_endpoint(means)
    sems9 = duplicate_endpoint(sems)

    fig, ax = plt.subplots(figsize=(2.45, 1.9))
    for values in animal_probabilities:
        ax.plot(X9, to9(values), color="0.25", lw=0.35, alpha=0.32, zorder=1)
    ax.bar(
        X9,
        means9,
        width=0.86,
        color=BAR_GREY,
        edgecolor="black",
        linewidth=0.65,
        zorder=2,
    )
    ax.errorbar(
        X9,
        means9,
        yerr=sems9,
        fmt="none",
        ecolor="black",
        elinewidth=0.75,
        capsize=1.7,
        capthick=0.75,
        zorder=3,
    )
    ax.axhline(CHANCE, color="black", ls=(0, (3, 3)), lw=0.55, zorder=0)
    zero_position = 4
    ax.text(
        zero_position,
        means9[zero_position] + sems9[zero_position] + 0.012,
        significance_text(corrected[zero_position]),
        ha="center",
        va="bottom",
        fontsize=8,
    )
    ax.set_ylim(0, 0.325)
    style_histogram_axis(ax)
    panel_label(ax, "h")
    save(fig, "figure_03h")


def render_extended_data_10k() -> None:
    rows = read_csv("independence_null_one_sided.csv")
    excess = np.array([float(row["excess_mean"]) for row in rows])
    excess_sem = np.array([float(row["excess_sem"]) for row in rows])
    corrected = np.array([
        float(row["wilcoxon_one_sided_greater_p_bonferroni_8"])
        for row in rows
    ])
    excess9 = duplicate_endpoint(excess)
    sem9 = duplicate_endpoint(excess_sem)

    fig, ax = plt.subplots(figsize=(2.55, 1.9))
    ax.bar(
        X9,
        excess9,
        width=0.86,
        color=BAR_GREY,
        edgecolor="black",
        linewidth=0.65,
    )
    ax.errorbar(
        X9,
        excess9,
        yerr=sem9,
        fmt="none",
        ecolor="black",
        elinewidth=0.75,
        capsize=1.7,
        capthick=0.75,
    )
    ax.axhline(0, color="black", lw=0.65)
    ax.set_ylim(-0.032, 0.092)
    for index, p_value in enumerate(corrected):
        marker = significance_text(p_value)
        if marker:
            ax.text(
                index,
                excess9[index] + sem9[index] + 0.006,
                marker,
                ha="center",
                va="bottom",
                fontsize=8,
            )
    style_histogram_axis(ax, ylabel="Δ Probability")
    panel_label(ax, "k")
    save(fig, "extended_data_figure_10k")


def main() -> None:
    with np.load(SOURCE / "figure_03h_values.npz", allow_pickle=False) as archive:
        animal_probabilities = np.asarray(archive["animal_probabilities"], dtype=float)
    render_figure_3h(animal_probabilities)
    render_extended_data_10k()


if __name__ == "__main__":
    main()
