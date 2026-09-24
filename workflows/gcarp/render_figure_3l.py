"""Render Figure 3l from the GCARP component-by-session matrices."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPO_ROOT / "source_data" / "gcarp"
OUT_PANELS = REPO_ROOT / "outputs" / "gcarp" / "panels"
OUT_RESULTS = REPO_ROOT / "outputs" / "gcarp" / "results"

SOURCES = {
    "Movement": ("movement.mat", "delta_R2_axis_plane"),
    "Tilt": ("ground_tilt.mat", "delta_degScore_axis_plane"),
    "Posture": ("eigen_posture_1.mat", "delta_R2_axis_plane"),
}

# Significance annotations shown in Figure 3l.
SIGNIFICANCE_MARKERS = {
    "Movement": {1: "**"},
    "Tilt": {1: "*", 2: "**", 3: "*"},
    "Posture": {1: "*", 2: "*", 3: "*"},
}

YLIMS = {
    "Movement": (-0.015, 0.43),
    "Tilt": (-0.012, 0.065),
    "Posture": (-0.004, 0.065),
}

YTICKS = {
    "Movement": [0.0, 0.1, 0.2, 0.3, 0.4],
    "Tilt": [0.0, 0.02, 0.04, 0.06],
    "Posture": [0.0, 0.02, 0.04, 0.06],
}

YMINOR_TICKS = {
    "Movement": [0.05, 0.15, 0.25, 0.35],
    "Tilt": [0.01, 0.03, 0.05],
    "Posture": [0.01, 0.03, 0.05],
}


def load_values() -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for title, (filename, variable) in SOURCES.items():
        values = np.asarray(loadmat(SOURCE_ROOT / filename)[variable], dtype=float)
        if values.shape != (10, 20):
            raise ValueError(f"{filename}: expected shape (10, 20), found {values.shape}")
        result[title] = values
    return result


def plot_panel(
    ax: plt.Axes,
    title: str,
    values: np.ndarray,
    stars: dict[str, dict[int, str]],
) -> None:
    x = np.arange(1, 11)
    mean = values.mean(axis=1)
    sem = values.std(axis=1, ddof=1) / np.sqrt(values.shape[1])
    ax.errorbar(
        x,
        mean,
        yerr=sem,
        color="black",
        linewidth=1.0,
        elinewidth=0.8,
        capsize=2.0,
        capthick=0.8,
    )
    lower, upper = YLIMS[title]
    ax.set_xlim(0.55, 10.45)
    ax.set_ylim(lower, upper)
    ax.set_xticks([1, 4, 7, 10])
    ax.set_xticks([2, 3, 5, 6, 8, 9], minor=True)
    ax.set_yticks(YTICKS[title])
    ax.set_yticks(YMINOR_TICKS[title], minor=True)
    ax.set_xlabel("GCCA component")
    ax.set_ylabel("GCARP contribution")
    ax.set_title(title, pad=5)
    for component, marker in stars[title].items():
        index = component - 1
        offset = (upper - lower) * 0.025
        ax.text(component, mean[index] + sem[index] + offset, marker, ha="center", va="bottom", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=2.5, width=0.7, labelsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)


def render_panel(
    title: str,
    values: np.ndarray,
    stem: str,
    stars: dict[str, dict[int, str]],
) -> None:
    fig, ax = plt.subplots(figsize=(2.3, 1.9))
    plot_panel(ax, title, values, stars)
    fig.tight_layout()
    OUT_PANELS.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PANELS / f"{stem}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(OUT_PANELS / f"{stem}.png", bbox_inches="tight", dpi=300, transparent=True)
    plt.close(fig)


def render_combined(
    values: dict[str, np.ndarray],
    stars: dict[str, dict[int, str]],
    stem: str,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(6.9, 1.9))
    for ax, title in zip(axes, ("Movement", "Tilt", "Posture")):
        plot_panel(ax, title, values[title], stars)
    fig.tight_layout(w_pad=1.0)
    OUT_PANELS.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PANELS / f"{stem}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(OUT_PANELS / f"{stem}.png", bbox_inches="tight", dpi=300, transparent=True)
    plt.close(fig)


def main() -> int:
    values = load_values()
    stems = {
        "Movement": "figure_03l_movement",
        "Tilt": "figure_03l_tilt",
        "Posture": "figure_03l_posture",
    }
    for title, matrix in values.items():
        render_panel(title, matrix, stems[title], SIGNIFICANCE_MARKERS)
    render_combined(values, SIGNIFICANCE_MARKERS, "figure_03l")

    OUT_RESULTS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT_RESULTS / "figure_03l_values.npz",
        movement=values["Movement"],
        ground_tilt=values["Tilt"],
        eigenposture_1=values["Posture"],
    )
    print(f"Wrote Figure 3l panels to {OUT_PANELS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
