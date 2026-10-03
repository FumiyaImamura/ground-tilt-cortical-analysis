"""Render Figure 3l from the GCARP component-by-session matrices."""

from __future__ import annotations

from pathlib import Path
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.stats import wilcoxon


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPO_ROOT / "source_data" / "gcarp"
OUT_PANELS = REPO_ROOT / "outputs" / "gcarp" / "panels"
OUT_RESULTS = REPO_ROOT / "outputs" / "gcarp" / "results"

SOURCES = {
    "Movement": ("movement.mat", "delta_R2_axis_plane"),
    "Tilt": ("ground_tilt.mat", "delta_degScore_axis_plane"),
    "Posture": ("eigen_posture_1.mat", "delta_R2_axis_plane"),
}

BONFERRONI_FAMILY = 10

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


def significance_symbol(p_value: float) -> str:
    if p_value < 0.001:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def matlab_default_two_sided_signrank(values: np.ndarray) -> float:
    """Match MATLAB's default two-sided signrank method at each sample size."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values != 0)]
    method = "exact" if values.size <= 15 else "approx"
    return float(
        wilcoxon(
            values,
            alternative="two-sided",
            method=method,
            correction=False,
        ).pvalue
    )


def component_tests(
    values: dict[str, np.ndarray],
) -> tuple[dict[str, dict[int, str]], list[dict[str, object]]]:
    markers: dict[str, dict[int, str]] = {}
    rows: list[dict[str, object]] = []
    for title, matrix in values.items():
        markers[title] = {}
        for component, component_values in enumerate(matrix, start=1):
            raw_p = matlab_default_two_sided_signrank(component_values)
            corrected_p = min(1.0, raw_p * BONFERRONI_FAMILY)
            marker = significance_symbol(corrected_p)
            if marker:
                markers[title][component] = marker
            rows.append(
                {
                    "target": title,
                    "component": component,
                    "n_sessions": int(np.isfinite(component_values).sum()),
                    "p_raw_two_sided_signed_rank": raw_p,
                    "p_bonferroni_10": corrected_p,
                    "significance_symbol": marker or "n.s.",
                }
            )
    return markers, rows


def write_tests(rows: list[dict[str, object]]) -> None:
    OUT_RESULTS.mkdir(parents=True, exist_ok=True)
    with (OUT_RESULTS / "figure_03l_component_tests.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


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
    markers, test_rows = component_tests(values)
    stems = {
        "Movement": "figure_03l_movement",
        "Tilt": "figure_03l_tilt",
        "Posture": "figure_03l_posture",
    }
    for title, matrix in values.items():
        render_panel(title, matrix, stems[title], markers)
    render_combined(values, markers, "figure_03l")

    OUT_RESULTS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT_RESULTS / "figure_03l_values.npz",
        movement=values["Movement"],
        ground_tilt=values["Tilt"],
        eigenposture_1=values["Posture"],
    )
    write_tests(test_rows)
    print(f"Wrote Figure 3l panels to {OUT_PANELS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
