"""Render Extended Data Figure 9 from the bundled mouse-level source data.

The workflow reads the analysis tables under ``source_data/sex_indicator`` and
does not require the original session-level analysis workspace.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = REPO_ROOT / "source_data" / "sex_indicator"
OUTPUT_DIR = REPO_ROOT / "outputs" / "sex_indicator"


PANEL_SPECS = (
    {
        "outcome": "direction_selective_fraction",
        "title": "Direction selectivity",
        "ylabel": "Direction-selective\nneurons (%)",
        "scale": 100.0,
        "ylim": (0.0, 8.5),
        "yticks": (0.0, 2.5, 5.0, 7.5),
    },
    {
        "outcome": "decoding_p8",
        "title": "8DT",
        "ylabel": "\u0394 Decoding accuracy",
        "scale": 1.0,
        "ylim": (0.0, 0.085),
        "yticks": (0.0, 0.025, 0.05, 0.075),
    },
    {
        "outcome": "decoding_rot",
        "title": "T-Rot",
        "ylabel": "\u0394 Decoding accuracy",
        "scale": 1.0,
        "ylim": (0.0, 0.34),
        "yticks": (0.0, 0.1, 0.2, 0.3),
    },
    {
        "outcome": "decoding_tilt",
        "title": "Roll",
        "ylabel": "\u0394 Decoding accuracy",
        "scale": 1.0,
        "ylim": (0.0, 0.135),
        "yticks": (0.0, 0.04, 0.08, 0.12),
    },
    {
        "outcome": "decoding_yaw",
        "title": "Yaw",
        "ylabel": "\u0394 Decoding accuracy",
        "scale": 1.0,
        "ylim": (0.0, 0.28),
        "yticks": (0.0, 0.08, 0.16, 0.24),
    },
    {
        "outcome": "decoding_zvert",
        "title": "Vert",
        "ylabel": "\u0394 Decoding accuracy",
        "scale": 1.0,
        "ylim": (0.0, 0.28),
        "yticks": (0.0, 0.08, 0.16, 0.24),
    },
)


GROUP_SPECS = (
    {
        "column": "sex",
        "table": "paper_table_sex.csv",
        "order": ("male", "female"),
        "labels": ("Male\n($n$=7)", "Female\n($n$=3)"),
        "panel_letter": "a",
    },
    {
        "column": "indicator",
        "table": "paper_table_indicator.csv",
        "order": ("GCaMP8s", "GCaMP6s"),
        "labels": ("GCaMP8s\n($n$=8)", "GCaMP6s\n($n$=2)"),
        "panel_letter": "b",
    },
)


def compact_number(value: float, _position: float | None = None) -> str:
    """Match the compact numeric tick labels in the manuscript composite."""

    if abs(value) < 1e-12:
        return "0"
    return f"{value:g}"


def draw_panel(
    ax: plt.Axes,
    wide: pd.DataFrame,
    table: pd.DataFrame,
    group_spec: dict,
    panel_spec: dict,
) -> None:
    outcome = panel_spec["outcome"]
    scale = float(panel_spec["scale"])
    groups = group_spec["order"]
    table_row = table.loc[table["outcome"] == outcome].iloc[0]

    for x, group in enumerate(groups):
        values = (
            wide.loc[wide[group_spec["column"]] == group, outcome]
            .dropna()
            .to_numpy(dtype=float)
            * scale
        )
        offsets = (
            np.linspace(-0.065, 0.065, values.size)
            if values.size > 1
            else np.array([0.0])
        )
        ax.scatter(
            x + offsets,
            values,
            s=21,
            marker="o",
            facecolor="#A6A6A6",
            edgecolor="black",
            linewidth=0.55,
            zorder=3,
        )

        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=1))
        ax.errorbar(
            x,
            mean,
            yerr=sd,
            fmt="none",
            ecolor="black",
            elinewidth=0.8,
            capsize=2.7,
            capthick=0.8,
            zorder=4,
        )
        ax.hlines(mean, x - 0.16, x + 0.16, color="black", linewidth=1.35, zorder=5)

        if group == table_row["target_group"]:
            expected_n = int(table_row["target_n"])
            expected_mean = float(table_row["target_mean"])
            expected_sd = float(table_row["target_sd"])
        elif group == table_row["reference_group"]:
            expected_n = int(table_row["reference_n"])
            expected_mean = float(table_row["reference_mean"])
            expected_sd = float(table_row["reference_sd"])
        else:
            raise ValueError(f"{group!r} is absent from the summary table for {outcome}")

        if values.size != expected_n:
            raise ValueError(f"Unexpected sample size for {group}, {outcome}")
        if not np.isclose(mean, expected_mean, rtol=0.0, atol=1e-12):
            raise ValueError(f"Mean does not match the summary table for {group}, {outcome}")
        if not np.isclose(sd, expected_sd, rtol=0.0, atol=1e-12):
            raise ValueError(f"SD does not match the summary table for {group}, {outcome}")

    p_value = float(table_row["exact_two_sided_permutation_p"])
    p_label = f"p={p_value:.3f}"
    ax.plot(
        [0, 0, 1, 1],
        [0.82, 0.85, 0.85, 0.82],
        color="black",
        linewidth=0.7,
        clip_on=False,
        transform=ax.get_xaxis_transform(),
    )
    ax.text(
        0.5,
        0.87,
        p_label,
        ha="center",
        va="bottom",
        fontsize=7.2,
        transform=ax.get_xaxis_transform(),
    )

    ax.set_xlim(-0.42, 1.42)
    ax.set_ylim(*panel_spec["ylim"])
    ax.set_yticks(panel_spec["yticks"])
    ax.yaxis.set_major_formatter(FuncFormatter(compact_number))
    ax.set_xticks((0, 1), group_spec["labels"])
    ax.set_ylabel(panel_spec["ylabel"], fontsize=8.0)
    ax.set_title(panel_spec["title"], fontsize=9.5, pad=3.5)
    ax.tick_params(axis="both", labelsize=7.1, width=0.75, length=3, pad=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_linewidth(0.8)



def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    wide = pd.read_csv(SOURCE_DIR / "mouse_level_outcomes.csv")

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.linewidth": 0.8,
            "svg.fonttype": "none",
            "mathtext.default": "regular",
        }
    )

    fig = plt.figure(figsize=(8.2, 9.0), facecolor="white")
    grid = fig.add_gridspec(
        nrows=4,
        ncols=3,
        left=0.105,
        right=0.975,
        bottom=0.055,
        top=0.965,
        wspace=0.67,
        hspace=0.62,
    )

    for group_index, group_spec in enumerate(GROUP_SPECS):
        table = pd.read_csv(SOURCE_DIR / group_spec["table"])
        for panel_index, panel_spec in enumerate(PANEL_SPECS):
            row = group_index * 2 + panel_index // 3
            col = panel_index % 3
            ax = fig.add_subplot(grid[row, col])
            draw_panel(ax, wide, table, group_spec, panel_spec)

    fig.text(0.010, 0.985, "a", fontsize=11.5, fontweight="bold", ha="left", va="top")
    fig.text(0.010, 0.498, "b", fontsize=11.5, fontweight="bold", ha="left", va="top")

    svg_path = OUTPUT_DIR / "extended_data_figure_9_sex_indicator.svg"
    png_path = OUTPUT_DIR / "extended_data_figure_9_sex_indicator.png"
    fig.savefig(svg_path, facecolor="white")
    fig.savefig(png_path, dpi=300, facecolor="white")
    plt.close(fig)

    svg_text = svg_path.read_text(encoding="utf-8")
    expected_p_labels = [
        f"p={float(value):.3f}"
        for table_name in ("paper_table_sex.csv", "paper_table_indicator.csv")
        for value in pd.read_csv(SOURCE_DIR / table_name)[
            "exact_two_sided_permutation_p"
        ]
    ]
    missing_labels = [label for label in expected_p_labels if label not in svg_text]
    if missing_labels:
        raise RuntimeError(f"Missing p-value labels in SVG: {missing_labels}")

    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
