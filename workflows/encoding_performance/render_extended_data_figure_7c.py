"""Render Extended Data Figure 7c from the bundled histogram arrays."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from scipy.io import loadmat


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = REPO_ROOT / "source_data" / "encoding_performance"
OUTPUT_DIR = REPO_ROOT / "outputs" / "encoding_performance"
SOURCE_FILE = SOURCE_DIR / "ExtFig_Rev3Q2_panel_a_data.mat"

X_LABEL = "Fold-averaged cross-validated Pearson r"
Y_LABEL = "Fraction of neurons"

COLORS = np.asarray(
    [
        [0.84, 0.19, 0.15],
        [0.20, 0.45, 0.72],
        [0.25, 0.62, 0.28],
        [0.55, 0.27, 0.60],
        [0.95, 0.55, 0.05],
        [0.60, 0.32, 0.16],
        [0.90, 0.45, 0.70],
        [0.40, 0.40, 0.40],
        [0.00, 0.65, 0.68],
        [0.72, 0.63, 0.00],
    ]
)


def compact_number(value: float, _position: float | None = None) -> str:
    if abs(value) < 1e-12:
        return "0"
    return f"{value:g}"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    data = loadmat(SOURCE_FILE, squeeze_me=True, struct_as_record=False)
    centres = np.asarray(data["centres"], dtype=float).reshape(-1)
    edges = np.asarray(data["edges"], dtype=float).reshape(-1)
    hist_probability = np.asarray(data["histProbability"], dtype=float)

    if hist_probability.shape != (20, 96):
        raise ValueError(f"Expected a 20 x 96 histogram array, got {hist_probability.shape}")
    if centres.size != hist_probability.shape[1] or edges.size != centres.size + 1:
        raise ValueError("Histogram centres/edges do not match the bundled probability array")

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.linewidth": 0.8,
            "svg.fonttype": "none",
        }
    )
    fig, ax = plt.subplots(figsize=(6.3, 4.0), facecolor="white")
    fig.subplots_adjust(left=0.145, right=0.98, bottom=0.19, top=0.97)

    # The bundled rows are the 20 sessions in the paper order, two
    # hemispheres per animal. Colour therefore repeats for each adjacent pair;
    # solid and dashed lines distinguish the paired hemispheres.
    for session_index, probability in enumerate(hist_probability):
        animal_index = session_index // 2
        line_style = "-" if session_index % 2 == 0 else "--"
        ax.plot(
            centres,
            probability,
            color=COLORS[animal_index],
            linestyle=line_style,
            linewidth=1.15,
        )

    ax.set_xlim(-0.12, 0.84)
    ax.set_ylim(0.0, 0.18)
    ax.set_xticks((-0.1, 0.0, 0.2, 0.4, 0.6, 0.8))
    ax.set_yticks((0.0, 0.04, 0.08, 0.12, 0.16))
    ax.xaxis.set_major_formatter(FuncFormatter(compact_number))
    ax.yaxis.set_major_formatter(FuncFormatter(compact_number))
    ax.set_xlabel(X_LABEL, fontsize=10.0)
    ax.set_ylabel(Y_LABEL, fontsize=10.0)
    ax.tick_params(axis="both", labelsize=8.5, width=0.8, length=4, pad=2)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_linewidth(0.8)

    svg_path = OUTPUT_DIR / "extended_data_figure_7c_encoding_performance.svg"
    png_path = OUTPUT_DIR / "extended_data_figure_7c_encoding_performance.png"
    fig.savefig(svg_path, facecolor="white")
    fig.savefig(png_path, dpi=600, facecolor="white")
    plt.close(fig)

    print(f"Wrote {svg_path}")
    print(f"Wrote {png_path}")


if __name__ == "__main__":
    main()
