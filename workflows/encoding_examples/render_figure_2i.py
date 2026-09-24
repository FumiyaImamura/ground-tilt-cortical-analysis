"""Render the observed and predicted activity traces in Figure 2i."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source_data" / "encoding_examples" / "figure_02i_model_output.npz"
OUT = ROOT / "outputs" / "encoding_examples"
PANELS = OUT / "panels"
ORANGE = "#e66101"


def normalized_pair(data: np.lib.npyio.NpzFile, index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    start = int(data["display_start_index"][index])
    count = int(data["display_sample_count"])
    observed = data["observed"][index, start : start + count]
    reconstructed = data["reconstructed"][index, start : start + count]
    observed_time = np.arange(count, dtype=float) / 10.0
    center = np.mean(np.r_[observed, reconstructed])
    scale = np.ptp(np.r_[observed, reconstructed])
    return observed_time, (observed - center) / scale, (reconstructed - center) / scale


def render(data: np.lib.npyio.NpzFile, stem: str, score_label: str) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(5.25, 3.8), sharex=True, gridspec_kw={"hspace": 0.03})
    for index, axis in enumerate(axes):
        time, observed, reconstructed = normalized_pair(data, index)
        axis.plot(time, observed, color="black", lw=1.15)
        axis.plot(time, reconstructed, color=ORANGE, lw=1.15)
        axis.text(
            -5.5,
            0.40,
            f"{score_label} = {data['displayed_score'][index]:.2f}",
            ha="left",
            va="top",
            fontsize=10,
            style="italic" if score_label == "r" else "normal",
        )
        axis.set_ylim(-0.60, 0.58)
        axis.axis("off")
    axes[0].plot([], [], color="black", lw=1.15, label="Data")
    axes[0].plot([], [], color=ORANGE, lw=1.15, label="Predicted")
    axes[0].legend(frameon=False, loc="upper right", fontsize=10, handlelength=1.0, handletextpad=0.4)
    axes[0].plot([36, 46], [0.48, 0.48], color="black", lw=1.2)
    axes[0].text(41, 0.51, "10 s", ha="center", va="bottom", fontsize=9)
    axes[-1].set_xlim(-6, 70)
    fig.text(0.01, 0.98, "i", ha="left", va="top", fontsize=14, fontweight="bold")
    fig.savefig(PANELS / f"{stem}.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(PANELS / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    PANELS.mkdir(parents=True, exist_ok=True)
    data = np.load(SOURCE)
    render(data, "figure_02i", "r")
    print(f"Wrote Figure 2i to {PANELS}")


if __name__ == "__main__":
    main()
