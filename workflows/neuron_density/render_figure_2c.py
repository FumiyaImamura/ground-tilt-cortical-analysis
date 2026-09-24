"""Render the recorded-neuron locations and density portions of Figure 2c."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.signal import convolve2d


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE = REPO_ROOT / "source_data" / "neuron_density" / "figure_02c_source.mat"
OUT_PANELS = REPO_ROOT / "outputs" / "neuron_density" / "panels"

# Match the physical cortical-image calibration applied by CSP.plotScaleBar.
# MATLAB uses daspect([(911/18.89)/(904/15.02), 1, 1]); Matplotlib's
# corresponding numeric aspect is the same displayed y-unit/x-unit ratio.
CORTICAL_DATA_ASPECT = (911.0 / 18.89) / (904.0 / 15.02)


def load_source() -> dict[str, np.ndarray]:
    source = loadmat(SOURCE, squeeze_me=True)
    required = {"all_roi", "session_index", "cell_density", "atlas_mask", "kernel"}
    missing = required.difference(source)
    if missing:
        raise KeyError(f"Missing source variables: {', '.join(sorted(missing))}")
    return source


def recompute_density(roi: np.ndarray, kernel: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    rounded = np.rint(roi).astype(int) - 1
    accumulated = np.zeros(shape, dtype=np.int64)
    np.add.at(accumulated, (rounded[:, 0], rounded[:, 1]), 1)
    return convolve2d(accumulated, kernel, mode="same")


def crop_limits(roi: np.ndarray) -> tuple[tuple[float, float], tuple[float, float]]:
    y_min, x_min = np.floor(roi.min(axis=0))
    y_max, x_max = np.ceil(roi.max(axis=0))
    return (max(0, x_min - 18), x_max + 18), (y_max + 22, max(0, y_min - 22))


def overlay_atlas(ax: plt.Axes, atlas_mask: np.ndarray, color: str, alpha: float) -> None:
    rgba = np.zeros((*atlas_mask.shape, 4), dtype=float)
    rgb = matplotlib.colors.to_rgb(color)
    rgba[..., :3] = rgb
    rgba[..., 3] = atlas_mask.astype(float) * alpha
    ax.imshow(rgba, interpolation="nearest")


def format_axis(ax: plt.Axes, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect(CORTICAL_DATA_ASPECT)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def render(source: dict[str, np.ndarray]) -> None:
    roi = np.asarray(source["all_roi"], dtype=float)
    density = np.asarray(source["cell_density"], dtype=float)
    atlas_mask = np.asarray(source["atlas_mask"], dtype=bool)
    xlim, ylim = crop_limits(roi)

    OUT_PANELS.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(2.0, 3.4))
    ax.scatter(
        roi[:, 1],
        roi[:, 0],
        s=0.28,
        c="black",
        alpha=0.16,
        linewidths=0,
        rasterized=True,
    )
    overlay_atlas(ax, atlas_mask, color="#555555", alpha=0.8)
    format_axis(ax, xlim, ylim)
    fig.savefig(OUT_PANELS / "figure_02c_locations.svg", bbox_inches="tight", transparent=True)
    fig.savefig(OUT_PANELS / "figure_02c_locations.png", bbox_inches="tight", dpi=300, transparent=True)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(2.4, 3.4))
    image = ax.imshow(density, cmap="Greys", vmin=0, vmax=1100, interpolation="nearest")
    overlay_atlas(ax, atlas_mask, color="#aaaaaa", alpha=0.75)
    format_axis(ax, xlim, ylim)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.055, pad=0.035, ticks=[0, 500, 1000])
    colorbar.set_label("Number of neurons within a\n0.3-mm-diameter circle", rotation=270, labelpad=23)
    colorbar.ax.tick_params(labelsize=8, length=2.5, width=0.7)
    colorbar.outline.set_linewidth(0.7)
    fig.savefig(OUT_PANELS / "figure_02c_density.svg", bbox_inches="tight", transparent=True)
    fig.savefig(OUT_PANELS / "figure_02c_density.png", bbox_inches="tight", dpi=300, transparent=True)
    plt.close(fig)


def main() -> int:
    source = load_source()
    roi = np.asarray(source["all_roi"], dtype=float)
    kernel = np.asarray(source["kernel"], dtype=np.int64)
    density_shape = np.asarray(source["cell_density"]).shape
    source["cell_density"] = recompute_density(roi, kernel, density_shape)
    render(source)
    print(f"Wrote Figure 2c panels to {OUT_PANELS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
