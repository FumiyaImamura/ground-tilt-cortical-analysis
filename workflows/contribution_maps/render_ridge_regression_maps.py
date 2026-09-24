"""Render ridge-regression contribution and neuron-density maps.

The source bundle contains the saved full-cohort map arrays, display masks,
atlas overlay, color limits, and MATLAB colormaps used for the manuscript.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.cm import ScalarMappable
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE = REPO_ROOT / "source_data" / "contribution_maps" / "ridge_regression_maps.npz"
OUT_PANELS = REPO_ROOT / "outputs" / "contribution_maps" / "panels"

# CSP.plotScaleBar calibrates the atlas image against its physical dimensions:
#   x pixels/mm = 911 * (1.77 / 18.89)
#   y pixels/mm = 904 * (1.77 / 15.02)
# MATLAB then applies daspect([x_pixels_per_mm / y_pixels_per_mm, 1, 1]).
# Matplotlib's numeric aspect is the displayed y-unit/x-unit ratio, so the
# same value reproduces the manuscript cortical geometry.
CORTICAL_DATA_ASPECT = (911.0 / 18.89) / (904.0 / 15.02)


def load_source() -> dict[str, np.ndarray]:
    with np.load(SOURCE) as archive:
        source = {key: archive[key] for key in archive.files}
    required = {
        "map_names",
        "relative_maps",
        "relative_masks",
        "relative_high",
        "normalized_maps",
        "normalized_masks",
        "normalized_high",
        "density_maps",
        "density_masks",
        "density_high",
        "atlas_image",
        "atlas_mask",
        "jet_colormap",
        "parula_colormap",
    }
    missing = required.difference(source)
    if missing:
        raise KeyError(f"Missing source variables: {', '.join(sorted(missing))}")
    return source


def compose_map(
    values: np.ndarray,
    mask: np.ndarray,
    high: float,
    atlas_image: np.ndarray,
    atlas_mask: np.ndarray,
    colors: np.ndarray,
) -> np.ndarray:
    scaled = np.nan_to_num(values / high, nan=0.0, posinf=1.0, neginf=0.0)
    indices = np.rint(np.clip(scaled, 0.0, 1.0) * (len(colors) - 1)).astype(int)
    rgb = colors[indices]
    rgb = np.where(mask[..., None], rgb, 0.0)
    alpha = np.clip(atlas_mask, 0.0, 1.0)[..., None]
    return np.clip(rgb * (1.0 - alpha) + atlas_image * alpha, 0.0, 1.0)


def draw_map(
    fig: plt.Figure,
    ax: plt.Axes,
    values: np.ndarray,
    mask: np.ndarray,
    high: float,
    title: str,
    colorbar_label: str,
    source: dict[str, np.ndarray],
    colors: np.ndarray,
) -> None:
    rgb = compose_map(
        values,
        mask,
        high,
        source["atlas_image"],
        source["atlas_mask"],
        colors,
    )
    ax.imshow(rgb, interpolation="nearest", aspect=CORTICAL_DATA_ASPECT)
    ax.set_title(title, fontsize=9, fontweight="normal", pad=3)
    ax.set_axis_off()
    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(0.0, high), cmap=ListedColormap(colors)),
        ax=ax,
        fraction=0.040,
        pad=0.025,
    )
    colorbar.set_label(colorbar_label, rotation=270, labelpad=10, fontsize=7)
    colorbar.ax.tick_params(labelsize=6, length=2, width=0.5, pad=1.5)
    colorbar.outline.set_linewidth(0.5)


def save(fig: plt.Figure, stem: str) -> None:
    OUT_PANELS.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PANELS / f"{stem}.svg", bbox_inches="tight", transparent=False)
    fig.savefig(OUT_PANELS / f"{stem}.png", bbox_inches="tight", dpi=300, transparent=False)
    plt.close(fig)


def render_row(
    source: dict[str, np.ndarray],
    indices: list[int],
    values_key: str,
    masks_key: str,
    high_key: str,
    label: str,
    colors_key: str,
    stem: str,
    titles: list[str] | None = None,
) -> None:
    fig, axes = plt.subplots(1, len(indices), figsize=(2.15 * len(indices), 3.7))
    axes = np.atleast_1d(axes)
    colors = source[colors_key]
    names = source["map_names"] if titles is None else np.asarray(titles)
    for ax, index, title in zip(axes, indices, names):
        draw_map(
            fig,
            ax,
            source[values_key][..., index],
            source[masks_key][..., index],
            float(source[high_key][index]),
            str(title),
            label,
            source,
            colors,
        )
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.92, wspace=0.28)
    save(fig, stem)


def render_extended_data(
    source: dict[str, np.ndarray],
    values_key: str,
    masks_key: str,
    high_key: str,
    label: str,
    stem: str,
) -> None:
    fig = plt.figure(figsize=(10.8, 7.2))
    grid = fig.add_gridspec(2, 5, left=0.015, right=0.985, bottom=0.02, top=0.98, wspace=0.30, hspace=0.16)
    colors = source["parula_colormap"]
    names = [str(value) for value in source["map_names"]]
    for index in range(9):
        row, column = (0, index) if index < 5 else (1, index - 5)
        ax = fig.add_subplot(grid[row, column])
        draw_map(
            fig,
            ax,
            source[values_key][..., index],
            source[masks_key][..., index],
            float(source[high_key][index]),
            names[index],
            label,
            source,
            colors,
        )
    save(fig, stem)


def main() -> int:
    source = load_source()
    render_row(
        source,
        [0, 1, 2, 3],
        "relative_maps",
        "relative_masks",
        "relative_high",
        "Relative contribution",
        "jet_colormap",
        "figure_02k",
        ["Forelimb", "Hindlimb", "Trunk", "Tail"],
    )
    render_row(
        source,
        [4],
        "relative_maps",
        "relative_masks",
        "relative_high",
        "Relative contribution",
        "jet_colormap",
        "figure_02l",
        ["Ground tilt"],
    )
    render_row(
        source,
        [5, 6, 7, 8],
        "relative_maps",
        "relative_masks",
        "relative_high",
        "Relative contribution",
        "jet_colormap",
        "figure_02n",
        ["Eigenposture 1", "Eigenposture 2", "Eigenposture 3", "Eigenposture 4"],
    )
    render_extended_data(
        source,
        "normalized_maps",
        "normalized_masks",
        "normalized_high",
        "Normalized neuron density",
        "extended_data_figure_08a",
    )
    render_extended_data(
        source,
        "density_maps",
        "density_masks",
        "density_high",
        "Neuron density",
        "extended_data_figure_08b",
    )
    print(f"Wrote ridge-regression map panels to {OUT_PANELS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
