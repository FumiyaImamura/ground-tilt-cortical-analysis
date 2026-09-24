"""Render the Figure 2d example neural, platform, and right-toe traces."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import numpy as np
from PIL import Image
from scipy.ndimage import label, median_filter


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source_data" / "example_activity" / "figure_02d_source.npz"
FIELD = ROOT / "source_data" / "example_activity" / "figure_02d_imaging_field.png"
ATLAS = ROOT / "source_data" / "example_activity" / "figure_02d_atlas_registration.png"
OUT = ROOT / "outputs" / "example_activity"
PANELS = OUT / "panels"


def atlas_overlay(size: tuple[int, int]) -> np.ndarray:
    """Return the original atlas registration as translucent fills and lines."""
    atlas = np.asarray(Image.open(ATLAS).convert("RGB"))
    white = np.all(atlas > 220, axis=2)
    components, count = label(white, np.ones((3, 3), dtype=int))
    component_sizes = [(int(np.sum(components == i)), i) for i in range(1, count + 1)]
    # The second-largest white component is the connected cortical boundary network.
    boundary_component = sorted(component_sizes, reverse=True)[1][1]
    boundaries = components == boundary_component
    # Remove atlas text while retaining the regional colour fields. The labels
    # are smaller than this local median window and are therefore replaced by
    # the surrounding registration colour.
    fill = median_filter(atlas, size=(19, 19, 1))
    colored = (fill.max(axis=2) > 55) & ~np.all(fill > 220, axis=2)

    rgba = np.zeros((*atlas.shape[:2], 4), dtype=np.uint8)
    rgba[colored, :3] = fill[colored]
    rgba[colored, 3] = 58
    rgba[boundaries, :3] = 255
    rgba[boundaries, 3] = 255

    target_width, target_height = size
    scaled_width = round(target_height * atlas.shape[1] / atlas.shape[0])
    overlay = Image.fromarray(rgba).resize((scaled_width, target_height), Image.Resampling.LANCZOS)
    left = max(0, (scaled_width - target_width) // 2)
    return np.asarray(overlay.crop((left, 0, left + target_width, target_height)))


def draw_trace_panel(axis: plt.Axes, data: np.lib.npyio.NpzFile) -> None:
    activity = data["activity"]
    platform = data["platform"]
    right_toe = data["right_toe"]
    x = np.arange(1, activity.shape[0] + 1)
    interval = 10.0
    gray = "#555555"
    direction_colors = hsv_to_rgb(np.c_[np.arange(8) / 8, np.ones(8), np.ones(8)])

    for motion, start, end in zip(data["trial_type"], data["trial_start"], data["trial_end"]):
        color = direction_colors[int(motion) - 5]
        axis.axvspan(start, end, facecolor=color, edgecolor=color, alpha=0.05, linewidth=0.8)

    axis.plot(x, activity - interval * np.arange(6), color=gray, linewidth=1.25)
    axis.plot(x, platform - interval * (6 + np.arange(2)), color=gray, linewidth=1.25)
    axis.plot(x, right_toe - interval * (8 + np.arange(3)), color=gray, linewidth=1.25)

    for row in range(6):
        axis.text(-20, -interval * row, str(row + 1), ha="right", va="center", fontsize=10)
    axis.text(-20, -interval * 6, r"$\theta$", ha="right", va="center", fontsize=11)
    axis.text(-20, -interval * 7, r"$\phi$", ha="right", va="center", fontsize=11)
    axis.text(-20, -interval * 8, "X", ha="right", va="center", fontsize=10)
    axis.text(-20, -interval * 9, "Y", ha="right", va="center", fontsize=10)
    axis.text(-20, -interval * 10, "Z", ha="right", va="center", fontsize=10)
    axis.text(-105, -interval * 6.5, "Platform", ha="right", va="center", fontsize=10)
    axis.text(-105, -interval * 9, "Right toe", ha="right", va="center", fontsize=10)
    axis.text(activity.shape[0] + 12, 3, r"$\Delta F/F$", fontsize=11)

    vectors = {
        5: (0, 1), 6: (-1, 1), 7: (-1, 0), 8: (-1, -1),
        9: (0, -1), 10: (1, -1), 11: (1, 0), 12: (1, 1),
    }
    for motion, start, end in zip(data["trial_type"], data["trial_start"], data["trial_end"]):
        dx, dy = vectors[int(motion)]
        center = (float(start) + float(end)) / 2
        color = direction_colors[int(motion) - 5]
        axis.annotate(
            "", xy=(center + 13 * dx, 14 + 7 * dy), xytext=(center - 13 * dx, 14 - 7 * dy),
            arrowprops={"arrowstyle": "-|>", "color": color, "lw": 2.2, "mutation_scale": 12},
            annotation_clip=False,
        )

    # Time scale: 4 s at 10 Hz.
    y_scale = -108
    axis.plot([activity.shape[0] - 40, activity.shape[0]], [y_scale, y_scale], color="black", lw=1.5)
    axis.text(activity.shape[0] - 20, y_scale - 4, "4 s", ha="center", va="top", fontsize=10)
    axis.set_xlim(-10, activity.shape[0] + 70)
    axis.set_ylim(-114, 18)
    axis.axis("off")


def main() -> None:
    PANELS.mkdir(parents=True, exist_ok=True)
    data = np.load(SOURCE)

    fig = plt.figure(figsize=(12, 6.15))
    grid = fig.add_gridspec(1, 2, width_ratios=(0.35, 0.65), wspace=0.34)
    image_axis = fig.add_subplot(grid[0, 0])
    trace_axis = fig.add_subplot(grid[0, 1])

    field = Image.open(FIELD).convert("RGB").resize((330, 480), Image.Resampling.LANCZOS)
    image_axis.imshow(field)
    image_axis.imshow(atlas_overlay((330, 480)))
    image_axis.axis("off")
    image_axis.plot([230, 329], [510, 510], color="black", lw=2.0, clip_on=False)
    image_axis.text(280, 502, "1 mm", ha="center", va="bottom", fontsize=10)

    draw_trace_panel(trace_axis, data)
    fig.text(0.01, 0.98, "d", ha="left", va="top", fontsize=14, fontweight="bold")
    fig.savefig(PANELS / "figure_02d.svg", bbox_inches="tight", facecolor="white")
    fig.savefig(PANELS / "figure_02d.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"Wrote Figure 2d panels to {PANELS}")


if __name__ == "__main__":
    main()
