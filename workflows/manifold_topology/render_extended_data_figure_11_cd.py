"""Render Extended Data Figure 11c-d from bundled source values."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import wilcoxon


MOTION_KEYS = ("rot", "p8", "zvert")
MOTION_LABELS = {"rot": "T-Rot", "p8": "8DT", "zvert": "Vert"}
MOTION_COLORS = {"rot": "#7E57C2", "p8": "#5C8FC9", "zvert": "#F5A623"}
DISPLAY_DIMS = (3, 8)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7.5,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "svg.hashsalt": "ground-tilt-ed11",
        "figure.dpi": 150,
    }
)


def sem(values: np.ndarray, axis: int) -> np.ndarray:
    """Sample SEM along ``axis``, ignoring non-finite values."""

    n = np.sum(np.isfinite(values), axis=axis)
    return np.nanstd(values, axis=axis, ddof=1) / np.sqrt(np.maximum(n, 1))


def star(p_value: float) -> str:
    if p_value < 1e-3:
        return "***"
    if p_value < 1e-2:
        return "**"
    if p_value < 5e-2:
        return "*"
    return ""


def corrected_tests(
    values: np.ndarray,
    motion_indices: dict[str, int],
    tested_w2: tuple[int, ...],
) -> dict[tuple[str, int], float]:
    """One-sided paired Wilcoxon tests versus Vert with panel-wise Bonferroni."""

    raw: list[float] = []
    keys: list[tuple[str, int]] = []
    vert = motion_indices["zvert"]
    for motion in ("rot", "p8"):
        mi = motion_indices[motion]
        for wi in tested_w2:
            x = values[wi, :, mi]
            y = values[wi, :, vert]
            valid = np.isfinite(x) & np.isfinite(y)
            raw.append(float(wilcoxon(x[valid], y[valid], alternative="greater").pvalue))
            keys.append((motion, wi))
    adjusted = np.minimum(np.asarray(raw) * len(raw), 1.0)
    return dict(zip(keys, adjusted, strict=True))


def load_source(source: Path) -> dict[str, np.ndarray]:
    with np.load(source, allow_pickle=False) as archive:
        required = {"dims", "w2", "motions", "lri", "maxlife_normalized"}
        missing = required.difference(archive.files)
        if missing:
            raise ValueError(f"missing arrays in {source}: {sorted(missing)}")
        data = {name: np.asarray(archive[name]) for name in required}

    dims = tuple(int(value) for value in data["dims"])
    motions = tuple(str(value) for value in data["motions"])
    if not all(dimension in dims for dimension in DISPLAY_DIMS):
        raise ValueError(f"expected dimensions {DISPLAY_DIMS}, found {dims}")
    if motions != MOTION_KEYS:
        raise ValueError(f"expected motion order {MOTION_KEYS}, found {motions}")
    if not np.array_equal(data["w2"], np.asarray([0.0, 0.5, 1.0, 3.0])):
        raise ValueError(f"unexpected w2 values: {data['w2']}")
    expected_shape = (4, len(dims), 20, 3)
    for key in ("lri", "maxlife_normalized"):
        if data[key].shape != expected_shape:
            raise ValueError(f"{key} has shape {data[key].shape}, expected {expected_shape}")
        if not np.all(np.isfinite(data[key])):
            raise ValueError(f"{key} contains non-finite values")
    return data


def render_panel(
    values: np.ndarray,
    dims: tuple[int, ...],
    w2: np.ndarray,
    motions: tuple[str, ...],
    *,
    ylabel: str,
    tested_w2: tuple[int, ...],
    ylim: tuple[float, float],
    yticks: tuple[float, ...],
    legend_location: str,
    output_svg: Path,
    output_png: Path,
) -> dict[int, dict[tuple[str, int], float]]:
    dim_indices = [dims.index(dimension) for dimension in DISPLAY_DIMS]
    motion_indices = {motion: motions.index(motion) for motion in motions}
    selected = values[:, dim_indices, :, :]
    panel_tests: dict[int, dict[tuple[str, int], float]] = {}

    fig, axes = plt.subplots(1, 2, figsize=(5.15, 2.35), sharey=True)
    fig.subplots_adjust(left=0.12, right=0.985, bottom=0.22, top=0.88, wspace=0.25)
    y_span = ylim[1] - ylim[0]

    for column, dimension in enumerate(DISPLAY_DIMS):
        ax = axes[column]
        panel = selected[:, column, :, :]
        adjusted = corrected_tests(panel, motion_indices, tested_w2)
        panel_tests[dimension] = adjusted
        means: dict[str, np.ndarray] = {}
        errors: dict[str, np.ndarray] = {}

        for motion in motions:
            mi = motion_indices[motion]
            means[motion] = np.mean(panel[:, :, mi], axis=1)
            errors[motion] = sem(panel[:, :, mi], axis=1)
            ax.plot(
                w2,
                means[motion],
                color=MOTION_COLORS[motion],
                linewidth=1.3,
                label=MOTION_LABELS[motion],
            )
            ax.fill_between(
                w2,
                means[motion] - errors[motion],
                means[motion] + errors[motion],
                color=MOTION_COLORS[motion],
                alpha=0.25,
                linewidth=0,
            )

        # Place labels just above their corresponding curves.  When the T-Rot
        # and 8DT labels would collide, lift the higher label by one text row.
        for wi in tested_w2:
            positions: list[tuple[str, float, str]] = []
            for motion in ("rot", "p8"):
                mark = star(adjusted[(motion, wi)])
                if mark:
                    positions.append(
                        (
                            motion,
                            float(means[motion][wi] + errors[motion][wi] + 0.025 * y_span),
                            mark,
                        )
                    )
            if len(positions) == 2 and abs(positions[0][1] - positions[1][1]) < 0.06 * y_span:
                high = int(positions[1][1] > positions[0][1])
                motion, y_value, mark = positions[high]
                positions[high] = (motion, y_value + 0.055 * y_span, mark)
            for _, y_value, mark in positions:
                ax.text(
                    w2[wi],
                    min(y_value, ylim[1] - 0.035 * y_span),
                    mark,
                    color="black",
                    fontsize=6.2,
                    fontweight="bold",
                    ha="center",
                    va="bottom",
                )

        ax.set_title(f"{dimension} dimensions")
        ax.set_xticks(w2)
        ax.set_xticklabels(("0", "0.5", "1", "3"))
        ax.set_xlabel(r"$w_2$")
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        if column == 0:
            ax.set_ylabel(ylabel)

    axes[1].legend(frameon=False, loc=legend_location, handlelength=1.7, handletextpad=0.35)
    output_svg.parent.mkdir(parents=True, exist_ok=True)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_svg, metadata={"Date": None})
    fig.savefig(output_png, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return panel_tests


def main() -> None:
    repository = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=repository
        / "source_data"
        / "manifold_topology"
        / "fig5df_dimensions_3_5_8_data.npz",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository / "outputs" / "manifold_topology" / "panels",
    )
    parser.add_argument(
        "--preview-dir",
        type=Path,
        default=repository / "outputs" / "manifold_topology" / "panels",
    )
    args = parser.parse_args()

    data = load_source(args.source)
    dims = tuple(int(value) for value in data["dims"])
    motions = tuple(str(value) for value in data["motions"])

    jobs = (
        (
            "c",
            data["maxlife_normalized"],
            "Normalized max H1 lifetime",
            (1, 2, 3),
            (0.23, 1.30),
            (0.4, 0.8, 1.2),
            "center right",
        ),
        (
            "d",
            data["lri"],
            "LRI",
            (0, 1, 2, 3),
            (0.16, 0.87),
            (0.2, 0.4, 0.6, 0.8),
            "upper right",
        ),
    )
    for letter, values, ylabel, tested_w2, ylim, yticks, legend_location in jobs:
        svg = args.output_dir / f"extended_data_figure_11{letter}.svg"
        png = args.preview_dir / f"extended_data_figure_11{letter}.png"
        tests = render_panel(
            values,
            dims,
            data["w2"],
            motions,
            ylabel=ylabel,
            tested_w2=tested_w2,
            ylim=ylim,
            yticks=yticks,
            legend_location=legend_location,
            output_svg=svg,
            output_png=png,
        )
        print(f"wrote {svg}")
        print(f"wrote {png}")
        for dimension in DISPLAY_DIMS:
            entries = ", ".join(
                f"{MOTION_LABELS[motion]} w2={data['w2'][wi]:g}: p_adj={p_value:.8g}"
                for (motion, wi), p_value in tests[dimension].items()
            )
            print(f"  {letter}, {dimension} dimensions: {entries}")


if __name__ == "__main__":
    main()
