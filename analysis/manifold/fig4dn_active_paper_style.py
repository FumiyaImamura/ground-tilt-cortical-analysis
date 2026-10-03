"""cPCA and gcPCA population trajectories for Fig. 4d-h and 4j-n.

For one session:

* only neurons with activity skewness > 0.5 are used;
* activity is projected onto the first 1,000 principal components refitted on these
  neurons, and then onto cPCA (w2 = 0) or gcPCA (w2 = 0.9) axes, with 8DT as
  condition 1, Vert as condition 2 and w1 = 0.5;
* trajectories are averaged over interleaved groups of approximately five trials;
* cPCA and gcPCA panels share identical limits within each motion, and PC1, PC2 and
  PC3 share one span (equal scaling);
* eigenvector signs are arbitrary. If a previously saved projection of the same
  session is available (saved_file), the displayed axes are oriented to it. gcPCA is
  then oriented to cPCA. This changes only the display orientation.

Usage
-----
python fig4dn_active_paper_style.py [plane]

plane is zero-based. The default, 17 (Animal_09 right), is the session shown in
Fig. 4d-n.
"""

from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIGS = ROOT / "figures"
RESULTS = ROOT / "results"
sys.path.insert(0, str(HERE))

import trajectory_helpers


PLANE = int(sys.argv[1]) if len(sys.argv) > 1 else 17
W1 = 0.5
W2_C = 0.0
W2_G = 0.9
NTRI = 5
# Camera view of the Fig. 4d-n panels.
ELEV = 10
AZIM = 80
ROLL = 0
MOTIONS = [
    ("p8", "8DT", "d", "j"),
    ("rot", "T-Rot", "e", "k"),
    ("tilt", "Roll", "f", "l"),
    ("zvert", "Vert", "g", "m"),
    ("yaw", "Yaw", "h", "n"),
]


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "axes.linewidth": 0.65,
        "svg.fonttype": "none",
        "figure.dpi": 180,
        "savefig.dpi": 300,
    }
)


def draw_trajectory(ax, csp, projection, field):
    """Plot one motion with the LoadDataFromMatlab 3-D trajectory helpers."""
    values = projection[field][..., :3]
    linewidth = 1.05
    if field == "p8":
        csp.plot3d_dktc_p8(
            values,
            csp.ktp_mtypes["p8"],
            colormap="hsv",
            linewidth=linewidth,
            alpha=1,
            trials_to_combine=NTRI,
            ax=ax,
        )
    elif field == "rot":
        csp.plot3d_ktc_rot(
            values,
            csp.ktp_mtypes["rot"],
            ax=ax,
            trials_to_combine=NTRI,
            alpha=1,
            linewidth=linewidth,
        )
    elif field == "tilt":
        csp.plot3d_ktc(
            values,
            csp.ktp_mtypes["tilt"][:, :, 0],
            ax=ax,
            trials_to_combine=NTRI,
            code_direction=False,
            param_min_max=[0, 2 * np.pi],
            alpha=1,
            linewidth=linewidth,
        )
    elif field == "zvert":
        csp.plot3d_ktc(
            values,
            csp.ktp_mtypes["zvert"][:, :, 5],
            ax=ax,
            trials_to_combine=NTRI,
            code_direction=True,
            alpha=1,
            linewidth=linewidth,
        )
    elif field == "yaw":
        csp.plot3d_ktc(
            values,
            csp.ktp_mtypes["yaw"][:, :, 2],
            ax=ax,
            trials_to_combine=NTRI,
            code_direction=True,
            alpha=1,
            linewidth=linewidth,
        )
    else:
        raise KeyError(field)


def limits_and_ticks(points):
    """Paper-style limits and ticks, with EQUAL scaling on PC1, PC2 and PC3.

    All three coordinates are given the same span -- the largest of the three --
    centred on each coordinate's own midpoint. Combined with a cubic box aspect
    this makes one unit along PC1 the same visual length as one unit along PC2 or
    PC3, so the plotted trajectory keeps its true shape. Fitting each coordinate
    to its own range would stretch the narrow coordinates and distort the geometry.
    """
    values_all = np.asarray(points, dtype=float)
    spans = []
    for dim in range(3):
        v = values_all[:, dim]
        s = float(np.nanmax(v) - np.nanmin(v))
        spans.append(s if np.isfinite(s) and s > 0 else 1.0)
    common = max(spans)

    output = []
    for dim in range(3):
        values = np.asarray(points[:, dim], dtype=float)
        mid = 0.5 * (float(np.nanmin(values)) + float(np.nanmax(values)))
        lo = mid - 0.5 * common
        hi = mid + 0.5 * common
        span = common
        pad = 0.02 * span
        # Keep the padded data bounds instead of expanding the limits to the outer
        # locator ticks, so that the trajectory fills the plotting volume.
        locator = MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10])
        axis_lo = lo - pad
        axis_hi = hi + pad
        ticks = np.asarray(locator.tick_values(axis_lo, axis_hi), dtype=float)
        ticks = ticks[(ticks >= axis_lo) & (ticks <= axis_hi)]
        if ticks.size < 2:
            ticks = np.asarray([axis_lo, axis_hi], dtype=float)
        output.append((float(axis_lo), float(axis_hi), ticks))
    return output


def neuron_space_axes(pca_weights, model, rows=None, n_dim=3):
    """Map gcPCA axes back to ROI coordinates for sign comparisons."""
    weights = pca_weights if rows is None else pca_weights[rows]
    axes = np.asarray(weights @ model.components_[:, :n_dim], dtype=float)
    norms = np.linalg.norm(axes, axis=0, keepdims=True)
    return axes / np.maximum(norms, np.finfo(float).eps)


def align_model_signs(model, target_axes, pca_weights, label):
    """Resolve each model axis sign against a target in common ROI coordinates."""
    source_axes = neuron_space_axes(pca_weights, model)
    dots = np.sum(target_axes * source_axes, axis=0)
    signs = np.where(dots < 0, -1.0, 1.0)
    model.components_[:, : signs.size] *= signs[None, :]
    return {"label": label, "dots_before": dots, "applied_signs": signs}


def flattened_projection_axes(projection, n_dim=3):
    """Concatenate matching motion/time points and standardize each display axis."""
    values = np.concatenate(
        [np.asarray(projection[field][..., :n_dim]).reshape(-1, n_dim)
         for field, _, _, _ in MOTIONS],
        axis=0,
    )
    values = values - values.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(values, axis=0, keepdims=True)
    return values / np.maximum(norms, np.finfo(float).eps)


def align_model_to_saved_projection(model, active_pca, saved_projection, label):
    """Match axis signs to a previously saved projection of the same session."""
    active_projection = trajectory_helpers.mc.project(active_pca, model, dim=3)
    source = flattened_projection_axes(active_projection)
    target = flattened_projection_axes(saved_projection)
    dots = np.sum(target * source, axis=0)
    signs = np.where(dots < 0, -1.0, 1.0)
    model.components_[:, : signs.size] *= signs[None, :]
    return {"label": label, "dots_before": dots, "applied_signs": signs}


def style_axis(ax, title, panel_letter, limits, show_row_label=None):
    for setter, (lo, hi, ticks) in zip(
        (ax.set_xlim3d, ax.set_ylim3d, ax.set_zlim3d), limits
    ):
        setter(lo, hi)
    ax.set_xticks(limits[0][2])
    ax.set_yticks(limits[1][2])
    ax.set_zticks(limits[2][2])
    def endpoint_zero_formatter(ticks):
        ticks = np.asarray(ticks, dtype=float)

        def format_tick(value, _position):
            keep = np.isclose(value, ticks[0]) or np.isclose(value, ticks[-1])
            if np.any(np.isclose(ticks, 0.0)):
                keep = keep or np.isclose(value, 0.0)
            return f"{value:g}" if keep else ""

        return FuncFormatter(format_tick)

    ax.xaxis.set_major_formatter(endpoint_zero_formatter(limits[0][2]))
    ax.yaxis.set_major_formatter(endpoint_zero_formatter(limits[1][2]))
    ax.zaxis.set_major_formatter(endpoint_zero_formatter(limits[2][2]))

    ax.set_xlabel("PC1", labelpad=0)
    ax.set_ylabel("PC2", labelpad=0)
    ax.set_zlabel("PC3", labelpad=0)
    ax.tick_params(axis="x", pad=0, length=2, width=0.5)
    ax.tick_params(axis="y", pad=5, length=2, width=0.5)
    ax.tick_params(axis="z", pad=2, length=2, width=0.5)
    ax.view_init(elev=ELEV, azim=AZIM, roll=ROLL)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title(title, pad=1.5, fontweight="normal")
    ax.text2D(0.01, 0.99, panel_letter, transform=ax.transAxes,
              ha="left", va="top", fontsize=11, fontweight="bold")
    if show_row_label is not None:
        ax.text2D(-0.14, 0.52, show_row_label, transform=ax.transAxes,
                  ha="right", va="center", fontsize=11, fontstyle="italic")

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"].update(
            {"linewidth": 0.45, "linestyle": "-", "color": (0.68, 0.68, 0.68, 1)}
        )
        axis._axinfo["axisline"].update({"linewidth": 0.65, "color": (0, 0, 0, 1)})
        axis.pane.set_facecolor((1, 1, 1, 1))
        axis.pane.set_edgecolor((0.7, 0.7, 0.7, 1))


def render_single(csp, projection, field, title, panel_letter, limits, row_label=None):
    fig = plt.figure(figsize=(3.65, 3.45))
    ax = fig.add_subplot(111, projection="3d")
    draw_trajectory(ax, csp, projection, field)
    style_axis(ax, title, panel_letter, limits, None)
    # Matplotlib's 3-D tight bounding box underestimates the projected x-label
    # extent. Pull the individual-panel label inward so the exported SVG/PNG does
    # not clip "PC1" at the lower edge.
    ax.set_xlabel("PC1", labelpad=-8)
    if row_label is not None:
        fig.text(0.015, 0.52, row_label, ha="left", va="center",
                 fontsize=11, fontstyle="italic")
    fig.subplots_adjust(left=0.16, right=0.93, bottom=0.16, top=0.91)
    stem = FIGS / f"fig4{panel_letter}_active_paper_style"
    fig.savefig(stem.with_suffix(".svg"), transparent=False,
                bbox_inches="tight", pad_inches=0.08)
    fig.savefig(stem.with_suffix(".png"), transparent=False,
                bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)

    csp = trajectory_helpers.LoadDataFromMatlab.load_cell_SP_posture(
        dir_BrainRegions=r"/path/to/data/BrainRegions.mat",
        dir_AllMiceDataFile=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat",
        dir_ListOfDirectoriesForAllMice=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2_for_python.mat",
    )
    csp.load_data(PLANE)
    csp.preprocessData(
        cutoff_freq=10,
        SamplingRate_Hz=10,
        order=4,
        triError_threshold_in_mm=4,
    )
    trajectory_helpers.csp = csp

    n_stored = int(csp.cellac.shape[1])
    active_indices = np.flatnonzero(trajectory_helpers.mc.active_neurons(csp.skewness))
    n_active = int(active_indices.size)
    active_pca, _, active_pca_weights = trajectory_helpers.mc.reduce_pca(
        csp.ktc_mtypes, active_indices
    )
    if active_pca_weights is None:
        active_pca_weights = np.eye(n_active)

    active_models = {
        w2: trajectory_helpers.fit_gc(active_pca, "p8", w2) for w2 in (W2_C, W2_G)
    }

    # Eigenvector signs are arbitrary, so a trajectory can appear mirrored along an
    # axis. If a saved projection of the same session is available, orient the
    # displayed axes to it (display only).
    saved_file = (
        Path(r"/path/to/results/gcPCA_visualization_PCA1000_20250902")
        / csp.PlaneNames[PLANE]
        / "p8 vs zvert.pkl"
    )
    sign_audit = []
    if saved_file.exists():
        with saved_file.open("rb") as handle:
            saved_data = pickle.load(handle)
        saved_projection = {
            W2_C: saved_data["ktc_gcPCA"][1][0],
            W2_G: saved_data["ktc_gcPCA"][1][3],
        }
        for w2, name in ((W2_C, "cPCA to saved reference projection"),
                         (W2_G, "gcPCA to saved reference projection")):
            sign_audit.append(
                align_model_to_saved_projection(
                    active_models[w2], active_pca, saved_projection[w2], name
                )
            )

    # Then orient gcPCA to cPCA so that both rows share the same directions.
    active_c_axes = neuron_space_axes(active_pca_weights, active_models[W2_C])
    sign_audit.append(
        align_model_signs(
            active_models[W2_G], active_c_axes, active_pca_weights,
            "gcPCA to cPCA",
        )
    )

    projection = {
        w2: trajectory_helpers.mc.project(active_pca, active_models[w2], dim=3)
        for w2 in (W2_C, W2_G)
    }

    pair_limits = {}
    for field, _, _, _ in MOTIONS:
        points = np.concatenate(
            [trajectory_helpers.traj_points(projection[w2], field)
             for w2 in (W2_C, W2_G)],
            axis=0,
        )
        pair_limits[field] = limits_and_ticks(points)

    fig = plt.figure(figsize=(13.4, 6.55))
    axes = np.empty((2, 5), dtype=object)
    for row, (w2, row_label) in enumerate(((W2_C, "cPCA"), (W2_G, "gcPCA"))):
        for col, (field, title, letter_c, letter_g) in enumerate(MOTIONS):
            letter = letter_c if row == 0 else letter_g
            ax = fig.add_subplot(2, 5, row * 5 + col + 1, projection="3d")
            axes[row, col] = ax
            draw_trajectory(ax, csp, projection[w2], field)
            style_axis(
                ax,
                title,
                letter,
                pair_limits[field],
                None,
            )

    fig.text(0.012, 0.735, "cPCA", ha="left", va="center",
             fontsize=11, fontstyle="italic")
    fig.text(0.012, 0.275, "gcPCA", ha="left", va="center",
             fontsize=11, fontstyle="italic")
    fig.subplots_adjust(left=0.045, right=0.995, bottom=0.055, top=0.975,
                        wspace=0.00, hspace=0.13)
    combined = FIGS / "fig4dn_active_paper_style"
    fig.savefig(combined.with_suffix(".svg"), transparent=False)
    fig.savefig(combined.with_suffix(".png"), transparent=False)
    plt.close(fig)

    for col, (field, title, letter_c, letter_g) in enumerate(MOTIONS):
        render_single(csp, projection[W2_C], field, title, letter_c,
                      pair_limits[field], "cPCA" if col == 0 else None)
        render_single(csp, projection[W2_G], field, title, letter_g,
                      pair_limits[field], "gcPCA" if col == 0 else None)

    log = [
        "Figure 4d-h and 4j-n, active neurons",
        f"plane index: {PLANE}",
        f"session: {csp.PlaneNames[PLANE]}",
        f"stored ROIs: {n_stored}",
        f"active ROIs (activity skewness > 0.5): {n_active}",
        f"w1: {W1}",
        f"cPCA w2: {W2_C}",
        f"gcPCA w2: {W2_G}",
        f"trial-average parameter: {NTRI}",
        f"view: elevation {ELEV}, azimuth {AZIM}, roll {ROLL}",
        "axis rule: equal PC1/PC2/PC3 scaling, identical limits within each cPCA/gcPCA motion pair",
        f"orientation reference: {saved_file if saved_file.exists() else 'none'}",
        "axis signs: aligned to the saved reference projection, then gcPCA aligned to cPCA",
        "",
    ]
    for audit in sign_audit:
        log.append(audit["label"])
        log.append("  dot products before sign alignment: " + ", ".join(
            f"{value:+.6f}" for value in audit["dots_before"]
        ))
        log.append("  applied signs: " + ", ".join(
            f"{value:+.0f}" for value in audit["applied_signs"]
        ))
    log.append("")
    for field, title, _, _ in MOTIONS:
        log.append(f"matched cPCA/gcPCA limits: {title}")
        for axis_name, (lo, hi, ticks) in zip(
            ("PC1", "PC2", "PC3"), pair_limits[field]
        ):
            log.append(
                f"  {axis_name}: limits [{lo:g}, {hi:g}], ticks "
                + ", ".join(f"{v:g}" for v in ticks)
            )
    (RESULTS / "fig4dn_active_paper_style.txt").write_text("\n".join(log) + "\n", encoding="utf-8")

    print(f"session={csp.PlaneNames[PLANE]} stored={n_stored} active={n_active}")
    print(combined.with_suffix(".svg"))
    print(combined.with_suffix(".png"))


if __name__ == "__main__":
    main()
