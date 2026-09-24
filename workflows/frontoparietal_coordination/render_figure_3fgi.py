"""Render Figure 3f, 3g, and 3i from the bundled decoder predictions.

The input is the compact MATLAB result saved by
``TB_classi_cPCA_20250706_analysis.m``.  It contains the true and decoded
directions for the frontal and parietal populations at four cPCA weights.
The manuscript uses weight 0.5.

Figure 3f and 3g follow the original plotting block. The bundled predictions
have already been circularly averaged over 1.0-1.5 s within each trial.  A
column-normalized confusion matrix is therefore calculated directly from the
saved values for each of 20 sessions, and the session matrices are averaged
with equal weight.

Figure 3i uses the same saved trial-level predictions.  The joint distribution
of frontal and parietal decoding offsets is shown after
subtracting the product of its marginals.  One-sided binomial tests compare
each observed joint-bin count with its independence probability, followed by
Bonferroni correction across the 64 bins.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.ticker import ScalarFormatter
import numpy as np
from scipy.io import loadmat
from scipy.stats import binomtest


REPO_ROOT = Path(__file__).resolve().parents[2]
SOURCE = REPO_ROOT / "source_data" / "frontoparietal_coordination" / "Reg_comparison.mat"
PARULA_SOURCE = REPO_ROOT / "source_data" / "frontoparietal_coordination" / "parula_256.csv"
OUT_PANELS = REPO_ROOT / "outputs" / "frontoparietal_coordination" / "panels"
OUT_RESULTS = REPO_ROOT / "outputs" / "frontoparietal_coordination" / "results"

WEIGHT_INDEX = 1
PLOT_ORDER = np.array([4, 5, 6, 7, 0, 1, 2, 3])
CONFUSION_PLOT_ORDER = np.array([1, 2, 3, 4, 5, 6, 7, 0])
DEGREES_8 = np.array([-180, -135, -90, -45, 0, 45, 90, 135])
ARROWS = ["↗", "→", "↘", "↓", "↙", "←", "↖", "↑"]
ARROW_COLORS = ["#4b00d1", "#1455d9", "#12cbd4", "#10a52b", "#9ee63a", "#ff9b32", "#ef1e10", "#df13d4"]


PARULA = ListedColormap(np.loadtxt(PARULA_SOURCE, delimiter=","), name="matlab_parula")


def wrap_bin(angle: np.ndarray) -> np.ndarray:
    wrapped = np.angle(np.exp(1j * angle))
    return np.mod(np.rint(np.rad2deg(wrapped) / 45.0), 8).astype(int)


def circular_mean(values: np.ndarray, axis: int = 0) -> np.ndarray:
    return np.angle(np.mean(np.exp(1j * values), axis=axis))


def load_predictions() -> tuple[dict[str, list[np.ndarray]], np.ndarray]:
    data = loadmat(SOURCE)
    rec = data["recdata"][0, 0]
    resp = data["respdata"][0, 0]
    weights = data["cPCA_weights"].ravel()
    predictions: dict[str, list[np.ndarray]] = {
        "frontal": [],
        "parietal": [],
        "true": [],
    }
    for plane in range(rec["M1_M2"].shape[0]):
        frontal = rec["M1_M2"][plane, WEIGHT_INDEX].ravel()
        parietal = rec["PPC_S1"][plane, WEIGHT_INDEX].ravel()
        true_frontal = resp["M1_M2"][plane, WEIGHT_INDEX].ravel()
        true_parietal = resp["PPC_S1"][plane, WEIGHT_INDEX].ravel()
        if not np.allclose(true_frontal, true_parietal):
            raise ValueError(f"True directions differ between regions in session {plane + 1}")
        predictions["frontal"].append(frontal)
        predictions["parietal"].append(parietal)
        predictions["true"].append(true_frontal)
    return predictions, weights


def session_mean_confusion(predictions: dict[str, list[np.ndarray]], region: str) -> np.ndarray:
    matrices: list[np.ndarray] = []
    for decoded, true in zip(predictions[region], predictions["true"]):
        decoded_bin = wrap_bin(decoded)
        true_bin = wrap_bin(true)
        counts = np.zeros((8, 8), dtype=float)
        np.add.at(counts, (decoded_bin, true_bin), 1)
        matrices.append(counts / counts.sum(axis=0, keepdims=True))
    return np.mean(matrices, axis=0)


def style_direction_ticks(ax: plt.Axes) -> None:
    ax.set_xticks(np.arange(8), ARROWS)
    ax.set_yticks(np.arange(8), ARROWS)
    for label, color in zip(ax.get_xticklabels(), ARROW_COLORS):
        label.set_color(color)
        label.set_fontsize(10)
        label.set_fontweight("bold")
    for label, color in zip(ax.get_yticklabels(), ARROW_COLORS):
        label.set_color(color)
        label.set_fontsize(10)
        label.set_fontweight("bold")


def save_figure(fig: plt.Figure, stem: str) -> None:
    OUT_PANELS.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PANELS / f"{stem}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(OUT_PANELS / f"{stem}.png", bbox_inches="tight", dpi=300, transparent=True)
    plt.close(fig)


def render_confusion(matrix: np.ndarray, title: str, stem: str, vmax: float) -> None:
    matrix = matrix[np.ix_(CONFUSION_PLOT_ORDER, CONFUSION_PLOT_ORDER)]
    fig, ax = plt.subplots(figsize=(2.15, 2.15))
    image = ax.imshow(matrix, cmap="gray", vmin=0.05, vmax=vmax, origin="upper", interpolation="nearest")
    style_direction_ticks(ax)
    ax.set_xlabel("True direction", labelpad=3)
    ax.set_ylabel("Predicted direction", labelpad=3)
    ax.set_title(title, pad=6)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
    colorbar.set_label("Probability", rotation=270, labelpad=12)
    colorbar.outline.set_linewidth(0.6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    save_figure(fig, stem)


def duplicate_circular_endpoint(matrix: np.ndarray) -> np.ndarray:
    return np.block([[matrix, matrix[:, :1]], [matrix[:1, :], matrix[:1, :1]]])


def joint_excess_and_pvalues(
    predictions: dict[str, list[np.ndarray]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    frontal_offset = []
    parietal_offset = []
    for frontal, parietal, true in zip(
        predictions["frontal"], predictions["parietal"], predictions["true"]
    ):
        frontal_offset.append(wrap_bin(frontal - true).ravel())
        parietal_offset.append(wrap_bin(parietal - true).ravel())
    frontal_bin = np.concatenate(frontal_offset)
    parietal_bin = np.concatenate(parietal_offset)
    counts = np.zeros((8, 8), dtype=int)
    np.add.at(counts, (parietal_bin, frontal_bin), 1)
    total = int(counts.sum())
    joint = counts / total
    expected = np.outer(joint.sum(axis=1), joint.sum(axis=0))
    excess = joint - expected
    p_raw = np.empty((8, 8), dtype=float)
    for row in range(8):
        for column in range(8):
            p_raw[row, column] = binomtest(
                int(counts[row, column]),
                total,
                float(expected[row, column]),
                alternative="greater",
            ).pvalue
    p_bonferroni = np.minimum(1.0, p_raw * 64)
    return counts, joint, expected, excess, p_bonferroni, total


def render_joint_panel(excess: np.ndarray, p_bonferroni: np.ndarray) -> None:
    ordered_excess = excess[np.ix_(PLOT_ORDER, PLOT_ORDER)]
    ordered_sig = p_bonferroni[np.ix_(PLOT_ORDER, PLOT_ORDER)] < 1e-3
    excess9 = duplicate_circular_endpoint(ordered_excess)
    sig9 = duplicate_circular_endpoint(ordered_sig.astype(int))

    fig, axes = plt.subplots(1, 2, figsize=(4.4, 2.05), gridspec_kw={"width_ratios": [1.18, 1]})
    ax = axes[0]
    image = ax.imshow(
        excess9,
        cmap=PARULA,
        origin="upper",
        extent=(-202.5, 202.5, 202.5, -202.5),
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_xlabel("Offset (Frontal)")
    ax.set_ylabel("Offset (Parietal)")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, -2))
    colorbar.formatter = formatter
    colorbar.update_ticks()
    colorbar.set_label("P(X,Y)-P(X)P(Y)", rotation=270, labelpad=14)
    colorbar.outline.set_linewidth(0.6)

    ax = axes[1]
    ax.imshow(
        sig9,
        cmap=ListedColormap(["white", "black"]),
        vmin=0,
        vmax=1,
        origin="upper",
        extent=(-202.5, 202.5, 202.5, -202.5),
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-180, -90, 0, 90, 180])
    ax.set_xlabel("Offset (Frontal)")
    ax.set_ylabel("Offset (Parietal)")
    ax.set_title(r"$p<0.001$", fontsize=9, fontstyle="italic", pad=5)
    for panel_ax in axes:
        panel_ax.tick_params(length=2, width=0.6, labelsize=8)
        for spine in panel_ax.spines.values():
            spine.set_linewidth(0.6)
    fig.tight_layout(w_pad=1.0)
    save_figure(fig, "figure_03i")


def main() -> int:
    predictions, weights = load_predictions()
    if not np.isclose(weights[WEIGHT_INDEX], 0.5):
        raise ValueError(f"Expected cPCA weight 0.5, found {weights[WEIGHT_INDEX]}")

    parietal_confusion = session_mean_confusion(predictions, "parietal")
    frontal_confusion = session_mean_confusion(predictions, "frontal")
    render_confusion(parietal_confusion, "Parietal cortex", "figure_03f", vmax=0.42)
    render_confusion(frontal_confusion, "Frontal cortex", "figure_03g", vmax=0.30)

    counts, joint, expected, excess, p_bonferroni, sample_count = joint_excess_and_pvalues(predictions)
    render_joint_panel(excess, p_bonferroni)

    OUT_RESULTS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUT_RESULTS / "figure_03fgi_values.npz",
        parietal_confusion=parietal_confusion,
        frontal_confusion=frontal_confusion,
        parietal_confusion_display=parietal_confusion[
            np.ix_(CONFUSION_PLOT_ORDER, CONFUSION_PLOT_ORDER)
        ],
        frontal_confusion_display=frontal_confusion[
            np.ix_(CONFUSION_PLOT_ORDER, CONFUSION_PLOT_ORDER)
        ],
        joint_counts=counts,
        joint_probability=joint,
        independence_probability=expected,
        excess_probability=excess,
        binomial_p_bonferroni_64=p_bonferroni,
        cPCA_weight=weights[WEIGHT_INDEX],
        sample_count=sample_count,
    )

    print(f"Wrote Figure 3f, 3g, and 3i panels to {OUT_PANELS}")
    print(f"Wrote plotted values to {OUT_RESULTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
