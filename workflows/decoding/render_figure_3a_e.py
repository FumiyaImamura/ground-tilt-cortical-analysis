"""Render Figure 3a-e from session-level decoding results.

Classification panels are centered on the analytic probability 3/8. Error
bars show the sample standard deviation of the finite session values divided
by sqrt(20).

The significance symbols are calculated from the bundled session values.
For each motion and area, the area value is paired by session with the mean of
the other available areas and tested with a one-sided Wilcoxon signed-rank
test.  The eight area tests within a motion are Bonferroni-corrected.  This
procedure is used to place the significance symbols.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "Arial"
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.stats import wilcoxon


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "source_data" / "decoding" / "basedata.mat"
OUT = ROOT / "outputs" / "decoding"
PANELS = OUT / "panels"
RESULTS = OUT / "results"

MOTIONS = (
    ("p8", "a", "8DT", True),
    ("rot", "b", "T-Rot", True),
    ("tilt", "c", "Roll", False),
    ("yaw", "d", "Yaw", False),
    ("zvert", "e", "Vert", False),
)
REGIONS = ("M2", "M1", "Whisker", "UL", "LL", "Trunk", "PPC", "RSC")
LABELS = ("M2", "M1", "Whisker", "FL", "HL", "Trunk", "PPC", "RSC")
CHANCE = 3.0 / 8.0
def values_by_region(basedata, motion: str, classification: bool) -> dict[str, np.ndarray]:
    sub = getattr(basedata, motion)
    values = {}
    for region in REGIONS:
        raw = np.asarray(getattr(sub, region), dtype=float).ravel()
        if classification:
            raw = raw - CHANCE
        values[region] = raw
    return values


def summary(values: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    means, sems, ns = [], [], []
    for region in REGIONS:
        raw = values[region]
        finite = raw[np.isfinite(raw)]
        means.append(float(np.mean(finite)))
        sems.append(float(np.std(finite, ddof=1) / np.sqrt(raw.size)))
        ns.append(int(finite.size))
    return np.asarray(means), np.asarray(sems), np.asarray(ns)


def symbol(corrected_p: float) -> str:
    if corrected_p < 0.001:
        return "***"
    if corrected_p < 0.01:
        return "**"
    if corrected_p < 0.05:
        return "*"
    return ""


def paired_region_tests(values: dict[str, np.ndarray]) -> list[dict[str, object]]:
    rows = []
    for region in REGIONS:
        target = values[region]
        others = np.column_stack([values[other] for other in REGIONS if other != region])
        rest = np.nanmean(others, axis=1)
        valid = np.isfinite(target) & np.isfinite(rest)
        test = wilcoxon(
            target[valid],
            rest[valid],
            alternative="greater",
            method="approx" if int(valid.sum()) >= 20 else "exact",
            correction=True if int(valid.sum()) >= 20 else False,
        )
        raw_p = float(test.pvalue)
        corrected_p = min(1.0, 8.0 * raw_p)
        rows.append(
            {
                "region": region,
                "n_pairs": int(valid.sum()),
                "raw_p": raw_p,
                "bonferroni_factor": 8,
                "corrected_p": corrected_p,
                "symbol": symbol(corrected_p),
            }
        )
    return rows


def render_panel(
    motion: str,
    letter: str,
    title: str,
    values: dict[str, np.ndarray],
    tests: list[dict[str, object]],
) -> None:
    means, sems, _ = summary(values)
    x = np.arange(len(REGIONS))
    fig, ax = plt.subplots(figsize=(3.15, 3.25))
    ax.bar(
        x,
        means,
        0.68,
        yerr=sems,
        color="0.68",
        edgecolor="0.2",
        linewidth=0.65,
        error_kw={"ecolor": "0.2", "lw": 0.75, "capsize": 2},
    )
    for xi, mean, sem, row in zip(x, means, sems, tests):
        if row["symbol"]:
            ax.text(
                xi,
                mean + sem + 0.035 * means.max(),
                str(row["symbol"]),
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
            )
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS, rotation=45, ha="right", fontsize=8.5)
    ax.set_title(title, fontsize=14, pad=7)
    ax.set_ylabel("$\\Delta$ Decoding score", fontsize=10.5)
    ax.set_ylim(bottom=0)
    ax.spines[["top", "right"]].set_visible(False)
    y0 = -0.26 * means.max()
    ax.plot([2 - 0.35, 5 + 0.35], [y0, y0], color="black", lw=0.9, clip_on=False)
    ax.text(3.5, y0 - 0.05 * means.max(), "S1", ha="center", va="top", fontsize=9, clip_on=False)
    fig.subplots_adjust(left=0.25, right=0.98, top=0.86, bottom=0.34)
    for suffix in ("svg", "png"):
        fig.savefig(PANELS / f"figure_03{letter}.{suffix}", dpi=300, transparent=False)
    plt.close(fig)


def main() -> None:
    PANELS.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    basedata = loadmat(SOURCE, squeeze_me=True, struct_as_record=False)["basedata"]

    summary_rows = []
    test_rows = []
    for motion, letter, title, classification in MOTIONS:
        values = values_by_region(basedata, motion, classification)
        means, sems, ns = summary(values)
        tests = paired_region_tests(values)
        for region, label, mean, sem, n, test in zip(REGIONS, LABELS, means, sems, ns, tests):
            summary_rows.append(
                {
                    "panel": f"3{letter}",
                    "motion": motion,
                    "title": title,
                    "region": region,
                    "display_label": label,
                    "mean": float(mean),
                    "sem_sqrt20": float(sem),
                    "n_valid_sessions": int(n),
                }
            )
            test_rows.append(
                {
                    "panel": f"3{letter}",
                    "motion": motion,
                    **test,
                }
            )
        render_panel(motion, letter, title, values, tests)

    for path, rows in (
        (RESULTS / "figure_03a_e_values.csv", summary_rows),
        (RESULTS / "figure_03a_e_region_tests.csv", test_rows),
    ):
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    print(f"Wrote Figure 3a-e panels to {PANELS}")
    print(f"Wrote summary values and statistical tests to {RESULTS}")


if __name__ == "__main__":
    main()
