"""Replot real and shuffled decoding relative to the analytic null.

The same
analytic constant (3/8 for 8DT and T-Rot; zero for Roll, Yaw, and Vert) is
subtracted from real and shuffled scores in every area and mouse. The pooled
empirical mean is not subtracted. Negative
analytic-null-centered shuffle bars are retained.

Use `--output-dir` to select a destination other than the default output directory.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import scipy.io as sio
from scipy.stats import wilcoxon


REPO = Path(__file__).resolve().parents[2]
SOURCE = REPO / "source_data" / "decoding"
RESULTS = REPO / "outputs" / "decoding" / "results"
MOTIONS = [
    ("p8", "8DT", True, 3 / 8),
    ("rot", "T-Rot", True, 3 / 8),
    ("tilt", "Roll", False, 0.0),
    ("yaw", "Yaw", False, 0.0),
    ("zvert", "Vert", False, 0.0),
]
REGIONS = ["M2", "M1", "Whisker", "FL", "HL", "Trunk", "PPC", "RSC"]
SOURCE_REGION = {region: region for region in REGIONS} | {"FL": "UL", "HL": "LL"}
REAL_COLOR = "#5E5E5E"
SHUFFLE_COLOR = "#C7C7C7"


def read_csv(name):
    with (RESULTS / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(name, rows):
    with (RESULTS / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mouse_scores(perrep, motion, region):
    scores = np.asarray(getattr(getattr(perrep, motion), SOURCE_REGION[region]), dtype=float)
    if scores.shape != (20, 10):
        raise ValueError(f"Unexpected source shape: {motion}/{region}: {scores.shape}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        session = np.nanmean(scores, axis=1)
        return np.nanmean(session.reshape(10, 2), axis=1)


def sem(values):
    return float(np.std(values, ddof=1) / np.sqrt(values.size))


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def prepare_results():
    raw_rows = read_csv("empirical_shuffle_stats.csv")
    raw_lookup = {(row["motion_key"], row["region"]): row for row in raw_rows}
    raw_summary = {row["motion_key"]: row for row in read_csv("empirical_shuffle_summary.csv")}
    real = sio.loadmat(SOURCE / "perrep_scores.mat", struct_as_record=False, squeeze_me=True)["perrep"]
    shuffle = sio.loadmat(SOURCE / "perrep_scores_shuffle.mat", struct_as_record=False, squeeze_me=True)["perrep"]
    centered_rows, summary_rows, mouse_rows = [], [], []
    raw_errors, pooled_errors, sem_errors, contrast_errors, p_errors, saved_p_errors = [], [], [], [], [], []
    baseline_errors, centered_mean_errors = [], []

    for motion, title, classification, analytic in MOTIONS:
        by_region = {}
        for region in REGIONS:
            r = mouse_scores(real, motion, region)
            s = mouse_scores(shuffle, motion, region)
            valid = np.isfinite(r) & np.isfinite(s)
            by_region[region] = (r[valid], s[valid], np.flatnonzero(valid) + 1)

        pooled = np.concatenate([values[1] for values in by_region.values()])
        pooled_mean = float(np.mean(pooled))
        offset = analytic
        old_summary = raw_summary[motion]
        assert pooled.size == int(old_summary["pooled_region_mouse_values"])
        pooled_errors.append(abs(pooled_mean - float(old_summary["pooled_shuffle_mean"])))
        baseline_errors.append(abs(offset - float(old_summary["analytic_reference"])))
        motion_centered_shuffle = []

        for region in REGIONS:
            r, s, mice = by_region[region]
            original = raw_lookup[(motion, region)]
            baseline_errors.append(abs(offset - float(original["analytic_reference"])))
            assert r.size == int(original["n_mice"])
            assert int(np.sum(r > s)) == int(original["positive_mice"])
            for value, key in [(np.mean(r), "real_mean"), (sem(r), "real_mouse_sem"),
                               (np.mean(s), "shuffle_mean"), (sem(s), "shuffle_mouse_sem")]:
                raw_errors.append(abs(float(value) - float(original[key])))

            r_centered, s_centered = r - offset, s - offset
            centered_mean_errors.extend([
                abs(float(np.mean(r_centered)) - (float(np.mean(r)) - analytic)),
                abs(float(np.mean(s_centered)) - (float(np.mean(s)) - analytic)),
            ])
            contrast_errors.extend(np.abs((r_centered - s_centered) - (r - s)))
            sem_errors.extend([abs(sem(r_centered) - sem(r)), abs(sem(s_centered) - sem(s))])
            p_raw = float(wilcoxon(r, s, alternative="greater").pvalue)
            p_centered = float(wilcoxon(r_centered, s_centered, alternative="greater").pvalue)
            p_errors.append(abs(p_centered - p_raw))
            saved_p_errors.append(abs(p_raw - float(original["p_raw_paired_one_sided"])))
            p8 = min(1.0, 8 * p_centered)
            p40 = min(1.0, 40 * p_centered)
            assert abs(p8 - float(original["p_bonferroni_8"])) < 5.1e-11
            assert abs(p40 - float(original["p_bonferroni_40"])) < 5.1e-11
            symbol = stars(p8)
            assert symbol == original["star_bonferroni_8"]
            motion_centered_shuffle.append(float(np.mean(s_centered)))
            centered_rows.append({
                "motion_key": motion, "motion": title, "region": region,
                "n_mice": r.size, "positive_mice": int(np.sum(r > s)),
                "real_mean_raw": float(np.mean(r)), "shuffle_mean_raw": float(np.mean(s)),
                "pooled_empirical_null": pooled_mean, "subtracted_analytic_null": offset,
                "real_mean_centered": float(np.mean(r_centered)),
                "real_mouse_sem": sem(r_centered),
                "shuffle_mean_centered": float(np.mean(s_centered)),
                "shuffle_mouse_sem": sem(s_centered),
                "analytic_reference_raw": analytic, "analytic_reference_centered": analytic - offset,
                "p_raw_paired_one_sided": p_centered, "p_bonferroni_8": p8,
                "p_bonferroni_40": p40, "star_bonferroni_8": symbol,
            })
            for mouse, rv, sv, rc, sc in zip(mice, r, s, r_centered, s_centered):
                mouse_rows.append({
                    "motion_key": motion, "motion": title, "region": region, "mouse_index": int(mouse),
                    "real_raw": float(rv), "shuffle_raw": float(sv), "subtracted_analytic_null": offset,
                    "real_centered": float(rc), "shuffle_centered": float(sc),
                    "paired_difference": float(rc - sc),
                })
        summary_rows.append({
            "motion_key": motion, "motion": title, "pooled_region_mouse_values": pooled.size,
            "pooled_empirical_null": pooled_mean, "subtracted_analytic_null": offset,
            "analytic_reference_centered": analytic - offset,
            "pooled_centered_shuffle_mean": float(np.mean(pooled - offset)),
            "negative_centered_shuffle_bars": int(np.sum(np.asarray(motion_centered_shuffle) < 0)),
            "smallest_centered_shuffle_mean": min(motion_centered_shuffle),
            "largest_centered_shuffle_mean": max(motion_centered_shuffle),
        })

    # The stored raw CSV has eight decimal places; all residuals must be rounding only.
    assert max(raw_errors) <= 5.1e-9
    assert max(pooled_errors) <= 5.1e-9
    assert max(baseline_errors) == 0
    assert max(centered_mean_errors) <= 1e-14
    assert max(saved_p_errors) <= 5.1e-11
    assert max(contrast_errors) <= 1e-14
    assert max(sem_errors) <= 1e-14
    assert max(p_errors) == 0
    write_csv("empirical_shuffle_centered_stats.csv", centered_rows)
    write_csv("empirical_shuffle_centered_summary.csv", summary_rows)
    write_csv("empirical_shuffle_centered_mouse_values.csv", mouse_rows)
    return centered_rows


def plot_panel(ax, rows, title, classification, compact=False):
    x = np.arange(len(REGIONS))
    width = 0.36
    r = np.array([float(row["real_mean_centered"]) for row in rows])
    s = np.array([float(row["shuffle_mean_centered"]) for row in rows])
    re = np.array([float(row["real_mouse_sem"]) for row in rows])
    se = np.array([float(row["shuffle_mouse_sem"]) for row in rows])
    reference = float(rows[0]["analytic_reference_centered"])
    assert reference == 0
    assert all(float(row["analytic_reference_centered"]) == 0 for row in rows)
    assert all(float(row["subtracted_analytic_null"]) == (3 / 8 if classification else 0)
               for row in rows)
    bounds = np.concatenate([r - re, r + re, s - se, s + se, [0, reference]])
    span = float(np.ptp(bounds))
    lower = min(float(np.min(bounds)) - 0.035 * span, -0.045 * span)
    upper = float(np.max(bounds)) + (0.14 if compact else 0.30) * span
    bars = []
    for positions, means, errors, color, label in [
        (x - width / 2, r, re, REAL_COLOR, "Real"),
        (x + width / 2, s, se, SHUFFLE_COLOR, "Shuffle"),
    ]:
        group = ax.bar(positions, means, width, bottom=0.0, yerr=errors,
                       color=color, edgecolor="0.2", linewidth=0.45 if compact else 0.55,
                       label=label, error_kw=dict(ecolor="0.2", lw=0.55 if compact else 0.65,
                                                  capsize=1.2 if compact else 1.5))
        bars.extend(group)
        for region, bar in zip(REGIONS, group):
            bar.set_gid(f"{label.lower().replace(' ', '-')}-{region}")
    ax.axhline(reference, color="black", lw=0.8, zorder=3)
    for xi, rv, error, row in zip(x, r, re, rows):
        symbol = row["star_bonferroni_8"]
        if symbol != "n.s.":
            ax.text(xi, rv + error + 0.014 * span, symbol, ha="center", va="bottom",
                    fontsize=6.5 if compact else 8, fontweight="bold")
    ax.set_title(title, fontsize=10 if compact else 12, pad=6 if compact else 8)
    ax.set_ylabel(r"$\Delta$ Decoding accuracy", fontsize=8 if compact else 10)
    ax.set_xticks(x, REGIONS, rotation=45, ha="right", rotation_mode="anchor")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out", length=2.5 if compact else 3, labelsize=7 if compact else 9)
    ax.set_ylim(lower, upper)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=4))
    assert ax.get_ylim()[0] < 0 < ax.get_ylim()[1]
    assert all(bar.get_y() == 0 for bar in bars)
    assert lower <= np.min(bounds) and upper >= np.max(bounds)
    return len(bars)


def build_figures(rows, output_dir):
    lookup = {(row["motion_key"], row["region"]): row for row in rows}
    output_dir.mkdir(parents=True, exist_ok=True)
    for motion, title, classification, _ in MOTIONS:
        with plt.rc_context({"svg.fonttype": "none", "font.family": "sans-serif", "font.size": 9,
                             "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8}):
            fig, ax = plt.subplots(figsize=(4.45, 3.75))
            panel_rows = [lookup[(motion, region)] for region in REGIONS]
            plot_panel(ax, panel_rows, title, classification)
            ax.legend(handles=[
                Patch(facecolor=REAL_COLOR, edgecolor="0.2", label="Real"),
                Patch(facecolor=SHUFFLE_COLOR, edgecolor="0.2", label="Shuffle"),
            ], frameon=False, fontsize=7.5, loc="upper left")
            fig.subplots_adjust(left=0.17, right=0.98, top=0.90, bottom=0.25)
            for extension in ("svg", "png"):
                fig.savefig(output_dir / f"empirical_null_{motion}.{extension}", bbox_inches="tight",
                            **({"dpi": 240} if extension == "png" else {"metadata": {"Date": None}}))
            plt.close(fig)

    with plt.rc_context({"svg.fonttype": "path", "font.family": "Arial", "font.size": 8,
                         "axes.linewidth": 0.75, "xtick.major.width": 0.75, "ytick.major.width": 0.75}):
        fig, axes = plt.subplots(3, 2, figsize=(7.4, 9.0))
        axes = axes.ravel()
        for index, (motion, title, classification, _) in enumerate(MOTIONS):
            ax = axes[index]
            plot_panel(ax, [lookup[(motion, region)] for region in REGIONS],
                       title, classification, compact=True)
            ax.text(-0.15, 1.08, chr(ord("A") + index), transform=ax.transAxes,
                    fontsize=11, fontweight="bold", ha="left", va="top")
        axes[5].axis("off")
        axes[5].legend(handles=[
            Patch(facecolor=REAL_COLOR, edgecolor="0.2", label="Real"),
            Patch(facecolor=SHUFFLE_COLOR, edgecolor="0.2", label="Shuffle"),
        ], frameon=False, loc="upper center", fontsize=9)
        axes[5].text(0.5, 0.35,
                     "Asterisks: real > shuffle,\none-sided paired Wilcoxon\nsigned-rank tests; Bonferroni\nx8 within motion",
                     ha="center", va="center", fontsize=8.5, linespacing=1.4, transform=axes[5].transAxes)
        fig.subplots_adjust(left=0.10, right=0.98, top=0.96, bottom=0.08, hspace=0.55, wspace=0.32)
        for extension in ("svg", "png"):
            fig.savefig(output_dir / f"empirical_null_all_motions.{extension}", bbox_inches="tight",
                        **({"dpi": 300} if extension == "png" else {"metadata": {"Date": None}}))
        plt.close(fig)
    print(f"Wrote centered empirical-null figures to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO / "outputs" / "decoding" / "panels",
    )
    args = parser.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    prepare_results()
    # Draw from saved results, so the exact plotted values are traceable on disk.
    build_figures(read_csv("empirical_shuffle_centered_stats.csv"), args.output_dir)


if __name__ == "__main__":
    main()
