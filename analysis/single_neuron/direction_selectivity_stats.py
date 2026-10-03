"""Statistics for Fig. 2f,g from the saved direction-selectivity results.

Reads the file written by direction_selectivity.m and selects the
direction-selective neurons (activity skewness > 0.5, ANOVA p < 0.05 and direction
selectivity score > 0.1).

Fig. 2f: chi-square goodness-of-fit test of the preferred-direction counts against a
uniform distribution over the eight directions, followed by a two-sided exact binomial
test of each direction against 1/8, Bonferroni-corrected over the eight directions.

Fig. 2g: Pearson's chi-square test of independence between area and selectivity on the
counts pooled over sessions, followed by a two-sided Fisher's exact test of each area
against the pooled other seven areas, Bonferroni-corrected over the eight areas. An area
is left out of a session when it holds fewer than 50 neurons in that session.

display_symbol gives the asterisks drawn on the panels, which mark only directions or
areas whose fraction is above the uniform or the pooled-other fraction.

Usage:  python direction_selectivity_stats.py path/to/direction_selectivity_<time>.mat
Writes figure_02f_direction_binomial_tests.csv and figure_02g_region_fisher_tests.csv
into ../results.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.stats import binomtest, chi2_contingency, chisquare, fisher_exact


OUT = Path(__file__).resolve().parent.parent / "results"

# Order of the area masks in RegLog_plane (RegionNames_type in BrainRegions_with_RSC.mat).
# Whisker, UL, LL and Trunk are S1b, S1fl, S1hl and S1tr.
REGION_NAMES_SAVED = ("PPC", "M2", "M1", "Trunk", "Whisker", "V1", "UL", "LL", "RSC")
DISPLAY_REGIONS = ("M2", "M1", "Whisker", "UL", "LL", "Trunk", "PPC", "RSC")
DIRECTION_DEGREES = np.arange(45, 361, 45)

SKEW_THRESHOLD = 0.5
P_THRESHOLD = 0.05
SCORE_THRESHOLD = 0.1
MIN_REGION_NEURONS = 50


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def symbol(p_value: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "n.s."


def main(source_file: Path) -> None:
    source = loadmat(source_file, simplify_cells=True)
    direction_results = source["dir_sel"]
    skew_by_session = source["skew_plane"]
    region_masks_by_session = source["RegLog_plane"]

    direction_counts = np.zeros(8, dtype=int)
    valid_counts = np.zeros((len(direction_results), len(REGION_NAMES_SAVED)), dtype=float)
    total_counts = np.zeros_like(valid_counts)

    for session, (result, skew, region_masks) in enumerate(
        zip(direction_results, skew_by_session, region_masks_by_session, strict=True)
    ):
        skew = np.asarray(skew, dtype=float).ravel()
        keep_active = skew > SKEW_THRESHOLD
        p_values = np.asarray(result["pval"], dtype=float).ravel()[keep_active]
        scores = np.asarray(result["score"], dtype=float).ravel()[keep_active]
        directions = np.asarray(result["direction"], dtype=int).ravel()[keep_active]
        selected = (p_values < P_THRESHOLD) & (scores > SCORE_THRESHOLD)
        direction_counts += np.bincount(directions[selected] - 1, minlength=8)

        for region, mask in enumerate(region_masks):
            mask = np.asarray(mask, dtype=bool).ravel()[keep_active]
            n_region = int(mask.sum())
            if n_region >= MIN_REGION_NEURONS:
                valid_counts[session, region] = int(np.sum(selected[mask]))
                total_counts[session, region] = n_region
            else:
                valid_counts[session, region] = np.nan
                total_counts[session, region] = np.nan

    # Fig. 2f
    n_direction_selective = int(direction_counts.sum())
    goodness_of_fit = chisquare(direction_counts)
    direction_rows: list[dict[str, object]] = []
    for degree, count in zip(DIRECTION_DEGREES, direction_counts, strict=True):
        p_two_sided = float(
            binomtest(int(count), n_direction_selective, 1.0 / 8.0, alternative="two-sided").pvalue
        )
        two_sided_corrected = min(1.0, 8.0 * p_two_sided)
        two_sided_symbol = symbol(two_sided_corrected)
        direction_rows.append(
            {
                "direction_deg": int(degree),
                "count": int(count),
                "total_direction_selective_neurons": n_direction_selective,
                "p_raw_two_sided": p_two_sided,
                "p_bonferroni_8_two_sided": two_sided_corrected,
                "two_sided_symbol": two_sided_symbol,
                "display_symbol": two_sided_symbol
                if count > n_direction_selective / 8.0
                else "n.s.",
                "omnibus_chi_square": float(goodness_of_fit.statistic),
                "omnibus_df": len(direction_counts) - 1,
                "omnibus_p": float(goodness_of_fit.pvalue),
            }
        )

    # Fig. 2g
    region_order = [REGION_NAMES_SAVED.index(name) for name in DISPLAY_REGIONS]
    valid_pool = np.nansum(valid_counts[:, region_order], axis=0).astype(int)
    total_pool = np.nansum(total_counts[:, region_order], axis=0).astype(int)
    invalid_pool = total_pool - valid_pool
    contingency = np.column_stack((valid_pool, invalid_pool))
    chi2, chi2_p, chi2_df, _ = chi2_contingency(contingency, correction=False)

    region_rows: list[dict[str, object]] = []
    for index, region in enumerate(DISPLAY_REGIONS):
        table = np.asarray(
            [
                [valid_pool[index], invalid_pool[index]],
                [valid_pool.sum() - valid_pool[index], invalid_pool.sum() - invalid_pool[index]],
            ]
        )
        p_two_sided = float(fisher_exact(table, alternative="two-sided").pvalue)
        two_sided_corrected = min(1.0, 8.0 * p_two_sided)
        two_sided_symbol = symbol(two_sided_corrected)
        other_fraction = float(
            (valid_pool.sum() - valid_pool[index])
            / (total_pool.sum() - total_pool[index])
        )
        region_rows.append(
            {
                "region": region,
                "direction_selective_neurons": int(valid_pool[index]),
                "all_eligible_neurons": int(total_pool[index]),
                "pooled_fraction": float(valid_pool[index] / total_pool[index]),
                "other_areas_fraction": other_fraction,
                "p_raw_two_sided": p_two_sided,
                "p_bonferroni_8_two_sided": two_sided_corrected,
                "two_sided_symbol": two_sided_symbol,
                "display_symbol": two_sided_symbol
                if valid_pool[index] / total_pool[index] > other_fraction
                else "n.s.",
                "omnibus_chi_square": float(chi2),
                "omnibus_df": int(chi2_df),
                "omnibus_p": float(chi2_p),
            }
        )

    write_csv(OUT / "figure_02f_direction_binomial_tests.csv", direction_rows)
    write_csv(OUT / "figure_02g_region_fisher_tests.csv", region_rows)
    print(
        f"Figure 2f: {n_direction_selective} direction-selective neurons, "
        f"chi-square={goodness_of_fit.statistic:.8g}, df={len(direction_counts) - 1}, "
        f"p={goodness_of_fit.pvalue:.8g}"
    )
    print(f"Figure 2g: chi-square={chi2:.8g}, df={chi2_df}, p={chi2_p:.8g}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]))
