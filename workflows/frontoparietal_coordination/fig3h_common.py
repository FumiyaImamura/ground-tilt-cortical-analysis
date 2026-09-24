"""Shared circular-bin definitions for Figure 3h and Extended Data Figure 10k."""

import numpy as np


NBIN = 8
BIN_CENTERS = np.array([-180, -135, -90, -45, 0, 45, 90, 135])
PLOT_ORDER = np.array([4, 5, 6, 7, 0, 1, 2, 3])


def as_plot_vector(values: np.ndarray) -> np.ndarray:
    """Return the eight circular bins in plotting order with the endpoint repeated."""

    ordered = values[PLOT_ORDER]
    return np.concatenate([ordered, ordered[:1]])
