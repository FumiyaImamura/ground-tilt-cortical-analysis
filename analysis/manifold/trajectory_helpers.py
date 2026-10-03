"""Helper functions for the 3-D cPCA/gcPCA trajectories in Fig. 4d-h and 4j-n.

Imported by population_trajectories.py. gcPCA uses 8DT as condition 1 and Vert as
condition 2 with w1 = 0.5. Trajectories are averaged over interleaved groups of
approximately NTRI trials.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import manifold_core as mc

CFDIR = r"/path/to/python_helpers"
if CFDIR in sys.path: sys.path.remove(CFDIR)
sys.path.insert(0, CFDIR)
import LoadDataFromMatlab, utility_functions as utf

W1 = 0.5
NTRI = 5

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 7, 'axes.titlesize': 7, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 6,
    'axes.linewidth': 0.6, 'svg.fonttype': 'none', 'figure.dpi': 150,
})


# --------------------------------------------------------------------- helpers
def fit_gc(sub, cond1, w2):
    """gcPCA with `cond1` as condition 1 and Vert as condition 2."""
    g = mc.gcPCA()
    a_in = g.CovarianceAndNormalize_intraTrial(_flat_if_p8(sub[cond1], cond1))
    a_ac = g.CovarianceAndNormalize_interTrial(_flat_if_p8(sub[cond1], cond1))
    z_in = g.CovarianceAndNormalize_intraTrial(sub['zvert'])
    g.fit((1 - W1) * a_in - W1 * a_ac - w2 * z_in)
    return g


def _flat_if_p8(v, fld):
    """8DT is (dir, trial, time, cell). The eight directions are concatenated into
    one 24 s trajectory per trial before the covariances are taken."""
    if fld != 'p8':
        return v
    return utf.multiBlockCat(v, dim2squish=0, dim2concat=2).squeeze(axis=0)


def traj_points(proj, fld):
    v = proj[fld][..., :3]
    if fld == 'p8':
        k = utf.jumping_mean(v, dim=1, num_to_combine=NTRI)
        return utf.multiBlockCat(k, dim2squish=[1, 0], dim2concat=2).squeeze().reshape(-1, 3)
    k = utf.jumping_mean(v, dim=0, num_to_combine=NTRI)
    return utf.multiBlockCat(k, dim2squish=0, dim2concat=1).squeeze().reshape(-1, 3)
