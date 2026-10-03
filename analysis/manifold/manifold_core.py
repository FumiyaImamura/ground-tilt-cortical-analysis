"""Core routines for the cPCA/gcPCA, persistent-homology and LRI analyses.

Conventions
-----------
* motion keys: 'rot' = T-Rot, 'p8' = 8DT, 'zvert' = Vert, 'tilt' = Roll, 'yaw' = Yaw
* neurons: active neurons, i.e. activity skewness > 0.5 (active_neurons())
* gcPCA objective, T-Rot as condition 1 and Vert as condition 2:
      S(w1,w2) = (1-w1) Sigma_within(rot) - w1 Sigma_across(rot) - w2 Sigma_within(zvert)
  w1 = 0.5 throughout. w1 = w2 = 0 reduces to standard PCA.
* point cloud for persistent homology: utility_functions.jumping_mean averages
  interleaved, non-overlapping groups of approximately 5 trials, using all time
  points and the first 5 gcPCA dimensions.
* LRI = integral(lambda_1) / sum_{i<=10} integral(lambda_i) of the H1 persistence
  landscape.
"""
import os, sys
import numpy as np

CFDIR = r"/path/to/python_helpers"
if CFDIR in sys.path:
    sys.path.remove(CFDIR)
sys.path.insert(0, CFDIR)
import utility_functions as utf                      # noqa: E402
from generalized_cPCA import gcPCA                   # noqa: E402
from ripser import ripser                            # noqa: E402

SKEW_THRESHOLD = 0.5    # active neurons: activity skewness > 0.5
W1 = 0.5
W2_LIST = (0.0, 0.5, 1.0, 3.0)
PH_DIM = 5
TRI_MEAN = 5
K_MAX = 10
N_GRID = 400
PCA_CAP = 1000          # paper: project onto the first 1000 PCs (only if n > 1000)
P8_INDICE = [14]        # 8DT uses the single time bin at peak tilt, so its point
                        # cloud has 8 points per trial


# ------------------------------------------------------------- neuron selection
def active_neurons(skewness, threshold=SKEW_THRESHOLD):
    """Boolean mask of the neurons analysed in the manuscript (skewness > threshold)."""
    return np.asarray(skewness, dtype=float).ravel() > threshold


# ----------------------------------------------------------------- point clouds
def point_cloud(ktc, fld, tri_mean=TRI_MEAN, dim=PH_DIM):
    """Trial-averaged point cloud for persistent homology.

    For finite tri_mean this delegates to utility_functions.jumping_mean, which
    averages interleaved, non-overlapping groups of approximately tri_mean trials.
    """
    k = ktc[..., :dim]
    if fld == 'p8':
        k = k[..., P8_INDICE, :]
        k = utf.jumping_mean(k, dim=1, num_to_combine=tri_mean)
        return utf.multiBlockCat(k, dim2squish=[1, 0], dim2concat=2).squeeze()
    k = utf.jumping_mean(k, dim=0, num_to_combine=tri_mean)
    return utf.multiBlockCat(k, dim2squish=0, dim2concat=1).squeeze()


# ---------------------------------------------------------- persistent homology
def _finite(D):
    D = np.asarray(D, dtype=float)
    return D[np.isfinite(D).all(axis=1)]


def persistence_landscape(D, k_max=K_MAX, n_grid=N_GRID):
    """H1 persistence landscape: the k_max largest tent functions on an n_grid grid."""
    D = _finite(D)
    if D.size == 0:
        return np.linspace(0.0, 1.0, n_grid), np.zeros((k_max, n_grid))
    R = float(np.max(D[:, 1]))
    T = np.linspace(0.0, R, int(n_grid))
    H = np.maximum(np.minimum(T[None, :] - D[:, 0][:, None],
                              D[:, 1][:, None] - T[None, :]), 0.0)
    Hs = np.sort(H, axis=0)[::-1, :]
    if Hs.shape[0] < k_max:
        Hs = np.vstack([Hs, np.zeros((k_max - Hs.shape[0], Hs.shape[1]))])
    return T, Hs[:k_max, :]


def lri(dgm_h1, k_max=K_MAX):
    """Landscape ringness index; scale invariant, so no diagram rescaling needed."""
    T, L = persistence_landscape(dgm_h1, k_max=k_max)
    den = float(np.trapz(L.sum(axis=0), T))
    if den <= 0:
        return np.nan
    return float(np.trapz(L[0], T)) / den


def max_h1_lifetime(dgm_h1):
    d = _finite(dgm_h1)
    if d.size == 0:
        return np.nan
    return float(np.max(d[:, 1] - d[:, 0]))


# ------------------------------------------------------------------ gcPCA fits
def slice_neurons(ktc_mtypes, idx, fields=None):
    it = ktc_mtypes.items() if fields is None else (
        (f, ktc_mtypes[f]) for f in fields)
    return {f: np.ascontiguousarray(v[..., idx]) for f, v in it}


def reduce_pca(ktc_mtypes, idx, cap=PCA_CAP, seed=0, fields=None):
    """Optionally pre-project a neuron subset onto its first `cap` PCs.

    Follows the paper's rule for subset analyses: the initial PCA step is applied
    only when the subset contains more than `cap` neurons; smaller subsets are
    analysed directly (truncating to `cap` components would be a full-rank
    rotation anyway, which leaves the gcPCA subspace unchanged).

    Returns (projected_dict, n_neurons, W) with W the (n, cap) projection matrix,
    or None when no PCA was applied.
    """
    raw = slice_neurons(ktc_mtypes, idx, fields)
    n = raw['rot'].shape[-1]
    if n <= cap:
        return raw, n, None
    from sklearn.decomposition import PCA
    flat = np.concatenate([v.reshape(-1, n) for v in raw.values()], axis=0)
    p = PCA(n_components=cap, svd_solver='randomized', random_state=seed)
    p.fit(flat)
    W = p.components_.T                      # (n, cap)
    return {f: v @ W for f, v in raw.items()}, n, W


def fit_gcpca(sub, w1=W1, w2=0.0):
    g = gcPCA()
    rot_in = g.CovarianceAndNormalize_intraTrial(sub['rot'])
    if w1 == 0.0 and w2 == 0.0:
        g.fit(rot_in)                        # standard PCA
        return g
    rot_ac = g.CovarianceAndNormalize_interTrial(sub['rot'])
    zv_in = g.CovarianceAndNormalize_intraTrial(sub['zvert'])
    g.fit((1 - w1) * rot_in - w1 * rot_ac - w2 * zv_in)
    return g


def project(sub, g, dim=PH_DIM):
    out = {}
    for f, v in sub.items():
        out[f] = v @ g.components_[:, :dim]
    return out


def split_halves(d, p8_axis=1):
    """Contiguous train/test split of the trial axis (first half, second half).
    Returns (train, test)."""
    tr, te = {}, {}
    for f, v in d.items():
        ax = p8_axis if f.startswith('p8') else 0
        mid = v.shape[ax] // 2
        tr[f] = np.take(v, np.arange(mid), axis=ax)
        te[f] = np.take(v, np.arange(mid, v.shape[ax]), axis=ax)
    return tr, te


def ring_metrics(sub, w1=W1, w2=0.0, motions=('rot', 'p8', 'zvert'), dim=PH_DIM):
    """gcPCA -> first `dim` dimensions -> Ripser H1 -> {LRI, max H1 lifetime}."""
    g = fit_gcpca(sub, w1=w1, w2=w2)
    proj = project(sub, g, dim=dim)
    res = {}
    for f in motions:
        cloud = point_cloud(proj[f], f, dim=dim)
        dgms = ripser(cloud)['dgms']
        res[f] = dict(lri=lri(dgms[1]), maxlife=max_h1_lifetime(dgms[1]))
    return res, g, proj
