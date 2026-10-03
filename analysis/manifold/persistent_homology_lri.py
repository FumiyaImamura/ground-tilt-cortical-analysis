"""Persistent homology and LRI of gcPCA embeddings for every session (Fig. 5d,f).

For each session, the trial-arranged activity of active neurons (activity skewness
> 0.5) is projected onto its first 1,000 principal components. For each w2 in W2_LIST,
gcPCA is fitted with T-Rot as condition 1, Vert as condition 2 and w1 = 0.5. The T-Rot,
8DT and Vert activity is projected onto the first 5 gcPCA dimensions and averaged over
interleaved groups of approximately 5 trials, and the H1 persistence diagram (Ripser)
gives the LRI and the maximum H1 lifetime.

Usage:  python persistent_homology_lri.py [first_plane] [last_plane]
Writes one .npz per session into ../results. A_is has the dimensions
subset x w2 x motion x metric.
"""
import os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import manifold_core as mc

CFDIR = r"/path/to/python_helpers"
if CFDIR in sys.path: sys.path.remove(CFDIR)
sys.path.insert(0, CFDIR)
import LoadDataFromMatlab

RESULTS = os.path.abspath(os.path.join(HERE, '..', 'results'))
os.makedirs(RESULTS, exist_ok=True)

SUBSETS = ['active']
W2_LIST = list(mc.W2_LIST)
MOTIONS = ['rot', 'p8', 'zvert']
METRICS = ['lri', 'maxlife']


def run_plane(csp, plane, planename):
    t0 = time.perf_counter()
    csp.load_data(plane)
    csp.preprocessData(cutoff_freq=10, SamplingRate_Hz=10, order=4,
                       triError_threshold_in_mm=4)
    n_cells = csp.cellac.shape[1]
    active = mc.active_neurons(csp.skewness)
    if active.size != n_cells:
        raise ValueError(f'skewness length {active.size} != n_cells {n_cells}')
    sub, n_active, _ = mc.reduce_pca(csp.ktc_mtypes, np.flatnonzero(active))
    A_is = np.full((len(SUBSETS), len(W2_LIST), len(MOTIONS), len(METRICS)), np.nan)
    for wi, w2 in enumerate(W2_LIST):
        res, _, _ = mc.ring_metrics(sub, w1=mc.W1, w2=w2, motions=MOTIONS)
        for j, m in enumerate(MOTIONS):
            for k, met in enumerate(METRICS):
                A_is[0, wi, j, k] = res[m][met]
    out = os.path.join(RESULTS, f'plane{plane:02d}_manifold.npz')
    np.savez_compressed(
        out, A_is=A_is, subsets=np.asarray(SUBSETS),
        cond_A_is=np.asarray([str(w) for w in W2_LIST]), w2_list=np.asarray(W2_LIST),
        motions=np.asarray(MOTIONS), metrics=np.asarray(METRICS),
        n_active=n_active, planename=planename)
    print(f'  active neurons={n_active}  T-Rot LRI at w2={W2_LIST}: '
          + ', '.join(f'{v:.3f}' for v in A_is[0, :, 0, 0])
          + f'  [{time.perf_counter() - t0:.0f}s] -> {out}', flush=True)


if __name__ == '__main__':
    first = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    last = int(sys.argv[2]) if len(sys.argv) > 2 else 19
    csp = LoadDataFromMatlab.load_cell_SP_posture(
        dir_BrainRegions=r'/path/to/data/BrainRegions.mat',
        dir_AllMiceDataFile=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat",
        dir_ListOfDirectoriesForAllMice=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2_for_python.mat")
    T0 = time.perf_counter()
    for plane in range(first, last + 1):
        pn = csp.PlaneNames[plane]
        print(f'=== plane {plane}: {pn} ===', flush=True)
        run_plane(csp, plane, pn)
    print(f'ALL DONE in {(time.perf_counter() - T0) / 60:.1f} min', flush=True)
