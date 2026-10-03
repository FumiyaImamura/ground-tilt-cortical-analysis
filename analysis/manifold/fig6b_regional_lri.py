"""Cross-validated LRI within cortical subdivisions (Fig. 6b).

For each session and each group (all areas, PPC, S1 and M1+M2), the active neurons
(activity skewness > 0.5) inside the registered region masks are analysed with two arms
evaluated on the same held-out trials:

  gcPCA          PCA(n_components=1000) fit on the continuous activity of the selected
                 neurons and applied to the trial-arranged data, which are split into
                 contiguous trial halves. gcPCA (T-Rot as condition 1, Vert as
                 condition 2, w1 = 0.5, w2 = 3) is fit on the first half and applied to
                 the held-out half.
  standard PCA   the trial-arranged data are split into halves first, then
                 PCA(n_components=1000) is fit on the flattened training T-Rot data and
                 applied to the held-out half.

The held-out T-Rot trials are averaged together, reduced to 5 dimensions, and the H1 LRI
is taken with k_max = 10. Fig. 6b shows gcPCA LRI minus standard-PCA LRI. Groups with
fewer than 50 neurons are skipped. The region masks are written by export_region_masks.m.

Usage:  python fig6b_regional_lri.py [--first 0] [--last 19] [--force]
Writes ../results/fig6b_m1m2_plane??.npz
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from scipy.io import loadmat

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manifold_core as mc  # noqa: E402

CF = r"/path/to/python_helpers"
if CF in sys.path:
    sys.path.remove(CF)
sys.path.insert(0, CF)
import LoadDataFromMatlab, utility_functions as utf  # noqa: E402
from PCA import PCA  # noqa: E402

RESULTS = HERE.parent / "results"
RESULTS.mkdir(parents=True, exist_ok=True)
# Written by export_region_masks.m.
REGION_MASK_FILE = RESULTS / "region_masks_with_RSC.mat"

PCADIM = 1000
W1, W2, DIM = 0.5, 3.0, 5
REGLIST = ["M1", "M2", "PPC", "LL", "UL", "Whisker", "Trunk", "RSC"]
REGION_KEYS = ("AllAreas", "PPC", "S1", "M1M2")


def observed_lri(csp, mask):
    """Observed arm: gcPCA LRI on the held-out trials."""
    pca = PCA(n_components=PCADIM, data=csp.cellac_cont[..., mask])
    ktc = {f: utf.contig(pca.transform(v[..., mask]))
           for f, v in csp.ktc_mtypes.items() if f in ('rot', 'zvert')}
    tr, te = mc.split_halves(ktc)
    g = mc.fit_gcpca(tr, w1=W1, w2=W2)
    cloud = mc.point_cloud(te['rot'] @ g.components_[:, :DIM], 'rot',
                           tri_mean=np.inf, dim=DIM)
    return mc.lri(mc.ripser(cloud)['dgms'][1])


def null_lri(csp, mask):
    """Null arm: split first, then standard PCA on the flattened training T-Rot."""
    raw = {f: np.ascontiguousarray(v[..., mask])
           for f, v in csp.ktc_mtypes.items() if f in ('rot', 'zvert')}
    tr, te = mc.split_halves(raw)
    rot_flat = utf.multiBlockCat(tr['rot'], dim2squish=0, dim2concat=1).squeeze(axis=0)
    pca = PCA(n_components=PCADIM, data=rot_flat)
    cloud = mc.point_cloud(pca.transform(te['rot'])[..., :DIM], 'rot',
                           tri_mean=np.inf, dim=DIM)
    return mc.lri(mc.ripser(cloud)['dgms'][1])


def load_region_masks():
    """(plane -> (n_roi, 8) bool array in REGLIST order, plane names)."""
    m = loadmat(str(REGION_MASK_FILE), squeeze_me=True, simplify_cells=True)
    names = [str(x) for x in np.atleast_1d(m['regNamesChar'])]
    planes = [str(x) for x in np.atleast_1d(m['planeNamesChar'])]
    order = [names.index(r) for r in REGLIST]
    masks = [np.asarray(M).astype(bool)[:, order] for M in m['maskCell']]
    return masks, planes


def load_registered_masks() -> tuple[list[np.ndarray], list[str]]:
    """Return masks in REGLIST order and the 20 session names."""
    masks, plane_names = load_region_masks()
    return [np.asarray(mask, dtype=bool) for mask in masks], list(plane_names)


def masks_for_plane(plane: int, n_cells: int | None = None) -> dict[str, np.ndarray]:
    """Construct the four groups used in Fig. 6b.

    S1 is the union of LL, UL, Whisker and Trunk. M1M2 is the union of the registered
    M1 and M2 masks.
    """
    masks, _ = load_registered_masks()
    matrix = masks[plane]
    if n_cells is not None and matrix.shape[0] != n_cells:
        raise ValueError(f"region-mask length {matrix.shape[0]} != n_cells {n_cells}")
    index = {name: i for i, name in enumerate(REGLIST)}
    return {
        "AllAreas": np.ones(matrix.shape[0], dtype=bool),
        "PPC": matrix[:, index["PPC"]].copy(),
        "S1": matrix[:, [index[x] for x in ("LL", "UL", "Whisker", "Trunk")]].any(axis=1),
        "M1M2": matrix[:, [index["M1"], index["M2"]]].any(axis=1),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--last", type=int, default=19)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    csp = LoadDataFromMatlab.load_cell_SP_posture(
        dir_BrainRegions=r"/path/to/data/BrainRegions.mat",
        dir_AllMiceDataFile=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice.mat",
        dir_ListOfDirectoriesForAllMice=r"/path/to/data/AllMiceData/ListOfDirectoriesForAllMice2_for_python.mat",
    )
    total_start = time.perf_counter()
    for plane in range(args.first, args.last + 1):
        output_path = RESULTS / f"fig6b_m1m2_plane{plane:02d}.npz"
        if output_path.exists() and not args.force:
            print(f"skip existing {output_path}", flush=True)
            continue
        session = str(csp.PlaneNames[plane])
        plane_start = time.perf_counter()
        print(f"=== plane {plane}: {session} ===", flush=True)
        csp.load_data(plane)
        csp.preprocessData(
            cutoff_freq=10, SamplingRate_Hz=10, order=4, triError_threshold_in_mm=4
        )
        csp.preprocessData_continuous(
            cutoff_freq=10,
            SamplingRate_Hz=10,
            order=4,
            triError_threshold_in_mm=4,
            UseCachedData=True,
        )
        n_cells = csp.cellac.shape[1]
        active = mc.active_neurons(csp.skewness)
        registered = masks_for_plane(plane, n_cells)
        values = np.full((len(REGION_KEYS), 2), np.nan, dtype=float)
        counts = np.zeros(len(REGION_KEYS), dtype=int)
        for region_index, region in enumerate(REGION_KEYS):
            mask = registered[region] & active
            counts[region_index] = int(mask.sum())
            if counts[region_index] < 50:
                continue
            start = time.perf_counter()
            values[region_index, 0] = observed_lri(csp, mask)
            values[region_index, 1] = null_lri(csp, mask)
            print(
                f"  {region:8s} n={counts[region_index]:5d} "
                f"gcPCA={values[region_index, 0]:.6f} "
                f"standardPCA={values[region_index, 1]:.6f} "
                f"delta={values[region_index, 0] - values[region_index, 1]:+.6f} "
                f"[{time.perf_counter() - start:.0f}s]",
                flush=True,
            )
        np.savez_compressed(
            output_path,
            lri=values,
            n_active=counts,
            regions=np.asarray(REGION_KEYS),
            planename=np.asarray(session),
            plane=np.asarray(plane),
            active_threshold=np.asarray(mc.SKEW_THRESHOLD),
            w1=np.asarray(W1),
            w2=np.asarray(W2),
            ph_dim=np.asarray(DIM),
        )
        print(f"wrote {output_path} [{time.perf_counter() - plane_start:.0f}s]", flush=True)
    print(f"DONE in {(time.perf_counter() - total_start) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
