# Analysis scripts

MATLAB and Python scripts for single-neuron direction selectivity, the encoding models,
regional and frontal–parietal decoding, GCARP, cPCA/gcPCA population trajectories,
persistent homology and the landscape ringness index (LRI). The MATLAB batch analyses
include the job scripts used on the TSUBAME 4.0 cluster. Figure panels are generated from
the compact source data by the code in [`../workflows`](../workflows).

## Contents

| Analysis | Figure | Source files |
| --- | --- | --- |
| Direction selectivity during 8DT, and the comparison of 8DT and T-Rot responses | Fig. 2f,g, Fig. 4b | [`single_neuron/DirectionSelectivity_20260327.m`](single_neuron/DirectionSelectivity_20260327.m), [`single_neuron/direction_selectivity_stats.py`](single_neuron/direction_selectivity_stats.py) |
| Full encoding model, and models with only platform or only posture and movement predictors | Extended Data Fig. 7c | [`encoding/RidgeRegression_continuous_20250621`](encoding/RidgeRegression_continuous_20250621) |
| Tilt-predictor ablation, fitted together with the full model | Fig. 2l, Extended Data Fig. 8a,b | [`encoding/tilt_ablation/fit_theta_phi_ablation.m`](encoding/tilt_ablation/fit_theta_phi_ablation.m) |
| Forelimb, hindlimb and trunk predictor ablations | Fig. 2k, Extended Data Fig. 8a,b | [`encoding/Ridge_reductionChunk_20250705`](encoding/Ridge_reductionChunk_20250705) |
| Tail predictor ablation | Fig. 2k, Extended Data Fig. 8a,b | [`encoding/Ridge_reductionTail_20250708`](encoding/Ridge_reductionTail_20250708) |
| Eigenposture ablations | Fig. 2n, Extended Data Fig. 8a,b | [`encoding/Ridge_cont_eigpos_20250628`](encoding/Ridge_cont_eigpos_20250628) |
| Regional decoding of 8DT, T-Rot, Roll, Yaw and Vert | Fig. 3a–e, Extended Data Fig. 10a–e | [`decoding`](decoding), one MATLAB script and one job script per condition |
| Tilt-direction decoding from frontal and from parietal neurons | Fig. 3f–i, Extended Data Fig. 10k | [`frontoparietal/TB_cPCA_20250706`](frontoparietal/TB_cPCA_20250706), [`frontoparietal/TB_classi_cPCA_20250706_analysis.m`](frontoparietal/TB_classi_cPCA_20250706_analysis.m) |
| GCARP for movement, the first eigenposture and ground tilt | Fig. 3l | [`gcarp/TBregress_multi_allreg_anchor_20250617_CSP2`](gcarp/TBregress_multi_allreg_anchor_20250617_CSP2), [`gcarp/TB_regress_multireg_anchor_20250617_analysis.m`](gcarp/TB_regress_multireg_anchor_20250617_analysis.m), [`gcarp/TB_classi_anchor_20250712`](gcarp/TB_classi_anchor_20250712), [`gcarp/TB_classi_anchor_20250712_analysis.m`](gcarp/TB_classi_anchor_20250712_analysis.m) |
| Example single-neuron responses to 8DT and T-Rot | Fig. 4a | [`single_neuron/DirectionSelectivity_p8_and_rot_20250621.m`](single_neuron/DirectionSelectivity_p8_and_rot_20250621.m) |
| cPCA and gcPCA population trajectories | Fig. 4d–h, 4j–n | [`manifold/fig4dn_active_paper_style.py`](manifold/fig4dn_active_paper_style.py), [`manifold/trajectory_helpers.py`](manifold/trajectory_helpers.py) |
| Persistent homology and LRI | Fig. 5d,f | [`manifold/run_manifold_analysis.py`](manifold/run_manifold_analysis.py) |
| LRI within cortical subdivisions | Fig. 6b | [`manifold/fig6b_regional_lri.py`](manifold/fig6b_regional_lri.py), [`manifold/export_region_masks.m`](manifold/export_region_masks.m) |
| Shared functions | | [`manifold/manifold_core.py`](manifold/manifold_core.py), [`helpers`](helpers) |

Paths are placeholders (`/path/to/...`) and must be set before use. Animal identifiers
follow the anonymized identifiers of the data deposit (`Animal_01` to `Animal_10`). The
MATLAB scripts are organized in `%%` sections, some of which only load saved results or
plot them. Run the sections you need rather than the whole file.

## Software

MATLAB R2023a or later (the scripts use `name=value` arguments, `anynan` and the
`"omitmissing"` flag) with the Statistics and Machine Learning, Signal Processing, Image
Processing and Parallel Computing toolboxes. The decoding scripts and
`gcarp/TB_classi_anchor_20250712_analysis.m` also call `wrapToPi` or `wrapTo2Pi`, which
are installed with the Mapping Toolbox. The CSP helper class also contains functions
that use the Control System Toolbox.

Python 3.11 with the packages in [`requirements.txt`](requirements.txt).

## MATLAB helper functions

Set `path2functions` to **one** of the two helper folders:

- `helpers/matlab_20250622`, a snapshot of the helper functions from June 2025, for
  `encoding/RidgeRegression_continuous_20250621`, the `decoding` scripts and
  `single_neuron/DirectionSelectivity_p8_and_rot_20250621.m`.
- `helpers/matlab_current` for all other MATLAB scripts. They call helper functions
  that are not in the June 2025 snapshot.

Do not add both folders to the MATLAB path at once. Also set the default data and atlas
paths near the top of `Matlab_classdef/CSP.m` in the chosen folder, because some
sections construct `CSP()` without path arguments.

## Input files

The scripts start from the per-session files of the original acquisition and
preprocessing pipeline, listed in a session manifest. These files are not included in
this repository. The table lists the variables the scripts read from them.

| Input | Fields read |
| --- | --- |
| Neural activity | `cellac`, time by neuron |
| Platform and synchronization | `org`, including the acquisition-to-time mappings, platform parameters, trial boundaries and motion labels |
| Tracked posture | `Q_bp` in `DLCdata_3DpointsOnly.mat` |
| Posture metadata | `rec3d` in `DLCdata_info.mat` |
| Registered neuron coordinates | `reg.roi` |
| Activity skewness | `skew`. Neurons with skewness > 0.5 are analysed. |
| Cortical atlas | `BrainRegions.mat` or `BrainRegions_with_RSC.mat`, with region names, contours and atlas image |

The MATLAB manifest contains `newallmicedir` and `allmiceNames`. `newallmicedir` is a
cell array of per-session structs whose fields `cellac`, `SyncAndSP`,
`DLCdata_3DpointsOnly`, `DLCdata_info`, `reg` and `skewness` hold the corresponding file
paths. CSP takes the text before the first underscore of a session name as the animal
and reads `left` or `right` from the rest, so name the sessions, for example,
`Animal01_left` and `Animal01_right`. The Python loader
(`helpers/python/LoadDataFromMatlab.py`) reads the same information from the
`char_allmicedir` and `char_allmiceNames` fields of a SciPy-readable manifest.
`PlaneNames` in the loader lists the sessions in manifest order.

The recordings deposited on Figshare are the same sessions after synchronization and
preprocessing, in a different layout (`Neural_activity` at 10 Hz, `Platform_param`,
`Posture` and trial annotations), restricted to neurons with activity skewness > 0.5.
The scripts here do not read that layout directly.

## Running the MATLAB scripts

Each batch script analyses one session, selected by `SGE_TASK_ID` (1 to 20). After
setting the paths, run one session locally from the script's folder:

```matlab
setenv('SGE_TASK_ID','1');
TB_regress_p8_20250719
```

or submit the job script, which runs all 20 sessions as an array job on an SGE-type
scheduler:

```bash
cd analysis/decoding/p8
qsub -g YOUR_GROUP TB_p8.sh
```

Add `-t 1-1` to run only the first session. Update the log paths and resource requests
in the job script for your system. The jobs load MATLAB with `module load matlab`, and
the decoding, frontal–parietal and GCARP scripts start a 12-worker parallel pool
(`parpool('local', 12)`). Lower the number on a machine whose local cluster profile
allows fewer workers.

Each decoding job samples 200 neurons from each of eight cortical areas, ten times,
reduces them to 20 principal components within each of ten cross-validation folds, and
decodes with 100 bagged trees (Fig. 3a–e). Extended Data Fig. 10a–e repeat the decoding
with `numcell2sample` set to 100 or 300.

The frontal–parietal job (`frontoparietal/TB_cPCA_20250706`) decodes the tilt direction
of each session in eight classes, once from all M1 and M2 neurons and once from all S1
(trunk, whisker, forelimb and hindlimb) and PPC neurons. Within each of ten
cross-validation folds, the activity is reduced to 1,000 principal components and then
to 100 cPCA dimensions (cPCA weights 0.3, 0.5, 0.7 and 0.9) and decoded with 100 bagged
trees. `TB_classi_cPCA_20250706_analysis.m` collects the decoded and true directions of
the decoders that use all 100 dimensions into `Reg_comparison.mat`, the file in
`../source_data/frontoparietal_coordination`. Fig. 3f–i and Extended Data Fig. 10k use
the cPCA weight 0.5.

The GCARP jobs (`gcarp/TBregress_multi_allreg_anchor_20250617_CSP2` for movement and the
eigenpostures, `gcarp/TB_classi_anchor_20250712` for the tilt direction) sample 200
neurons from each of M1, M2, S1tr, S1b, S1fl, S1hl, PPC and RSC, five times per session.
Areas with fewer than 200 neurons in a session are left out. Within each of ten
cross-validation folds, each area is reduced to 100 principal components, GCCA is fitted
across the areas, and 100 bagged trees decode from the shared GCCA components, either all
of them (`reg_data_base`) or with one of components 1 to 10 removed (`reg_data`).
`TB_regress_multireg_anchor_20250617_analysis.m` and
`TB_classi_anchor_20250712_analysis.m` compute, for each session and component, the drop
in cross-validated R² (movement and the first eigenposture) or in the fraction of tilt
directions decoded within 45° (ground tilt) when the component is removed. They save
`tDXYZbp_magnitude_mean.mat`, `eigen_posture_1.mat` and `theta.mat`, which are
`movement.mat`, `eigen_posture_1.mat` and `ground_tilt.mat` in `../source_data/gcarp`.

`DirectionSelectivity_20260327.m` computes, for every neuron of every session, the
direction selectivity score, the preferred direction and the p-value of a one-way ANOVA
across the eight 8DT directions (`assess_direction_selectivity` in CSP), and saves them
with the activity skewness, the area masks and the mean T-Rot responses. Its later
sections compute the Fig. 2f direction counts, the Fig. 2g percentages and the Fig. 4b
correlations from the saved file. Direction-selective neurons have an activity skewness
above 0.5, a score above 0.1 and an ANOVA p-value below 0.05.

`fit_theta_phi_ablation.m` is a function that fits, for each session, the full model
and the model without the 82 tilt predictors (tilt direction and tilt angle) through
the same code. Fig. 2l uses these paired scores. For example,
`fit_theta_phi_ablation(1:20, true)` fits both models for all sessions.

## Python scripts

Run the Python scripts from `analysis/manifold` with `helpers/python` on the Python
path, or set `CFDIR` (`CF` in `fig6b_regional_lri.py`) near the top of each script. Set
the manifest and atlas paths passed to `load_cell_SP_posture()` in each script. Results
are written to `analysis/results` and figures to `analysis/figures`. All three scripts
use the neurons with activity skewness > 0.5 (`manifold_core.active_neurons()`).

- `fig4dn_active_paper_style.py [plane]` computes Fig. 4d–h and 4j–n for one session
  (zero-based `plane`, default 17). Activity is projected onto the first 1,000
  principal components of the session's active neurons and then onto cPCA (w2 = 0) or
  gcPCA (w2 = 0.9) axes. Because the sign of each axis is arbitrary, a trajectory can
  come out mirrored along an axis. If a previously saved projection of the same session
  is available (`saved_file`), the displayed axes are oriented to match it.
- `run_manifold_analysis.py [first] [last]` computes the LRI and the maximum H1
  lifetime of the T-Rot, 8DT and Vert point clouds at w2 = 0, 0.5, 1 and 3 for each
  session (Fig. 5d,f). The per-session outputs have the format of the
  `plane??_manifold.npz` files in `../source_data/manifold_topology`.
- `fig6b_regional_lri.py` computes the cross-validated LRI of gcPCA and of standard PCA
  for all areas, PPC, S1 and M1+M2 (Fig. 6b). Run `export_region_masks.m` first. The
  per-session outputs have the format of the `fig6b_m1m2_plane??.npz` files in
  `../source_data/regional_topology`.

`single_neuron/direction_selectivity_stats.py` reads the file saved by
`DirectionSelectivity_20260327.m` and computes the Fig. 2f,g tests (chi-square tests,
then two-sided exact binomial tests per direction and two-sided Fisher's exact tests per
area, Bonferroni-corrected). It writes two tables to `analysis/results`:

```bash
python single_neuron/direction_selectivity_stats.py path/to/direction_selectivity_<time>.mat
```

The LRI of an H1 persistence diagram can also be computed directly:

```python
from ripser import ripser
import manifold_core as mc

h1 = ripser(points)['dgms'][1]   # points: samples by dimensions
value = mc.lri(h1)               # k_max = 10, as in the paper
```

## Optional exporters

Some MATLAB plotting sections save figures through `ImageSaver`, which calls the
external `plot2svg` exporter. `plot2svg` is not bundled and is not needed for model
fitting.
