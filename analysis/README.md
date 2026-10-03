# Analysis scripts

MATLAB and Python scripts for single-neuron direction selectivity, the encoding models,
regional and frontal–parietal decoding, GCARP, cPCA/gcPCA population trajectories,
persistent homology and the landscape ringness index (LRI). The MATLAB batch analyses
include the job scripts used on the TSUBAME 4.0 cluster. Figure panels are generated from
the compact source data by the code in [`../workflows`](../workflows).

## Contents

| Analysis | Figure | Scripts |
| --- | --- | --- |
| Direction selectivity during 8DT, and the comparison of 8DT and T-Rot responses | Fig. 2f,g, Fig. 4b | [`single_neuron/direction_selectivity.m`](single_neuron/direction_selectivity.m), [`single_neuron/direction_selectivity_stats.py`](single_neuron/direction_selectivity_stats.py) |
| Full encoding model, and models with only platform or only posture and movement predictors | Extended Data Fig. 7c | [`encoding/ridge_full_model.m`](encoding/ridge_full_model.m) |
| Tilt-predictor ablation, fitted together with the full model | Fig. 2l, Extended Data Fig. 8a,b | [`encoding/ridge_tilt_ablation.m`](encoding/ridge_tilt_ablation.m) |
| Forelimb, hindlimb and trunk predictor ablations | Fig. 2k, Extended Data Fig. 8a,b | [`encoding/ridge_body_part_ablation.m`](encoding/ridge_body_part_ablation.m) |
| Tail predictor ablation | Fig. 2k, Extended Data Fig. 8a,b | [`encoding/ridge_tail_ablation.m`](encoding/ridge_tail_ablation.m) |
| Eigenposture ablations | Fig. 2n, Extended Data Fig. 8a,b | [`encoding/ridge_eigenposture_ablation.m`](encoding/ridge_eigenposture_ablation.m) |
| Regional decoding of 8DT, T-Rot, Roll, Yaw and Vert | Fig. 3a–e, Extended Data Fig. 10a–e | [`decoding/decode_8dt.m`](decoding/decode_8dt.m), [`decoding/decode_trot.m`](decoding/decode_trot.m), [`decoding/decode_roll.m`](decoding/decode_roll.m), [`decoding/decode_yaw.m`](decoding/decode_yaw.m), [`decoding/decode_vert.m`](decoding/decode_vert.m) |
| Tilt-direction decoding from frontal and from parietal neurons | Fig. 3f–i, Extended Data Fig. 10k | [`frontoparietal/decode_frontal_parietal.m`](frontoparietal/decode_frontal_parietal.m), [`frontoparietal/collect_decoded_directions.m`](frontoparietal/collect_decoded_directions.m) |
| GCARP for movement, the first eigenposture and ground tilt | Fig. 3l | [`gcarp/gcarp_decode_movement_posture.m`](gcarp/gcarp_decode_movement_posture.m), [`gcarp/gcarp_decode_tilt.m`](gcarp/gcarp_decode_tilt.m), [`gcarp/gcarp_contribution_movement_posture.m`](gcarp/gcarp_contribution_movement_posture.m), [`gcarp/gcarp_contribution_tilt.m`](gcarp/gcarp_contribution_tilt.m) |
| Example single-neuron responses to 8DT and T-Rot | Fig. 4a | [`single_neuron/example_neurons_8dt_trot.m`](single_neuron/example_neurons_8dt_trot.m) |
| cPCA and gcPCA population trajectories | Fig. 4d–h, 4j–n | [`manifold/population_trajectories.py`](manifold/population_trajectories.py), [`manifold/trajectory_helpers.py`](manifold/trajectory_helpers.py) |
| Persistent homology and LRI | Fig. 5d,f | [`manifold/persistent_homology_lri.py`](manifold/persistent_homology_lri.py) |
| LRI within cortical subdivisions | Fig. 6b | [`manifold/regional_lri.py`](manifold/regional_lri.py), [`manifold/export_region_masks.m`](manifold/export_region_masks.m) |
| Shared functions | | [`manifold/manifold_core.py`](manifold/manifold_core.py), [`helpers`](helpers) |

Each MATLAB batch script has a job script with the same name (`.sh`) for the cluster.

Paths are placeholders (`/path/to/...`) and must be set before use. Animal identifiers
follow the anonymized identifiers of the data deposit (`Animal_01` to `Animal_10`). The
MATLAB scripts are organized in `%%` sections, some of which only load saved results or
plot them. Run the sections you need rather than the whole file.

## Software

MATLAB R2023a or later (the scripts use `name=value` arguments, `anynan` and the
`"omitmissing"` flag) with the Statistics and Machine Learning, Signal Processing, Image
Processing and Parallel Computing toolboxes. The decoding scripts and
`gcarp/gcarp_contribution_tilt.m` also call `wrapToPi` or `wrapTo2Pi`, which
are installed with the Mapping Toolbox. The CSP helper class also contains functions
that use the Control System Toolbox.

Python 3.11 with the packages in [`requirements.txt`](requirements.txt).

## MATLAB helper functions

Set `path2functions` to **one** of the two helper folders:

- `helpers/matlab_v1`, a snapshot of the helper functions from June 2025, for
  `encoding/ridge_full_model.m`, the `decoding` scripts and
  `single_neuron/example_neurons_8dt_trot.m`.
- `helpers/matlab_v2` for all other MATLAB scripts. They call helper functions that are
  not in `matlab_v1`.

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
decode_8dt
```

or submit the job script, which runs all 20 sessions as an array job on an SGE-type
scheduler:

```bash
cd analysis/decoding
qsub -g YOUR_GROUP decode_8dt.sh
```

Add `-t 1-1` to run only the first session. Update the log paths (`/path/to/logs`) and
resource requests in the job script for your system. Each batch script writes its results
to the folder set in `MasterDir` (`/path/to/results/<script name>`), which is where the
follow-up scripts read them. The jobs load MATLAB with `module load matlab`, and
the decoding, frontal–parietal and GCARP scripts start a 12-worker parallel pool
(`parpool('local', 12)`). Lower the number on a machine whose local cluster profile
allows fewer workers.

Each decoding job samples 200 neurons from each of eight cortical areas, ten times,
reduces them to 20 principal components within each of ten cross-validation folds, and
decodes with 100 bagged trees (Fig. 3a–e). Extended Data Fig. 10a–e repeat the decoding
with `numcell2sample` set to 100 or 300.

The frontal–parietal job (`decode_frontal_parietal.m`) decodes the tilt direction
of each session in eight classes, once from all M1 and M2 neurons and once from all S1
(trunk, whisker, forelimb and hindlimb) and PPC neurons. Within each of ten
cross-validation folds, the activity is reduced to 1,000 principal components and then
to 100 cPCA dimensions (cPCA weights 0.3, 0.5, 0.7 and 0.9) and decoded with 100 bagged
trees. `collect_decoded_directions.m` collects the decoded and true directions of
the decoders that use all 100 dimensions into `Reg_comparison.mat`, the file in
`../source_data/frontoparietal_coordination`. Fig. 3f–i and Extended Data Fig. 10k use
the cPCA weight 0.5.

The GCARP jobs (`gcarp_decode_movement_posture.m` for movement and the eigenpostures,
`gcarp_decode_tilt.m` for the tilt direction) sample 200
neurons from each of M1, M2, S1tr, S1b, S1fl, S1hl, PPC and RSC, five times per session.
Areas with fewer than 200 neurons in a session are left out. Within each of ten
cross-validation folds, each area is reduced to 100 principal components, GCCA is fitted
across the areas, and 100 bagged trees decode from the shared GCCA components, either all
of them (`reg_data_base`) or with one of components 1 to 10 removed (`reg_data`).
`gcarp_contribution_movement_posture.m` and `gcarp_contribution_tilt.m` compute, for each
session and component, the drop
in cross-validated R² (movement and the first eigenposture) or in the fraction of tilt
directions decoded within 45° (ground tilt) when the component is removed. They save
`tDXYZbp_magnitude_mean.mat`, `eigen_posture_1.mat` and `theta.mat`, which are
`movement.mat`, `eigen_posture_1.mat` and `ground_tilt.mat` in `../source_data/gcarp`.

`direction_selectivity.m` computes, for every neuron of every session, the
direction selectivity score, the preferred direction and the p-value of a one-way ANOVA
across the eight 8DT directions (`assess_direction_selectivity` in CSP), and saves them
with the activity skewness, the area masks and the mean T-Rot responses. Its later
sections compute the Fig. 2f direction counts, the Fig. 2g percentages and the Fig. 4b
correlations from the saved file. Direction-selective neurons have an activity skewness
above 0.5, a score above 0.1 and an ANOVA p-value below 0.05.

`ridge_tilt_ablation.m` is a function that fits, for each session, the full model
and the model without the 82 tilt predictors (tilt direction and tilt angle) through
the same code. Fig. 2l uses these paired scores. For example,
`ridge_tilt_ablation(1:20, true)` fits both models for all sessions.

## Python scripts

Run the Python scripts from `analysis/manifold` with `helpers/python` on the Python
path, or set `CFDIR` (`CF` in `regional_lri.py`) near the top of each script. Set
the manifest and atlas paths passed to `load_cell_SP_posture()` in each script. Results
are written to `analysis/results` and figures to `analysis/figures`. All three scripts
use the neurons with activity skewness > 0.5 (`manifold_core.active_neurons()`).

- `population_trajectories.py [plane]` computes Fig. 4d–h and 4j–n for one session
  (zero-based `plane`, default 17). Activity is projected onto the first 1,000
  principal components of the session's active neurons and then onto cPCA (w2 = 0) or
  gcPCA (w2 = 0.9) axes. Because the sign of each axis is arbitrary, a trajectory can
  come out mirrored along an axis. If a previously saved projection of the same session
  is available (`saved_file`), the displayed axes are oriented to match it.
- `persistent_homology_lri.py [first] [last]` computes the LRI and the maximum H1
  lifetime of the T-Rot, 8DT and Vert point clouds at w2 = 0, 0.5, 1 and 3 for each
  session (Fig. 5d,f). The per-session outputs have the format of the
  `plane??_manifold.npz` files in `../source_data/manifold_topology`.
- `regional_lri.py` computes the cross-validated LRI of gcPCA and of standard PCA
  for all areas, PPC, S1 and M1+M2 (Fig. 6b). Run `export_region_masks.m` first. The
  per-session outputs have the format of the `fig6b_m1m2_plane??.npz` files in
  `../source_data/regional_topology`.

`single_neuron/direction_selectivity_stats.py` reads the file saved by
`direction_selectivity.m` and computes the Fig. 2f,g tests (chi-square tests,
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
