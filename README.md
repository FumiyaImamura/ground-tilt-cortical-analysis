# Ground tilt cortical analysis

This repository contains figure-generation code and compact source data for selected
quantitative panels in *A cortex-wide self-consistent manifold for a body–environment
reference frame*.

The included workflows generate panels from Figures 2–6 and Extended Data Figures
7–11. They cover example neural activity, encoding, ridge-regression contribution
maps, regional decoding,
frontal–parietal coordination, shared population components, manifold geometry,
persistent homology, and region-restricted analyses.

## Repository contents

- `source_data/`: compact numerical data used by the included figure workflows
- `workflows/`: figure-generation code organized by analysis family
- `scripts/reproduce_all.py`: entry point for generating the included panels
- `environment/`: Conda environment definitions
- `docs/FIGURE_CODE_MAP.md`: panel-to-code and source-data map
- `docs/METHODS_TO_CODE.md`: analysis-to-code crosswalk
- `reference_figures/`: manuscript figure exports for orientation
- `outputs/`: locally generated panels and result tables, ignored by Git

The six-degree-of-freedom platform-control software and engineering files are
maintained separately.

## Installation

Clone the repository and create the standard environment:

```powershell
git clone git@github.com:FumiyaImamura/ground-tilt-cortical-analysis.git
cd ground-tilt-cortical-analysis
conda env create -f environment/environment.yml
conda activate ground-tilt-cortical-analysis
python -m pip install -e .
```

Figures 3h and 4–6, together with Extended Data Figures 10k and 11c,d, use the
topology environment:

```powershell
conda env create -f environment/topology.yml
conda activate ground-tilt-cortical-topology
python -m pip install -e .
```

## Generate figures

Generate the standard panels:

```powershell
python scripts/reproduce_all.py --profile standard
```

Generate the manifold and topology panels:

```powershell
python scripts/reproduce_all.py --profile topology
```

List the available tasks or generate one selected figure group:

```powershell
python scripts/reproduce_all.py --list
python scripts/reproduce_all.py --task figure_3f_g_i
```

Generated files are written under `outputs/`. See
[`docs/FIGURE_CODE_MAP.md`](docs/FIGURE_CODE_MAP.md) for the code and source-data
directory associated with each included panel.

## Data

Compact figure source data are included in `source_data/`. Raw imaging movies,
behavior videos, and larger processed datasets are distributed separately as
described in the paper's Data Availability statement.

## Citation and license

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Repository code is
distributed under the [`MIT License`](LICENSE). Python dependencies retain their
respective licenses.
