# Figure-to-code map

This document maps the quantitative manuscript panels included in this repository to
their figure-generation code and compact source-data directories. Schematics,
histology, photographs, and other non-computational artwork are not listed.

## Main figures

| Panel | Analysis | Figure-generation code | Source data |
|---|---|---|---|
| Fig. 2c | Recorded-neuron locations and density | `workflows/neuron_density/render_figure_2c.py` | `source_data/neuron_density/` |
| Fig. 2d | Example neuronal, platform, and toe trajectories | `workflows/example_activity/render_figure_2d.py` | `source_data/example_activity/` |
| Fig. 2i | Example encoding-model fits | `workflows/encoding_examples/render_figure_2i.py` | `source_data/encoding_examples/` |
| Fig. 3a–e | Regional population decoding | `workflows/decoding/render_figure_3a_e.py` | `source_data/decoding/` |
| Fig. 3f,g,i | Frontal and parietal decoding-error structure | `workflows/frontoparietal_coordination/render_figure_3fgi.py` | `source_data/frontoparietal_coordination/` |
| Fig. 3h | Frontal–parietal error coordination | `workflows/frontoparietal_coordination/step5_paper_ready_figures.py` | `source_data/frontoparietal_coordination/` |
| Fig. 3l | GCARP component contributions | `workflows/gcarp/render_figure_3l.py` | `source_data/gcarp/` |
| Fig. 4a,b,d–h,j–q | Neuronal responses and population manifolds | `workflows/manifold_topology/render_figures_4_6.py` | `source_data/manifold_topology/` |
| Fig. 5b–d,f | Persistent homology and landscape ringiness | `workflows/manifold_topology/render_figures_4_6.py` | `source_data/manifold_topology/` |
| Fig. 6a–c | Region-restricted manifold analyses | `workflows/manifold_topology/render_figures_4_6.py` | `source_data/regional_topology/` and `source_data/manifold_topology/` |

## Extended Data figures

| Panel | Analysis | Figure-generation code | Source data |
|---|---|---|---|
| ED7c | Encoding-performance distribution | `workflows/encoding_performance/render_extended_data_figure_7c.py` | `source_data/encoding_performance/` |
| ED9 | Sex and calcium-indicator comparisons | `workflows/sex_indicator/render_extended_data_figure_9.py` | `source_data/sex_indicator/` |
| ED10a–e | Neuron-count sensitivity | `workflows/decoding/build_paper_ready_assets.py` | `source_data/decoding/` |
| ED10f–j | Shuffled-target controls | `workflows/decoding/recenter_empirical_null_figures.py` | `source_data/decoding/` |
| ED10k | Excess over independence | `workflows/frontoparietal_coordination/step5_paper_ready_figures.py` | `source_data/frontoparietal_coordination/` |
| ED11c,d | Embedding-dimension controls | `workflows/manifold_topology/render_extended_data_figure_11_cd.py` | `source_data/manifold_topology/` |

## Figure-generation commands

The standard profile generates Fig. 2c,d,i, Fig. 3a–g,i,l, ED7c, ED9, and ED10a–j:

```powershell
python scripts/reproduce_all.py --profile standard
```

The topology profile generates Fig. 3h, the computational panels in Figures 4–6,
ED10k, and ED11c,d:

```powershell
python scripts/reproduce_all.py --profile topology
```

Use `python scripts/reproduce_all.py --list` to display individual task names and
`--task TASK_NAME` to run a selected figure group.
