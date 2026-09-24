# Data dictionary

The `source_data/` directory is organized by figure-analysis family:

- `neuron_density/`: registered neuron coordinates and density-map inputs for Fig. 2c
- `example_activity/`: neural, platform, and toe traces plus imaging assets for Fig. 2d
- `encoding_examples/`: observed and model-predicted traces for Fig. 2i
- `encoding_performance/`: encoding-performance values for Extended Data Fig. 7c
- `decoding/`: regional decoding, neuron-count, and shuffled-target data for Fig. 3a–e
  and Extended Data Fig. 10a–j
- `frontoparietal_coordination/`: frontal and parietal decoder outputs and summary
  distributions for Fig. 3f–i and Extended Data Fig. 10k
- `gcarp/`: component-by-session movement, tilt, and posture contributions for Fig. 3l
- `sex_indicator/`: animal-level values and statistical summaries for Extended Data
  Fig. 9
- `manifold_topology/`: response summaries, manifold trajectories, ellipse values,
  persistent-homology values, and session-level topology results for Figs. 4–6 and
  Extended Data Fig. 11c,d
- `regional_topology/`: regional trajectory and ringiness data for Fig. 6a,b

The individual renderers document the array names and table columns they read.
