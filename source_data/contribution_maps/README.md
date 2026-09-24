# Ridge-regression contribution maps

`ridge_regression_maps.npz` contains the full-cohort numerical maps used for
Figure 2k,l,n and Extended Data Figure 8a,b.

The archive contains three map stacks in the order given by `map_names`:

- `relative_maps`: local mean relative contribution among model-encoding neurons
- `normalized_maps`: density of high-contribution neurons relative to model-encoding neurons
- `density_maps`: density of neurons that are both model-encoding and high-contribution relative to all neurons

Each stack has a corresponding display mask and upper color limit. The atlas image,
atlas alpha mask, and manuscript colormaps are included so the panel renderer does
not depend on external files.

Run:

```powershell
python workflows/contribution_maps/render_ridge_regression_maps.py
```

The generated SVG and PNG panels are written to
`outputs/contribution_maps/panels/`.
