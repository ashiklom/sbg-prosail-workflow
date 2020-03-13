# PROSAIL simulation workflow for SBG Uncertainty tasks

Currently, the main workhorse here is the `test-workflow.R` script.
This will do the following:

1. Generate a spatially-autocorrelated random LAI field.

2. Generate a purely random grid of PROSPECT and leaf angle distribution parameters.

3. Combine these two to generate a grid of PROSAIL-simulated hyperspectral reflectance at each grid cell.

4. Write the output to ENVI-format binary files (BSQ) that should be usable as input to Isofit or Hypertrace.
This can be read in directly by the Python `spectral` package -- the `read-bsq.py` script provides a barebones example.

5. As a simple example analysis, calculate the bi-directional reflectance NDVI and regress it against the true LAI.

## Future plans

- [ ] Combine with Hypertrace, Isofit, or an atmospheric RTM (e.g. LibRadTran, 6S, MODTRAN) to generate TOA reflectance
- [ ] Replace random inputs with vegetation model-driven (e.g. CLM, ED2) or observed (e.g. FIA) inputs
