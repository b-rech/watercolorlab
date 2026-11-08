# API reference

## Main classes and functions

- `watercolorlab.Spectra` — shared data container and processor composition.
- `watercolorlab.read_files` — reads instrument files and builds the metadata
  table consumed by `Spectra`.
- `Spectra.cdom.compute_absorption` — computes CDOM absorption coefficients.
- `Spectra.cdom.expfit` — fits the exponential CDOM model.
- `Spectra.particulate.compute_absorption` — computes particulate absorption
  products.

The repository also contains a legacy generated API page for
[`Spectra`](watercolorlab/spectra.html).
