# CDOM processing

CDOM processing is available through
`Spectra.cdom.compute_absorption()` and
`Spectra.cdom.expfit()`.

## Workflow

1. Read the measurement files with `read_files`.
2. Create a `Spectra` object.
3. Compute blank-corrected CDOM absorption.
4. Fit the exponential CDOM model when reference absorption and slope are
   required.

```python
import watercolorlab as wcl

data = wcl.read_files(
    path="data/cdom",
    instrument="wpi",
    format_string="{pattern}_{id}_{replicate}_{campaign}",
    sample_pattern="ponto",
    blank_pattern="milliq",
)

spectra = wcl.Spectra(data)
spectra.cdom.compute_absorption(
    wl_null_point_correction=(590, 600),
    pathlength=0.1,
)
spectra.cdom.expfit(wl_ref=443, wl_range=(350, 500))
```

The resulting tables are available as `spectra.absorption` and
`spectra.cdom_fitted`. For the complete mathematical description, see the
[CDOM ATBD](../cdom_atbd.md).
