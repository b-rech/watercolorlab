# Particulate processing

Particulate processing is available through
`Spectra.particulate.compute_absorption()`.

The workflow requires total and depigmented transmittance and reflectance
measurements, identified through the metadata returned by `read_files`.

```python
import watercolorlab as wcl

data = wcl.read_files(
    path="data/particulate",
    instrument="perkinelmer",
    format_string="{id}_{rmode}_{pigment}",
    sample_pattern="sample",
    blank_pattern="blank",
    depig_pattern="EXT",
)

spectra = wcl.Spectra(data)
spectra.particulate.compute_absorption(
    type_measurement="tr",
    trans_pattern="T",
    wl_range=(400, 750),
)
```

See the [particulate ATBD](../particulate_atbd.md) for the equations,
corrections, quality-control rules, and output products.
