'''
Module: spectra.py
Purpose: Shared spectral data container and processor composition.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Python: 3.11+

Dependencies: pandas, watercolorlab.cdom, watercolorlab.consistency, watercolorlab.correction,
watercolorlab.particulate

Public API: Spectra, spectra.cdom.compute_absorption(),
spectra.particulate.compute_absorption(), spectra.consistency.check()
'''

# %% Dependencies

import pandas as pd

from .cdom import CDOMProcessor
from .consistency import ConsistencyProcessor
from .correction import BreakMethods
from .particulate import ParticulateProcessor


# %% Class definition

class Spectra:
    """Class to handle laboratory spectral measurements.

    Processing attributes are initialized to ``None`` and populated by the
    corresponding processor methods. Declaring them here keeps the public
    result API visible to static type checkers.
    """

    raw_samples: pd.DataFrame | None
    absorbance: pd.DataFrame | None
    absorption: pd.DataFrame | None
    cdom_fitted: pd.DataFrame | None
    blank_curves: pd.DataFrame | None
    blank_salinity: pd.DataFrame | None
    npc_offset: pd.Series | None
    suspicious: list | None
    absorbance_total: pd.DataFrame | None
    absorbance_depig: pd.DataFrame | None
    part_abs_total: pd.DataFrame | None
    part_abs_depig: pd.DataFrame | None
    part_abs_phyto: pd.DataFrame | None

    def __init__(self, data: 'pd.DataFrame'):
        """Create a Spectra container and attach its processing modules.

        Parameters
        ----------
        data : pandas.DataFrame
            Table containing an ``id`` column, metadata, and wavelength
            columns whose names contain numeric values.
        """
        # Keep an independent table and normalize wavelength column names.
        self.raw_data = data.copy()
        wavelength_columns = self.raw_data.filter(regex=r'\d').columns
        wavelength_values = pd.Index(
            wavelength_columns.map(float).map(int),
            dtype='int64',
        )
        self.raw_data.rename(
            columns=dict(zip(wavelength_columns, wavelength_values)),
            inplace=True,
        )
        self.raw_spectra = self.raw_data[wavelength_values]
        self.raw_spectra.columns = wavelength_values
        self.wls = wavelength_values.to_numpy()

        # Processing outputs are created lazily by the domain-specific
        # processors, but are declared here so their public API is explicit.
        self.raw_samples = None
        self.absorbance = None
        self.absorption = None
        self.cdom_fitted = None
        self.blank_curves = None
        self.blank_salinity = None
        self.npc_offset = None
        self.suspicious = None
        self.absorbance_total = None
        self.absorbance_depig = None
        self.part_abs_total = None
        self.part_abs_depig = None
        self.part_abs_phyto = None

        # Expose metadata columns as attributes used by the processors.
        meta_cols = self.raw_data.columns.difference(self.raw_spectra.columns)
        for col in meta_cols:
            setattr(self, col, self.raw_data[col])

        # Domain-specific processors operate on this shared object.
        self.cdom = CDOMProcessor(self)
        self.consistency = ConsistencyProcessor(self)
        self.correction = BreakMethods(self)
        self.particulate = ParticulateProcessor(self)

    def __repr__(self):
        """Return a compact representation with observation dimensions."""
        return (f'Spectra(observations={len(self.raw_spectra)}, '
                f'wavelengths={len(self.wls)})')

    def __str__(self):
        """Return a human-readable description of the object."""
        return (f'Object Spectra with {len(self.raw_spectra)} observations and '
                f'{len(self.wls)} wavelengths')
