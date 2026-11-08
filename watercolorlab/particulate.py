'''
Module: particulate.py
Purpose: Particulate absorption retrieval from transmittance-reflectance data.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Python: 3.11+

Dependencies: numpy, pandas, matplotlib, seaborn
Public API: ParticulateProcessor.compute_absorption
'''

import numpy as np
import pandas as pd

from .plots import plot_part_absorption_raw


class ParticulateProcessor:
    """Compute particulate absorption from transmittance/reflectance spectra.

    The processor is attached to a :class:`~watercolorlab.spectra.Spectra`
    instance and stores intermediate and final results on that shared object.
    """

    def __init__(self, spectra):
        """Create a processor associated with ``spectra``.

        Parameters
        ----------
        spectra : Spectra
            Container holding raw spectra and the metadata needed to classify
            transmission, reflection, total, depigmented, sample, and blank
            measurements.
        """
        self.spectra = spectra

    def compute_absorption(
            self,
            vol_diameter,
            percentual=False,
            type_measurement='absorbance',
            wl_offset=800,
            use_tau=True,
            trans_pattern='T',
            wl_range=None,
            plot_raw=False,
            split_total_depig_blank=True,
            custom_blanks=None):
        """Compute total and depigmented particulate absorption.

        The method converts the input to decimal transmittance/reflectance,
        separates the four required sample measurement classes, applies blank
        and reference-wavelength corrections, optionally applies the tau
        correction, and converts corrected absorbance to particulate
        absorption using the supplied filter dimensions. Depigmented
        normalized transmittance and reflectance values slightly above one
        are clipped to one to accommodate measurement uncertainty around a
        blank-equivalent curve. The associated
        ``Spectra`` object must be built from ``read_files`` data containing
        ``id`` and ``rmode`` metadata, sample/blank classification columns,
        and ``is_total``/``is_depig`` classification from ``depig_pattern``.

        Parameters
        ----------
        vol_diameter : tuple or pandas.DataFrame
            Either a ``(volume, diameter)`` pair shared by every sample or an
            indexed table with ``volume`` and ``diameter`` columns. Values are
            interpreted in cubic metres and metres, respectively, and must be
            finite and strictly positive.
        percentual : bool, optional
            If ``True``, divide input transmittance and reflectance values by
            100 before processing. It is invalid for absorbance input.
        type_measurement : {'absorbance', 'tr'}, optional
            Representation of the raw spectral columns. ``'absorbance'``
            converts values with ``10 ** -absorbance``; ``'tr'`` treats them
            as transmittance/reflectance values already.
        wl_offset : int, optional
            Wavelength where total and depigmented curves are aligned. It
            must be present in the selected wavelength columns.
        use_tau : bool, optional
            If ``True``, apply the tau correction to transmittance optical
            depth. Wavelength 750 nm is required in this mode.
        split_total_depig_blank : bool, optional
            If ``True``, calculate separate total and depigmented blanks. If
            ``False``, combine the two pigmentation classes within each
            measurement mode and reuse that blank for both paths.
        wl_range : tuple, list, or None, optional
            Inclusive wavelength interval to process. Both bounds must be
            present in the available wavelength columns, finite, and ordered
            from lower to higher wavelength.
        plot_raw : bool, optional
            If ``True``, display raw and offset-corrected diagnostic plots.
        custom_blanks : pandas.DataFrame or None, optional
            Optional replacement blanks in the same representation specified
            by ``type_measurement``. The index must contain exactly one row
            labelled by ``trans_pattern`` and one row labelled by a reflectance pattern,
            and the columns must be numeric integer wavelengths, including
            labels that can be converted to integers. Wavelengths must remain
            unique after normalization. This option requires
            ``split_total_depig_blank=False``.

        Returns
        -------
        Spectra
            The associated object, with ``absorbance_total``,
            ``absorbance_depig``, ``part_abs_total``, ``part_abs_depig``, and
            ``part_abs_phyto`` result attributes.

        Raises
        ------
        TypeError
            If ``custom_blanks`` or ``vol_diameter`` has an unsupported type.
        ValueError
            If wavelengths, measurement pairings, blanks, filter dimensions,
            or other required inputs are invalid.
        """
        spectra = self.spectra

        # A new computation invalidates results from any previous particulate
        # processing run, including when the new run fails validation.
        for attribute in (
                'absorbance_total', 'absorbance_depig',
                'part_abs_total', 'part_abs_depig', 'part_abs_phyto'):
            setattr(spectra, attribute, None)

        if not np.isfinite(spectra.raw_spectra.to_numpy()).all():
            raise ValueError(
                'Raw spectral measurements contain NaN or infinite values.'
            )

        if percentual and type_measurement == 'absorbance':
            raise ValueError(
                'percentual=True is only valid when type_measurement="tr".'
            )

        # Normalize the input representation first. The correction equations
        # operate on decimal T/R values, regardless of how the instrument
        # exported the original data.
        if type_measurement == 'absorbance':
            data = np.power(10, -spectra.raw_spectra.copy())
        elif type_measurement == 'tr':
            data = spectra.raw_spectra.copy()
        else:
            raise ValueError('type_measurement must be either "absorbance" or "tr"')

        # Restrict all subsequent calculations to the requested inclusive
        # wavelength interval while preserving the original column order.
        if wl_range is not None:
            if not isinstance(wl_range, (tuple, list, np.ndarray)) \
                    or len(wl_range) != 2:
                raise ValueError(
                    'wl_range must contain exactly two wavelength values.'
                )
            wli, wlf = wl_range
            if any(
                    not isinstance(wavelength, (int, float, np.number))
                    or not np.isfinite(wavelength)
                    for wavelength in (wli, wlf)):
                raise ValueError(
                    'wl_range values must be finite numeric wavelengths.'
                )
            if wli > wlf:
                raise ValueError('wl_range must satisfy wli <= wlf.')
            if wli not in spectra.wls or wlf not in spectra.wls:
                raise ValueError('wl_range values out of bounds.')
            selected_wavelengths = [
                wavelength for wavelength in data.columns
                if wli <= wavelength <= wlf
            ]
            if not selected_wavelengths:
                raise ValueError(
                    'wl_range does not contain any available wavelengths.'
                )
            data = data.loc[:, selected_wavelengths]

        # Percentage input must be converted before blank subtraction and
        # division; otherwise every correction would be off by a factor of 100.
        if percentual:
            data /= 100

        # Keep the selected wavelengths available for validation of offsets,
        # custom blanks, and later tau calculations.
        wls = data.columns

        # Convert metadata into boolean masks used to build the four sample
        # tables and the corresponding blank reference curves.
        required_metadata = {
            'id', 'rmode', 'is_blank', 'is_sample', 'is_total'
        }
        missing_metadata = required_metadata.difference(
            spectra.raw_data.columns)
        if missing_metadata:
            raise ValueError(
                'Particulate processing requires metadata columns: '
                f'{sorted(required_metadata)}. Missing: '
                f'{sorted(missing_metadata)}.'
            )
        fblank = spectra.is_blank.to_numpy()
        fsample = spectra.is_sample.to_numpy()
        ftotal = spectra.is_total.to_numpy()
        sample_data = data.loc[fsample]
        if not np.isfinite(sample_data.to_numpy()).all():
            raise ValueError(
                'Sample measurements contain NaN or infinite values.'
            )
        relevant_rows = fsample | fblank
        rmode_values = pd.Index(
            spectra.rmode.loc[relevant_rows].dropna().unique())
        trans_label = str(trans_pattern).casefold()
        trans_modes = [
            mode for mode in rmode_values
            if str(mode).casefold() == trans_label
        ]
        if len(trans_modes) != 1:
            raise ValueError(
                f'trans_pattern {trans_pattern!r} must match exactly one '
                'mode in spectra.rmode, ignoring case.'
            )
        trans_mode = trans_modes[0]
        reflectance_modes = [
            mode for mode in rmode_values
            if str(mode).casefold() != trans_label
        ]
        if len(reflectance_modes) != 1:
            raise ValueError(
                'The sample and blank rows must contain exactly one '
                'reflectance pattern besides trans_pattern '
                f'{trans_pattern!r}; found {reflectance_modes!r}.'
            )
        reflect_pattern = reflectance_modes[0]
        ftrans = spectra.rmode.eq(trans_mode).to_numpy()
        frefl = spectra.rmode.eq(reflect_pattern).to_numpy()

        # Validate and normalize caller-supplied blanks before combining them
        # with dataset measurements. This keeps custom and automatic blank
        # paths in the same decimal T/R representation.
        custom_blank_data = None
        if custom_blanks is not None:
            if split_total_depig_blank:
                raise ValueError(
                    'split_total_depig_blank must be False when '
                    'custom_blanks is provided.')
            if not isinstance(custom_blanks, pd.DataFrame):
                raise TypeError('custom_blanks must be a pandas DataFrame.')
            if custom_blanks.empty:
                raise ValueError('custom_blanks must not be empty.')

            blank_modes = pd.Index(custom_blanks.index)
            normalized_modes = blank_modes.map(lambda mode: str(mode).casefold())
            reflect_label = str(reflect_pattern).casefold()
            if not normalized_modes.isin([trans_label, reflect_label]).all():
                raise ValueError(
                    'custom_blanks index must contain only transmittance and '
                    'reflectance pattern labels.')
            if ((normalized_modes == trans_label).sum() != 1
                    or (normalized_modes == reflect_label).sum() != 1):
                raise ValueError(
                    'custom_blanks must contain exactly one row for each '
                    'measurement mode.')

            custom_blank_data = custom_blanks.copy()
            try:
                numeric_wavelengths = np.asarray(
                    [float(wavelength) for wavelength in custom_blank_data.columns],
                    dtype=float)
                if not np.equal(numeric_wavelengths,
                                numeric_wavelengths.astype(int)).all():
                    raise ValueError
                custom_blank_data.columns = pd.Index(
                    numeric_wavelengths.astype(int), dtype='int64')
                if custom_blank_data.columns.has_duplicates:
                    raise ValueError
            except (TypeError, ValueError, OverflowError) as error:
                raise ValueError(
                    'custom_blanks columns must be unique integer '
                    'wavelengths.') from error

            missing_wavelengths = wls.difference(custom_blank_data.columns)
            if not missing_wavelengths.empty:
                raise ValueError(
                    'custom_blanks is missing wavelengths: '
                    f'{missing_wavelengths.tolist()}')
            custom_blank_data = custom_blank_data.loc[:, wls]
            custom_blank_data.index = pd.Index([
                trans_pattern if mode == trans_label else reflect_pattern
                for mode in normalized_modes
            ])
            custom_blank_data = custom_blank_data.apply(pd.to_numeric,
                                                        errors='coerce')
            if custom_blank_data.isna().any().any():
                raise ValueError(
                    'custom_blanks must contain only numeric measurements.')
            if type_measurement == 'absorbance':
                custom_blank_data = (-custom_blank_data).rpow(10)
            if percentual:
                custom_blank_data /= 100

        # Extract the four required sample classes. Their original file keys
        # are temporary; the tables are aligned by parsed sample ID below.
        trans_total_raw = data.loc[fsample & ftrans & ftotal].sort_index()
        refl_total_raw = data.loc[fsample & frefl & ftotal].sort_index()
        trans_depig_raw = data.loc[fsample & ftrans & ~ftotal].sort_index()
        refl_depig_raw = data.loc[fsample & frefl & ~ftotal].sort_index()

        # Pair the four measurements by their parsed sample ID, not by their
        # current row position. Sorting the source keys is not sufficient:
        # one missing or differently ordered file would otherwise make a
        # measurement inherit another sample's correction and label.
        raw_curves = {
            'transmittance total': trans_total_raw,
            'reflectance total': refl_total_raw,
            'transmittance depigmented': trans_depig_raw,
            'reflectance depigmented': refl_depig_raw,
        }
        aligned_curves = {}
        for curve_name, curves in raw_curves.items():
            sample_ids = spectra.id.reindex(curves.index)
            if sample_ids.isna().any():
                raise ValueError(
                    f'Could not determine sample IDs for {curve_name} curves.'
                )
            if sample_ids.duplicated().any():
                raise ValueError(
                    f'{curve_name.capitalize()} curves contain duplicate '
                    'measurements for a sample ID.'
                )
            aligned = curves.copy()
            aligned.index = pd.Index(sample_ids.to_numpy(), name='id')
            aligned_curves[curve_name] = aligned.sort_index()

        expected_ids = None
        for curve_name, curves in aligned_curves.items():
            curve_ids = curves.index
            if expected_ids is None:
                expected_ids = curve_ids
                continue
            if not curve_ids.equals(expected_ids):
                raise ValueError(
                    'Each sample must have exactly one total and one '
                    'depigmented measurement for both transmittance and '
                    f'reflectance; {curve_name} IDs do not match.'
                )

        (trans_total_raw, refl_total_raw, trans_depig_raw, refl_depig_raw) = (
            aligned_curves['transmittance total'],
            aligned_curves['reflectance total'],
            aligned_curves['transmittance depigmented'],
            aligned_curves['reflectance depigmented'],
        )

        # Calculate one reference curve for each measurement mode and
        # pigmentation class, either from custom blanks or from the dataset.
        if custom_blank_data is not None:
            Tblank_total = custom_blank_data.loc[trans_pattern]
            Rblank_total = custom_blank_data.loc[reflect_pattern]
            Tblank_depig = Tblank_total.copy()
            Rblank_depig = Rblank_total.copy()
            print('Using custom blanks for transmittance and reflectance.')
        elif split_total_depig_blank:
            Tblank_total = data.loc[fblank & ftrans & ftotal].mean()
            Rblank_total = data.loc[fblank & frefl & ftotal].mean()
            Tblank_depig = data.loc[fblank & ftrans & ~ftotal].mean()
            Rblank_depig = data.loc[fblank & frefl & ~ftotal].mean()
            print('Using separate total and depigmented blanks from the dataset.')
        else:
            Tblank_total = data.loc[fblank & ftrans].mean()
            Rblank_total = data.loc[fblank & frefl].mean()
            Tblank_depig = Tblank_total.copy()
            Rblank_depig = Rblank_total.copy()
            print('Using combined total and depigmented blanks from the dataset.')

        # If only one pigmentation class has a valid blank, reuse it for the
        # other class and make that fallback explicit to the caller.
        if Tblank_total.isna().all() and not Tblank_depig.isna().all():
            Tblank_total = Tblank_depig.copy()
            print('Warning: Transmittance blank total is missing; '
                  'using depigmented blank instead.')
        elif Tblank_depig.isna().all() and not Tblank_total.isna().all():
            Tblank_depig = Tblank_total.copy()
            print('Warning: Transmittance blank depigmented is missing; '
                  'using total blank instead.')

        if Rblank_total.isna().all() and not Rblank_depig.isna().all():
            Rblank_total = Rblank_depig.copy()
            print('Warning: Reflectance blank total is missing; '
                  'using depigmented blank instead.')
        elif Rblank_depig.isna().all() and not Rblank_total.isna().all():
            Rblank_depig = Rblank_total.copy()
            print('Warning: Reflectance blank depigmented is missing; '
                  'using total blank instead.')

        # An empty row selection still produces a non-empty Series after
        # pandas calculates its mean; every value in that Series is NaN.
        # Reject references containing any non-finite value before arithmetic can
        # propagate missing values into the particulate calculation.
        blank_curves = {
            'transmittance total': Tblank_total,
            'transmittance depigmented': Tblank_depig,
            'reflectance total': Rblank_total,
            'reflectance depigmented': Rblank_depig,
        }
        invalid_blanks = [
            name for name, curve in blank_curves.items()
            if curve.empty or not np.isfinite(curve.to_numpy()).all()
        ]
        if invalid_blanks:
            raise ValueError(
                'Missing or invalid blank measurements for: '
                f'{", ".join(invalid_blanks)}.'
            )

        # The total and depigmented curves must share this wavelength because
        # the next step anchors each depigmented curve to its total partner.
        if wl_offset not in wls:
            raise ValueError('wl_offset out of bounds.')

        # Shift each total transmittance curve downward when it exceeds its
        # blank. A single per-sample offset is chosen from the largest excess.
        trans_blank_offset = (trans_total_raw.sub(Tblank_total, axis=1)
                              .max(axis=1)
                              .clip(lower=0))
        trans_total = trans_total_raw.sub(trans_blank_offset, axis=0)

        # Apply the same blank-boundary correction independently to total
        # reflectance curves.
        refl_blank_offset = (refl_total_raw.sub(Rblank_total, axis=1)
                             .max(axis=1)
                             .clip(lower=0))
        refl_total = refl_total_raw.sub(refl_blank_offset, axis=0)

        # Shift each depigmented curve so it equals its total counterpart at
        # the reference wavelength while preserving its spectral shape.
        trans_offset = trans_depig_raw[wl_offset] - trans_total[wl_offset].to_numpy()
        trans_depig = trans_depig_raw.sub(trans_offset, axis=0)

        refl_offset = refl_depig_raw[wl_offset] - refl_total[wl_offset].to_numpy()
        refl_depig = refl_depig_raw.sub(refl_offset, axis=0)

        # Diagnostic plotting receives already ID-aligned tables so displayed
        # curves cannot be labeled with another sample's metadata.
        if plot_raw:
            plot_ids = pd.Series(
                trans_total_raw.index,
                index=trans_total_raw.index,
                name='id')
            plot_part_absorption_raw(
                trans_total_raw, trans_depig_raw,
                refl_total_raw, refl_depig_raw,
                trans_total, trans_depig, refl_total, refl_depig,
                Tblank_total, Tblank_depig, Rblank_total, Rblank_depig,
                plot_ids, trans_mode, reflect_pattern)

        # Normalize every corrected sample curve by its matching blank curve.
        # All four tables have identical sample-ID indexes at this point.
        Tt = trans_total.div(Tblank_total, axis=1)
        Td = trans_depig.div(Tblank_depig, axis=1)
        Rt = refl_total.div(Rblank_total, axis=1)
        Rd = refl_depig.div(Rblank_depig, axis=1)

        normalized_curves = {
            'total transmittance': Tt,
            'depigmented transmittance': Td,
            'total reflectance': Rt,
            'depigmented reflectance': Rd,
        }
        depigmented_above_one = (Td > 1).sum().sum() + (Rd > 1).sum().sum()
        if depigmented_above_one:
            print(
                'Warning: '
                f'{depigmented_above_one} depigmented normalized '
                'transmittance/reflectance values exceeded one and were '
                'set to one.'
            )
            Td = Td.clip(upper=1)
            Rd = Rd.clip(upper=1)
            normalized_curves['depigmented transmittance'] = Td
            normalized_curves['depigmented reflectance'] = Rd

        for curve_name, curve in normalized_curves.items():
            if not np.isfinite(curve.to_numpy()).all():
                raise ValueError(
                    f'Non-finite values were calculated for {curve_name}.'
                )
            if (curve <= 0).any().any() or (curve > 1).any().any():
                raise ValueError(
                    f'{curve_name.capitalize()} values must be greater than '
                    'zero and no greater than one.'
                )

        # Convert transmittance to optical depth and apply the empirical tau
        # adjustment when requested. Values outside the calibrated range use
        # tau=1 rather than the fitted expression.
        if use_tau:
            if 750 not in data.columns:
                raise ValueError('Wavelength 750 is required when use_tau=True')
            odt_total = np.log10(1 / Tt)
            odt_depig = np.log10(1 / Td)
            odts_total = odt_total.sub(.5 * odt_total[750], axis=0)
            odts_depig = odt_depig.sub(.5 * odt_depig[750], axis=0)
            tau_total = 1.15 - 0.17 * odts_total
            tau_total[(odts_total <= .02) | (odts_total >= .7)] = 1
            tau_depig = 1.15 - 0.17 * odts_depig
            tau_depig[(odts_depig <= .02) | (odts_depig >= .7)] = 1
        else:
            tau_total = tau_depig = 1

        # Apply the particulate T/R relationship separately to total and
        # depigmented measurements.
        absorbance_total = ((1 - Tt + Rblank_total * (Tt - Rt)) /
                            (1 + Rblank_total * Tt * tau_total))
        absorbance_depig = ((1 - Td + Rblank_depig * (Td - Rd)) /
                            (1 + Rblank_depig * Td * tau_depig))
        if (not np.isfinite(absorbance_total.to_numpy()).all()
                or not np.isfinite(absorbance_depig.to_numpy()).all()):
            raise ValueError(
                'Non-finite particulate absorbance was calculated from the '
                'provided measurements.'
            )
        negative_total = (absorbance_total < 0).any().any()
        negative_depig = (absorbance_depig < 0).any().any()
        if negative_total or negative_depig:
            clipped_curves = []
            if negative_total:
                clipped_curves.append('total')
            if negative_depig:
                clipped_curves.append('depigmented')
            print(
                'Warning: negative particulate absorbance values were '
                f'clipped to zero for {", ".join(clipped_curves)} curves.'
            )
        absorbance_total = absorbance_total.clip(lower=0)
        absorbance_depig = absorbance_depig.clip(lower=0)
        if ((absorbance_total >= 1).any().any()
                or (absorbance_depig >= 1).any().any()):
            raise ValueError(
                'Calculated particulate absorbance must be less than one '
                'before optical-density conversion.'
            )

        spectra.absorbance_total = absorbance_total
        spectra.absorbance_depig = absorbance_depig

        # Convert absorbance to optical density before applying the
        # filter-area/filtered-volume normalization.
        od_total = np.log10(1 / (1 - absorbance_total))
        od_depig = np.log10(1 / (1 - absorbance_depig))

        # Build the geometry normalization. A tuple supplies one factor for
        # every sample; a DataFrame supplies a factor indexed by sample ID.
        if isinstance(vol_diameter, tuple) and len(vol_diameter) == 2:
            volume, diameter = vol_diameter
            # Coerce tuple values through pandas so strings and non-finite
            # values fail with a domain-specific error instead of arithmetic
            # producing NaN or a low-level TypeError.
            dimensions = pd.to_numeric(
                pd.Series([volume, diameter]),
                errors='coerce')
            if dimensions.isna().any() or not np.isfinite(
                    dimensions.to_numpy(dtype=float)).all():
                raise ValueError(
                    'volume and diameter must be finite numeric values.')
            volume, diameter = dimensions
            if volume <= 0 or diameter <= 0:
                raise ValueError(
                    'volume and diameter must be greater than zero.')
            normalization = np.pi * diameter ** 2 / (4 * volume)
        elif isinstance(vol_diameter, pd.DataFrame):
            required = {'volume', 'diameter'}
            if not required.issubset(vol_diameter.columns):
                raise ValueError(
                    'vol_diameter must contain volume and diameter columns.')
            if vol_diameter.index.has_duplicates:
                raise ValueError(
                    'vol_diameter must have unique sample IDs in its index.')
            # Reindex by the computed sample IDs so dimensions cannot be
            # accidentally paired by row position.
            dimensions = vol_diameter.reindex(od_total.index)
            if dimensions[['volume', 'diameter']].isna().any().any():
                raise ValueError(
                    'vol_diameter must have volume and diameter values for '
                    'every sample id.')
            volume = pd.to_numeric(dimensions['volume'], errors='coerce')
            diameter = pd.to_numeric(dimensions['diameter'], errors='coerce')
            if (volume.isna().any() or diameter.isna().any()
                    or not np.isfinite(volume.to_numpy(dtype=float)).all()
                    or not np.isfinite(diameter.to_numpy(dtype=float)).all()):
                raise ValueError(
                    'volume and diameter must be finite numeric values.')
            if (volume <= 0).any() or (diameter <= 0).any():
                raise ValueError(
                    'volume and diameter must be greater than zero.')
            normalization = np.pi * diameter ** 2 / (4 * volume)
        else:
            raise TypeError(
                'vol_diameter must be a (volume, diameter) tuple or a '
                'DataFrame with volume and diameter columns.')

        # Apply the empirical particulate absorption relationship and the
        # geometry factor to obtain the final total and depigmented results.
        spectra.part_abs_total = np.log(10) * .719 * (od_total ** 1.2287)
        spectra.part_abs_depig = np.log(10) * .719 * (od_depig ** 1.2287)
        spectra.part_abs_total = spectra.part_abs_total.mul(normalization, axis=0)
        spectra.part_abs_depig = spectra.part_abs_depig.mul(normalization, axis=0)
        spectra.part_abs_phyto = spectra.part_abs_total - spectra.part_abs_depig
        return spectra
