"""
Module: cdom.py
Purpose: CDOM absorption correction and exponential curve fitting.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Last update: 2026-11-07
Python: 3.11+

Dependencies: numpy, pandas, scipy
Public API: CDOMProcessor.compute_absorption, CDOMProcessor.expfit
"""

import numpy as np
import pandas as pd
import scipy.optimize


class CDOMProcessor:
    """Compute CDOM absorption and exponential-fit parameters.

    The processor is attached to a :class:`~watercolorlab.spectra.Spectra`
    instance. Processing results are therefore stored on that shared object,
    while this class provides the CDOM-specific operations.
    """

    def __init__(self, spectra):
        """Create a processor associated with ``spectra``.

        Parameters
        ----------
        spectra : Spectra
            Object containing raw spectral measurements and their metadata.
            The processor writes intermediate and final CDOM results to this
            object.
        """
        self.spectra = spectra

    def __getattr__(self, name):
        """Delegate unknown attribute reads to the associated spectra object.

        This makes attributes such as ``raw_data``, ``wls``, and
        ``absorption`` available through ``spectra.cdom`` without duplicating
        the data on the processor.
        """
        return getattr(self.spectra, name)

    def __setattr__(self, name, value):
        """Store processor state on the associated spectra object.

        The ``spectra`` reference itself must remain on the processor.
        Every other attribute is deliberately written to the shared
        :class:`Spectra` object so all processing modules see the same state.
        """
        if name == 'spectra':
            object.__setattr__(self, name, value)
        else:
            setattr(self.spectra, name, value)

    def compute_absorption(
            self,
            blank=None,
            wl_null_point_correction=600,
            pathlength=.1,
            mean_replicates=True):
        """
        Convert absorbance measurements into CDOM absorption coefficients.

        The calculation optionally subtracts a blank curve, applies a
        null-point correction, and converts absorbance to an absorption
        coefficient using ``2.303 * absorbance / pathlength``. When multiple
        blank curves are supplied, the blank whose long-wavelength response
        is closest to each sample is selected independently.

        Parameters
        ----------
        blank : bool, pd.DataFrame, or None, optional
            If `False`, no blank correction is applied. If `True`, the mean
            blank curves in the dataset are used. If `None`, the mean blank
            curves are used when available; otherwise a zero blank is assumed.
            A DataFrame uses its wavelength columns as the blank reference.
        wl_null_point_correction : int, tuple, or None, optional
            Wavelength or inclusive wavelength interval used to remove the
            null-point offset. A single wavelength must be present in the
            absorbance columns. Interval bounds must be ordered and within
            the available wavelength range. ``None`` disables this correction.
        pathlength : float, optional
            Instrument optical path length in meters. Must be finite and
            greater than zero.
        mean_replicates : bool, optional
            If ``True``, average sample replicates by their ``id`` metadata
            before blank correction. If ``False``, retain individual samples.

        Returns
        -------
        Spectra
            The associated object, with ``raw_samples``, ``absorbance``, and
            ``absorption`` attributes populated.

        Raises
        ------
        ValueError
            If the pathlength, sample data, blank data, or null-point
            correction is invalid.
        TypeError
            If ``blank`` is not ``None``, a boolean, or a DataFrame.
        """
        # Absorption is inversely proportional to pathlength, so a non-positive
        # path would be physically invalid and would make the conversion
        # meaningless.
        if (not isinstance(pathlength, (int, float, np.number))
                or not np.isfinite(pathlength) or pathlength <= 0):
            raise ValueError(
                'pathlength must be a finite number greater than zero'
            )

        # A new absorption computation invalidates results and metadata from
        # previous CDOM processing modes.
        for attribute in (
                'raw_samples', 'absorbance', 'absorption', 'cdom_fitted',
                'blank_curves', 'blank_salinity', 'npc_offset'):
            setattr(self.spectra, attribute, None)

        # Resolve the blank source. ``None`` means "use existing blanks when
        # available"; False creates a zero blank; True averages dataset
        # blanks; a DataFrame supplies caller-selected reference curves.
        if pd.Index(self.wls).has_duplicates:
            raise ValueError(
                'The sample data contain duplicate wavelength columns.'
            )

        if blank is None:
            blank = bool(self.is_blank.any())

        if blank is False:
            print('No blank correction will be applied')
            blank = pd.DataFrame({wl: 0 for wl in self.wls}, index=['blank'])

        elif blank is True:
            print('Blank correction will be applied with existing curves (averaged)')
            if not self.is_blank.any():
                raise ValueError('No blank curves were found in the data.')
            blank = (self.raw_spectra.loc[self.is_blank.to_numpy()]
                     .mean()
                     .rename('blank')
                     .to_frame()
                     .transpose())

        elif isinstance(blank, pd.DataFrame):
            print('Blank correction will be applied with provided salinity curves')
            if blank.empty:
                raise ValueError('The provided blank DataFrame is empty.')
            blank = blank.copy()
            blank.columns = blank.columns.map(int)
            if blank.columns.has_duplicates:
                raise ValueError(
                    'The blank data contain duplicate wavelength columns.'
                )

        else:
            raise TypeError(
                'blank must be one of: None, True, False, or a pandas.DataFrame.'
            )

        # A blank without wavelength columns cannot be aligned with samples.
        if blank.columns.empty:
            raise ValueError('The blank dataset does not contain wavelength columns.')

        # Restrict both datasets to their common integer wavelength range.
        # Using a common range prevents subtraction from silently introducing
        # missing values where one source has no measurement.
        wlmin = max(self.wls.min(), blank.columns.min())
        wlmax = min(self.wls.max(), blank.columns.max())
        if wlmin > wlmax:
            raise ValueError(
                'The sample and blank wavelength ranges do not overlap.'
            )
        wl_range = range(wlmin, wlmax + 1)
        print(f'Using wavelength range {wlmin}-{wlmax} nm for CDOM processing')

        # Build the sample table. In averaged mode, numeric wavelength columns
        # are grouped by station/sample ID; otherwise each original row is
        # retained for later calculations.
        if mean_replicates:
            samples = (self.raw_data.loc[self.is_sample]
                       .groupby('id')
                       .mean(numeric_only=True)
                       .loc[:, self.wls])
        else:
            samples = self.raw_spectra.loc[self.is_sample.to_numpy()]

        if samples.empty:
            raise ValueError('No sample curves were found in the data.')

        samples.columns = samples.columns.map(int)
        self.raw_samples = samples[wl_range]
        blank = blank[wl_range]

        # Apply one common blank directly. With multiple blanks, use the
        # wavelengths at or above 685 nm to choose the reference whose shape
        # is closest to each sample's long-wavelength response.
        if len(blank) == 1:
            self.absorbance = self.raw_samples.sub(blank.squeeze(), axis=1)
        else:
            blank_ref = blank.loc[:, blank.columns.to_numpy() >= 685]
            if blank_ref.empty:
                raise ValueError(
                    'The blank DataFrame must include wavelengths >= 685 nm '
                    'to select the closest blank curve.'
                )

            # Compare samples with candidate blanks only in the long-
            # wavelength region used to identify the closest reference.
            ref_candidates = blank_ref.loc[:, blank_ref.columns.to_numpy() >= 685]
            corr_dict = {}
            blank_dict = {}
            blank_sal_dict = {}

            for station, curve in self.raw_samples.iterrows():
                # Compare every candidate blank against this sample only over
                # the long-wavelength reference region.
                sample_ref = curve.loc[curve.index >= 685]
                blank_sel = (
                    abs(ref_candidates.sub(sample_ref, axis=1))
                    .mean(axis=1)
                    .sort_values()
                    .index[0])
                corr_dict[station] = curve - blank.loc[blank_sel]
                blank_dict[station] = blank.loc[blank_sel]
                blank_sal_dict[station] = blank_sel

            # Preserve the selected blank and its identifying value for users
            # who need to inspect or reproduce the correction.
            self.absorbance = pd.DataFrame(corr_dict).transpose()
            self.blank_curves = pd.DataFrame(blank_dict).transpose()
            self.blank_salinity = pd.DataFrame(blank_sal_dict, index=[0])

        # Estimate and remove the null-point offset. The offset is calculated
        # separately for every sample so one curve cannot affect another.
        check_offset = True
        if isinstance(wl_null_point_correction, (int, float)):
            if wl_null_point_correction not in self.absorbance.columns:
                raise ValueError(
                    f'Wavelength {wl_null_point_correction} '
                    'is not available in the absorbance data.'
                )
            self.npc_offset = self.absorbance[wl_null_point_correction]
            absorbance_corr = self.absorbance.sub(self.npc_offset, axis=0)

        elif isinstance(wl_null_point_correction, (list, tuple, np.ndarray)):
            if len(wl_null_point_correction) != 2:
                raise ValueError(
                    'wl_null_point_correction must be a pair (wli, wlf) '
                    'or a single wavelength.'
                )
            wli = wl_null_point_correction[0]
            wlf = wl_null_point_correction[1]
            if wli > wlf:
                raise ValueError(
                    'The null point interval must satisfy wli <= wlf.'
                )
            available_wavelengths = self.absorbance.columns
            if (wli < available_wavelengths.min()
                    or wlf > available_wavelengths.max()):
                raise ValueError(
                    'The null point interval must be within the available '
                    f'wavelength range ({available_wavelengths.min()}, '
                    f'{available_wavelengths.max()}) nm.'
                )
            correction_interval = self.absorbance.loc[:, wli:wlf]
            if correction_interval.empty or correction_interval.shape[1] == 0:
                raise ValueError(
                    'The null point interval does not contain any available '
                    'wavelengths.'
                )
            self.npc_offset = correction_interval.mean(axis=1)
            absorbance_corr = self.absorbance.sub(self.npc_offset, axis=0)

        elif wl_null_point_correction is None:
            # Keep the blank-corrected absorbance unchanged when the caller
            # explicitly disables null-point correction.
            absorbance_corr = self.absorbance
            check_offset = False
            print('Null-point correction is disabled')

        else:
            raise ValueError('Set a proper parameter to wl_null_point_correction')

        # Report unusually large offsets while still returning the computed
        # result. The threshold is a quality-control warning, not a failure.
        if check_offset and (self.npc_offset > .0015).any():
            print('Warning: there are offsets > 0.0015 AU in the null '
                  'point correction; check the ``npc_offset`` attribute.')

        # Convert the corrected absorbance to absorption coefficients
        self.absorption = np.log(10) * absorbance_corr / pathlength
        return self

    def expfit(self, wl_ref=440, wl_range=(350, 500)):
        """Fit an exponential CDOM model for each sample.

        The fitted model is ``a(lambda) = a(ref) * exp(-slope * (lambda -
        ref))``. Each row is fitted independently, and the resulting
        reference absorption and slope are stored in ``cdom_fitted``.

        Parameters
        ----------
        wl_ref : int, optional
            Reference wavelength used in the model and in the output column
            name ``a<wl_ref>``. It must be present in the spectral columns.
        wl_range : tuple, list, array-like, or None, optional
            Inclusive fitting interval. Both endpoints must be present in the
            spectral columns. ``None`` uses all available wavelengths.

        Returns
        -------
        Spectra
            The associated object, with a ``cdom_fitted`` DataFrame indexed by
            sample ID.

        Raises
        ------
        AttributeError
            If ``compute_absorption`` has not been called.
        TypeError
            If ``wl_range`` is neither ``None`` nor an array-like pair.
        ValueError
            If the requested wavelengths are unavailable, fewer than two
            wavelengths are available for fitting, or a sample contains
            non-finite values or cannot be fitted.
        """
        # Fitting requires the absorption coefficients generated by the
        # previous processing stage.
        if self.absorption is None:
            raise AttributeError(
                'The absorption data are missing. Run `compute_absorption` '
                'before calling `expfit`.'
            )

        # Resolve the fitting interval. ``None`` means the complete available
        # spectrum; an explicit interval is checked before slicing so a
        # misspelled or unavailable wavelength cannot produce an empty fit.
        if wl_range is None:
            wli, wlf = self.wls.min(), self.wls.max()
        elif isinstance(wl_range, (list, tuple, np.ndarray)):
            wli, wlf = map(int, wl_range)
            if wli > wlf:
                raise ValueError(
                    'The wavelength interval must satisfy wli <= wlf.'
                )
            if wli not in self.wls or wlf not in self.wls:
                raise ValueError(
                    'The requested wavelength range is outside bounds.'
                )
        else:
            raise TypeError(
                'wl_range must be a tuple, list, array-like, or None.'
            )

        # The reference wavelength is part of the model, so it must be
        # available even when it lies outside the selected fitting interval.
        if wl_ref not in self.wls:
            raise ValueError(
                f'Wavelength {wl_ref} is not available in the absorbance data.'
            )

        # Define the model passed to scipy.optimize.curve_fit. ``abs_ref`` is
        # the fitted coefficient at the reference wavelength and ``slope`` is
        # the exponential decay constant.
        def exp_curve(wl, abs_ref, slope):
            return abs_ref * np.exp(-slope * (wl - wl_ref))

        # Select the requested columns and prepare output names. The output
        # uses one row per station and one pair of fitted parameters.
        if wl_range is None:
            ydata = self.absorption
        else:
            fit_columns = [
                wavelength for wavelength in self.absorption.columns
                if wli <= wavelength <= wlf
            ]
            ydata = self.absorption.loc[:, fit_columns]
        if ydata.shape[1] < 2:
            raise ValueError(
                'At least two wavelengths are required for exponential fitting.'
            )
        x = ydata.columns.map(int).to_numpy(dtype=float)
        abs_key = f'a{wl_ref:.0f}'
        slope_key = f'slope_{wli:.0f}_{wlf:.0f}'
        fit_dict = {'id': [], abs_key: [], slope_key: []}

        # Fit each station independently. Coercion makes malformed values
        # explicit NaNs; non-finite rows are rejected before optimization.
        for station, y in ydata.iterrows():
            y = pd.to_numeric(y, errors='coerce').to_numpy(dtype=float)
            if not np.all(np.isfinite(y)):
                raise ValueError(
                    f'Non-finite absorbance values were found for sample '
                    f'{station}.'
                )
            try:
                coeffs, _ = scipy.optimize.curve_fit(
                    exp_curve, x, y, (0.1, 0.02), maxfev=20000)
            except RuntimeError as exc:
                raise ValueError(
                    f'Curve fitting failed for sample {station} in the interval '
                    f'[{wli}, {wlf}] nm.'
                ) from exc
            if not np.all(np.isfinite(coeffs)):
                raise ValueError(
                    f'Non-finite fit parameters were returned for sample {station}.'
                )
            fit_dict['id'].append(station)
            fit_dict[abs_key].append(coeffs[0])
            fit_dict[slope_key].append(coeffs[1])

        # Convert the collected parameter lists to the public result table.
        self.cdom_fitted = pd.DataFrame(fit_dict).set_index('id')
        return self
