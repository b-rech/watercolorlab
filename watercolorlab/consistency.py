'''
Module: consistency.py
Purpose: Replicate consistency analysis and suspicious-curve filtering.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Updated: 2026-09-27
Python: 3.11+

Dependencies: pandas, matplotlib, seaborn
Public API: ConsistencyProcessor.check
'''

from math import isfinite
from numbers import Real
from pathlib import Path
import hashlib
import re
import warnings

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


class ConsistencyProcessor:
    '''
    Consistency analysis associated with a Spectra instance.
    '''

    def __init__(self, spectra):
        '''
        Create a consistency processor attached to a Spectra object.
        '''
        self.spectra = spectra

    def check(self,
            std_threshold=.1,
            measurements='all',
            groupby='id',
            remove_suspicious=False,
            plot_suspicious=False,
            recursive=False,
            legend=True,
            plot_path=None,
            inplace=False,
            dpi=200):
        '''
        Check replicate spread and optionally remove suspicious curves.

        Parameters
        ----------
        std_threshold : float, optional
            Maximum accepted standard deviation at each wavelength. All
            selected spectral values must be finite. A group is suspicious
            when the sample standard deviation exceeds this threshold at any
            wavelength; no normalization is applied across wavelengths.
        measurements : {'all', 'sample', 'blank'}, optional
            Subset of measurements to analyze.
        groupby : str or list, optional
            Metadata column(s) defining a station or measurement group.
        remove_suspicious : bool, optional
            Remove the highest-deviation curve from each suspicious group.
        plot_suspicious : bool, optional
            Display plots for suspicious groups.
        recursive : bool, optional
            Repeat removal until no suspicious groups remain.
        legend : bool, optional
            Show a legend on suspicious-curve plots.
        plot_path : pathlib.Path, optional
            Directory where suspicious-curve plots are saved.
        inplace : bool, optional
            When ``True``, remove suspicious curves from the associated
            :class:`Spectra` object and return that same object. When
            ``False``, return a filtered replacement instead.
        dpi : int, optional
            Resolution used for suspicious-group plots. Must be positive and
            finite.

        Notes
        -----
        ``spectra.suspicious`` contains the raw DataFrame index labels of
        selected suspicious curves, not their group labels or sample IDs.
        With ``inplace=True``, filtering rebuilds the associated
        :class:`Spectra` state, so previously calculated processor outputs
        are reset and external references to the old processor instances are
        no longer attached to the current object.
        The sample standard deviation (`ddof=1`) is evaluated independently
        at each wavelength. A group is suspicious when any wavelength exceeds
        ``std_threshold``. Exactly one curve with the largest mean absolute
        deviation from the group median is selected per suspicious group and
        pass; ties select the first curve in row order.

        When ``recursive=True``, all selected raw index labels are preserved
        in ``spectra.suspicious`` across passes. If ``plot_path`` is given,
        it is created when necessary. In-place filtering rebuilds
        :class:`Spectra` and resets processor outputs.

        Returns
        -------
        Spectra
            The associated object, or a filtered replacement when removal is
            requested.
        '''

        spectra = self.spectra

        # Fail early so invalid settings do not produce misleading results.
        if not spectra.raw_data.index.is_unique:
            raise ValueError(
                'Consistency analysis requires unique raw_data index labels.'
            )
        if (not isinstance(std_threshold, Real)
            or not isfinite(std_threshold)):
            raise ValueError('std_threshold must be a finite non-negative number')
        if std_threshold < 0:
            raise ValueError('std_threshold must be a finite non-negative number')
        if (isinstance(dpi, bool) or not isinstance(dpi, Real)
            or not isfinite(dpi) or dpi <= 0):
            raise ValueError('dpi must be a finite positive number')
        boolean_options = {
            'remove_suspicious': remove_suspicious,
            'plot_suspicious': plot_suspicious,
            'recursive': recursive,
            'legend': legend,
            'inplace': inplace,
        }
        invalid_boolean_options = [
            name for name, value in boolean_options.items()
            if not isinstance(value, bool)
        ]
        if invalid_boolean_options:
            raise TypeError(
                'Boolean options must be bool: '
                + ', '.join(invalid_boolean_options)
            )
        if recursive and not remove_suspicious:
            raise ValueError(
                'recursive=True requires remove_suspicious=True'
            )
        if measurements not in {'all', 'sample', 'blank'}:
            raise ValueError("measurements must be 'all', 'sample', or 'blank'")
        required_mask = {
            'sample': 'is_sample',
            'blank': 'is_blank',
        }.get(measurements)
        if required_mask and not hasattr(spectra, required_mask):
            raise ValueError(
                f"Spectra input must provide the '{required_mask}' mask "
                f"for measurements='{measurements}'."
            )
        if required_mask:
            mask = getattr(spectra, required_mask)
            if not isinstance(mask, pd.Series):
                raise TypeError(
                    f"Spectra '{required_mask}' must be a pandas Series."
                )
            if not pd.api.types.is_bool_dtype(mask.dtype):
                raise TypeError(
                    f"Spectra '{required_mask}' must contain boolean values."
                )
            if mask.isna().any():
                raise ValueError(
                    f"Spectra '{required_mask}' cannot contain missing values."
                )
            if not mask.index.equals(spectra.raw_data.index):
                raise ValueError(
                    f"Spectra '{required_mask}' must align with raw_data."
                )
        if plot_path is not None:
            plot_path = Path(plot_path)

        # Normalize one grouping column and multiple grouping columns to the
        # same list form. Pandas returns scalar labels for one key and tuples
        # for multiple keys; the rest of this method handles both forms.
        if isinstance(groupby, str):
            group_keys = [groupby]
        else:
            try:
                group_keys = list(groupby)
            except TypeError as error:
                raise TypeError(
                    'groupby must be a column name or an iterable of column names'
                ) from error
        if not group_keys:
            raise ValueError('groupby must contain at least one column')
        missing = set(group_keys).difference(spectra.raw_data.columns)
        if missing:
            raise KeyError(f'groupby columns not found: {sorted(missing)}')

        # Use Spectra's canonical wavelength list instead of selecting columns
        # with a regex, which could accidentally include metadata columns.
        wavelength_columns = [wl for wl in spectra.wls
                              if wl in spectra.raw_data.columns]
        if not wavelength_columns:
            raise ValueError('No wavelength columns are available for analysis')

        # Restrict the analysis to samples, blanks, or the complete dataset.
        # The filtered frame keeps the original row index so suspicious rows
        # can later be removed from the complete raw dataset.
        if measurements == 'sample':
            data = spectra.raw_data.loc[spectra.is_sample].copy()
        elif measurements == 'blank':
            data = spectra.raw_data.loc[spectra.is_blank].copy()
        else:
            data = spectra.raw_data.copy()
        if data.empty:
            raise ValueError(
                f"No rows are available for measurements='{measurements}'."
            )

        # Pandas skips missing values when calculating standard deviation.
        # Reject incomplete spectral curves explicitly so a group cannot be
        # classified as consistent from only a subset of its wavelengths.
        spectral_data = data[wavelength_columns].apply(
            pd.to_numeric, errors='coerce')
        if not np.isfinite(spectral_data.to_numpy(dtype=float)).all():
            raise ValueError(
                'Consistency analysis requires finite values in every '
                'selected spectral curve.'
            )

        # A group is suspicious when the replicate standard deviation exceeds
        # the threshold at least at one wavelength. The sample standard
        # deviation naturally leaves groups with one replicate non-suspicious.
        grouped = data.groupby(group_keys, sort=False, dropna=False)
        insufficient_groups = [
            label for label, group_data in grouped
            if len(group_data) < 2
        ]
        if insufficient_groups:
            warnings.warn(
                'Consistency analysis cannot estimate sample standard '
                f'deviation for groups with fewer than two curves: '
                f'{insufficient_groups}',
                UserWarning,
                stacklevel=2,
            )
        suspicious_groups = []
        for label, group_data in grouped:
            spread = group_data[wavelength_columns].std()
            if spread.gt(std_threshold).any():
                suspicious_groups.append((label, group_data))

        if suspicious_groups:
            # Plot only groups that failed the spread check. Plotting is kept
            # inside this branch so normal analyses do not create figures.
            if plot_suspicious and plot_path:
                if plot_path.exists() and not plot_path.is_dir():
                    raise ValueError(
                        f'plot_path is not a directory: {plot_path}'
                    )
                plot_path.mkdir(parents=True, exist_ok=True)
            for label, group_data in suspicious_groups:
                if plot_suspicious:
                    figure, axis = plt.subplots(figsize=(4, 3), dpi=dpi)
                    plot_data = group_data[wavelength_columns].transpose()
                    plot_data.columns = [
                        f'index={index}' for index in group_data.index
                    ]
                    plot_data.plot(
                        legend=legend, lw=.5, ax=axis)
                    label_text = label if len(group_keys) > 1 else str(label)
                    axis.set(xlabel='Wavelength (nm)',
                             ylabel='Spectral measurement',
                             title=f'Suspicious curve: ID {label_text}')
                    axis.grid(which='both', lw=.25, color='black')
                    if legend:
                        axis.legend(frameon=False)
                    sns.despine(ax=axis)
                    if plot_path:
                        label_text = str(label)
                        safe_label = re.sub(r'[^A-Za-z0-9._-]+', '_',
                                            label_text).strip('._')
                        safe_label = safe_label or 'group'
                        label_hash = hashlib.sha1(
                            label_text.encode('utf-8')
                        ).hexdigest()[:10]
                        filename = (
                            f'suspicious_{safe_label}_{label_hash}.jpg'
                        )
                        output_path = plot_path / filename
                        try:
                            figure.savefig(output_path)
                        except (OSError, ValueError) as error:
                            raise RuntimeError(
                                f'Could not save suspicious-group plot for '
                                f'{label_text!r} to {output_path}.'
                            ) from error
                    plt.show()
                    plt.close(figure)

            # Compare every curve with its group median. The resulting score
            # is the mean absolute deviation across wavelengths, avoiding a
            # dependence on the number of wavelengths. When scores tie,
            # select the first curve in the existing row order.
            median = grouped[wavelength_columns].transform('median')
            deviation_score = (
                data[wavelength_columns].sub(median).abs().mean(axis=1)
            )
            suspicious = []
            for _, group_data in suspicious_groups:
                scores = deviation_score.loc[group_data.index]
                suspicious.append(scores.idxmax())
            spectra.suspicious = suspicious

            # Optionally remove the worst replicate from the original dataset.
            # A recursive run is useful when a second-outlier group becomes
            # visible only after the first suspicious curve is removed.
            if remove_suspicious:
                filtered = spectra.raw_data.drop(
                    index=spectra.suspicious
                )
                if inplace:
                    # Drop rows on the original table, then rebuild the
                    # derived spectra, metadata, and processor attributes.
                    removed = list(spectra.suspicious)
                    spectra.raw_data.drop(index=removed,
                                          inplace=True)
                    spectra.__init__(spectra.raw_data)
                    target = spectra
                else:
                    removed = list(spectra.suspicious)
                    target = type(spectra)(filtered)
                target.suspicious = removed
                if recursive:
                    prior_removed = list(target.suspicious or [])
                    result = target.consistency.check(
                        std_threshold, measurements, groupby,
                        remove_suspicious, plot_suspicious, recursive,
                        legend, plot_path, inplace, dpi)
                    current_removed = list(result.suspicious or [])
                    result.suspicious = list(dict.fromkeys(
                        prior_removed + current_removed
                    ))
                    return result
                return target
            return spectra

        spectra.suspicious = None
        return spectra
