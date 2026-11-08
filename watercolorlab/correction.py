"""
Module: correction.py
Purpose: Correct spectral breaks in absorbance measurements.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Python: 3.11+

Dependencies: numpy, pandas, matplotlib, seaborn
Public API: BreakMethods.treat_breaks
"""

import numpy as np
import pandas as pd
import seaborn as sns


class BreakMethods:
    '''
    Break correction methods associated with a :class:`Spectra` object.
    '''

    def __init__(self, spectra):
        '''
        Create break correction methods attached to a Spectra instance.
        '''
        self.spectra = spectra

    def treat_breaks(self,
                     break_centers=(795, 684, 558, 422),
                     break_width=(1, 1),
                     edge_width=2,
                     correct_775_785=True,
                     inline=False,
                     dpi=200):
        '''
        Correct spectral breaks and return a new Spectra object.

        Parameters
        ----------
        break_centers : iterable of int, optional
            Approximate center wavelength for each break.
        break_width : tuple of int, optional
            Width of the interpolation interval as ``(above, below)`` the
            center wavelength. For example, ``break_width=(1, 1)`` around
            795 nm interpolates 796, 795, and 794 nm. The default is
            ``(1, 1)``.
        edge_width : int, optional
            Number of wavelengths used on each side to estimate the offset.
            With ``edge_width=2`` around a break from 796 to 794 nm, the
            lower edge uses 792-794 nm and the upper edge uses 796-798 nm.
            Pandas includes both endpoints in these ranges. The default is
            ``2``.
        correct_775_785 : bool, optional
            Whether to interpolate the fixed 775-785 nm correction range.
        inline : bool, optional
            Display plots inline, such as in a Jupyter notebook, instead of
            opening blocking matplotlib windows. The default is ``False``.
        dpi : int, optional
            Dots per inch for the generated figures. The default is ``200``.
        '''

        plt = self._get_pyplot(inline)

        upper_width, lower_width = break_width
        if upper_width < 0 or lower_width < 0:
            raise ValueError('break_width values must not be negative')
        if edge_width < 0:
            raise ValueError('edge_width must not be negative')

        break_labels = {
            chr(ord('A') + index): center
            for index, center in enumerate(break_centers)
        }

        # Work on a copy so the source Spectra remains unchanged.
        raw_data = self.spectra.raw_data.copy(deep=True)
        wavelength_columns = []
        wavelength_values = []

        # Keep original column names while using numeric wavelengths internally.
        for column in raw_data.columns:
            try:
                wavelength = float(column)
            except (TypeError, ValueError):
                continue
            if wavelength.is_integer():
                wavelength_columns.append(column)
                wavelength_values.append(int(wavelength))

        if not wavelength_columns:
            raise ValueError('No wavelength columns were found')

        wavelength_index = pd.Index(wavelength_values)
        corrected_positions = [
            raw_data.columns.get_loc(column) for column in wavelength_columns
        ]

        for row_position in range(len(raw_data)):
            values = raw_data.iloc[row_position, corrected_positions].to_numpy()
            spectrum = pd.Series(values, index=wavelength_index, dtype=float)
            original_spectrum = spectrum.copy()
            station_key = raw_data.index[row_position]

            if correct_775_785:
                self._interpolate_gap(spectrum, (774, 786), (776, 784))

            figure, axis = plt.subplots(figsize=(5.75, 5.75 / 2), dpi=dpi)
            axis.plot(original_spectrum.index.to_numpy(),
                      original_spectrum.to_numpy(),
                      lw=.5, color='red', label='Original')
            self._mark_breaks(axis, original_spectrum, break_labels, upper_width,
                              lower_width)
            axis.set_xlabel('Wavelength (nm)')
            axis.set_ylabel('Absorbance (AU)')
            axis.set_title(str(station_key))
            axis.legend(frameon=False, facecolor='white', framealpha=1)
            sns.despine(ax=axis)
            self._show_figure(figure, inline)

            selected = input(
                'Indicate which break to correct '
                f'({", ".join(break_labels)}, [None]): '
            )

            for label in selected.upper():
                if label not in break_labels:
                    continue

                center = break_labels[label]
                # Build the wavelength interval containing the instrument break.
                break_wavelengths = np.arange(
                    center + upper_width,
                    center - lower_width - 1,
                    -1
                )
                missing = pd.Index(break_wavelengths).difference(
                    wavelength_index
                )
                if len(missing):
                    raise ValueError(
                        f'Break near {center} nm is outside the available '
                        f'wavelength range'
                    )

                wavelength_high = int(break_wavelengths[0])
                wavelength_low = int(break_wavelengths[-1])
                low_edge = spectrum.loc[
                    wavelength_low - edge_width:wavelength_low
                ].mean()
                high_edge = spectrum.loc[
                    wavelength_high:wavelength_high + edge_width
                ].mean()
                # Shift the lower-wavelength segment to align both sides.
                offset = (2 * spectrum.loc[wavelength_high]
                          - low_edge - high_edge)
                spectrum.loc[:wavelength_low] += offset

                # Replace the interior break values with a straight line.
                interpolation = spectrum.loc[
                    sorted(break_wavelengths)
                ].copy()
                interpolation.iloc[1:-1] = np.nan
                interpolation.interpolate(method='linear', inplace=True)
                spectrum.loc[interpolation.index] = interpolation

            raw_data.iloc[row_position, corrected_positions] = spectrum.to_numpy()

            figure, axis = plt.subplots(figsize=(5.75, 5.75 / 2), dpi=300)
            axis.plot(original_spectrum.index.to_numpy(),
                      original_spectrum.to_numpy(),
                      lw=.5, color='red', label='Original')
            axis.plot(spectrum.index.to_numpy(), spectrum.to_numpy(),
                      lw=.5, color='green', label='Corrected')
            self._mark_breaks(axis, spectrum, break_labels, upper_width,
                              lower_width)
            axis.set_xlabel('Wavelength (nm)')
            axis.set_ylabel('Absorbance (AU)')
            axis.set_title(str(station_key))
            axis.legend(frameon=False, facecolor='white', framealpha=1)
            sns.despine(ax=axis)
            self._show_figure(figure, inline)

        return type(self.spectra)(raw_data)

    @staticmethod
    def _interpolate_gap(spectrum, interpolation_range, gap):
        """Interpolate a fixed interior gap between two known endpoints."""
        range_start, range_end = interpolation_range
        gap_start, gap_end = gap
        interval = spectrum.loc[range_start:range_end].copy()
        interval.loc[gap_start:gap_end] = np.nan
        interval.interpolate(method='linear', inplace=True)
        spectrum.loc[interval.index] = interval

    @staticmethod
    def _mark_breaks(axis, spectrum, break_labels, upper_width, lower_width):
        """Shade and label possible break intervals on a spectrum plot."""
        for label, center in break_labels.items():
            high = center + upper_width
            low = center - lower_width
            axis.axvspan(low, high, alpha=.5, facecolor='gray',
                         edgecolor='none', linewidth=0)
            axis.text((low + high) / 2, spectrum.max(), label,
                      ha='center', va='bottom',
                      bbox={'facecolor': 'white', 'edgecolor': 'none',
                            'pad': 1})

    @staticmethod
    def _show_figure(figure, inline):
        """Display a figure inline or in a blocking window."""
        if inline:
            from IPython.display import display
            display(figure)
        else:
            import matplotlib.pyplot as plt
            plt.show(block=True)
            plt.close(figure)

    @staticmethod
    def _get_pyplot(inline):
        """Return pyplot using an inline or interactive GUI backend."""
        import matplotlib

        if not inline and matplotlib.get_backend().lower() != 'qtagg':
            try:
                matplotlib.use('QtAgg', force=True)
            except ImportError:
                # Keep the configured backend when Qt bindings are unavailable.
                pass

        import matplotlib.pyplot as plt
        return plt
