"""
Module: plots.py
Purpose: Plot raw, corrected, and particulate absorption spectra.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Python: 3.11+

Dependencies: pandas, matplotlib, seaborn
Public API: plot_part_absorption_raw, plot_part_absorption_absorbance
"""

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def _add_station_index(table: pd.DataFrame,
                       station_ids: pd.Series) -> pd.DataFrame:
    """Add station IDs as a temporary index level without changing keys."""
    stations = station_ids.reindex(table.index)
    if stations.isna().any():
        raise ValueError('station_ids must contain every curve index.')

    indexed = table.copy()
    indexed.index = pd.MultiIndex.from_arrays(
        [indexed.index, stations.to_numpy()],
        names=['key', 'station'])
    return indexed


def plot_part_absorption_raw(
        trans_total_raw: pd.DataFrame,
        trans_depig_raw: pd.DataFrame,
        refle_total_raw: pd.DataFrame,
        refle_depig_raw: pd.DataFrame,
        trans_total: pd.DataFrame,
        trans_depig: pd.DataFrame,
        refle_total: pd.DataFrame,
        refle_depig: pd.DataFrame,
        Tblank_total: pd.Series,
        Tblank_depig: pd.Series,
        Rblank_total: pd.Series,
        Rblank_depig: pd.Series,
        station_ids: pd.Series,
        trans_label='T',
        reflect_label='R'):
    """Plot raw and offset-corrected particulate curves per station.

    The first panel shows the input curves; the second shows curves after
    transmittance and pigmentation offsets have been applied.
    """
    raw_curves = tuple(_add_station_index(table, station_ids) for table in (
        trans_total_raw, trans_depig_raw, refle_total_raw, refle_depig_raw))
    corrected_curves = tuple(_add_station_index(table, station_ids) for table in (
        trans_total, trans_depig, refle_total, refle_depig))
    sample_ids = pd.Index([
        station
        for curve in raw_curves
        for station in curve.index.get_level_values('station')
    ]).drop_duplicates()

    for sample_id in sample_ids:
        _, axes = plt.subplots(1, 2, figsize=(5.75, 3), dpi=200,
                               sharey=False, constrained_layout=True)

        raw_plot_curves = (
            (raw_curves[0].xs(sample_id, level='station'),
             'blue', '-', f'{trans_label} total'),
            (raw_curves[1].xs(sample_id, level='station'),
             'blue', '--', f'{trans_label} depigmented'),
            (raw_curves[2].xs(sample_id, level='station'),
             'red', '-', f'{reflect_label} total'),
            (raw_curves[3].xs(sample_id, level='station'),
             'red', '--', f'{reflect_label} depigmented'))
        corrected_plot_curves = (
            (corrected_curves[0].xs(sample_id, level='station'),
             'blue', '-', f'{trans_label} total'),
            (corrected_curves[1].xs(sample_id, level='station'),
             'blue', '--', f'{trans_label} depigmented'),
            (corrected_curves[2].xs(sample_id, level='station'),
             'red', '-', f'{reflect_label} total'),
            (corrected_curves[3].xs(sample_id, level='station'),
             'red', '--', f'{reflect_label} depigmented'))

        for curves, axis in ((raw_plot_curves, axes[0]),
                     (corrected_plot_curves, axes[1])):
            for curve, color, linestyle, label in curves:
                curve.transpose().plot(ax=axis, color=color, lw=.5,
                                       linestyle=linestyle,
                                       legend=False, label=label)

        for axis in axes:
            Tblank_total.plot(ax=axis, color='green', lw=.7,
                              label=f'{trans_label} blank')
            Tblank_depig.plot(ax=axis, color='green', lw=.7,
                              ls='--',
                              label=f'{trans_label} depigmented blank')
            Rblank_total.plot(ax=axis, color='black', lw=.7,
                              label=f'{reflect_label} blank')
            Rblank_depig.plot(ax=axis, color='black', lw=.7,
                              ls='--',
                              label=f'{reflect_label} depigmented blank')

        axes[0].set(title=f'ID {sample_id} - raw')
        axes[1].set(title=f'ID {sample_id} - after offsets')
        axes[0].legend(fontsize=6, frameon=False)
        sns.despine()
        plt.show()
        plt.close(axes[0].figure)


def plot_part_absorption_absorbance(
        absorbance_total: pd.DataFrame,
    absorbance_depig: pd.DataFrame,
    station_ids: pd.Series):
    """Plot total and depigmented absorbance curves per station."""
    total = _add_station_index(absorbance_total, station_ids)
    depig = _add_station_index(absorbance_depig, station_ids)
    sample_ids = pd.Index(total.index.get_level_values('station')).drop_duplicates()

    for sample_id in sample_ids:
        _, axis = plt.subplots(figsize=(3, 2), dpi=300,
                               constrained_layout=True)

        total.xs(sample_id, level='station').transpose().plot(
            ax=axis, color='blue', lw=.5, legend=False, label='Total')
        depig.xs(sample_id, level='station').transpose().plot(
            ax=axis, color='blue', lw=.5, linestyle='--',
            legend=False, label='Depigmented')

        axis.set(title=f'ID {sample_id}', xlabel='Wavelength (nm)',
                 ylabel='Corrected absorbance (AU)')
        axis.legend(frameon=False)
        axis.grid(which='both', color='black', lw=.5)
        sns.despine()
        plt.show()
        plt.close(axis.figure)
