'''
Module: io.py
Purpose: Read laboratory spectral measurement files.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Python: 3.11+

Dependencies: pandas, parse

Public API: read_files
'''

# %% Dependencies

from pathlib import Path
import pandas as pd
from parse import parse


# %% Function to read files ------------------------------------------------------------

def read_files(path,
               instrument,
               format_string,
               sample_pattern,
               blank_pattern,
               depig_pattern=None,
               decimal='.'):
    '''
    Read and combine spectrophotometric measurement files.

    Parameters
    ----------
    path : str or pathlib.Path
        Directory containing the measurement files. The directory must contain
        only files from the dataset to be read.

    instrument : str
        Instrument that generated the files. Supported values are ``'wpi'``,
        ``'perkinelmer'``, ``'shimadzu'``, and ``'drawell'``.

    format_string : str
        ``parse``-compatible template used to extract metadata from each file
        name. It must contain the fields ``pattern`` and ``id``.

        For example, ``'{pattern}_{id}_{replicate}_{campaign}'`` parses a
        filename such as ``'station_01_a_c01'`` into named fields.

        See the ``parse`` package documentation for template syntax.

    sample_pattern : str
        Value of the parsed ``pattern`` field used to identify sample
        measurements.

    blank_pattern : str
        Value of the parsed ``pattern`` field used to identify blank
        measurements.

    depig_pattern : str, optional
        Substring used to classify depigmented measurements. When provided,
        the returned DataFrame includes ``is_depig`` and ``is_total`` columns.

    decimal : str, default='.'
        Decimal separator used when reading the measurement files.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by ``key``. It contains the metadata parsed from
        each filename, the boolean columns ``is_blank`` and ``is_sample``,
        optional depigmentation flags, and one column per wavelength.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    NotADirectoryError
        If ``path`` is not a directory.
    ValueError
        If the instrument is unsupported, the format is missing ``pattern``
        or ``id``, or a filename cannot be parsed.
    '''

    # -------------------------------------------------------------------------
    # Initial checks
    # -------------------------------------------------------------------------

    # Check format string consistency
    if 'id' not in format_string:
        raise ValueError('The format string must contain an "id" string')
    if 'pattern' not in format_string:
        raise ValueError('The format string must contain a "pattern" string')

    # Update path
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'Path does not exist: {path}')
    if not path.is_dir():
        raise NotADirectoryError(f'Expected a directory, got: {path}')

    supported_instruments = {'wpi', 'shimadzu', 'perkinelmer', 'drawell'}
    if instrument not in supported_instruments:
        raise ValueError(
            'Unsupported instrument. Supported values are: '
            f'{sorted(supported_instruments)}'
        )

    # -------------------------------------------------------------------------
    # Specific processing of WPI data
    # -------------------------------------------------------------------------

    if instrument == 'wpi':

        # Create dict to receive curves
        curve_dict = {}

        # Iterate over files
        for file_path in path.iterdir():

            # Skip .csv files
            if file_path.suffix.lower() == '.csv':
                continue

            # Open file
            raw = pd.read_table(filepath_or_buffer=file_path,
                                decimal=decimal,
                                index_col=0,
                                skiprows=44,
                                header=None,
                                encoding='latin1',
                                engine='python')

            # Add to list
            curve_dict[file_path.stem] = raw.dropna().mean(axis=1)

        # Create dataframe
        curves = pd.DataFrame(curve_dict).transpose()

    # -------------------------------------------------------------------------
    # Specific processing of Shimadzu data
    # -------------------------------------------------------------------------

    elif instrument == 'shimadzu':

        # Create list to store tables
        table_list = []

        # Iterate over files
        for file_path in path.iterdir():

            # Skip .csv files
            if file_path.suffix.lower() == '.csv':
                continue

            # Open and store
            raw = pd.read_table(file_path, decimal=decimal, index_col=0)
            raw.index.rename('wl', inplace=True)
            table_list.append(raw)
            print(f'Reading {file_path}')

        # Concatenate and transpose
        if not table_list:
            raise ValueError(
                f'No files were found in {path} for instrument "{instrument}"'
                )

        curves = pd.concat(table_list, axis=1).transpose()

        # Format wavelengths to integer
        curves.columns = curves.columns.map(int)

    # -------------------------------------------------------------------------
    # Specific processing of PerkinElmer data
    # -------------------------------------------------------------------------
    elif instrument == 'perkinelmer':

        # Create dict to receive curves
        curve_dict = {}

        # Iterate over files
        for file_path in path.iterdir():

            # Open file
            raw = pd.read_csv(filepath_or_buffer=file_path,
                              sep=';',
                              decimal=decimal,
                              index_col=0,
                              skiprows=1)

            # Reverse if wavelengths are in decreasing order
            if raw.index[0] > raw.index[-1]:
                raw = raw[::-1]

            # Add to list
            curve_dict[file_path.stem] = raw.dropna().mean(axis=1)

        # Create dataframe
        curves = pd.DataFrame(curve_dict).transpose()

    # -------------------------------------------------------------------------
    # Specific processing of Drawell data
    # -------------------------------------------------------------------------

    elif instrument == 'drawell':

        # Create dict to receive curves
        curve_dict = {}

        # Iterate over files
        for file_path in path.iterdir():

            # Skip .csv files
            if file_path.suffix.lower() == '.csv':
                continue

            # Open file
            raw = pd.read_table(filepath_or_buffer=file_path,
                                decimal=decimal,
                                sep=r'\s+',
                                skiprows=6,
                                index_col=1,
                                encoding='latin1',
                                engine='python')

            # Add to list
            curve_dict[file_path.stem] = raw['Abs'].dropna()

        # Create dataframe
        curves = pd.DataFrame(curve_dict).transpose()

    # -------------------------------------------------------------------------
    # Further processing of data
    # -------------------------------------------------------------------------

    # Check for NaN
    if curves.isna().any().any():
        print('Warning: there are NaNs in the data')

    # Parse names
    try:
        parsed = curves.index.map(lambda x: parse(format_string, x).named)
    except Exception as exc:
        raise ValueError(
            '\nCould not parse metadata. Check "format_string" '
            'and the consistency of the naming convention'
        ) from exc

    # Create dataframe with metadata
    meta = pd.DataFrame(parsed.tolist(), index=curves.index)

    missing_fields = {'pattern', 'id'} - set(meta.columns)
    if missing_fields:
        raise ValueError(
            'Parsed metadata is missing required fields: '
            f'{sorted(missing_fields)}'
        )

    # Create column to identify sample and blank measurements
    meta['is_blank'] = meta['pattern'] == blank_pattern
    meta['is_sample'] = meta['pattern'] == sample_pattern

    if depig_pattern is not None:
        meta['is_depig'] = [depig_pattern in x for x in curves.index]
        meta['is_total'] = [depig_pattern not in x for x in curves.index]

    # Final dataframe
    df = meta.merge(curves, left_index=True, right_index=True)
    df.index.name = 'key'

    # Return data
    return df

# %%
