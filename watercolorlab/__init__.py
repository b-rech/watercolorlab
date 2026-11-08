"""
Module: __init__.py
Purpose: Expose the public watercolorlab package API.
Author: Bruno Rech
Institution: INPE
Created: 2026-03-12
Python: 3.11+

Dependencies: watercolorlab.spectra, watercolorlab.io
Public API: Spectra, read_files
"""

from .spectra import Spectra
from .io import read_files

__all__ = [
    'Spectra',
    'read_files',
]

__version__ = '0.1.0'