"""
Script: __init__.py
Objective: Expose the maintained Python ubESTARFM public API.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Python package modules.
Outputs: Public package names and version metadata.
Usage: Imported by the ubestarfm package.
Dependencies: Python standard library.
"""

from __future__ import annotations

from .io import RasterData, read_raster, validate_geometry
from .model import UBESTARFMModel

__all__ = [
    "RasterData",
    "UBESTARFMModel",
    "read_raster",
    "validate_geometry",
]
__version__ = "3.0.0.dev0"
