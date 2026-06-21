"""
Script: io.py
Objective: Read, validate, and write raster data for Python ubESTARFM.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Raster paths, Rasterio datasets, NumPy arrays, and raster profiles.
Outputs: Validated RasterData values and profiles.
Usage: Imported by the ubestarfm package.
Dependencies: numpy and rasterio.
"""

from __future__ import annotations

from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
from rasterio.io import DatasetReader


@dataclass(frozen=True, slots=True)
class RasterData:
    """A two-dimensional float64 raster array and its spatial profile."""

    values: np.ndarray
    profile: dict[str, Any]

    def __post_init__(self) -> None:
        values = np.ascontiguousarray(self.values, dtype=np.float64)
        if values.ndim != 2:
            raise ValueError("Raster values must be a two-dimensional array.")

        profile = dict(self.profile)
        if int(profile.get("height", -1)) != values.shape[0]:
            raise ValueError("Raster profile height does not match its values.")
        if int(profile.get("width", -1)) != values.shape[1]:
            raise ValueError("Raster profile width does not match its values.")
        if profile.get("transform") is None or profile.get("crs") is None:
            raise ValueError("Raster profile must contain a transform and CRS.")

        object.__setattr__(self, "values", values)
        object.__setattr__(self, "profile", profile)


def _read_dataset(dataset: DatasetReader) -> RasterData:
    if dataset.count != 1:
        raise ValueError("Raster inputs must contain exactly one layer.")
    masked = dataset.read(1, masked=True)
    values = np.ascontiguousarray(masked.filled(np.nan), dtype=np.float64)
    values[~np.isfinite(values)] = np.nan
    profile = dataset.profile.copy()
    profile.update(height=dataset.height, width=dataset.width, count=1)
    return RasterData(values=values, profile=profile)


def read_raster(
    source: str | PathLike[str] | DatasetReader | np.ndarray,
    *,
    profile: dict[str, Any] | None = None,
) -> RasterData:
    """Read one raster layer or validate an in-memory array."""
    if isinstance(source, np.ndarray):
        if profile is None:
            raise ValueError("An explicit raster profile is required for arrays.")
        values = np.ascontiguousarray(source, dtype=np.float64)
        values[~np.isfinite(values)] = np.nan
        return RasterData(values=values, profile=profile)

    if isinstance(source, DatasetReader):
        return _read_dataset(source)

    path = Path(source)
    with rasterio.open(path) as dataset:
        return _read_dataset(dataset)


def validate_geometry(rasters: list[RasterData]) -> None:
    """Require identical shape, affine transform, and coordinate system."""
    if len(rasters) < 2:
        return

    first = rasters[0]
    for raster in rasters[1:]:
        shape_matches = raster.values.shape == first.values.shape
        transform_matches = first.profile["transform"].almost_equals(
            raster.profile["transform"],
            precision=1e-12,
        )
        crs_matches = raster.profile["crs"] == first.profile["crs"]
        if not (shape_matches and transform_matches and crs_matches):
            raise ValueError(
                "All rasters must have the same dimensions, transform, and CRS."
            )
