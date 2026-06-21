"""
Script: test_io.py
Objective: Verify Python raster loading and geometry validation.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Bundled GeoTIFFs and in-memory RasterData objects.
Outputs: Pytest assertions for Python raster I/O behavior.
Usage: python3 -m pytest python/tests/test_io.py -v
Dependencies: numpy, pytest, rasterio, and ubestarfm.
"""

from __future__ import annotations

import numpy as np
import pytest
from rasterio.crs import CRS
from rasterio.transform import from_origin

from ubestarfm.io import RasterData, read_raster, validate_geometry


def test_read_raster_returns_float64_row_major_array(example_path) -> None:
    raster = read_raster(example_path("MOD11A1_LST_cloudrm_20160218.tif"))

    assert raster.values.shape == (400, 400)
    assert raster.values.dtype == np.float64
    assert raster.values.flags.c_contiguous
    assert raster.values[0, 0] == pytest.approx(312.74, abs=1e-4)
    assert np.isfinite(raster.values).all()


def test_geometry_mismatch_is_rejected() -> None:
    first = RasterData(
        values=np.ones((3, 3), dtype=np.float64),
        profile={
            "height": 3,
            "width": 3,
            "transform": from_origin(0, 3, 1, 1),
            "crs": CRS.from_epsg(4326),
            "nodata": np.nan,
            "dtype": "float64",
            "count": 1,
            "driver": "GTiff",
        },
    )
    second = RasterData(
        values=np.ones((4, 3), dtype=np.float64),
        profile={
            **first.profile,
            "height": 4,
            "transform": from_origin(0, 4, 1, 1),
        },
    )

    with pytest.raises(ValueError, match="same dimensions, transform, and CRS"):
        validate_geometry([first, second])


def test_geometry_accepts_floating_point_equivalent_transforms() -> None:
    first = RasterData(
        values=np.ones((2, 2), dtype=np.float64),
        profile={
            "height": 2,
            "width": 2,
            "transform": from_origin(146, -34.8, 0.001, 0.001),
            "crs": CRS.from_epsg(4326),
            "nodata": np.nan,
            "dtype": "float64",
            "count": 1,
            "driver": "GTiff",
        },
    )
    second = RasterData(
        values=np.ones((2, 2), dtype=np.float64),
        profile={
            **first.profile,
            "transform": from_origin(
                146,
                -34.8,
                0.001,
                0.001 + 2e-17,
            ),
        },
    )

    validate_geometry([first, second])
