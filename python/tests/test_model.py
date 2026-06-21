"""
Script: test_model.py
Objective: Verify immutable and versioned Python model objects.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Small valid and malformed in-memory model arrays.
Outputs: Pytest assertions for Python model validation.
Usage: python3 -m pytest python/tests/test_model.py -v
Dependencies: numpy, pytest, and ubestarfm.
"""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest
from rasterio.crs import CRS
from rasterio.transform import from_origin

from ubestarfm.model import UBESTARFMModel


def valid_model() -> UBESTARFMModel:
    return UBESTARFMModel(
        schema_version=1,
        package_version="3.0.0.dev0",
        reference_values=np.zeros((4, 2, 2), dtype=np.float64),
        reference_valid=np.ones((2, 2), dtype=np.bool_),
        candidate_masks=np.zeros((4, 2), dtype=np.uint8),
        patch_ids=np.ones((2, 2), dtype=np.int32),
        patch_thresholds=np.ones((1, 2), dtype=np.float64),
        profile={
            "height": 2,
            "width": 2,
            "transform": from_origin(0, 2, 1, 1),
            "crs": CRS.from_epsg(4326),
            "nodata": np.nan,
            "dtype": "float32",
            "count": 1,
            "driver": "GTiff",
        },
        window_radius=1,
        patch_size=2,
        method="zero_bias",
    )


def test_model_rejects_unsupported_schema() -> None:
    with pytest.raises(ValueError, match="Unsupported model schema"):
        replace(valid_model(), schema_version=99)


def test_model_is_frozen() -> None:
    model = valid_model()
    with pytest.raises(FrozenInstanceError):
        model.method = "baseline"
