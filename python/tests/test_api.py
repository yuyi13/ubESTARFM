"""
Script: test_api.py
Objective: Verify public Python training and raster prediction APIs.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Bundled GeoTIFF fixtures and temporary output paths.
Outputs: Pytest assertions for public Python API behavior.
Usage: python3 -m pytest python/tests/test_api.py -v
Dependencies: numpy, pytest, rasterio, and ubestarfm.
"""

from __future__ import annotations

import numpy as np
import pytest
import rasterio

from ubestarfm import predict, predict_batch, train


def test_public_prediction_preserves_geometry(example_path, tmp_path) -> None:
    model = train(
        example_path("Landsat_LST_cloudrm_20160205.tif"),
        example_path("Landsat_LST_cloudrm_20160308.tif"),
        example_path("MOD11A1_LST_cloudrm_20160205.tif"),
        example_path("MOD11A1_LST_cloudrm_20160308.tif"),
        window_radius=2,
        patch_size=100,
    )
    target = example_path("MOD11A1_LST_cloudrm_20160218.tif")
    output = tmp_path / "prediction.tif"
    single = predict(model, target, output_path=output)
    batch = predict_batch(model, [target, target])

    assert output.exists()
    np.testing.assert_allclose(batch[0].values, batch[1].values, equal_nan=True)
    with rasterio.open(output) as dataset:
        assert dataset.shape == single.values.shape
        assert dataset.transform == single.profile["transform"]
        assert dataset.crs == single.profile["crs"]


def test_prediction_refuses_unrequested_overwrite(example_path, tmp_path) -> None:
    model = train(
        example_path("Landsat_LST_cloudrm_20160205.tif"),
        example_path("Landsat_LST_cloudrm_20160308.tif"),
        example_path("MOD11A1_LST_cloudrm_20160205.tif"),
        example_path("MOD11A1_LST_cloudrm_20160308.tif"),
        window_radius=1,
        patch_size=100,
    )
    output = tmp_path / "prediction.tif"
    output.touch()

    with pytest.raises(FileExistsError):
        predict(
            model,
            example_path("MOD11A1_LST_cloudrm_20160218.tif"),
            output_path=output,
        )
