"""
Script: test_cross_language.py
Objective: Enforce R, Python, and published-output numerical compatibility.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Bundled reference rasters, R fixture runner, and legacy fused GeoTIFF.
Outputs: Pytest assertions and comparison diagnostics.
Usage: python3 -m pytest tests/cross_language/test_cross_language.py -v
Dependencies: numpy, pytest, rasterio, Rscript, and both ubESTARFM packages.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import rasterio

from ubestarfm import predict, train


ROOT = Path(__file__).resolve().parents[2]


def read_values(path: Path):
    with rasterio.open(path) as dataset:
        values = dataset.read(1, masked=True).filled(np.nan).astype(np.float64)
        return values, dataset.transform, dataset.crs


def assert_compatible(first, second, tolerance: float = 1e-4) -> None:
    first_values, first_transform, first_crs = first
    second_values, second_transform, second_crs = second
    np.testing.assert_allclose(
        tuple(first_transform),
        tuple(second_transform),
        rtol=0,
        atol=1e-12,
    )
    assert first_crs == second_crs
    np.testing.assert_array_equal(np.isnan(first_values), np.isnan(second_values))
    finite = np.isfinite(first_values) & np.isfinite(second_values)
    assert np.max(np.abs(first_values[finite] - second_values[finite])) <= tolerance


def test_r_python_and_published_outputs_agree(tmp_path) -> None:
    subprocess.run(
        [
            "Rscript",
            str(ROOT / "tests" / "cross_language" / "run_r_fixture.R"),
            str(tmp_path),
        ],
        cwd=ROOT,
        check=True,
    )

    data = ROOT / "inst" / "extdata"
    model = train(
        data / "Landsat_LST_cloudrm_20160205.tif",
        data / "Landsat_LST_cloudrm_20160308.tif",
        data / "MOD11A1_LST_cloudrm_20160205.tif",
        data / "MOD11A1_LST_cloudrm_20160308.tif",
        window_radius=25,
        patch_size=200,
        workers=4,
    )
    python_prediction = predict(
        model,
        data / "MOD11A1_LST_cloudrm_20160218.tif",
        workers=4,
    )

    r_output = read_values(tmp_path / "r_prediction.tif")
    r_reloaded = read_values(tmp_path / "r_prediction_reloaded.tif")
    legacy = read_values(ROOT / "tests" / "fixtures" / "legacy" / "fused_result.tif")
    python_output = (
        python_prediction.values,
        python_prediction.profile["transform"],
        python_prediction.profile["crs"],
    )

    assert_compatible(r_output, r_reloaded, tolerance=0)
    assert_compatible(r_output, legacy)
    assert_compatible(python_output, r_output)
    assert int((tmp_path / "r_model_size_bytes.txt").read_text().strip()) < 80 * 1024**2
