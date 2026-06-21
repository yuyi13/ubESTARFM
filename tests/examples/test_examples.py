"""
Script: test_examples.py
Objective: Smoke-test all maintained R and Python examples.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Example scripts, bundled GeoTIFFs, and temporary output directories.
Outputs: Pytest assertions for generated example products.
Usage: python3 -m pytest tests/examples/test_examples.py -v
Dependencies: pytest, rasterio, Rscript, Python, and both ubESTARFM packages.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import rasterio


ROOT = Path(__file__).resolve().parents[2]


def run_example(command: list[str], output_directory: Path) -> None:
    subprocess.run(
        [*command, "--output-dir", str(output_directory), "--quick"],
        cwd=ROOT,
        check=True,
    )


def assert_raster(path: Path) -> None:
    assert path.exists()
    with rasterio.open(path) as dataset:
        assert dataset.shape == (400, 400)
        assert dataset.crs.to_epsg() == 4326


def test_r_examples(tmp_path) -> None:
    single = tmp_path / "r-single"
    batch = tmp_path / "r-batch"
    single.mkdir()
    batch.mkdir()
    run_example(["Rscript", "examples/R/single_target.R"], single)
    run_example(["Rscript", "examples/R/batch_targets.R"], batch)

    assert_raster(single / "r_single_20160218.tif")
    assert_raster(batch / "r_batch_20160218_a.tif")
    assert_raster(batch / "r_batch_20160218_b.tif")


def test_python_examples(tmp_path) -> None:
    single = tmp_path / "python-single"
    batch = tmp_path / "python-batch"
    single.mkdir()
    batch.mkdir()
    run_example(["python3", "examples/python/single_target.py"], single)
    run_example(["python3", "examples/python/batch_targets.py"], batch)

    assert_raster(single / "python_single_20160218.tif")
    assert_raster(batch / "python_batch_20160218_a.tif")
    assert_raster(batch / "python_batch_20160218_b.tif")
