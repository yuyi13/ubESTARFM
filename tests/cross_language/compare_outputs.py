#!/usr/bin/env python3
"""
Script: compare_outputs.py
Objective: Report numerical and spatial differences between two raster outputs.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Two GeoTIFF paths and an optional absolute tolerance.
Outputs: Comparison diagnostics on stdout and a nonzero status on failure.
Usage: python3 tests/cross_language/compare_outputs.py first.tif second.tif
Dependencies: numpy and rasterio.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import rasterio


def read(path: Path):
    with rasterio.open(path) as dataset:
        values = dataset.read(1, masked=True).filled(np.nan).astype(np.float64)
        return values, dataset.transform, dataset.crs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("--tolerance", type=float, default=1e-4)
    args = parser.parse_args()

    first, first_transform, first_crs = read(args.first)
    second, second_transform, second_crs = read(args.second)
    mask_difference = int(np.count_nonzero(np.isnan(first) != np.isnan(second)))
    finite = np.isfinite(first) & np.isfinite(second)
    absolute = np.abs(first[finite] - second[finite])
    maximum = float(absolute.max()) if absolute.size else 0.0
    mean = float(absolute.mean()) if absolute.size else 0.0
    transform_equal = bool(
        np.allclose(
            tuple(first_transform),
            tuple(second_transform),
            rtol=0,
            atol=1e-12,
        )
    )
    print(f"pixels={first.size}")
    print(f"mask_differences={mask_difference}")
    print(f"maximum_absolute_difference={maximum:.12g}")
    print(f"mean_absolute_difference={mean:.12g}")
    print(f"transform_equal={transform_equal}")
    print(f"crs_equal={first_crs == second_crs}")

    return int(
        mask_difference != 0
        or maximum > args.tolerance
        or not transform_equal
        or first_crs != second_crs
    )


if __name__ == "__main__":
    raise SystemExit(main())
