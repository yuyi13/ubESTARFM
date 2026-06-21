#!/usr/bin/env python3
"""
Script: single_target.py
Objective: Demonstrate one Python ubESTARFM prediction from a reusable model.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Bundled GeoTIFFs under inst/extdata and --output-dir.
Outputs: One fused GeoTIFF.
Usage: python3 examples/python/single_target.py --output-dir examples/outputs
Dependencies: Python package ubestarfm.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ubestarfm import predict, train


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    data = root / "inst" / "extdata"
    radius = 2 if args.quick else 25
    model = train(
        data / "Landsat_LST_cloudrm_20160205.tif",
        data / "Landsat_LST_cloudrm_20160308.tif",
        data / "MOD11A1_LST_cloudrm_20160205.tif",
        data / "MOD11A1_LST_cloudrm_20160308.tif",
        window_radius=radius,
        patch_size=200,
        workers=4,
    )
    predict(
        model,
        data / "MOD11A1_LST_cloudrm_20160218.tif",
        output_path=args.output_dir / "python_single_20160218.tif",
        overwrite=True,
        workers=4,
    )


if __name__ == "__main__":
    main()
