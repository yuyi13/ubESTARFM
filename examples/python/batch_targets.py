#!/usr/bin/env python3
"""
Script: batch_targets.py
Objective: Demonstrate Python batch prediction after one reference-pair training.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Bundled GeoTIFFs under inst/extdata and --output-dir.
Outputs: Two fused GeoTIFFs and one reusable NPZ model.
Usage: python3 examples/python/batch_targets.py --output-dir examples/outputs
Dependencies: Python package ubestarfm.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ubestarfm import predict_batch, save_model, train


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
    target = data / "MOD11A1_LST_cloudrm_20160218.tif"
    predict_batch(
        model,
        [target, target],
        output_paths=[
            args.output_dir / "python_batch_20160218_a.tif",
            args.output_dir / "python_batch_20160218_b.tif",
        ],
        overwrite=True,
        workers=4,
    )
    save_model(model, args.output_dir / "python_reference_pair_model.npz")
    print(
        "The bundled repository contains one target date, so this tutorial "
        "repeats 20160218 to demonstrate model reuse."
    )


if __name__ == "__main__":
    main()
