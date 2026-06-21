#!/usr/bin/env python3
"""
Script: benchmark.py
Objective: Benchmark Python ubESTARFM training, prediction, and model memory.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Bundled GeoTIFFs and --output-dir.
Outputs: Python benchmark JSON and Markdown reports.
Usage: python3 benchmarks/benchmark.py --output-dir benchmarks/results
Dependencies: Python package ubestarfm.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from ubestarfm import predict, predict_batch, train


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = Path(__file__).resolve().parents[1]
    data = root / "inst" / "extdata"
    started = time.perf_counter()
    model = train(
        data / "Landsat_LST_cloudrm_20160205.tif",
        data / "Landsat_LST_cloudrm_20160308.tif",
        data / "MOD11A1_LST_cloudrm_20160205.tif",
        data / "MOD11A1_LST_cloudrm_20160308.tif",
        workers=4,
    )
    training_seconds = time.perf_counter() - started

    target = data / "MOD11A1_LST_cloudrm_20160218.tif"
    started = time.perf_counter()
    predict(model, target, workers=4)
    single_seconds = time.perf_counter() - started
    started = time.perf_counter()
    predict_batch(model, [target, target], workers=4)
    batch_seconds = time.perf_counter() - started

    metrics = {
        "training_seconds": training_seconds,
        "single_prediction_seconds": single_seconds,
        "two_target_batch_seconds": batch_seconds,
        "model_size_bytes": model.nbytes,
    }
    (args.output_dir / "python_benchmark.json").write_text(
        json.dumps(metrics, indent=2) + "\n"
    )
    (args.output_dir / "python_benchmark.md").write_text(
        "\n".join(
            [
                "# Python benchmark",
                "",
                f"- Training: {training_seconds:.3f} seconds",
                f"- Warm single prediction: {single_seconds:.3f} seconds",
                f"- Two-target batch: {batch_seconds:.3f} seconds",
                f"- Model size: {model.nbytes / 1024**2:.2f} MiB",
                "",
            ]
        )
    )


if __name__ == "__main__":
    main()
