"""
Script: test_candidate_masks.py
Objective: Verify Python packed candidate membership and deterministic training.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Small in-memory reference arrays and bundled GeoTIFFs.
Outputs: Pytest assertions for Numba candidate masks and model memory.
Usage: python3 -m pytest python/tests/test_candidate_masks.py -v
Dependencies: numpy, pytest, and ubestarfm.
"""

from __future__ import annotations

import numpy as np

from ubestarfm.api import train, train_arrays, unpack_candidates


def test_candidate_bits_use_column_major_window_order(small_reference_fixture) -> None:
    model = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
    )
    rows, cols = unpack_candidates(model, row=2, col=2)

    np.testing.assert_array_equal(rows, np.array([2, 3, 1, 2, 3, 1, 2, 3]))
    np.testing.assert_array_equal(cols, np.array([1, 1, 2, 2, 2, 3, 3, 3]))


def test_masks_are_method_independent(small_reference_fixture) -> None:
    zero_bias = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
        method="zero_bias",
    )
    baseline = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
        method="baseline",
    )

    np.testing.assert_array_equal(zero_bias.candidate_masks, baseline.candidate_masks)


def test_parallel_training_is_deterministic(small_reference_fixture) -> None:
    sequential = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=2,
        workers=1,
    )
    parallel = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=2,
        workers=2,
    )

    np.testing.assert_array_equal(sequential.candidate_masks, parallel.candidate_masks)
    np.testing.assert_allclose(sequential.patch_thresholds, parallel.patch_thresholds)


def test_bundled_model_stays_below_80_mib(example_path) -> None:
    model = train(
        example_path("Landsat_LST_cloudrm_20160205.tif"),
        example_path("Landsat_LST_cloudrm_20160308.tif"),
        example_path("MOD11A1_LST_cloudrm_20160205.tif"),
        example_path("MOD11A1_LST_cloudrm_20160308.tif"),
        window_radius=25,
        patch_size=200,
        workers=4,
    )

    assert model.nbytes < 80 * 1024**2
    assert model.candidate_masks.shape == (400 * 400, 326)
