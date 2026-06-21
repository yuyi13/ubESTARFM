"""
Script: test_predict_kernel.py
Objective: Verify Python predictions against independent ubESTARFM equations.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Small trained models and synthetic coarse target arrays.
Outputs: Pytest assertions for detailed, fallback, and batch predictions.
Usage: python3 -m pytest python/tests/test_predict_kernel.py -v
Dependencies: numpy, pytest, and ubestarfm.
"""

from __future__ import annotations

import numpy as np
import pytest

from ubestarfm.api import predict_arrays, train_arrays, unpack_candidates


def reference_predict_pixel(model, target, row, col, value_range):
    radius = model.window_radius
    rows = slice(max(0, row - radius), min(target.shape[0], row + radius + 1))
    cols = slice(max(0, col - radius), min(target.shape[1], col + radius + 1))
    fine_1, fine_2, coarse_1, coarse_2 = model.reference_values
    window_valid = model.reference_valid[rows, cols] & np.isfinite(target[rows, cols])

    candidate_rows, candidate_cols = unpack_candidates(model, row=row, col=col)
    keep = np.isfinite(target[candidate_rows, candidate_cols])
    candidate_rows = candidate_rows[keep]
    candidate_cols = candidate_cols[keep]
    index = (candidate_rows, candidate_cols)

    target_window = target[rows, cols][window_valid]
    coarse_1_window = coarse_1[rows, cols][window_valid]
    coarse_2_window = coarse_2[rows, cols][window_valid]

    if candidate_rows.size > 5:
        fine_candidates_1 = fine_1[index].copy()
        fine_candidates_2 = fine_2[index].copy()
        coarse_candidates_1 = coarse_1[index]
        coarse_candidates_2 = coarse_2[index]
        fine_pixel_1 = fine_1[row, col]
        fine_pixel_2 = fine_2[row, col]
        if model.method == "zero_bias":
            bias_1 = -fine_candidates_1.mean() + coarse_candidates_1.mean()
            bias_2 = -fine_candidates_2.mean() + coarse_candidates_2.mean()
            fine_candidates_1 += bias_1
            fine_candidates_2 += bias_2
            fine_pixel_1 += bias_1
            fine_pixel_2 += bias_2

        spectral = 1 - 0.5 * (
            np.abs(
                (fine_candidates_1 - coarse_candidates_1)
                / (fine_candidates_1 + coarse_candidates_1)
            )
            + np.abs(
                (fine_candidates_2 - coarse_candidates_2)
                / (fine_candidates_2 + coarse_candidates_2)
            )
        )
        spectral[(spectral > 1) | (spectral < -1)] = 0.5
        spatial = 1 + np.sqrt((col - candidate_cols) ** 2 + (row - candidate_rows) ** 2) / radius
        inverse = 1 / ((1 - spectral) * spatial + 1e-7)
        weight = inverse / inverse.sum()

        difference_1 = abs(np.mean(target_window - coarse_1_window)) + 1e-10
        difference_2 = abs(np.mean(target_window - coarse_2_window)) + 1e-10
        temporal_1 = (1 / difference_1) / (1 / difference_1 + 1 / difference_2)
        temporal_2 = 1 - temporal_1
        target_candidates = target[index]
        prediction_1 = fine_pixel_1 + np.sum(weight * (target_candidates - coarse_candidates_1))
        prediction_2 = fine_pixel_2 + np.sum(weight * (target_candidates - coarse_candidates_2))
        prediction = temporal_1 * prediction_1 + temporal_2 * prediction_2
        if prediction <= value_range[0] or prediction >= value_range[1]:
            prediction = temporal_1 * np.sum(weight * fine_candidates_1) + temporal_2 * np.sum(
                weight * fine_candidates_2
            )
        return prediction

    fine_pixel_1 = fine_1[row, col]
    fine_pixel_2 = fine_2[row, col]
    if model.method == "zero_bias":
        fine_pixel_1 = (
            fine_pixel_1 - fine_1[rows, cols][window_valid].mean() + coarse_1_window.mean()
        )
        fine_pixel_2 = (
            fine_pixel_2 - fine_2[rows, cols][window_valid].mean() + coarse_2_window.mean()
        )
    difference_1 = np.mean(target_window - coarse_1_window) + 1e-10
    difference_2 = np.mean(target_window - coarse_2_window) + 1e-10
    temporal_1 = (1 / abs(difference_1)) / (1 / abs(difference_1) + 1 / abs(difference_2))
    return temporal_1 * (fine_pixel_1 + difference_1) + (1 - temporal_1) * (
        fine_pixel_2 + difference_2
    )


@pytest.mark.parametrize("method", ["zero_bias", "baseline"])
def test_prediction_matches_independent_equations(small_reference_fixture, method) -> None:
    model = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
        method=method,
    )
    target = small_reference_fixture["coarse_1"] + 0.5
    expected = reference_predict_pixel(model, target, 2, 2, (-np.inf, np.inf))
    actual = predict_arrays(model, [target], (-np.inf, np.inf))[0]

    assert actual[2, 2] == pytest.approx(expected, abs=1e-10)


def test_missing_targets_trigger_candidate_filtering_and_fallback(
    small_reference_fixture,
) -> None:
    model = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
    )
    target = small_reference_fixture["coarse_1"] + 0.5
    target[1:4, 1:4] = np.nan
    target[2, 2] = small_reference_fixture["coarse_1"][2, 2] + 0.5
    expected = reference_predict_pixel(model, target, 2, 2, (-np.inf, np.inf))
    actual = predict_arrays(model, [target], (-np.inf, np.inf))[0]

    assert actual[2, 2] == pytest.approx(expected, abs=1e-10)


def test_single_and_batch_predictions_are_identical(small_reference_fixture) -> None:
    model = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
    )
    first = small_reference_fixture["coarse_1"] + 0.5
    second = small_reference_fixture["coarse_1"] + 1
    batch = predict_arrays(model, [first, second], (-np.inf, np.inf))
    single = predict_arrays(model, [first], (-np.inf, np.inf))

    np.testing.assert_allclose(batch[0], single[0], equal_nan=True)
