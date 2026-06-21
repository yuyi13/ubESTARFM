"""
Script: kernels.py
Objective: Compile candidate-mask and prediction kernels for Python ubESTARFM.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Contiguous reference arrays, packed masks, target arrays, and patches.
Outputs: Packed candidate masks and fine-resolution prediction arrays.
Usage: Imported by the ubestarfm package.
Dependencies: numba and numpy.
"""

from __future__ import annotations

import math

import numba
import numpy as np


@numba.njit(cache=True, nogil=True)
def train_patch(
    fine_1: np.ndarray,
    fine_2: np.ndarray,
    coarse_1: np.ndarray,
    coarse_2: np.ndarray,
    row_start: int,
    row_end: int,
    col_start: int,
    col_end: int,
    window_radius: int,
    threshold_1: float,
    threshold_2: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Create packed candidate masks for one half-open patch."""
    rows, cols = fine_1.shape
    side = 2 * window_radius + 1
    mask_bytes = (side * side + 7) // 8
    patch_cells = (row_end - row_start) * (col_end - col_start)
    cell_ids = np.empty(patch_cells, dtype=np.int64)
    masks = np.zeros((patch_cells, mask_bytes), dtype=np.uint8)

    local_cell = 0
    for row in range(row_start, row_end):
        for col in range(col_start, col_end):
            cell_ids[local_cell] = row * cols + col
            center_valid = (
                np.isfinite(fine_1[row, col])
                and np.isfinite(fine_2[row, col])
                and np.isfinite(coarse_1[row, col])
                and np.isfinite(coarse_2[row, col])
            )
            if center_valid:
                candidate_row_start = max(0, row - window_radius)
                candidate_row_end = min(rows, row + window_radius + 1)
                candidate_col_start = max(0, col - window_radius)
                candidate_col_end = min(cols, col + window_radius + 1)
                for candidate_col in range(candidate_col_start, candidate_col_end):
                    for candidate_row in range(candidate_row_start, candidate_row_end):
                        candidate_valid = (
                            np.isfinite(fine_1[candidate_row, candidate_col])
                            and np.isfinite(fine_2[candidate_row, candidate_col])
                            and np.isfinite(coarse_1[candidate_row, candidate_col])
                            and np.isfinite(coarse_2[candidate_row, candidate_col])
                        )
                        if not candidate_valid:
                            continue
                        if (
                            abs(fine_1[candidate_row, candidate_col] - fine_1[row, col])
                            < threshold_1
                            and abs(fine_2[candidate_row, candidate_col] - fine_2[row, col])
                            < threshold_2
                        ):
                            delta_row = candidate_row - row
                            delta_col = candidate_col - col
                            bit = (delta_col + window_radius) * side + (
                                delta_row + window_radius
                            )
                            masks[local_cell, bit // 8] |= np.uint8(1 << (bit % 8))
            local_cell += 1

    return cell_ids, masks


@numba.njit(cache=True, nogil=True)
def predict_patch(
    reference_values: np.ndarray,
    candidate_masks: np.ndarray,
    targets: np.ndarray,
    row_start: int,
    row_end: int,
    col_start: int,
    col_end: int,
    window_radius: int,
    zero_bias: bool,
    value_min: float,
    value_max: float,
) -> np.ndarray:
    """Predict every target for one half-open patch."""
    fine_1 = reference_values[0]
    fine_2 = reference_values[1]
    coarse_1 = reference_values[2]
    coarse_2 = reference_values[3]
    target_count, rows, cols = targets.shape
    side = 2 * window_radius + 1
    maximum_candidates = side * side
    output = np.full(
        (target_count, row_end - row_start, col_end - col_start),
        np.nan,
        dtype=np.float64,
    )

    for row in range(row_start, row_end):
        for col in range(col_start, col_end):
            center_valid = (
                np.isfinite(fine_1[row, col])
                and np.isfinite(fine_2[row, col])
                and np.isfinite(coarse_1[row, col])
                and np.isfinite(coarse_2[row, col])
            )
            if not center_valid:
                continue

            cell = row * cols + col
            window_row_start = max(0, row - window_radius)
            window_row_end = min(rows, row + window_radius + 1)
            window_col_start = max(0, col - window_radius)
            window_col_end = min(cols, col + window_radius + 1)

            for target_index in range(target_count):
                if not np.isfinite(targets[target_index, row, col]):
                    continue

                target_window_sum = 0.0
                coarse_1_window_sum = 0.0
                coarse_2_window_sum = 0.0
                fine_1_window_sum = 0.0
                fine_2_window_sum = 0.0
                window_count = 0
                for window_col in range(window_col_start, window_col_end):
                    for window_row in range(window_row_start, window_row_end):
                        if not (
                            np.isfinite(fine_1[window_row, window_col])
                            and np.isfinite(fine_2[window_row, window_col])
                            and np.isfinite(coarse_1[window_row, window_col])
                            and np.isfinite(coarse_2[window_row, window_col])
                            and np.isfinite(targets[target_index, window_row, window_col])
                        ):
                            continue
                        target_window_sum += targets[target_index, window_row, window_col]
                        coarse_1_window_sum += coarse_1[window_row, window_col]
                        coarse_2_window_sum += coarse_2[window_row, window_col]
                        fine_1_window_sum += fine_1[window_row, window_col]
                        fine_2_window_sum += fine_2[window_row, window_col]
                        window_count += 1

                candidate_rows = np.empty(maximum_candidates, dtype=np.int32)
                candidate_cols = np.empty(maximum_candidates, dtype=np.int32)
                candidate_count = 0
                for delta_col in range(-window_radius, window_radius + 1):
                    for delta_row in range(-window_radius, window_radius + 1):
                        bit = (delta_col + window_radius) * side + (
                            delta_row + window_radius
                        )
                        byte = candidate_masks[cell, bit // 8]
                        if byte & np.uint8(1 << (bit % 8)) == 0:
                            continue
                        candidate_row = row + delta_row
                        candidate_col = col + delta_col
                        if not np.isfinite(
                            targets[target_index, candidate_row, candidate_col]
                        ):
                            continue
                        candidate_rows[candidate_count] = candidate_row
                        candidate_cols[candidate_count] = candidate_col
                        candidate_count += 1

                mean_target = target_window_sum / window_count
                mean_coarse_1 = coarse_1_window_sum / window_count
                mean_coarse_2 = coarse_2_window_sum / window_count

                if candidate_count > 5:
                    fine_mean_1 = 0.0
                    fine_mean_2 = 0.0
                    coarse_mean_1 = 0.0
                    coarse_mean_2 = 0.0
                    for candidate in range(candidate_count):
                        candidate_row = candidate_rows[candidate]
                        candidate_col = candidate_cols[candidate]
                        fine_mean_1 += fine_1[candidate_row, candidate_col]
                        fine_mean_2 += fine_2[candidate_row, candidate_col]
                        coarse_mean_1 += coarse_1[candidate_row, candidate_col]
                        coarse_mean_2 += coarse_2[candidate_row, candidate_col]
                    fine_mean_1 /= candidate_count
                    fine_mean_2 /= candidate_count
                    coarse_mean_1 /= candidate_count
                    coarse_mean_2 /= candidate_count

                    bias_1 = -fine_mean_1 + coarse_mean_1 if zero_bias else 0.0
                    bias_2 = -fine_mean_2 + coarse_mean_2 if zero_bias else 0.0
                    fine_pixel_1 = fine_1[row, col] + bias_1
                    fine_pixel_2 = fine_2[row, col] + bias_2

                    inverse_distance = np.empty(maximum_candidates, dtype=np.float64)
                    inverse_distance_sum = 0.0
                    for candidate in range(candidate_count):
                        candidate_row = candidate_rows[candidate]
                        candidate_col = candidate_cols[candidate]
                        fine_candidate_1 = fine_1[candidate_row, candidate_col] + bias_1
                        fine_candidate_2 = fine_2[candidate_row, candidate_col] + bias_2
                        coarse_candidate_1 = coarse_1[candidate_row, candidate_col]
                        coarse_candidate_2 = coarse_2[candidate_row, candidate_col]
                        spectral = 1.0 - 0.5 * (
                            abs(
                                (fine_candidate_1 - coarse_candidate_1)
                                / (fine_candidate_1 + coarse_candidate_1)
                            )
                            + abs(
                                (fine_candidate_2 - coarse_candidate_2)
                                / (fine_candidate_2 + coarse_candidate_2)
                            )
                        )
                        if spectral > 1.0 or spectral < -1.0:
                            spectral = 0.5
                        spatial = 1.0 + math.sqrt(
                            (col - candidate_col) ** 2 + (row - candidate_row) ** 2
                        ) / window_radius
                        combined = (1.0 - spectral) * spatial + 1e-7
                        inverse_distance[candidate] = 1.0 / combined
                        inverse_distance_sum += inverse_distance[candidate]

                    difference_1 = abs(mean_target - mean_coarse_1) + 1e-10
                    difference_2 = abs(mean_target - mean_coarse_2) + 1e-10
                    inverse_temporal_1 = 1.0 / difference_1
                    inverse_temporal_2 = 1.0 / difference_2
                    temporal_1 = inverse_temporal_1 / (
                        inverse_temporal_1 + inverse_temporal_2
                    )
                    temporal_2 = 1.0 - temporal_1

                    change_1 = 0.0
                    change_2 = 0.0
                    fallback_1 = 0.0
                    fallback_2 = 0.0
                    for candidate in range(candidate_count):
                        candidate_row = candidate_rows[candidate]
                        candidate_col = candidate_cols[candidate]
                        weight = inverse_distance[candidate] / inverse_distance_sum
                        target_candidate = targets[
                            target_index, candidate_row, candidate_col
                        ]
                        change_1 += weight * (
                            target_candidate - coarse_1[candidate_row, candidate_col]
                        )
                        change_2 += weight * (
                            target_candidate - coarse_2[candidate_row, candidate_col]
                        )
                        fallback_1 += weight * (
                            fine_1[candidate_row, candidate_col] + bias_1
                        )
                        fallback_2 += weight * (
                            fine_2[candidate_row, candidate_col] + bias_2
                        )

                    prediction = temporal_1 * (fine_pixel_1 + change_1) + temporal_2 * (
                        fine_pixel_2 + change_2
                    )
                    if prediction <= value_min or prediction >= value_max:
                        prediction = temporal_1 * fallback_1 + temporal_2 * fallback_2
                else:
                    fine_pixel_1 = fine_1[row, col]
                    fine_pixel_2 = fine_2[row, col]
                    if zero_bias:
                        fine_pixel_1 = (
                            fine_pixel_1
                            - fine_1_window_sum / window_count
                            + mean_coarse_1
                        )
                        fine_pixel_2 = (
                            fine_pixel_2
                            - fine_2_window_sum / window_count
                            + mean_coarse_2
                        )
                    difference_1 = mean_target - mean_coarse_1 + 1e-10
                    difference_2 = mean_target - mean_coarse_2 + 1e-10
                    inverse_temporal_1 = 1.0 / abs(difference_1)
                    inverse_temporal_2 = 1.0 / abs(difference_2)
                    temporal_1 = inverse_temporal_1 / (
                        inverse_temporal_1 + inverse_temporal_2
                    )
                    prediction = temporal_1 * (
                        fine_pixel_1 + difference_1
                    ) + (1.0 - temporal_1) * (fine_pixel_2 + difference_2)

                output[target_index, row - row_start, col - col_start] = prediction

    return output
