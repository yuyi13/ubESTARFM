# Script: helper-fixtures.R
# Objective: Provide shared raster paths and small numerical fixtures for R tests.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Installed package example data.
# Outputs: Helper functions and in-memory test fixtures.
# Usage: Loaded automatically by testthat.
# Dependencies: R package ubESTARFM.

example_path <- function(filename) {
  installed_path <- system.file(
    "extdata",
    filename,
    package = "ubESTARFM",
    mustWork = FALSE
  )
  if (nzchar(installed_path)) {
    return(installed_path)
  }
  file.path("inst", "extdata", filename)
}

small_reference_fixture <- function() {
  fine_1   <- matrix(seq_len(25), nrow = 5L, byrow = TRUE)
  fine_2   <- fine_1 + 100
  coarse_1 <- fine_1 + 200
  coarse_2 <- fine_1 + 300
  coarse_1[2, 2] <- NA_real_
  list(
    fine_1   = fine_1,
    fine_2   = fine_2,
    coarse_1 = coarse_1,
    coarse_2 = coarse_2
  )
}

train_example_model <- function(
  window_radius = 25L,
  patch_size = 200L,
  method = "zero_bias",
  workers = 1L
) {
  ubestarfm_train(
    fine_1 = example_path("Landsat_LST_cloudrm_20160205.tif"),
    fine_2 = example_path("Landsat_LST_cloudrm_20160308.tif"),
    coarse_1 = example_path("MOD11A1_LST_cloudrm_20160205.tif"),
    coarse_2 = example_path("MOD11A1_LST_cloudrm_20160308.tif"),
    window_radius = window_radius,
    patch_size = patch_size,
    method = method,
    workers = workers
  )
}

reference_predict_pixel <- function(model, coarse_target, row, col, value_range) {
  radius <- model$parameters$window_radius
  nrow   <- model$metadata$nrow
  ncol   <- model$metadata$ncol
  rows   <- max(1L, row - radius):min(nrow, row + radius)
  cols   <- max(1L, col - radius):min(ncol, col + radius)

  fine_1   <- model$reference_values$fine_1
  fine_2   <- model$reference_values$fine_2
  coarse_1 <- model$reference_values$coarse_1
  coarse_2 <- model$reference_values$coarse_2
  target_valid <- model$reference_valid & is.finite(coarse_target)
  window_valid <- target_valid[rows, cols, drop = FALSE]

  candidates <- ubESTARFM:::ubestarfm_unpack_candidates(
    model,
    row = row,
    col = col
  )
  keep <- is.finite(coarse_target[cbind(candidates$row, candidates$col)])
  candidates <- candidates[keep, , drop = FALSE]

  target_window   <- coarse_target[rows, cols, drop = FALSE][window_valid]
  coarse_1_window <- coarse_1[rows, cols, drop = FALSE][window_valid]
  coarse_2_window <- coarse_2[rows, cols, drop = FALSE][window_valid]

  if (nrow(candidates) > 5L) {
    index <- cbind(candidates$row, candidates$col)
    fine_candidates_1   <- fine_1[index]
    fine_candidates_2   <- fine_2[index]
    coarse_candidates_1 <- coarse_1[index]
    coarse_candidates_2 <- coarse_2[index]

    if (model$parameters$method == "zero_bias") {
      bias_1 <- -mean(fine_candidates_1) + mean(coarse_candidates_1)
      bias_2 <- -mean(fine_candidates_2) + mean(coarse_candidates_2)
      fine_candidates_1 <- fine_candidates_1 + bias_1
      fine_candidates_2 <- fine_candidates_2 + bias_2
      fine_pixel_1 <- fine_1[row, col] + bias_1
      fine_pixel_2 <- fine_2[row, col] + bias_2
    } else {
      fine_pixel_1 <- fine_1[row, col]
      fine_pixel_2 <- fine_2[row, col]
    }

    spectral_distance <- 1 - 0.5 * (
      abs(
        (fine_candidates_1 - coarse_candidates_1) /
          (fine_candidates_1 + coarse_candidates_1)
      ) +
        abs(
          (fine_candidates_2 - coarse_candidates_2) /
            (fine_candidates_2 + coarse_candidates_2)
        )
    )
    abnormal <- spectral_distance > 1 | spectral_distance < -1
    spectral_distance[abnormal] <- 0.5
    spatial_distance <- 1 + sqrt(
      (col - candidates$col)^2 + (row - candidates$row)^2
    ) / radius
    combined_distance <- (1 - spectral_distance) * spatial_distance + 1e-7
    weight <- (1 / combined_distance) / sum(1 / combined_distance)

    temporal_difference_1 <- abs(mean(target_window - coarse_1_window)) + 1e-10
    temporal_difference_2 <- abs(mean(target_window - coarse_2_window)) + 1e-10
    temporal_weight_1 <- (1 / temporal_difference_1) /
      (1 / temporal_difference_1 + 1 / temporal_difference_2)
    temporal_weight_2 <- 1 - temporal_weight_1

    target_candidates <- coarse_target[index]
    prediction_1 <- fine_pixel_1 + sum(
      weight * (target_candidates - coarse_candidates_1)
    )
    prediction_2 <- fine_pixel_2 + sum(
      weight * (target_candidates - coarse_candidates_2)
    )
    prediction <- temporal_weight_1 * prediction_1 +
      temporal_weight_2 * prediction_2

    if (prediction <= value_range[1] || prediction >= value_range[2]) {
      prediction <- temporal_weight_1 * sum(weight * fine_candidates_1) +
        temporal_weight_2 * sum(weight * fine_candidates_2)
    }
    return(prediction)
  }

  if (model$parameters$method == "zero_bias") {
    fine_pixel_1 <- fine_1[row, col] -
      mean(fine_1[rows, cols, drop = FALSE][window_valid]) +
      mean(coarse_1_window)
    fine_pixel_2 <- fine_2[row, col] -
      mean(fine_2[rows, cols, drop = FALSE][window_valid]) +
      mean(coarse_2_window)
  } else {
    fine_pixel_1 <- fine_1[row, col]
    fine_pixel_2 <- fine_2[row, col]
  }

  temporal_difference_1 <- mean(target_window - coarse_1_window) + 1e-10
  temporal_difference_2 <- mean(target_window - coarse_2_window) + 1e-10
  temporal_weight_1 <- (1 / abs(temporal_difference_1)) /
    (1 / abs(temporal_difference_1) + 1 / abs(temporal_difference_2))
  temporal_weight_2 <- 1 - temporal_weight_1
  temporal_weight_1 * (fine_pixel_1 + temporal_difference_1) +
    temporal_weight_2 * (fine_pixel_2 + temporal_difference_2)
}
