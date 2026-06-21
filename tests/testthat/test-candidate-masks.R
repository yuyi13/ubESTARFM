# Script: test-candidate-masks.R
# Objective: Verify packed candidate membership and ordering in R training.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Small in-memory reference arrays.
# Outputs: testthat assertions for candidate masks.
# Usage: Rscript -e 'library(ubestarfm); testthat::test_file("tests/testthat/test-candidate-masks.R")'
# Dependencies: R packages testthat and ubestarfm.

test_that("candidate bits use column-major local-window order", {
  fixture <- small_reference_fixture()
  result <- ubestarfm:::ubestarfm_train_arrays(
    fixture$fine_1,
    fixture$fine_2,
    fixture$coarse_1,
    fixture$coarse_2,
    window_radius = 1L,
    patch_size = 5L,
    method = "zero_bias",
    workers = 1L
  )

  center_bits <- ubestarfm:::ubestarfm_unpack_candidates(
    result,
    row = 3L,
    col = 3L
  )
  expect_equal(center_bits$row, c(3L, 4L, 2L, 3L, 4L, 2L, 3L, 4L))
  expect_equal(center_bits$col, c(2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L))
})

test_that("candidate membership uses strict thresholds and valid references", {
  fine_1 <- matrix(c(0, 1, 2, 3), nrow = 2L, byrow = TRUE)
  fine_2 <- fine_1
  coarse_1 <- fine_1 + 10
  coarse_2 <- fine_1 + 20
  coarse_1[1, 2] <- NA_real_

  model <- ubestarfm:::ubestarfm_train_arrays(
    fine_1,
    fine_2,
    coarse_1,
    coarse_2,
    window_radius = 1L,
    patch_size = 2L
  )
  candidates <- ubestarfm:::ubestarfm_unpack_candidates(
    model,
    row = 1L,
    col = 1L
  )

  expect_false(any(candidates$row == 1L & candidates$col == 2L))
  expect_true(all(candidates$row >= 1L & candidates$col >= 1L))
})

test_that("candidate masks do not depend on correction method", {
  fixture <- small_reference_fixture()
  zero_bias <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L, method = "zero_bias"))
  )
  baseline <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L, method = "baseline"))
  )

  expect_identical(zero_bias$candidate_masks, baseline$candidate_masks)
})

test_that("training is independent of target dates", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )

  expect_false("coarse_target" %in% names(model))
})
