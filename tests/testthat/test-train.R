# Script: test-train.R
# Objective: Verify reusable R model training, parallel determinism, and memory.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled GeoTIFFs and small in-memory arrays.
# Outputs: testthat assertions for the training API.
# Usage: Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-train.R")'
# Dependencies: R packages testthat and ubESTARFM.

test_that("parallel and sequential training are deterministic", {
  fixture <- small_reference_fixture()
  sequential <- do.call(
    ubESTARFM:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 2L, workers = 1L))
  )
  parallel <- do.call(
    ubESTARFM:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 2L, workers = 2L))
  )

  expect_identical(sequential$candidate_masks, parallel$candidate_masks)
  expect_equal(sequential$patch_thresholds, parallel$patch_thresholds)
})

test_that("bundled candidate-cache model stays below 80 MiB", {
  model <- train_example_model(window_radius = 25L, patch_size = 200L)

  expect_lt(as.numeric(object.size(model)), 80 * 1024^2)
  expect_type(model$candidate_masks, "raw")
  expect_equal(length(model$candidate_masks), 400L * 400L * 326L)
})
