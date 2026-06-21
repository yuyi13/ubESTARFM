# Script: test-model-persistence.R
# Objective: Verify R model save and load round trips.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Small reusable models and temporary RDS paths.
# Outputs: testthat assertions for model persistence.
# Usage: Rscript -e 'library(ubestarfm); testthat::test_file("tests/testthat/test-model-persistence.R")'
# Dependencies: R packages testthat and ubestarfm.

test_that("models survive compressed and uncompressed round trips", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )

  for (compress in c(TRUE, FALSE)) {
    path <- tempfile(fileext = ".rds")
    ubestarfm_save_model(model, path, compress = compress)
    loaded <- ubestarfm_load_model(path)
    expect_identical(loaded$candidate_masks, model$candidate_masks)
    expect_equal(loaded$reference_values, model$reference_values)
  }
})

test_that("unsupported model schemas are rejected on load", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm:::ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  model$schema_version <- 99L
  path <- tempfile(fileext = ".rds")
  saveRDS(model, path)

  expect_error(ubestarfm_load_model(path), "Unsupported")
})
