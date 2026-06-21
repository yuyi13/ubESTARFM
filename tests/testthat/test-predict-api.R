# Script: test-predict-api.R
# Objective: Verify public R single-target and batch prediction APIs.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled GeoTIFF fixtures and reusable models.
# Outputs: testthat assertions for prediction objects and files.
# Usage: Rscript -e 'library(ubestarfm); testthat::test_file("tests/testthat/test-predict-api.R")'
# Dependencies: R packages testthat, terra, and ubestarfm.

test_that("public single and batch predictions preserve geometry", {
  model <- train_example_model(window_radius = 2L, patch_size = 100L)
  target_path <- example_path("MOD11A1_LST_cloudrm_20160218.tif")
  output_path <- tempfile(fileext = ".tif")

  single <- ubestarfm_predict(
    model,
    target_path,
    value_range = c(250, 350),
    output_path = output_path
  )
  batch <- ubestarfm_predict_batch(
    model,
    list(target_path, target_path),
    value_range = c(250, 350)
  )

  expect_s4_class(single, "SpatRaster")
  expect_true(file.exists(output_path))
  expect_true(terra::compareGeom(single, batch[[1L]], stopOnError = FALSE))
  expect_equal(terra::values(batch[[1L]]), terra::values(batch[[2L]]))
})

test_that("prediction refuses to overwrite unless requested", {
  model <- train_example_model(window_radius = 1L, patch_size = 100L)
  output_path <- tempfile(fileext = ".tif")
  file.create(output_path)

  expect_error(
    ubestarfm_predict(
      model,
      example_path("MOD11A1_LST_cloudrm_20160218.tif"),
      output_path = output_path
    ),
    "already exists"
  )
})
