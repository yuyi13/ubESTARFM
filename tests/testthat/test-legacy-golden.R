# Script: test-legacy-golden.R
# Objective: Compare modern R predictions with the published tracked output.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled reference rasters, target raster, and legacy fused GeoTIFF.
# Outputs: testthat assertions for numerical and spatial compatibility.
# Usage: Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-legacy-golden.R")'
# Dependencies: R packages testthat, terra, and ubESTARFM.

test_that("modern R matches the published fused result", {
  model <- train_example_model(
    window_radius = 25L,
    patch_size = 200L,
    workers = 4L
  )
  modern <- ubestarfm_predict(
    model,
    example_path("MOD11A1_LST_cloudrm_20160218.tif"),
    value_range = c(250, 350),
    workers = 4L
  )
  legacy <- terra::rast(
    testthat::test_path("..", "fixtures", "legacy", "fused_result.tif")
  )

  expect_true(terra::compareGeom(modern, legacy, stopOnError = FALSE))
  modern_values <- terra::values(modern, mat = FALSE)
  legacy_values <- terra::values(legacy, mat = FALSE)
  expect_identical(is.na(modern_values), is.na(legacy_values))
  finite <- is.finite(modern_values) & is.finite(legacy_values)
  expect_lte(max(abs(modern_values[finite] - legacy_values[finite])), 1e-4)
})
