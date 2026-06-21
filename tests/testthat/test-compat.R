# Script: test-compat.R
# Objective: Verify compatibility with the published ubESTARFM function signature.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled reference and target GeoTIFFs.
# Outputs: testthat assertions for the compatibility wrapper.
# Usage: Rscript -e 'library(ubestarfm); testthat::test_file("tests/testthat/test-compat.R")'
# Dependencies: R packages testthat, terra, and ubestarfm.

test_that("published wrapper matches explicit train and predict", {
  temporary_directory <- tempfile()
  dir.create(temporary_directory)
  unrelated <- file.path(temporary_directory, "unrelated.tif")
  terra::writeRaster(
    terra::rast(nrows = 1, ncols = 1, vals = 1),
    unrelated
  )
  output_path <- file.path(temporary_directory, "compatibility.tif")

  wrapped <- ubESTARFM(
    w = 2L,
    DN_min = 250,
    DN_max = 350,
    patch_long = 100L,
    tmp_path = temporary_directory,
    out_path = output_path,
    method = "zero bias",
    rst_fine1 = example_path("Landsat_LST_cloudrm_20160205.tif"),
    rst_fine2 = example_path("Landsat_LST_cloudrm_20160308.tif"),
    rst_coarse1 = example_path("MOD11A1_LST_cloudrm_20160205.tif"),
    rst_coarse2 = example_path("MOD11A1_LST_cloudrm_20160308.tif"),
    rst_coarse0 = example_path("MOD11A1_LST_cloudrm_20160218.tif")
  )
  model <- train_example_model(window_radius = 2L, patch_size = 100L)
  explicit <- ubestarfm_predict(
    model,
    example_path("MOD11A1_LST_cloudrm_20160218.tif")
  )

  expect_s4_class(wrapped, "SpatRaster")
  expect_true(file.exists(output_path))
  expect_true(file.exists(unrelated))
  expect_equal(terra::values(wrapped), terra::values(explicit))
})

test_that("published wrapper validates method names", {
  expect_error(
    ubESTARFM(
      DN_min = 250,
      DN_max = 350,
      out_path = tempfile(fileext = ".tif"),
      method = "unknown",
      rst_fine1 = example_path("Landsat_LST_cloudrm_20160205.tif"),
      rst_fine2 = example_path("Landsat_LST_cloudrm_20160308.tif"),
      rst_coarse1 = example_path("MOD11A1_LST_cloudrm_20160205.tif"),
      rst_coarse2 = example_path("MOD11A1_LST_cloudrm_20160308.tif"),
      rst_coarse0 = example_path("MOD11A1_LST_cloudrm_20160218.tif")
    ),
    "zero bias.*baseline"
  )
})
