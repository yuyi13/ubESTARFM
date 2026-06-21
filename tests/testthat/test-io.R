# Script: test-io.R
# Objective: Verify R raster loading and geometry validation.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled GeoTIFF fixtures and in-memory SpatRaster objects.
# Outputs: testthat assertions for raster I/O behavior.
# Usage: Rscript -e 'testthat::test_file("tests/testthat/test-io.R")'
# Dependencies: R packages testthat, terra, and ubestarfm.

test_that("single-layer GeoTIFFs load as row-column matrices", {
  input <- ubestarfm_read_raster(
    example_path("MOD11A1_LST_cloudrm_20160218.tif")
  )

  expect_s4_class(input$raster, "SpatRaster")
  expect_equal(dim(input$values), c(400L, 400L))
  expect_equal(input$values[1, 1], 312.74, tolerance = 1e-4)
  expect_true(all(is.finite(input$values)))
})

test_that("reference rasters must share geometry", {
  first <- terra::rast(
    nrows = 3,
    ncols = 3,
    xmin = 0,
    xmax = 3,
    ymin = 0,
    ymax = 3
  )
  second <- terra::rast(
    nrows = 4,
    ncols = 3,
    xmin = 0,
    xmax = 3,
    ymin = 0,
    ymax = 4
  )

  expect_error(
    ubestarfm_validate_geometry(list(first, second)),
    "same dimensions, extent, resolution, and CRS"
  )
})
