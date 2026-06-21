# Script: test-model.R
# Objective: Verify validation of versioned R model objects.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: In-memory valid and malformed model objects.
# Outputs: testthat assertions for model validation.
# Usage: Rscript -e 'testthat::test_file("tests/testthat/test-model.R")'
# Dependencies: R packages testthat and ubESTARFM.

test_that("model validation rejects malformed models", {
  expect_error(
    ubestarfm_validate_model(list(schema_version = 1L)),
    "ubestarfm_model"
  )
})
