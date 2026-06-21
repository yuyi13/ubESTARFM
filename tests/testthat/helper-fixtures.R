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
