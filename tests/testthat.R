# Script: testthat.R
# Objective: Run the ubESTARFM R package test suite.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Installed ubESTARFM package and tests under tests/testthat.
# Outputs: testthat test results.
# Usage: Rscript tests/testthat.R
# Dependencies: R packages ubESTARFM and testthat.

library(testthat)
library(ubESTARFM)

test_check("ubESTARFM")
