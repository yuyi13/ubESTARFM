# Script: package.R
# Objective: Register compiled code and package-level imports for ubESTARFM.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Compiled Rcpp symbols.
# Outputs: Package namespace registration.
# Usage: Internal package module; loaded through the ubESTARFM namespace.
# Dependencies: R package Rcpp.

#' ubESTARFM package
#'
#' @useDynLib ubESTARFM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @keywords internal
"_PACKAGE"
