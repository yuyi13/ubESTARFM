# Script: package.R
# Objective: Register compiled code and package-level imports for ubestarfm.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Compiled Rcpp symbols.
# Outputs: Package namespace registration.
# Usage: Internal package module; loaded through the ubestarfm namespace.
# Dependencies: R package Rcpp.

#' ubestarfm package
#'
#' @useDynLib ubestarfm, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @keywords internal
"_PACKAGE"
