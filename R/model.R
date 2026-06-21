# Script: model.R
# Objective: Define and validate versioned reusable ubESTARFM model objects.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Reference arrays, candidate masks, spatial metadata, and parameters.
# Outputs: Validated objects of class ubestarfm_model.
# Usage: Internal package module; loaded through the ubestarfm namespace.
# Dependencies: Base R and utils.

new_ubestarfm_model <- function(
  reference_values,
  reference_valid,
  candidate_masks,
  mask_bytes,
  patch_ids,
  patch_thresholds,
  metadata,
  parameters
) {
  structure(
    list(
      schema_version    = 1L,
      package_version   = as.character(utils::packageVersion("ubestarfm")),
      reference_values  = reference_values,
      reference_valid   = reference_valid,
      candidate_masks   = candidate_masks,
      mask_bytes        = as.integer(mask_bytes),
      patch_ids         = patch_ids,
      patch_thresholds  = patch_thresholds,
      metadata          = metadata,
      parameters        = parameters
    ),
    class = "ubestarfm_model"
  )
}

#' Validate a reusable ubESTARFM model
#'
#' @param model An object expected to inherit from `ubestarfm_model`.
#'
#' @return The model, invisibly, when valid.
#' @export
ubestarfm_validate_model <- function(model) {
  if (!inherits(model, "ubestarfm_model")) {
    stop("model must inherit from ubestarfm_model.", call. = FALSE)
  }
  if (!identical(model$schema_version, 1L)) {
    stop("Unsupported ubestarfm_model schema version.", call. = FALSE)
  }

  required_references <- c("fine_1", "fine_2", "coarse_1", "coarse_2")
  if (
    !is.list(model$reference_values) ||
      !all(required_references %in% names(model$reference_values))
  ) {
    stop("model must contain four named reference arrays.", call. = FALSE)
  }

  dimensions <- vapply(model$reference_values[required_references], length, integer(1))
  if (length(unique(dimensions)) != 1L) {
    stop("All reference arrays must have equal lengths.", call. = FALSE)
  }

  if (!is.logical(model$reference_valid) || length(model$reference_valid) != dimensions[1]) {
    stop("reference_valid must match the reference arrays.", call. = FALSE)
  }
  if (!is.raw(model$candidate_masks)) {
    stop("candidate_masks must be stored as a raw vector.", call. = FALSE)
  }
  if (!is.numeric(model$mask_bytes) || length(model$mask_bytes) != 1L) {
    stop("mask_bytes must be a positive scalar.", call. = FALSE)
  }

  required_metadata <- c("nrow", "ncol")
  if (
    !is.list(model$metadata) ||
      !all(required_metadata %in% names(model$metadata))
  ) {
    stop("model metadata must contain nrow and ncol.", call. = FALSE)
  }
  if (model$metadata$nrow * model$metadata$ncol != dimensions[1]) {
    stop("model dimensions do not match its reference arrays.", call. = FALSE)
  }
  expected_mask_length <- dimensions[1] * as.integer(model$mask_bytes)
  if (length(model$candidate_masks) != expected_mask_length) {
    stop("candidate_masks length does not match model dimensions.", call. = FALSE)
  }

  if (
    !is.list(model$parameters) ||
      !model$parameters$method %in% c("zero_bias", "baseline")
  ) {
    stop("model parameters contain an unsupported method.", call. = FALSE)
  }

  invisible(model)
}

#' Save a reusable ubESTARFM model
#'
#' @param model A reusable `ubestarfm_model`.
#' @param path Destination RDS path.
#' @param compress Compression passed to [saveRDS()].
#'
#' @return The normalized output path, invisibly.
#' @export
ubestarfm_save_model <- function(model, path, compress = TRUE) {
  ubestarfm_validate_model(model)
  saveRDS(model, path, compress = compress)
  invisible(normalizePath(path, mustWork = TRUE))
}

#' Load a reusable ubESTARFM model
#'
#' @param path An RDS model path.
#'
#' @return A validated `ubestarfm_model`.
#' @export
ubestarfm_load_model <- function(path) {
  model <- readRDS(path)
  ubestarfm_validate_model(model)
  model
}
