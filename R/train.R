# Script: train.R
# Objective: Train compact reusable ubESTARFM candidate-membership models.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Two fine and two coarse reference rasters or matrices.
# Outputs: Versioned ubestarfm_model objects.
# Usage: Internal package module; loaded through the ubestarfm namespace.
# Dependencies: R packages terra and Rcpp plus compiled ubESTARFM kernels.

ubestarfm_make_patches <- function(nrow, ncol, patch_size) {
  row_starts <- seq.int(1L, nrow, by = patch_size)
  col_starts <- seq.int(1L, ncol, by = patch_size)
  patches <- vector("list", length(row_starts) * length(col_starts))

  patch_id <- 1L
  for (row_start in row_starts) {
    for (col_start in col_starts) {
      patches[[patch_id]] <- list(
        id        = patch_id,
        row_start = row_start,
        row_end   = min(nrow, row_start + patch_size - 1L),
        col_start = col_start,
        col_end   = min(ncol, col_start + patch_size - 1L)
      )
      patch_id <- patch_id + 1L
    }
  }

  patches
}

ubestarfm_validate_training_parameters <- function(
  window_radius,
  patch_size,
  method,
  workers
) {
  window_radius <- as.integer(window_radius)
  patch_size    <- as.integer(patch_size)
  workers       <- as.integer(workers)
  method        <- match.arg(method, c("zero_bias", "baseline"))

  if (length(window_radius) != 1L || is.na(window_radius) || window_radius < 1L) {
    stop("window_radius must be a positive integer.", call. = FALSE)
  }
  if (length(patch_size) != 1L || is.na(patch_size) || patch_size < 1L) {
    stop("patch_size must be a positive integer.", call. = FALSE)
  }
  if (length(workers) != 1L || is.na(workers) || workers < 1L) {
    stop("workers must be a positive integer.", call. = FALSE)
  }

  list(
    window_radius = window_radius,
    patch_size    = patch_size,
    method        = method,
    workers       = workers
  )
}

ubestarfm_train_arrays <- function(
  fine_1,
  fine_2,
  coarse_1,
  coarse_2,
  window_radius = 25L,
  patch_size = 200L,
  method = "zero_bias",
  workers = 1L,
  metadata = NULL
) {
  references <- list(
    fine_1   = fine_1,
    fine_2   = fine_2,
    coarse_1 = coarse_1,
    coarse_2 = coarse_2
  )
  if (!all(vapply(references, is.matrix, logical(1)))) {
    stop("All reference inputs must be matrices.", call. = FALSE)
  }
  dimensions <- lapply(references, dim)
  if (!all(vapply(dimensions, identical, logical(1), dimensions[[1L]]))) {
    stop("All reference matrices must have identical dimensions.", call. = FALSE)
  }

  parameters <- ubestarfm_validate_training_parameters(
    window_radius,
    patch_size,
    method,
    workers
  )
  nrow <- nrow(fine_1)
  ncol <- ncol(fine_1)
  side <- 2L * parameters$window_radius + 1L
  mask_bytes <- as.integer(ceiling(side^2 / 8))

  references <- lapply(references, function(values) {
    storage.mode(values) <- "double"
    values[!is.finite(values)] <- NA_real_
    values
  })
  reference_valid <- Reduce(
    `&`,
    lapply(references, is.finite)
  )
  if (!any(reference_valid)) {
    stop("The reference rasters contain no jointly valid pixels.", call. = FALSE)
  }

  patches <- ubestarfm_make_patches(nrow, ncol, parameters$patch_size)
  patch_thresholds <- matrix(
    NA_real_,
    nrow = length(patches),
    ncol = 2L,
    dimnames = list(NULL, c("fine_1", "fine_2"))
  )
  patch_ids <- matrix(NA_integer_, nrow = nrow, ncol = ncol)

  for (patch in patches) {
    rows <- patch$row_start:patch$row_end
    cols <- patch$col_start:patch$col_end
    patch_thresholds[patch$id, ] <- c(
      stats::sd(references$fine_1[rows, cols], na.rm = TRUE),
      stats::sd(references$fine_2[rows, cols], na.rm = TRUE)
    )
    patch_ids[rows, cols] <- patch$id
  }

  worker <- function(patch) {
    thresholds <- patch_thresholds[patch$id, ]
    ubestarfm_train_patch_cpp(
      references$fine_1,
      references$fine_2,
      references$coarse_1,
      references$coarse_2,
      patch$row_start - 1L,
      patch$row_end - 1L,
      patch$col_start - 1L,
      patch$col_end - 1L,
      parameters$window_radius,
      thresholds[1L],
      thresholds[2L]
    )
  }
  patch_results <- ubestarfm_lapply(
    patches,
    worker,
    workers = parameters$workers
  )
  candidate_masks <- ubestarfm_assemble_masks_cpp(
    patch_results,
    nrow * ncol,
    mask_bytes
  )

  if (is.null(metadata)) {
    metadata <- list(nrow = nrow, ncol = ncol)
  } else {
    metadata$nrow <- nrow
    metadata$ncol <- ncol
  }
  model <- new_ubestarfm_model(
    reference_values = references,
    reference_valid = reference_valid,
    candidate_masks = candidate_masks,
    mask_bytes = mask_bytes,
    patch_ids = patch_ids,
    patch_thresholds = patch_thresholds,
    metadata = metadata,
    parameters = list(
      window_radius = parameters$window_radius,
      patch_size    = parameters$patch_size,
      method        = parameters$method,
      mask_side     = side,
      mask_bytes    = mask_bytes
    )
  )
  ubestarfm_validate_model(model)
  model
}

#' Train a reusable ubESTARFM reference-pair model
#'
#' @param fine_1 Fine-resolution raster for the first reference date.
#' @param fine_2 Fine-resolution raster for the second reference date.
#' @param coarse_1 Coarse-resolution raster for the first reference date.
#' @param coarse_2 Coarse-resolution raster for the second reference date.
#' @param window_radius Moving-window radius in pixels.
#' @param patch_size Patch side length in pixels.
#' @param method Either `"zero_bias"` or `"baseline"`.
#' @param cache Cache mode. Reusable training supports `"candidates"`.
#' @param workers Number of patch workers.
#'
#' @return An object of class `ubestarfm_model`.
#' @export
ubestarfm_train <- function(
  fine_1,
  fine_2,
  coarse_1,
  coarse_2,
  window_radius = 25L,
  patch_size = 200L,
  method = "zero_bias",
  cache = "candidates",
  workers = 1L
) {
  if (!identical(cache, "candidates")) {
    stop("ubestarfm_train() supports cache = \"candidates\".", call. = FALSE)
  }

  inputs <- list(
    fine_1   = ubestarfm_read_raster(fine_1, "fine_1"),
    fine_2   = ubestarfm_read_raster(fine_2, "fine_2"),
    coarse_1 = ubestarfm_read_raster(coarse_1, "coarse_1"),
    coarse_2 = ubestarfm_read_raster(coarse_2, "coarse_2")
  )
  ubestarfm_validate_geometry(lapply(inputs, `[[`, "raster"))

  template <- inputs$fine_1$raster
  extent <- as.vector(terra::ext(template))
  metadata <- list(
    nrow       = terra::nrow(template),
    ncol       = terra::ncol(template),
    extent     = unname(extent),
    resolution = unname(terra::res(template)),
    crs        = terra::crs(template, proj = TRUE),
    nodata     = terra::NAflag(template)
  )

  ubestarfm_train_arrays(
    inputs$fine_1$values,
    inputs$fine_2$values,
    inputs$coarse_1$values,
    inputs$coarse_2$values,
    window_radius = window_radius,
    patch_size = patch_size,
    method = method,
    workers = workers,
    metadata = metadata
  )
}

ubestarfm_unpack_candidates <- function(model, row, col) {
  ubestarfm_validate_model(model)
  row <- as.integer(row)
  col <- as.integer(col)
  radius <- model$parameters$window_radius
  side <- model$parameters$mask_side
  nrow <- model$metadata$nrow
  ncol <- model$metadata$ncol

  if (row < 1L || row > nrow || col < 1L || col > ncol) {
    stop("row and col must identify a model cell.", call. = FALSE)
  }

  cell <- (row - 1L) * ncol + (col - 1L)
  base <- cell * model$mask_bytes
  candidate_rows <- integer()
  candidate_cols <- integer()

  for (bit in 0:(side^2 - 1L)) {
    byte <- as.integer(model$candidate_masks[base + bit %/% 8L + 1L])
    is_set <- bitwAnd(byte, bitwShiftL(1L, bit %% 8L)) != 0L
    if (!is_set) {
      next
    }
    delta_col <- bit %/% side - radius
    delta_row <- bit %% side - radius
    candidate_rows <- c(candidate_rows, row + delta_row)
    candidate_cols <- c(candidate_cols, col + delta_col)
  }

  data.frame(row = candidate_rows, col = candidate_cols)
}
