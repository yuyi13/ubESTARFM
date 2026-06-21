# Script: io.R
# Objective: Read, validate, and write single-layer rasters for ubESTARFM.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Raster paths, SpatRaster objects, matrices, and raster templates.
# Outputs: Validated raster data and SpatRaster prediction objects.
# Usage: Internal package module; loaded through the ubESTARFM namespace.
# Dependencies: R package terra.

#' Read one raster layer
#'
#' @param x A file path or a `terra::SpatRaster`.
#' @param name Name used in validation messages.
#'
#' @return A list containing the source raster and a row-column matrix.
#' @export
ubestarfm_read_raster <- function(x, name = deparse(substitute(x))) {
  raster <- if (inherits(x, "SpatRaster")) {
    x
  } else {
    terra::rast(x)
  }

  if (terra::nlyr(raster) != 1L) {
    stop(name, " must contain exactly one raster layer.", call. = FALSE)
  }

  values <- terra::values(raster, mat = FALSE)
  values[!is.finite(values)] <- NA_real_
  matrix_values <- matrix(
    values,
    nrow = terra::nrow(raster),
    ncol = terra::ncol(raster),
    byrow = TRUE
  )

  list(raster = raster, values = matrix_values)
}

#' Validate raster geometry
#'
#' @param rasters A list of `terra::SpatRaster` objects.
#'
#' @return `TRUE`, invisibly, when all geometries match.
#' @export
ubestarfm_validate_geometry <- function(rasters) {
  if (length(rasters) < 2L) {
    return(invisible(TRUE))
  }

  reference <- rasters[[1L]]
  compatible <- vapply(
    rasters[-1L],
    function(candidate) {
      terra::compareGeom(
        reference,
        candidate,
        lyrs = FALSE,
        crs = TRUE,
        ext = TRUE,
        rowcol = TRUE,
        res = TRUE,
        stopOnError = FALSE
      )
    },
    logical(1)
  )

  if (!all(compatible)) {
    stop(
      "All rasters must have the same dimensions, extent, resolution, and CRS.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

ubestarfm_matrix_to_raster <- function(values, template) {
  if (!is.matrix(values)) {
    stop("values must be a matrix.", call. = FALSE)
  }
  expected_dimensions <- c(terra::nrow(template), terra::ncol(template))
  if (!all(dim(values) == expected_dimensions)) {
    stop("Matrix dimensions must match the raster template.", call. = FALSE)
  }

  terra::setValues(template, as.vector(t(values)))
}
