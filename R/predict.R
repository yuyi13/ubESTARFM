# Script: predict.R
# Objective: Predict one or more target dates from a reusable ubESTARFM model.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: A reusable model, coarse target rasters, ranges, and output paths.
# Outputs: Fine-resolution prediction matrices and SpatRaster objects.
# Usage: Internal package module; loaded through the ubESTARFM namespace.
# Dependencies: R package terra and compiled ubESTARFM kernels.

ubestarfm_validate_value_range <- function(value_range) {
  if (
    !is.numeric(value_range) ||
      length(value_range) != 2L ||
      anyNA(value_range) ||
      value_range[1L] >= value_range[2L]
  ) {
    stop(
      "value_range must contain an increasing lower and upper bound.",
      call. = FALSE
    )
  }
  as.numeric(value_range)
}

ubestarfm_predict_arrays <- function(
  model,
  coarse_targets,
  value_range = c(250, 350),
  workers = 1L
) {
  ubestarfm_validate_model(model)
  value_range <- ubestarfm_validate_value_range(value_range)
  if (!is.list(coarse_targets) || length(coarse_targets) < 1L) {
    stop("coarse_targets must be a non-empty list of matrices.", call. = FALSE)
  }
  expected_dimensions <- c(model$metadata$nrow, model$metadata$ncol)
  if (!all(vapply(
    coarse_targets,
    function(values) is.matrix(values) && identical(dim(values), expected_dimensions),
    logical(1)
  ))) {
    stop("All target matrices must match the model dimensions.", call. = FALSE)
  }

  coarse_targets <- lapply(coarse_targets, function(values) {
    storage.mode(values) <- "double"
    values[!is.finite(values)] <- NA_real_
    values
  })
  target_values <- unlist(
    lapply(coarse_targets, function(values) as.vector(t(values))),
    use.names = FALSE
  )

  patches <- ubestarfm_make_patches(
    model$metadata$nrow,
    model$metadata$ncol,
    model$parameters$patch_size
  )
  worker <- function(patch) {
    ubestarfm_predict_patch_cpp(
      model$reference_values$fine_1,
      model$reference_values$fine_2,
      model$reference_values$coarse_1,
      model$reference_values$coarse_2,
      model$candidate_masks,
      model$mask_bytes,
      target_values,
      length(coarse_targets),
      patch$row_start - 1L,
      patch$row_end - 1L,
      patch$col_start - 1L,
      patch$col_end - 1L,
      model$parameters$window_radius,
      model$parameters$method,
      value_range[1L],
      value_range[2L]
    )
  }
  patch_results <- ubestarfm_lapply(patches, worker, workers = workers)

  predictions <- lapply(
    seq_along(coarse_targets),
    function(target_index) {
      matrix(
        NA_real_,
        nrow = model$metadata$nrow,
        ncol = model$metadata$ncol
      )
    }
  )
  for (patch_index in seq_along(patches)) {
    patch <- patches[[patch_index]]
    result <- patch_results[[patch_index]]
    patch_nrow <- patch$row_end - patch$row_start + 1L
    patch_ncol <- patch$col_end - patch$col_start + 1L
    rows <- patch$row_start:patch$row_end
    cols <- patch$col_start:patch$col_end
    for (target_index in seq_along(predictions)) {
      predictions[[target_index]][rows, cols] <- matrix(
        result[, target_index],
        nrow = patch_nrow,
        ncol = patch_ncol,
        byrow = TRUE
      )
    }
  }

  predictions
}

ubestarfm_template_from_model <- function(model) {
  extent <- model$metadata$extent
  terra::rast(
    nrows = model$metadata$nrow,
    ncols = model$metadata$ncol,
    xmin = extent[1L],
    xmax = extent[2L],
    ymin = extent[3L],
    ymax = extent[4L],
    crs = model$metadata$crs
  )
}

ubestarfm_write_prediction <- function(
  prediction,
  output_path,
  overwrite,
  nodata
) {
  output_path <- normalizePath(
    output_path,
    winslash = "/",
    mustWork = FALSE
  )
  if (file.exists(output_path) && !overwrite) {
    stop("Output already exists: ", output_path, call. = FALSE)
  }

  output_directory <- dirname(output_path)
  if (!dir.exists(output_directory)) {
    created <- dir.create(output_directory, recursive = TRUE)
    if (!created) {
      stop("Could not create output directory: ", output_directory, call. = FALSE)
    }
  }

  temporary_path <- tempfile(
    pattern = paste0(".", basename(output_path), "-"),
    tmpdir = output_directory,
    fileext = ".tif"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  na_flag <- if (is.numeric(nodata) && length(nodata) == 1L && is.finite(nodata)) {
    nodata
  } else {
    -3.4e38
  }
  terra::writeRaster(
    prediction,
    temporary_path,
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = na_flag
  )

  if (file.exists(output_path)) {
    unlink(output_path)
  }
  if (!file.rename(temporary_path, output_path)) {
    stop("Could not move completed prediction to: ", output_path, call. = FALSE)
  }
  invisible(output_path)
}

#' Predict multiple target dates
#'
#' @param model A reusable `ubestarfm_model`.
#' @param coarse_targets A list of target raster paths or `SpatRaster` objects.
#' @param output_paths Optional output paths, one per target.
#' @param value_range Allowed output range before the published fallback.
#' @param overwrite Whether existing output files may be replaced.
#' @param workers Number of patch workers.
#'
#' @return A list of predicted `terra::SpatRaster` objects.
#' @export
ubestarfm_predict_batch <- function(
  model,
  coarse_targets,
  output_paths = NULL,
  value_range = c(250, 350),
  overwrite = FALSE,
  workers = 1L
) {
  ubestarfm_validate_model(model)
  if (!is.list(coarse_targets)) {
    coarse_targets <- list(coarse_targets)
  }
  if (length(coarse_targets) < 1L) {
    stop("coarse_targets must contain at least one raster.", call. = FALSE)
  }

  inputs <- lapply(
    seq_along(coarse_targets),
    function(index) {
      ubestarfm_read_raster(
        coarse_targets[[index]],
        paste0("coarse_targets[[", index, "]]")
      )
    }
  )
  template <- ubestarfm_template_from_model(model)
  ubestarfm_validate_geometry(c(
    list(template),
    lapply(inputs, `[[`, "raster")
  ))

  if (!is.null(output_paths) && length(output_paths) != length(inputs)) {
    stop("output_paths must have one entry per target.", call. = FALSE)
  }
  matrices <- ubestarfm_predict_arrays(
    model,
    lapply(inputs, `[[`, "values"),
    value_range = value_range,
    workers = workers
  )
  predictions <- lapply(
    matrices,
    ubestarfm_matrix_to_raster,
    template = template
  )

  if (!is.null(output_paths)) {
    for (index in seq_along(predictions)) {
      ubestarfm_write_prediction(
        predictions[[index]],
        output_paths[[index]],
        overwrite,
        model$metadata$nodata
      )
    }
  }

  predictions
}

#' Predict one target date
#'
#' @param model A reusable `ubestarfm_model`.
#' @param coarse_target A target raster path or `SpatRaster`.
#' @param value_range Allowed output range before the published fallback.
#' @param output_path Optional output GeoTIFF path.
#' @param overwrite Whether an existing output may be replaced.
#' @param workers Number of patch workers.
#'
#' @return A predicted `terra::SpatRaster`.
#' @export
ubestarfm_predict <- function(
  model,
  coarse_target,
  value_range = c(250, 350),
  output_path = NULL,
  overwrite = FALSE,
  workers = 1L
) {
  output_paths <- if (is.null(output_path)) NULL else list(output_path)
  ubestarfm_predict_batch(
    model,
    list(coarse_target),
    output_paths = output_paths,
    value_range = value_range,
    overwrite = overwrite,
    workers = workers
  )[[1L]]
}

#' Train once and predict a target batch
#'
#' @inheritParams ubestarfm_train
#' @param coarse_targets A list of target rasters.
#' @param output_paths Optional output paths.
#' @param value_range Allowed output range.
#' @param overwrite Whether existing outputs may be replaced.
#'
#' @return A list of predicted `terra::SpatRaster` objects.
#' @export
ubestarfm_fit_predict_batch <- function(
  fine_1,
  fine_2,
  coarse_1,
  coarse_2,
  coarse_targets,
  output_paths = NULL,
  window_radius = 25L,
  patch_size = 200L,
  method = "zero_bias",
  value_range = c(250, 350),
  overwrite = FALSE,
  workers = 1L
) {
  model <- ubestarfm_train(
    fine_1,
    fine_2,
    coarse_1,
    coarse_2,
    window_radius = window_radius,
    patch_size = patch_size,
    method = method,
    workers = workers
  )
  ubestarfm_predict_batch(
    model,
    coarse_targets,
    output_paths = output_paths,
    value_range = value_range,
    overwrite = overwrite,
    workers = workers
  )
}
