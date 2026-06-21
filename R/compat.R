# Script: compat.R
# Objective: Preserve the published single-target ubESTARFM function signature.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Published function arguments and five reference or target rasters.
# Outputs: A fine-resolution SpatRaster and output GeoTIFF.
# Usage: Internal package module; loaded through the ubESTARFM namespace.
# Dependencies: R package ubESTARFM.

#' Published ubESTARFM compatibility wrapper
#'
#' @param w Moving-window radius in pixels.
#' @param DN_min Lower plausible output value.
#' @param DN_max Upper plausible output value.
#' @param patch_long Patch side length in pixels.
#' @param tmp_path Accepted for source compatibility; no longer used.
#' @param out_path Output GeoTIFF path.
#' @param method Either `"zero bias"` or `"baseline"`.
#' @param rst_fine1 Fine raster for reference date one.
#' @param rst_fine2 Fine raster for reference date two.
#' @param rst_coarse1 Coarse raster for reference date one.
#' @param rst_coarse2 Coarse raster for reference date two.
#' @param rst_coarse0 Coarse raster for the prediction date.
#'
#' @return A predicted `terra::SpatRaster`.
#' @export
ubESTARFM <- function(
  w = 25,
  DN_min,
  DN_max,
  patch_long = 200,
  tmp_path = tempdir(),
  out_path,
  method = "zero bias",
  rst_fine1,
  rst_fine2,
  rst_coarse1,
  rst_coarse2,
  rst_coarse0
) {
  modern_method <- switch(
    method,
    "zero bias" = "zero_bias",
    "baseline" = "baseline",
    stop("method must be 'zero bias' or 'baseline'.", call. = FALSE)
  )
  invisible(tmp_path)

  model <- ubestarfm_train(
    fine_1 = rst_fine1,
    fine_2 = rst_fine2,
    coarse_1 = rst_coarse1,
    coarse_2 = rst_coarse2,
    window_radius = as.integer(w),
    patch_size = as.integer(patch_long),
    method = modern_method
  )
  ubestarfm_predict(
    model,
    coarse_target = rst_coarse0,
    value_range = c(DN_min, DN_max),
    output_path = out_path,
    overwrite = TRUE
  )
}
