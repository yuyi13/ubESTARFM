#!/usr/bin/env Rscript
# Script: batch_targets.R
# Objective: Demonstrate R batch prediction after one reference-pair training.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled GeoTIFFs under inst/extdata and --output-dir.
# Outputs: Two fused GeoTIFFs and one reusable RDS model.
# Usage: Rscript examples/R/batch_targets.R --output-dir examples/outputs
# Dependencies: R package ubESTARFM.

args <- commandArgs(trailingOnly = TRUE)
output_index <- match("--output-dir", args)
if (is.na(output_index) || output_index == length(args)) {
  stop("--output-dir must be followed by a directory.", call. = FALSE)
}
output_directory <- normalizePath(
  args[output_index + 1L],
  mustWork = TRUE
)
quick <- "--quick" %in% args
window_radius <- if (quick) 2L else 25L

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
data_directory <- file.path(root, "inst", "extdata")

library(ubESTARFM)

model <- ubestarfm_train(
  fine_1 = file.path(data_directory, "Landsat_LST_cloudrm_20160205.tif"),
  fine_2 = file.path(data_directory, "Landsat_LST_cloudrm_20160308.tif"),
  coarse_1 = file.path(data_directory, "MOD11A1_LST_cloudrm_20160205.tif"),
  coarse_2 = file.path(data_directory, "MOD11A1_LST_cloudrm_20160308.tif"),
  window_radius = window_radius,
  patch_size = 200L,
  workers = 4L
)
target <- file.path(data_directory, "MOD11A1_LST_cloudrm_20160218.tif")
predictions <- ubestarfm_predict_batch(
  model,
  list(target, target),
  output_paths = file.path(
    output_directory,
    c("r_batch_20160218_a.tif", "r_batch_20160218_b.tif")
  ),
  overwrite = TRUE,
  workers = 4L
)
ubestarfm_save_model(
  model,
  file.path(output_directory, "r_reference_pair_model.rds")
)
invisible(predictions)

message(
  "The bundled repository contains one target date, so this tutorial repeats ",
  "20160218 to demonstrate that a batch reuses one trained model."
)
