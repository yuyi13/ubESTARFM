#!/usr/bin/env Rscript
# Script: run_r_fixture.R
# Objective: Generate R predictions and model-size evidence for parity tests.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled GeoTIFFs and one output-directory command-line argument.
# Outputs: Two prediction GeoTIFFs, one RDS model, and a model-size text file.
# Usage: Rscript tests/cross_language/run_r_fixture.R <output-directory>
# Dependencies: R package ubESTARFM.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_r_fixture.R <output-directory>", call. = FALSE)
}

root <- normalizePath(".", mustWork = TRUE)
output_directory <- normalizePath(args[1L], mustWork = TRUE)
data_directory <- file.path(root, "inst", "extdata")

library(ubESTARFM)

model <- ubestarfm_train(
  fine_1 = file.path(
    data_directory,
    "Landsat_LST_cloudrm_20160205.tif"
  ),
  fine_2 = file.path(
    data_directory,
    "Landsat_LST_cloudrm_20160308.tif"
  ),
  coarse_1 = file.path(
    data_directory,
    "MOD11A1_LST_cloudrm_20160205.tif"
  ),
  coarse_2 = file.path(
    data_directory,
    "MOD11A1_LST_cloudrm_20160308.tif"
  ),
  window_radius = 25L,
  patch_size = 200L,
  workers = 4L
)

target <- file.path(
  data_directory,
  "MOD11A1_LST_cloudrm_20160218.tif"
)
ubestarfm_predict(
  model,
  target,
  output_path = file.path(output_directory, "r_prediction.tif"),
  workers = 4L
)

model_path <- file.path(output_directory, "r_model.rds")
ubestarfm_save_model(model, model_path, compress = TRUE)
reloaded <- ubestarfm_load_model(model_path)
ubestarfm_predict(
  reloaded,
  target,
  output_path = file.path(
    output_directory,
    "r_prediction_reloaded.tif"
  ),
  workers = 4L
)

writeLines(
  as.character(as.numeric(object.size(model))),
  file.path(output_directory, "r_model_size_bytes.txt")
)
