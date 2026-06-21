#!/usr/bin/env Rscript
# Script: benchmark.R
# Objective: Benchmark modern R ubESTARFM training, prediction, and model memory.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled GeoTIFFs and --output-dir.
# Outputs: R benchmark JSON and Markdown reports.
# Usage: Rscript benchmarks/benchmark.R --output-dir benchmarks/results
# Dependencies: R package ubestarfm.

args <- commandArgs(trailingOnly = TRUE)
output_index <- match("--output-dir", args)
if (is.na(output_index) || output_index == length(args)) {
  stop("--output-dir must be followed by a directory.", call. = FALSE)
}
output_directory <- args[output_index + 1L]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
data_directory <- file.path(root, "inst", "extdata")

library(ubestarfm)

training <- system.time({
  model <- ubestarfm_train(
    file.path(data_directory, "Landsat_LST_cloudrm_20160205.tif"),
    file.path(data_directory, "Landsat_LST_cloudrm_20160308.tif"),
    file.path(data_directory, "MOD11A1_LST_cloudrm_20160205.tif"),
    file.path(data_directory, "MOD11A1_LST_cloudrm_20160308.tif"),
    workers = 4L
  )
})
target <- file.path(data_directory, "MOD11A1_LST_cloudrm_20160218.tif")
single <- system.time({
  prediction <- ubestarfm_predict(model, target, workers = 4L)
})
batch <- system.time({
  predictions <- ubestarfm_predict_batch(
    model,
    list(target, target),
    workers = 4L
  )
})
invisible(prediction)
invisible(predictions)

metrics <- c(
  training_seconds = unname(training["elapsed"]),
  single_prediction_seconds = unname(single["elapsed"]),
  two_target_batch_seconds = unname(batch["elapsed"]),
  model_size_bytes = as.numeric(object.size(model))
)
json <- sprintf(
  paste0(
    "{\n",
    "  \"training_seconds\": %.6f,\n",
    "  \"single_prediction_seconds\": %.6f,\n",
    "  \"two_target_batch_seconds\": %.6f,\n",
    "  \"model_size_bytes\": %.0f\n",
    "}\n"
  ),
  metrics[1L],
  metrics[2L],
  metrics[3L],
  metrics[4L]
)
writeLines(json, file.path(output_directory, "r_benchmark.json"))
writeLines(
  c(
    "# R benchmark",
    "",
    sprintf("- Training: %.3f seconds", metrics[1L]),
    sprintf("- Warm single prediction: %.3f seconds", metrics[2L]),
    sprintf("- Two-target batch: %.3f seconds", metrics[3L]),
    sprintf("- Model size: %.2f MiB", metrics[4L] / 1024^2)
  ),
  file.path(output_directory, "r_benchmark.md")
)
