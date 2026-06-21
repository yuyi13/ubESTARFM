# Script: parallel.R
# Objective: Execute deterministic patch work sequentially or in parallel.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Work items, a worker function, and a worker count.
# Outputs: Ordered lists of patch results.
# Usage: Internal package module; loaded through the ubestarfm namespace.
# Dependencies: Base R parallel package.

ubestarfm_lapply <- function(x, fun, workers = 1L) {
  workers <- as.integer(workers)
  if (
    length(workers) != 1L ||
      is.na(workers) ||
      workers < 1L
  ) {
    stop("workers must be a positive integer.", call. = FALSE)
  }

  if (workers == 1L || length(x) <= 1L) {
    return(lapply(x, fun))
  }

  if (.Platform$OS.type != "windows") {
    return(parallel::mclapply(
      x,
      fun,
      mc.cores = min(workers, length(x)),
      mc.preschedule = TRUE
    ))
  }

  cluster <- parallel::makeCluster(min(workers, length(x)))
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, library(ubestarfm))
  parallel::parLapply(cluster, x, fun)
}
