# Changelog

## 3.0.0

- Added reusable train, single-predict, and batch-predict APIs in R and Python.
- Replaced maintained `raster`, `foreach`, and `doParallel` use with `terra`,
  Rcpp, Rasterio, NumPy, and Numba.
- Added compact packed candidate caches, model persistence, tests, examples,
  tutorials, benchmarks, and cross-language validation.
- Preserved the published R implementation under `legacy/`.
- Replaced numbered repository directories with descriptive paths.
- Corrected target-date missing-value handling by filtering candidates and
  renormalizing weights during prediction.
- Standardized the maintained R package identity as lowercase `ubestarfm`
  while retaining the published `ubESTARFM()` compatibility function.
