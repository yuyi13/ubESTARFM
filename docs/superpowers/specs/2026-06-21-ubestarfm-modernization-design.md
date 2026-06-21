# ubESTARFM Modernization Design

**Date:** 2026-06-21  
**Status:** Pending written-spec review  
**Scientific reference:** `legacy/ubESTARFM.R`, moved without content changes from the published `0_algorithm/ubESTARFM.R`

## 1. Purpose

Modernize ubESTARFM so a pair of fine/coarse reference dates is prepared once and then reused to predict one or more target dates. Deliver maintained R and Python implementations, examples, tests, benchmarks, and a repository layout that clearly separates current software from archival paper-processing code.

The published `ubESTARFM.R` implementation is the sole scientific reference. The untracked and previously unsuccessful `ubESTARFM_fast.R` is not a reference implementation and will not be repaired or incorporated.

## 2. Scope

The implementation is optimized for site-sized rasters that fit in memory. The bundled 400 by 400 GeoTIFFs are the primary functional and performance fixtures. Continental and out-of-core processing are outside this release.

The work includes:

- A modern R package using `terra`, Rcpp, and optional base-R parallel execution.
- An equivalent Python package using Rasterio, NumPy, Numba, and optional standard-library parallel execution.
- Explicit train, predict, and batch-predict APIs in both languages.
- A compatibility `ubESTARFM()` R wrapper for existing one-target workflows.
- Compact reusable model storage.
- Cross-language scientific validation.
- R and Python tutorials and runnable examples.
- Repository restructuring, automated checks, migration documentation, and reproducible benchmarks.

The work does not include:

- Rewriting all archival paper-processing scripts to modern libraries.
- Supporting rasters that cannot fit in memory.
- GPU execution.
- A shared binary model format between R and Python.
- Persisting every candidate weight for every output pixel.

## 3. Repository Layout

The repository will use descriptive directories rather than numbered workflow prefixes:

```text
R/                              Modern R public API and internal helpers
src/                            Rcpp computational kernels
python/
  pyproject.toml                Python package and development metadata
  src/ubestarfm/                Python public API and computational kernels
  tests/                        Python unit and integration tests
tests/
  testthat/                     R unit and integration tests
  fixtures/                     Small language-neutral numerical fixtures
inst/
  extdata/                      Bundled example GeoTIFFs
examples/
  R/                            R single-target and batch examples
  python/                       Python single-target and batch examples
  outputs/                      Generated outputs; ignored except README
legacy/
  ubESTARFM.R                   Byte-for-byte published R implementation
  examples/                     Former published example and visualization
  paper-processing/             Former 4_lst_processing_scripts
data/
  ozflux/                       Former 5_ozflux_lst
docs/
  tutorials/                    R and Python user tutorials
  migration/                    Old-to-new path and API guidance
  superpowers/specs/            Approved design specifications
  superpowers/plans/            Implementation plans
benchmarks/                     Reproducible timing and memory benchmarks
figures/                        Documentation figures
```

Path migration:

| Current path | New path |
|---|---|
| `0_algorithm/ubESTARFM.R` | `legacy/ubESTARFM.R` |
| `0_algorithm/example.R` | `legacy/examples/example.R` |
| `0_algorithm/visualise.R` | `legacy/examples/visualise.R` |
| `1_test_data/` | `inst/extdata/` |
| `2_tmp_path/` | Removed; safe temporary directories are used |
| `3_output/` | `examples/outputs/` |
| `4_lst_processing_scripts/` | `legacy/paper-processing/` |
| `5_ozflux_lst/` | `data/ozflux/` |

The last commit with the old layout, `4721ac0`, will receive an annotated preservation tag before the modernization release. The modernization release will use a new major version because it changes package structure and public APIs.

## 4. Scientific Compatibility

### 4.1 Reference behavior

The modern engines retain the published algorithm's:

- Patch-specific similar-pixel thresholds.
- Similar-pixel selection using both fine-resolution reference dates.
- Zero-bias correction and baseline modes.
- Spectral and spatial distance formulas.
- Temporal weighting formulas.
- Five-candidate cutoff and fallback calculation.
- Lower and upper output-value correction.
- Edge-window behavior.

R's sample standard deviation is used for patch thresholds. Python uses the equivalent `ddof = 1`. Candidate ordering is defined canonically as column-major within each local window to match R matrix vectorization. Both implementations use double precision for calculations and preserve GeoTIFF spatial metadata.

### 4.2 Deliberate corrections

The modern implementation corrects behavior that is demonstrably defective or prevents reusable training:

- Target-date values are excluded from reference-pair training.
- At prediction time, the stored reference candidate mask is intersected with finite target-date pixels.
- Candidate weights are renormalized after target-date filtering.
- If target filtering leaves five or fewer candidates, the published fallback branch is used with the target-valid window.
- Invalid method names, incompatible dimensions, grids, coordinate reference systems, and value ranges raise explicit errors.
- Temporary output discovery is not based on every `.tif` in a shared directory.
- Patch outputs are assembled by location rather than averaged through a raster stack.

The target coarse raster does not change reference-date similarity membership when it is complete. However, its missing-value pattern changes the effective valid candidate set in the published algorithm. Therefore, fully reusable coefficients are guaranteed only for targets sharing the same validity mask. The modern batch engine uses a fast path for shared or complete target masks and performs target-specific filtering when masks differ.

### 4.3 Numerical acceptance

For finite pixels where the published behavior is defined:

- Modern R versus the tracked published example output: maximum absolute difference no greater than `1e-4` K.
- Python versus modern R: maximum absolute difference no greater than `1e-4` K.
- Unit-level kernels in R and Python: maximum absolute difference no greater than `1e-8`.
- Missing-value masks and output geometry must match exactly.

Differences caused solely by corrected target-date missing-value handling are tested separately and documented in migration notes.

## 5. Public APIs

### 5.1 R

```r
model <- ubestarfm_train(
  fine_1,
  fine_2,
  coarse_1,
  coarse_2,
  window_radius = 25L,
  patch_size = 200L,
  method = "zero_bias",
  cache = "candidates",
  workers = 1L
)

prediction <- ubestarfm_predict(
  model,
  coarse_target,
  value_range = c(250, 350),
  output_path = NULL,
  overwrite = FALSE,
  workers = 1L
)

predictions <- ubestarfm_predict_batch(
  model,
  coarse_targets,
  output_paths = NULL,
  value_range = c(250, 350),
  overwrite = FALSE,
  workers = 1L
)
```

`fine_1`, `fine_2`, `coarse_1`, `coarse_2`, and target inputs accept `terra::SpatRaster` objects or file paths. Batch targets accept a list of rasters or paths. Results are `SpatRaster` objects when retained in memory and are also written when output paths are supplied.

The compatibility wrapper retains the published argument names:

```r
ubESTARFM(
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
)
```

It translates `"zero bias"` to `"zero_bias"`, trains one model, predicts one target, writes the requested output, and returns the result. `tmp_path` is accepted for source compatibility but is not used for unsafe file discovery.

Model persistence:

```r
ubestarfm_save_model(model, "model.rds", compress = TRUE)
model <- ubestarfm_load_model("model.rds")
```

### 5.2 Python

```python
from ubestarfm import (
    train,
    predict,
    predict_batch,
    save_model,
    load_model,
)

model = train(
    fine_1,
    fine_2,
    coarse_1,
    coarse_2,
    window_radius=25,
    patch_size=200,
    method="zero_bias",
    cache="candidates",
    workers=1,
)

prediction = predict(
    model,
    coarse_target,
    value_range=(250.0, 350.0),
    output_path=None,
    overwrite=False,
    workers=1,
)

predictions = predict_batch(
    model,
    coarse_targets,
    output_paths=None,
    value_range=(250.0, 350.0),
    overwrite=False,
    workers=1,
)
```

Inputs accept paths, Rasterio dataset readers, or two-dimensional NumPy arrays accompanied by an explicit raster profile. File-based results preserve the reference fine raster's transform, coordinate reference system, dimensions, and nodata convention.

Python models are saved as compressed NumPy archives plus JSON metadata. R and Python model files are intentionally language-specific; parity is enforced at the prediction level.

## 6. Model Representation and Memory

### 6.1 Why exact weight storage is rejected

For the bundled 400 by 400 reference rasters with `window_radius = 25`, a measured sample contains approximately 1,380 similar pixels per valid target pixel on average. Storing a 32-bit candidate index and 64-bit weight for every relationship would require approximately 2.5 GiB before general object overhead.

The failed fast implementation additionally stored nested lists, string keys, coordinates, duplicated fine and coarse values, bounds, and validity indices. This explains previously observed model sizes of hundreds of megabytes or more.

### 6.2 Candidate-bitset model

For each output pixel, the model stores a fixed 51 by 51 candidate-membership bitset. A set bit means the corresponding local-window position:

1. Is valid in all four reference rasters.
2. Meets the fine-date-one similarity threshold.
3. Meets the fine-date-two similarity threshold.

The model also stores:

- Four reference arrays.
- The reference validity mask.
- Patch identifiers and patch thresholds.
- Raster dimensions, transform, coordinate reference system, nodata value, and data type metadata.
- Window radius, patch size, method, model schema version, and package version.

The model does not store candidate weights, duplicated candidate values, coordinate vectors, or per-pixel nested objects. Bias terms and fallback statistics are calculated from the compact mask because target-date validity can alter the effective candidate set.

For a 400 by 400 raster and a 51 by 51 window, candidate bitsets require approximately 50 MiB. Reference arrays and metadata add approximately 10 to 20 MiB. The acceptance ceiling for the uncompressed in-memory candidate-cache model on the bundled fixture is 80 MiB. The compressed serialized size is measured and reported but is not constrained to a machine-independent byte ceiling.

### 6.3 Cache modes

- `cache = "candidates"` is the default. It creates the reusable bitset model and supports saving, loading, and later prediction.
- `cache = "none"` is available through a combined fit-and-predict helper used for one-off batch processing. It processes candidate masks patch by patch, applies them to all supplied targets, and discards them. It is not saveable.

The explicit `train()` API always produces a reusable candidate-cache model. The one-off helper is named `ubestarfm_fit_predict_batch()` in R and `fit_predict_batch()` in Python.

## 7. Computational Design

### 7.1 Training

Training:

1. Loads and validates the four reference rasters.
2. Converts nodata values to missing values.
3. Divides the raster into deterministic non-overlapping patches.
4. Computes the two sample-standard-deviation thresholds for each patch.
5. For every valid reference pixel, evaluates the local window and writes candidate membership into its fixed bitset.
6. Returns a versioned model with reference data and spatial metadata.

Training contains no target-date raster and is executed once per reference pair.

### 7.2 Prediction

For each pixel and target:

1. Skip the pixel when its target coarse value is missing.
2. Intersect the stored candidate bitset with finite target pixels in the local window.
3. Use the detailed branch when more than five candidates remain.
4. Reconstruct candidate reference values and coordinates directly from flat arrays.
5. Apply zero-bias correction when requested.
6. Compute spectral/spatial weights and temporal weights.
7. Predict from both reference pairs and combine them.
8. Apply the published abnormal-value fallback.
9. Use the published low-candidate fallback when five or fewer candidates remain.

Batch prediction loops over pixels and reconstructs pair-dependent candidate information once when target validity masks permit reuse. Target values and temporal weights remain target-specific.

### 7.3 R engine

- `terra` handles GeoTIFF input, output, geometry validation, nodata conversion, and raster metadata.
- Rcpp implements candidate-bitset creation and prediction kernels.
- Base `parallel` provides optional patch-level execution. Sequential execution is the default to avoid unnecessary process and memory overhead for small sites.
- Worker cleanup is guaranteed with `on.exit()` where clusters are used.
- No `foreach`, `doParallel`, or `raster` dependency is introduced into the modern package.

### 7.4 Python engine

- Rasterio handles GeoTIFF input, output, masks, and metadata.
- NumPy owns contiguous numeric arrays and model serialization.
- Numba compiles candidate-bitset creation and prediction kernels.
- Standard-library concurrent execution is optional at the patch level. Sequential execution is the default.
- The public API does not require Xarray, Dask, or Joblib.

## 8. Testing Strategy

### 8.1 R tests

`testthat` covers:

- Input and geometry validation.
- Patch threshold calculation.
- Candidate bit packing and unpacking.
- Detailed and fallback prediction branches.
- Zero-bias and baseline modes.
- Edge windows.
- Missing target pixels and candidate renormalization.
- Value-range correction.
- Model save/load round trips.
- Single versus batch prediction equality.
- Compatibility-wrapper behavior.
- Cache-size ceiling.

### 8.2 Python tests

`pytest` covers the equivalent Python behavior, serialization, raster metadata, and public API.

### 8.3 Cross-language fixtures

Small hand-checked arrays in `tests/fixtures/` exercise exact kernel behavior. R and Python write prediction arrays to a shared fixture format, and a comparison script enforces the stated tolerances.

The tracked published `fused_result.tif` becomes an immutable legacy golden output under `tests/fixtures/legacy/`. Its checksum is recorded. A documented regeneration script uses the archived published implementation and its legacy dependencies, but normal CI does not depend on deprecated packages.

### 8.4 End-to-end tests

The bundled Yanco example verifies:

- One training run followed by one prediction.
- One training run followed by at least two target predictions.
- Equality between a batch result and the corresponding individual prediction.
- Modern R agreement with the published golden output.
- Python agreement with modern R.
- Preservation of dimensions, transform, coordinate reference system, nodata mask, and output range.

## 9. Documentation and Examples

The README becomes a concise entry point with:

- Scientific overview and citation.
- Installation for R and Python.
- Single-target and multi-target quick starts.
- Links to tutorials, migration notes, benchmarks, and legacy code.

Tutorials:

- `docs/tutorials/r-single-and-batch.md`
- `docs/tutorials/python-single-and-batch.md`
- `docs/tutorials/model-caching.md`
- `docs/migration/from-published-r-api.md`
- `docs/migration/repository-layout.md`

Examples use project-root-relative paths, create output directories explicitly, never call `setwd()`, and clean up resources. Revised R and Python scripts follow the workspace script-header standard. Revised plotting code explicitly requests Arial where supported, and R plotting uses base graphics unless an existing archival file already depends on another plotting system.

## 10. Repository Quality

The modernization adds:

- R `DESCRIPTION`, `NAMESPACE`, generated documentation, `testthat`, and `.Rbuildignore`.
- Python `pyproject.toml` with build, runtime, test, and formatting configuration.
- `.gitignore` entries for generated outputs, caches, temporary files, compiled objects, and local environments.
- GitHub Actions for R checks, Python tests, syntax checks, and cross-language fixture comparison.
- `CONTRIBUTING.md` with test and style commands.
- `CHANGELOG.md` documenting the modernization and intentional numerical differences.
- License metadata in both package systems consistent with the repository MIT license.

The existing syntax error in the archival `10_figures.R` is corrected before it is moved. Archival processing scripts otherwise retain their scientific content; they receive compliant headers only if they are revised during the migration.

Permission-only changes from the original dirty working tree are not included. Data, documentation, and raster files remain non-executable. Only scripts with valid shebangs may be executable.

## 11. Benchmarks and Acceptance

Benchmarks run on the bundled 400 by 400 rasters and report:

- Published one-target elapsed time and peak memory where legacy dependencies are available.
- Modern training time.
- Modern one-target prediction time after training.
- Modern two-target and multi-target batch time.
- Candidate-cache in-memory and serialized sizes.
- R/Python output differences.

Timing results are descriptive rather than enforced in CI because shared runners are variable. Structural efficiency is enforced:

- Candidate search runs once per reference pair.
- Batch prediction does not retrain.
- Single and batch outputs are numerically equivalent.
- The bundled reusable model remains at or below 80 MiB uncompressed in memory.
- No per-pixel nested relationship lists or exact-weight cache are created.

## 12. Error Handling

Both APIs fail early with actionable messages for:

- Non-single-layer rasters.
- Different dimensions, extents, resolutions, transforms, or coordinate reference systems.
- No valid reference pixels.
- Patch dimensions or window radii outside valid bounds.
- Unsupported methods or cache modes.
- Invalid output ranges.
- Existing outputs when `overwrite = FALSE`.
- Missing output-parent directories that cannot be created.
- Model schema or package-version incompatibility.

Interrupted writes use temporary files in the destination directory followed by atomic replacement when supported. Partial outputs are removed on failure.

## 13. Delivery Sequence

Implementation is divided into independently verifiable stages:

1. Restructure the repository and preserve the published implementation and fixtures.
2. Establish R and Python package skeletons, test infrastructure, and shared fixtures.
3. Implement and test candidate-bitset storage.
4. Implement the R training and prediction kernels and compatibility wrapper.
5. Implement the equivalent Python engine.
6. Add cross-language and published-output validation.
7. Add examples, tutorials, migration documentation, and benchmarks.
8. Add CI and complete repository-quality checks.

Each behavioral change follows test-driven development. Final verification includes clean Git status, R package checks, Python tests, cross-language comparisons, syntax checks, example execution, cache-size measurement, and benchmark generation.
