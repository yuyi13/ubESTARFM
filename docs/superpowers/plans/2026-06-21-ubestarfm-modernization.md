# ubESTARFM Modernization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace repeated per-target ubESTARFM fitting with reusable reference-pair training, maintained R and Python engines, modern raster I/O, compact candidate caching, reproducible validation, and a descriptive repository layout.

**Architecture:** Preserve the published R implementation unchanged under `legacy/` and implement independent modern engines against its scientific behavior. Both engines store fixed-window candidate membership as packed bits, reconstruct weights during prediction, and expose train, single-predict, batch-predict, and one-off fit/predict APIs. R uses `terra` and Rcpp; Python uses Rasterio, NumPy, and Numba.

**Tech Stack:** R 4.3+, terra 1.9+, Rcpp 1.1+, testthat; Python 3.12+, NumPy 2+, Rasterio 1.5+, Numba 0.65+, pytest; GitHub Actions.

**Approved design:** `docs/superpowers/specs/2026-06-21-ubestarfm-modernization-design.md`

---

## File Map

Every new `.R` and `.py` file listed below must begin with the workspace-standard
header. Use `2026-06-21` for `Created` and `Last updated`, `Yi Yu` for `Author`,
the file's stated responsibility for `Objective`, its public arguments or
fixtures for `Inputs`, and its return values/files for `Outputs`. Package modules
use `Usage: Internal package module; loaded through the ubESTARFM namespace.`
or `Usage: Imported by the ubestarfm package.` Tests use the exact focused test
command from their task as `Usage`. Copy the language-specific structure from
`~/.codex/templates/script_headers/`; do not abbreviate or omit fields.

### Preserved and migrated material

- Move unchanged: `0_algorithm/ubESTARFM.R` → `legacy/ubESTARFM.R`
- Move unchanged: `0_algorithm/example.R` → `legacy/examples/example.R`
- Move unchanged: `0_algorithm/visualise.R` → `legacy/examples/visualise.R`
- Move and repair syntax/header: `4_lst_processing_scripts/10_figures.R` → `legacy/paper-processing/10_figures.R`
- Move unchanged: remaining `4_lst_processing_scripts/*` → `legacy/paper-processing/`
- Move unchanged: `1_test_data/*` → `inst/extdata/`
- Move unchanged: `5_ozflux_lst/*` → `data/ozflux/`
- Move golden result: `3_output/fused_result.tif` → `tests/fixtures/legacy/fused_result.tif`
- Move rendered example: `3_output/visualisation.png` → `examples/outputs/legacy_visualisation.png`
- Move archival PDF: `README.pdf` → `legacy/README.pdf`
- Delete obsolete: `2_tmp_path/README.md`

### R package

- Create: `DESCRIPTION`
- Create: `NAMESPACE`
- Create: `.Rbuildignore`
- Create: `R/io.R`
- Create: `R/model.R`
- Create: `R/parallel.R`
- Create: `R/train.R`
- Create: `R/predict.R`
- Create: `R/compat.R`
- Create/generated: `R/RcppExports.R`
- Create: `src/candidate_masks.cpp`
- Create: `src/predict.cpp`
- Create/generated: `src/RcppExports.cpp`
- Create: `src/Makevars`
- Create: `src/Makevars.win`
- Create: `tests/testthat.R`
- Create: `tests/testthat/helper-fixtures.R`
- Create: focused `tests/testthat/test-*.R` files

### Python package

- Create: `python/pyproject.toml`
- Create: `python/README.md`
- Create: `python/src/ubestarfm/__init__.py`
- Create: `python/src/ubestarfm/io.py`
- Create: `python/src/ubestarfm/model.py`
- Create: `python/src/ubestarfm/kernels.py`
- Create: `python/src/ubestarfm/api.py`
- Create: `python/tests/conftest.py`
- Create: focused `python/tests/test_*.py` files

### Shared verification and documentation

- Create: `tests/fixtures/legacy/SHA256SUMS`
- Create: `tests/fixtures/small/*.csv`
- Create: `tests/cross_language/run_r_fixture.R`
- Create: `tests/cross_language/compare_outputs.py`
- Create: `examples/R/single_target.R`
- Create: `examples/R/batch_targets.R`
- Create: `examples/python/single_target.py`
- Create: `examples/python/batch_targets.py`
- Create: `benchmarks/benchmark.R`
- Create: `benchmarks/benchmark.py`
- Create: `docs/tutorials/*.md`
- Create: `docs/migration/*.md`
- Create: `.github/workflows/r-check.yaml`
- Create: `.github/workflows/python-test.yaml`
- Create: `.github/workflows/cross-language.yaml`
- Rewrite: `README.md`
- Create: `CONTRIBUTING.md`
- Create: `CHANGELOG.md`
- Create: `.gitignore`

---

### Task 1: Preserve the Published Version and Migrate the Repository Layout

**Files:**

- Move the paths listed under “Preserved and migrated material”
- Create: `tests/fixtures/legacy/SHA256SUMS`
- Create: `tests/layout/test_layout.py`
- Create: `.gitignore`
- Create: `docs/migration/repository-layout.md`

- [ ] **Step 1: Write the failing layout and preservation test**

Create `tests/layout/test_layout.py`:

```python
from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def test_descriptive_layout_exists() -> None:
    expected = [
        ROOT / "legacy" / "ubESTARFM.R",
        ROOT / "inst" / "extdata" / "MOD11A1_LST_cloudrm_20160218.tif",
        ROOT / "data" / "ozflux" / "Yanco_lst.csv",
        ROOT / "tests" / "fixtures" / "legacy" / "fused_result.tif",
        ROOT / "examples" / "outputs" / "legacy_visualisation.png",
    ]
    assert all(path.exists() for path in expected)


def test_numbered_layout_is_removed() -> None:
    assert not any((ROOT / name).exists() for name in (
        "0_algorithm",
        "1_test_data",
        "2_tmp_path",
        "3_output",
        "4_lst_processing_scripts",
        "5_ozflux_lst",
    ))


def test_published_algorithm_is_byte_for_byte_preserved() -> None:
    assert sha256(ROOT / "legacy" / "ubESTARFM.R") == (
        "131f15efab127deccbe411cf800ea7967e621e8c916de677430299448f47e67b"
    )


def test_legacy_golden_output_is_preserved() -> None:
    assert sha256(ROOT / "tests" / "fixtures" / "legacy" / "fused_result.tif") == (
        "c4d216842fe9624f5355feac88ba6cb1b8dd35afd0e66272578561bc2c92a490"
    )
```

- [ ] **Step 2: Run the test and verify it fails because the new paths do not exist**

Run:

```bash
python3 -m pytest tests/layout/test_layout.py -v
```

Expected: FAIL in `test_descriptive_layout_exists` and `test_numbered_layout_is_removed`.

- [ ] **Step 3: Add an annotated local preservation tag**

Run:

```bash
git tag -a legacy-layout-v2.0.1 4721ac0 \
  -m "Last published repository layout before the ubESTARFM modernization"
```

Do not push the tag in this task.

- [ ] **Step 4: Move preserved content with Git-aware operations**

Run:

```bash
mkdir -p legacy/examples legacy/paper-processing inst/extdata data/ozflux
mkdir -p tests/fixtures/legacy examples/outputs docs/migration
git mv 0_algorithm/ubESTARFM.R legacy/ubESTARFM.R
git mv 0_algorithm/example.R legacy/examples/example.R
git mv 0_algorithm/visualise.R legacy/examples/visualise.R
git mv 4_lst_processing_scripts/* legacy/paper-processing/
rmdir 0_algorithm 4_lst_processing_scripts
git mv 1_test_data/* inst/extdata/
rmdir 1_test_data
git mv 5_ozflux_lst/* data/ozflux/
rmdir 5_ozflux_lst
git mv 3_output/fused_result.tif tests/fixtures/legacy/fused_result.tif
git mv 3_output/visualisation.png examples/outputs/legacy_visualisation.png
rmdir 3_output
git mv README.pdf legacy/README.pdf
git rm 2_tmp_path/README.md
rmdir 2_tmp_path
```

- [ ] **Step 5: Repair only the known archival syntax error and add its required header**

At the top of `legacy/paper-processing/10_figures.R`, add:

```r
#!/usr/bin/env Rscript
# Script: 10_figures.R
# Objective: Reproduce the figures used in the published ubESTARFM study.
# Author: Yi Yu
# Created: 2023-01-01
# Last updated: 2026-06-21
# Inputs: Archived study metrics and raster products at the configured paths.
# Outputs: Publication figure files.
# Usage: Rscript legacy/paper-processing/10_figures.R
# Dependencies: R packages stringr, ggsci, raster, and fields; archived study data.
```

Delete the invalid standalone line:

```r
xlim=c(1,6), ylim=c(1,3)
```

Do not otherwise refactor this archival script. Set all data and documentation files to mode `0644`. Set scripts to `0755` only when they have a valid shebang; otherwise use `0644`.

- [ ] **Step 6: Add checksum, ignore, and migration files**

Create `tests/fixtures/legacy/SHA256SUMS`:

```text
131f15efab127deccbe411cf800ea7967e621e8c916de677430299448f47e67b  legacy/ubESTARFM.R
c4d216842fe9624f5355feac88ba6cb1b8dd35afd0e66272578561bc2c92a490  tests/fixtures/legacy/fused_result.tif
3884534784c6851c1dde90e37fe0e16cffedea4f62a637e86fab4f3a3eee9f69  examples/outputs/legacy_visualisation.png
```

Create `.gitignore`:

```gitignore
.Rcheck/
.Rhistory
.Rproj.user/
*.Rproj
src/*.o
src/*.so
src/*.dll
python/.venv/
python/.pytest_cache/
python/.ruff_cache/
python/**/__pycache__/
python/**/*.py[cod]
python/build/
python/dist/
python/*.egg-info/
examples/outputs/*
!examples/outputs/README.md
benchmarks/results/
*.rds
*.npz
```

Document every old-to-new path in `docs/migration/repository-layout.md`, including the golden-output exception from `3_output/`.

- [ ] **Step 7: Run layout, checksum, and syntax verification**

Run:

```bash
python3 -m pytest tests/layout/test_layout.py -v
sha256sum -c tests/fixtures/legacy/SHA256SUMS
Rscript -e 'files <- list.files("legacy", pattern = "[.]R$", recursive = TRUE, full.names = TRUE); for (f in files) parse(file = f); cat(length(files), "R files parsed\n")'
python3 -m py_compile legacy/paper-processing/02_process_landsat_lst.py
git diff --check
```

Expected: all commands exit 0. The R parse command reports 14 files.

- [ ] **Step 8: Commit the migration**

```bash
git add -A
git commit -m "refactor: adopt descriptive repository layout"
```

---

### Task 2: Establish the R Package, Raster I/O, and Model Contract

**Files:**

- Create: `DESCRIPTION`
- Create: `NAMESPACE`
- Create: `.Rbuildignore`
- Create: `R/io.R`
- Create: `R/model.R`
- Create: `tests/testthat.R`
- Create: `tests/testthat/helper-fixtures.R`
- Create: `tests/testthat/test-io.R`
- Create: `tests/testthat/test-model.R`

- [ ] **Step 1: Add package-level failing tests**

Create `tests/testthat.R`:

```r
library(testthat)
library(ubESTARFM)

test_check("ubESTARFM")
```

Create `tests/testthat/helper-fixtures.R`:

```r
example_path <- function(filename) {
  system.file("extdata", filename, package = "ubESTARFM", mustWork = TRUE)
}

small_reference_fixture <- function() {
  fine_1   <- matrix(seq_len(25), nrow = 5L, byrow = TRUE)
  fine_2   <- fine_1 + 100
  coarse_1 <- fine_1 + 200
  coarse_2 <- fine_1 + 300
  coarse_1[2, 2] <- NA_real_
  list(
    fine_1   = fine_1,
    fine_2   = fine_2,
    coarse_1 = coarse_1,
    coarse_2 = coarse_2
  )
}
```

Create `tests/testthat/test-io.R`:

```r
test_that("single-layer GeoTIFFs load as row-column matrices", {
  input <- ubestarfm_read_raster(example_path("MOD11A1_LST_cloudrm_20160218.tif"))

  expect_s4_class(input$raster, "SpatRaster")
  expect_equal(dim(input$values), c(400L, 400L))
  expect_equal(input$values[1, 1], 312.74, tolerance = 1e-4)
  expect_true(all(is.finite(input$values)))
})

test_that("reference rasters must share geometry", {
  first <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  second <- terra::rast(nrows = 4, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 4)

  expect_error(
    ubestarfm_validate_geometry(list(first, second)),
    "same dimensions, extent, resolution, and CRS"
  )
})
```

Create `tests/testthat/test-model.R`:

```r
test_that("model validation rejects malformed models", {
  expect_error(
    ubestarfm_validate_model(list(schema_version = 1L)),
    "ubestarfm_model"
  )
})
```

- [ ] **Step 2: Run the tests and verify missing-package/API failures**

Run:

```bash
Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'
```

Expected: FAIL because package functions do not exist.

- [ ] **Step 3: Create R package metadata**

Create `DESCRIPTION` with:

```text
Package: ubESTARFM
Type: Package
Title: Unbiased Enhanced Spatial and Temporal Adaptive Reflectance Fusion
Version: 3.0.0.9000
Authors@R: person("Yi", "Yu", email = "yi.yu1@anu.edu.au", role = c("aut", "cre"))
Description: Trains reusable reference-pair models and predicts fine-resolution
    land surface temperature using the unbiased ESTARFM algorithm.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.2
Depends: R (>= 4.3.0)
Imports:
    methods,
    Rcpp (>= 1.1.1),
    terra (>= 1.9.0)
LinkingTo:
    Rcpp
Suggests:
    testthat (>= 3.2.0)
Config/testthat/edition: 3
SystemRequirements: C++17
```

Create `.Rbuildignore` for `python`, `legacy`, `benchmarks`, `docs/superpowers`, and local build artifacts.

Before adding function code, add compliant headers to every new file under
`R/` and `tests/testthat/`. The shebang is required for directly executable
examples and runners; package modules and test files omit the shebang but place
the complete comment header at line 1.

- [ ] **Step 4: Implement raster input and geometry validation**

In `R/io.R`, implement and export:

```r
ubestarfm_read_raster <- function(x, name = deparse(substitute(x))) {
  raster <- if (inherits(x, "SpatRaster")) x else terra::rast(x)
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

ubestarfm_validate_geometry <- function(rasters) {
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
```

Add a private matrix-to-raster helper using `terra::setValues(template, as.vector(t(values)))`.

- [ ] **Step 5: Implement the versioned model contract**

In `R/model.R`, define:

```r
new_ubestarfm_model <- function(
  reference_values,
  reference_valid,
  candidate_masks,
  mask_bytes,
  patch_ids,
  patch_thresholds,
  metadata,
  parameters
) {
  structure(
    list(
      schema_version = 1L,
      package_version = as.character(utils::packageVersion("ubESTARFM")),
      reference_values = reference_values,
      reference_valid = reference_valid,
      candidate_masks = candidate_masks,
      mask_bytes = as.integer(mask_bytes),
      patch_ids = patch_ids,
      patch_thresholds = patch_thresholds,
      metadata = metadata,
      parameters = parameters
    ),
    class = "ubestarfm_model"
  )
}
```

`ubestarfm_validate_model()` must verify class, schema version, four equal-length reference arrays, raw mask length, dimensions, method, and required metadata.

- [ ] **Step 6: Generate namespace stubs and run tests**

Add roxygen export/import declarations, then run:

```bash
Rscript -e 'Rcpp::compileAttributes(".")'
Rscript -e 'roxygen2::roxygenise(".")'
R CMD INSTALL .
Rscript -e 'testthat::test_local(".", reporter = "summary")'
```

Expected: all I/O and model tests pass.

- [ ] **Step 7: Commit the R package foundation**

```bash
git add DESCRIPTION NAMESPACE .Rbuildignore R/ tests/testthat.R tests/testthat
git commit -m "feat(r): establish package and raster model contracts"
```

---

### Task 3: Implement Candidate-Bitset Training in R

**Files:**

- Create: `src/candidate_masks.cpp`
- Create: `src/Makevars`
- Create: `src/Makevars.win`
- Create/generated: `src/RcppExports.cpp`
- Create/generated: `R/RcppExports.R`
- Create: `R/parallel.R`
- Create: `R/train.R`
- Create: `tests/testthat/test-candidate-masks.R`
- Create: `tests/testthat/test-train.R`

- [ ] **Step 1: Write failing unit tests for bit ordering and candidate selection**

Use a 5 by 5 fixture with `window_radius = 1L`. Assert:

```r
test_that("candidate bits use column-major local-window order", {
  fixture <- small_reference_fixture()
  result <- ubestarfm_train_arrays(
    fixture$fine_1,
    fixture$fine_2,
    fixture$coarse_1,
    fixture$coarse_2,
    window_radius = 1L,
    patch_size = 5L,
    method = "zero_bias",
    workers = 1L
  )

  center_bits <- ubestarfm_unpack_candidates(result, row = 3L, col = 3L)
  expect_equal(center_bits$row, c(3L, 4L, 2L, 3L, 4L, 2L, 3L, 4L))
  expect_equal(center_bits$col, c(2L, 2L, 3L, 3L, 3L, 4L, 4L, 4L))
})

test_that("training is independent of all target dates", {
  model <- train_example_model(window_radius = 2L, patch_size = 20L)
  expect_false("coarse_target" %in% names(model))
})
```

Also test strict `< threshold`, edge clipping, missing reference cells, baseline/zero-bias membership equality, and deterministic results for `workers = 1L` versus `2L`.

- [ ] **Step 2: Run candidate tests and verify they fail because training is absent**

```bash
Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-candidate-masks.R")'
```

Expected: FAIL with missing `ubestarfm_train_arrays`.

- [ ] **Step 3: Implement packed-bit helpers and the C++ training kernel**

In `src/candidate_masks.cpp`, use C++17 and Rcpp. Implement:

```cpp
inline std::size_t cell_index(int row, int col, int ncol) {
  return static_cast<std::size_t>(row) * ncol + col;
}

inline int local_bit_index(int delta_row, int delta_col, int radius) {
  const int side = 2 * radius + 1;
  return (delta_col + radius) * side + (delta_row + radius);
}

inline void set_bit(Rcpp::RawVector& masks, std::size_t base, int bit) {
  const std::size_t byte = base + static_cast<std::size_t>(bit / 8);
  masks[byte] = static_cast<Rbyte>(masks[byte] | (1u << (bit % 8)));
}
```

Export `ubestarfm_train_patch_cpp()` with these inputs:

```cpp
Rcpp::List ubestarfm_train_patch_cpp(
    const Rcpp::NumericMatrix& fine_1,
    const Rcpp::NumericMatrix& fine_2,
    const Rcpp::NumericMatrix& coarse_1,
    const Rcpp::NumericMatrix& coarse_2,
    int row_start,
    int row_end,
    int col_start,
    int col_end,
    int window_radius,
    double threshold_1,
    double threshold_2);
```

Return zero-based `cell_ids` and a raw vector containing exactly `number_of_patch_cells * ceiling(side^2 / 8)` bytes. Set a bit only when all four reference values are finite and both fine differences are strictly below their patch thresholds.

- [ ] **Step 4: Implement deterministic patch creation and optional execution**

In `R/parallel.R`, implement `ubestarfm_lapply(x, fun, workers)`:

- Validate `workers >= 1`.
- Use `lapply()` for one worker.
- Use `parallel::mclapply()` on non-Windows systems.
- Use a PSOCK cluster with `on.exit(parallel::stopCluster(cluster), add = TRUE)` on Windows.
- Preserve input order.

In `R/train.R`, implement:

- `ubestarfm_make_patches(nrow, ncol, patch_size)`.
- R sample-SD thresholds with `stats::sd(..., na.rm = TRUE)`.
- Assembly of patch raw vectors into one fixed-size global raw vector.
- `ubestarfm_train_arrays()` for tests.
- Public `ubestarfm_train()` for raster/path inputs.

Model parameters must include `window_radius`, `patch_size`, `method`, `mask_side`, and `mask_bytes`.

- [ ] **Step 5: Generate Rcpp bindings and run focused tests**

```bash
Rscript -e 'Rcpp::compileAttributes(".")'
R CMD INSTALL .
Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-candidate-masks.R")'
Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-train.R")'
```

Expected: all training tests pass.

- [ ] **Step 6: Measure the bundled model and enforce the memory ceiling**

Add:

```r
test_that("bundled candidate-cache model stays below 80 MiB", {
  model <- train_example_model(window_radius = 25L, patch_size = 200L)
  expect_lt(as.numeric(object.size(model)), 80 * 1024^2)
  expect_type(model$candidate_masks, "raw")
  expect_equal(length(model$candidate_masks), 400L * 400L * 326L)
})
```

Run:

```bash
Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-train.R", reporter = "summary")'
```

Expected: PASS and model size below 80 MiB.

- [ ] **Step 7: Commit R training**

```bash
git add R src tests/testthat
git commit -m "feat(r): train compact reusable candidate models"
```

---

### Task 4: Implement R Single and Batch Prediction

**Files:**

- Create: `src/predict.cpp`
- Create/modify: `R/predict.R`
- Create/modify: `R/model.R`
- Create: `tests/testthat/test-predict-kernel.R`
- Create: `tests/testthat/test-predict-api.R`
- Create: `tests/testthat/test-model-persistence.R`

- [ ] **Step 1: Write failing prediction tests from a pure-R oracle**

In `tests/testthat/helper-fixtures.R`, add a small pure-R reference function that directly follows lines 87–203 of `legacy/ubESTARFM.R` for one pixel. It must calculate:

- Target-valid window indices.
- Detailed branch for more than five candidates.
- Zero-bias correction.
- Spectral/spatial and temporal weights.
- Abnormal-value fallback.
- Five-or-fewer-candidate fallback.

Use this independent oracle:

```r
reference_predict_pixel <- function(model, coarse_target, row, col, value_range) {
  radius <- model$parameters$window_radius
  nrow   <- model$metadata$nrow
  ncol   <- model$metadata$ncol
  rows   <- max(1L, row - radius):min(nrow, row + radius)
  cols   <- max(1L, col - radius):min(ncol, col + radius)

  fine_1   <- model$reference_values$fine_1
  fine_2   <- model$reference_values$fine_2
  coarse_1 <- model$reference_values$coarse_1
  coarse_2 <- model$reference_values$coarse_2
  target_valid <- model$reference_valid & is.finite(coarse_target)
  window_valid <- target_valid[rows, cols, drop = FALSE]

  candidates <- ubestarfm_unpack_candidates(model, row = row, col = col)
  keep <- is.finite(coarse_target[cbind(candidates$row, candidates$col)])
  candidates <- candidates[keep, , drop = FALSE]

  target_window  <- coarse_target[rows, cols, drop = FALSE][window_valid]
  coarse_1_window <- coarse_1[rows, cols, drop = FALSE][window_valid]
  coarse_2_window <- coarse_2[rows, cols, drop = FALSE][window_valid]

  if (nrow(candidates) > 5L) {
    index <- cbind(candidates$row, candidates$col)
    fine_candidates_1   <- fine_1[index]
    fine_candidates_2   <- fine_2[index]
    coarse_candidates_1 <- coarse_1[index]
    coarse_candidates_2 <- coarse_2[index]

    if (model$parameters$method == "zero_bias") {
      bias_1 <- -mean(fine_candidates_1) + mean(coarse_candidates_1)
      bias_2 <- -mean(fine_candidates_2) + mean(coarse_candidates_2)
      fine_candidates_1 <- fine_candidates_1 + bias_1
      fine_candidates_2 <- fine_candidates_2 + bias_2
      fine_pixel_1 <- fine_1[row, col] + bias_1
      fine_pixel_2 <- fine_2[row, col] + bias_2
    } else {
      fine_pixel_1 <- fine_1[row, col]
      fine_pixel_2 <- fine_2[row, col]
    }

    spectral_distance <- 1 - 0.5 * (
      abs(
        (fine_candidates_1 - coarse_candidates_1) /
          (fine_candidates_1 + coarse_candidates_1)
      ) +
        abs(
          (fine_candidates_2 - coarse_candidates_2) /
            (fine_candidates_2 + coarse_candidates_2)
        )
    )
    abnormal <- spectral_distance > 1 | spectral_distance < -1
    spectral_distance[abnormal] <- 0.5
    spatial_distance <- 1 + sqrt(
      (col - candidates$col)^2 + (row - candidates$row)^2
    ) / radius
    combined_distance <- (1 - spectral_distance) * spatial_distance + 1e-7
    weight <- (1 / combined_distance) / sum(1 / combined_distance)

    temporal_difference_1 <- abs(mean(target_window - coarse_1_window)) + 1e-10
    temporal_difference_2 <- abs(mean(target_window - coarse_2_window)) + 1e-10
    temporal_weight_1 <- (1 / temporal_difference_1) /
      (1 / temporal_difference_1 + 1 / temporal_difference_2)
    temporal_weight_2 <- 1 - temporal_weight_1

    target_candidates <- coarse_target[index]
    prediction_1 <- fine_pixel_1 + sum(
      weight * (target_candidates - coarse_candidates_1)
    )
    prediction_2 <- fine_pixel_2 + sum(
      weight * (target_candidates - coarse_candidates_2)
    )
    prediction <- temporal_weight_1 * prediction_1 +
      temporal_weight_2 * prediction_2

    if (prediction <= value_range[1] || prediction >= value_range[2]) {
      prediction <- temporal_weight_1 * sum(weight * fine_candidates_1) +
        temporal_weight_2 * sum(weight * fine_candidates_2)
    }
    return(prediction)
  }

  if (model$parameters$method == "zero_bias") {
    fine_pixel_1 <- fine_1[row, col] -
      mean(fine_1[rows, cols, drop = FALSE][window_valid]) +
      mean(coarse_1_window)
    fine_pixel_2 <- fine_2[row, col] -
      mean(fine_2[rows, cols, drop = FALSE][window_valid]) +
      mean(coarse_2_window)
  } else {
    fine_pixel_1 <- fine_1[row, col]
    fine_pixel_2 <- fine_2[row, col]
  }

  temporal_difference_1 <- mean(target_window - coarse_1_window) + 1e-10
  temporal_difference_2 <- mean(target_window - coarse_2_window) + 1e-10
  temporal_weight_1 <- (1 / abs(temporal_difference_1)) /
    (1 / abs(temporal_difference_1) + 1 / abs(temporal_difference_2))
  temporal_weight_2 <- 1 - temporal_weight_1
  temporal_weight_1 * (fine_pixel_1 + temporal_difference_1) +
    temporal_weight_2 * (fine_pixel_2 + temporal_difference_2)
}
```

Add these complete tests:

```r
test_that("detailed prediction matches the published equations", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 0.5
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
  actual <- ubestarfm_predict_arrays(model, list(target), c(-Inf, Inf))[[1L]]
  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})

test_that("fallback prediction matches the published equations", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 0.5
  target[2:4, 2:4] <- NA_real_
  target[3, 3] <- fixture$coarse_1[3, 3] + 0.5
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
  actual <- ubestarfm_predict_arrays(model, list(target), c(-Inf, Inf))[[1L]]
  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})

test_that("target missing values filter and renormalize candidates", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 0.5
  target[4, 4] <- NA_real_
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
  actual <- ubestarfm_predict_arrays(model, list(target), c(-Inf, Inf))[[1L]]
  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})

test_that("single and batch prediction are identical", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  first <- fixture$coarse_1 + 0.5
  second <- fixture$coarse_1 + 1
  batch <- ubestarfm_predict_arrays(model, list(first, second), c(-Inf, Inf))
  expect_equal(
    batch[[1L]],
    ubestarfm_predict_arrays(model, list(first), c(-Inf, Inf))[[1L]]
  )
})

test_that("baseline and zero-bias modes match their independent oracles", {
  fixture <- small_reference_fixture()
  target <- fixture$coarse_1 + 0.5
  for (method in c("baseline", "zero_bias")) {
    model <- do.call(
      ubestarfm_train_arrays,
      c(
        fixture,
        list(
          window_radius = 1L,
          patch_size = 5L,
          method = method
        )
      )
    )
    expected <- reference_predict_pixel(model, target, 3L, 3L, c(-Inf, Inf))
    actual <- ubestarfm_predict_arrays(model, list(target), c(-Inf, Inf))[[1L]]
    expect_equal(actual[3, 3], expected, tolerance = 1e-10)
  }
})

test_that("value-range correction uses corrected fine candidates", {
  fixture <- small_reference_fixture()
  model <- do.call(
    ubestarfm_train_arrays,
    c(fixture, list(window_radius = 1L, patch_size = 5L))
  )
  target <- fixture$coarse_1 + 1000
  expected <- reference_predict_pixel(model, target, 3L, 3L, c(0, 1))
  actual <- ubestarfm_predict_arrays(model, list(target), c(0, 1))[[1L]]
  expect_equal(actual[3, 3], expected, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run focused tests and verify missing prediction failures**

```bash
Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-predict-kernel.R")'
```

Expected: FAIL because the prediction kernel is absent.

- [ ] **Step 3: Implement C++ mask decoding and prediction**

In `src/predict.cpp`, implement:

```cpp
inline bool test_bit(
    const Rcpp::RawVector& masks,
    std::size_t base,
    int bit);

Rcpp::NumericMatrix ubestarfm_predict_patch_cpp(
    const Rcpp::NumericMatrix& fine_1,
    const Rcpp::NumericMatrix& fine_2,
    const Rcpp::NumericMatrix& coarse_1,
    const Rcpp::NumericMatrix& coarse_2,
    const Rcpp::RawVector& candidate_masks,
    int mask_bytes,
    const Rcpp::NumericVector& targets,
    int target_count,
    int row_start,
    int row_end,
    int col_start,
    int col_end,
    int window_radius,
    std::string method,
    double value_min,
    double value_max);
```

`targets` is a contiguous row-major array with layout `[target, row, col]`. The returned matrix has one row per patch cell and one column per target.

Use this exact local candidate order:

```cpp
for (int delta_col = -radius; delta_col <= radius; ++delta_col) {
  for (int delta_row = -radius; delta_row <= radius; ++delta_row) {
    const int bit = local_bit_index(delta_row, delta_col, radius);
    // Decode the bit, reject out-of-grid cells, and append finite candidates.
  }
}
```

For each target, intersect candidate bits with finite target values. Compute the temporal-window means from all reference-valid and target-valid window cells, not only candidates. Add `1e-10`, matching `0.01^5`, before temporal inversion.

- [ ] **Step 4: Implement R prediction APIs and atomic output**

In `R/predict.R`, implement:

- `ubestarfm_predict_arrays()` for tests.
- `ubestarfm_predict()`.
- `ubestarfm_predict_batch()`.
- `ubestarfm_fit_predict_batch()` for `cache = "none"`.
- Target geometry validation.
- `value_range` validation.
- Temporary-file write in the destination directory followed by `file.rename()`.
- Output existence checks unless `overwrite = TRUE`.

Batch prediction must stack target matrices into one contiguous row-major numeric vector and invoke patch prediction once per patch. It must not call `ubestarfm_train()` internally.

- [ ] **Step 5: Implement model persistence**

In `R/model.R`:

```r
ubestarfm_save_model <- function(model, path, compress = TRUE) {
  ubestarfm_validate_model(model)
  saveRDS(model, path, compress = compress)
  invisible(normalizePath(path, mustWork = TRUE))
}

ubestarfm_load_model <- function(path) {
  model <- readRDS(path)
  ubestarfm_validate_model(model)
  model
}
```

Test compressed and uncompressed round trips and reject unsupported schema versions.

- [ ] **Step 6: Run prediction, persistence, and package tests**

```bash
Rscript -e 'Rcpp::compileAttributes(".")'
R CMD INSTALL .
Rscript -e 'testthat::test_local(".", reporter = "summary")'
```

Expected: all tests pass.

- [ ] **Step 7: Commit R prediction**

```bash
git add R src tests/testthat
git commit -m "feat(r): predict single and batch target dates"
```

---

### Task 5: Add the R Compatibility Wrapper and Legacy Golden Validation

**Files:**

- Create: `R/compat.R`
- Create: `tests/testthat/test-compat.R`
- Create: `tests/testthat/test-legacy-golden.R`
- Create: `docs/migration/from-published-r-api.md`

- [ ] **Step 1: Write failing compatibility tests**

Test that `ubESTARFM()`:

- Accepts the published argument names and `"zero bias"` spelling.
- Writes `out_path`.
- Returns a `SpatRaster`.
- Ignores unrelated `.tif` files in `tmp_path`.
- Produces the same values as explicit `ubestarfm_train()` plus `ubestarfm_predict()`.

Add a golden test comparing the modern R output to `tests/fixtures/legacy/fused_result.tif` with maximum finite-pixel absolute difference `<= 1e-4` K and identical missing mask and geometry.

- [ ] **Step 2: Run the tests and verify wrapper failures**

```bash
Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-compat.R")'
```

Expected: FAIL because `ubESTARFM()` is absent from the package.

- [ ] **Step 3: Implement the compatibility wrapper**

Implement the approved signature exactly:

```r
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
```

Retain `tmp_path` only for source compatibility. Document that it is intentionally unused.

- [ ] **Step 4: Diagnose any golden mismatch before changing tolerance**

Run:

```bash
Rscript -e 'library(ubESTARFM); testthat::test_file("tests/testthat/test-legacy-golden.R", reporter = "progress")'
```

If it fails, write a targeted test for the first differing pixel and compare candidate order, patch threshold, target-valid window, and summation order. Do not raise the `1e-4` tolerance without documenting a scientifically justified difference in the approved design and migration guide.

- [ ] **Step 5: Run complete R verification**

```bash
R CMD build .
R CMD check --no-manual ubESTARFM_3.0.0.9000.tar.gz
Rscript -e 'testthat::test_local(".", reporter = "summary")'
```

Expected: `R CMD check` has no ERROR or WARNING, and all tests pass.

- [ ] **Step 6: Commit the compatibility layer**

```bash
git add R tests/testthat docs/migration
git commit -m "feat(r): preserve published single-target API"
```

---

### Task 6: Establish the Python Package, Raster I/O, and Model Contract

**Files:**

- Create: `python/pyproject.toml`
- Create: `python/README.md`
- Create: `python/src/ubestarfm/__init__.py`
- Create: `python/src/ubestarfm/io.py`
- Create: `python/src/ubestarfm/model.py`
- Create: `python/tests/conftest.py`
- Create: `python/tests/test_io.py`
- Create: `python/tests/test_model.py`

- [ ] **Step 1: Write failing Python I/O and model tests**

Create tests equivalent to the R contract:

```python
def test_read_raster_returns_float64_row_major_array(example_path):
    raster = read_raster(example_path("MOD11A1_LST_cloudrm_20160218.tif"))
    assert raster.values.shape == (400, 400)
    assert raster.values.dtype == np.float64
    assert raster.values.flags.c_contiguous
    assert raster.values[0, 0] == pytest.approx(312.74, abs=1e-4)


def test_geometry_mismatch_is_rejected():
    first = RasterData(
        values=np.ones((3, 3), dtype=np.float64),
        profile={
            "height": 3,
            "width": 3,
            "transform": from_origin(0, 3, 1, 1),
            "crs": CRS.from_epsg(4326),
            "nodata": np.nan,
            "dtype": "float64",
            "count": 1,
            "driver": "GTiff",
        },
    )
    second = RasterData(
        values=np.ones((4, 3), dtype=np.float64),
        profile={
            **first.profile,
            "height": 4,
            "transform": from_origin(0, 4, 1, 1),
        },
    )
    with pytest.raises(ValueError, match="same dimensions, transform, and CRS"):
        validate_geometry([first, second])
```

Import `CRS` from `rasterio.crs`, `from_origin` from
`rasterio.transform`, and `RasterData`, `read_raster`, and
`validate_geometry` from `ubestarfm.io`. Add model tests that construct a
valid 2 by 2 `UBESTARFMModel`, replace `schema_version` with `99`, and assert
`ValueError("Unsupported model schema")`; assigning to a frozen field must
raise `dataclasses.FrozenInstanceError`.

- [ ] **Step 2: Run tests and verify import failures**

```bash
python3 -m pytest python/tests/test_io.py python/tests/test_model.py -v
```

Expected: FAIL because `ubestarfm` is not installed.

- [ ] **Step 3: Create Python metadata**

Create `python/pyproject.toml`:

```toml
[build-system]
requires = ["hatchling>=1.27"]
build-backend = "hatchling.build"

[project]
name = "ubestarfm"
version = "3.0.0.dev0"
description = "Reusable unbiased ESTARFM training and prediction"
readme = "README.md"
requires-python = ">=3.12"
license = "MIT"
authors = [{ name = "Yi Yu", email = "yi.yu1@anu.edu.au" }]
dependencies = [
  "numpy>=2.0",
  "rasterio>=1.5",
  "numba>=0.65",
]

[project.optional-dependencies]
test = ["pytest>=9", "ruff"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-ra"

[tool.hatch.build.targets.wheel]
packages = ["src/ubestarfm"]

[tool.ruff]
line-length = 100
target-version = "py312"
```

Create `python/README.md` with the package name, a two-sentence summary, and
links to the root README and `docs/tutorials/python-single-and-batch.md`.

Before adding module code, add compliant docstring headers to every new `.py`
file under `python/src/ubestarfm/` and `python/tests/`. Package modules omit the
shebang; directly executable scripts added in later tasks include it.

- [ ] **Step 4: Implement typed raster and model objects**

In `io.py`, define frozen `RasterData(values, profile)` and support paths, open Rasterio datasets, and arrays with explicit profiles. Convert masked/nodata values to `np.nan`.

In `model.py`, define:

```python
@dataclass(frozen=True, slots=True)
class UBESTARFMModel:
    schema_version: int
    package_version: str
    reference_values: np.ndarray
    reference_valid: np.ndarray
    candidate_masks: np.ndarray
    patch_ids: np.ndarray
    patch_thresholds: np.ndarray
    profile: dict[str, object]
    window_radius: int
    patch_size: int
    method: str
```

Validate dtypes, C contiguity, shapes, mask byte count, schema version, and method in `__post_init__`.

- [ ] **Step 5: Install editable package and run tests**

```bash
python3 -m pip install -e 'python[test]'
python3 -m pytest python/tests/test_io.py python/tests/test_model.py -v
python3 -m ruff check python
```

Expected: all tests and lint pass.

- [ ] **Step 6: Commit Python foundation**

```bash
git add python
git commit -m "feat(python): establish package and raster model contracts"
```

---

### Task 7: Implement Python Numba Training and Prediction

**Files:**

- Create: `python/src/ubestarfm/kernels.py`
- Create: `python/src/ubestarfm/api.py`
- Modify: `python/src/ubestarfm/model.py`
- Modify: `python/src/ubestarfm/__init__.py`
- Create: `python/tests/test_candidate_masks.py`
- Create: `python/tests/test_predict_kernel.py`
- Create: `python/tests/test_api.py`
- Create: `python/tests/test_serialization.py`

- [ ] **Step 1: Write failing tests equivalent to R**

Port the small-array R fixtures and assertions exactly:

- Candidate bit positions and ordering.
- Strict threshold comparison.
- Detailed and fallback equations.
- Missing-target filtering.
- Baseline and zero-bias behavior.
- Single/batch equality.
- One/two-worker equality.
- Save/load equality.
- Model memory below 80 MiB on bundled inputs.

- [ ] **Step 2: Run tests and verify missing-kernel failures**

```bash
python3 -m pytest \
  python/tests/test_candidate_masks.py \
  python/tests/test_predict_kernel.py -v
```

Expected: FAIL because training and prediction functions are absent.

- [ ] **Step 3: Implement Numba candidate-mask training**

In `kernels.py`, implement `@numba.njit(cache=True, nogil=True)` helpers:

```python
def _cell_index(row: int, col: int, ncols: int) -> int:
    return row * ncols + col


def _local_bit_index(delta_row: int, delta_col: int, radius: int) -> int:
    side = 2 * radius + 1
    return (delta_col + radius) * side + (delta_row + radius)
```

Implement `train_patch()` returning `uint8` masks for the specified patch. Use `np.std(..., ddof=1)` thresholds calculated outside the kernel to match R. `train()` assembles deterministic patches and uses `ThreadPoolExecutor` only when `workers > 1`; Numba kernels must release the GIL.

- [ ] **Step 4: Implement Numba prediction**

Implement `predict_patch()` with the same loop nesting, finite checks, temporal epsilon, and branch logic as the Rcpp kernel. Accept targets with shape `(target_count, rows, cols)` and return `(target_count, patch_rows, patch_cols)`.

The API layer implements:

- `train()`
- `predict()`
- `predict_batch()`
- `fit_predict_batch()`
- `save_model()`
- `load_model()`

Serialization uses one compressed `.npz` containing arrays and a JSON string
scalar named `metadata_json`. Convert the affine transform to its first six
coefficients and the CRS to WKT before JSON encoding; reconstruct Rasterio
objects when loading. Write through a destination-directory temporary file and
atomically replace.

- [ ] **Step 5: Run Python tests and enforce memory**

```bash
python3 -m pytest python/tests -v
python3 -m ruff check python
```

Expected: all tests pass. The model arrays’ `nbytes` sum is below 80 MiB.

- [ ] **Step 6: Commit Python engine**

```bash
git add python
git commit -m "feat(python): add reusable Numba training and prediction"
```

---

### Task 8: Add Cross-Language and Golden-Output Acceptance

**Files:**

- Create: `tests/fixtures/small/reference_arrays.csv`
- Create: `tests/cross_language/run_r_fixture.R`
- Create: `tests/cross_language/compare_outputs.py`
- Create: `tests/cross_language/test_cross_language.py`
- Create: `tests/cross_language/README.md`

- [ ] **Step 1: Write the failing cross-language test**

`test_cross_language.py` must:

1. Run `Rscript tests/cross_language/run_r_fixture.R <tempdir>`.
2. Run Python train/batch prediction on the same bundled inputs.
3. Compare modern R, modern Python, and the legacy golden GeoTIFF.
4. Assert geometry and masks exactly.
5. Assert maximum finite difference `<= 1e-4`.

- [ ] **Step 2: Run and verify expected failure before runner implementation**

```bash
python3 -m pytest tests/cross_language/test_cross_language.py -v
```

Expected: FAIL because the R runner and comparison outputs do not exist.

- [ ] **Step 3: Implement the R runner with a compliant script header**

`run_r_fixture.R` must train once, predict the bundled February 18 target, save the model, reload it, and predict again. It writes:

- `r_prediction.tif`
- `r_prediction_reloaded.tif`
- `r_model_size_bytes.txt`

Use `on.exit()` for temporary cleanup and no `setwd()`.

- [ ] **Step 4: Implement Python comparison diagnostics**

`compare_outputs.py` reports:

- Pixel count.
- Missing-mask differences.
- Maximum and mean absolute differences.
- Location and values of the largest difference.
- Transform and CRS equality.

Exit nonzero when tolerance or geometry checks fail.

- [ ] **Step 5: Run full cross-language acceptance**

```bash
python3 -m pytest tests/cross_language/test_cross_language.py -v
Rscript -e 'testthat::test_local(".", reporter = "summary")'
python3 -m pytest python/tests tests/layout tests/cross_language -v
```

Expected: all checks pass.

- [ ] **Step 6: Commit acceptance fixtures**

```bash
git add tests
git commit -m "test: enforce R Python and published-output parity"
```

---

### Task 9: Add Examples, Tutorials, Migration Guidance, and Benchmarks

**Files:**

- Create: `examples/outputs/README.md`
- Create: `examples/R/single_target.R`
- Create: `examples/R/batch_targets.R`
- Create: `examples/python/single_target.py`
- Create: `examples/python/batch_targets.py`
- Create: `docs/tutorials/r-single-and-batch.md`
- Create: `docs/tutorials/python-single-and-batch.md`
- Create: `docs/tutorials/model-caching.md`
- Create: `benchmarks/benchmark.R`
- Create: `benchmarks/benchmark.py`
- Rewrite: `README.md`
- Create: `CHANGELOG.md`
- Create: `CONTRIBUTING.md`

- [ ] **Step 1: Write smoke tests for all examples**

Create a test that runs each example with a temporary output directory supplied through `--output-dir`. Assert expected files exist and can be opened.

Run:

```bash
python3 -m pytest tests/examples/test_examples.py -v
```

Expected: FAIL because modern examples do not exist.

- [ ] **Step 2: Implement R examples with standard headers**

Each R script starts with:

```r
#!/usr/bin/env Rscript
# Script: <name>.R
# Objective: Demonstrate <single or batch> ubESTARFM prediction.
# Author: Yi Yu
# Created: 2026-06-21
# Last updated: 2026-06-21
# Inputs: Bundled GeoTIFFs under inst/extdata and --output-dir.
# Outputs: Fused GeoTIFFs and an optional saved model.
# Usage: Rscript examples/R/<name>.R --output-dir examples/outputs
# Dependencies: R package ubESTARFM.
```

The batch example uses the bundled February 18 target twice under distinct labels only to demonstrate API equivalence until a second target fixture is added. It clearly states this in output and documentation; tests verify the two outputs match.

- [ ] **Step 3: Implement Python examples with standard headers**

Use the corresponding required Python docstring header, `argparse`, project-root paths derived from `Path(__file__)`, and no current-directory assumptions.

- [ ] **Step 4: Implement benchmarks**

Benchmarks report JSON and Markdown containing:

- R training time and model size.
- R warm single-target prediction time.
- R two-target batch time.
- Python equivalents.
- R/Python maximum output difference.
- Legacy timing only when `raster`, `foreach`, and `doParallel` are installed.

Benchmark scripts must include standard headers and place generated results under ignored `benchmarks/results/`.

- [ ] **Step 5: Rewrite user documentation**

README order:

1. Scientific purpose and citation.
2. R installation and quick start.
3. Python installation and quick start.
4. Train-once/batch-predict explanation.
5. Model-cache memory tradeoff.
6. Links to tutorials, benchmarks, migration, and legacy code.

Use the updated social handle already present in the user’s dirty working tree only if explicitly carried into this branch during integration; do not infer or overwrite unrelated local changes.

- [ ] **Step 6: Run example and documentation verification**

```bash
python3 -m pytest tests/examples/test_examples.py -v
Rscript examples/R/single_target.R --output-dir "$(mktemp -d)"
Rscript examples/R/batch_targets.R --output-dir "$(mktemp -d)"
python3 examples/python/single_target.py --output-dir "$(mktemp -d)"
python3 examples/python/batch_targets.py --output-dir "$(mktemp -d)"
```

Expected: all commands exit 0.

- [ ] **Step 7: Commit documentation and examples**

```bash
git add README.md CHANGELOG.md CONTRIBUTING.md examples docs/tutorials docs/migration benchmarks tests/examples
git commit -m "docs: add R and Python tutorials and benchmarks"
```

---

### Task 10: Add CI and Perform Final Repository Verification

**Files:**

- Create: `.github/workflows/r-check.yaml`
- Create: `.github/workflows/python-test.yaml`
- Create: `.github/workflows/cross-language.yaml`
- Modify as required: package metadata, tests, documentation

- [ ] **Step 1: Add failing CI-configuration checks**

Create `tests/layout/test_ci.py` asserting:

- Three workflow files exist.
- Workflow YAML contains R check, Python tests, cross-language tests, and archival syntax checks.
- Python matrices include 3.12 and 3.13.
- R workflow uses release R.

Run:

```bash
python3 -m pytest tests/layout/test_ci.py -v
```

Expected: FAIL because workflows do not exist.

- [ ] **Step 2: Add R workflow**

Use `r-lib/actions` to:

- Set up release R.
- Install system dependencies and R dependencies.
- Run `R CMD check --no-manual`.
- Parse all archival R scripts.
- Verify legacy checksums.

- [ ] **Step 3: Add Python workflow**

Use Python 3.12 and 3.13 to:

- Install `python[test]`.
- Run Ruff.
- Run Python unit tests and layout tests.
- Compile the archival Python script.

- [ ] **Step 4: Add cross-language workflow**

On Ubuntu:

- Install release R and Python 3.12.
- Install the R and Python packages.
- Run the shared cross-language test.
- Upload comparison diagnostics on failure.

- [ ] **Step 5: Run fresh full verification**

Run all commands from a clean shell:

```bash
git diff --check
sha256sum -c tests/fixtures/legacy/SHA256SUMS
Rscript -e 'files <- list.files("legacy", pattern = "[.]R$", recursive = TRUE, full.names = TRUE); for (f in files) parse(file = f)'
python3 -m py_compile legacy/paper-processing/02_process_landsat_lst.py
Rscript -e 'Rcpp::compileAttributes(".")'
Rscript -e 'roxygen2::roxygenise(".")'
R CMD build .
R CMD check --no-manual ubESTARFM_3.0.0.9000.tar.gz
python3 -m pip install -e 'python[test]'
python3 -m ruff check python tests
python3 -m pytest python/tests tests -v
Rscript -e 'testthat::test_local(".", reporter = "summary")'
```

Expected:

- No checksum failures.
- No archival syntax failures.
- `R CMD check`: 0 ERROR, 0 WARNING.
- Ruff: clean.
- All R, Python, layout, example, and cross-language tests pass.

- [ ] **Step 6: Generate final benchmark evidence**

```bash
mkdir -p benchmarks/results
Rscript benchmarks/benchmark.R --output-dir benchmarks/results
python3 benchmarks/benchmark.py --output-dir benchmarks/results
```

Review:

- Candidate search occurs once per trained model.
- Single and batch outputs match.
- R/Python maximum difference is within `1e-4`.
- In-memory model size is below 80 MiB.

- [ ] **Step 7: Request code review before completion**

Invoke `superpowers:requesting-code-review`. Address actionable findings with `superpowers:receiving-code-review`, writing a failing test before each behavioral correction.

- [ ] **Step 8: Commit CI and final quality fixes**

```bash
git add .github tests DESCRIPTION NAMESPACE R src python README.md CHANGELOG.md CONTRIBUTING.md
git commit -m "ci: verify modern R and Python ubESTARFM engines"
```

- [ ] **Step 9: Re-run final verification and inspect Git state**

Repeat Step 5 after the final commit, then run:

```bash
git status --short --branch
git log --oneline --decorate --max-count=15
git tag --list 'legacy-layout-v2.0.1' --format='%(refname:short) %(objectname:short)'
```

Expected: clean feature branch, focused commit history, and the local preservation tag pointing to `4721ac0`.

---

## Execution Notes

- Do not copy the untracked `0_algorithm/ubESTARFM_fast.R` from the user’s original dirty worktree.
- Do not alter `legacy/ubESTARFM.R`; checksum verification is mandatory after every repository-wide formatting operation.
- Do not use the deprecated `raster`, `foreach`, or `doParallel` packages in maintained code.
- Do not run formatters over `legacy/`.
- Keep R assignment and formatting consistent with the workspace R coding standard.
- Use Arial explicitly in maintained plotting code.
- Treat any golden-output mismatch as a debugging task, not a reason to weaken tolerance.
- Preserve unrelated user changes when this branch is eventually integrated.
