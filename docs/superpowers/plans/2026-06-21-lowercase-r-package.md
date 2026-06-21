# Lowercase R Package Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rename the maintained R package to `ubestarfm`, document direct GitHub installation, and retain the ubESTARFM algorithm name, published compatibility function, repository name, and legacy artifacts.

**Architecture:** Treat the lowercase string as package identity rather than a global textual rename. Change package metadata first, regenerate Rcpp and roxygen outputs from their sources, then update package consumers and documentation while regression tests enforce the boundary between lowercase package identity and mixed-case scientific provenance.

**Tech Stack:** R 4.3+, Rcpp, roxygen2, testthat, terra, GitHub Actions, Python 3.12+, pytest.

---

### Task 1: Add a package-identity regression guard

**Files:**
- Create: `tests/layout/test_r_package_identity.py`

- [ ] **Step 1: Add a test that defines the lowercase package-name contract**

Create `tests/layout/test_r_package_identity.py` with:

```python
"""
Script: test_r_package_identity.py
Objective: Verify lowercase R package identity while preserving ubESTARFM provenance.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: R package metadata, generated bindings, tests, examples, CI, and README.
Outputs: Pytest assertions for package naming and compatibility boundaries.
Usage: python3 -m pytest tests/layout/test_r_package_identity.py -v
Dependencies: Python standard library and pytest.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

PACKAGE_IDENTITY_FILES = (
    "DESCRIPTION",
    "NAMESPACE",
    "src/RcppExports.cpp",
    "tests/testthat.R",
    ".github/workflows/r-check.yaml",
    "CONTRIBUTING.md",
)

PACKAGE_IDENTITY_GLOBS = (
    "R/*.R",
    "tests/testthat/*.R",
    "examples/R/*.R",
    "benchmarks/*.R",
    "tests/cross_language/*.R",
)

FORBIDDEN_PACKAGE_IDENTITIES = (
    "Package: ubESTARFM",
    "library(ubESTARFM)",
    "ubESTARFM:::",
    'package = "ubESTARFM"',
    'packageVersion("ubESTARFM")',
    'test_check("ubESTARFM")',
    "useDynLib(ubESTARFM",
    "_ubESTARFM_ubestarfm_",
    "R_init_ubESTARFM",
    "ubESTARFM_*.tar.gz",
    "ubESTARFM_3.0.0.9000.tar.gz",
)


def test_r_package_identity_is_lowercase() -> None:
    description = (ROOT / "DESCRIPTION").read_text()
    namespace = (ROOT / "NAMESPACE").read_text()

    assert description.startswith("Package: ubestarfm\n")
    assert "URL: https://github.com/yuyi13/ubESTARFM" in description
    assert "BugReports: https://github.com/yuyi13/ubESTARFM/issues" in description
    assert "useDynLib(ubestarfm, .registration = TRUE)" in namespace


def test_readme_uses_github_installation() -> None:
    readme = (ROOT / "README.md").read_text()

    assert 'remotes::install_github("yuyi13/ubESTARFM")' in readme
    assert "library(ubestarfm)" in readme


def test_maintained_files_do_not_use_the_old_package_identity() -> None:
    stale = {}
    paths = [ROOT / relative_path for relative_path in PACKAGE_IDENTITY_FILES]
    for pattern in PACKAGE_IDENTITY_GLOBS:
        paths.extend(ROOT.glob(pattern))

    for path in sorted(set(paths)):
        text = path.read_text()
        matches = [
            pattern
            for pattern in FORBIDDEN_PACKAGE_IDENTITIES
            if pattern in text
        ]
        if matches:
            stale[str(path.relative_to(ROOT))] = matches
    assert stale == {}


def test_published_names_and_artifacts_are_retained() -> None:
    namespace = (ROOT / "NAMESPACE").read_text()
    compatibility = (ROOT / "R" / "compat.R").read_text()

    assert "export(ubESTARFM)" in namespace
    assert "ubESTARFM <- function(" in compatibility
    assert (ROOT / "legacy" / "ubESTARFM.R").exists()
```

- [ ] **Step 2: Run the new test and confirm it fails against the old package identity**

Run:

```bash
python3 -m pytest tests/layout/test_r_package_identity.py -v
```

Expected: failures report `Package: ubESTARFM`, mixed-case native registration,
and missing GitHub-installation metadata.

- [ ] **Step 3: Commit the failing regression guard**

```bash
git add tests/layout/test_r_package_identity.py
git commit -m "test: define lowercase R package identity"
```

Expected: one new test file is committed while the test still fails.

### Task 2: Rename package metadata and regenerate native bindings

**Files:**
- Modify: `DESCRIPTION`
- Modify: `R/package.R`
- Modify: `R/model.R`
- Modify: `R/parallel.R`
- Regenerate: `R/RcppExports.R`
- Regenerate: `src/RcppExports.cpp`
- Regenerate: `NAMESPACE`
- Create: `man/ubestarfm-package.Rd`
- Delete: `man/ubESTARFM-package.Rd`

- [ ] **Step 1: Change the package metadata and add repository metadata**

Update the relevant `DESCRIPTION` fields to:

```text
Package: ubestarfm
Type: Package
Title: Unbiased Enhanced Spatial and Temporal Adaptive Reflectance Fusion
Version: 3.0.0.9000
Authors@R: person("Yi", "Yu", email = "yi.yu1@anu.edu.au", role = c("aut", "cre"))
Description: Trains reusable reference-pair models and predicts fine-resolution
    land surface temperature using the unbiased ESTARFM algorithm.
License: MIT + file LICENSE
URL: https://github.com/yuyi13/ubESTARFM
BugReports: https://github.com/yuyi13/ubESTARFM/issues
```

Leave the remaining encoding, dependency, roxygen, and system-requirement
fields unchanged.

- [ ] **Step 2: Change package-level source references**

In `R/package.R`, preserve the script header and use:

```r
# Objective: Register compiled code and package-level imports for ubestarfm.
# Usage: Internal package module; loaded through the ubestarfm namespace.

#' ubestarfm package
#'
#' @useDynLib ubestarfm, .registration = TRUE
```

In `R/model.R`, change the package-version lookup to:

```r
package_version = as.character(utils::packageVersion("ubestarfm"))
```

Retain the existing aligned assignment style around this line.

In `R/parallel.R`, change the Windows worker package load to:

```r
parallel::clusterEvalQ(cluster, library(ubestarfm))
```

Update the `Objective`, `Usage`, or `Dependencies` header fields in these files
where they refer to the old package identity. Do not lowercase the ubESTARFM
algorithm name in function documentation.

- [ ] **Step 3: Regenerate Rcpp bindings from the lowercase package metadata**

Run:

```bash
Rscript -e 'Rcpp::compileAttributes(".")'
```

Expected:

- `R/RcppExports.R` calls `_ubestarfm_ubestarfm_*` symbols.
- `src/RcppExports.cpp` defines `_ubestarfm_ubestarfm_*` wrappers.
- The native initializer is `R_init_ubestarfm`.

- [ ] **Step 4: Regenerate namespace and package documentation**

Run:

```bash
Rscript -e 'roxygen2::roxygenise(".")'
```

Expected:

- `NAMESPACE` contains `useDynLib(ubestarfm, .registration = TRUE)`.
- `man/ubestarfm-package.Rd` is created.

Delete the obsolete generated file:

```text
man/ubESTARFM-package.Rd
```

- [ ] **Step 5: Verify generated package identity**

Run:

```bash
rg -n \
  'Package:|useDynLib|_ubestarfm_ubestarfm_|R_init_ubestarfm|packageVersion|library\\(' \
  DESCRIPTION NAMESPACE R/package.R R/model.R R/parallel.R \
  R/RcppExports.R src/RcppExports.cpp man/ubestarfm-package.Rd
```

Expected: package and native identities are lowercase; `ubESTARFM()` is not
renamed.

- [ ] **Step 6: Commit metadata and generated bindings**

```bash
git add \
  DESCRIPTION NAMESPACE R/package.R R/model.R R/parallel.R \
  R/RcppExports.R src/RcppExports.cpp \
  man/ubestarfm-package.Rd man/ubESTARFM-package.Rd
git commit -m "refactor: rename R package to ubestarfm"
```

Expected: package metadata, generated bindings, and generated package
documentation are committed together.

### Task 3: Update tests, examples, benchmarks, and cross-language fixtures

**Files:**
- Modify: `tests/testthat.R`
- Modify: `tests/testthat/helper-fixtures.R`
- Modify: `tests/testthat/test-candidate-masks.R`
- Modify: `tests/testthat/test-model-persistence.R`
- Modify: `tests/testthat/test-predict-kernel.R`
- Modify: `tests/testthat/test-train.R`
- Modify: `tests/testthat/test-compat.R`
- Modify: `tests/testthat/test-io.R`
- Modify: `tests/testthat/test-legacy-golden.R`
- Modify: `tests/testthat/test-model.R`
- Modify: `tests/testthat/test-predict-api.R`
- Modify: `examples/R/single_target.R`
- Modify: `examples/R/batch_targets.R`
- Modify: `benchmarks/benchmark.R`
- Modify: `tests/cross_language/run_r_fixture.R`
- Modify: `R/compat.R`
- Modify: `R/io.R`
- Modify: `R/train.R`
- Modify: `R/predict.R`

- [ ] **Step 1: Update the test runner and installed-data lookup**

In `tests/testthat.R`, use:

```r
library(testthat)
library(ubestarfm)

test_check("ubestarfm")
```

Update its header fields to refer to the `ubestarfm` package identity.

In `tests/testthat/helper-fixtures.R`, use:

```r
package = "ubestarfm"
```

and:

```r
candidates <- ubestarfm:::ubestarfm_unpack_candidates(
```

- [ ] **Step 2: Update internal namespace access in R tests**

Replace package-qualified internal calls only:

```text
ubESTARFM:::ubestarfm_train_arrays
ubESTARFM:::ubestarfm_predict_arrays
ubESTARFM:::ubestarfm_unpack_candidates
```

with:

```text
ubestarfm:::ubestarfm_train_arrays
ubestarfm:::ubestarfm_predict_arrays
ubestarfm:::ubestarfm_unpack_candidates
```

Retain calls to the exported compatibility function `ubESTARFM(...)`.
Update script-header `Usage` and `Dependencies` fields from
`library(ubESTARFM)` or package `ubESTARFM` to lowercase package identity.

- [ ] **Step 3: Update maintained R entry points**

In both R example scripts, the R benchmark, and the cross-language fixture,
change:

```r
library(ubESTARFM)
```

to:

```r
library(ubestarfm)
```

Update their dependency header fields to `R package ubestarfm`. Preserve
scientific descriptions such as “R ubESTARFM prediction”.

- [ ] **Step 4: Update package-identity comments in maintained R modules**

In `R/compat.R`, `R/io.R`, `R/train.R`, and `R/predict.R`, update only header
text that says the module is loaded through or depends on the old package
namespace. Keep algorithm descriptions and the `ubESTARFM()` function unchanged.

- [ ] **Step 5: Verify core package consumers use the lowercase identity**

Run:

```bash
rg -n \
  'library\\(ubESTARFM\\)|ubESTARFM:::|package = "ubESTARFM"|test_check\\("ubESTARFM"\\)' \
  R tests/testthat.R tests/testthat examples/R benchmarks/benchmark.R \
  tests/cross_language/run_r_fixture.R
```

Expected: no matches. Documentation and CI are updated in Task 4.

- [ ] **Step 6: Run the R test suite under the lowercase package name**

Run:

```bash
Rscript -e 'testthat::test_local(reporter = "summary")'
```

Expected: the full R suite passes, including the unchanged `ubESTARFM()`
compatibility tests and the legacy golden-output test.

- [ ] **Step 7: Commit package consumers**

```bash
git add R tests/testthat.R tests/testthat examples/R benchmarks/benchmark.R \
  tests/cross_language/run_r_fixture.R
git commit -m "refactor: use lowercase R package identity"
```

Expected: maintained R consumers use `ubestarfm`; published function calls
remain mixed-case.

### Task 4: Update installation, migration, CI, and contributor documentation

**Files:**
- Modify: `README.md`
- Modify: `docs/tutorials/r-single-and-batch.md`
- Modify: `docs/migration/from-published-r-api.md`
- Modify: `CONTRIBUTING.md`
- Modify: `CHANGELOG.md`
- Modify: `.github/workflows/r-check.yaml`
- Modify: `.Rbuildignore`
- Modify: `tests/layout/test_ci.py`

- [ ] **Step 1: Make direct GitHub installation the primary README path**

Replace the R installation subsection in `README.md` with:

````markdown
Install the R package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("yuyi13/ubESTARFM")
library(ubestarfm)
```

The R package contains compiled C++ code. Installation from source therefore
requires an R build toolchain, such as
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows or the
Command Line Tools for Xcode on macOS.

For development from a local clone:

```bash
R CMD INSTALL .
```
````

Change the R quick-start code to:

```r
library(ubestarfm)
```

Keep the title, algorithm prose, `ubESTARFM()` compatibility section,
repository URL, and legacy paths mixed-case.

- [ ] **Step 2: Update the R tutorial and migration guide**

At the start of `docs/tutorials/r-single-and-batch.md`, use:

````markdown
Install directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("yuyi13/ubESTARFM")
library(ubestarfm)
```

For development from a local clone, run `R CMD INSTALL .`.
````

In `docs/migration/from-published-r-api.md`, add this package-name migration
before the API example:

````markdown
The maintained R package is named `ubestarfm`. Replace
`library(ubESTARFM)` with:

```r
library(ubestarfm)
```

If the earlier development package is installed, remove it in a clean R
session before installing the lowercase package:

```r
remove.packages("ubESTARFM")
install.packages("remotes")
remotes::install_github("yuyi13/ubESTARFM")
```

The published compatibility function remains `ubESTARFM()`.
````

- [ ] **Step 3: Update build and CI artifact names**

In `.Rbuildignore`, replace:

```text
^ubESTARFM_.*\.tar\.gz$
```

with:

```text
^ubestarfm_.*\.tar\.gz$
```

In `.github/workflows/r-check.yaml`, use:

```yaml
R CMD check --no-manual ubestarfm_*.tar.gz
```

In `CONTRIBUTING.md`, use:

```bash
R CMD check --no-manual ubestarfm_3.0.0.9000.tar.gz
```

Add a `CHANGELOG.md` bullet under 3.0.0:

```markdown
- Standardized the maintained R package identity as lowercase `ubestarfm`
  while retaining the published `ubESTARFM()` compatibility function.
```

- [ ] **Step 4: Strengthen the CI layout assertion**

In `tests/layout/test_ci.py`, update the script header’s `Last updated` date
only if it differs from 2026-06-21, then add this assertion inside
`test_r_workflow_covers_package_and_archive()`:

```python
assert "ubestarfm_*.tar.gz" in text
assert "ubESTARFM_*.tar.gz" not in text
```

- [ ] **Step 5: Run documentation and CI regression tests**

Run:

```bash
python3 -m pytest \
  tests/layout/test_r_package_identity.py \
  tests/layout/test_ci.py \
  -v
```

Expected: all package-identity and CI tests pass.

- [ ] **Step 6: Scan for stale package-identity forms**

Run:

```bash
rg -n \
  'library\\(ubESTARFM\\)|ubESTARFM:::|package = "ubESTARFM"|packageVersion\\("ubESTARFM"\\)|test_check\\("ubESTARFM"\\)|useDynLib\\(ubESTARFM|_ubESTARFM_ubestarfm_|R_init_ubESTARFM|ubESTARFM_\\*\\.tar\\.gz|ubESTARFM_3\\.0\\.0\\.9000\\.tar\\.gz' \
  DESCRIPTION NAMESPACE R tests/testthat.R tests/testthat examples/R \
  benchmarks tests/cross_language/run_r_fixture.R \
  .github/workflows/r-check.yaml CONTRIBUTING.md .Rbuildignore
```

Expected: no matches.

- [ ] **Step 7: Commit documentation and CI**

```bash
git add \
  README.md docs/tutorials/r-single-and-batch.md \
  docs/migration/from-published-r-api.md CONTRIBUTING.md CHANGELOG.md \
  .github/workflows/r-check.yaml .Rbuildignore tests/layout/test_ci.py
git commit -m "docs: add lowercase R GitHub installation"
```

Expected: installation and migration docs consistently use `ubestarfm`.

### Task 5: Verify installation, parallel execution, build, and cross-language behavior

**Files:**
- Verify: all changed package, test, and documentation files

- [ ] **Step 1: Confirm generated files are reproducible**

Run:

```bash
Rscript -e 'Rcpp::compileAttributes(".")'
Rscript -e 'roxygen2::roxygenise(".")'
git diff --exit-code
```

Expected: both generators complete and `git diff --exit-code` reports no
generated-file drift.

- [ ] **Step 2: Install and load the lowercase package in an isolated library**

Run:

```bash
test_library="${TMPDIR:-/tmp}/ubestarfm-test-library-$(id -u)"
user_library="$(Rscript -e 'cat(.libPaths()[1L])')"
rm -rf "$test_library"
mkdir -p "$test_library"
R CMD INSTALL --library="$test_library" .
R_LIBS_USER="$test_library:$user_library" Rscript -e \
  'library(ubestarfm); stopifnot("package:ubestarfm" %in% search())'
```

Expected: installation succeeds and `library(ubestarfm)` loads the installed
package.

- [ ] **Step 3: Verify a multi-worker prediction from the isolated install**

Run:

```bash
test_library="${TMPDIR:-/tmp}/ubestarfm-test-library-$(id -u)"
user_library="$(Rscript -e 'cat(.libPaths()[1L])')"
R_LIBS_USER="$test_library:$user_library" Rscript - <<'RS'
library(ubestarfm)

cluster <- parallel::makeCluster(2L)
loaded <- parallel::clusterEvalQ(cluster, {
  library(ubestarfm)
  "package:ubestarfm" %in% search()
})
parallel::stopCluster(cluster)
stopifnot(all(unlist(loaded)))

data_directory <- "inst/extdata"
model <- ubestarfm_train(
  fine_1 = file.path(data_directory, "Landsat_LST_cloudrm_20160205.tif"),
  fine_2 = file.path(data_directory, "Landsat_LST_cloudrm_20160308.tif"),
  coarse_1 = file.path(data_directory, "MOD11A1_LST_cloudrm_20160205.tif"),
  coarse_2 = file.path(data_directory, "MOD11A1_LST_cloudrm_20160308.tif"),
  window_radius = 2L,
  patch_size = 200L,
  workers = 2L
)
prediction <- ubestarfm_predict(
  model,
  file.path(data_directory, "MOD11A1_LST_cloudrm_20160218.tif"),
  workers = 2L
)
stopifnot(inherits(prediction, "SpatRaster"))
RS
```

Expected: training and prediction finish with two workers and return a
`SpatRaster`.

- [ ] **Step 4: Run the full R tests**

Run:

```bash
Rscript -e 'testthat::test_local(reporter = "summary")'
```

Expected: all R tests pass.

- [ ] **Step 5: Build and check the lowercase package tarball**

Run:

```bash
R CMD build .
R CMD check --no-manual ubestarfm_3.0.0.9000.tar.gz
```

Expected: `R CMD check` ends with `Status: OK`.

- [ ] **Step 6: Run the full Python and repository test suite**

Run:

```bash
test_library="${TMPDIR:-/tmp}/ubestarfm-test-library-$(id -u)"
user_library="$(Rscript -e 'cat(.libPaths()[1L])')"
rm -rf "$test_library"
mkdir -p "$test_library"
R CMD INSTALL --library="$test_library" ubestarfm_3.0.0.9000.tar.gz
R_LIBS_USER="$test_library:$user_library" PYTHONPATH=python/src \
  python3 -m pytest tests python/tests -q
```

Expected: all Python, layout, example, and cross-language tests pass.

- [ ] **Step 7: Verify preserved artifacts and README links**

Run:

```bash
sha256sum -c tests/fixtures/legacy/SHA256SUMS
python3 - <<'PY'
import re
from pathlib import Path

root = Path.cwd()
text = Path("README.md").read_text()
targets = re.findall(r"!?\[[^\]]*\]\(([^)]+)\)", text)
missing = []
for target in targets:
    local_target = target.split("#", 1)[0]
    if not local_target or local_target.startswith(("http://", "https://", "mailto:")):
        continue
    if not (root / local_target).exists():
        missing.append(local_target)

assert not missing, sorted(set(missing))
assert 'remotes::install_github("yuyi13/ubESTARFM")' in text
assert "library(ubestarfm)" in text
PY
```

Expected: checksums pass and all local README links resolve.

- [ ] **Step 8: Remove generated verification artifacts**

Remove only artifacts produced by this verification:

```bash
test_library="${TMPDIR:-/tmp}/ubestarfm-test-library-$(id -u)"
rm -rf ubestarfm.Rcheck ubestarfm_3.0.0.9000.tar.gz "$test_library"
```

Expected: source files and committed fixtures remain untouched.

- [ ] **Step 9: Run final repository checks**

Run:

```bash
git diff --check master...HEAD
git status --short --branch
git log --oneline master..HEAD
```

Expected:

- `git diff --check` produces no output.
- The worktree is clean on `chore/lowercase-r-package`.
- The branch contains the design, plan, regression-test, package-rename,
  consumer-update, and documentation commits.
