# Lowercase R Package Design

**Date:** 2026-06-21
**Status:** Pending written-spec review

## Purpose

Rename the maintained R package from `ubESTARFM` to `ubestarfm` so its package
identity matches the lowercase Python package, and allow users to install the R
package directly from GitHub without cloning the repository.

## Naming Boundary

Use `ubestarfm` for the maintained R package identity:

- `Package:` in `DESCRIPTION`.
- `library()` and namespace-qualified calls.
- Native-library and Rcpp registration names.
- Package lookup, version lookup, tests, examples, benchmarks, CI artifacts,
  package documentation, and dependency comments.

Retain `ubESTARFM` for the scientific algorithm and published provenance:

- Algorithm name in prose and citations.
- Published compatibility function `ubESTARFM()`.
- GitHub repository `yuyi13/ubESTARFM`.
- Archived file `legacy/ubESTARFM.R`.
- Historical paths, DOI records, titles, and figures.

The Python package remains `ubestarfm` and requires no API or package-identity
change.

## GitHub Installation

The primary R installation instructions will be:

```r
install.packages("remotes")
remotes::install_github("yuyi13/ubESTARFM")
library(ubestarfm)
```

The `remotes` package is an installation helper and will not become a runtime
dependency of `ubestarfm`. Because `ubestarfm` contains C++ code, the
documentation will note that source-package compilation requires a suitable R
build toolchain, such as Rtools on Windows or the command-line developer tools
on macOS.

Local installation from a clone will remain documented as a developer
alternative:

```bash
R CMD INSTALL .
```

## Package Metadata and Native Code

The package rename will update:

- `DESCRIPTION`, including repository and issue-tracker metadata.
- The package-level roxygen declaration and generated package documentation.
- `NAMESPACE` and `useDynLib()` registration.
- Rcpp-generated R and C++ registration symbols.
- Internal calls to `packageVersion()`, `system.file()`, and package loading on
  parallel workers.

Rcpp bindings will be regenerated with `Rcpp::compileAttributes()` after the
package metadata changes. Documentation and namespace files will be regenerated
with `roxygen2` where their source is roxygen-managed.

## Compatibility and Migration

The exported `ubESTARFM()` compatibility function remains unchanged, so code
using the published function can continue to call it after loading the new
lowercase package:

```r
library(ubestarfm)
result <- ubESTARFM(...)
```

The package-name change means existing `library(ubESTARFM)` calls must become
`library(ubestarfm)`. Migration documentation will recommend removing an
installed development copy named `ubESTARFM`, restarting R, and then installing
the lowercase package:

```r
remove.packages("ubESTARFM")
install.packages("remotes")
remotes::install_github("yuyi13/ubESTARFM")
```

Saved `ubestarfm_model` objects remain structurally compatible. Their embedded
package-version metadata is informational and does not depend on the old
package name when loading or predicting.

## Repository Updates

Update all package-identity references in:

- R source and generated Rcpp bindings.
- R tests and shared fixtures.
- Maintained R examples and benchmark scripts.
- Cross-language R fixtures.
- GitHub Actions R build/check commands.
- README installation and R usage sections.
- R tutorial, migration guide, contribution guide, and generated manual pages.

Do not alter the published algorithm file or its checksum. Do not mechanically
lowercase scientific prose, the compatibility function, the repository name,
or legacy paths.

## Validation

The completed rename must pass:

1. A repository scan that permits `ubESTARFM` only for the algorithm,
   compatibility function, repository, and legacy provenance.
2. `Rcpp::compileAttributes()` and roxygen generation without uncommitted
   generated-file drift.
3. Installation of the package under the name `ubestarfm`.
4. Loading with `library(ubestarfm)`.
5. A multi-worker prediction to verify parallel worker package loading.
6. The full R test suite.
7. `R CMD build` and `R CMD check --no-manual` using the lowercase tarball
   name.
8. The full Python, layout, example, and cross-language test suite.
9. The published algorithm and golden-output checksum checks.
10. README local-link and installation-snippet validation.

## Scope

This change does not rename the GitHub repository, algorithm, compatibility
function, Python API, model class, public `ubestarfm_*` functions, legacy
artifacts, or DOI records. It does not publish a release or push changes to
GitHub.
