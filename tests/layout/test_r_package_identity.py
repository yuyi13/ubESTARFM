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
