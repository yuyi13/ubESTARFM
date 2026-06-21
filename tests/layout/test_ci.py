"""
Script: test_ci.py
Objective: Verify CI coverage for R, Python, cross-language, and archival checks.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: GitHub Actions workflow YAML files.
Outputs: Pytest assertions for required CI jobs and version matrices.
Usage: python3 -m pytest tests/layout/test_ci.py -v
Dependencies: Python standard library and pytest.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
WORKFLOWS = ROOT / ".github" / "workflows"


def test_required_workflows_exist() -> None:
    assert (WORKFLOWS / "r-check.yaml").exists()
    assert (WORKFLOWS / "python-test.yaml").exists()
    assert (WORKFLOWS / "cross-language.yaml").exists()


def test_r_workflow_covers_package_and_archive() -> None:
    text = (WORKFLOWS / "r-check.yaml").read_text()
    assert "r-version: release" in text
    assert "R CMD check" in text
    assert "sha256sum -c" in text
    assert "list.files(\"legacy\"" in text


def test_python_workflow_covers_supported_versions_and_archive() -> None:
    text = (WORKFLOWS / "python-test.yaml").read_text()
    assert '"3.12"' in text
    assert '"3.13"' in text
    assert "ruff check" in text
    assert "py_compile" in text


def test_cross_language_workflow_runs_shared_acceptance() -> None:
    text = (WORKFLOWS / "cross-language.yaml").read_text()
    assert "test_cross_language.py" in text
