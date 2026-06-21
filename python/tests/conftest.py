"""
Script: conftest.py
Objective: Provide Python test paths and bundled raster fixtures.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Repository root and files under inst/extdata.
Outputs: Pytest fixtures and local package import configuration.
Usage: python3 -m pytest python/tests -v
Dependencies: Python standard library and pytest.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "python" / "src"
sys.path.insert(0, str(SOURCE))


@pytest.fixture
def example_path():
    """Return a function resolving a bundled raster filename."""

    def resolve(filename: str) -> Path:
        return ROOT / "inst" / "extdata" / filename

    return resolve


@pytest.fixture
def small_reference_fixture() -> dict[str, np.ndarray]:
    """Return a five-by-five reference-pair fixture."""
    fine_1 = np.arange(1, 26, dtype=np.float64).reshape(5, 5)
    fine_2 = fine_1 + 100
    coarse_1 = fine_1 + 200
    coarse_2 = fine_1 + 300
    coarse_1[1, 1] = np.nan
    return {
        "fine_1": fine_1,
        "fine_2": fine_2,
        "coarse_1": coarse_1,
        "coarse_2": coarse_2,
    }
