"""
Script: test_layout.py
Objective: Verify the modern repository layout and preserve published artifacts.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Repository files and published artifact checksums.
Outputs: Pytest assertions for layout and artifact integrity.
Usage: python3 -m pytest tests/layout/test_layout.py -v
Dependencies: Python standard library and pytest.
"""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    """Return the SHA-256 digest for a file."""
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
    numbered_paths = (
        "0_algorithm",
        "1_test_data",
        "2_tmp_path",
        "3_output",
        "4_lst_processing_scripts",
        "5_ozflux_lst",
    )
    assert not any((ROOT / name).exists() for name in numbered_paths)


def test_published_algorithm_is_byte_for_byte_preserved() -> None:
    assert sha256(ROOT / "legacy" / "ubESTARFM.R") == (
        "131f15efab127deccbe411cf800ea7967e621e8c916de677430299448f47e67b"
    )


def test_legacy_golden_output_is_preserved() -> None:
    assert sha256(ROOT / "tests" / "fixtures" / "legacy" / "fused_result.tif") == (
        "c4d216842fe9624f5355feac88ba6cb1b8dd35afd0e66272578561bc2c92a490"
    )
