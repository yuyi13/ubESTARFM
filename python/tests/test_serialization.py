"""
Script: test_serialization.py
Objective: Verify Python model serialization and schema validation.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Small reusable models and temporary NPZ paths.
Outputs: Pytest assertions for save and load behavior.
Usage: python3 -m pytest python/tests/test_serialization.py -v
Dependencies: numpy, pytest, and ubestarfm.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from ubestarfm import load_model, save_model
from ubestarfm.api import train_arrays


def test_model_round_trip(small_reference_fixture, tmp_path) -> None:
    model = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
    )
    path = tmp_path / "model.npz"
    save_model(model, path)
    loaded = load_model(path)

    np.testing.assert_array_equal(loaded.candidate_masks, model.candidate_masks)
    np.testing.assert_allclose(loaded.reference_values, model.reference_values)
    assert loaded.profile["transform"] == model.profile["transform"]
    assert loaded.profile["crs"] == model.profile["crs"]


def test_load_rejects_unsupported_schema(small_reference_fixture, tmp_path) -> None:
    model = train_arrays(
        **small_reference_fixture,
        window_radius=1,
        patch_size=5,
    )
    path = tmp_path / "model.npz"
    save_model(model, path)
    with np.load(path, allow_pickle=False) as stored:
        arrays = {name: stored[name] for name in stored.files if name != "metadata_json"}
        metadata = json.loads(str(stored["metadata_json"]))
    metadata["schema_version"] = 99
    np.savez_compressed(path, **arrays, metadata_json=np.array(json.dumps(metadata)))

    with pytest.raises(ValueError, match="Unsupported model schema"):
        load_model(path)
