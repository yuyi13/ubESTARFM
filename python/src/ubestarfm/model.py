"""
Script: model.py
Objective: Define and validate immutable reusable Python ubESTARFM models.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Reference arrays, candidate masks, spatial metadata, and parameters.
Outputs: Validated UBESTARFMModel instances.
Usage: Imported by the ubestarfm package.
Dependencies: numpy.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass(frozen=True, slots=True)
class UBESTARFMModel:
    """A compact, reusable reference-pair model."""

    schema_version: int
    package_version: str
    reference_values: np.ndarray
    reference_valid: np.ndarray
    candidate_masks: np.ndarray
    patch_ids: np.ndarray
    patch_thresholds: np.ndarray
    profile: dict[str, Any]
    window_radius: int
    patch_size: int
    method: str

    def __post_init__(self) -> None:
        if self.schema_version != 1:
            raise ValueError("Unsupported model schema version.")
        if self.method not in {"zero_bias", "baseline"}:
            raise ValueError("method must be 'zero_bias' or 'baseline'.")
        if self.window_radius < 1 or self.patch_size < 1:
            raise ValueError("window_radius and patch_size must be positive.")

        references = np.ascontiguousarray(self.reference_values, dtype=np.float64)
        valid = np.ascontiguousarray(self.reference_valid, dtype=np.bool_)
        masks = np.ascontiguousarray(self.candidate_masks, dtype=np.uint8)
        patch_ids = np.ascontiguousarray(self.patch_ids, dtype=np.int32)
        thresholds = np.ascontiguousarray(self.patch_thresholds, dtype=np.float64)
        profile = dict(self.profile)

        if references.ndim != 3 or references.shape[0] != 4:
            raise ValueError("reference_values must have shape (4, rows, cols).")
        rows, cols = references.shape[1:]
        if valid.shape != (rows, cols):
            raise ValueError("reference_valid shape does not match references.")
        if patch_ids.shape != (rows, cols):
            raise ValueError("patch_ids shape does not match references.")
        if thresholds.ndim != 2 or thresholds.shape[1] != 2:
            raise ValueError("patch_thresholds must have two columns.")

        side = 2 * self.window_radius + 1
        mask_bytes = (side * side + 7) // 8
        if masks.shape != (rows * cols, mask_bytes):
            raise ValueError("candidate_masks shape does not match model parameters.")
        if int(profile.get("height", -1)) != rows:
            raise ValueError("profile height does not match model rows.")
        if int(profile.get("width", -1)) != cols:
            raise ValueError("profile width does not match model columns.")
        if profile.get("transform") is None or profile.get("crs") is None:
            raise ValueError("profile must contain a transform and CRS.")

        for array in (references, valid, masks, patch_ids, thresholds):
            array.setflags(write=False)
        object.__setattr__(self, "reference_values", references)
        object.__setattr__(self, "reference_valid", valid)
        object.__setattr__(self, "candidate_masks", masks)
        object.__setattr__(self, "patch_ids", patch_ids)
        object.__setattr__(self, "patch_thresholds", thresholds)
        object.__setattr__(self, "profile", profile)

    @property
    def mask_bytes(self) -> int:
        """Number of packed candidate bytes stored per output cell."""
        return self.candidate_masks.shape[1]

    @property
    def nbytes(self) -> int:
        """Total bytes used by the model's NumPy arrays."""
        return sum(
            array.nbytes
            for array in (
                self.reference_values,
                self.reference_valid,
                self.candidate_masks,
                self.patch_ids,
                self.patch_thresholds,
            )
        )
