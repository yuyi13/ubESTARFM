"""
Script: api.py
Objective: Train, predict, serialize, and write Python ubESTARFM products.
Author: Yi Yu
Created: 2026-06-21
Last updated: 2026-06-21
Inputs: Reference and target rasters, arrays, model files, and parameters.
Outputs: Reusable models, RasterData predictions, NPZ models, and GeoTIFFs.
Usage: Imported by the ubestarfm package.
Dependencies: numpy, rasterio, and ubestarfm kernels.
"""

from __future__ import annotations

import json
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import rasterio
from affine import Affine
from rasterio.crs import CRS
from rasterio.transform import from_origin

from . import __version__
from .io import RasterData, read_raster, validate_geometry
from .kernels import predict_patch, train_patch
from .model import UBESTARFMModel


def _patches(rows: int, cols: int, patch_size: int) -> list[tuple[int, int, int, int]]:
    return [
        (
            row_start,
            min(rows, row_start + patch_size),
            col_start,
            min(cols, col_start + patch_size),
        )
        for row_start in range(0, rows, patch_size)
        for col_start in range(0, cols, patch_size)
    ]


def _default_profile(rows: int, cols: int) -> dict[str, Any]:
    return {
        "driver": "GTiff",
        "height": rows,
        "width": cols,
        "count": 1,
        "dtype": "float64",
        "crs": CRS.from_epsg(4326),
        "transform": from_origin(0, rows, 1, 1),
        "nodata": np.nan,
    }


def _sample_sd(values: np.ndarray) -> float:
    finite = values[np.isfinite(values)]
    if finite.size < 2:
        return float("nan")
    return float(np.std(finite, ddof=1))


def _validate_training_parameters(
    window_radius: int,
    patch_size: int,
    method: str,
    workers: int,
) -> tuple[int, int, str, int]:
    window_radius = int(window_radius)
    patch_size = int(patch_size)
    workers = int(workers)
    if window_radius < 1 or patch_size < 1 or workers < 1:
        raise ValueError("window_radius, patch_size, and workers must be positive.")
    if method not in {"zero_bias", "baseline"}:
        raise ValueError("method must be 'zero_bias' or 'baseline'.")
    return window_radius, patch_size, method, workers


def train_arrays(
    fine_1: np.ndarray,
    fine_2: np.ndarray,
    coarse_1: np.ndarray,
    coarse_2: np.ndarray,
    *,
    window_radius: int = 25,
    patch_size: int = 200,
    method: str = "zero_bias",
    workers: int = 1,
    profile: dict[str, Any] | None = None,
) -> UBESTARFMModel:
    """Train a reusable model from aligned two-dimensional arrays."""
    window_radius, patch_size, method, workers = _validate_training_parameters(
        window_radius,
        patch_size,
        method,
        workers,
    )
    references = np.ascontiguousarray(
        np.stack([fine_1, fine_2, coarse_1, coarse_2]),
        dtype=np.float64,
    )
    references[~np.isfinite(references)] = np.nan
    if references.ndim != 3 or references.shape[0] != 4:
        raise ValueError("All reference arrays must be aligned and two-dimensional.")
    rows, cols = references.shape[1:]
    reference_valid = np.all(np.isfinite(references), axis=0)
    if not reference_valid.any():
        raise ValueError("The reference arrays contain no jointly valid pixels.")

    patches = _patches(rows, cols, patch_size)
    thresholds = np.empty((len(patches), 2), dtype=np.float64)
    patch_ids = np.empty((rows, cols), dtype=np.int32)
    for patch_id, (row_start, row_end, col_start, col_end) in enumerate(patches):
        thresholds[patch_id] = (
            _sample_sd(references[0, row_start:row_end, col_start:col_end]),
            _sample_sd(references[1, row_start:row_end, col_start:col_end]),
        )
        patch_ids[row_start:row_end, col_start:col_end] = patch_id + 1

    def execute(patch_id: int):
        row_start, row_end, col_start, col_end = patches[patch_id]
        return train_patch(
            references[0],
            references[1],
            references[2],
            references[3],
            row_start,
            row_end,
            col_start,
            col_end,
            window_radius,
            thresholds[patch_id, 0],
            thresholds[patch_id, 1],
        )

    first_result = execute(0)
    if len(patches) == 1:
        patch_results = [first_result]
    elif workers == 1:
        patch_results = [first_result, *(execute(index) for index in range(1, len(patches)))]
    else:
        with ThreadPoolExecutor(max_workers=min(workers, len(patches))) as executor:
            remaining = list(executor.map(execute, range(1, len(patches))))
        patch_results = [first_result, *remaining]

    side = 2 * window_radius + 1
    mask_bytes = (side * side + 7) // 8
    candidate_masks = np.zeros((rows * cols, mask_bytes), dtype=np.uint8)
    for cell_ids, masks in patch_results:
        candidate_masks[cell_ids] = masks

    return UBESTARFMModel(
        schema_version=1,
        package_version=__version__,
        reference_values=references,
        reference_valid=reference_valid,
        candidate_masks=candidate_masks,
        patch_ids=patch_ids,
        patch_thresholds=thresholds,
        profile=_default_profile(rows, cols) if profile is None else dict(profile),
        window_radius=window_radius,
        patch_size=patch_size,
        method=method,
    )


def train(
    fine_1,
    fine_2,
    coarse_1,
    coarse_2,
    *,
    window_radius: int = 25,
    patch_size: int = 200,
    method: str = "zero_bias",
    cache: str = "candidates",
    workers: int = 1,
) -> UBESTARFMModel:
    """Train a reusable model from raster paths or RasterData inputs."""
    if cache != "candidates":
        raise ValueError("train() supports cache='candidates'.")
    inputs = [read_raster(source) for source in (fine_1, fine_2, coarse_1, coarse_2)]
    validate_geometry(inputs)
    return train_arrays(
        *(raster.values for raster in inputs),
        window_radius=window_radius,
        patch_size=patch_size,
        method=method,
        workers=workers,
        profile=inputs[0].profile,
    )


def unpack_candidates(
    model: UBESTARFMModel,
    *,
    row: int,
    col: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Decode candidate coordinates for one zero-based output cell."""
    rows, cols = model.reference_valid.shape
    if row < 0 or row >= rows or col < 0 or col >= cols:
        raise ValueError("row and col must identify a model cell.")
    side = 2 * model.window_radius + 1
    cell = row * cols + col
    candidate_rows: list[int] = []
    candidate_cols: list[int] = []
    for bit in range(side * side):
        byte = int(model.candidate_masks[cell, bit // 8])
        if byte & (1 << (bit % 8)) == 0:
            continue
        delta_col = bit // side - model.window_radius
        delta_row = bit % side - model.window_radius
        candidate_rows.append(row + delta_row)
        candidate_cols.append(col + delta_col)
    return (
        np.asarray(candidate_rows, dtype=np.int32),
        np.asarray(candidate_cols, dtype=np.int32),
    )


def _validate_value_range(value_range: Sequence[float]) -> tuple[float, float]:
    if len(value_range) != 2:
        raise ValueError("value_range must contain two bounds.")
    value_min, value_max = float(value_range[0]), float(value_range[1])
    if np.isnan(value_min) or np.isnan(value_max) or value_min >= value_max:
        raise ValueError("value_range must contain increasing bounds.")
    return value_min, value_max


def predict_arrays(
    model: UBESTARFMModel,
    coarse_targets: Sequence[np.ndarray],
    value_range: Sequence[float] = (250.0, 350.0),
    *,
    workers: int = 1,
) -> list[np.ndarray]:
    """Predict one or more aligned target arrays."""
    value_min, value_max = _validate_value_range(value_range)
    workers = int(workers)
    if workers < 1:
        raise ValueError("workers must be positive.")
    targets = np.ascontiguousarray(np.stack(coarse_targets), dtype=np.float64)
    targets[~np.isfinite(targets)] = np.nan
    expected_shape = model.reference_valid.shape
    if targets.ndim != 3 or targets.shape[1:] != expected_shape:
        raise ValueError("All target arrays must match the model dimensions.")

    patches = _patches(*expected_shape, model.patch_size)

    def execute(patch):
        row_start, row_end, col_start, col_end = patch
        return predict_patch(
            model.reference_values,
            model.candidate_masks,
            targets,
            row_start,
            row_end,
            col_start,
            col_end,
            model.window_radius,
            model.method == "zero_bias",
            value_min,
            value_max,
        )

    first_result = execute(patches[0])
    if len(patches) == 1:
        patch_results = [first_result]
    elif workers == 1:
        patch_results = [first_result, *(execute(patch) for patch in patches[1:])]
    else:
        with ThreadPoolExecutor(max_workers=min(workers, len(patches))) as executor:
            remaining = list(executor.map(execute, patches[1:]))
        patch_results = [first_result, *remaining]

    predictions = np.full(targets.shape, np.nan, dtype=np.float64)
    for patch, result in zip(patches, patch_results, strict=True):
        row_start, row_end, col_start, col_end = patch
        predictions[:, row_start:row_end, col_start:col_end] = result
    return [predictions[index] for index in range(predictions.shape[0])]


def _model_raster(model: UBESTARFMModel) -> RasterData:
    values = np.empty(model.reference_valid.shape, dtype=np.float64)
    return RasterData(values=values, profile=model.profile)


def _write_prediction(
    prediction: RasterData,
    output_path: str | os.PathLike[str],
    *,
    overwrite: bool,
) -> None:
    path = Path(output_path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    profile = dict(prediction.profile)
    nodata = profile.get("nodata")
    if nodata is None or not np.isfinite(nodata):
        nodata = -3.4e38
    profile.update(
        driver="GTiff",
        height=prediction.values.shape[0],
        width=prediction.values.shape[1],
        count=1,
        dtype="float32",
        nodata=nodata,
    )
    with tempfile.NamedTemporaryFile(
        prefix=f".{path.name}-",
        suffix=".tif",
        dir=path.parent,
        delete=False,
    ) as stream:
        temporary_path = Path(stream.name)
    try:
        with rasterio.open(temporary_path, "w", **profile) as dataset:
            dataset.write(prediction.values.astype(np.float32), 1)
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)


def predict_batch(
    model: UBESTARFMModel,
    coarse_targets: Iterable,
    *,
    output_paths: Sequence[str | os.PathLike[str]] | None = None,
    value_range: Sequence[float] = (250.0, 350.0),
    overwrite: bool = False,
    workers: int = 1,
) -> list[RasterData]:
    """Predict multiple target rasters without retraining."""
    inputs = [read_raster(source) for source in coarse_targets]
    if not inputs:
        raise ValueError("coarse_targets must contain at least one raster.")
    validate_geometry([_model_raster(model), *inputs])
    if output_paths is not None and len(output_paths) != len(inputs):
        raise ValueError("output_paths must contain one path per target.")

    matrices = predict_arrays(
        model,
        [raster.values for raster in inputs],
        value_range,
        workers=workers,
    )
    predictions = [
        RasterData(values=matrix, profile=dict(model.profile)) for matrix in matrices
    ]
    if output_paths is not None:
        for prediction, output_path in zip(predictions, output_paths, strict=True):
            _write_prediction(prediction, output_path, overwrite=overwrite)
    return predictions


def predict(
    model: UBESTARFMModel,
    coarse_target,
    *,
    output_path: str | os.PathLike[str] | None = None,
    value_range: Sequence[float] = (250.0, 350.0),
    overwrite: bool = False,
    workers: int = 1,
) -> RasterData:
    """Predict one target raster."""
    paths = None if output_path is None else [output_path]
    return predict_batch(
        model,
        [coarse_target],
        output_paths=paths,
        value_range=value_range,
        overwrite=overwrite,
        workers=workers,
    )[0]


def fit_predict_batch(
    fine_1,
    fine_2,
    coarse_1,
    coarse_2,
    coarse_targets: Iterable,
    **kwargs,
) -> list[RasterData]:
    """Train a reference pair once and predict a supplied target batch."""
    prediction_keys = {"output_paths", "value_range", "overwrite"}
    prediction_options = {
        key: kwargs.pop(key) for key in list(kwargs) if key in prediction_keys
    }
    model = train(fine_1, fine_2, coarse_1, coarse_2, **kwargs)
    return predict_batch(model, coarse_targets, **prediction_options)


def _profile_to_json(profile: dict[str, Any]) -> dict[str, Any]:
    transform = profile["transform"]
    crs = profile["crs"]
    nodata = profile.get("nodata")
    return {
        "driver": profile.get("driver", "GTiff"),
        "height": int(profile["height"]),
        "width": int(profile["width"]),
        "count": 1,
        "dtype": str(profile.get("dtype", "float32")),
        "transform": list(transform)[:6],
        "crs_wkt": crs.to_wkt(),
        "nodata": None if nodata is None or not np.isfinite(nodata) else float(nodata),
    }


def _profile_from_json(profile: dict[str, Any]) -> dict[str, Any]:
    nodata = profile["nodata"]
    return {
        "driver": profile["driver"],
        "height": int(profile["height"]),
        "width": int(profile["width"]),
        "count": 1,
        "dtype": profile["dtype"],
        "transform": Affine(*profile["transform"]),
        "crs": CRS.from_wkt(profile["crs_wkt"]),
        "nodata": np.nan if nodata is None else float(nodata),
    }


def save_model(
    model: UBESTARFMModel,
    path: str | os.PathLike[str],
) -> Path:
    """Save a compressed, language-specific model archive."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    metadata = {
        "schema_version": model.schema_version,
        "package_version": model.package_version,
        "profile": _profile_to_json(model.profile),
        "window_radius": model.window_radius,
        "patch_size": model.patch_size,
        "method": model.method,
    }
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}-",
        suffix=".npz",
        dir=destination.parent,
        delete=False,
    ) as stream:
        temporary_path = Path(stream.name)
    try:
        np.savez_compressed(
            temporary_path,
            reference_values=model.reference_values,
            reference_valid=model.reference_valid,
            candidate_masks=model.candidate_masks,
            patch_ids=model.patch_ids,
            patch_thresholds=model.patch_thresholds,
            metadata_json=np.array(json.dumps(metadata)),
        )
        os.replace(temporary_path, destination)
    finally:
        temporary_path.unlink(missing_ok=True)
    return destination.resolve()


def load_model(path: str | os.PathLike[str]) -> UBESTARFMModel:
    """Load and validate a compressed Python model archive."""
    with np.load(path, allow_pickle=False) as stored:
        metadata = json.loads(str(stored["metadata_json"].item()))
        return UBESTARFMModel(
            schema_version=int(metadata["schema_version"]),
            package_version=str(metadata["package_version"]),
            reference_values=stored["reference_values"],
            reference_valid=stored["reference_valid"],
            candidate_masks=stored["candidate_masks"],
            patch_ids=stored["patch_ids"],
            patch_thresholds=stored["patch_thresholds"],
            profile=_profile_from_json(metadata["profile"]),
            window_radius=int(metadata["window_radius"]),
            patch_size=int(metadata["patch_size"]),
            method=str(metadata["method"]),
        )
