#!/usr/bin/env python3
"""Benchmark single-read storage layouts used by dimelo-toolkit.

This module is intentionally independent of the production I/O path.  It compares
the current HDF5 representation with flattened HDF5, NetCDF4, Zarr v3, and Parquet
without changing any public dimelo API or dependency.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gc
import gzip
import importlib.metadata
import json
import os
import platform
import shutil
import statistics
import sys
import threading
import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

import h5py
import numpy as np
import psutil

try:  # benchmark-only optional dependency
    import netCDF4
except ImportError:  # pragma: no cover - exercised by availability checks
    netCDF4 = None

try:  # benchmark-only optional dependency
    import pyarrow as pa
    import pyarrow.dataset as pads
    import pyarrow.parquet as pq
except ImportError:  # pragma: no cover - exercised by availability checks
    pa = None
    pads = None
    pq = None

try:  # benchmark-only optional dependency
    import zarr
except ImportError:  # pragma: no cover - exercised by availability checks
    zarr = None


SCHEMA_VERSION = "single-read-benchmark-v1"
DEFAULT_FIXTURE = Path(
    "dimelo/test/data/test_targets/megalodon_peaks_nothresh/"
    "reads.combined_basemods.h5"
)
DEFAULT_BACKENDS = (
    "legacy-hdf5",
    "flat-hdf5",
    "netcdf-vlen",
    "netcdf-flat",
    "zarr-flat",
    "parquet",
)
STRING_COLUMNS = ("read_name", "chromosome", "strand", "motif")
INTEGER_COLUMNS = ("read_start", "read_end")
VECTOR_COLUMNS = ("mod_vector", "val_vector")
PRIMARY_METRICS = (
    "artifact_bytes",
    "write_seconds",
    "region_read_seconds",
    "random_read_seconds",
)


@dataclass(frozen=True)
class Query:
    chromosome: str
    start: int
    end: int
    motif: str | None = None


@dataclass
class ReadData:
    read_name: np.ndarray
    chromosome: np.ndarray
    read_start: np.ndarray
    read_end: np.ndarray
    strand: np.ndarray
    motif: np.ndarray
    mod_vector: list[np.ndarray]
    val_vector: list[np.ndarray]
    threshold_by_motif: dict[str, float | None]

    def __len__(self) -> int:
        return len(self.read_name)

    def take(self, indices: Iterable[int]) -> ReadData:
        idx = np.asarray(list(indices), dtype=np.int64)
        return ReadData(
            read_name=self.read_name[idx],
            chromosome=self.chromosome[idx],
            read_start=self.read_start[idx],
            read_end=self.read_end[idx],
            strand=self.strand[idx],
            motif=self.motif[idx],
            mod_vector=[self.mod_vector[int(i)] for i in idx],
            val_vector=[self.val_vector[int(i)] for i in idx],
            threshold_by_motif=dict(self.threshold_by_motif),
        )

    def validate(self) -> None:
        lengths = {
            len(self.read_name),
            len(self.chromosome),
            len(self.read_start),
            len(self.read_end),
            len(self.strand),
            len(self.motif),
            len(self.mod_vector),
            len(self.val_vector),
        }
        if len(lengths) != 1:
            raise ValueError(f"Canonical columns have inconsistent lengths: {lengths}")
        for index, (mod, val) in enumerate(
            zip(self.mod_vector, self.val_vector, strict=True)
        ):
            if mod.dtype != np.uint8 or val.dtype != np.uint8:
                raise ValueError(f"Row {index} vectors must use uint8")
            if len(mod) != len(val):
                raise ValueError(f"Row {index} vector lengths differ")


class Backend(Protocol):
    name: str
    suffix: str

    def available(self) -> tuple[bool, str | None]: ...

    def write(self, path: Path, data: ReadData, chunk_rows: int) -> None: ...

    def read_all(self, path: Path) -> ReadData: ...

    def read_region(self, path: Path, query: Query) -> ReadData: ...

    def read_indices(self, path: Path, indices: np.ndarray) -> ReadData: ...

    def open_metadata(self, path: Path) -> dict[str, Any]: ...

    def capabilities(self) -> dict[str, Any]: ...


def _decode_strings(values: Any) -> np.ndarray:
    return np.asarray(
        [
            value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
            for value in values
        ],
        dtype=object,
    )


def _json_threshold_map(data: ReadData) -> str:
    return json.dumps(data.threshold_by_motif, sort_keys=True, separators=(",", ":"))


def _threshold_map(raw: str | bytes | None, motifs: np.ndarray) -> dict[str, float | None]:
    if raw is not None:
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        return {str(key): value for key, value in json.loads(raw).items()}
    return {str(motif): None for motif in np.unique(motifs)}


def _query_indices(data: ReadData, query: Query) -> np.ndarray:
    mask = (
        (data.chromosome == query.chromosome)
        & (data.read_end > query.start)
        & (data.read_start < query.end)
    )
    if query.motif is not None:
        mask &= data.motif == query.motif
    return np.flatnonzero(mask)


def _flatten_vectors(vectors: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    offsets = np.zeros(len(vectors) + 1, dtype=np.int64)
    if vectors:
        offsets[1:] = np.cumsum([len(vector) for vector in vectors], dtype=np.int64)
        values = np.concatenate(vectors).astype(np.uint8, copy=False)
    else:
        values = np.empty(0, dtype=np.uint8)
    return offsets, values


def _inflate_vectors(
    offsets: np.ndarray, values: np.ndarray, *, copy: bool = True
) -> list[np.ndarray]:
    return [
        (
            np.asarray(values[int(start) : int(end)], dtype=np.uint8).copy()
            if copy
            else np.asarray(values[int(start) : int(end)], dtype=np.uint8)
        )
        for start, end in zip(offsets[:-1], offsets[1:], strict=True)
    ]


def _select_flat_vectors(
    offsets: np.ndarray, values: Any, indices: np.ndarray
) -> list[np.ndarray]:
    if len(indices) == len(offsets) - 1 and np.array_equal(
        indices, np.arange(len(indices), dtype=np.int64)
    ):
        return _inflate_vectors(offsets, np.asarray(values[:], dtype=np.uint8))
    return [
        np.asarray(values[int(offsets[index]) : int(offsets[index + 1])], dtype=np.uint8)
        for index in indices
    ]


def _hdf5_take(dataset: Any, indices: np.ndarray) -> Any:
    indices = np.asarray(indices, dtype=np.int64)
    if len(indices) == 0:
        return dataset[:0]
    if len(indices) == 1 or np.all(indices[1:] >= indices[:-1]):
        return dataset[list(indices)]
    return [dataset[int(index)] for index in indices]


def _slice_bounds(total: int, chunk_rows: int) -> Iterable[tuple[int, int]]:
    for start in range(0, total, chunk_rows):
        yield start, min(start + chunk_rows, total)


def load_legacy_fixture(path: Path) -> ReadData:
    """Decode a production HDF5 fixture into the canonical logical schema."""
    with h5py.File(path, "r") as handle:
        motifs = _decode_strings(handle["motif"][:])
        threshold_raw = None
        if "threshold_by_motif_json" in handle:
            threshold_raw = handle["threshold_by_motif_json"][()]
        elif "threshold" in handle:
            threshold = float(handle["threshold"][()])
            threshold_raw = json.dumps(
                {
                    str(motif): None if np.isnan(threshold) else threshold
                    for motif in np.unique(motifs)
                }
            )

        def decode_vector(value: Any) -> np.ndarray:
            payload = np.asarray(value, dtype=np.uint8).tobytes()
            return np.frombuffer(gzip.decompress(payload), dtype=np.uint8).copy()

        data = ReadData(
            read_name=_decode_strings(handle["read_name"][:]),
            chromosome=_decode_strings(handle["chromosome"][:]),
            read_start=np.asarray(handle["read_start"][:], dtype=np.int64),
            read_end=np.asarray(handle["read_end"][:], dtype=np.int64),
            strand=_decode_strings(handle["strand"][:]),
            motif=motifs,
            mod_vector=[decode_vector(value) for value in handle["mod_vector"][:]],
            val_vector=[decode_vector(value) for value in handle["val_vector"][:]],
            threshold_by_motif=_threshold_map(threshold_raw, motifs),
        )
    data.validate()
    return data


def generate_synthetic_data(template: ReadData, rows: int, seed: int) -> ReadData:
    """Sample a fixture deterministically while preserving ragged-vector distributions."""
    if rows <= 0:
        raise ValueError("rows must be positive")
    rng = np.random.default_rng(seed)
    sampled = rng.integers(0, len(template), size=rows, dtype=np.int64)
    vectors_mod = [template.mod_vector[int(index)] for index in sampled]
    vectors_val = [template.val_vector[int(index)] for index in sampled]
    lengths = np.asarray([len(vector) for vector in vectors_mod], dtype=np.int64)
    chromosomes = np.asarray(
        [f"chr{1 + (index // 50_000) % 3}" for index in range(rows)], dtype=object
    )
    starts = np.arange(rows, dtype=np.int64) * 75
    # Reset genomic positions per chromosome-sized block so every chromosome has useful queries.
    starts %= 50_000 * 75
    motifs = template.motif[sampled]
    data = ReadData(
        read_name=np.asarray([f"synthetic_{index:09d}" for index in range(rows)], dtype=object),
        chromosome=chromosomes,
        read_start=starts,
        read_end=starts + np.maximum(lengths, 1),
        strand=template.strand[sampled],
        motif=motifs,
        mod_vector=vectors_mod,
        val_vector=vectors_val,
        threshold_by_motif={
            str(motif): template.threshold_by_motif.get(str(motif))
            for motif in np.unique(motifs)
        },
    )
    data.validate()
    return data


def make_edge_case_data(thresholded: bool = False) -> ReadData:
    """Small canonical dataset used by tests and smoke runs."""
    threshold = 0.5 if thresholded else None
    data = ReadData(
        read_name=np.asarray(["empty", "short", "long", "motif-2"], dtype=object),
        chromosome=np.asarray(["chr1", "chr1", "chr2", "chr1"], dtype=object),
        read_start=np.asarray([0, 10, 100, 20], dtype=np.int64),
        read_end=np.asarray([1, 13, 4196, 24], dtype=np.int64),
        strand=np.asarray(["+", "-", "+", "+"], dtype=object),
        motif=np.asarray(["A,0", "A,0", "A,0", "CG,0"], dtype=object),
        mod_vector=[
            np.empty(0, dtype=np.uint8),
            np.asarray([0, 127, 255], dtype=np.uint8),
            np.arange(4096, dtype=np.uint16).astype(np.uint8),
            np.asarray([1, 0, 1, 0], dtype=np.uint8),
        ],
        val_vector=[
            np.empty(0, dtype=np.uint8),
            np.asarray([1, 1, 1], dtype=np.uint8),
            np.ones(4096, dtype=np.uint8),
            np.asarray([1, 0, 1, 0], dtype=np.uint8),
        ],
        threshold_by_motif={"A,0": threshold, "CG,0": threshold},
    )
    data.validate()
    return data


def assert_logical_equal(expected: ReadData, actual: ReadData) -> None:
    expected.validate()
    actual.validate()
    for column in (*STRING_COLUMNS, *INTEGER_COLUMNS):
        left = np.asarray(getattr(expected, column))
        right = np.asarray(getattr(actual, column))
        if not np.array_equal(left, right):
            raise AssertionError(f"Column {column!r} did not round-trip")
    if expected.threshold_by_motif != actual.threshold_by_motif:
        raise AssertionError("Threshold metadata did not round-trip")
    for column in VECTOR_COLUMNS:
        left_vectors = getattr(expected, column)
        right_vectors = getattr(actual, column)
        if len(left_vectors) != len(right_vectors):
            raise AssertionError(f"Vector column {column!r} row count differs")
        for index, (left, right) in enumerate(
            zip(left_vectors, right_vectors, strict=True)
        ):
            if not np.array_equal(left, right):
                raise AssertionError(f"{column}[{index}] did not round-trip")


def _prepare_path(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()
    path.parent.mkdir(parents=True, exist_ok=True)


class LegacyHDF5Backend:
    name = "legacy-hdf5"
    suffix = ".h5"

    def available(self) -> tuple[bool, str | None]:
        return True, None

    def capabilities(self) -> dict[str, Any]:
        return {"reopen_append": True, "layout": "gzip bytes in HDF5 VLEN uint8"}

    def write(self, path: Path, data: ReadData, chunk_rows: int) -> None:
        _prepare_path(path)
        string_type = h5py.string_dtype("utf-8")
        vlen_type = h5py.vlen_dtype(np.dtype("uint8"))
        with h5py.File(path, "w") as handle:
            handle.create_dataset("threshold", data=_legacy_threshold(data))
            handle.create_dataset("threshold_by_motif_json", data=_json_threshold_map(data))
            for column in STRING_COLUMNS:
                handle.create_dataset(
                    column,
                    shape=(0,),
                    maxshape=(None,),
                    dtype=string_type,
                    compression="gzip",
                    compression_opts=9,
                )
            for column in INTEGER_COLUMNS:
                handle.create_dataset(
                    column,
                    shape=(0,),
                    maxshape=(None,),
                    dtype="i8",
                    compression="gzip",
                    compression_opts=9,
                )
            for column in VECTOR_COLUMNS:
                handle.create_dataset(column, shape=(0,), maxshape=(None,), dtype=vlen_type)

            for start, end in _slice_bounds(len(data), chunk_rows):
                for column in (*STRING_COLUMNS, *INTEGER_COLUMNS, *VECTOR_COLUMNS):
                    dataset = handle[column]
                    dataset.resize((end,))
                    if column in VECTOR_COLUMNS:
                        values = [
                            np.frombuffer(
                                gzip.compress(vector.tobytes(), compresslevel=1),
                                dtype=np.uint8,
                            )
                            for vector in getattr(data, column)[start:end]
                        ]
                        dataset[start:end] = values
                    else:
                        dataset[start:end] = getattr(data, column)[start:end]

    def _read(self, path: Path, indices: np.ndarray | None) -> ReadData:
        with h5py.File(path, "r") as handle:
            if indices is None:
                indices = np.arange(len(handle["read_name"]), dtype=np.int64)
            selection = np.asarray(indices, dtype=np.int64)
            motifs = _decode_strings(_hdf5_take(handle["motif"], selection))

            def vectors(column: str) -> list[np.ndarray]:
                return [
                    np.frombuffer(
                        gzip.decompress(np.asarray(value, dtype=np.uint8).tobytes()),
                        dtype=np.uint8,
                    ).copy()
                    for value in _hdf5_take(handle[column], selection)
                ]

            return ReadData(
                read_name=_decode_strings(_hdf5_take(handle["read_name"], selection)),
                chromosome=_decode_strings(_hdf5_take(handle["chromosome"], selection)),
                read_start=np.asarray(_hdf5_take(handle["read_start"], selection), dtype=np.int64),
                read_end=np.asarray(_hdf5_take(handle["read_end"], selection), dtype=np.int64),
                strand=_decode_strings(_hdf5_take(handle["strand"], selection)),
                motif=motifs,
                mod_vector=vectors("mod_vector"),
                val_vector=vectors("val_vector"),
                threshold_by_motif=_threshold_map(
                    handle["threshold_by_motif_json"][()], motifs
                ),
            )

    def read_all(self, path: Path) -> ReadData:
        return self._read(path, None)

    def read_indices(self, path: Path, indices: np.ndarray) -> ReadData:
        return self._read(path, indices)

    def read_region(self, path: Path, query: Query) -> ReadData:
        with h5py.File(path, "r") as handle:
            metadata = ReadData(
                read_name=_decode_strings(handle["read_name"][:]),
                chromosome=_decode_strings(handle["chromosome"][:]),
                read_start=np.asarray(handle["read_start"][:], dtype=np.int64),
                read_end=np.asarray(handle["read_end"][:], dtype=np.int64),
                strand=_decode_strings(handle["strand"][:]),
                motif=_decode_strings(handle["motif"][:]),
                mod_vector=[],
                val_vector=[],
                threshold_by_motif={},
            )
        return self._read(path, _query_indices_metadata(metadata, query))

    def open_metadata(self, path: Path) -> dict[str, Any]:
        with h5py.File(path, "r") as handle:
            return {"rows": len(handle["read_name"]), "keys": sorted(handle.keys())}


def _legacy_threshold(data: ReadData) -> float:
    values = list(data.threshold_by_motif.values())
    if not values or any(value is None for value in values):
        return float("nan")
    unique = {float(value) for value in values if value is not None}
    return unique.pop() if len(unique) == 1 else 0.0


def _query_indices_metadata(data: ReadData, query: Query) -> np.ndarray:
    mask = (
        (data.chromosome == query.chromosome)
        & (data.read_end > query.start)
        & (data.read_start < query.end)
    )
    if query.motif is not None:
        mask &= data.motif == query.motif
    return np.flatnonzero(mask)


class FlatHDF5Backend:
    name = "flat-hdf5"
    suffix = ".flat.h5"

    def available(self) -> tuple[bool, str | None]:
        return True, None

    def capabilities(self) -> dict[str, Any]:
        return {"reopen_append": True, "layout": "contiguous values plus offsets"}

    def write(self, path: Path, data: ReadData, chunk_rows: int) -> None:
        _prepare_path(path)
        string_type = h5py.string_dtype("utf-8")
        with h5py.File(path, "w") as handle:
            handle.attrs["schema_version"] = SCHEMA_VERSION
            handle.attrs["threshold_by_motif_json"] = _json_threshold_map(data)
            for column in STRING_COLUMNS:
                handle.create_dataset(
                    column,
                    shape=(0,),
                    maxshape=(None,),
                    dtype=string_type,
                    chunks=(max(1, min(chunk_rows, 4096)),),
                    compression="gzip",
                    compression_opts=1,
                )
            for column in INTEGER_COLUMNS:
                handle.create_dataset(
                    column,
                    shape=(0,),
                    maxshape=(None,),
                    dtype="i8",
                    chunks=(max(1, min(chunk_rows, 4096)),),
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                )
            for prefix in ("mod", "val"):
                handle.create_dataset(
                    f"{prefix}_values",
                    shape=(0,),
                    maxshape=(None,),
                    dtype="u1",
                    chunks=(262_144,),
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                )
                handle.create_dataset(
                    f"{prefix}_offsets",
                    data=np.asarray([0], dtype=np.int64),
                    maxshape=(None,),
                    chunks=(max(2, min(chunk_rows + 1, 4096)),),
                    compression="gzip",
                    compression_opts=1,
                    shuffle=True,
                )

            for start, end in _slice_bounds(len(data), chunk_rows):
                for column in (*STRING_COLUMNS, *INTEGER_COLUMNS):
                    dataset = handle[column]
                    dataset.resize((end,))
                    dataset[start:end] = getattr(data, column)[start:end]
                for prefix, column in (("mod", "mod_vector"), ("val", "val_vector")):
                    offsets, values = _flatten_vectors(getattr(data, column)[start:end])
                    value_dataset = handle[f"{prefix}_values"]
                    old_value_count = len(value_dataset)
                    value_dataset.resize((old_value_count + len(values),))
                    value_dataset[old_value_count:] = values
                    offset_dataset = handle[f"{prefix}_offsets"]
                    old_offset_count = len(offset_dataset)
                    offset_dataset.resize((old_offset_count + end - start,))
                    offset_dataset[old_offset_count:] = offsets[1:] + old_value_count

    def _read(self, path: Path, indices: np.ndarray | None) -> ReadData:
        with h5py.File(path, "r") as handle:
            if indices is None:
                indices = np.arange(len(handle["read_name"]), dtype=np.int64)
            indices = np.asarray(indices, dtype=np.int64)
            motifs = _decode_strings(_hdf5_take(handle["motif"], indices))
            return ReadData(
                read_name=_decode_strings(_hdf5_take(handle["read_name"], indices)),
                chromosome=_decode_strings(_hdf5_take(handle["chromosome"], indices)),
                read_start=np.asarray(_hdf5_take(handle["read_start"], indices), dtype=np.int64),
                read_end=np.asarray(_hdf5_take(handle["read_end"], indices), dtype=np.int64),
                strand=_decode_strings(_hdf5_take(handle["strand"], indices)),
                motif=motifs,
                mod_vector=_select_flat_vectors(
                    np.asarray(handle["mod_offsets"][:], dtype=np.int64),
                    handle["mod_values"],
                    indices,
                ),
                val_vector=_select_flat_vectors(
                    np.asarray(handle["val_offsets"][:], dtype=np.int64),
                    handle["val_values"],
                    indices,
                ),
                threshold_by_motif=_threshold_map(
                    handle.attrs.get("threshold_by_motif_json"), motifs
                ),
            )

    def read_all(self, path: Path) -> ReadData:
        return self._read(path, None)

    def read_indices(self, path: Path, indices: np.ndarray) -> ReadData:
        return self._read(path, indices)

    def read_region(self, path: Path, query: Query) -> ReadData:
        with h5py.File(path, "r") as handle:
            chromosome = _decode_strings(handle["chromosome"][:])
            starts = np.asarray(handle["read_start"][:], dtype=np.int64)
            ends = np.asarray(handle["read_end"][:], dtype=np.int64)
            motifs = _decode_strings(handle["motif"][:])
        mask = (chromosome == query.chromosome) & (ends > query.start) & (starts < query.end)
        if query.motif is not None:
            mask &= motifs == query.motif
        return self._read(path, np.flatnonzero(mask))

    def open_metadata(self, path: Path) -> dict[str, Any]:
        with h5py.File(path, "r") as handle:
            return {"rows": len(handle["read_name"]), "keys": sorted(handle.keys())}


class NetCDFBackend:
    suffix = ".nc"

    def __init__(self, *, flat: bool):
        self.flat = flat
        self.name = "netcdf-flat" if flat else "netcdf-vlen"

    def available(self) -> tuple[bool, str | None]:
        if netCDF4 is None:
            return False, "netCDF4 is not installed"
        return True, None

    def capabilities(self) -> dict[str, Any]:
        return {
            "reopen_append": True,
            "layout": "contiguous values plus offsets" if self.flat else "native NetCDF VLEN",
        }

    def write(self, path: Path, data: ReadData, chunk_rows: int) -> None:
        assert netCDF4 is not None
        _prepare_path(path)
        with netCDF4.Dataset(path, "w", format="NETCDF4") as handle:
            handle.setncattr("schema_version", SCHEMA_VERSION)
            handle.setncattr("threshold_by_motif_json", _json_threshold_map(data))
            handle.createDimension("read", None)
            for column in STRING_COLUMNS:
                handle.createVariable(column, str, ("read",))
            for column in INTEGER_COLUMNS:
                handle.createVariable(
                    column, "i8", ("read",), compression="zlib", complevel=1, shuffle=True
                )
            if self.flat:
                handle.createDimension("offset", None)
                handle.createDimension("mod_value", None)
                handle.createDimension("val_value", None)
                for prefix in ("mod", "val"):
                    handle.createVariable(
                        f"{prefix}_offsets",
                        "i8",
                        ("offset",),
                        compression="zlib",
                        complevel=1,
                        shuffle=True,
                    )
                    handle.createVariable(
                        f"{prefix}_values",
                        "u1",
                        (f"{prefix}_value",),
                        compression="zlib",
                        complevel=1,
                        shuffle=True,
                    )
                    handle.variables[f"{prefix}_offsets"][0] = 0
            else:
                vlen_type = handle.createVLType(np.uint8, "uint8_vlen")
                for column in VECTOR_COLUMNS:
                    handle.createVariable(column, vlen_type, ("read",))

            value_counts = {"mod": 0, "val": 0}
            for start, end in _slice_bounds(len(data), chunk_rows):
                for column in (*STRING_COLUMNS, *INTEGER_COLUMNS):
                    handle.variables[column][start:end] = getattr(data, column)[start:end]
                if self.flat:
                    for prefix, column in (("mod", "mod_vector"), ("val", "val_vector")):
                        offsets, values = _flatten_vectors(getattr(data, column)[start:end])
                        value_start = value_counts[prefix]
                        handle.variables[f"{prefix}_values"][
                            value_start : value_start + len(values)
                        ] = values
                        handle.variables[f"{prefix}_offsets"][start + 1 : end + 1] = (
                            offsets[1:] + value_start
                        )
                        value_counts[prefix] += len(values)
                else:
                    for column in VECTOR_COLUMNS:
                        values = np.empty(end - start, dtype=object)
                        values[:] = getattr(data, column)[start:end]
                        handle.variables[column][start:end] = values

    def _read(self, path: Path, indices: np.ndarray | None) -> ReadData:
        assert netCDF4 is not None
        with netCDF4.Dataset(path, "r") as handle:
            if indices is None:
                indices = np.arange(len(handle.dimensions["read"]), dtype=np.int64)
            indices = np.asarray(indices, dtype=np.int64)
            motifs = _decode_strings(handle.variables["motif"][indices])
            if self.flat:
                mod = _select_flat_vectors(
                    np.asarray(handle.variables["mod_offsets"][:], dtype=np.int64),
                    handle.variables["mod_values"],
                    indices,
                )
                val = _select_flat_vectors(
                    np.asarray(handle.variables["val_offsets"][:], dtype=np.int64),
                    handle.variables["val_values"],
                    indices,
                )
            else:
                mod = [np.asarray(value, dtype=np.uint8) for value in handle.variables["mod_vector"][indices]]
                val = [np.asarray(value, dtype=np.uint8) for value in handle.variables["val_vector"][indices]]
            return ReadData(
                read_name=_decode_strings(handle.variables["read_name"][indices]),
                chromosome=_decode_strings(handle.variables["chromosome"][indices]),
                read_start=np.asarray(handle.variables["read_start"][indices], dtype=np.int64),
                read_end=np.asarray(handle.variables["read_end"][indices], dtype=np.int64),
                strand=_decode_strings(handle.variables["strand"][indices]),
                motif=motifs,
                mod_vector=mod,
                val_vector=val,
                threshold_by_motif=_threshold_map(
                    handle.getncattr("threshold_by_motif_json"), motifs
                ),
            )

    def read_all(self, path: Path) -> ReadData:
        return self._read(path, None)

    def read_indices(self, path: Path, indices: np.ndarray) -> ReadData:
        return self._read(path, indices)

    def read_region(self, path: Path, query: Query) -> ReadData:
        assert netCDF4 is not None
        with netCDF4.Dataset(path, "r") as handle:
            chromosome = _decode_strings(handle.variables["chromosome"][:])
            starts = np.asarray(handle.variables["read_start"][:], dtype=np.int64)
            ends = np.asarray(handle.variables["read_end"][:], dtype=np.int64)
            motifs = _decode_strings(handle.variables["motif"][:])
        mask = (chromosome == query.chromosome) & (ends > query.start) & (starts < query.end)
        if query.motif is not None:
            mask &= motifs == query.motif
        return self._read(path, np.flatnonzero(mask))

    def open_metadata(self, path: Path) -> dict[str, Any]:
        assert netCDF4 is not None
        with netCDF4.Dataset(path, "r") as handle:
            return {
                "rows": len(handle.dimensions["read"]),
                "keys": sorted(handle.variables),
            }


class ZarrBackend:
    name = "zarr-flat"
    suffix = ".zarr"

    def available(self) -> tuple[bool, str | None]:
        if zarr is None:
            return False, "zarr is not installed"
        major = int(str(zarr.__version__).split(".", maxsplit=1)[0])
        if major < 3:
            return False, f"Zarr >=3 is required, found {zarr.__version__}"
        return True, None

    def capabilities(self) -> dict[str, Any]:
        return {
            "reopen_append": True,
            "layout": "Zarr v3 sharded contiguous values plus offsets",
        }

    def _create_array(
        self,
        group: Any,
        name: str,
        *,
        shape: tuple[int, ...],
        dtype: Any,
        chunks: tuple[int, ...],
        shards: tuple[int, ...],
    ) -> Any:
        assert zarr is not None
        compressor = zarr.codecs.ZstdCodec(level=3)
        return group.create_array(
            name,
            shape=shape,
            dtype=dtype,
            chunks=chunks,
            shards=shards,
            compressors=[compressor],
        )

    def write(self, path: Path, data: ReadData, chunk_rows: int) -> None:
        assert zarr is not None
        _prepare_path(path)
        group = zarr.create_group(store=str(path), overwrite=True, zarr_format=3)
        group.attrs["schema_version"] = SCHEMA_VERSION
        group.attrs["threshold_by_motif_json"] = _json_threshold_map(data)
        row_chunk = max(1, min(chunk_rows, 4096))
        row_shard = max(row_chunk, min(max(row_chunk * 8, row_chunk), max(len(data), row_chunk)))
        encoded_strings: dict[str, np.ndarray] = {}
        for column in STRING_COLUMNS:
            dictionary = sorted({str(value) for value in getattr(data, column)})
            codes = {value: index for index, value in enumerate(dictionary)}
            encoded_strings[column] = np.asarray(
                [codes[str(value)] for value in getattr(data, column)], dtype=np.int32
            )
            group.attrs[f"{column}_dictionary"] = dictionary
            self._create_array(
                group,
                column,
                shape=(0,),
                dtype="i4",
                chunks=(row_chunk,),
                shards=(row_shard,),
            )
        for column in INTEGER_COLUMNS:
            self._create_array(
                group,
                column,
                shape=(0,),
                dtype="i8",
                chunks=(row_chunk,),
                shards=(row_shard,),
            )
        value_chunk = 262_144
        value_shard = value_chunk * 8
        for prefix in ("mod", "val"):
            self._create_array(
                group,
                f"{prefix}_values",
                shape=(0,),
                dtype="u1",
                chunks=(value_chunk,),
                shards=(value_shard,),
            )
            offsets = self._create_array(
                group,
                f"{prefix}_offsets",
                shape=(1,),
                dtype="i8",
                chunks=(row_chunk,),
                shards=(row_shard,),
            )
            offsets[0] = 0

        value_counts = {"mod": 0, "val": 0}
        for start, end in _slice_bounds(len(data), chunk_rows):
            for column in (*STRING_COLUMNS, *INTEGER_COLUMNS):
                array = group[column]
                array.resize((end,))
                values = encoded_strings[column][start:end] if column in STRING_COLUMNS else getattr(data, column)[start:end]
                array[start:end] = values
            for prefix, column in (("mod", "mod_vector"), ("val", "val_vector")):
                offsets, values = _flatten_vectors(getattr(data, column)[start:end])
                value_array = group[f"{prefix}_values"]
                value_start = value_counts[prefix]
                value_array.resize((value_start + len(values),))
                value_array[value_start:] = values
                offset_array = group[f"{prefix}_offsets"]
                offset_array.resize((end + 1,))
                offset_array[start + 1 : end + 1] = offsets[1:] + value_start
                value_counts[prefix] += len(values)

    def _decode_column(self, group: Any, column: str, selection: Any) -> np.ndarray:
        dictionary = np.asarray(group.attrs[f"{column}_dictionary"], dtype=object)
        codes = np.asarray(selection, dtype=np.int64)
        return dictionary[codes]

    def _read(self, path: Path, indices: np.ndarray | None) -> ReadData:
        assert zarr is not None
        group = zarr.open_group(store=str(path), mode="r")
        if indices is None:
            indices = np.arange(group["read_name"].shape[0], dtype=np.int64)
        indices = np.asarray(indices, dtype=np.int64)
        motifs = self._decode_column(group, "motif", group["motif"].oindex[indices])
        return ReadData(
            read_name=self._decode_column(
                group, "read_name", group["read_name"].oindex[indices]
            ),
            chromosome=self._decode_column(
                group, "chromosome", group["chromosome"].oindex[indices]
            ),
            read_start=np.asarray(group["read_start"].oindex[indices], dtype=np.int64),
            read_end=np.asarray(group["read_end"].oindex[indices], dtype=np.int64),
            strand=self._decode_column(
                group, "strand", group["strand"].oindex[indices]
            ),
            motif=motifs,
            mod_vector=_select_flat_vectors(
                np.asarray(group["mod_offsets"][:], dtype=np.int64),
                group["mod_values"],
                indices,
            ),
            val_vector=_select_flat_vectors(
                np.asarray(group["val_offsets"][:], dtype=np.int64),
                group["val_values"],
                indices,
            ),
            threshold_by_motif=_threshold_map(
                group.attrs.get("threshold_by_motif_json"), motifs
            ),
        )

    def read_all(self, path: Path) -> ReadData:
        return self._read(path, None)

    def read_indices(self, path: Path, indices: np.ndarray) -> ReadData:
        return self._read(path, indices)

    def read_region(self, path: Path, query: Query) -> ReadData:
        assert zarr is not None
        group = zarr.open_group(store=str(path), mode="r")
        chromosome = self._decode_column(group, "chromosome", group["chromosome"][:])
        starts = np.asarray(group["read_start"][:], dtype=np.int64)
        ends = np.asarray(group["read_end"][:], dtype=np.int64)
        motifs = self._decode_column(group, "motif", group["motif"][:])
        mask = (chromosome == query.chromosome) & (ends > query.start) & (starts < query.end)
        if query.motif is not None:
            mask &= motifs == query.motif
        return self._read(path, np.flatnonzero(mask))

    def open_metadata(self, path: Path) -> dict[str, Any]:
        assert zarr is not None
        group = zarr.open_group(store=str(path), mode="r")
        return {
            "rows": group["read_name"].shape[0],
            "keys": sorted(group.array_keys()),
        }


class ParquetBackend:
    name = "parquet"
    suffix = ".parquet.d"

    def available(self) -> tuple[bool, str | None]:
        if pa is None or pads is None or pq is None:
            return False, "pyarrow is not installed"
        return True, None

    def capabilities(self) -> dict[str, Any]:
        return {
            "reopen_append": True,
            "layout": "partition fragments with list<uint8> vectors",
        }

    def _table(self, data: ReadData, row_offset: int = 0) -> Any:
        assert pa is not None
        schema = pa.schema(
            [
                pa.field("_row_index", pa.int64()),
                pa.field("read_name", pa.string()),
                pa.field("chromosome", pa.string()),
                pa.field("read_start", pa.int64()),
                pa.field("read_end", pa.int64()),
                pa.field("strand", pa.string()),
                pa.field("motif", pa.string()),
                pa.field("mod_vector", pa.large_list(pa.uint8())),
                pa.field("val_vector", pa.large_list(pa.uint8())),
            ],
            metadata={
                b"schema_version": SCHEMA_VERSION.encode(),
                b"threshold_by_motif_json": _json_threshold_map(data).encode(),
            },
        )
        arrays = [pa.array(np.arange(row_offset, row_offset + len(data)), type=pa.int64())]
        for column in STRING_COLUMNS:
            arrays.append(pa.array(getattr(data, column).tolist(), type=pa.string()))
        for column in INTEGER_COLUMNS:
            arrays.append(pa.array(getattr(data, column), type=pa.int64()))
        for column in VECTOR_COLUMNS:
            arrays.append(
                pa.array(getattr(data, column), type=pa.large_list(pa.uint8()))
            )
        # Arrays above follow canonical grouping, not schema order; reorder explicitly.
        by_name = {
            "_row_index": arrays[0],
            "read_name": arrays[1],
            "chromosome": arrays[2],
            "strand": arrays[3],
            "motif": arrays[4],
            "read_start": arrays[5],
            "read_end": arrays[6],
            "mod_vector": arrays[7],
            "val_vector": arrays[8],
        }
        return pa.Table.from_arrays([by_name[field.name] for field in schema], schema=schema)

    def write(self, path: Path, data: ReadData, chunk_rows: int) -> None:
        assert pq is not None
        _prepare_path(path)
        path.mkdir(parents=True)
        (path / "_dimelo_metadata.json").write_text(
            json.dumps(
                {
                    "schema_version": SCHEMA_VERSION,
                    "threshold_by_motif": data.threshold_by_motif,
                },
                indent=2,
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        for fragment, (start, end) in enumerate(_slice_bounds(len(data), chunk_rows)):
            table = self._table(data.take(range(start, end)), row_offset=start)
            pq.write_table(
                table,
                path / f"part-{fragment:06d}.parquet",
                compression="zstd",
                use_dictionary=list(STRING_COLUMNS),
                row_group_size=chunk_rows,
                write_statistics=True,
            )

    def _dataset(self, path: Path) -> Any:
        assert pads is not None
        return pads.dataset(path, format="parquet", exclude_invalid_files=True)

    def _from_table(
        self, path: Path, table: Any, requested_order: np.ndarray | None = None
    ) -> ReadData:
        metadata = json.loads((path / "_dimelo_metadata.json").read_text(encoding="utf-8"))
        if table.num_rows:
            row_indices = np.asarray(table["_row_index"].to_numpy(), dtype=np.int64)
            if requested_order is not None:
                position_by_row = {int(row): offset for offset, row in enumerate(row_indices)}
                order = [position_by_row[int(row)] for row in requested_order]
                table = table.take(pa.array(order, type=pa.int64()))
            elif np.any(row_indices[1:] < row_indices[:-1]):
                order = np.argsort(row_indices)
                table = table.take(pa.array(order, type=pa.int64()))
        motifs = np.asarray(table["motif"].to_pylist(), dtype=object)

        def vectors(column: str) -> list[np.ndarray]:
            output = []
            for array in table[column].chunks:
                offsets = np.asarray(array.offsets.to_numpy(), dtype=np.int64)
                values = np.asarray(array.values.to_numpy(), dtype=np.uint8)
                output.extend(_inflate_vectors(offsets, values, copy=False))
            return output

        return ReadData(
            read_name=np.asarray(table["read_name"].to_pylist(), dtype=object),
            chromosome=np.asarray(table["chromosome"].to_pylist(), dtype=object),
            read_start=np.asarray(table["read_start"].to_numpy(), dtype=np.int64),
            read_end=np.asarray(table["read_end"].to_numpy(), dtype=np.int64),
            strand=np.asarray(table["strand"].to_pylist(), dtype=object),
            motif=motifs,
            mod_vector=vectors("mod_vector"),
            val_vector=vectors("val_vector"),
            threshold_by_motif={
                str(key): value for key, value in metadata["threshold_by_motif"].items()
            },
        )

    def read_all(self, path: Path) -> ReadData:
        return self._from_table(path, self._dataset(path).to_table())

    def read_indices(self, path: Path, indices: np.ndarray) -> ReadData:
        assert pads is not None
        requested_order = np.asarray(indices, dtype=np.int64)
        expression = pads.field("_row_index").isin(requested_order.tolist())
        return self._from_table(
            path,
            self._dataset(path).to_table(filter=expression),
            requested_order=requested_order,
        )

    def read_region(self, path: Path, query: Query) -> ReadData:
        assert pads is not None
        expression = (
            (pads.field("chromosome") == query.chromosome)
            & (pads.field("read_end") > query.start)
            & (pads.field("read_start") < query.end)
        )
        if query.motif is not None:
            expression &= pads.field("motif") == query.motif
        return self._from_table(path, self._dataset(path).to_table(filter=expression))

    def open_metadata(self, path: Path) -> dict[str, Any]:
        dataset = self._dataset(path)
        return {
            "rows": dataset.count_rows(),
            "keys": dataset.schema.names,
            "fragments": len(list(dataset.get_fragments())),
        }


BACKENDS: dict[str, Backend] = {
    "legacy-hdf5": LegacyHDF5Backend(),
    "flat-hdf5": FlatHDF5Backend(),
    "netcdf-vlen": NetCDFBackend(flat=False),
    "netcdf-flat": NetCDFBackend(flat=True),
    "zarr-flat": ZarrBackend(),
    "parquet": ParquetBackend(),
}


def artifact_stats(path: Path) -> dict[str, int]:
    if path.is_file():
        return {"artifact_bytes": path.stat().st_size, "artifact_files": 1}
    files = [candidate for candidate in path.rglob("*") if candidate.is_file()]
    return {
        "artifact_bytes": sum(candidate.stat().st_size for candidate in files),
        "artifact_files": len(files),
    }


def _measure(operation: Callable[[], Any]) -> tuple[Any, float, int]:
    process = psutil.Process()
    peak_rss = process.memory_info().rss
    stop_sampling = threading.Event()

    def sample_rss() -> None:
        nonlocal peak_rss
        while not stop_sampling.wait(0.01):
            peak_rss = max(peak_rss, process.memory_info().rss)

    sampler = threading.Thread(target=sample_rss, daemon=True)
    sampler.start()
    start = time.perf_counter()
    try:
        result = operation()
        elapsed = time.perf_counter() - start
    finally:
        stop_sampling.set()
        sampler.join()
        peak_rss = max(peak_rss, process.memory_info().rss)
    return result, elapsed, int(peak_rss)


def _summary(values: list[float]) -> dict[str, float]:
    return {
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def _release_native_memory() -> None:
    gc.collect()
    if pa is not None:
        pa.default_memory_pool().release_unused()


def _concurrent_region_worker(
    backend_name: str, path: str, query: tuple[str, int, int, str | None]
) -> int:
    backend = BACKENDS[backend_name]
    selected = backend.read_region(Path(path), Query(*query))
    return len(selected)


def _concurrent_region_read(
    backend: Backend, path: Path, query: Query, workers: int
) -> list[int]:
    if workers == 1:
        return [len(backend.read_region(path, query))]
    packed_query = (query.chromosome, query.start, query.end, query.motif)
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        return list(
            executor.map(
                _concurrent_region_worker,
                [backend.name] * workers,
                [str(path)] * workers,
                [packed_query] * workers,
            )
        )


def benchmark_backend(
    backend: Backend,
    data: ReadData,
    output_path: Path,
    *,
    chunk_rows: int,
    warmups: int,
    repeats: int,
    random_rows: int,
    concurrency: tuple[int, ...],
    seed: int,
) -> dict[str, Any]:
    available, reason = backend.available()
    if not available:
        return {"backend": backend.name, "status": "skipped", "reason": reason}

    rng = np.random.default_rng(seed)
    random_indices = np.sort(
        rng.choice(len(data), size=min(random_rows, len(data)), replace=False)
    )
    midpoint = int(np.median(data.read_start))
    query = Query(
        chromosome=str(data.chromosome[len(data) // 2]),
        start=max(0, midpoint - 5_000),
        end=midpoint + 5_000,
        motif=str(data.motif[len(data) // 2]),
    )
    expected_region = data.take(_query_indices(data, query))
    expected_random = data.take(random_indices)

    timings: dict[str, list[float]] = {
        "write_seconds": [],
        "open_seconds": [],
        "full_read_seconds": [],
        "region_read_seconds": [],
        "random_read_seconds": [],
    }
    peaks: dict[str, list[int]] = {key: [] for key in timings}
    total_runs = warmups + repeats
    for run in range(total_runs):
        _, write_time, write_peak = _measure(
            lambda: backend.write(output_path, data, chunk_rows)
        )
        _release_native_memory()
        _, open_time, open_peak = _measure(lambda: backend.open_metadata(output_path))
        _release_native_memory()
        full, full_time, full_peak = _measure(lambda: backend.read_all(output_path))
        assert_logical_equal(data, full)
        del full
        _release_native_memory()
        region, region_time, region_peak = _measure(
            lambda: backend.read_region(output_path, query)
        )
        assert_logical_equal(expected_region, region)
        del region
        _release_native_memory()
        random_selection, random_time, random_peak = _measure(
            lambda: backend.read_indices(output_path, random_indices)
        )
        assert_logical_equal(expected_random, random_selection)
        del random_selection
        _release_native_memory()
        if run >= warmups:
            for key, elapsed, peak in (
                ("write_seconds", write_time, write_peak),
                ("open_seconds", open_time, open_peak),
                ("full_read_seconds", full_time, full_peak),
                ("region_read_seconds", region_time, region_peak),
                ("random_read_seconds", random_time, random_peak),
            ):
                timings[key].append(elapsed)
                peaks[key].append(peak)

    concurrency_results: dict[str, dict[str, Any]] = {}
    for workers in concurrency:
        counts, elapsed, peak = _measure(
            lambda workers=workers: _concurrent_region_read(
                backend, output_path, query, workers
            )
        )
        if counts != [len(expected_region)] * workers:
            raise AssertionError(f"Concurrent reads returned incorrect counts for {backend.name}")
        concurrency_results[str(workers)] = {
            "seconds": elapsed,
            "parent_peak_rss_bytes": peak,
            "total_rows": sum(counts),
        }

    result: dict[str, Any] = {
        "backend": backend.name,
        "status": "success",
        "rows": len(data),
        "query_rows": len(expected_region),
        "random_rows": len(random_indices),
        "capabilities": backend.capabilities(),
        "timings": {key: _summary(values) for key, values in timings.items()},
        "peak_rss_bytes": {
            key: max(values) for key, values in peaks.items()
        },
        "concurrent_region_reads": concurrency_results,
        **artifact_stats(output_path),
    }
    return result


def environment_metadata() -> dict[str, Any]:
    packages = {}
    for package in ("h5py", "netCDF4", "numpy", "pyarrow", "zarr"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "packages": packages,
    }


def recommend(results: list[dict[str, Any]]) -> dict[str, Any]:
    successful = [result for result in results if result.get("status") == "success"]
    if not successful:
        return {"recommendation": "no decision", "reason": "No backend completed."}

    def metric(result: dict[str, Any], name: str) -> float:
        if name == "artifact_bytes":
            return float(result[name])
        return float(result["timings"][name]["median"])

    baseline = next(
        (result for result in successful if result["backend"] == "legacy-hdf5"),
        successful[0],
    )
    improvements: dict[str, dict[str, float]] = {}
    eligible: list[dict[str, Any]] = []
    for result in successful:
        changes = {
            name: (metric(baseline, name) - metric(result, name)) / metric(baseline, name)
            for name in PRIMARY_METRICS
            if metric(baseline, name) > 0
        }
        improvements[result["backend"]] = changes
        improved = sum(value >= 0.20 for value in changes.values())
        regressed = any(value < -0.10 for value in changes.values())
        if improved >= 2 and not regressed:
            eligible.append(result)

    best_by_metric = {
        name: min(successful, key=lambda result: metric(result, name))["backend"]
        for name in PRIMARY_METRICS
    }
    flat_hdf5 = next(
        (result for result in successful if result["backend"] == "flat-hdf5"), None
    )
    if flat_hdf5 is not None:
        within_ten_percent = all(
            metric(flat_hdf5, name)
            <= min(metric(result, name) for result in successful) * 1.10
            for name in PRIMARY_METRICS
        )
        if within_ten_percent:
            return {
                "recommendation": "adopt flattened HDF5 v2",
                "reason": "Flattened HDF5 is within 10% of the best backend on every primary metric.",
                "best_by_metric": best_by_metric,
                "improvements_vs_legacy": improvements,
            }

    if eligible:
        winner = min(
            eligible,
            key=lambda result: sum(
                metric(result, name) / metric(baseline, name) for name in PRIMARY_METRICS
            ),
        )
        mapping = {
            "flat-hdf5": "adopt flattened HDF5 v2",
            "zarr-flat": "prototype Zarr",
            "parquet": "prototype Parquet",
            "netcdf-flat": "prototype NetCDF4",
            "netcdf-vlen": "prototype NetCDF4",
            "legacy-hdf5": "retain current HDF5",
        }
        return {
            "recommendation": mapping[winner["backend"]],
            "reason": "Winner satisfies the 20% improvement and 10% regression gates.",
            "best_by_metric": best_by_metric,
            "improvements_vs_legacy": improvements,
        }
    return {
        "recommendation": "retain current HDF5",
        "reason": "No alternative satisfied the production follow-up gates.",
        "best_by_metric": best_by_metric,
        "improvements_vs_legacy": improvements,
    }


def _flatten_result_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for dataset in payload["datasets"]:
        for result in dataset["results"]:
            row = {
                "dataset": dataset["name"],
                "rows": dataset["rows"],
                "backend": result["backend"],
                "status": result["status"],
                "artifact_bytes": result.get("artifact_bytes"),
                "artifact_files": result.get("artifact_files"),
                "reason": result.get("reason"),
            }
            for metric in (
                "write_seconds",
                "open_seconds",
                "full_read_seconds",
                "region_read_seconds",
                "random_read_seconds",
            ):
                row[metric] = result.get("timings", {}).get(metric, {}).get("median")
            rows.append(row)
    return rows


def _primary_metric_value(result: dict[str, Any], metric: str) -> float | None:
    if result.get("status") != "success":
        return None
    if metric == "artifact_bytes":
        value = result.get(metric)
        return float(value) if value is not None else None
    timing = result.get("timings", {}).get(metric, {}).get("median")
    return float(timing) if timing is not None else None


def _format_metric(metric: str, value: float | None) -> str:
    if value is None:
        return "n/a"
    if metric == "artifact_bytes":
        return f"{value / (1024**2):.1f} MiB"
    return f"{value:.4f} s"


def _format_percent(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value * 100:+.1f}%"


def _write_report(payload: dict[str, Any], output_dir: Path) -> None:
    datasets = payload.get("datasets", [])
    if not datasets:
        return
    largest = max(datasets, key=lambda dataset: int(dataset["rows"]))
    decision = payload.get("decision", {})
    improvements = decision.get("improvements_vs_legacy", {})
    environment = payload.get("environment", {})
    configuration = payload.get("configuration", {})

    lines = [
        "# Single-Read Storage Benchmark Report",
        "",
        f"Generated: `{payload.get('created_at', 'unknown')}`",
        "",
        "## Recommendation",
        "",
        f"- Recommendation: `{decision.get('recommendation', 'no decision')}`",
        f"- Reason: {decision.get('reason', 'No reason recorded.')}",
        f"- Largest measured dataset: `{largest['name']}` ({largest['rows']} rows)",
        f"- Backends: {', '.join(result['backend'] for result in largest['results'])}",
        "",
        "## Primary Metrics On Largest Dataset",
        "",
        "| Backend | Status | Size | Write | Region Read | Random Read | Files |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    for result in largest["results"]:
        lines.append(
            "| {backend} | {status} | {size} | {write} | {region} | {random} | {files} |".format(
                backend=result["backend"],
                status=result["status"],
                size=_format_metric("artifact_bytes", _primary_metric_value(result, "artifact_bytes")),
                write=_format_metric("write_seconds", _primary_metric_value(result, "write_seconds")),
                region=_format_metric(
                    "region_read_seconds", _primary_metric_value(result, "region_read_seconds")
                ),
                random=_format_metric(
                    "random_read_seconds", _primary_metric_value(result, "random_read_seconds")
                ),
                files=result.get("artifact_files", "n/a"),
            )
        )

    lines.extend(
        [
            "",
            "## Improvement Vs Legacy HDF5",
            "",
            "| Backend | Size | Write | Region Read | Random Read |",
            "| --- | --- | --- | --- | --- |",
        ]
    )
    for backend in sorted(improvements):
        change = improvements[backend]
        lines.append(
            "| {backend} | {size} | {write} | {region} | {random} |".format(
                backend=backend,
                size=_format_percent(change.get("artifact_bytes")),
                write=_format_percent(change.get("write_seconds")),
                region=_format_percent(change.get("region_read_seconds")),
                random=_format_percent(change.get("random_read_seconds")),
            )
        )

    lines.extend(
        [
            "",
            "## Benchmark Configuration",
            "",
            f"- Concurrency settings: `{configuration.get('concurrency')}`",
            f"- Warmups per benchmark: `{configuration.get('warmups')}`",
            f"- Measured repetitions: `{configuration.get('repeats')}`",
            f"- Chunk rows: `{configuration.get('chunk_rows')}`",
            f"- Random-row sample size: `{configuration.get('random_rows')}`",
            "",
            "## Environment",
            "",
            f"- Python: `{environment.get('python', 'unknown')}`",
            f"- Platform: `{environment.get('platform', 'unknown')}`",
            f"- Libraries: `{json.dumps(environment.get('packages', {}), sort_keys=True)}`",
            "",
        ]
    )
    (output_dir / "storage_benchmark_report.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def write_results(payload: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "storage_benchmark_results.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )
    rows = _flatten_result_rows(payload)
    with (output_dir / "storage_benchmark_results.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["status"])
        writer.writeheader()
        writer.writerows(rows)
    _write_report(payload, output_dir)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE)
    parser.add_argument(
        "--sizes",
        default="fixture,20000,200000",
        help="Comma-separated dataset sizes; use 'fixture' for the unexpanded fixture.",
    )
    parser.add_argument(
        "--backends",
        default=",".join(DEFAULT_BACKENDS),
        help=f"Comma-separated backends: {', '.join(DEFAULT_BACKENDS)}",
    )
    parser.add_argument("--work-dir", type=Path, default=Path("docs/benchmarks/tmp_bench"))
    parser.add_argument("--results-dir", type=Path, default=Path("docs/benchmarks/results"))
    parser.add_argument("--chunk-rows", type=int, default=2_000)
    parser.add_argument("--random-rows", type=int, default=256)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--concurrency", default="1,4,8")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> dict[str, Any]:
    template = load_legacy_fixture(args.fixture)
    backend_names = [name.strip() for name in args.backends.split(",") if name.strip()]
    unknown = sorted(set(backend_names) - set(BACKENDS))
    if unknown:
        raise ValueError(f"Unknown backends: {', '.join(unknown)}")
    concurrency = tuple(int(value) for value in args.concurrency.split(",") if value)
    size_specs = [value.strip() for value in args.sizes.split(",") if value.strip()]
    datasets: list[dict[str, Any]] = []
    all_results: list[dict[str, Any]] = []
    for size_spec in size_specs:
        if size_spec == "fixture":
            data = template
            dataset_name = "fixture"
        else:
            rows = int(size_spec)
            data = generate_synthetic_data(template, rows, args.seed)
            dataset_name = f"synthetic-{rows}"
        dataset_results = []
        for backend_name in backend_names:
            backend = BACKENDS[backend_name]
            output_path = args.work_dir / dataset_name / f"data{backend.suffix}"
            print(f"[{dataset_name}] {backend_name}", flush=True)
            _release_native_memory()
            try:
                result = benchmark_backend(
                    backend,
                    data,
                    output_path,
                    chunk_rows=args.chunk_rows,
                    warmups=args.warmups,
                    repeats=args.repeats,
                    random_rows=args.random_rows,
                    concurrency=concurrency,
                    seed=args.seed,
                )
            except Exception as error:
                result = {
                    "backend": backend_name,
                    "status": "failed",
                    "reason": f"{type(error).__name__}: {error}",
                }
            dataset_results.append(result)
            all_results.append(result)
            print(f"  {result['status']}: {result.get('reason', '')}", flush=True)
        datasets.append(
            {"name": dataset_name, "rows": len(data), "results": dataset_results}
        )

    largest_dataset = datasets[-1]["results"] if datasets else []
    payload = {
        "schema_version": SCHEMA_VERSION,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "environment": environment_metadata(),
        "configuration": {
            "fixture": str(args.fixture),
            "sizes": size_specs,
            "backends": backend_names,
            "chunk_rows": args.chunk_rows,
            "random_rows": args.random_rows,
            "warmups": args.warmups,
            "repeats": args.repeats,
            "concurrency": concurrency,
            "seed": args.seed,
        },
        "datasets": datasets,
        "decision": recommend(largest_dataset),
    }
    write_results(payload, args.results_dir)
    return payload


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    payload = run(args)
    print(json.dumps(payload["decision"], indent=2, sort_keys=True))
    failed = [
        result
        for dataset in payload["datasets"]
        for result in dataset["results"]
        if result["status"] == "failed"
    ]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
