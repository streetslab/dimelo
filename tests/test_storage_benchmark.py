from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

MODULE_PATH = (
    Path(__file__).parents[1]
    / "docs"
    / "benchmarks"
    / "single_read_storage_benchmark.py"
)
SPEC = importlib.util.spec_from_file_location("single_read_storage_benchmark", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
storage = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = storage
SPEC.loader.exec_module(storage)


@pytest.mark.parametrize("thresholded", [False, True])
@pytest.mark.parametrize("backend_name", list(storage.BACKENDS))
def test_storage_backend_round_trip_and_selection(tmp_path, backend_name, thresholded):
    backend = storage.BACKENDS[backend_name]
    available, reason = backend.available()
    if not available:
        pytest.skip(reason)

    data = storage.make_edge_case_data(thresholded=thresholded)
    path = tmp_path / f"data{backend.suffix}"
    backend.write(path, data, chunk_rows=2)

    storage.assert_logical_equal(data, backend.read_all(path))
    storage.assert_logical_equal(data.take([1, 3]), backend.read_indices(path, np.array([1, 3])))

    query = storage.Query("chr1", 9, 30, "A,0")
    expected = data.take(storage._query_indices(data, query))
    storage.assert_logical_equal(expected, backend.read_region(path, query))

    empty_query = storage.Query("chrMissing", 0, 100, None)
    expected_empty = data.take([])
    storage.assert_logical_equal(expected_empty, backend.read_region(path, empty_query))


def test_generate_synthetic_data_is_deterministic():
    template = storage.make_edge_case_data()
    first = storage.generate_synthetic_data(template, 20, seed=17)
    second = storage.generate_synthetic_data(template, 20, seed=17)
    storage.assert_logical_equal(first, second)


def test_recommend_prefers_flat_hdf5_within_ten_percent():
    def result(name, size, write, region, random_read):
        return {
            "backend": name,
            "status": "success",
            "artifact_bytes": size,
            "timings": {
                "write_seconds": {"median": write},
                "region_read_seconds": {"median": region},
                "random_read_seconds": {"median": random_read},
            },
        }

    decision = storage.recommend(
        [
            result("legacy-hdf5", 100, 10, 10, 10),
            result("flat-hdf5", 51, 5.1, 5.1, 5.1),
            result("zarr-flat", 50, 5, 5, 5),
        ]
    )
    assert decision["recommendation"] == "adopt flattened HDF5 v2"
