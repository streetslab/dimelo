#!/usr/bin/env python3
"""Benchmark HDF5 vs OpenZL for dimelo-toolkit intermediate read data.

This script is intended as a local/private benchmark harness. It compares the
existing HDF5 read/write path with an optional OpenZL-backed serialization path.
The OpenZL benchmark is only enabled if an importable `openzl` package is available.
"""

from __future__ import annotations

import argparse
import gzip
import json
import pickle
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from dimelo.load_processed import read_vectors_from_hdf5

try:
    import openzl
except ImportError:  # pragma: no cover
    openzl = None


DEFAULT_MOTIFS = ["A,0", "CG,0"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Benchmark HDF5 and OpenZL compression performance for dimelo intermediate data."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("docs/benchmarks/tmp_bench"),
        help="Directory to write benchmark outputs.",
    )
    parser.add_argument(
        "--benchmark-hdf5",
        action="store_true",
        help="Run the HDF5 benchmark.",
    )
    parser.add_argument(
        "--benchmark-openzl",
        action="store_true",
        help="Run the OpenZL benchmark if available.",
    )
    parser.add_argument(
        "--use-synthetic-data",
        action="store_true",
        help="Use synthetic sample data instead of trying to ingest a real sample file.",
    )
    parser.add_argument(
        "--num-reads",
        type=int,
        default=200,
        help="Number of synthetic reads to generate for the benchmark.",
    )
    parser.add_argument(
        "--save-results",
        action="store_true",
        help="Save benchmark results to benchmark_results.json.",
    )
    return parser.parse_args()


def ensure_directory(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def compress_uint8_vector(data: np.ndarray, compress_level: int = 1) -> bytes:
    return gzip.compress(np.asarray(data, dtype=np.uint8).tobytes(), compresslevel=compress_level)


def generate_synthetic_data(num_reads: int) -> dict[str, list[Any]]:
    data = {
        "read_name": [],
        "chromosome": [],
        "read_start": [],
        "read_end": [],
        "strand": [],
        "motif": [],
        "mod_vector": [],
        "val_vector": [],
    }
    start_position = 0
    for read_index in range(num_reads):
        read_length = int(40 + 20 * (read_index % 4))
        read_name = f"read_{read_index:06d}"
        motif = DEFAULT_MOTIFS[read_index % len(DEFAULT_MOTIFS)]
        strand = "+" if (read_index % 2 == 0) else "-"
        chromosome = "chr1"
        read_start = start_position
        read_end = read_start + read_length
        start_position += read_length // 2

        mod_vector = np.zeros(read_length, dtype=np.uint8)
        val_vector = np.zeros(read_length, dtype=np.uint8)
        positions = np.arange(0, read_length, 5)
        mod_vector[positions] = np.random.randint(0, 255, size=len(positions), dtype=np.uint8)
        val_vector[positions] = 1

        data["read_name"].append(read_name)
        data["chromosome"].append(chromosome)
        data["read_start"].append(read_start)
        data["read_end"].append(read_end)
        data["strand"].append(strand)
        data["motif"].append(motif)
        data["mod_vector"].append(compress_uint8_vector(mod_vector, compress_level=1))
        data["val_vector"].append(compress_uint8_vector(val_vector, compress_level=1))

    return data


def write_hdf5_direct(
    data: dict[str, list[Any]],
    output_h5: Path,
    compress_level: int = 1,
) -> dict[str, Any]:
    start = time.perf_counter()
    with h5py.File(output_h5, "w") as h5:
        dt_str = h5py.string_dtype(encoding="utf-8")
        dt_vlen = h5py.vlen_dtype(np.dtype("uint8"))

        h5.create_dataset("threshold", data=np.array(np.nan), dtype="f8")
        h5.create_dataset(
            "read_name",
            data=np.asarray(data["read_name"], dtype=dt_str),
            dtype=dt_str,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "chromosome",
            data=np.asarray(data["chromosome"], dtype=dt_str),
            dtype=dt_str,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "read_start",
            data=np.asarray(data["read_start"], dtype="i"),
            dtype="i",
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "read_end",
            data=np.asarray(data["read_end"], dtype="i"),
            dtype="i",
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "strand",
            data=np.asarray(data["strand"], dtype=dt_str),
            dtype=dt_str,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "motif",
            data=np.asarray(data["motif"], dtype=dt_str),
            dtype=dt_str,
            compression="gzip",
            compression_opts=9,
        )
        h5.create_dataset(
            "mod_vector",
            shape=(len(data["mod_vector"]),),
            maxshape=(None,),
            dtype=dt_vlen,
            data=None,
        )
        h5["mod_vector"][:] = [np.frombuffer(b, dtype=np.uint8) for b in data["mod_vector"]]
        h5.create_dataset(
            "val_vector",
            shape=(len(data["val_vector"]),),
            maxshape=(None,),
            dtype=dt_vlen,
            data=None,
        )
        h5["val_vector"][:] = [np.frombuffer(b, dtype=np.uint8) for b in data["val_vector"]]
    end = time.perf_counter()
    return {
        "write_time": end - start,
        "file_size": output_h5.stat().st_size,
        "path": str(output_h5),
    }


def read_hdf5_via_dimelo(output_h5: Path, motifs: list[str]) -> dict[str, Any]:
    start = time.perf_counter()
    records, datasets, _ = read_vectors_from_hdf5(
        file=output_h5,
        motifs=motifs,
        regions=None,
        calculate_mod_fractions=False,
        quiet=True,
    )
    elapsed = time.perf_counter() - start
    return {
        "read_time": elapsed,
        "records": len(records),
        "datasets": datasets,
    }


def openzl_available() -> bool:
    return openzl is not None


def compress_with_openzl(payload: bytes) -> bytes:
    if not openzl_available():
        raise RuntimeError("OpenZL is not installed")
    if hasattr(openzl, "compress"):
        return openzl.compress(payload)
    if hasattr(openzl, "compress_bytes"):
        return openzl.compress_bytes(payload)
    raise RuntimeError("Unsupported OpenZL API: no compress() or compress_bytes()")


def decompress_with_openzl(payload: bytes) -> bytes:
    if not openzl_available():
        raise RuntimeError("OpenZL is not installed")
    if hasattr(openzl, "decompress"):
        return openzl.decompress(payload)
    if hasattr(openzl, "decompress_bytes"):
        return openzl.decompress_bytes(payload)
    raise RuntimeError("Unsupported OpenZL API: no decompress() or decompress_bytes()")


def write_openzl_standalone(
    data: dict[str, list[Any]], output_file: Path
) -> dict[str, Any]:
    if not openzl_available():
        return {
            "status": "skipped",
            "error": "OpenZL package not installed",
        }
    payload = pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL)
    start = time.perf_counter()
    compressed = compress_with_openzl(payload)
    output_file.write_bytes(compressed)
    elapsed = time.perf_counter() - start
    return {
        "write_time": elapsed,
        "file_size": output_file.stat().st_size,
        "path": str(output_file),
    }


def read_openzl_standalone(output_file: Path) -> dict[str, Any]:
    if not openzl_available():
        return {
            "status": "skipped",
            "error": "OpenZL package not installed",
        }
    payload = output_file.read_bytes()
    start = time.perf_counter()
    decompressed = decompress_with_openzl(payload)
    data = pickle.loads(decompressed)
    elapsed = time.perf_counter() - start
    return {
        "read_time": elapsed,
        "record_count": len(data.get("read_name", [])),
        "dataset_keys": sorted(data.keys()),
    }


def compare_results(results: list[dict[str, Any]], output_dir: Path) -> None:
    print("\nBenchmark results")
    print("-----------------")
    for result in results:
        print(json.dumps(result, indent=2, sort_keys=True))
    results_path = output_dir / "benchmark_results.json"
    with results_path.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=2)
    print(f"Saved results to: {results_path}")


def main() -> None:
    args = parse_args()
    output_dir = ensure_directory(args.output_dir)
    data = generate_synthetic_data(args.num_reads)

    results: list[dict[str, Any]] = []
    if args.benchmark_hdf5:
        hdf5_path = output_dir / "benchmark_data.h5"
        write_metrics = write_hdf5_direct(data, hdf5_path, compress_level=1)
        read_metrics = read_hdf5_via_dimelo(hdf5_path, motifs=DEFAULT_MOTIFS)
        results.append(
            {
                "format": "HDF5+gzip",
                **write_metrics,
                **read_metrics,
                "status": "success",
            }
        )

    if args.benchmark_openzl:
        openzl_path = output_dir / "benchmark_data.openzl"
        write_metrics = write_openzl_standalone(data, openzl_path)
        read_metrics: dict[str, Any] = {}
        if write_metrics.get("status") != "skipped":
            read_metrics = read_openzl_standalone(openzl_path)
        results.append(
            {
                "format": "OpenZL standalone",
                **write_metrics,
                **read_metrics,
            }
        )

    if not results:
        print("No benchmark selected. Use --benchmark-hdf5, --benchmark-openzl, or both.")
        return

    if args.save_results:
        compare_results(results, output_dir)
    else:
        print("\nResults summary")
        for result in results:
            print(f"- {result['format']}: {result.get('write_time', 'n/a'):.4f}s write, {result.get('read_time', 'n/a'):.4f}s read, {result.get('file_size', 'n/a')} bytes")


if __name__ == "__main__":
    main()
