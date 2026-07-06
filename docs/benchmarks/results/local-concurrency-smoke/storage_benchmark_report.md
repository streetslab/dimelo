# Single-Read Storage Benchmark Report

Generated: `2026-07-03T21:15:38-0700`

## Recommendation

- Recommendation: `retain current HDF5`
- Reason: No alternative satisfied the production follow-up gates.
- Largest measured dataset: `fixture` (1926 rows)
- Backends: legacy-hdf5, flat-hdf5, netcdf-vlen, netcdf-flat, zarr-flat, parquet

## Primary Metrics On Largest Dataset

| Backend | Status | Size | Write | Region Read | Random Read | Files |
| --- | --- | --- | --- | --- | --- | --- |
| legacy-hdf5 | success | 4.6 MiB | 0.4280 s | 0.0082 s | 0.0222 s | 1 |
| flat-hdf5 | success | 4.4 MiB | 0.7034 s | 0.0054 s | 0.2862 s | 1 |
| netcdf-vlen | success | 85.7 MiB | 0.1309 s | 0.0112 s | 0.1529 s | 1 |
| netcdf-flat | success | 5.8 MiB | 0.5482 s | 0.0144 s | 0.2544 s | 1 |
| zarr-flat | success | 2.5 MiB | 0.4221 s | 0.0947 s | 1.0899 s | 61 |
| parquet | success | 4.1 MiB | 1.1248 s | 0.4832 s | 0.5439 s | 2 |

## Improvement Vs Legacy HDF5

| Backend | Size | Write | Region Read | Random Read |
| --- | --- | --- | --- | --- |
| flat-hdf5 | +5.8% | -64.4% | +34.8% | -1191.1% |
| legacy-hdf5 | +0.0% | +0.0% | +0.0% | +0.0% |
| netcdf-flat | -24.8% | -28.1% | -75.1% | -1047.5% |
| netcdf-vlen | -1746.1% | +69.4% | -36.5% | -589.9% |
| parquet | +10.6% | -162.8% | -5775.7% | -2353.8% |
| zarr-flat | +47.1% | +1.4% | -1051.0% | -4816.8% |

## Benchmark Configuration

- Concurrency settings: `(1, 4, 8)`
- Warmups per benchmark: `0`
- Measured repetitions: `1`
- Chunk rows: `2000`
- Random-row sample size: `256`

## Environment

- Python: `3.11.3 (main, May 21 2023, 16:08:19) [Clang 14.0.3 (clang-1403.0.22.14.1)]`
- Platform: `macOS-26.5.1-x86_64-i386-64bit`
- Libraries: `{"h5py": "3.11.0", "netCDF4": "1.7.4", "numpy": "1.26.4", "pyarrow": "24.0.0", "zarr": "3.0.10"}`
