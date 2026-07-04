# Single-Read Storage Benchmark Report

Generated: `2026-07-03T13:10:47-0700`

## Recommendation

- Recommendation: `retain current HDF5`
- Reason: No alternative satisfied the production follow-up gates.
- Largest measured dataset: `synthetic-200000` (200000 rows)
- Backends: legacy-hdf5, flat-hdf5, netcdf-vlen, netcdf-flat, zarr-flat, parquet

## Primary Metrics On Largest Dataset

| Backend | Status | Size | Write | Region Read | Random Read | Files |
| --- | --- | --- | --- | --- | --- | --- |
| legacy-hdf5 | success | 474.6 MiB | 29.4481 s | 0.5578 s | 0.1407 s | 1 |
| flat-hdf5 | success | 456.9 MiB | 69.4255 s | 0.6042 s | 0.8137 s | 1 |
| netcdf-vlen | success | 8878.7 MiB | 21.7056 s | 0.5581 s | 0.3558 s | 1 |
| netcdf-flat | success | 594.5 MiB | 50.5870 s | 0.5692 s | 0.3946 s | 1 |
| zarr-flat | success | 338.1 MiB | 42.0732 s | 1.1188 s | 1.3084 s | 4451 |
| parquet | success | 459.4 MiB | 76.7794 s | 2.5170 s | 66.5758 s | 101 |

## Improvement Vs Legacy HDF5

| Backend | Size | Write | Region Read | Random Read |
| --- | --- | --- | --- | --- |
| flat-hdf5 | +3.7% | -135.8% | -8.3% | -478.3% |
| legacy-hdf5 | +0.0% | +0.0% | +0.0% | +0.0% |
| netcdf-flat | -25.3% | -71.8% | -2.1% | -180.4% |
| netcdf-vlen | -1770.8% | +26.3% | -0.1% | -152.9% |
| parquet | +3.2% | -160.7% | -351.2% | -47215.6% |
| zarr-flat | +28.8% | -42.9% | -100.6% | -829.9% |

## Benchmark Configuration

- Concurrency settings: `[1]`
- Warmups per benchmark: `1`
- Measured repetitions: `5`
- Chunk rows: `2000`
- Random-row sample size: `256`

## Environment

- Python: `3.11.3 (main, May 21 2023, 16:08:19) [Clang 14.0.3 (clang-1403.0.22.14.1)]`
- Platform: `macOS-26.5.1-x86_64-i386-64bit`
- Libraries: `{"h5py": "3.11.0", "netCDF4": "1.7.4", "numpy": "1.26.4", "pyarrow": "24.0.0", "zarr": "3.0.10"}`
