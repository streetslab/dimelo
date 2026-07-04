# Single-read storage benchmark

This benchmark compares the current DiMeLo HDF5 representation with flattened
HDF5, NetCDF4, Zarr v3, and Parquet. Parquet uses `large_list<uint8>` vectors so
large datasets are not constrained by Arrow's 32-bit list offsets. The benchmark
does not change production dependencies or file formats.

## Environment

Use the project environment for its existing `numpy` and `h5py` installation,
then install the benchmark-only packages:

```bash
conda run -n dimelo-toolkit python -m pip install \
  -r docs/benchmarks/storage_benchmark_requirements.txt
```

The benchmark records exact library versions in its JSON output. OpenZL remains
available through `hdf5_openzl_compression_benchmark.py`, but is not part of the
required decision matrix because it does not currently have a stable dependency.

## Running

Run the complete local profile from the repository root:

```bash
conda run -n dimelo-toolkit python \
  docs/benchmarks/single_read_storage_benchmark.py
```

The defaults use the tracked 1,926-row fixture and deterministic 20K- and
200K-row expansions, one warm-up, five measured repetitions, and 1/4/8-process
region reads. Temporary artifacts are written beneath the ignored
`docs/benchmarks/tmp_bench/` directory. Results are written as JSON, CSV, and a
markdown summary report in `docs/benchmarks/results/`.

For a quick correctness and API smoke run:

```bash
conda run -n dimelo-toolkit python \
  docs/benchmarks/single_read_storage_benchmark.py \
  --sizes fixture --warmups 0 --repeats 1 --concurrency 1
```

Run the same command on a representative HPC POSIX filesystem by passing a work
directory on that filesystem. Keep the results directory in the repository so
the small result files can be compared and committed:

```bash
python docs/benchmarks/single_read_storage_benchmark.py \
  --work-dir /path/on/parallel/filesystem/dimelo-storage-benchmark \
  --results-dir docs/benchmarks/results/hpc
```

## Interpretation

Logical equality is a hard gate. A candidate must improve at least two of total
bytes, write time, region-read time, or random-read time by 20% relative to the
legacy HDF5 baseline without regressing another primary metric by more than 10%.
Flattened HDF5 is preferred if it is within 10% of the best candidate on every
primary metric. Zarr artifact count is reported because directory metadata can
dominate performance on shared filesystems.

`peak_rss_bytes` is sampled every 10 ms from the benchmark process and therefore
includes memory allocated by native HDF5, NetCDF, Arrow, and compression
libraries. Concurrent-read entries report the parent process RSS separately;
they do not sum child-process peaks.
