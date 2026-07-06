import h5py
import numpy as np
import pytest

from dimelo import load_processed
from dimelo import utils as dimelo_utils


def _single_region_chunks():
    return [
        {
            "chromosome": "chr1",
            "region_start": 100,
            "region_end": 104,
            "subregion_start": 100,
            "subregion_end": 104,
            "strand": "+",
        }
    ]


class _FakeTabix:
    def __init__(self, *_args, **_kwargs):
        self.contigs = {"chr1"}

    def fetch(self, _chromosome, _start, _end):
        return ["row_a", "row_b"]


def test_process_pileup_row_parses_modkit_18_column_format():
    parsed_motif = dimelo_utils.ParsedMotif("A,0")
    row = "chr1\t9167144\t9167145\tY\t12\t+\t9167144\t9167145\t255,0,0\t12\t0.00\t3\t9\t0\t0\t0\t0\t0"
    keep, coord, modified, valid = load_processed.process_pileup_row(
        row=row,
        parsed_motif=parsed_motif,
        region_strand="+",
        single_strand=False,
    )
    assert keep is True
    assert coord == 9167144
    assert modified == 3
    assert valid == 12


def test_process_pileup_row_skips_malformed_rows_without_raising():
    parsed_motif = dimelo_utils.ParsedMotif("A,0")
    malformed_row = "chr1\t1\t2\tY\tbad\t+\t1\t2\t255,0,0\tbad\t0.00\tbad"
    with pytest.warns(RuntimeWarning, match="Skipping malformed bedMethyl row"):
        keep, coord, modified, valid = load_processed.process_pileup_row(
            row=malformed_row,
            parsed_motif=parsed_motif,
            region_strand="+",
            single_strand=False,
        )
    assert (keep, coord, modified, valid) == (False, 0, 0, 0)


def test_resolve_cores_for_task_count_auto_small_batch_is_single_core(monkeypatch):
    monkeypatch.setattr(load_processed.utils, "cores_to_run", lambda _cores: 8)
    resolved = load_processed._resolve_cores_for_task_count(
        requested_cores=None,
        task_count=20,
    )
    assert resolved == 1


def test_resolve_cores_for_task_count_explicit_cores_are_preserved(monkeypatch):
    monkeypatch.setattr(
        load_processed.utils, "cores_to_run", lambda _cores: int(_cores)
    )
    resolved = load_processed._resolve_cores_for_task_count(
        requested_cores=4,
        task_count=2,
    )
    assert resolved == 4


def test_memory_limited_cores_caps_auto_parallelism(monkeypatch):
    monkeypatch.setattr(load_processed.utils, "cores_to_run", lambda _cores: 8)
    monkeypatch.setattr(
        load_processed, "_available_memory_bytes", lambda: 1024 * 1024 * 1024
    )
    monkeypatch.setattr(load_processed, "AUTO_PARALLEL_MEMORY_FRACTION", 0.5)
    monkeypatch.setattr(
        load_processed, "AUTO_PARALLEL_PROCESS_OVERHEAD_BYTES", 64 * 1024 * 1024
    )

    resolved = load_processed._memory_limited_cores(
        requested_cores=None,
        task_count=1024,
        bytes_per_task=128 * 1024 * 1024,
        extra_shared_bytes=0,
    )
    assert resolved == 2


def test_memory_limited_cores_keeps_explicit_override(monkeypatch):
    monkeypatch.setattr(
        load_processed.utils, "cores_to_run", lambda _cores: int(_cores)
    )
    monkeypatch.setattr(
        load_processed, "_available_memory_bytes", lambda: 128 * 1024 * 1024
    )

    resolved = load_processed._memory_limited_cores(
        requested_cores=6,
        task_count=1,
        bytes_per_task=1024 * 1024 * 1024,
        extra_shared_bytes=0,
    )
    assert resolved == 6


def test_memory_limited_cores_falls_back_when_memory_unknown(monkeypatch):
    monkeypatch.setattr(load_processed.utils, "cores_to_run", lambda _cores: 8)
    monkeypatch.setattr(load_processed, "_available_memory_bytes", lambda: None)

    resolved = load_processed._memory_limited_cores(
        requested_cores=None,
        task_count=160,
        bytes_per_task=256 * 1024 * 1024,
        extra_shared_bytes=0,
    )
    assert resolved == 8


def test_slurm_memory_limit_bytes_uses_mem_per_cpu(monkeypatch):
    monkeypatch.setenv("SLURM_MEM_PER_CPU", "2048")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    assert load_processed._slurm_memory_limit_bytes() == 2048 * 4 * 1024 * 1024


def test_slurm_memory_limit_bytes_uses_mem_per_node(monkeypatch):
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "8192")
    monkeypatch.delenv("SLURM_MEM_PER_CPU", raising=False)
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    assert load_processed._slurm_memory_limit_bytes() == 8192 * 1024 * 1024


def test_available_memory_bytes_respects_override(monkeypatch):
    monkeypatch.setenv("DIMELO_AVAILABLE_MEMORY_BYTES", "12345")
    monkeypatch.setattr(load_processed, "_cgroup_memory_limit_bytes", lambda: 999)
    monkeypatch.setattr(load_processed, "_slurm_memory_limit_bytes", lambda: 999)
    assert load_processed._available_memory_bytes() == 12345


def test_cores_to_run_respects_slurm_and_affinity(monkeypatch):
    monkeypatch.setattr(dimelo_utils.multiprocessing, "cpu_count", lambda: 64)
    monkeypatch.setattr(
        dimelo_utils.os,
        "sched_getaffinity",
        lambda _pid: set(range(8)),
        raising=False,
    )
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")
    assert dimelo_utils.cores_to_run(None) == 4
    assert dimelo_utils.cores_to_run(12) == 4


def test_get_tabix_file_caches_by_path(monkeypatch):
    class _FakeTabixFile:
        calls = 0

        def __init__(self, _path):
            _FakeTabixFile.calls += 1
            self.contigs = set()

    load_processed._clear_tabix_cache()
    monkeypatch.setattr(load_processed.pysam, "TabixFile", _FakeTabixFile)

    first = load_processed._get_tabix_file("same.bed.gz")
    second = load_processed._get_tabix_file("same.bed.gz")
    third = load_processed._get_tabix_file("other.bed.gz")

    assert first is second
    assert first is not third
    assert _FakeTabixFile.calls == 2


def test_pileup_vectors_from_bedmethyl_single_core_bypasses_process_pool(monkeypatch):
    def _raise_if_called(*_args, **_kwargs):
        raise AssertionError(
            "ProcessPoolExecutor should not be created for single-core execution."
        )

    def fake_process_pileup_row(row, **_kwargs):
        if row == "row_a":
            return True, 100, 2, 5
        if row == "row_b":
            return True, 102, 1, 3
        raise AssertionError(f"Unexpected row payload: {row}")

    monkeypatch.setattr(
        load_processed.utils,
        "regions_dict_from_input",
        lambda *_args, **_kwargs: {"chr1": [(100, 104, "+")]},
    )
    monkeypatch.setattr(
        load_processed.utils,
        "process_chunks_from_regions_dict",
        lambda *_args, **_kwargs: _single_region_chunks(),
    )
    monkeypatch.setattr(load_processed.utils, "cores_to_run", lambda _cores: 1)
    monkeypatch.setattr(
        load_processed.concurrent.futures, "ProcessPoolExecutor", _raise_if_called
    )
    monkeypatch.setattr(load_processed.pysam, "TabixFile", _FakeTabix)
    monkeypatch.setattr(load_processed, "process_pileup_row", fake_process_pileup_row)

    modified, valid = load_processed.pileup_vectors_from_bedmethyl(
        bedmethyl_file="fake.bed.gz",
        motif="A,0",
        regions="chr1:100-104,+",
        quiet=True,
    )

    np.testing.assert_array_equal(modified, np.array([2, 0, 1, 0], dtype=np.int32))
    np.testing.assert_array_equal(valid, np.array([5, 0, 3, 0], dtype=np.int32))


def test_regions_to_list_single_core_bypasses_process_pool(monkeypatch):
    def _raise_if_called(*_args, **_kwargs):
        raise AssertionError(
            "ProcessPoolExecutor should not be created for single-core execution."
        )

    monkeypatch.setattr(
        load_processed.utils,
        "regions_dict_from_input",
        lambda *_args, **_kwargs: {"chr1": [(100, 110, "+"), (120, 130, "-")]},
    )
    monkeypatch.setattr(load_processed.utils, "cores_to_run", lambda _cores: 1)
    monkeypatch.setattr(
        load_processed.concurrent.futures, "ProcessPoolExecutor", _raise_if_called
    )

    result = load_processed.regions_to_list(
        function_handle=lambda regions, **_kwargs: f"ok:{regions}",
        regions="ignored",
        quiet=True,
    )

    assert result == ["ok:chr1:100-110,+", "ok:chr1:120-130,-"]


def test_regions_to_list_auto_cores_prefers_single_core_for_small_batches(monkeypatch):
    def _raise_if_called(*_args, **_kwargs):
        raise AssertionError(
            "ProcessPoolExecutor should not be created for small auto-core batches."
        )

    monkeypatch.setattr(
        load_processed.utils,
        "regions_dict_from_input",
        lambda *_args, **_kwargs: {"chr1": [(100, 110, "+"), (120, 130, "-")]},
    )
    monkeypatch.setattr(load_processed.utils, "cores_to_run", lambda _cores: 8)
    monkeypatch.setattr(
        load_processed.concurrent.futures, "ProcessPoolExecutor", _raise_if_called
    )

    result = load_processed.regions_to_list(
        function_handle=lambda regions, **_kwargs: f"ok:{regions}",
        regions="ignored",
        quiet=True,
        cores=None,
    )

    assert result == ["ok:chr1:100-110,+", "ok:chr1:120-130,-"]


def test_regions_to_list_explicit_cores_still_uses_process_pool(monkeypatch):
    class _FakeExecutor:
        last_max_workers = None

        def __init__(self, *, max_workers):
            _FakeExecutor.last_max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def map(self, fn, iterable):
            return map(fn, iterable)

    monkeypatch.setattr(
        load_processed.utils,
        "regions_dict_from_input",
        lambda *_args, **_kwargs: {"chr1": [(100, 110, "+"), (120, 130, "-")]},
    )
    monkeypatch.setattr(
        load_processed.utils, "cores_to_run", lambda _cores: int(_cores)
    )
    monkeypatch.setattr(
        load_processed.concurrent.futures, "ProcessPoolExecutor", _FakeExecutor
    )

    result = load_processed.regions_to_list(
        function_handle=lambda regions, **_kwargs: f"ok:{regions}",
        regions="ignored",
        quiet=True,
        cores=2,
    )

    assert result == ["ok:chr1:100-110,+", "ok:chr1:120-130,-"]
    assert _FakeExecutor.last_max_workers == 2


def test_regions_to_list_parallel_batches_preserve_order(monkeypatch):
    class _FakeExecutor:
        observed_batch_count = None

        def __init__(self, *, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def map(self, fn, iterable):
            batches = list(iterable)
            _FakeExecutor.observed_batch_count = len(batches)
            return map(fn, batches)

    regions = [(idx * 10, idx * 10 + 5, "+") for idx in range(100)]
    monkeypatch.setattr(
        load_processed.utils,
        "regions_dict_from_input",
        lambda *_args, **_kwargs: {"chr1": regions},
    )
    monkeypatch.setattr(
        load_processed.utils, "cores_to_run", lambda _cores: int(_cores)
    )
    monkeypatch.setattr(
        load_processed.concurrent.futures, "ProcessPoolExecutor", _FakeExecutor
    )

    result = load_processed.regions_to_list(
        function_handle=lambda regions, **_kwargs: f"ok:{regions}",
        regions="ignored",
        quiet=True,
        cores=2,
    )

    expected = [f"ok:chr1:{start}-{end},+" for start, end, _strand in regions]
    assert result == expected
    assert _FakeExecutor.observed_batch_count is not None
    assert _FakeExecutor.observed_batch_count < len(regions)


def test_regions_to_list_uses_external_executor_when_provided(monkeypatch):
    class _ExternalExecutor:
        called = False

        def map(self, fn, iterable):
            _ExternalExecutor.called = True
            return map(fn, iterable)

    def _raise_if_called(*_args, **_kwargs):
        raise AssertionError(
            "ProcessPoolExecutor should not be created when external executor is provided."
        )

    monkeypatch.setattr(
        load_processed.utils,
        "regions_dict_from_input",
        lambda *_args, **_kwargs: {"chr1": [(100, 110, "+"), (120, 130, "-")]},
    )
    monkeypatch.setattr(
        load_processed.utils, "cores_to_run", lambda _cores: int(_cores)
    )
    monkeypatch.setattr(
        load_processed.concurrent.futures, "ProcessPoolExecutor", _raise_if_called
    )

    executor = _ExternalExecutor()
    result = load_processed.regions_to_list(
        function_handle=lambda regions, **_kwargs: f"ok:{regions}",
        regions="ignored",
        quiet=True,
        cores=4,
        executor=executor,
    )

    assert result == ["ok:chr1:100-110,+", "ok:chr1:120-130,-"]
    assert _ExternalExecutor.called


def test_pileup_counts_from_bedmethyl_single_core_bypasses_process_pool(monkeypatch):
    def _raise_if_called(*_args, **_kwargs):
        raise AssertionError(
            "ProcessPoolExecutor should not be created for single-core execution."
        )

    def fake_process_pileup_row(row, **_kwargs):
        if row == "row_a":
            return True, 0, 2, 5
        if row == "row_b":
            return True, 0, 1, 3
        raise AssertionError(f"Unexpected row payload: {row}")

    monkeypatch.setattr(
        load_processed.utils,
        "regions_dict_from_input",
        lambda *_args, **_kwargs: {"chr1": [(100, 104, "+")]},
    )
    monkeypatch.setattr(
        load_processed.utils,
        "process_chunks_from_regions_dict",
        lambda *_args, **_kwargs: _single_region_chunks(),
    )
    monkeypatch.setattr(load_processed.utils, "cores_to_run", lambda _cores: 1)
    monkeypatch.setattr(
        load_processed.concurrent.futures, "ProcessPoolExecutor", _raise_if_called
    )
    monkeypatch.setattr(load_processed.pysam, "TabixFile", _FakeTabix)
    monkeypatch.setattr(load_processed, "process_pileup_row", fake_process_pileup_row)

    modified, valid = load_processed.pileup_counts_from_bedmethyl(
        bedmethyl_file="fake.bed.gz",
        motif="A,0",
        regions="chr1:100-104,+",
        quiet=True,
    )

    assert modified == 3
    assert valid == 8


def test_read_vectors_from_hdf5_mod_fraction_order_and_missing_defaults(tmp_path):
    h5_path = tmp_path / "reads.h5"
    with h5py.File(h5_path, "w") as h5:
        h5.create_dataset("read_name", data=np.array(["r1", "r1", "r2"], dtype="S2"))
        h5.create_dataset(
            "chromosome", data=np.array(["chr1", "chr1", "chr1"], dtype="S4")
        )
        h5.create_dataset("read_start", data=np.array([0, 0, 10], dtype=np.int32))
        h5.create_dataset("read_end", data=np.array([4, 4, 14], dtype=np.int32))
        h5.create_dataset("motif", data=np.array(["A,0", "CG,0", "A,0"], dtype="S4"))
        h5.create_dataset("strand", data=np.array(["+", "+", "+"], dtype="S1"))
        h5.create_dataset(
            "mod_vector",
            data=np.array(
                [
                    [1, 0],
                    [1, 1],
                    [0, 1],
                ],
                dtype=np.uint8,
            ),
        )
        h5.create_dataset(
            "val_vector",
            data=np.array(
                [
                    [1, 1],
                    [1, 1],
                    [1, 1],
                ],
                dtype=np.uint8,
            ),
        )

    reads, datasets, _regions_dict = load_processed.read_vectors_from_hdf5(
        file=h5_path,
        motifs=["CG,0", "A,0"],
        regions=None,
        sort_by=["read_name", "motif"],
        calculate_mod_fractions=True,
    )

    assert datasets[-2:] == ["CG,0_mod_fraction", "A,0_mod_fraction"]

    dataset_idx = {name: i for i, name in enumerate(datasets)}
    r1_rows = [row for row in reads if row[dataset_idx["read_name"]] == "r1"]
    r2_rows = [row for row in reads if row[dataset_idx["read_name"]] == "r2"]

    assert len(r1_rows) == 2
    assert len(r2_rows) == 1

    for row in r1_rows:
        assert row[dataset_idx["CG,0_mod_fraction"]] == 1.0
        assert row[dataset_idx["A,0_mod_fraction"]] == 0.5

    r2_row = r2_rows[0]
    assert r2_row[dataset_idx["CG,0_mod_fraction"]] == 0.0
    assert r2_row[dataset_idx["A,0_mod_fraction"]] == 0.5
