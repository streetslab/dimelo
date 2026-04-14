import threading
import time
from pathlib import Path

import numpy as np

from dimelo import region_analysis
from dimelo.models import SampleSpec


def _samples() -> list[SampleSpec]:
    return [
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5", metadata={"pileup_path": "s1.bed.gz"}),
        SampleSpec(sample_id="s2", condition="TR", extract_h5="s2.h5", metadata={"pileup_path": "s2.bed.gz"}),
    ]


def test_build_region_feature_table_parallelizes_across_samples(monkeypatch):
    started = 0
    started_lock = threading.Lock()
    second_started = threading.Event()

    def fake_region_feature_matrix_from_pileup(*, bedmethyl_file, motif, regions, cores, **kwargs):
        nonlocal started
        with started_lock:
            started += 1
            local_started = started
            if started >= 2:
                second_started.set()

        if local_started == 1:
            # If sample processing is serial, this wait will timeout and inflate wall time.
            second_started.wait(timeout=1.0)
        time.sleep(0.05)

        if str(bedmethyl_file).endswith("s1.bed.gz"):
            return np.array([[1.0, 2.0]], dtype=float), [("chr1", 0, 10, "+")]
        return np.array([[3.0, 4.0]], dtype=float), [("chr1", 10, 20, "-")]

    monkeypatch.setattr(
        region_analysis.cluster,
        "region_feature_matrix_from_pileup",
        fake_region_feature_matrix_from_pileup,
    )

    start = time.perf_counter()
    matrix, metadata_rows = region_analysis.build_region_feature_table(
        samples=_samples(),
        motifs=["A,0"],
        matched_regions=Path("regions.bed"),
        cores=2,
    )
    elapsed = time.perf_counter() - start

    assert elapsed < 0.5
    assert matrix.shape == (2, 2)
    assert matrix.tolist() == [[1.0, 2.0], [3.0, 4.0]]
    assert metadata_rows[0]["sample_id"] == "s1"
    assert metadata_rows[1]["sample_id"] == "s2"


def test_build_region_feature_table_reuses_regions_executor_in_serial_sample_mode(monkeypatch):
    class _FakeProcessPoolExecutor:
        instances = []

        def __init__(self, *, max_workers):
            self.max_workers = max_workers
            self.shutdown_called = False
            _FakeProcessPoolExecutor.instances.append(self)

        def shutdown(self, wait=True, cancel_futures=False):
            self.shutdown_called = True

    observed_executor_ids = []

    def fake_region_feature_matrix_from_pileup(
        *, bedmethyl_file, motif, regions, cores, regions_executor=None, **kwargs
    ):
        observed_executor_ids.append(id(regions_executor))
        if str(bedmethyl_file).endswith("s1.bed.gz"):
            return np.array([[1.0, 2.0]], dtype=float), [("chr1", 0, 10, "+")]
        return np.array([[3.0, 4.0]], dtype=float), [("chr1", 10, 20, "-")]

    monkeypatch.setattr(
        region_analysis.concurrent.futures,
        "ProcessPoolExecutor",
        _FakeProcessPoolExecutor,
    )
    monkeypatch.setattr(region_analysis.utils, "cores_to_run", lambda _cores: 8)
    monkeypatch.setattr(
        region_analysis.cluster,
        "region_feature_matrix_from_pileup",
        fake_region_feature_matrix_from_pileup,
    )

    matrix, metadata_rows = region_analysis.build_region_feature_table(
        samples=_samples(),
        motifs=["A,0"],
        matched_regions=Path("regions.bed"),
        cores=None,
    )

    assert matrix.shape == (2, 2)
    assert len(metadata_rows) == 2
    assert len(_FakeProcessPoolExecutor.instances) == 1
    assert _FakeProcessPoolExecutor.instances[0].shutdown_called
    assert len(set(observed_executor_ids)) == 1
