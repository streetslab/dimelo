import numpy as np
import pandas as pd
import pytest

from dimelo import workflows
from dimelo.models import DatasetArtifact
from dimelo.models import SampleSpec


def test_shared_cluster_distribution_read_global(monkeypatch):
    fake_samples = [
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5"),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5"),
    ]

    def fake_extract(*args, **kwargs):
        hdf5_file = kwargs["hdf5_file"]

        class R:
            if hdf5_file == "s1.h5":
                data_matrix = np.array([[0.0, 0.0, 0.1, 0.1], [0.0, 0.1, 0.0, 0.1]])
                metadata = [
                    {"read_name": "s1-r1", "chromosome": "chr1", "region_start": 0, "region_end": 4},
                    {"read_name": "s1-r2", "chromosome": "chr1", "region_start": 10, "region_end": 14},
                ]
            else:
                data_matrix = np.array([[1.0, 0.9, 1.0, 0.9], [0.9, 1.0, 0.9, 1.0]])
                metadata = [
                    {"read_name": "s2-r1", "chromosome": "chr1", "region_start": 20, "region_end": 24},
                    {"read_name": "s2-r2", "chromosome": "chr1", "region_start": 30, "region_end": 34},
                ]
            val_matrix = np.ones((2, 4), dtype=float)
            datasets = []
            regions_dict = None

        return R()

    def fake_features(result, **kwargs):
        return result.data_matrix, ["f0", "f1", "f2", "f3"]

    monkeypatch.setattr(workflows.cluster, "extract_read_windows", fake_extract)
    monkeypatch.setattr(workflows.cluster, "read_window_feature_matrix", fake_features)

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="read_global",
        motifs=["A,0"],
        n_clusters=2,
        training_sample_per_dataset=2,
        make_plots=False,
    )

    assert not result.assignments.empty
    assert not result.cluster_distribution.empty
    assert not result.condition_distribution.empty
    assert not result.cluster_profiles.empty
    assert set(result.assignments["sample_id"]) == {"s1", "s2"}
    assert "replicate" in result.assignments.columns
    assert "cluster_distribution_bar" in result.plot_data
    assert "cluster_distribution_heatmap" in result.plot_data
    assert result.metadata["artifact_policy"] == "prefer_cached"
    assert set(result.metadata["cache_misses"]) == {"s1", "s2"}


def test_shared_cluster_distribution_rejects_unsupported_mode():
    with pytest.raises(NotImplementedError, match="read_global"):
        workflows.shared_cluster_distribution(
            samples=[SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5")],
            mode="region_anchored",
            motifs=["A,0"],
        )


def test_shared_cluster_distribution_requires_two_datasets():
    with pytest.raises(ValueError, match="at least two datasets"):
        workflows.shared_cluster_distribution(
            samples=[SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5")],
            mode="read_global",
            motifs=["A,0"],
        )


def test_shared_cluster_distribution_tracks_condition_replicates(monkeypatch):
    fake_samples = [
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5"),
        SampleSpec(sample_id="s2", condition="NS", extract_h5="s2.h5"),
        SampleSpec(sample_id="s3", condition="15min", extract_h5="s3.h5"),
    ]

    matrices = {
        "s1.h5": np.array([[0.0, 0.0], [0.0, 0.1]]),
        "s2.h5": np.array([[0.1, 0.0], [0.0, 0.0]]),
        "s3.h5": np.array([[1.0, 1.0], [0.9, 1.0]]),
    }

    def fake_extract(*args, **kwargs):
        data = matrices[kwargs["hdf5_file"]]

        class R:
            data_matrix = data
            val_matrix = np.ones_like(data)
            metadata = [
                {"read_name": f"{kwargs['hdf5_file']}-r1", "chromosome": "chr1", "region_start": 0, "region_end": 2},
                {"read_name": f"{kwargs['hdf5_file']}-r2", "chromosome": "chr1", "region_start": 2, "region_end": 4},
            ]
            datasets = []
            regions_dict = None

        return R()

    monkeypatch.setattr(workflows.cluster, "extract_read_windows", fake_extract)
    monkeypatch.setattr(
        workflows.cluster,
        "read_window_feature_matrix",
        lambda result, **kwargs: (result.data_matrix, ["f0", "f1"]),
    )

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="read_global",
        motifs=["A,0"],
        n_clusters=2,
        training_sample_per_dataset=2,
    )

    ns_rows = result.condition_distribution[result.condition_distribution["condition"] == "NS"]
    assert set(ns_rows["replicate_n"]) == {2}


def test_shared_cluster_distribution_prefers_cached_extract_artifact(monkeypatch):
    cached_artifact = DatasetArtifact(
        sample_id="s1",
        artifact_type="extract",
        path="cached-s1.h5",
        format="h5",
        params={
            "motifs": ["A,0"],
            "matched_regions": None,
            "window_size": None,
            "signal_normalization": "none",
            "feature_scaling": "robust_zscore",
            "cluster_basis": "shape_plus_level",
        },
        provenance={
            "pipeline": "parse_bam",
            "source_files": ["s1.h5"],
            "source_fingerprints": [{"path": "s1.h5", "missing": True}],
            "upstream_lineage": [],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.0.0"},
    )
    fake_samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"artifacts": [cached_artifact]},
        ),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5"),
    ]
    seen_paths = []

    def fake_extract(*args, **kwargs):
        seen_paths.append(kwargs["hdf5_file"])

        class R:
            data_matrix = np.array([[0.0, 0.0], [1.0, 1.0]])
            val_matrix = np.ones((2, 2), dtype=float)
            metadata = [
                {"read_name": "r1", "chromosome": "chr1", "region_start": 0, "region_end": 2},
                {"read_name": "r2", "chromosome": "chr1", "region_start": 2, "region_end": 4},
            ]
            datasets = []
            regions_dict = None

        return R()

    monkeypatch.setattr(workflows.cluster, "extract_read_windows", fake_extract)
    monkeypatch.setattr(
        workflows.cluster,
        "read_window_feature_matrix",
        lambda result, **kwargs: (result.data_matrix, ["f0", "f1"]),
    )

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="read_global",
        motifs=["A,0"],
        n_clusters=2,
        training_sample_per_dataset=2,
    )

    assert seen_paths[0] == "cached-s1.h5"
    assert result.metadata["cache_hits"]["s1"] == "cached-s1.h5"


def test_shared_cluster_distribution_rebuilds_when_source_fingerprint_mismatches(monkeypatch):
    stale_artifact = DatasetArtifact(
        sample_id="s1",
        artifact_type="extract",
        path="cached-s1.h5",
        format="h5",
        params={
            "motifs": ["A,0"],
            "matched_regions": None,
            "window_size": None,
            "signal_normalization": "none",
            "feature_scaling": "robust_zscore",
            "cluster_basis": "shape_plus_level",
        },
        provenance={
            "pipeline": "parse_bam",
            "source_files": ["s1.h5"],
            "source_fingerprints": [{"path": "s1.h5", "size": 123, "mtime_ns": 1}],
            "upstream_lineage": [],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.0.0"},
    )
    fake_samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"artifacts": [stale_artifact]},
        ),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5"),
    ]
    seen_paths = []

    def fake_extract(*args, **kwargs):
        seen_paths.append(kwargs["hdf5_file"])

        class R:
            data_matrix = np.array([[0.0, 0.0], [1.0, 1.0]])
            val_matrix = np.ones((2, 2), dtype=float)
            metadata = [
                {"read_name": "r1", "chromosome": "chr1", "region_start": 0, "region_end": 2},
                {"read_name": "r2", "chromosome": "chr1", "region_start": 2, "region_end": 4},
            ]
            datasets = []
            regions_dict = None

        return R()

    monkeypatch.setattr(workflows.cluster, "extract_read_windows", fake_extract)
    monkeypatch.setattr(
        workflows.cluster,
        "read_window_feature_matrix",
        lambda result, **kwargs: (result.data_matrix, ["f0", "f1"]),
    )

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="read_global",
        motifs=["A,0"],
        n_clusters=2,
        training_sample_per_dataset=2,
    )

    assert seen_paths[0] == "s1.h5"
    assert "s1" not in result.metadata["cache_hits"]
