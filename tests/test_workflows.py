import numpy as np
import pandas as pd
import pytest

from dimelo import workflows
from dimelo.models import ContrastSpec, DatasetArtifact
from dimelo.models import (
    RegionContrastResult,
    RegionDiscoveryResult,
    SampleSpec,
    SharedClusterModel,
    SharedClusterResult,
)


def _workflow_samples():
    return [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            regions_bed="control-s1.bed",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="s2.h5",
            regions_bed="control-s2.bed",
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]


def _mock_discovery_result(*args, **kwargs):
    hits = pd.DataFrame(
        [
            {
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": "+",
                "score_value": 0.9,
                "adjusted_p_value": 0.01,
            },
            {
                "window_id": "chr1:500-1000",
                "chromosome": "chr1",
                "start": 500,
                "end": 1000,
                "strand": "-",
                "score_value": 0.6,
                "adjusted_p_value": 0.02,
            },
            {
                "window_id": "chr1:1000-1500",
                "chromosome": "chr1",
                "start": 1000,
                "end": 1500,
                "strand": "+",
                "score_value": 0.2,
                "adjusted_p_value": 0.20,
            },
        ]
    )
    return RegionDiscoveryResult(
        hits=hits,
        windows=hits.copy(),
        contrast=None,
        plot_data={},
        metadata={"score": "effect_size_only"},
    )


def _mock_cluster_result(*args, **kwargs):
    assignments = pd.DataFrame(
        [
            {"sample_id": "s1", "condition": "NS", "cluster": "C0"},
            {"sample_id": "s2", "condition": "15min", "cluster": "C1"},
        ]
    )
    return SharedClusterResult(
        model=SharedClusterModel(
            mode=kwargs["mode"],
            motifs=list(kwargs["motifs"]),
            feature_names=["f0"],
            preprocessing={},
            estimator=object(),
            cluster_labels=["C0", "C1"],
            fit_metadata={},
        ),
        assignments=assignments,
        cluster_distribution=pd.DataFrame(
            [
                {"sample_id": "s1", "condition": "NS", "cluster": "C0", "count": 1, "fraction": 1.0},
                {"sample_id": "s2", "condition": "15min", "cluster": "C1", "count": 1, "fraction": 1.0},
            ]
        ),
        condition_distribution=pd.DataFrame(
            [
                {"condition": "NS", "cluster": "C0", "fraction": 1.0, "replicate_n": 1},
                {"condition": "15min", "cluster": "C1", "fraction": 1.0, "replicate_n": 1},
            ]
        ),
        distribution_change=None,
        cluster_profiles=pd.DataFrame([{"cluster": "C0", "count": 1, "f0": 0.0}]),
        region_summaries=None,
        plot_data={},
        metadata={"matched_regions": kwargs.get("matched_regions")},
    )


def _mock_region_contrast_result(*args, **kwargs):
    contrast = kwargs["contrast"]
    regions = pd.DataFrame(
        [
            {
                "region_id": "chr1:0-500,+",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": "+",
                "fraction": 0.8,
                "reference_fraction": 0.2,
                "delta_fraction": 0.6,
                "log2_fc": 2.0,
                "effect_size": 0.6,
                "rank": 1,
                "numerator_modified_count": 8,
                "numerator_valid_count": 10,
                "numerator_replicate_n": 1,
                "denominator_modified_count": 2,
                "denominator_valid_count": 10,
                "denominator_replicate_n": 1,
            },
            {
                "region_id": "chr1:500-1000,-",
                "chromosome": "chr1",
                "start": 500,
                "end": 1000,
                "strand": "-",
                "fraction": 0.6,
                "reference_fraction": 0.3,
                "delta_fraction": 0.3,
                "log2_fc": 1.0,
                "effect_size": 0.3,
                "rank": 2,
                "numerator_modified_count": 6,
                "numerator_valid_count": 10,
                "numerator_replicate_n": 1,
                "denominator_modified_count": 3,
                "denominator_valid_count": 10,
                "denominator_replicate_n": 1,
            },
        ]
    )
    summary = regions[
        [
            "region_id",
            "fraction",
            "reference_fraction",
            "delta_fraction",
            "log2_fc",
            "rank",
            "numerator_modified_count",
            "numerator_valid_count",
            "numerator_replicate_n",
            "denominator_modified_count",
            "denominator_valid_count",
            "denominator_replicate_n",
        ]
    ].copy()
    return RegionContrastResult(
        regions=regions,
        summary=summary,
        contrast=contrast,
        plot_data={"region_effect_sizes": summary.copy()},
        metadata={
            "analysis_unit": kwargs.get("analysis_unit", "ensemble_region"),
            "representation": kwargs.get("representation", "modified_fraction"),
            "signal_source": kwargs.get("signal_source", "pileup_counts"),
            "test": kwargs.get("test", "effect_size_only"),
            "multiple_testing": kwargs.get("multiple_testing", "fdr_bh"),
        },
    )


def _fake_region_anchored_extract(*args, **kwargs):
    regions = kwargs["regions"]
    region_list = list(regions) if isinstance(regions, list) else [regions]
    metadata = []
    for region in region_list:
        chrom_coords, strand = region.split(",")
        chrom, coords = chrom_coords.split(":")
        start, end = coords.split("-")
        metadata.append((chrom, int(start), int(end), strand))
    matrix = np.array([[0.2, 0.8], [0.8, 0.2]], dtype=np.float32)
    return matrix[: len(metadata)], metadata


def test_discovery_cluster_workflow_returns_both_results(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    result = workflows.discovery_cluster_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
    )

    assert result.discovery.hits.shape[0] == 3
    assert result.clustering.model.mode == "region_anchored"
    assert list(result.selected_regions.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    assert result.metadata["selection"]["mode"] == "top_n"
    assert result.metadata["selection"]["top_n"] == 250


def test_discovery_cluster_workflow_selects_top_n_hits(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    result = workflows.discovery_cluster_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        selection={"mode": "top_n", "top_n": 1},
    )

    assert result.selected_regions["name"].tolist() == ["chr1:0-500"]


def test_discovery_cluster_workflow_selects_all_hits(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    result = workflows.discovery_cluster_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        selection={"mode": "all"},
    )

    assert set(result.selected_regions["name"].tolist()) == {
        "chr1:0-500",
        "chr1:500-1000",
        "chr1:1000-1500",
    }
    assert result.metadata["selection"]["top_n"] is None


def test_discovery_cluster_workflow_passes_selected_regions_into_clustering(monkeypatch):
    captured = {}

    def fake_cluster(**kwargs):
        captured["matched_regions"] = kwargs["matched_regions"]
        return _mock_cluster_result(**kwargs)

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", fake_cluster)

    result = workflows.discovery_cluster_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        selection={"mode": "top_n", "top_n": 2},
    )

    assert captured["matched_regions"] == ["chr1:0-500,+", "chr1:500-1000,-"]
    assert list(result.selected_regions.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    assert list(result.selected_regions["name"]) == ["chr1:0-500", "chr1:500-1000"]


def test_discovery_cluster_workflow_errors_when_no_hits_survive_selection(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    with pytest.raises(ValueError, match="No discovery hits remained after selection"):
        workflows.discovery_cluster_workflow(
            samples=_workflow_samples(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1500},
            discovery={"window_size": 500, "step_size": 500},
            clustering={"mode": "region_anchored", "n_clusters": 2},
            selection={"mode": "top_n", "top_n": 0},
        )


def test_discovery_cluster_workflow_rejects_unknown_selection_mode(monkeypatch):
    called = {"scan_genome": False}

    def fake_scan_genome(*args, **kwargs):
        called["scan_genome"] = True
        raise AssertionError("scan_genome should not be called for invalid selection config")

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", fake_scan_genome)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    with pytest.raises(ValueError, match="Unsupported selection mode"):
        workflows.discovery_cluster_workflow(
            samples=_workflow_samples(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1500},
            discovery={"window_size": 500, "step_size": 500},
            clustering={"mode": "region_anchored", "n_clusters": 2},
            selection={"mode": "invalid"},
        )
    assert called["scan_genome"] is False


def test_discovery_cluster_workflow_rejects_invalid_clustering_config_before_scan(monkeypatch):
    called = {"scan_genome": False}

    def fake_scan_genome(*args, **kwargs):
        called["scan_genome"] = True
        raise AssertionError("scan_genome should not be called for invalid clustering config")

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", fake_scan_genome)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    with pytest.raises(NotImplementedError, match="clusterer='minibatch_kmeans' only"):
        workflows.discovery_cluster_workflow(
            samples=_workflow_samples(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1500},
            discovery={"window_size": 500, "step_size": 500},
            clustering={"mode": "region_anchored", "clusterer": "agglomerative", "n_clusters": 2},
        )
    assert called["scan_genome"] is False


def test_discovery_cluster_workflow_region_anchored_uses_serializable_matched_regions(monkeypatch):
    captured = {"resolve_params": []}

    def fake_resolve_artifact(requested_artifact, available_artifacts, artifact_policy):
        captured["resolve_params"].append(requested_artifact.params["matched_regions"])
        return None

    def fake_region_features(*args, **kwargs):
        regions = kwargs["regions"]
        captured["regions"] = regions
        return _fake_region_anchored_extract(*args, **kwargs)

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "resolve_artifact", fake_resolve_artifact)
    monkeypatch.setattr(
        workflows.cluster,
        "region_feature_matrix_from_pileup",
        fake_region_features,
    )

    result = workflows.discovery_cluster_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2, "make_plots": False},
        selection={"mode": "top_n", "top_n": 2},
    )

    assert captured["regions"] == ["chr1:0-500,+", "chr1:500-1000,-"]
    assert captured["resolve_params"] == [
        ["chr1:0-500,+", "chr1:500-1000,-"],
        ["chr1:0-500,+", "chr1:500-1000,-"],
    ]
    assert result.clustering.metadata["matched_regions"] == ["chr1:0-500,+", "chr1:500-1000,-"]


def test_discovery_cluster_workflow_artifact_params_keep_matched_regions_json_serializable(monkeypatch):
    captured = {}
    region_spec = ["chr1:0-500,+", "chr1:500-1000,-"]
    artifact_samples = []
    for sample in _workflow_samples():
        artifact_samples.append(
            SampleSpec(
                sample_id=sample.sample_id,
                condition=sample.condition,
                extract_h5=sample.extract_h5,
                regions_bed=sample.regions_bed,
                metadata={
                    "pileup_path": sample.metadata["pileup_path"],
                    "artifacts": [
                        DatasetArtifact(
                            sample_id=sample.sample_id,
                            artifact_type="pileup",
                            path=sample.metadata["pileup_path"],
                            format="bed.gz",
                            params={
                                "motifs": ["A,0"],
                                "matched_regions": region_spec,
                                "signal_normalization": "none",
                                "feature_scaling": "robust_zscore",
                                "cluster_basis": "shape_plus_level",
                            },
                            provenance={"pipeline": "parse_bam", "source_files": [], "source_fingerprints": []},
                        )
                    ],
                },
            )
        )

    def fake_resolve_artifact(requested_artifact, available_artifacts, artifact_policy):
        captured.setdefault("params", []).append(requested_artifact.params["matched_regions"])
        return available_artifacts[0]

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "resolve_artifact", fake_resolve_artifact)
    monkeypatch.setattr(
        workflows.cluster,
        "region_feature_matrix_from_pileup",
        _fake_region_anchored_extract,
    )

    result = workflows.discovery_cluster_workflow(
        samples=artifact_samples,
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2, "make_plots": False},
        selection={"mode": "top_n", "top_n": 2},
    )

    assert captured["params"] == [region_spec, region_spec]
    assert result.clustering.metadata["cache_hits"] == {"s1": "s1.bed.gz", "s2": "s2.bed.gz"}


def test_discovery_cluster_workflow_materializes_samples_iterable(monkeypatch):
    captured = {}

    def fake_discovery(*args, **kwargs):
        captured["discovery_sample_ids"] = [sample.sample_id for sample in kwargs["samples"]]
        return _mock_discovery_result(*args, **kwargs)

    def fake_cluster(**kwargs):
        captured["sample_ids"] = [sample.sample_id for sample in kwargs["samples"]]
        return _mock_cluster_result(**kwargs)

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", fake_discovery)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", fake_cluster)

    result = workflows.discovery_cluster_workflow(
        samples=(sample for sample in _workflow_samples()),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        selection={"mode": "top_n", "top_n": 2},
    )

    assert result.discovery.hits.shape[0] == 3
    assert result.clustering.model.mode == "region_anchored"
    assert captured["discovery_sample_ids"] == ["s1", "s2"]
    assert captured["sample_ids"] == ["s1", "s2"]
    assert list(result.selected_regions["name"]) == ["chr1:0-500", "chr1:500-1000"]


def test_discovery_cluster_workflow_materializes_motif_generator(monkeypatch):
    captured = {}

    def fake_discovery(*args, **kwargs):
        captured["discovery_motifs"] = list(kwargs["motifs"])
        return _mock_discovery_result(*args, **kwargs)

    def fake_cluster(**kwargs):
        captured["clustering_motifs"] = list(kwargs["motifs"])
        return _mock_cluster_result(**kwargs)

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", fake_discovery)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", fake_cluster)

    workflows.discovery_cluster_workflow(
        samples=_workflow_samples(),
        motifs=(motif for motif in ["A,0", "CG,1"]),
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        selection={"mode": "top_n", "top_n": 2},
    )

    assert captured["discovery_motifs"] == ["A,0", "CG,1"]
    assert captured["clustering_motifs"] == ["A,0", "CG,1"]


def test_discovery_cluster_contrast_workflow_returns_all_results(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)
    monkeypatch.setattr(workflows.region_contrasts, "score_regions", _mock_region_contrast_result)

    result = workflows.discovery_cluster_contrast_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        contrasts={
            "contrast": ContrastSpec(
                mode="pairwise",
                numerator=["15min"],
                denominator=["NS"],
            ),
            "test": "effect_size_only",
        },
    )

    assert result.discovery.hits.shape[0] == 3
    assert result.clustering.model.mode == "region_anchored"
    assert result.contrasts.metadata["test"] == "effect_size_only"
    assert result.metadata["contrast_scope"] == "selected"
    assert list(result.selected_regions.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    assert result.metadata["full_scan_windows"].equals(result.discovery.windows)


def test_discovery_cluster_contrast_workflow_scores_selected_regions_by_default(monkeypatch):
    captured = {}

    def fake_score_regions(**kwargs):
        captured["regions"] = kwargs["regions"]
        return _mock_region_contrast_result(**kwargs)

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)
    monkeypatch.setattr(workflows.region_contrasts, "score_regions", fake_score_regions)

    result = workflows.discovery_cluster_contrast_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        contrasts={
            "contrast": ContrastSpec(
                mode="pairwise",
                numerator=["15min"],
                denominator=["NS"],
            ),
        },
        selection={"mode": "top_n", "top_n": 2},
    )

    assert captured["regions"] == ["chr1:0-500,+", "chr1:500-1000,-"]
    assert result.selected_regions["name"].tolist() == ["chr1:0-500", "chr1:500-1000"]


def test_discovery_cluster_contrast_workflow_normalizes_clustering_region_ids_for_joining(
    monkeypatch,
):
    def fake_cluster(**kwargs):
        result = _mock_cluster_result(**kwargs)
        result.assignments = pd.DataFrame(
            [
                {
                    "region_id": "chr1:0-500",
                    "sample_id": "s1",
                    "condition": "NS",
                    "cluster": "C0",
                    "strand": "+",
                },
                {
                    "region_id": "chr1:500-1000",
                    "sample_id": "s2",
                    "condition": "15min",
                    "cluster": "C1",
                    "strand": "-",
                },
            ]
        )
        result.region_summaries = pd.DataFrame(
            [
                {
                    "region_id": "chr1:0-500",
                    "sample_id": "s1",
                    "condition": "NS",
                    "cluster": "C0",
                    "count": 1,
                    "fraction": 1.0,
                    "strand": "+",
                },
                {
                    "region_id": "chr1:500-1000",
                    "sample_id": "s2",
                    "condition": "15min",
                    "cluster": "C1",
                    "count": 1,
                    "fraction": 1.0,
                    "strand": "-",
                },
            ]
        )
        return result

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", fake_cluster)
    monkeypatch.setattr(workflows.region_contrasts, "score_regions", _mock_region_contrast_result)

    result = workflows.discovery_cluster_contrast_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        contrasts={
            "contrast": ContrastSpec(
                mode="pairwise",
                numerator=["15min"],
                denominator=["NS"],
            ),
        },
        selection={"mode": "top_n", "top_n": 2},
    )

    contrast_region_ids = set(result.contrasts.regions["region_id"])
    assert set(result.clustering.assignments["region_id"]) == contrast_region_ids
    assert set(result.clustering.region_summaries["region_id"]) == contrast_region_ids


def test_discovery_cluster_contrast_workflow_uses_custom_contrast_regions_when_provided(
    monkeypatch,
):
    captured = {}

    def fake_score_regions(**kwargs):
        captured["regions"] = kwargs["regions"]
        return _mock_region_contrast_result(**kwargs)

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)
    monkeypatch.setattr(workflows.region_contrasts, "score_regions", fake_score_regions)

    custom_regions = ["chr2:0-100,+", "chr2:100-200,-"]
    result = workflows.discovery_cluster_contrast_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        contrasts={
            "contrast": ContrastSpec(
                mode="pairwise",
                numerator=["15min"],
                denominator=["NS"],
            ),
            "regions": custom_regions,
        },
        selection={"mode": "top_n", "top_n": 2},
    )

    assert captured["regions"] == custom_regions
    assert result.metadata["contrast_scope"] == "custom"
    assert result.selected_regions["name"].tolist() == ["chr1:0-500", "chr1:500-1000"]


def test_discovery_cluster_contrast_workflow_preserves_full_scan_windows_context(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)
    monkeypatch.setattr(workflows.region_contrasts, "score_regions", _mock_region_contrast_result)

    result = workflows.discovery_cluster_contrast_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        contrasts={
            "contrast": ContrastSpec(
                mode="pairwise",
                numerator=["15min"],
                denominator=["NS"],
            ),
        },
    )

    assert result.metadata["full_scan_windows"].equals(result.discovery.windows)


def test_discovery_cluster_contrast_workflow_rejects_missing_contrast_config(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    with pytest.raises(ValueError, match=r"requires contrasts\['contrast'\]"):
        workflows.discovery_cluster_contrast_workflow(
            samples=_workflow_samples(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1500},
            discovery={"window_size": 500, "step_size": 500},
            clustering={"mode": "region_anchored", "n_clusters": 2},
            contrasts={"test": "effect_size_only"},
        )


def test_discovery_cluster_contrast_workflow_fast_fails_invalid_contrast_config_before_scan(
    monkeypatch,
):
    called = {"scan_genome": False}

    def fake_scan_genome(*args, **kwargs):
        called["scan_genome"] = True
        raise AssertionError("scan_genome should not be called for invalid contrast config")

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", fake_scan_genome)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    with pytest.raises(ValueError, match="analysis_unit='ensemble_region'"):
        workflows.discovery_cluster_contrast_workflow(
            samples=_workflow_samples(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1500},
            discovery={"window_size": 500, "step_size": 500},
            clustering={"mode": "region_anchored", "n_clusters": 2},
            contrasts={
                "contrast": ContrastSpec(
                    mode="pairwise",
                    numerator=["15min"],
                    denominator=["NS"],
                ),
                "analysis_unit": "single_read",
            },
        )

    assert called["scan_genome"] is False


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
    with pytest.raises(NotImplementedError, match="region_anchored"):
        workflows.shared_cluster_distribution(
            samples=[
                SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5"),
                SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5"),
            ],
            mode="unknown_mode",
            motifs=["A,0"],
        )


def test_shared_cluster_distribution_requires_two_datasets():
    with pytest.raises(ValueError, match="at least two datasets"):
        workflows.shared_cluster_distribution(
            samples=[SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5")],
            mode="read_global",
            motifs=["A,0"],
        )


def test_shared_cluster_distribution_region_anchored(monkeypatch):
    fake_samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            regions_bed="r1.bed",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="s2.h5",
            regions_bed="r2.bed",
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]

    def fake_region_table(*args, **kwargs):
        return np.array([[0.2, 0.8], [0.7, 0.3]]), [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "chromosome": "chr1",
                "start": 0,
                "end": 2,
                "strand": "+",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "15min",
                "replicate": None,
                "chromosome": "chr1",
                "start": 0,
                "end": 2,
                "strand": "+",
            },
        ]

    monkeypatch.setattr(workflows.region_analysis, "build_region_feature_table", fake_region_table)

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="region_anchored",
        motifs=["A,0"],
        matched_regions="matched.bed",
        n_clusters=2,
        make_plots=False,
    )

    assert not result.assignments.empty
    assert "region_id" in result.assignments.columns
    assert result.region_summaries is not None
    assert {"region_id", "sample_id", "condition", "cluster", "count", "fraction"} <= set(
        result.region_summaries.columns
    )


def test_shared_cluster_distribution_region_anchored_requires_matched_regions():
    with pytest.raises(ValueError, match="requires matched_regions"):
        workflows.shared_cluster_distribution(
            samples=_workflow_samples(),
            mode="region_anchored",
            motifs=["A,0"],
            make_plots=False,
        )


def test_shared_cluster_distribution_region_anchored_accepts_list_matched_regions(monkeypatch):
    captured = {}

    def fake_region_table(*args, **kwargs):
        captured["matched_regions"] = kwargs["matched_regions"]
        return np.array([[0.2, 0.8], [0.7, 0.3]]), [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "chromosome": "chr1",
                "start": 0,
                "end": 2,
                "strand": "+",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "15min",
                "replicate": None,
                "chromosome": "chr1",
                "start": 0,
                "end": 2,
                "strand": "+",
            },
        ]

    monkeypatch.setattr(workflows.region_analysis, "build_region_feature_table", fake_region_table)

    result = workflows.shared_cluster_distribution(
        samples=_workflow_samples(),
        mode="region_anchored",
        motifs=["A,0"],
        matched_regions=["chr1:0-2,+", "chr1:2-4,-"],
        n_clusters=2,
        make_plots=False,
    )

    assert captured["matched_regions"] == ["chr1:0-2,+", "chr1:2-4,-"]
    assert not result.assignments.empty


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


def test_shared_cluster_distribution_region_anchored_rejects_multi_motif(monkeypatch):
    fake_samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            regions_bed="r1.bed",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="s2.h5",
            regions_bed="r2.bed",
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]

    with pytest.raises(ValueError, match="exactly one motif"):
        workflows.shared_cluster_distribution(
            samples=fake_samples,
            mode="region_anchored",
            motifs=["A,0", "CG,0"],
            matched_regions="matched.bed",
        )


def test_shared_cluster_distribution_region_anchored_control_regions(monkeypatch):
    fake_samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            regions_bed="control-s1.bed",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="s2.h5",
            regions_bed="control-s2.bed",
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]

    def fake_region_table(*args, **kwargs):
        matched_regions = kwargs["matched_regions"]
        if matched_regions == "matched.bed":
            return np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32), [
                {
                    "region_id": "reg1",
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": None,
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 2,
                    "strand": "+",
                },
                {
                    "region_id": "reg1",
                    "sample_id": "s2",
                    "condition": "15min",
                    "replicate": None,
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 2,
                    "strand": "+",
                },
            ]
        if matched_regions == "control-s1.bed":
            return np.array([[0.5, 0.5]], dtype=np.float32), [
                {
                    "region_id": "c1",
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": None,
                    "chromosome": "chr1",
                    "start": 10,
                    "end": 12,
                    "strand": "+",
                }
            ]
        return np.array([[0.25, 0.25]], dtype=np.float32), [
            {
                "region_id": "c2",
                "sample_id": "s2",
                "condition": "15min",
                "replicate": None,
                "chromosome": "chr1",
                "start": 20,
                "end": 22,
                "strand": "+",
            }
        ]

    monkeypatch.setattr(workflows.region_analysis, "build_region_feature_table", fake_region_table)

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="region_anchored",
        motifs=["A,0"],
        matched_regions="matched.bed",
        signal_normalization="control_regions",
        n_clusters=2,
    )

    assert result.metadata["sample_normalization"]["s1"]["global_offset"] == 0.5
    assert result.metadata["sample_normalization"]["s2"]["global_offset"] == 0.25


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
