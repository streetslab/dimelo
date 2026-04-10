import numpy as np
import pandas as pd
import pytest

from dimelo.models import (
    ContrastSpec,
    GlobalAnalysisResult,
    RegionContrastResult,
    RegionDiscoveryResult,
    SharedClusterModel,
    SharedClusterResult,
)
from dimelo import plot_enrichment_profile, plot_reads, plotting


def test_axis_spec_accepts_fixed_window_region_5to3():
    spec = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=1000,
        downstream_bp=1000,
    )

    plotting.validate_axis_spec(spec, plot_family="aggregate_profile")


def test_axis_spec_rejects_segment_map_without_segments():
    spec = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=None,
    )

    with pytest.raises(ValueError, match="segment_map requires segments"):
        plotting.validate_axis_spec(spec, plot_family="aggregate_profile")


def test_single_read_rejects_scaled_segments():
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec(
                segment_id="body",
                label="Body",
                start_ref=100,
                end_ref=500,
                mode="scaled",
                bins=50,
            )
        ],
    )

    with pytest.raises(ValueError, match="single_read"):
        plotting.validate_axis_spec(axis, plot_family="single_read_raster")


def test_aggregation_spec_accepts_equal_region_default():
    spec = plotting.AggregationSpec()

    plotting.validate_aggregation_spec(spec)


def test_prepare_cluster_distribution_bar_data_handles_empty_input_shape():
    cluster_distribution = pd.DataFrame(
        columns=["sample_id", "condition", "cluster", "count", "fraction"]
    )

    result = plotting.prepare_cluster_distribution_bar_data(cluster_distribution)

    assert list(result.columns) == ["sample_id", "condition", "cluster", "count", "fraction"]
    assert result.empty


def test_prepare_cluster_distribution_bar_data_keeps_stable_sorted_output():
    cluster_distribution = pd.DataFrame(
        {
            "sample_id": ["s2", "s1", "s1", "s2"],
            "condition": ["15min", "NS", "NS", "15min"],
            "cluster": ["C0", "C0", "C0", "C1"],
            "count": [1, 2, 3, 4],
            "fraction": [0.25, 0.66, 0.34, 0.75],
        }
    )

    result = plotting.prepare_cluster_distribution_bar_data(cluster_distribution)

    assert list(result["sample_id"]) == ["s1", "s1", "s2", "s2"]
    assert list(result["condition"]) == ["NS", "NS", "15min", "15min"]
    assert list(result["cluster"]) == ["C0", "C0", "C0", "C1"]
    assert list(result["count"]) == [2, 3, 1, 4]


def test_prepare_cluster_distribution_heatmap_data_handles_empty_input_shape():
    condition_distribution = pd.DataFrame(columns=["condition", "cluster", "fraction"])

    result = plotting.prepare_cluster_distribution_heatmap_data(condition_distribution)

    assert list(result.columns) == ["condition"]
    assert result.empty


def test_prepare_cluster_distribution_heatmap_data_pivots_columns_in_sorted_order():
    condition_distribution = pd.DataFrame(
        {
            "condition": ["NS", "NS", "15min", "15min"],
            "cluster": ["C1", "C0", "C1", "C0"],
            "fraction": [0.25, 0.75, 0.75, 0.25],
        }
    )

    result = plotting.prepare_cluster_distribution_heatmap_data(condition_distribution)

    assert list(result.columns) == ["condition", "C0", "C1"]
    assert list(result["condition"]) == ["15min", "NS"]
    fifteen_min = result[result["condition"] == "15min"].iloc[0]
    ns_row = result[result["condition"] == "NS"].iloc[0]
    assert fifteen_min["C0"] == 0.25
    assert fifteen_min["C1"] == 0.75
    assert ns_row["C0"] == 0.75
    assert ns_row["C1"] == 0.25


def _make_shared_cluster_result() -> SharedClusterResult:
    model = SharedClusterModel(
        mode="region_anchored",
        motifs=["A,0"],
        feature_names=["f0", "f1"],
        preprocessing={"signal_normalization": "none"},
        estimator=object(),
        cluster_labels=["C1", "C0"],
        fit_metadata={"clusterer": "minibatch_kmeans", "n_clusters": 2},
    )
    cluster_distribution = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "cluster": "C0",
                "count": 2,
                "fraction": 2 / 3,
            },
            {
                "sample_id": "s1",
                "condition": "NS",
                "cluster": "C1",
                "count": 1,
                "fraction": 1 / 3,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "cluster": "C0",
                "count": 1,
                "fraction": 1 / 4,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "cluster": "C1",
                "count": 3,
                "fraction": 3 / 4,
            },
        ]
    )
    condition_distribution = pd.DataFrame(
        [
            {
                "condition": "NS",
                "cluster": "C0",
                "count": 2,
                "fraction": 2 / 3,
                "replicate_n": 1,
            },
            {
                "condition": "NS",
                "cluster": "C1",
                "count": 1,
                "fraction": 1 / 3,
                "replicate_n": 1,
            },
            {
                "condition": "treated",
                "cluster": "C0",
                "count": 1,
                "fraction": 1 / 4,
                "replicate_n": 1,
            },
            {
                "condition": "treated",
                "cluster": "C1",
                "count": 3,
                "fraction": 3 / 4,
                "replicate_n": 1,
            },
        ]
    )
    distribution_change = pd.DataFrame(
        [
            {
                "condition": "treated",
                "cluster": "C0",
                "count": 1,
                "fraction": 1 / 4,
                "replicate_n": 1,
                "reference_fraction": 2 / 3,
                "delta_fraction": -5 / 12,
                "log2_fc": -1.415037499278844,
            },
            {
                "condition": "treated",
                "cluster": "C1",
                "count": 3,
                "fraction": 3 / 4,
                "replicate_n": 1,
                "reference_fraction": 1 / 3,
                "delta_fraction": 5 / 12,
                "log2_fc": 1.1699250014423124,
            },
        ]
    )
    return SharedClusterResult(
        model=model,
        assignments=pd.DataFrame(columns=["sample_id", "condition", "cluster"]),
        cluster_distribution=cluster_distribution,
        condition_distribution=condition_distribution,
        distribution_change=distribution_change,
        cluster_profiles=pd.DataFrame(columns=["cluster", "count", "f0", "f1"]),
        region_summaries=None,
        plot_data={},
        metadata={"mode": "region_anchored"},
    )


def test_prepare_shared_cluster_distribution_data_returns_distribution_payload():
    result = _make_shared_cluster_result()

    payload = plotting.prepare_shared_cluster_distribution_data(result=result)
    sample_distribution = payload["sample_distribution"]
    condition_distribution = payload["condition_distribution"]
    distribution_change = payload["distribution_change"]

    assert set(payload) == {
        "sample_distribution",
        "condition_distribution",
        "distribution_change",
        "metadata",
    }
    assert sample_distribution["sample_id"].tolist() == ["s1", "s1", "s2", "s2"]
    assert sample_distribution["condition"].tolist() == ["NS", "NS", "treated", "treated"]
    assert sample_distribution["cluster"].tolist() == ["C0", "C1", "C0", "C1"]
    assert sample_distribution["count"].tolist() == [2, 1, 1, 3]
    assert sample_distribution["fraction"].tolist() == pytest.approx(
        [2 / 3, 1 / 3, 1 / 4, 3 / 4]
    )
    assert condition_distribution["condition"].tolist() == ["NS", "NS", "treated", "treated"]
    assert condition_distribution["cluster"].tolist() == ["C0", "C1", "C0", "C1"]
    assert condition_distribution["count"].tolist() == [2, 1, 1, 3]
    assert condition_distribution["fraction"].tolist() == pytest.approx(
        [2 / 3, 1 / 3, 1 / 4, 3 / 4]
    )
    assert condition_distribution["replicate_n"].tolist() == [1, 1, 1, 1]
    assert distribution_change["condition"].tolist() == ["treated", "treated"]
    assert distribution_change["cluster"].tolist() == ["C0", "C1"]
    assert distribution_change["count"].tolist() == [1, 3]
    assert distribution_change["fraction"].tolist() == pytest.approx([1 / 4, 3 / 4])
    assert distribution_change["reference_fraction"].tolist() == pytest.approx([2 / 3, 1 / 3])
    assert distribution_change["delta_fraction"].tolist() == pytest.approx([-5 / 12, 5 / 12])
    assert distribution_change["log2_fc"].tolist() == pytest.approx(
        [-1.415037499278844, 1.1699250014423124]
    )
    assert payload["metadata"]["mode"] == "region_anchored"


def test_plotting_matplotlib_module_exports_save_figure():
    from dimelo import plotting_matplotlib

    assert hasattr(plotting_matplotlib, "save_figure")


def test_save_figure_writes_png(tmp_path):
    from matplotlib import pyplot as plt
    from dimelo import plotting_matplotlib

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])

    output_path = tmp_path / "figure.png"
    plotting_matplotlib.save_figure(fig, output_path)

    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_prepare_shared_cluster_distribution_data_handles_missing_change_table():
    result = _make_shared_cluster_result()
    result.distribution_change = None

    payload = plotting.prepare_shared_cluster_distribution_data(result=result)
    sample_distribution = payload["sample_distribution"]
    condition_distribution = payload["condition_distribution"]

    assert list(payload["distribution_change"].columns) == [
        "condition",
        "cluster",
        "count",
        "fraction",
        "replicate_n",
        "reference_fraction",
        "delta_fraction",
        "log2_fc",
    ]
    assert payload["distribution_change"].empty
    assert sample_distribution["sample_id"].tolist() == ["s1", "s1", "s2", "s2"]
    assert sample_distribution["cluster"].tolist() == ["C0", "C1", "C0", "C1"]
    assert sample_distribution["count"].tolist() == [2, 1, 1, 3]
    assert sample_distribution["fraction"].tolist() == pytest.approx(
        [2 / 3, 1 / 3, 1 / 4, 3 / 4]
    )
    assert condition_distribution["condition"].tolist() == ["NS", "NS", "treated", "treated"]
    assert condition_distribution["cluster"].tolist() == ["C0", "C1", "C0", "C1"]
    assert condition_distribution["count"].tolist() == [2, 1, 1, 3]
    assert condition_distribution["fraction"].tolist() == pytest.approx(
        [2 / 3, 1 / 3, 1 / 4, 3 / 4]
    )
    assert condition_distribution["replicate_n"].tolist() == [1, 1, 1, 1]


def test_plot_shared_cluster_distribution_matplotlib_returns_figure_and_axis():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_distribution_data(result=_make_shared_cluster_result())
    payload = dict(payload)
    payload["metadata"] = {**payload["metadata"], "cluster_labels": ["C1"]}
    payload["sample_distribution"] = pd.DataFrame(
        [
            {
                "sample_id": "s2",
                "condition": "treated",
                "cluster": "C3",
                "count": 1,
                "fraction": 0.25,
            },
            {
                "sample_id": "s1",
                "condition": "NS",
                "cluster": "C20",
                "count": 2,
                "fraction": 0.75,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "cluster": "C20",
                "count": 3,
                "fraction": 0.75,
            },
            {
                "sample_id": "s1",
                "condition": "NS",
                "cluster": "C3",
                "count": 1,
                "fraction": 0.25,
            },
        ]
    )
    payload["condition_distribution"] = pd.DataFrame(
        [
            {
                "condition": "treated",
                "cluster": "C3",
                "count": 1,
                "fraction": 0.25,
                "replicate_n": 1,
            },
            {
                "condition": "treated",
                "cluster": "C20",
                "count": 3,
                "fraction": 0.75,
                "replicate_n": 1,
            },
            {
                "condition": "NS",
                "cluster": "C3",
                "count": 1,
                "fraction": 0.25,
                "replicate_n": 1,
            },
            {
                "condition": "NS",
                "cluster": "C20",
                "count": 2,
                "fraction": 0.75,
                "replicate_n": 1,
            },
        ]
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_distribution_matplotlib(payload, level="sample")

    assert fig is not None
    assert ax is not None
    assert [text.get_text() for text in ax.get_legend().texts] == ["C3", "C20"]
    assert [tick.get_text() for tick in ax.get_xticklabels() if tick.get_text()] == ["s2", "s1"]

    fig, ax = plotting_matplotlib.plot_shared_cluster_distribution_matplotlib(payload, level="condition")

    assert fig is not None
    assert ax is not None
    assert [text.get_text() for text in ax.get_legend().texts] == ["C3", "C20"]
    assert [tick.get_text() for tick in ax.get_xticklabels() if tick.get_text()] == ["treated", "NS"]

    duplicate_payload = dict(payload)
    duplicate_payload["condition_distribution"] = pd.concat(
        [payload["condition_distribution"], payload["condition_distribution"].iloc[[0]]],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="duplicate cluster fractions"):
        plotting_matplotlib.plot_shared_cluster_distribution_matplotlib(duplicate_payload)


def test_plot_shared_cluster_change_matplotlib_returns_figure_and_axis():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_distribution_data(result=_make_shared_cluster_result())
    payload = dict(payload)
    payload["metadata"] = {**payload["metadata"], "cluster_labels": ["C1"]}
    payload["distribution_change"] = pd.DataFrame(
        [
            {
                "condition": "treated",
                "cluster": "C3",
                "count": 1,
                "fraction": 0.25,
                "replicate_n": 1,
                "reference_fraction": 0.5,
                "delta_fraction": -0.25,
                "log2_fc": -1.0,
            },
            {
                "condition": "treated",
                "cluster": "C20",
                "count": 3,
                "fraction": 0.75,
                "replicate_n": 1,
                "reference_fraction": 0.5,
                "delta_fraction": 0.25,
                "log2_fc": 1.0,
            },
            {
                "condition": "baseline",
                "cluster": "C3",
                "count": 2,
                "fraction": 0.5,
                "replicate_n": 1,
                "reference_fraction": 0.5,
                "delta_fraction": 0.0,
                "log2_fc": 0.0,
            },
            {
                "condition": "baseline",
                "cluster": "C20",
                "count": 2,
                "fraction": 0.5,
                "replicate_n": 1,
                "reference_fraction": 0.5,
                "delta_fraction": 0.0,
                "log2_fc": 0.0,
            },
        ]
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_change_matplotlib(payload)

    assert fig is not None
    assert ax is not None
    assert [tick.get_text() for tick in ax.get_xticklabels()] == ["C3", "C20"]
    assert [tick.get_text() for tick in ax.get_yticklabels()] == ["treated", "baseline"]

    duplicate_payload = dict(payload)
    duplicate_payload["distribution_change"] = pd.concat(
        [payload["distribution_change"], payload["distribution_change"].iloc[[0]]],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="duplicate delta fractions"):
        plotting_matplotlib.plot_shared_cluster_change_matplotlib(duplicate_payload)


def _make_shared_cluster_profile_result() -> SharedClusterResult:
    return SharedClusterResult(
        model=SharedClusterModel(
            mode="region_anchored",
            motifs=["A,0"],
            feature_names=["f0", "f1"],
            preprocessing={"signal_normalization": "none"},
            estimator=object(),
            cluster_labels=["C0", "C1"],
            fit_metadata={"clusterer": "minibatch_kmeans", "n_clusters": 2},
        ),
        assignments=pd.DataFrame(columns=["sample_id", "condition", "cluster"]),
        cluster_distribution=pd.DataFrame(
            columns=["sample_id", "condition", "cluster", "count", "fraction"]
        ),
        condition_distribution=pd.DataFrame(
            columns=["condition", "cluster", "count", "fraction", "replicate_n"]
        ),
        distribution_change=None,
        cluster_profiles=pd.DataFrame(
            [
                {"cluster": "C0", "count": 3, "f0": 0.1, "f1": 0.2},
                {"cluster": "C1", "count": 4, "f0": 0.8, "f1": 0.9},
            ]
        ),
        region_summaries=None,
        plot_data={},
        metadata={"mode": "region_anchored"},
    )


def test_prepare_shared_cluster_profile_data_returns_long_form_profiles():
    result = _make_shared_cluster_profile_result()

    payload = plotting.prepare_shared_cluster_profile_data(result=result)
    profile_table = payload["profile_table"]

    assert set(payload) == {"profile_table", "metadata"}
    assert set(profile_table.columns) == {"cluster", "feature", "value", "count"}
    assert profile_table.to_dict("records") == [
        {"cluster": "C0", "feature": "f0", "value": pytest.approx(0.1), "count": 3},
        {"cluster": "C0", "feature": "f1", "value": pytest.approx(0.2), "count": 3},
        {"cluster": "C1", "feature": "f0", "value": pytest.approx(0.8), "count": 4},
        {"cluster": "C1", "feature": "f1", "value": pytest.approx(0.9), "count": 4},
    ]
    assert payload["metadata"]["feature_names"] == ["f0", "f1"]


def test_prepare_shared_cluster_profile_data_respects_feature_subset():
    result = _make_shared_cluster_profile_result()

    payload = plotting.prepare_shared_cluster_profile_data(result=result, features=["f1"])
    profile_table = payload["profile_table"]

    assert profile_table.to_dict("records") == [
        {"cluster": "C0", "feature": "f1", "value": pytest.approx(0.2), "count": 3},
        {"cluster": "C1", "feature": "f1", "value": pytest.approx(0.9), "count": 4},
    ]


def _make_shared_cluster_region_result() -> SharedClusterResult:
    return SharedClusterResult(
        model=SharedClusterModel(
            mode="region_anchored",
            motifs=["A,0"],
            feature_names=["f0", "f1"],
            preprocessing={"signal_normalization": "none"},
            estimator=object(),
            cluster_labels=["C0", "C1"],
            fit_metadata={"clusterer": "minibatch_kmeans", "n_clusters": 2},
        ),
        assignments=pd.DataFrame(columns=["sample_id", "condition", "cluster"]),
        cluster_distribution=pd.DataFrame(
            columns=["sample_id", "condition", "cluster", "count", "fraction"]
        ),
        condition_distribution=pd.DataFrame(
            columns=["condition", "cluster", "count", "fraction", "replicate_n"]
        ),
        distribution_change=None,
        cluster_profiles=pd.DataFrame(columns=["cluster", "count", "f0", "f1"]),
        region_summaries=pd.DataFrame(
            [
                {
                    "region_id": "reg1",
                    "sample_id": "s1",
                    "condition": "NS",
                    "cluster": "C0",
                    "count": 2,
                    "fraction": 2 / 3,
                },
                {
                    "region_id": "reg1",
                    "sample_id": "s1",
                    "condition": "NS",
                    "cluster": "C1",
                    "count": 1,
                    "fraction": 1 / 3,
                },
                {
                    "region_id": "reg1",
                    "sample_id": "s2",
                    "condition": "treated",
                    "cluster": "C0",
                    "count": 1,
                    "fraction": 1 / 4,
                },
                {
                    "region_id": "reg1",
                    "sample_id": "s2",
                    "condition": "treated",
                    "cluster": "C1",
                    "count": 3,
                    "fraction": 3 / 4,
                },
            ]
        ),
        plot_data={},
        metadata={"mode": "region_anchored"},
    )


def test_prepare_shared_cluster_region_data_returns_sample_and_condition_tables():
    result = _make_shared_cluster_region_result()

    payload = plotting.prepare_shared_cluster_region_data(result=result)

    assert set(payload) == {"region_table", "condition_region_table", "metadata"}
    region_table = payload["region_table"]
    condition_region_table = payload["condition_region_table"]
    assert region_table.to_dict("records") == [
        {
            "region_id": "reg1",
            "sample_id": "s1",
            "condition": "NS",
            "cluster": "C0",
            "count": 2,
            "fraction": pytest.approx(2 / 3),
        },
        {
            "region_id": "reg1",
            "sample_id": "s1",
            "condition": "NS",
            "cluster": "C1",
            "count": 1,
            "fraction": pytest.approx(1 / 3),
        },
        {
            "region_id": "reg1",
            "sample_id": "s2",
            "condition": "treated",
            "cluster": "C0",
            "count": 1,
            "fraction": pytest.approx(1 / 4),
        },
        {
            "region_id": "reg1",
            "sample_id": "s2",
            "condition": "treated",
            "cluster": "C1",
            "count": 3,
            "fraction": pytest.approx(3 / 4),
        },
    ]
    assert condition_region_table.to_dict("records") == [
        {
            "region_id": "reg1",
            "condition": "NS",
            "cluster": "C0",
            "count": 2,
            "fraction_mean": pytest.approx(2 / 3),
            "fraction_median": pytest.approx(2 / 3),
            "sample_n": 1,
        },
        {
            "region_id": "reg1",
            "condition": "NS",
            "cluster": "C1",
            "count": 1,
            "fraction_mean": pytest.approx(1 / 3),
            "fraction_median": pytest.approx(1 / 3),
            "sample_n": 1,
        },
        {
            "region_id": "reg1",
            "condition": "treated",
            "cluster": "C0",
            "count": 1,
            "fraction_mean": pytest.approx(1 / 4),
            "fraction_median": pytest.approx(1 / 4),
            "sample_n": 1,
        },
        {
            "region_id": "reg1",
            "condition": "treated",
            "cluster": "C1",
            "count": 3,
            "fraction_mean": pytest.approx(3 / 4),
            "fraction_median": pytest.approx(3 / 4),
            "sample_n": 1,
        },
    ]
    assert payload["metadata"]["mode"] == "region_anchored"


def test_prepare_shared_cluster_region_data_can_disable_condition_aggregation():
    result = _make_shared_cluster_region_result()

    payload = plotting.prepare_shared_cluster_region_data(
        result=result,
        aggregate_conditions=False,
    )

    assert payload["condition_region_table"].empty


def test_prepare_shared_cluster_region_data_rejects_missing_region_summaries():
    result = _make_shared_cluster_region_result()
    result.region_summaries = None

    with pytest.raises(ValueError, match="region_summaries"):
        plotting.prepare_shared_cluster_region_data(result=result)


def _minimal_region_contrast_result() -> RegionContrastResult:
    regions = pd.DataFrame(
        [
            {"region_id": "chr1:90-110,+", "condition": "NS", "fraction": 0.20, "rank": 2},
            {"region_id": "chr1:90-110,+", "condition": "15min", "fraction": 0.55, "rank": 2},
            {"region_id": "chr1:190-210,-", "condition": "NS", "fraction": 0.30, "rank": 1},
            {"region_id": "chr1:190-210,-", "condition": "15min", "fraction": 0.70, "rank": 1},
        ]
    )
    summary = pd.DataFrame(
        [
            {
                "region_id": "chr1:90-110,+",
                "fraction": 0.55,
                "reference_fraction": 0.20,
                "delta_fraction": 0.35,
                "rank": 2,
            },
            {
                "region_id": "chr1:190-210,-",
                "fraction": 0.70,
                "reference_fraction": 0.30,
                "delta_fraction": 0.40,
                "rank": 1,
            },
        ]
    )
    return RegionContrastResult(
        regions=regions,
        summary=summary,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["15min"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        metadata={
            "analysis_unit": "ensemble_region",
            "representation": "modified_fraction",
            "signal_source": "pileup_counts",
            "test": "effect_size_only",
        },
        plot_data={},
    )


def _group_vs_group_region_contrast_result() -> RegionContrastResult:
    summary = pd.DataFrame(
        [
            {
                "region_id": "chr1:90-110,+",
                "fraction": 0.60,
                "reference_fraction": 0.10,
                "delta_fraction": 0.50,
                "rank": 1,
            }
        ]
    )
    return RegionContrastResult(
        regions=pd.DataFrame(
            [
                {"region_id": "chr1:90-110,+", "condition": "NS", "fraction": 0.10, "rank": 1},
                {"region_id": "chr1:90-110,+", "condition": "15min", "fraction": 0.50, "rank": 1},
                {"region_id": "chr1:90-110,+", "condition": "30min", "fraction": 0.70, "rank": 1},
            ]
        ),
        summary=summary,
        contrast=ContrastSpec(
            mode="group_vs_group",
            numerator=["15min", "30min"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        metadata={
            "analysis_unit": "ensemble_region",
            "representation": "modified_fraction",
            "signal_source": "pileup_counts",
            "test": "effect_size_only",
        },
        plot_data={},
    )


def _region_contrast_plot_setup(position_rows: list[dict[str, object]]) -> tuple[
    RegionContrastResult,
    pd.DataFrame,
    plotting.AxisSpec,
    plotting.AggregationSpec,
]:
    return (
        _minimal_region_contrast_result(),
        pd.DataFrame(position_rows),
        plotting.AxisSpec(
            orientation="region_5to3",
            coordinate_mode="fixed_window",
            anchor="center",
            upstream_bp=20,
            downstream_bp=20,
        ),
        plotting.AggregationSpec(),
    )


def _region_contrast_position_rows(*, include_grouping_key: bool = True) -> list[dict[str, object]]:
    base_rows = [
        {
            "region_id": "chr1:90-110,+",
            "condition": "NS",
            "position": 95,
            "anchor": 100,
            "value": 0.1,
            "region_strand": "+",
        },
        {
            "region_id": "chr1:90-110,+",
            "condition": "15min",
            "position": 95,
            "anchor": 100,
            "value": 0.6,
            "region_strand": "+",
        },
        {
            "region_id": "chr1:90-110,+",
            "condition": "NS",
            "position": 105,
            "anchor": 100,
            "value": 0.2,
            "region_strand": "+",
        },
        {
            "region_id": "chr1:90-110,+",
            "condition": "15min",
            "position": 105,
            "anchor": 100,
            "value": 0.7,
            "region_strand": "+",
        },
        {
            "region_id": "chr1:190-210,-",
            "condition": "NS",
            "position": 195,
            "anchor": 200,
            "value": 0.2,
            "region_strand": "-",
        },
        {
            "region_id": "chr1:190-210,-",
            "condition": "15min",
            "position": 195,
            "anchor": 200,
            "value": 0.8,
            "region_strand": "-",
        },
        {
            "region_id": "chr1:190-210,-",
            "condition": "NS",
            "position": 205,
            "anchor": 200,
            "value": 0.3,
            "region_strand": "-",
        },
        {
            "region_id": "chr1:190-210,-",
            "condition": "15min",
            "position": 205,
            "anchor": 200,
            "value": 0.9,
            "region_strand": "-",
        },
    ]
    if include_grouping_key:
        return base_rows
    return [
        {key: value for key, value in row.items() if key != "condition"}
        for row in base_rows[:1]
    ]


def _make_region_discovery_result() -> RegionDiscoveryResult:
    windows = pd.DataFrame(
        [
            {
                "window_id": "chr1:0-100:+",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": "+",
                "score_value": 0.2,
                "rank": 3,
            },
            {
                "window_id": "chr1:100-200:+",
                "chromosome": "chr1",
                "start": 100,
                "end": 200,
                "strand": "+",
                "score_value": 0.9,
                "rank": 1,
            },
            {
                "window_id": "chr2:0-100:+",
                "chromosome": "chr2",
                "start": 0,
                "end": 100,
                "strand": "+",
                "score_value": 0.5,
                "rank": 2,
            },
        ]
    )
    hits = windows.loc[windows["rank"] <= 2].copy()
    return RegionDiscoveryResult(
        windows=windows,
        hits=hits,
        contrast=None,
        plot_data={
            "window_score_table": windows.copy(),
            "top_hits_table": hits.copy(),
        },
        metadata={
            "score": "effect_size_only",
            "contrast_mode": "pairwise",
            "merge_hits": False,
        },
    )


def _make_global_analysis_result() -> GlobalAnalysisResult:
    summary = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "modified_count": 10,
                "valid_count": 20,
                "global_fraction": 0.5,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "modified_count": 16,
                "valid_count": 20,
                "global_fraction": 0.8,
            },
            {
                "sample_id": "s3",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "modified_count": 12,
                "valid_count": 20,
                "global_fraction": 0.6,
            },
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "CG,0",
                "modified_count": 4,
                "valid_count": 20,
                "global_fraction": 0.2,
            },
        ]
    )
    normalization_factors = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "global_fraction": 0.5,
                "reference_fraction": 0.6333333333,
                "global_offset": -0.1333333333,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "global_fraction": 0.8,
                "reference_fraction": 0.6333333333,
                "global_offset": 0.1666666667,
            },
            {
                "sample_id": "s3",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "global_fraction": 0.6,
                "reference_fraction": 0.6333333333,
                "global_offset": -0.0333333333,
            },
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "CG,0",
                "global_fraction": 0.2,
                "reference_fraction": 0.2,
                "global_offset": 0.0,
            },
        ]
    )
    windows = pd.DataFrame(
        columns=[
            "sample_id",
            "condition",
            "replicate",
            "motif",
            "window_id",
            "chromosome",
            "start",
            "end",
            "strand",
            "modified_count",
            "valid_count",
            "window_fraction",
        ]
    )
    return GlobalAnalysisResult(
        summary=summary,
        windows=windows,
        normalization_factors=normalization_factors,
        plot_data={
            "global_fraction_bar": summary.copy(),
            "window_fraction_table": windows.copy(),
        },
        metadata={"window_size": 100000, "step_size": 50000},
    )


def _make_global_window_result() -> GlobalAnalysisResult:
    result = _make_global_analysis_result()
    result.windows = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-100",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 5,
                "valid_count": 10,
                "window_fraction": 0.5,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-100",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 8,
                "valid_count": 10,
                "window_fraction": 0.8,
            },
            {
                "sample_id": "s3",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "window_id": "chr1:0-100",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 10,
                "window_fraction": 0.4,
            },
            {
                "sample_id": "s3",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "window_id": "chr2:0-100",
                "chromosome": "chr2",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 6,
                "valid_count": 10,
                "window_fraction": 0.6,
            },
        ]
    )
    result.plot_data["window_fraction_table"] = result.windows.copy()
    return result


def test_prepare_region_discovery_scan_data_returns_expected_tables():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_scan_data(result=result)

    assert set(payload) == {"scan_table", "hit_table", "metadata"}
    assert list(payload["scan_table"]["contig"]) == ["chr1", "chr1", "chr2"]
    assert list(payload["hit_table"]["rank"]) == [1, 2]
    assert payload["scan_table"]["is_hit"].tolist() == [False, True, True]
    assert payload["metadata"]["contig_order"] == ["chr1", "chr2"]
    assert payload["metadata"]["score_column"] == "score_value"


def test_prepare_region_discovery_scan_data_filters_contigs_in_requested_order():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_scan_data(
        result=result,
        contigs=["chr2"],
    )

    assert payload["metadata"]["contig_order"] == ["chr2"]
    assert payload["scan_table"]["contig"].tolist() == ["chr2"]
    assert payload["hit_table"]["contig"].tolist() == ["chr2"]
    assert "chr1" not in payload["scan_table"]["contig"].tolist()
    assert "chr1" not in payload["hit_table"]["contig"].tolist()


def test_prepare_region_discovery_scan_data_limits_hit_overlay():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_scan_data(
        result=result,
        top_n_hits=1,
    )

    assert payload["hit_table"]["rank"].tolist() == [1]
    assert payload["scan_table"]["is_hit"].tolist() == [False, True, False]


def test_prepare_region_discovery_scan_data_rejects_negative_top_n_hits():
    result = _make_region_discovery_result()

    with pytest.raises(ValueError, match="top_n_hits must be non-negative"):
        plotting.prepare_region_discovery_scan_data(
            result=result,
            top_n_hits=-1,
        )


def test_prepare_region_discovery_scan_data_rejects_hits_missing_from_filtered_windows():
    result = _make_region_discovery_result()
    inconsistent_hit = pd.DataFrame(
        [
            {
                "window_id": "chr2:100-200:+",
                "chromosome": "chr2",
                "start": 100,
                "end": 200,
                "strand": "+",
                "score_value": 0.8,
                "rank": 1,
            }
        ]
    )
    inconsistent_result = RegionDiscoveryResult(
        windows=result.windows.copy(),
        hits=pd.concat([result.hits, inconsistent_hit], ignore_index=True),
        contrast=result.contrast,
        plot_data=result.plot_data,
        metadata=result.metadata,
        figures=result.figures,
    )

    with pytest.raises(ValueError, match="result.hits contains window_id values not present in the filtered windows table"):
        plotting.prepare_region_discovery_scan_data(
            result=inconsistent_result,
            contigs=["chr2"],
        )


def test_prepare_region_discovery_hit_context_data_uses_top_ranked_hits():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_hit_context_data(
        result=result,
        top_n=1,
        padding_windows=1,
    )

    assert set(payload) == {"context_table", "selected_hits", "metadata"}
    assert payload["selected_hits"]["window_id"].tolist() == ["chr1:100-200:+"]
    assert payload["selected_hits"]["rank"].tolist() == [1]
    assert payload["context_table"]["selected_hit_id"].tolist() == ["chr1:100-200:+", "chr1:100-200:+"]
    assert payload["context_table"]["window_id"].tolist() == ["chr1:0-100:+", "chr1:100-200:+"]
    assert payload["context_table"]["selected_hit_rank"].nunique() == 1
    assert payload["context_table"]["selected_hit_rank"].iloc[0] == 1
    assert payload["context_table"]["is_selected_hit"].sum() == 1


def test_prepare_region_discovery_hit_context_data_adds_relative_window_offsets():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_hit_context_data(
        result=result,
        top_n=1,
        padding_windows=1,
    )

    assert payload["context_table"]["window_id"].tolist() == ["chr1:0-100:+", "chr1:100-200:+"]
    assert payload["context_table"]["relative_window_offset"].tolist() == [-1, 0]


def test_prepare_region_discovery_hit_context_data_returns_empty_payload_for_no_hits():
    base_result = _make_region_discovery_result()
    result = RegionDiscoveryResult(
        windows=base_result.windows.copy(),
        hits=pd.DataFrame(columns=base_result.hits.columns),
        contrast=None,
        plot_data={},
        metadata={"score": "effect_size_only", "contrast_mode": "pairwise", "merge_hits": False},
    )

    payload = plotting.prepare_region_discovery_hit_context_data(result=result)

    assert payload["context_table"].empty
    assert payload["selected_hits"].empty
    assert payload["metadata"]["selection_mode"] == "top_n"


def test_plot_region_discovery_scan_matplotlib_returns_figure_and_axes():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_region_discovery_scan_data(result=_make_region_discovery_result())

    fig, axes = plotting_matplotlib.plot_region_discovery_scan_matplotlib(payload)

    assert fig is not None
    assert axes is not None


def test_plot_region_discovery_hit_context_matplotlib_returns_figure_and_axes():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_region_discovery_hit_context_data(
        result=_make_region_discovery_result(),
        top_n=2,
    )

    fig, axes = plotting_matplotlib.plot_region_discovery_hit_context_matplotlib(payload)

    assert fig is not None
    assert axes is not None


def test_prepare_global_analysis_summary_data_returns_expected_tables():
    result = _make_global_analysis_result()

    payload = plotting.prepare_global_analysis_summary_data(result=result)

    assert set(payload) == {"sample_summary", "condition_summary", "normalization_table", "metadata"}
    assert payload["sample_summary"]["sample_id"].tolist() == ["s1", "s2", "s3", "s1"]
    assert payload["normalization_table"]["sample_id"].tolist() == ["s1", "s2", "s3", "s1"]
    assert (
        payload["condition_summary"][["condition", "motif"]].apply(tuple, axis=1).tolist()
        == [("NS", "A,0"), ("treated", "A,0"), ("NS", "CG,0")]
    )
    assert payload["condition_summary"]["sample_n"].tolist() == [1, 2, 1]
    assert payload["metadata"]["motifs"] == ["A,0", "CG,0"]


def test_prepare_global_analysis_summary_data_computes_condition_means():
    result = _make_global_analysis_result()

    payload = plotting.prepare_global_analysis_summary_data(result=result)
    condition_summary = payload["condition_summary"].set_index(["condition", "motif"])

    assert condition_summary.loc[("NS", "A,0"), "global_fraction_mean"] == pytest.approx(0.5)
    assert condition_summary.loc[("treated", "A,0"), "global_fraction_mean"] == pytest.approx(0.7)
    assert condition_summary.loc[("treated", "A,0"), "global_fraction_median"] == pytest.approx(0.7)


def test_prepare_global_analysis_summary_data_filters_motifs():
    result = _make_global_analysis_result()

    payload = plotting.prepare_global_analysis_summary_data(
        result=result,
        motifs=["A,0"],
    )

    assert payload["sample_summary"]["motif"].unique().tolist() == ["A,0"]
    assert payload["normalization_table"]["motif"].unique().tolist() == ["A,0"]
    assert "CG,0" not in payload["sample_summary"]["motif"].tolist()
    assert "CG,0" not in payload["normalization_table"]["motif"].tolist()


def test_prepare_global_analysis_summary_data_skips_condition_aggregation_when_disabled():
    result = _make_global_analysis_result()

    payload = plotting.prepare_global_analysis_summary_data(
        result=result,
        aggregate_conditions=False,
    )

    assert payload["condition_summary"].empty
    assert payload["metadata"]["aggregate_conditions"] is False


def test_plot_global_analysis_summary_matplotlib_returns_figure_and_axis():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_global_analysis_summary_data(result=_make_global_analysis_result())

    fig, ax = plotting_matplotlib.plot_global_analysis_summary_matplotlib(payload)

    assert fig is not None
    assert ax is not None


def test_prepare_global_analysis_window_data_returns_expected_tables():
    result = _make_global_window_result()

    payload = plotting.prepare_global_analysis_window_data(result=result)

    assert set(payload) == {"window_table", "condition_window_table", "metadata"}
    assert payload["window_table"]["contig"].tolist() == ["chr1", "chr1", "chr1", "chr2"]
    assert payload["window_table"]["window_midpoint"].tolist() == [50.0, 50.0, 50.0, 50.0]
    assert payload["condition_window_table"].empty
    assert payload["metadata"]["contig_order"] == ["chr1", "chr2"]


def test_plot_global_analysis_window_matplotlib_returns_figure_and_axes():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_global_analysis_window_data(result=_make_global_window_result())

    fig, axes = plotting_matplotlib.plot_global_analysis_window_matplotlib(payload)

    assert fig is not None
    assert axes is not None


def test_prepare_global_analysis_window_data_filters_contigs_in_requested_order():
    result = _make_global_window_result()

    payload = plotting.prepare_global_analysis_window_data(
        result=result,
        contigs=["chr2"],
    )

    assert payload["metadata"]["contig_order"] == ["chr2"]
    assert payload["window_table"]["contig"].tolist() == ["chr2"]


def test_prepare_global_analysis_window_data_aggregates_conditions():
    result = _make_global_window_result()

    payload = plotting.prepare_global_analysis_window_data(
        result=result,
        aggregate_conditions=True,
    )

    condition_rows = payload["condition_window_table"]
    ns_row = condition_rows.loc[condition_rows["condition"] == "NS"].iloc[0]
    treated_chr1 = condition_rows.loc[
        (condition_rows["condition"] == "treated")
        & (condition_rows["contig"] == "chr1")
    ].iloc[0]

    assert ns_row["window_fraction_mean"] == pytest.approx(0.5)
    assert treated_chr1["window_fraction_mean"] == pytest.approx(0.6)
    assert treated_chr1["sample_n"] == 2


def test_prepare_global_analysis_window_data_aggregates_per_sample_windows_across_replicates_with_chrom_alias():
    result = _make_global_window_result()
    duplicate_sample_window = pd.DataFrame(
        [
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "window_id": "chr1:0-100",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 2,
                "valid_count": 10,
                "window_fraction": 0.2,
            }
        ]
    )
    result.windows = pd.concat([result.windows, duplicate_sample_window], ignore_index=True).rename(
        columns={"chromosome": "chrom"}
    )
    result.plot_data["window_fraction_table"] = result.windows.copy()

    payload = plotting.prepare_global_analysis_window_data(
        result=result,
        aggregate_conditions=True,
    )

    treated_chr1 = payload["condition_window_table"].loc[
        (payload["condition_window_table"]["condition"] == "treated")
        & (payload["condition_window_table"]["contig"] == "chr1")
    ].iloc[0]

    assert treated_chr1["window_fraction_mean"] == pytest.approx(0.45)
    assert treated_chr1["window_fraction_median"] == pytest.approx(0.45)
    assert treated_chr1["sample_n"] == 2


def test_prepare_region_discovery_hit_context_data_rejects_padding_bp():
    result = _make_region_discovery_result()

    with pytest.raises(ValueError, match="padding_bp is not supported"):
        plotting.prepare_region_discovery_hit_context_data(
            result=result,
            padding_bp=100,
        )


def test_prepare_region_discovery_hit_context_data_rejects_selected_hits_missing_from_windows():
    result = _make_region_discovery_result()
    inconsistent_hit = pd.DataFrame(
        [
            {
                "window_id": "chr1:200-300:+",
                "chromosome": "chr1",
                "start": 200,
                "end": 300,
                "strand": "+",
                "score_value": 0.7,
                "rank": 0,
            }
        ]
    )
    inconsistent_result = RegionDiscoveryResult(
        windows=result.windows.copy(),
        hits=pd.concat([result.hits, inconsistent_hit], ignore_index=True),
        contrast=result.contrast,
        plot_data=result.plot_data,
        metadata=result.metadata,
        figures=result.figures,
    )

    with pytest.raises(
        ValueError,
        match="selected hits contain window_id values not present in result.windows",
    ):
        plotting.prepare_region_discovery_hit_context_data(
            result=inconsistent_result,
            top_n=1,
        )


def test_prepare_region_contrast_profile_data_returns_all_value_modes():
    result, position_table, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows()
    )

    payload = plotting.prepare_region_contrast_profile_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )

    plot_table = payload["plot_table"]
    assert len(plot_table) == 12
    assert set(payload["plot_table"]["value_mode"]) == {"numerator", "denominator", "delta"}
    assert plot_table.groupby(["region_id", "position", "value_mode"]).size().to_dict() == {
        ("chr1:90-110,+", 95, "numerator"): 1,
        ("chr1:90-110,+", 95, "denominator"): 1,
        ("chr1:90-110,+", 95, "delta"): 1,
        ("chr1:90-110,+", 105, "numerator"): 1,
        ("chr1:90-110,+", 105, "denominator"): 1,
        ("chr1:90-110,+", 105, "delta"): 1,
        ("chr1:190-210,-", 195, "numerator"): 1,
        ("chr1:190-210,-", 195, "denominator"): 1,
        ("chr1:190-210,-", 195, "delta"): 1,
        ("chr1:190-210,-", 205, "numerator"): 1,
        ("chr1:190-210,-", 205, "denominator"): 1,
        ("chr1:190-210,-", 205, "delta"): 1,
    }
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:90-110,+")
        & (plot_table["position"] == 95)
        & (plot_table["value_mode"] == "numerator"),
        "value",
    ].iloc[0] == pytest.approx(0.6)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:90-110,+")
        & (plot_table["position"] == 95)
        & (plot_table["value_mode"] == "denominator"),
        "value",
    ].iloc[0] == pytest.approx(0.1)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:90-110,+")
        & (plot_table["position"] == 95)
        & (plot_table["value_mode"] == "delta"),
        "value",
    ].iloc[0] == pytest.approx(0.5)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:90-110,+")
        & (plot_table["position"] == 105)
        & (plot_table["value_mode"] == "numerator"),
        "value",
    ].iloc[0] == pytest.approx(0.7)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:90-110,+")
        & (plot_table["position"] == 105)
        & (plot_table["value_mode"] == "denominator"),
        "value",
    ].iloc[0] == pytest.approx(0.2)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:90-110,+")
        & (plot_table["position"] == 105)
        & (plot_table["value_mode"] == "delta"),
        "value",
    ].iloc[0] == pytest.approx(0.5)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:190-210,-")
        & (plot_table["position"] == 195)
        & (plot_table["value_mode"] == "numerator"),
        "value",
    ].iloc[0] == pytest.approx(0.8)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:190-210,-")
        & (plot_table["position"] == 195)
        & (plot_table["value_mode"] == "denominator"),
        "value",
    ].iloc[0] == pytest.approx(0.2)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:190-210,-")
        & (plot_table["position"] == 195)
        & (plot_table["value_mode"] == "delta"),
        "value",
    ].iloc[0] == pytest.approx(0.6)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:190-210,-")
        & (plot_table["position"] == 205)
        & (plot_table["value_mode"] == "numerator"),
        "value",
    ].iloc[0] == pytest.approx(0.9)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:190-210,-")
        & (plot_table["position"] == 205)
        & (plot_table["value_mode"] == "denominator"),
        "value",
    ].iloc[0] == pytest.approx(0.3)
    assert plot_table.loc[
        (plot_table["region_id"] == "chr1:190-210,-")
        & (plot_table["position"] == 205)
        & (plot_table["value_mode"] == "delta"),
        "value",
    ].iloc[0] == pytest.approx(0.6)
    assert payload["metadata"]["plot_family"] == "region_contrast_profile"


def test_plot_region_contrast_profile_matplotlib_defaults_to_delta_and_honors_ax_title():
    from dimelo import plotting_matplotlib
    from matplotlib import pyplot as plt

    result, position_table, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows()
    )

    payload = plotting.prepare_region_contrast_profile_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )

    x_column = "plot_x" if "plot_x" in payload["plot_table"].columns else "position"
    expected = (
        payload["plot_table"]
        .loc[lambda table: table["value_mode"] == "delta", [x_column, "value"]]
        .groupby(x_column, as_index=False, sort=True)
        .mean(numeric_only=True)
    )

    fig, provided_ax = plt.subplots()
    fig, ax = plotting_matplotlib.plot_region_contrast_profile_matplotlib(
        payload,
        ax=provided_ax,
        title="Custom region contrast profile",
    )

    assert fig is provided_ax.figure
    assert ax is provided_ax
    assert ax.get_title() == "Custom region contrast profile"
    assert len(ax.lines) == 1
    np.testing.assert_allclose(np.asarray(ax.lines[0].get_xdata()), expected[x_column].to_numpy())
    np.testing.assert_allclose(np.asarray(ax.lines[0].get_ydata()), expected["value"].to_numpy())
    plt.close(fig)


def test_prepare_region_contrast_profile_data_collapses_same_coordinate_labels_within_each_side():
    result = _group_vs_group_region_contrast_result()
    position_table = pd.DataFrame(
        [
            {
                "region_id": "chr1:90-110,+",
                "condition": "NS",
                "position": 95,
                "anchor": 100,
                "value": 0.1,
                "region_strand": "+",
            },
            {
                "region_id": "chr1:90-110,+",
                "condition": "15min",
                "position": 95,
                "anchor": 100,
                "value": 0.5,
                "region_strand": "+",
            },
            {
                "region_id": "chr1:90-110,+",
                "condition": "30min",
                "position": 95,
                "anchor": 100,
                "value": 0.7,
                "region_strand": "+",
            },
            {
                "region_id": "chr1:90-110,+",
                "condition": "NS",
                "position": 105,
                "anchor": 100,
                "value": 0.2,
                "region_strand": "+",
            },
            {
                "region_id": "chr1:90-110,+",
                "condition": "15min",
                "position": 105,
                "anchor": 100,
                "value": 0.6,
                "region_strand": "+",
            },
            {
                "region_id": "chr1:90-110,+",
                "condition": "30min",
                "position": 105,
                "anchor": 100,
                "value": 0.8,
                "region_strand": "+",
            },
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=20,
        downstream_bp=20,
    )
    aggregation = plotting.AggregationSpec()

    payload = plotting.prepare_region_contrast_profile_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )

    plot_table = payload["plot_table"]
    assert len(plot_table) == 6
    assert plot_table.groupby(["position", "value_mode"]).size().to_dict() == {
        (95, "numerator"): 1,
        (95, "denominator"): 1,
        (95, "delta"): 1,
        (105, "numerator"): 1,
        (105, "denominator"): 1,
        (105, "delta"): 1,
    }
    assert plot_table.loc[
        (plot_table["position"] == 95) & (plot_table["value_mode"] == "numerator"),
        "value",
    ].iloc[0] == pytest.approx(0.6)
    assert plot_table.loc[
        (plot_table["position"] == 95) & (plot_table["value_mode"] == "denominator"),
        "value",
    ].iloc[0] == pytest.approx(0.1)
    assert plot_table.loc[
        (plot_table["position"] == 95) & (plot_table["value_mode"] == "delta"),
        "value",
    ].iloc[0] == pytest.approx(0.5)
    assert plot_table.loc[
        (plot_table["position"] == 105) & (plot_table["value_mode"] == "numerator"),
        "value",
    ].iloc[0] == pytest.approx(0.7)
    assert plot_table.loc[
        (plot_table["position"] == 105) & (plot_table["value_mode"] == "denominator"),
        "value",
    ].iloc[0] == pytest.approx(0.2)
    assert plot_table.loc[
        (plot_table["position"] == 105) & (plot_table["value_mode"] == "delta"),
        "value",
    ].iloc[0] == pytest.approx(0.5)


def test_prepare_region_contrast_heatmap_data_orders_rows_by_rank():
    result, position_table, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows()
    )
    result.summary = pd.concat(
        [
            result.summary,
            pd.DataFrame(
                [
                    {
                        "region_id": "chr1:290-310,+",
                        "fraction": 0.15,
                        "reference_fraction": 0.05,
                        "delta_fraction": 0.10,
                        "rank": 3,
                    },
                    {
                        "region_id": "chr1:290-310,+",
                        "fraction": 0.16,
                        "reference_fraction": 0.05,
                        "delta_fraction": 0.11,
                        "rank": 8,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

    payload = plotting.prepare_region_contrast_heatmap_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )

    rank_rows = payload["plot_table"].loc[:, ["region_id", "row_order"]].drop_duplicates()
    assert list(rank_rows.sort_values("row_order")["region_id"]) == [
        "chr1:190-210,-",
        "chr1:90-110,+",
    ]
    assert list(payload["summary_table"]["region_id"]) == [
        "chr1:190-210,-",
        "chr1:90-110,+",
    ]
    assert payload["metadata"]["plot_family"] == "region_contrast_heatmap"


def test_plot_region_contrast_heatmap_matplotlib_defaults_to_delta_and_honors_ax_title():
    from dimelo import plotting_matplotlib
    from matplotlib import pyplot as plt

    result, position_table, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows()
    )

    payload = plotting.prepare_region_contrast_heatmap_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )

    x_column = "plot_x" if "plot_x" in payload["plot_table"].columns else "position"
    expected = (
        payload["plot_table"]
        .loc[lambda table: table["value_mode"] == "delta", [x_column, "row_order", "value"]]
        .pivot_table(index="row_order", columns=x_column, values="value", aggfunc="mean")
        .sort_index(axis=0)
        .sort_index(axis=1)
        .to_numpy()
    )

    fig, provided_ax = plt.subplots()
    fig, ax = plotting_matplotlib.plot_region_contrast_heatmap_matplotlib(
        payload,
        ax=provided_ax,
        title="Custom region contrast heatmap",
    )

    assert fig is provided_ax.figure
    assert ax is provided_ax
    assert ax.get_title() == "Custom region contrast heatmap"
    assert len(ax.images) == 1
    np.testing.assert_allclose(np.asarray(ax.images[0].get_array()), expected)
    plt.close(fig)


def test_plot_region_contrast_heatmap_matplotlib_accepts_minimal_payload_without_summary_table():
    from dimelo import plotting_matplotlib
    from matplotlib import pyplot as plt

    result, position_table, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows()
    )

    payload = plotting.prepare_region_contrast_heatmap_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )
    payload = {
        "plot_table": payload["plot_table"].copy(),
        "metadata": dict(payload["metadata"]),
    }

    fig, ax = plotting_matplotlib.plot_region_contrast_heatmap_matplotlib(payload)

    assert fig is not None
    assert ax is not None
    assert len(ax.images) == 1
    assert [tick.get_text() for tick in ax.get_yticklabels()] == [
        "chr1:190-210,-",
        "chr1:90-110,+",
    ]
    plt.close(fig)


def test_prepare_region_contrast_heatmap_data_rejects_conflicting_ranks_for_same_region():
    result, position_table, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows()
    )
    result.summary = pd.concat(
        [
            result.summary,
            pd.DataFrame(
                [
                    {
                        "region_id": "chr1:90-110,+",
                        "fraction": 0.55,
                        "reference_fraction": 0.20,
                        "delta_fraction": 0.35,
                        "rank": 7,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

    with pytest.raises(ValueError, match="exactly one rank value per plotted region"):
        plotting.prepare_region_contrast_heatmap_data(
            result=result,
            position_table=position_table,
            axis=axis,
            aggregation=aggregation,
            value_mode="all",
        )


def test_prepare_region_contrast_heatmap_data_requires_both_contrast_sides():
    result, _, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows()
    )
    position_table = pd.DataFrame(
        [
            {
                "region_id": "chr1:90-110,+",
                "condition": "NS",
                "position": 95,
                "anchor": 100,
                "value": 0.1,
                "region_strand": "+",
            }
        ]
    )

    with pytest.raises(ValueError, match="both contrast sides"):
        plotting.prepare_region_contrast_heatmap_data(
            result=result,
            position_table=position_table,
            axis=axis,
            aggregation=aggregation,
            value_mode="all",
        )


def test_prepare_region_contrast_profile_data_requires_grouping_key():
    result, position_table, axis, aggregation = _region_contrast_plot_setup(
        _region_contrast_position_rows(include_grouping_key=False)
    )

    with pytest.raises(ValueError, match="sample_id or condition"):
        plotting.prepare_region_contrast_profile_data(
            result=result,
            position_table=position_table,
            axis=axis,
            aggregation=aggregation,
            value_mode="all",
        )


@pytest.mark.parametrize(
    "position_rows",
    [
        [
            {
                "region_id": "chr1:90-110,+",
                "sample_id": "NS",
                "condition": "NS",
                "position": 95,
                "anchor": 100,
                "value": 0.1,
                "region_strand": "+",
            }
        ],
        [
            {
                "region_id": "chr1:90-110,+",
                "sample_id": "sample-a",
                "condition": "condition-a",
                "position": 95,
                "anchor": 100,
                "value": 0.1,
                "region_strand": "+",
            }
        ],
    ],
)
def test_prepare_region_contrast_profile_data_rejects_ambiguous_grouping_key(position_rows):
    result, _, _, _ = _region_contrast_plot_setup(_region_contrast_position_rows())
    position_table = pd.DataFrame(position_rows)
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=20,
        downstream_bp=20,
    )
    aggregation = plotting.AggregationSpec()

    with pytest.raises(ValueError, match="could not resolve a unique grouping key"):
        plotting.prepare_region_contrast_profile_data(
            result=result,
            position_table=position_table,
            axis=axis,
            aggregation=aggregation,
            value_mode="all",
        )


def test_prepare_single_read_plot_data_flips_negative_regions_to_5to3():
    reads = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "region_strand": "-",
                "event_pos": 110,
                "anchor": 100,
                "read_id": "r1",
            }
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="custom",
        upstream_bp=20,
        downstream_bp=20,
    )

    payload = plotting.prepare_single_read_plot_data(
        reads,
        plot_family="single_read_raster",
        axis=axis,
        position_column="event_pos",
        anchor_column="anchor",
        region_strand_column="region_strand",
    )

    assert payload["plot_table"].loc[0, "plot_x"] == -10


def test_prepare_single_read_plot_data_genomic_fixed_window_does_not_require_strand_column():
    reads = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "event_pos": 110,
                "anchor": 100,
                "read_id": "r1",
            }
        ]
    )
    axis = plotting.AxisSpec(
        orientation="genomic",
        coordinate_mode="fixed_window",
        anchor="custom",
        upstream_bp=20,
        downstream_bp=20,
    )

    payload = plotting.prepare_single_read_plot_data(
        reads,
        plot_family="single_read_raster",
        axis=axis,
        position_column="event_pos",
        anchor_column="anchor",
    )

    assert list(payload["plot_table"]["plot_x"]) == [10.0]


def test_prepare_single_read_plot_data_filters_rows_outside_fixed_window_bounds():
    reads = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "region_strand": "+",
                "event_pos": 89,
                "anchor": 100,
                "read_id": "r0",
            },
            {
                "region_id": "reg1",
                "region_strand": "+",
                "event_pos": 95,
                "anchor": 100,
                "read_id": "r1",
            },
            {
                "region_id": "reg1",
                "region_strand": "+",
                "event_pos": 105,
                "anchor": 100,
                "read_id": "r2",
            },
            {
                "region_id": "reg1",
                "region_strand": "+",
                "event_pos": 111,
                "anchor": 100,
                "read_id": "r3",
            },
        ]
    )
    axis = plotting.AxisSpec(
        orientation="genomic",
        coordinate_mode="fixed_window",
        anchor="custom",
        upstream_bp=10,
        downstream_bp=10,
    )

    payload = plotting.prepare_single_read_plot_data(
        reads,
        plot_family="single_read_raster",
        axis=axis,
        position_column="event_pos",
        anchor_column="anchor",
    )

    assert list(payload["plot_table"]["read_id"]) == ["r1", "r2"]
    assert list(payload["plot_table"]["plot_x"]) == [-5.0, 5.0]


@pytest.mark.parametrize("upstream_bp, downstream_bp", [(-1, 10), (10, -1)])
def test_validate_axis_spec_rejects_negative_fixed_window_bounds(upstream_bp, downstream_bp):
    axis = plotting.AxisSpec(
        orientation="genomic",
        coordinate_mode="fixed_window",
        anchor="custom",
        upstream_bp=upstream_bp,
        downstream_bp=downstream_bp,
    )

    with pytest.raises(ValueError, match="non-negative"):
        plotting.validate_axis_spec(axis, plot_family="aggregate_profile")


def test_prepare_single_read_plot_data_rejects_invalid_region_strand_in_region_5to3():
    reads = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "region_strand": "?",
                "event_pos": 110,
                "anchor": 100,
                "read_id": "r1",
            }
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="custom",
        upstream_bp=20,
        downstream_bp=20,
    )

    with pytest.raises(ValueError, match="region_strand values"):
        plotting.prepare_single_read_plot_data(
            reads,
            plot_family="single_read_raster",
            axis=axis,
            position_column="event_pos",
            anchor_column="anchor",
            region_strand_column="region_strand",
        )


def test_prepare_aggregate_plot_data_retains_metadata_for_fixed_window():
    table = pd.DataFrame(
        [
            {"region_id": "reg1", "region_strand": "+", "event_pos": 95, "anchor": 100, "signal": 1.0},
            {"region_id": "reg1", "region_strand": "+", "event_pos": 105, "anchor": 100, "signal": 3.0},
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=10,
        downstream_bp=10,
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="faceted",
    )

    payload = plotting.prepare_aggregate_plot_data(
        table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="signal",
        position_column="event_pos",
        anchor_column="anchor",
        region_strand_column="region_strand",
    )

    assert {"plot_table", "axis_table", "metadata"} <= set(payload)
    assert payload["metadata"]["orientation"] == "region_5to3"


@pytest.mark.parametrize(
    "relative, regions_dict, expected_anchor, expected_offset_center",
    [
        (True, None, "center", 0),
        (False, {"regions.bed": [(100, 140)]}, "absolute", 120),
    ],
)
def test_legacy_enrichment_profile_routes_regions_5to3prime_through_shared_prep(
    monkeypatch,
    relative,
    regions_dict,
    expected_anchor,
    expected_offset_center,
):
    called = {}
    captured = {}

    real_prepare_aggregate_plot_data = plotting.prepare_aggregate_plot_data

    def spy_prepare_aggregate_plot_data(table, *, plot_family, axis, aggregation, **kwargs):
        called["table"] = table
        called["plot_family"] = plot_family
        called["axis"] = axis
        called["aggregation"] = aggregation
        return real_prepare_aggregate_plot_data(
            table,
            plot_family=plot_family,
            axis=axis,
            aggregation=aggregation,
            **kwargs,
        )

    monkeypatch.setattr(plotting, "prepare_aggregate_plot_data", spy_prepare_aggregate_plot_data)
    monkeypatch.setattr(
        "dimelo.plot_enrichment_profile.get_enrichment_profiles",
        lambda **kwargs: [np.array([0.25, 0.75])],
    )
    if regions_dict is not None:
        monkeypatch.setattr(
            "dimelo.plot_enrichment_profile.utils.regions_dict_from_input",
            lambda *args, **kwargs: regions_dict,
        )

    def fake_make_enrichment_profile_plot(*, trace_vectors, sample_names, offset_center=0, **kwargs):
        captured["trace_vectors"] = trace_vectors
        captured["sample_names"] = sample_names
        captured["offset_center"] = offset_center
        return object()

    monkeypatch.setattr(
        "dimelo.plot_enrichment_profile.make_enrichment_profile_plot",
        fake_make_enrichment_profile_plot,
    )

    result = plot_enrichment_profile.plot_enrichment_profile(
        mod_file_names=["sample.fake"],
        regions_list=["regions.bed"],
        motifs=["A"],
        sample_names=["sample"],
        window_size=10,
        relative=relative,
        regions_5to3prime=True,
    )

    assert result is not None
    assert called["plot_family"] == "aggregate_profile"
    assert called["axis"].orientation == "region_5to3"
    assert called["axis"].coordinate_mode == "fixed_window"
    assert called["axis"].anchor == expected_anchor
    assert called["aggregation"].layout == "faceted"
    assert list(called["table"]["sample_name"]) == ["sample", "sample"]
    assert list(called["table"]["value"]) == [0.25, 0.75]
    assert list(called["table"]["region_strand"]) == ["+", "+"]
    assert len(captured["trace_vectors"]) == 1
    assert np.array_equal(captured["trace_vectors"][0], np.array([0.25, 0.75]))
    assert captured["sample_names"] == ["sample"]
    assert captured["offset_center"] == expected_offset_center


def test_prepare_single_read_plot_data_rejects_segment_map_axes():
    reads = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "region_strand": "+",
                "event_pos": 110,
                "anchor": 100,
                "read_id": "r1",
            }
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec(
                segment_id="body",
                label="Body",
                start_ref=100,
                end_ref=200,
            )
        ],
    )

    with pytest.raises(ValueError, match="segment_map is aggregate-only"):
        plotting.prepare_single_read_plot_data(
            reads,
            plot_family="single_read_raster",
            axis=axis,
            position_column="event_pos",
            anchor_column="anchor",
            region_strand_column="region_strand",
        )


def test_prepare_single_read_plot_data_rejects_scaled_segment_axes():
    reads = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "region_strand": "+",
                "event_pos": 110,
                "anchor": 100,
                "read_id": "r1",
            }
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec(
                segment_id="body",
                label="Body",
                start_ref=100,
                end_ref=200,
                mode="scaled",
                bins=25,
            )
        ],
    )

    with pytest.raises(ValueError, match="segment_map is aggregate-only"):
        plotting.prepare_single_read_plot_data(
            reads,
            plot_family="single_read_raster",
            axis=axis,
            position_column="event_pos",
            anchor_column="anchor",
            region_strand_column="region_strand",
        )


def test_plot_reads_routes_region_5to3prime_into_shared_single_read_prep(monkeypatch):
    captured = {}

    def fake_loader(**kwargs):
        return (
            [np.array([110.0])],
            np.array(["read-1"]),
            np.array(["A,0"]),
            {"chr1": [(100, 120, "-")]},
        )

    real_prepare = plotting.prepare_single_read_plot_data

    def spy_prepare(table, *, plot_family, axis, **kwargs):
        captured["table"] = table
        captured["plot_family"] = plot_family
        captured["axis"] = axis
        captured["kwargs"] = kwargs
        return real_prepare(table, plot_family=plot_family, axis=axis, **kwargs)

    class FakeAxes:
        def __init__(self):
            self.legend_ = None
            self.xlim = None

        def legend(self, *args, **kwargs):
            self.legend_args = (args, kwargs)

        def get_legend_handles_labels(self):
            return ([], [])

        def set_xlim(self, limits):
            self.xlim = limits

    monkeypatch.setattr(
        "dimelo.plot_reads.load_processed.readwise_binary_modification_arrays",
        fake_loader,
    )
    monkeypatch.setattr(
        "dimelo.plot_reads.plotting.prepare_single_read_plot_data",
        spy_prepare,
    )
    monkeypatch.setattr("dimelo.plot_reads.sns.scatterplot", lambda **kwargs: FakeAxes())

    axes = plot_reads.plot_reads(
        mod_file_name="sample.h5",
        regions="regions.bed",
        motifs=["A,0"],
        window_size=20,
        relative=True,
        regions_5to3prime=True,
    )

    assert isinstance(axes, FakeAxes)
    assert captured["plot_family"] == "single_read_raster"
    assert captured["axis"].orientation == "region_5to3"
    assert captured["axis"].coordinate_mode == "fixed_window"
    assert captured["axis"].anchor == "center"
    assert captured["axis"].upstream_bp == 20
    assert captured["axis"].downstream_bp == 20
    assert list(captured["table"]["region_strand"]) == ["+"]


def test_prepare_aggregate_plot_data_builds_concatenated_segment_axis():
    table = pd.DataFrame(
        [
            {"region_id": "reg1", "segment_id": "upstream", "segment_pos": 0, "signal": 1.0},
            {"region_id": "reg1", "segment_id": "body", "segment_pos": 10, "signal": 2.0},
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream", 0, 100, "raw", bins=20),
            plotting.SegmentSpec("body", "Body", 100, 400, "scaled", bins=50, contiguous_with_previous=True),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    payload = plotting.prepare_aggregate_plot_data(
        table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="signal",
        segment_id_column="segment_id",
        segment_position_column="segment_pos",
    )

    assert list(payload["axis_table"]["segment_id"]) == ["upstream", "body"]
    assert list(payload["axis_table"]["plot_start"]) == [0, 20]
    assert list(payload["axis_table"]["plot_end"]) == [20, 70]
    assert payload["axis_table"]["plot_start"].is_monotonic_increasing
    assert list(payload["plot_table"]["plot_x"]) == [0.0, 30.0]


def test_prepare_aggregate_plot_data_rejects_duplicate_segment_ids_in_axis_segments():
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "upstream", "segment_pos": 1, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream A", 0, 100, "raw", bins=20),
            plotting.SegmentSpec("upstream", "Upstream B", 100, 200, "raw", bins=20),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    with pytest.raises(ValueError, match="duplicate segment_id values"):
        plotting.prepare_aggregate_plot_data(
            table,
            plot_family="aggregate_profile",
            axis=axis,
            aggregation=aggregation,
            value_column="signal",
            segment_id_column="segment_id",
            segment_position_column="segment_pos",
        )


def test_prepare_aggregate_plot_data_rejects_unknown_segment_ids():
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "unknown", "segment_pos": 1, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream", 0, 100, "raw", bins=20),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    with pytest.raises(ValueError, match="unknown segment_id values"):
        plotting.prepare_aggregate_plot_data(
            table,
            plot_family="aggregate_profile",
            axis=axis,
            aggregation=aggregation,
            value_column="signal",
            segment_id_column="segment_id",
            segment_position_column="segment_pos",
        )


def test_prepare_aggregate_plot_data_rejects_missing_segment_positions():
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "upstream", "segment_pos": float("nan"), "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream", 0, 100, "raw", bins=20),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    with pytest.raises(ValueError, match="missing or NaN values"):
        plotting.prepare_aggregate_plot_data(
            table,
            plot_family="aggregate_profile",
            axis=axis,
            aggregation=aggregation,
            value_column="signal",
            segment_id_column="segment_id",
            segment_position_column="segment_pos",
        )


@pytest.mark.parametrize("segment_pos", [-1, 20])
def test_prepare_aggregate_plot_data_rejects_segment_positions_outside_declared_span(segment_pos):
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "upstream", "segment_pos": segment_pos, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream", 0, 100, "raw", bins=20),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    with pytest.raises(ValueError, match="declared segment span"):
        plotting.prepare_aggregate_plot_data(
            table,
            plot_family="aggregate_profile",
            axis=axis,
            aggregation=aggregation,
            value_column="signal",
            segment_id_column="segment_id",
            segment_position_column="segment_pos",
        )


def test_prepare_aggregate_plot_data_marks_non_contiguous_segment_breaks():
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "exon1", "segment_pos": 5, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("exon1", "Exon 1", 100, 200, "scaled", bins=20),
            plotting.SegmentSpec(
                "exon3",
                "Exon 3",
                500,
                650,
                "scaled",
                bins=20,
                contiguous_with_previous=False,
                plot_gap_after=True,
            ),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="faceted",
    )

    payload = plotting.prepare_aggregate_plot_data(
        table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="signal",
        segment_id_column="segment_id",
        segment_position_column="segment_pos",
    )

    assert list(payload["axis_table"]["contiguous_with_previous"]) == [True, False]
    assert list(payload["axis_table"]["plot_gap_after"]) == [False, True]


@pytest.mark.parametrize(
    "segment",
    [
        plotting.SegmentSpec("bad_bins", "Bad bins", 0, 100, "raw", bins=0),
        plotting.SegmentSpec("bad_bins_neg", "Bad bins negative", 0, 100, "scaled", bins=-5),
        plotting.SegmentSpec("bad_raw", "Bad raw", 100, 100, "raw"),
        plotting.SegmentSpec("bad_raw_reverse", "Bad raw reverse", 200, 100, "raw", bins=20),
    ],
)
def test_prepare_aggregate_plot_data_rejects_malformed_segment_spans(segment):
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "bad_raw", "segment_pos": 1, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[segment],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    with pytest.raises(ValueError, match="invalid"):
        plotting.prepare_aggregate_plot_data(
            table,
            plot_family="aggregate_profile",
            axis=axis,
            aggregation=aggregation,
            value_column="signal",
            segment_id_column="segment_id",
            segment_position_column="segment_pos",
        )


def test_prepare_aggregate_plot_data_rejects_invalid_segment_mode():
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "upstream", "segment_pos": 1, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream", 0, 100, "bogus", bins=20),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    with pytest.raises(ValueError, match="mode must be 'raw' or 'scaled'"):
        plotting.prepare_aggregate_plot_data(
            table,
            plot_family="aggregate_profile",
            axis=axis,
            aggregation=aggregation,
            value_column="signal",
            segment_id_column="segment_id",
            segment_position_column="segment_pos",
        )


@pytest.mark.parametrize("start_ref, end_ref", [(100, 100), (200, 100)])
def test_prepare_aggregate_plot_data_rejects_zero_or_negative_width_scaled_segment_without_bins(
    start_ref, end_ref
):
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "body", "segment_pos": 1, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("body", "Body", start_ref, end_ref, "scaled"),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    with pytest.raises(ValueError, match="invalid scaled span"):
        plotting.prepare_aggregate_plot_data(
            table,
            plot_family="aggregate_profile",
            axis=axis,
            aggregation=aggregation,
            value_column="signal",
            segment_id_column="segment_id",
            segment_position_column="segment_pos",
        )


def test_prepare_aggregate_plot_data_segment_map_keeps_user_plot_columns():
    table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "segment_id": "upstream",
                "segment_pos": 3,
                "plot_start": "user-start",
                "plot_end": "user-end",
                "signal": 1.0,
            }
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream", 0, 100, "raw", bins=20),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    payload = plotting.prepare_aggregate_plot_data(
        table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="signal",
        segment_id_column="segment_id",
        segment_position_column="segment_pos",
    )

    assert list(payload["plot_table"]["plot_start"]) == ["user-start"]
    assert list(payload["plot_table"]["plot_end"]) == ["user-end"]
    assert list(payload["plot_table"]["plot_x"]) == [3.0]
