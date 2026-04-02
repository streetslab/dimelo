import pytest
import pandas as pd

from dimelo import plotting


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
