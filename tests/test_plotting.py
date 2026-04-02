import numpy as np
import pandas as pd
import pytest

from dimelo.models import ContrastSpec, RegionContrastResult
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
