from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from dimelo.distribution import _require_columns


@dataclass(frozen=True)
class SegmentSpec:
    segment_id: str
    label: str
    start_ref: int
    end_ref: int
    mode: str = "raw"
    bins: int | None = None
    plot_gap_after: bool = False
    contiguous_with_previous: bool = True


@dataclass(frozen=True)
class AxisSpec:
    orientation: str
    coordinate_mode: str
    anchor: str | None = None
    upstream_bp: int | None = None
    downstream_bp: int | None = None
    segments: list[SegmentSpec] | None = None


@dataclass(frozen=True)
class AggregationSpec:
    weighting: str = "equal_region"
    within_region_summary: str = "mean"
    signal_normalization: str = "none"
    layout: str = "faceted"


def validate_axis_spec(axis: AxisSpec, *, plot_family: str) -> None:
    if axis.orientation not in {"genomic", "region_5to3"}:
        raise ValueError("AxisSpec.orientation must be 'genomic' or 'region_5to3'.")

    if axis.coordinate_mode not in {"fixed_window", "segment_map"}:
        raise ValueError("AxisSpec.coordinate_mode must be 'fixed_window' or 'segment_map'.")

    if axis.coordinate_mode == "segment_map" and not axis.segments:
        raise ValueError("segment_map requires segments.")

    if axis.coordinate_mode == "fixed_window":
        if axis.upstream_bp is None or axis.downstream_bp is None:
            raise ValueError("fixed_window requires upstream_bp and downstream_bp.")
        if axis.upstream_bp < 0 or axis.downstream_bp < 0:
            raise ValueError("fixed_window upstream_bp and downstream_bp must be non-negative.")

    if plot_family == "single_read_raster" and axis.segments:
        if any(segment.mode == "scaled" for segment in axis.segments):
            raise ValueError(
                "single_read_raster plots must preserve coordinates and cannot use scaled segments."
            )


def validate_aggregation_spec(spec: AggregationSpec) -> None:
    if spec.weighting not in {"equal_region", "equal_observation", "coverage_weighted"}:
        raise ValueError("Unsupported weighting mode.")

    if spec.within_region_summary not in {"mean", "fraction", "density"}:
        raise ValueError("Unsupported within_region_summary.")

    if spec.signal_normalization not in {"none", "per_region", "global", "control_regions"}:
        raise ValueError("Unsupported signal_normalization.")

    if spec.layout not in {"concatenated", "faceted"}:
        raise ValueError("Unsupported layout mode.")


def _relative_position(position: float, anchor: float) -> float:
    return float(position) - float(anchor)


def _orient_position(relative_position: float, region_strand: str, orientation: str) -> float:
    if orientation != "region_5to3":
        return relative_position

    if region_strand not in {"+", "-"}:
        raise ValueError("region_5to3 fixed_window prep requires region_strand values of '+' or '-'.")

    if region_strand == "-":
        return -relative_position
    return relative_position


def _prepare_fixed_window_plot_data(
    table: pd.DataFrame,
    *,
    axis: AxisSpec,
    position_column: str,
    anchor_column: str,
    region_strand_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if axis.coordinate_mode != "fixed_window":
        raise ValueError("fixed_window prep only supports coordinate_mode='fixed_window'.")

    _require_columns(table, (position_column, anchor_column), "plot_table")
    if axis.orientation == "region_5to3":
        if region_strand_column is None:
            raise ValueError("region_5to3 fixed_window prep requires region_strand_column.")
        _require_columns(table, (region_strand_column,), "plot_table")

    plot_table = table.copy()
    relative_position = plot_table[position_column].astype(float) - plot_table[anchor_column].astype(float)
    if axis.orientation == "region_5to3":
        plot_table["plot_x"] = [
            _orient_position(relative_position=rel, region_strand=strand, orientation=axis.orientation)
            for rel, strand in zip(relative_position, plot_table[region_strand_column], strict=True)
        ]
    else:
        plot_table["plot_x"] = relative_position

    window_mask = plot_table["plot_x"].between(-axis.upstream_bp, axis.downstream_bp)
    plot_table = plot_table.loc[window_mask].copy().reset_index(drop=True)

    axis_table = pd.DataFrame(
        [
            {
                "segment_id": "window",
                "axis_min": -axis.upstream_bp,
                "axis_max": axis.downstream_bp,
            }
        ]
    )
    return plot_table, axis_table


def prepare_single_read_plot_data(
    table: pd.DataFrame,
    *,
    plot_family: str,
    axis: AxisSpec,
    position_column: str,
    anchor_column: str,
    region_strand_column: str | None = None,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family=plot_family)
    plot_table, axis_table = _prepare_fixed_window_plot_data(
        table,
        axis=axis,
        position_column=position_column,
        anchor_column=anchor_column,
        region_strand_column=region_strand_column,
    )
    metadata = {
        "plot_family": plot_family,
        "orientation": axis.orientation,
        "coordinate_mode": axis.coordinate_mode,
    }
    return {"plot_table": plot_table, "axis_table": axis_table, "metadata": metadata}


def prepare_aggregate_plot_data(
    table: pd.DataFrame,
    *,
    plot_family: str,
    axis: AxisSpec,
    aggregation: AggregationSpec,
    value_column: str,
    position_column: str,
    anchor_column: str,
    region_strand_column: str | None = None,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family=plot_family)
    validate_aggregation_spec(aggregation)
    _require_columns(table, (value_column,), "plot_table")

    plot_table, axis_table = _prepare_fixed_window_plot_data(
        table,
        axis=axis,
        position_column=position_column,
        anchor_column=anchor_column,
        region_strand_column=region_strand_column,
    )
    metadata = {
        "plot_family": plot_family,
        "orientation": axis.orientation,
        "coordinate_mode": axis.coordinate_mode,
        "weighting": aggregation.weighting,
        "within_region_summary": aggregation.within_region_summary,
        "signal_normalization": aggregation.signal_normalization,
        "layout": aggregation.layout,
    }
    return {"plot_table": plot_table, "axis_table": axis_table, "metadata": metadata}


def prepare_cluster_distribution_bar_data(cluster_distribution: pd.DataFrame) -> pd.DataFrame:
    _require_columns(
        cluster_distribution,
        ("sample_id", "condition", "cluster", "count", "fraction"),
        "cluster_distribution",
    )
    if cluster_distribution.empty:
        return cluster_distribution.loc[:, ["sample_id", "condition", "cluster", "count", "fraction"]].copy()

    return (
        cluster_distribution.loc[:, ["sample_id", "condition", "cluster", "count", "fraction"]]
        .sort_values(["sample_id", "condition", "cluster"], kind="stable")
        .reset_index(drop=True)
    )


def prepare_cluster_distribution_heatmap_data(condition_distribution: pd.DataFrame) -> pd.DataFrame:
    _require_columns(
        condition_distribution,
        ("condition", "cluster", "fraction"),
        "condition_distribution",
    )
    if condition_distribution.empty:
        return pd.DataFrame(columns=["condition"])

    heatmap = (
        condition_distribution.pivot_table(
            index="condition",
            columns="cluster",
            values="fraction",
            aggfunc="sum",
            fill_value=0.0,
        )
        .sort_index(axis=0)
        .reindex(sorted(condition_distribution["cluster"].unique()), axis=1, fill_value=0.0)
        .reset_index()
    )
    heatmap.columns.name = None
    return heatmap
