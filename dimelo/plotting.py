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
