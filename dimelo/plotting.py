from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import numpy as np
import pandas as pd

from dimelo.distribution import _require_columns

_SEGMENT_AXIS_INTERNAL_SEGMENT_ID = "__dimelo_segment_axis_segment_id"
_SEGMENT_AXIS_INTERNAL_PLOT_START = "__dimelo_segment_axis_plot_start"
_SEGMENT_AXIS_INTERNAL_PLOT_END = "__dimelo_segment_axis_plot_end"
_SEGMENT_AXIS_INTERNAL_SEGMENT_SPAN = "__dimelo_segment_axis_segment_span"
_SEGMENT_AXIS_INTERNAL_COLUMNS = {
    _SEGMENT_AXIS_INTERNAL_SEGMENT_ID,
    _SEGMENT_AXIS_INTERNAL_PLOT_START,
    _SEGMENT_AXIS_INTERNAL_PLOT_END,
    _SEGMENT_AXIS_INTERNAL_SEGMENT_SPAN,
}


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
        raise ValueError(
            "AxisSpec.coordinate_mode must be 'fixed_window' or 'segment_map'."
        )

    if axis.coordinate_mode == "segment_map" and not axis.segments:
        raise ValueError("segment_map requires segments.")

    if plot_family == "single_read_raster" and axis.coordinate_mode == "segment_map":
        raise ValueError("segment_map is aggregate-only for single_read_raster plots.")

    if axis.coordinate_mode == "fixed_window":
        if axis.upstream_bp is None or axis.downstream_bp is None:
            raise ValueError("fixed_window requires upstream_bp and downstream_bp.")
        if axis.upstream_bp < 0 or axis.downstream_bp < 0:
            raise ValueError(
                "fixed_window upstream_bp and downstream_bp must be non-negative."
            )

    if (
        plot_family == "single_read_raster"
        and axis.segments
        and any(segment.mode == "scaled" for segment in axis.segments)
    ):
        raise ValueError(
            "single_read_raster plots must preserve coordinates and cannot use scaled segments."
        )


def validate_aggregation_spec(spec: AggregationSpec) -> None:
    if spec.weighting not in {"equal_region", "equal_observation", "coverage_weighted"}:
        raise ValueError("Unsupported weighting mode.")

    if spec.within_region_summary not in {"mean", "fraction", "density"}:
        raise ValueError("Unsupported within_region_summary.")

    if spec.signal_normalization not in {
        "none",
        "per_region",
        "global",
        "control_regions",
    }:
        raise ValueError("Unsupported signal_normalization.")

    if spec.layout not in {"concatenated", "faceted"}:
        raise ValueError("Unsupported layout mode.")


def _region_contrast_grouping_key(result, position_table: pd.DataFrame) -> str:
    available_keys = [
        column
        for column in ("sample_id", "condition")
        if column in position_table.columns
    ]
    if not available_keys:
        raise ValueError(
            "region contrast plotting requires position_table to include sample_id or condition."
        )
    if len(available_keys) == 1:
        return available_keys[0]

    contrast_values: set[str] = set()
    contrast = getattr(result, "contrast", None)
    if contrast is not None:
        contrast_values.update(str(value) for value in (contrast.numerator or []))
        contrast_values.update(str(value) for value in (contrast.denominator or []))

    sample_ids = {str(value) for value in position_table["sample_id"].dropna().unique()}
    conditions = {str(value) for value in position_table["condition"].dropna().unique()}
    sample_match = bool(contrast_values) and contrast_values.issubset(sample_ids)
    condition_match = bool(contrast_values) and contrast_values.issubset(conditions)

    if sample_match and not condition_match:
        return "sample_id"
    if condition_match and not sample_match:
        return "condition"

    raise ValueError(
        "region contrast plotting could not resolve a unique grouping key from "
        "position_table columns sample_id and condition."
    )


def _region_contrast_metadata(result) -> dict[str, object]:
    metadata = getattr(result, "metadata", {}) or {}
    return {
        "analysis_unit": metadata.get("analysis_unit"),
        "representation": metadata.get("representation"),
        "signal_source": metadata.get("signal_source"),
        "test": metadata.get("test"),
    }


def _validate_global_analysis_result(result) -> None:
    if result is None:
        raise ValueError("plotting helpers require a GlobalAnalysisResult.")
    required_attrs = ("summary", "windows", "normalization_factors", "metadata")
    if not all(hasattr(result, attr) for attr in required_attrs):
        raise TypeError("plotting helpers require a GlobalAnalysisResult-like object.")


def _validate_shared_cluster_result(result) -> None:
    if result is None:
        raise ValueError("plotting helpers require a SharedClusterResult.")
    required_attrs = (
        "model",
        "cluster_distribution",
        "condition_distribution",
        "cluster_profiles",
        "plot_data",
        "metadata",
    )
    if not all(hasattr(result, attr) for attr in required_attrs):
        raise TypeError("plotting helpers require a SharedClusterResult-like object.")
    if not hasattr(result.model, "mode") or not hasattr(result.model, "cluster_labels"):
        raise TypeError("plotting helpers require a SharedClusterResult-like object.")


def _empty_distribution_change_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "condition",
            "cluster",
            "count",
            "fraction",
            "replicate_n",
            "reference_fraction",
            "delta_fraction",
            "log2_fc",
        ]
    )


def _ordered_non_null_values(table: pd.DataFrame, column: str) -> list[object]:
    if column not in table.columns:
        return []
    return table.loc[table[column].notna(), column].drop_duplicates().tolist()


def _filter_motif_table(
    table: pd.DataFrame, motifs: list[str] | None, *, owner: str
) -> pd.DataFrame:
    _require_columns(table, ("motif",), owner)
    if motifs is None:
        return table.copy()

    filtered = table.loc[table["motif"].isin(motifs)].copy()
    if filtered.empty:
        raise ValueError(f"Requested motifs are not present in {owner}.")
    return filtered


def _validate_region_discovery_result(result) -> None:
    if result is None:
        raise ValueError("plotting helpers require a RegionDiscoveryResult.")
    if (
        not hasattr(result, "windows")
        or not hasattr(result, "hits")
        or not hasattr(result, "metadata")
    ):
        raise TypeError("plotting helpers require a RegionDiscoveryResult-like object.")


def _select_discovery_score_column(
    windows: pd.DataFrame, score_column: str | None
) -> str:
    if score_column is not None:
        if score_column not in windows.columns:
            raise ValueError(f"Unknown discovery score column: {score_column}")
        return score_column
    if "score_value" in windows.columns or windows.empty:
        return "score_value"
    raise ValueError(
        "Could not infer a discovery score column from RegionDiscoveryResult.windows."
    )


def _empty_discovery_scan_table(score_column: str) -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "window_id",
            "contig",
            "start",
            "end",
            "strand",
            score_column,
            "window_midpoint",
            "is_hit",
        ]
    )


def _discovery_contig_column(table: pd.DataFrame, *, owner: str) -> str:
    if "chromosome" in table.columns:
        return "chromosome"
    if "chrom" in table.columns:
        return "chrom"
    raise ValueError(f"{owner} must include either 'chromosome' or 'chrom'.")


def _sort_discovery_table(
    table: pd.DataFrame,
    *,
    contig_order: list[str],
    sort_columns: list[str],
) -> pd.DataFrame:
    if table.empty:
        return table.reset_index(drop=True)

    sorted_table = table.copy()
    sorted_table["__dimelo_contig_order"] = pd.Categorical(
        sorted_table["contig"],
        categories=contig_order,
        ordered=True,
    )
    sorted_table = sorted_table.sort_values(
        ["__dimelo_contig_order", *sort_columns],
        kind="stable",
    ).drop(columns="__dimelo_contig_order")
    return sorted_table.reset_index(drop=True)


def _select_region_discovery_hits(
    hits: pd.DataFrame,
    *,
    top_n: int | None,
) -> tuple[pd.DataFrame, str]:
    if top_n is not None and top_n < 0:
        raise ValueError("top_n must be non-negative.")
    if hits.empty:
        return hits.copy(), "top_n"

    if "rank" in hits.columns:
        selected = hits.sort_values(["rank"], kind="stable")
    else:
        selected = hits.copy()
    if top_n is not None:
        selected = selected.head(top_n)
    return selected.copy().reset_index(drop=True), "top_n"


def _relative_position(position: float, anchor: float) -> float:
    return float(position) - float(anchor)


def _orient_position(
    relative_position: float, region_strand: str, orientation: str
) -> float:
    if orientation != "region_5to3":
        return relative_position

    if region_strand not in {"+", "-"}:
        raise ValueError(
            "region_5to3 fixed_window prep requires region_strand values of '+' or '-'."
        )

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
        raise ValueError(
            "fixed_window prep only supports coordinate_mode='fixed_window'."
        )

    _require_columns(table, (position_column, anchor_column), "plot_table")
    if axis.orientation == "region_5to3":
        if region_strand_column is None:
            raise ValueError(
                "region_5to3 fixed_window prep requires region_strand_column."
            )
        _require_columns(table, (region_strand_column,), "plot_table")

    plot_table = table.copy()
    relative_position = plot_table[position_column].astype(float) - plot_table[
        anchor_column
    ].astype(float)
    if axis.orientation == "region_5to3":
        plot_table["plot_x"] = [
            _orient_position(
                relative_position=rel,
                region_strand=strand,
                orientation=axis.orientation,
            )
            for rel, strand in zip(
                relative_position, plot_table[region_strand_column], strict=True
            )
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


def _build_segment_axis_table(segments: list[SegmentSpec]) -> pd.DataFrame:
    segment_id_counts = Counter(segment.segment_id for segment in segments)
    duplicate_segment_ids = sorted(
        segment_id for segment_id, count in segment_id_counts.items() if count > 1
    )
    if duplicate_segment_ids:
        raise ValueError(
            "segment_map axis.segments contains duplicate segment_id values: "
            f"{', '.join(duplicate_segment_ids)}."
        )

    rows: list[dict[str, object]] = []
    running_start = 0

    for segment in segments:
        if segment.mode not in {"raw", "scaled"}:
            raise ValueError(
                "segment_map axis.segments contains invalid mode values: "
                f"{segment.segment_id} has mode={segment.mode!r}; mode must be 'raw' or 'scaled'."
            )
        if segment.bins is not None and segment.bins <= 0:
            raise ValueError(
                "segment_map axis.segments contains invalid bins values: "
                f"{segment.segment_id} has bins={segment.bins!r}."
            )
        if (
            segment.mode == "scaled"
            and segment.bins is None
            and segment.end_ref <= segment.start_ref
        ):
            raise ValueError(
                "segment_map axis.segments contains invalid scaled span values: "
                f"{segment.segment_id} has start_ref={segment.start_ref!r} and end_ref={segment.end_ref!r}."
            )
        if segment.mode == "raw" and segment.end_ref <= segment.start_ref:
            raise ValueError(
                "segment_map axis.segments contains invalid raw span values: "
                f"{segment.segment_id} has start_ref={segment.start_ref!r} and end_ref={segment.end_ref!r}."
            )

        span = (
            segment.bins
            if segment.bins is not None
            else segment.end_ref - segment.start_ref
        )
        rows.append(
            {
                "segment_id": segment.segment_id,
                "label": segment.label,
                "start_ref": segment.start_ref,
                "end_ref": segment.end_ref,
                "mode": segment.mode,
                "bins": segment.bins,
                "plot_start": running_start,
                "plot_end": running_start + span,
                "contiguous_with_previous": segment.contiguous_with_previous,
                "plot_gap_after": segment.plot_gap_after,
                "discontinuity_before": not segment.contiguous_with_previous,
                "discontinuity_after": segment.plot_gap_after,
            }
        )
        running_start += span

    return pd.DataFrame(rows)


def _prepare_segment_map_plot_data(
    table: pd.DataFrame,
    *,
    axis: AxisSpec,
    segment_id_column: str,
    segment_position_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not axis.segments:
        raise ValueError("segment_map requires segments.")

    _require_columns(table, (segment_id_column, segment_position_column), "plot_table")
    reserved_columns = _SEGMENT_AXIS_INTERNAL_COLUMNS.intersection(table.columns)
    if reserved_columns:
        raise ValueError(
            "segment_map plotting received plot_table columns reserved for internal use: "
            f"{', '.join(sorted(reserved_columns))}."
        )

    axis_table = _build_segment_axis_table(axis.segments)
    plot_table = table.copy()
    axis_lookup = axis_table.loc[:, ["segment_id", "plot_start", "plot_end"]].rename(
        columns={
            "segment_id": _SEGMENT_AXIS_INTERNAL_SEGMENT_ID,
            "plot_start": _SEGMENT_AXIS_INTERNAL_PLOT_START,
            "plot_end": _SEGMENT_AXIS_INTERNAL_PLOT_END,
        }
    )
    axis_lookup[_SEGMENT_AXIS_INTERNAL_SEGMENT_SPAN] = (
        axis_lookup[_SEGMENT_AXIS_INTERNAL_PLOT_END]
        - axis_lookup[_SEGMENT_AXIS_INTERNAL_PLOT_START]
    )
    plot_table = plot_table.merge(
        axis_lookup,
        left_on=segment_id_column,
        right_on=_SEGMENT_AXIS_INTERNAL_SEGMENT_ID,
        how="left",
    )
    unknown_segment_mask = plot_table[_SEGMENT_AXIS_INTERNAL_PLOT_START].isna()
    if unknown_segment_mask.any():
        unknown_segment_ids = sorted(
            plot_table.loc[unknown_segment_mask, segment_id_column].astype(str).unique()
        )
        raise ValueError(
            "segment_map plotting received unknown segment_id values: "
            f"{', '.join(unknown_segment_ids)}."
        )

    segment_positions = pd.to_numeric(
        plot_table[segment_position_column], errors="raise"
    )
    missing_position_mask = segment_positions.isna()
    if missing_position_mask.any():
        invalid_rows = plot_table.loc[
            missing_position_mask, [segment_id_column, segment_position_column]
        ]
        raise ValueError(
            "segment_position_column contains missing or NaN values. "
            f"Invalid rows: {invalid_rows.to_dict(orient='records')}."
        )

    invalid_position_mask = (segment_positions < 0) | (
        segment_positions >= plot_table[_SEGMENT_AXIS_INTERNAL_SEGMENT_SPAN]
    )
    if invalid_position_mask.any():
        invalid_rows = plot_table.loc[
            invalid_position_mask, [segment_id_column, segment_position_column]
        ].copy()
        invalid_rows["segment_span"] = plot_table.loc[
            invalid_position_mask, _SEGMENT_AXIS_INTERNAL_SEGMENT_SPAN
        ].values
        raise ValueError(
            "segment_position_column values must stay within the declared segment span "
            "for each segment."
            f" Invalid rows: {invalid_rows.to_dict(orient='records')}."
        )

    plot_table["plot_x"] = (
        plot_table[_SEGMENT_AXIS_INTERNAL_PLOT_START] + segment_positions
    )
    plot_table = plot_table.drop(
        columns=[
            _SEGMENT_AXIS_INTERNAL_SEGMENT_ID,
            _SEGMENT_AXIS_INTERNAL_PLOT_START,
            _SEGMENT_AXIS_INTERNAL_PLOT_END,
            _SEGMENT_AXIS_INTERNAL_SEGMENT_SPAN,
        ]
    )

    return plot_table, axis_table


def _prepare_region_contrast_position_table(
    *,
    result,
    position_table: pd.DataFrame,
    grouping_key: str,
) -> pd.DataFrame:
    _require_columns(
        position_table,
        ("region_id", grouping_key, "position", "anchor", "value", "region_strand"),
        "position_table",
    )
    summary_region_ids = result.summary["region_id"].dropna().astype(str)
    filtered = position_table[
        position_table["region_id"].astype(str).isin(summary_region_ids)
    ].copy()
    if filtered.empty:
        raise ValueError(
            "position_table does not contain any region_id values present in the contrast result."
        )
    return filtered.reset_index(drop=True)


def _prepare_region_contrast_value_modes(
    *,
    result,
    position_table: pd.DataFrame,
    grouping_key: str,
) -> pd.DataFrame:
    contrast = result.contrast
    numerator_mask = position_table[grouping_key].isin(contrast.numerator or [])
    denominator_mask = position_table[grouping_key].isin(contrast.denominator or [])

    numerator = position_table.loc[numerator_mask].copy()
    denominator = position_table.loc[denominator_mask].copy()
    if numerator.empty or denominator.empty:
        raise ValueError(
            "position_table does not contain rows for both contrast sides."
        )

    join_keys = ["region_id", "position", "anchor", "region_strand"]
    numerator = (
        numerator.loc[:, join_keys + ["value"]]
        .groupby(join_keys, as_index=False, sort=False)
        .mean(numeric_only=True)
    )
    denominator = (
        denominator.loc[:, join_keys + ["value"]]
        .groupby(join_keys, as_index=False, sort=False)
        .mean(numeric_only=True)
    )

    coordinate_match = numerator.loc[:, join_keys].merge(
        denominator.loc[:, join_keys],
        on=join_keys,
        how="outer",
        indicator=True,
    )
    if not (coordinate_match["_merge"] == "both").all():
        mismatched_rows = coordinate_match.loc[
            coordinate_match["_merge"] != "both", join_keys + ["_merge"]
        ].copy()
        raise ValueError(
            "position_table contains mismatched coordinates between contrast sides. "
            f"Mismatched rows: {mismatched_rows.to_dict(orient='records')}."
        )

    paired = numerator.merge(
        denominator,
        on=join_keys,
        suffixes=("_numerator", "_denominator"),
        how="inner",
    )
    if paired.empty:
        raise ValueError(
            "Unable to compute delta because numerator and denominator positions do not align."
        )

    if "rank" in result.summary.columns:
        _require_columns(result.summary, ("region_id", "rank"), "result.summary")
        rank_table = result.summary.loc[:, ["region_id", "rank"]].drop_duplicates(
            subset=["region_id"]
        )
    else:
        rank_table = None

    def _attach_rank(table: pd.DataFrame) -> pd.DataFrame:
        if rank_table is None:
            return table
        return table.merge(rank_table, on="region_id", how="left")

    numerator_table = numerator.loc[
        :, ["region_id", "position", "anchor", "value", "region_strand"]
    ].copy()
    numerator_table["value_mode"] = "numerator"
    numerator_table[grouping_key] = "numerator"
    numerator_table = _attach_rank(numerator_table)

    denominator_table = denominator.loc[
        :, ["region_id", "position", "anchor", "value", "region_strand"]
    ].copy()
    denominator_table["value_mode"] = "denominator"
    denominator_table[grouping_key] = "denominator"
    denominator_table = _attach_rank(denominator_table)

    delta_table = paired.loc[
        :,
        [
            "region_id",
            "position",
            "anchor",
            "region_strand",
            "value_numerator",
            "value_denominator",
        ],
    ].copy()
    delta_table["value"] = (
        delta_table["value_numerator"] - delta_table["value_denominator"]
    )
    delta_table[grouping_key] = "delta"
    delta_table["value_mode"] = "delta"
    delta_table = delta_table.loc[
        :,
        [
            "region_id",
            grouping_key,
            "position",
            "anchor",
            "value",
            "region_strand",
            "value_mode",
        ],
    ]
    delta_table = _attach_rank(delta_table)

    ordered_columns = [
        "region_id",
        grouping_key,
        "position",
        "anchor",
        "value",
        "region_strand",
        "value_mode",
    ]
    if rank_table is not None:
        ordered_columns.append("rank")

    return pd.concat(
        [
            numerator_table.loc[:, ordered_columns],
            denominator_table.loc[:, ordered_columns],
            delta_table.loc[:, ordered_columns],
        ],
        ignore_index=True,
    )


def prepare_region_contrast_profile_data(
    *,
    result,
    position_table: pd.DataFrame,
    axis: AxisSpec,
    aggregation: AggregationSpec,
    value_mode: str = "all",
) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family="region_contrast_profile")
    validate_aggregation_spec(aggregation)

    grouping_key = _region_contrast_grouping_key(result, position_table)
    prepared_position_table = _prepare_region_contrast_position_table(
        result=result,
        position_table=position_table,
        grouping_key=grouping_key,
    )
    contrast_table = _prepare_region_contrast_value_modes(
        result=result,
        position_table=prepared_position_table,
        grouping_key=grouping_key,
    )

    if value_mode != "all":
        if value_mode not in {"numerator", "denominator", "delta"}:
            raise ValueError("Unsupported region contrast value_mode.")
        contrast_table = contrast_table.loc[
            contrast_table["value_mode"] == value_mode
        ].copy()

    prepared = prepare_aggregate_plot_data(
        contrast_table,
        plot_family="region_contrast_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="value",
        position_column="position",
        anchor_column="anchor",
        region_strand_column="region_strand",
    )
    metadata = dict(prepared["metadata"])
    metadata.update(_region_contrast_metadata(result))
    metadata["value_mode"] = value_mode
    prepared["metadata"] = metadata
    return prepared


def prepare_region_contrast_heatmap_data(
    *,
    result,
    position_table: pd.DataFrame,
    axis: AxisSpec,
    aggregation: AggregationSpec,
    value_mode: str = "all",
) -> dict[str, pd.DataFrame | dict[str, object]]:
    profile_payload = prepare_region_contrast_profile_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode=value_mode,
    )

    _require_columns(result.summary, ("region_id", "rank"), "result.summary")
    plot_region_ids = (
        profile_payload["plot_table"]
        .loc[:, ["region_id"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    summary_ranks = plot_region_ids.merge(
        result.summary.loc[:, ["region_id", "rank"]],
        on="region_id",
        how="left",
    )
    summary_ranks = summary_ranks.groupby("region_id", as_index=False, sort=False).agg(
        rank=(
            "rank",
            lambda values: (
                lambda non_na: non_na.iloc[0] if not non_na.empty else pd.NA
            )(pd.Series(values).dropna()),
        ),
        rank_count=("rank", lambda values: pd.Series(values).dropna().nunique()),
    )
    conflicting_rank_ids = sorted(
        summary_ranks.loc[summary_ranks["rank_count"] > 1, "region_id"]
        .astype(str)
        .unique()
    )
    if conflicting_rank_ids:
        raise ValueError(
            "result.summary must provide exactly one rank value per plotted region. "
            f"Conflicting region_id values: {', '.join(conflicting_rank_ids)}."
        )

    row_order = summary_ranks.loc[:, ["region_id", "rank"]].copy()
    missing_rank_ids = sorted(
        summary_ranks.loc[summary_ranks["rank_count"] == 0, "region_id"]
        .astype(str)
        .unique()
    )
    if missing_rank_ids:
        raise ValueError(
            "result.summary does not provide rank values for all plotted regions. "
            f"Missing region_id values: {', '.join(missing_rank_ids)}."
        )
    row_order = row_order.sort_values(["rank", "region_id"], kind="stable").reset_index(
        drop=True
    )
    row_order["row_order"] = range(len(row_order))

    plot_table = profile_payload["plot_table"].merge(
        row_order.loc[:, ["region_id", "row_order"]],
        on="region_id",
        how="left",
    )
    plot_table = plot_table.sort_values(
        ["row_order", "value_mode", "position"], kind="stable"
    ).reset_index(drop=True)

    profile_payload["plot_table"] = plot_table
    profile_payload["metadata"] = {
        **profile_payload["metadata"],
        "plot_family": "region_contrast_heatmap",
    }
    profile_payload["summary_table"] = row_order
    return profile_payload


def prepare_global_analysis_summary_data(
    *,
    result,
    motifs: list[str] | None = None,
    aggregate_conditions: bool = True,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_global_analysis_result(result)
    _require_columns(
        result.summary,
        (
            "sample_id",
            "condition",
            "motif",
            "modified_count",
            "valid_count",
            "global_fraction",
        ),
        "result.summary",
    )
    _require_columns(
        result.normalization_factors,
        (
            "sample_id",
            "condition",
            "motif",
            "global_fraction",
            "reference_fraction",
            "global_offset",
        ),
        "result.normalization_factors",
    )

    sample_summary = _filter_motif_table(result.summary, motifs, owner="result.summary")
    normalization_table = _filter_motif_table(
        result.normalization_factors,
        motifs,
        owner="result.normalization_factors",
    )

    if sample_summary.empty:
        condition_summary = pd.DataFrame(
            columns=[
                "condition",
                "motif",
                "global_fraction_mean",
                "global_fraction_median",
                "sample_n",
            ]
        )
        motif_values = list(motifs) if motifs is not None else []
    else:
        motif_values = sample_summary["motif"].drop_duplicates().tolist()
        if aggregate_conditions:
            condition_summary = (
                sample_summary.groupby(
                    ["condition", "motif"], as_index=False, sort=False
                )
                .agg(
                    global_fraction_mean=("global_fraction", "mean"),
                    global_fraction_median=("global_fraction", "median"),
                    sample_n=("sample_id", "nunique"),
                )
                .copy()
            )
        else:
            condition_summary = pd.DataFrame(
                columns=[
                    "condition",
                    "motif",
                    "global_fraction_mean",
                    "global_fraction_median",
                    "sample_n",
                ]
            )

    return {
        "sample_summary": sample_summary.reset_index(drop=True),
        "condition_summary": condition_summary.reset_index(drop=True),
        "normalization_table": normalization_table.reset_index(drop=True),
        "metadata": {
            "motifs": motif_values,
            "window_size": result.metadata.get("window_size"),
            "step_size": result.metadata.get("step_size"),
            "aggregate_conditions": aggregate_conditions,
        },
    }


def prepare_global_analysis_window_data(
    *,
    result,
    contigs: list[str] | None = None,
    motifs: list[str] | None = None,
    aggregate_conditions: bool = False,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_global_analysis_result(result)

    condition_window_columns = [
        "condition",
        "motif",
        "contig",
        "start",
        "end",
        "window_midpoint",
        "window_fraction_mean",
        "window_fraction_median",
        "sample_n",
    ]
    empty_window_columns = [
        "sample_id",
        "condition",
        "replicate",
        "motif",
        "window_id",
        "contig",
        "start",
        "end",
        "strand",
        "window_fraction",
        "window_midpoint",
    ]

    if result.windows is None:
        return {
            "window_table": pd.DataFrame(columns=empty_window_columns),
            "condition_window_table": pd.DataFrame(columns=condition_window_columns),
            "metadata": {
                "contig_order": [],
                "motifs": [] if motifs is None else list(motifs),
                "aggregate_conditions": aggregate_conditions,
            },
        }

    _require_columns(
        result.windows,
        (
            "sample_id",
            "condition",
            "replicate",
            "motif",
            "window_id",
            "start",
            "end",
            "strand",
            "window_fraction",
        ),
        "result.windows",
    )
    window_table = _filter_motif_table(result.windows, motifs, owner="result.windows")

    if window_table.empty:
        return {
            "window_table": pd.DataFrame(columns=empty_window_columns),
            "condition_window_table": pd.DataFrame(columns=condition_window_columns),
            "metadata": {
                "contig_order": [],
                "motifs": [] if motifs is None else list(motifs),
                "aggregate_conditions": aggregate_conditions,
            },
        }

    window_table = window_table.copy()
    windows_contig_column = _discovery_contig_column(
        window_table, owner="result.windows"
    )
    window_table["contig"] = window_table[windows_contig_column]

    if contigs is not None:
        window_table = window_table.loc[window_table["contig"].isin(contigs)].copy()
        if window_table.empty:
            raise ValueError("Requested contigs are not present in result.windows.")
        contig_order = list(contigs)
    else:
        contig_order = window_table["contig"].drop_duplicates().tolist()

    window_table["window_midpoint"] = (
        window_table["start"] + window_table["end"]
    ) / 2.0
    window_table = _sort_discovery_table(
        window_table,
        contig_order=contig_order,
        sort_columns=["start", "end", "sample_id"],
    )

    if aggregate_conditions:
        sample_window_table = (
            window_table.groupby(
                [
                    "sample_id",
                    "condition",
                    "motif",
                    "window_id",
                    "contig",
                    "start",
                    "end",
                    "window_midpoint",
                ],
                as_index=False,
                sort=False,
            )
            .agg(window_fraction=("window_fraction", "mean"))
            .copy()
        )
        condition_window_table = (
            sample_window_table.groupby(
                ["condition", "motif", "contig", "start", "end", "window_midpoint"],
                as_index=False,
                sort=False,
            )
            .agg(
                window_fraction_mean=("window_fraction", "mean"),
                window_fraction_median=("window_fraction", "median"),
                sample_n=("sample_id", "nunique"),
            )
            .copy()
        )
        condition_window_table = _sort_discovery_table(
            condition_window_table,
            contig_order=contig_order,
            sort_columns=["start", "end", "condition"],
        )
    else:
        condition_window_table = pd.DataFrame(columns=condition_window_columns)

    return {
        "window_table": window_table.reset_index(drop=True),
        "condition_window_table": condition_window_table.reset_index(drop=True),
        "metadata": {
            "contig_order": contig_order,
            "motifs": window_table["motif"].drop_duplicates().tolist(),
            "aggregate_conditions": aggregate_conditions,
        },
    }


def prepare_region_discovery_scan_data(
    *,
    result,
    contigs: list[str] | None = None,
    top_n_hits: int | None = 100,
    score_column: str | None = None,
    include_all_windows: bool = True,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_region_discovery_result(result)
    if top_n_hits is not None and top_n_hits < 0:
        raise ValueError("top_n_hits must be non-negative.")

    windows = result.windows.copy()
    hits = result.hits.copy()
    active_score_column = _select_discovery_score_column(windows, score_column)

    if windows.empty:
        empty_scan = _empty_discovery_scan_table(active_score_column)
        empty_hits = pd.DataFrame(
            columns=[
                "window_id",
                "contig",
                "start",
                "end",
                "strand",
                active_score_column,
            ]
        )
        return {
            "scan_table": empty_scan,
            "hit_table": empty_hits,
            "metadata": {
                "contig_order": [],
                "score_column": active_score_column,
                "score_mode": result.metadata.get("score"),
                "contrast_mode": result.metadata.get("contrast_mode"),
                "merge_hits": result.metadata.get("merge_hits"),
                "top_n_hits": top_n_hits,
            },
        }

    _require_columns(windows, ("window_id", "start", "end", "strand"), "result.windows")
    windows_contig_column = _discovery_contig_column(windows, owner="result.windows")
    windows["contig"] = windows[windows_contig_column]

    if hits.empty:
        hits = pd.DataFrame(columns=windows.columns)
    else:
        _require_columns(hits, ("window_id", "start", "end", "strand"), "result.hits")
        hits_contig_column = _discovery_contig_column(hits, owner="result.hits")
        hits["contig"] = hits[hits_contig_column]

    if contigs is not None:
        windows = windows.loc[windows["contig"].isin(contigs)].copy()
        hits = hits.loc[hits["contig"].isin(contigs)].copy()
        if windows.empty:
            raise ValueError(
                "Requested contigs are not present in RegionDiscoveryResult.windows."
            )
        contig_order = list(contigs)
    else:
        contig_order = windows["contig"].drop_duplicates().tolist()

    if top_n_hits is not None and not hits.empty and "rank" in hits.columns:
        hits = hits.sort_values(["rank"], kind="stable").head(top_n_hits).copy()

    filtered_window_ids = set(windows["window_id"].tolist())
    missing_hit_window_ids = sorted(
        str(window_id)
        for window_id in hits["window_id"].dropna().unique()
        if window_id not in filtered_window_ids
    )
    if missing_hit_window_ids:
        raise ValueError(
            "result.hits contains window_id values not present in the filtered windows table: "
            f"{', '.join(missing_hit_window_ids)}."
        )

    hit_window_ids = set(hits.get("window_id", pd.Series(dtype="object")).tolist())
    windows["window_midpoint"] = (windows["start"] + windows["end"]) / 2.0
    windows["is_hit"] = windows["window_id"].isin(hit_window_ids)

    scan_table = (
        windows if include_all_windows else windows.loc[windows["is_hit"]].copy()
    )
    scan_table = _sort_discovery_table(
        scan_table,
        contig_order=contig_order,
        sort_columns=["start", "end"],
    )

    hit_sort_columns = ["start", "end"]
    if "rank" in hits.columns:
        hit_sort_columns = ["rank", *hit_sort_columns]
    hit_table = _sort_discovery_table(
        hits,
        contig_order=contig_order,
        sort_columns=hit_sort_columns,
    )

    return {
        "scan_table": scan_table,
        "hit_table": hit_table,
        "metadata": {
            "contig_order": contig_order,
            "score_column": active_score_column,
            "score_mode": result.metadata.get("score"),
            "contrast_mode": result.metadata.get("contrast_mode"),
            "merge_hits": result.metadata.get("merge_hits"),
            "top_n_hits": top_n_hits,
        },
    }


def prepare_region_discovery_hit_context_data(
    *,
    result,
    top_n: int | None = 12,
    hit_ids: list[str] | None = None,
    padding_windows: int | None = 5,
    padding_bp: int | None = None,
    score_column: str | None = None,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_region_discovery_result(result)
    if padding_bp is not None:
        raise ValueError(
            "padding_bp is not supported for region discovery hit context prep."
        )
    if padding_windows is not None and padding_windows < 0:
        raise ValueError("padding_windows must be non-negative.")

    windows = result.windows.copy()
    hits = result.hits.copy()
    active_score_column = _select_discovery_score_column(windows, score_column)

    context_columns = [
        "window_id",
        "contig",
        "start",
        "end",
        "strand",
        active_score_column,
        "window_midpoint",
        "selected_hit_id",
        "selected_hit_rank",
        "relative_window_offset",
        "is_selected_hit",
    ]
    selection_mode = "top_n" if hit_ids is None else "explicit_ids"

    if hits.empty or windows.empty:
        return {
            "context_table": pd.DataFrame(columns=context_columns),
            "selected_hits": hits.copy().reset_index(drop=True),
            "metadata": {
                "selection_mode": selection_mode,
                "top_n": top_n,
                "padding_windows": padding_windows,
                "padding_bp": padding_bp,
                "score_column": active_score_column,
            },
        }

    _require_columns(windows, ("window_id", "start", "end", "strand"), "result.windows")
    _require_columns(hits, ("window_id", "start", "end", "strand"), "result.hits")
    windows_contig_column = _discovery_contig_column(windows, owner="result.windows")
    hits_contig_column = _discovery_contig_column(hits, owner="result.hits")
    windows["contig"] = windows[windows_contig_column]
    hits["contig"] = hits[hits_contig_column]
    windows["window_midpoint"] = (windows["start"] + windows["end"]) / 2.0

    contig_order = windows["contig"].drop_duplicates().tolist()
    windows = _sort_discovery_table(
        windows,
        contig_order=contig_order,
        sort_columns=["start", "end"],
    )

    if hit_ids is not None:
        selected_hits = (
            hits.loc[hits["window_id"].isin(hit_ids)].copy().reset_index(drop=True)
        )
    else:
        selected_hits, selection_mode = _select_region_discovery_hits(hits, top_n=top_n)

    if selected_hits.empty:
        return {
            "context_table": pd.DataFrame(columns=context_columns),
            "selected_hits": selected_hits,
            "metadata": {
                "selection_mode": selection_mode,
                "top_n": top_n,
                "padding_windows": padding_windows,
                "padding_bp": padding_bp,
                "score_column": active_score_column,
            },
        }

    available_window_ids = set(windows["window_id"].dropna().tolist())
    missing_selected_hit_ids = sorted(
        str(window_id)
        for window_id in selected_hits["window_id"].dropna().unique()
        if window_id not in available_window_ids
    )
    if missing_selected_hit_ids:
        raise ValueError(
            "selected hits contain window_id values not present in result.windows: "
            f"{', '.join(missing_selected_hit_ids)}."
        )

    windows_by_contig: dict[str, pd.DataFrame] = {}
    window_position_by_contig: dict[str, dict[str, int]] = {}
    for contig, contig_windows in windows.groupby("contig", sort=False):
        contig_table = contig_windows.reset_index(drop=True)
        windows_by_contig[contig] = contig_table
        window_position_by_contig[contig] = dict(
            zip(contig_table["window_id"], contig_table.index, strict=False)
        )

    context_frames: list[pd.DataFrame] = []
    padding = int(padding_windows or 0)
    has_rank = "rank" in selected_hits.columns

    for hit in selected_hits.itertuples(index=False):
        hit_contig = hit.contig
        hit_window_id = hit.window_id
        contig_windows = windows_by_contig.get(hit_contig)
        if contig_windows is None:
            continue

        hit_position = window_position_by_contig[hit_contig].get(hit_window_id)
        if hit_position is None:
            continue

        start_index = max(0, hit_position - padding)
        end_index = min(len(contig_windows), hit_position + padding + 1)
        context = contig_windows.iloc[start_index:end_index].copy()
        context["selected_hit_id"] = hit_window_id
        context["selected_hit_rank"] = hit.rank if has_rank else np.nan
        context["relative_window_offset"] = np.arange(
            start_index - hit_position, end_index - hit_position
        )
        context["is_selected_hit"] = context["window_id"] == hit_window_id
        context_frames.append(context)

    if context_frames:
        context_table = pd.concat(context_frames, ignore_index=True)
    else:
        context_table = pd.DataFrame(columns=context_columns)

    return {
        "context_table": context_table,
        "selected_hits": selected_hits,
        "metadata": {
            "selection_mode": selection_mode,
            "top_n": top_n,
            "padding_windows": padding_windows,
            "padding_bp": padding_bp,
            "score_column": active_score_column,
        },
    }


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
    position_column: str | None = None,
    anchor_column: str | None = None,
    region_strand_column: str | None = None,
    segment_id_column: str | None = None,
    segment_position_column: str | None = None,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family=plot_family)
    validate_aggregation_spec(aggregation)
    _require_columns(table, (value_column,), "plot_table")

    if axis.coordinate_mode == "segment_map":
        if segment_id_column is None or segment_position_column is None:
            raise ValueError(
                "segment_map plotting requires segment_id_column and segment_position_column."
            )
        plot_table, axis_table = _prepare_segment_map_plot_data(
            table,
            axis=axis,
            segment_id_column=segment_id_column,
            segment_position_column=segment_position_column,
        )
    else:
        if position_column is None or anchor_column is None:
            raise ValueError(
                "fixed_window plotting requires position_column and anchor_column."
            )
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
    if axis.coordinate_mode == "segment_map":
        metadata["segment_count"] = len(axis_table)
    return {"plot_table": plot_table, "axis_table": axis_table, "metadata": metadata}


def prepare_cluster_distribution_bar_data(
    cluster_distribution: pd.DataFrame,
) -> pd.DataFrame:
    _require_columns(
        cluster_distribution,
        ("sample_id", "condition", "cluster", "count", "fraction"),
        "cluster_distribution",
    )
    if cluster_distribution.empty:
        return cluster_distribution.loc[
            :, ["sample_id", "condition", "cluster", "count", "fraction"]
        ].copy()

    return (
        cluster_distribution.loc[
            :, ["sample_id", "condition", "cluster", "count", "fraction"]
        ]
        .sort_values(["sample_id", "condition", "cluster"], kind="stable")
        .reset_index(drop=True)
    )


def prepare_cluster_distribution_heatmap_data(
    condition_distribution: pd.DataFrame,
) -> pd.DataFrame:
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
        .reindex(
            sorted(condition_distribution["cluster"].unique()), axis=1, fill_value=0.0
        )
        .reset_index()
    )
    heatmap.columns.name = None
    return heatmap


def prepare_shared_cluster_distribution_data(
    *,
    result,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_shared_cluster_result(result)

    sample_distribution = prepare_cluster_distribution_bar_data(
        result.cluster_distribution
    )

    _require_columns(
        result.condition_distribution,
        ("condition", "cluster", "count", "fraction", "replicate_n"),
        "result.condition_distribution",
    )
    condition_distribution = (
        result.condition_distribution.loc[
            :, ["condition", "cluster", "count", "fraction", "replicate_n"]
        ]
        .sort_values(["condition", "cluster"], kind="stable")
        .reset_index(drop=True)
    )

    if result.distribution_change is None:
        distribution_change = _empty_distribution_change_table()
    else:
        _require_columns(
            result.distribution_change,
            (
                "condition",
                "cluster",
                "count",
                "fraction",
                "replicate_n",
                "reference_fraction",
                "delta_fraction",
                "log2_fc",
            ),
            "result.distribution_change",
        )
        distribution_change = (
            result.distribution_change.loc[
                :,
                [
                    "condition",
                    "cluster",
                    "count",
                    "fraction",
                    "replicate_n",
                    "reference_fraction",
                    "delta_fraction",
                    "log2_fc",
                ],
            ]
            .sort_values(["condition", "cluster"], kind="stable")
            .reset_index(drop=True)
        )

    metadata = {
        "mode": result.model.mode,
        "cluster_labels": list(result.model.cluster_labels),
        "has_distribution_change": not distribution_change.empty,
        "change_condition_order": distribution_change["condition"]
        .drop_duplicates()
        .tolist(),
    }
    return {
        "sample_distribution": sample_distribution,
        "condition_distribution": condition_distribution,
        "distribution_change": distribution_change,
        "metadata": metadata,
    }


def prepare_shared_cluster_profile_data(
    *,
    result,
    features: list[str] | None = None,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_shared_cluster_result(result)
    if not hasattr(result.model, "feature_names"):
        raise TypeError("plotting helpers require a SharedClusterResult-like object.")

    cluster_profiles = result.cluster_profiles.copy()
    feature_names = list(result.model.feature_names)

    requested_features = feature_names
    if features is not None:
        missing = [feature for feature in features if feature not in feature_names]
        if missing:
            raise ValueError(
                "Requested features are not present in cluster_profiles: "
                f"{', '.join(missing)}"
            )
        requested_features = [
            feature for feature in feature_names if feature in features
        ]

    _require_columns(
        cluster_profiles,
        ("cluster", "count", *feature_names),
        "result.cluster_profiles",
    )

    if cluster_profiles.empty:
        return {
            "profile_table": pd.DataFrame(
                columns=["cluster", "feature", "value", "count"]
            ),
            "metadata": {
                "feature_names": requested_features,
                "cluster_labels": list(result.model.cluster_labels),
            },
        }

    profile_table = (
        cluster_profiles.loc[:, ["cluster", "count", *requested_features]]
        .melt(id_vars=["cluster", "count"], var_name="feature", value_name="value")
        .sort_values(["cluster", "feature"], kind="stable")
        .reset_index(drop=True)
    )
    return {
        "profile_table": profile_table,
        "metadata": {
            "feature_names": requested_features,
            "cluster_labels": list(result.model.cluster_labels),
        },
    }


def prepare_shared_cluster_region_data(
    *,
    result,
    aggregate_conditions: bool = True,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_shared_cluster_result(result)
    if result.region_summaries is None:
        raise ValueError(
            "SharedClusterResult.region_summaries is required for region plotting."
        )

    _require_columns(
        result.region_summaries,
        ("region_id", "sample_id", "condition", "cluster", "count", "fraction"),
        "result.region_summaries",
    )
    region_table = (
        result.region_summaries.loc[
            :, ["region_id", "sample_id", "condition", "cluster", "count", "fraction"]
        ]
        .sort_values(["region_id", "sample_id", "cluster"], kind="stable")
        .reset_index(drop=True)
    )

    if aggregate_conditions and not region_table.empty:
        condition_region_table = (
            region_table.groupby(
                ["region_id", "condition", "cluster"], as_index=False, sort=False
            )
            .agg(
                count=("count", "sum"),
                fraction_mean=("fraction", "mean"),
                fraction_median=("fraction", "median"),
                sample_n=("sample_id", "nunique"),
            )
            .sort_values(["region_id", "condition", "cluster"], kind="stable")
            .reset_index(drop=True)
        )
    else:
        condition_region_table = pd.DataFrame(
            columns=[
                "region_id",
                "condition",
                "cluster",
                "count",
                "fraction_mean",
                "fraction_median",
                "sample_n",
            ]
        )

    return {
        "region_table": region_table,
        "condition_region_table": condition_region_table,
        "metadata": {
            "mode": result.model.mode,
            "cluster_labels": list(result.model.cluster_labels),
            "has_condition_aggregation": aggregate_conditions,
        },
    }


_READ_CLUSTER_ASSOCIATION_VALUE_MODES = {"fraction", "log2_enrichment"}
_READ_CLUSTER_ASSOCIATION_REGION_COLUMNS = (
    "region_id",
    "chrom",
    "chromosome",
    "start",
    "end",
    "strand",
)
_READ_CLUSTER_ASSOCIATION_LEGACY_EXCLUDE_COLUMNS = {
    "count",
    "fraction",
    "log2_enrichment",
    "total_reads",
    "dominant_cluster",
    "dominant_fraction",
    "entropy",
    "value",
    "value_mode",
}


def _format_region_coordinate(value) -> str:
    if pd.isna(value):
        raise ValueError("region coordinates must not be missing.")
    if isinstance(value, (np.integer, int)):
        return str(int(value))
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _synthesize_read_cluster_region_id(table: pd.DataFrame, *, owner: str) -> pd.Series:
    if "chrom" not in table.columns and "chromosome" not in table.columns:
        raise ValueError(
            f"{owner} requires region_id or chrom/chromosome, start, end, and strand columns."
        )
    _require_columns(table, ("start", "end", "strand"), owner)

    chrom_column = "chrom" if "chrom" in table.columns else "chromosome"
    if table[[chrom_column, "start", "end", "strand"]].isna().any().any():
        raise ValueError(
            f"{owner} requires region_id or chrom/chromosome, start, end, and strand columns."
        )

    return table.loc[:, [chrom_column, "start", "end", "strand"]].apply(
        lambda row: (
            f"{_format_region_coordinate(row.iloc[0])}:"
            f"{_format_region_coordinate(row.iloc[1])}-"
            f"{_format_region_coordinate(row.iloc[2])}:"
            f"{_format_region_coordinate(row.iloc[3])}"
        ),
        axis=1,
    )


def _ensure_read_cluster_region_ids(table: pd.DataFrame, *, owner: str) -> pd.DataFrame:
    normalized = table.copy()
    if "region_id" not in normalized.columns:
        normalized["region_id"] = _synthesize_read_cluster_region_id(
            normalized, owner=owner
        )
        return normalized

    missing_region_ids = normalized["region_id"].isna()
    if missing_region_ids.any():
        replacement_ids = _synthesize_read_cluster_region_id(
            normalized.loc[missing_region_ids], owner=owner
        )
        normalized.loc[missing_region_ids, "region_id"] = replacement_ids
    normalized["region_id"] = normalized["region_id"].astype(str)
    return normalized


def _normalize_long_form_read_cluster_region_association(
    table: pd.DataFrame,
    *,
    value_mode: str,
) -> pd.DataFrame:
    _require_columns(table, ("cluster",), "association_table")
    if value_mode not in _READ_CLUSTER_ASSOCIATION_VALUE_MODES:
        raise ValueError("Unsupported read-cluster association value_mode.")
    if value_mode not in table.columns:
        raise ValueError(
            f"association_table does not include the requested {value_mode} column."
        )

    normalized = _ensure_read_cluster_region_ids(table, owner="association_table")
    if normalized["region_id"].isna().any():
        raise ValueError(
            "association_table requires region_id or chrom/chromosome, start, end, and strand columns."
        )

    normalized = normalized.copy()
    normalized["value"] = normalized[value_mode]
    normalized["value_mode"] = value_mode
    normalized["source_format"] = "long_form"
    return normalized


def _infer_read_cluster_columns(table: pd.DataFrame) -> list[object]:
    excluded_columns = {
        column
        for column in _READ_CLUSTER_ASSOCIATION_REGION_COLUMNS
        if column in table.columns
    }
    excluded_columns.update(
        col
        for col in _READ_CLUSTER_ASSOCIATION_LEGACY_EXCLUDE_COLUMNS
        if col in table.columns
    )
    cluster_columns = [
        column for column in table.columns if column not in excluded_columns
    ]
    if "region_id" in cluster_columns:
        cluster_columns.remove("region_id")
    return cluster_columns


def _normalize_legacy_read_cluster_region_association(
    table: pd.DataFrame,
    *,
    value_mode: str,
) -> pd.DataFrame:
    if value_mode != "fraction":
        raise ValueError(
            "Legacy read-cluster association tables only support value_mode='fraction'."
        )

    cluster_columns = _infer_read_cluster_columns(table)
    if not cluster_columns:
        raise ValueError(
            "Could not infer cluster columns from the legacy association table."
        )

    normalized = _ensure_read_cluster_region_ids(table, owner="association_table")
    if normalized["region_id"].isna().any():
        raise ValueError(
            "association_table requires region_id or chrom/chromosome, start, end, and strand columns."
        )

    if "total_reads" not in normalized.columns:
        normalized["total_reads"] = normalized.loc[:, cluster_columns].sum(axis=1)
    if (normalized["total_reads"] <= 0).any():
        raise ValueError(
            "Legacy association rows must have positive total_reads or cluster counts."
        )

    id_vars = [column for column in normalized.columns if column not in cluster_columns]
    melted = normalized.melt(
        id_vars=id_vars,
        value_vars=cluster_columns,
        var_name="cluster",
        value_name="count",
    )
    melted["fraction"] = melted["count"] / melted["total_reads"]
    melted["value"] = melted["fraction"]
    melted["value_mode"] = "fraction"
    melted["source_format"] = "legacy_wide"
    return melted


def _read_cluster_region_association_cluster_order(table: pd.DataFrame) -> list[object]:
    return table.loc[:, "cluster"].drop_duplicates().tolist()


def _parse_region_id_components(region_id: object) -> tuple[str, int, int, str]:
    text = str(region_id)
    strand = "."
    core = text
    if "," in text:
        candidate_core, candidate_strand = text.rsplit(",", 1)
        if candidate_strand in {"+", "-", "."}:
            core = candidate_core
            strand = candidate_strand
    if ":" not in core or "-" not in core:
        return text, -1, -1, strand
    chrom, coords = core.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    try:
        start = int(start_s)
        end = int(end_s)
    except ValueError:
        start = -1
        end = -1
    return chrom, start, end, strand


def _region_coordinate_table_for_association(table: pd.DataFrame) -> pd.DataFrame:
    regions = table.loc[:, ["region_id"]].drop_duplicates().copy()
    has_chrom = "chrom" in table.columns or "chromosome" in table.columns
    has_start_end = "start" in table.columns and "end" in table.columns
    if has_chrom and has_start_end:
        chrom_col = "chrom" if "chrom" in table.columns else "chromosome"
        coords = (
            table.loc[
                :,
                ["region_id", chrom_col, "start", "end"]
                + (["strand"] if "strand" in table.columns else []),
            ]
            .drop_duplicates(subset=["region_id"])
            .copy()
        )
        coords = coords.rename(columns={chrom_col: "chrom"})
        if "strand" not in coords.columns:
            coords["strand"] = "."
        coords["start"] = (
            pd.to_numeric(coords["start"], errors="coerce").fillna(-1).astype(int)
        )
        coords["end"] = (
            pd.to_numeric(coords["end"], errors="coerce").fillna(-1).astype(int)
        )
        return coords.loc[:, ["region_id", "chrom", "start", "end", "strand"]]

    parsed = regions["region_id"].map(_parse_region_id_components)
    regions["chrom"] = parsed.map(lambda x: x[0])
    regions["start"] = parsed.map(lambda x: x[1])
    regions["end"] = parsed.map(lambda x: x[2])
    regions["strand"] = parsed.map(lambda x: x[3])
    return regions


def _chromosome_sort_key(chrom: object) -> tuple[int, object]:
    text = str(chrom)
    if text.lower().startswith("chr"):
        suffix = text[3:]
        if suffix.isdigit():
            return (0, int(suffix))
        if suffix in {"X", "x"}:
            return (1, 23)
        if suffix in {"Y", "y"}:
            return (1, 24)
        if suffix in {"M", "MT", "m", "mt"}:
            return (1, 25)
        return (2, suffix)
    return (3, text)


def _read_cluster_region_association_region_order(
    table: pd.DataFrame,
    *,
    region_sort: str,
    strength_aggregate: str,
) -> list[object]:
    if region_sort == "input":
        return table.loc[:, "region_id"].drop_duplicates().tolist()

    coords = _region_coordinate_table_for_association(table)
    coords = coords.copy()
    coords["chrom_sort"] = coords["chrom"].map(_chromosome_sort_key)
    coords = coords.sort_values(
        ["chrom_sort", "start", "end", "region_id"], kind="stable"
    )

    if region_sort == "genomic":
        return coords["region_id"].tolist()

    if region_sort == "association_strength":
        grouped = table.groupby("region_id", sort=False)["value"]
        strength = grouped.mean() if strength_aggregate == "mean" else grouped.max()
        coords = coords.merge(
            strength.rename("association_strength"),
            left_on="region_id",
            right_index=True,
            how="left",
        )
        coords = coords.sort_values(
            ["association_strength", "chrom_sort", "start", "end", "region_id"],
            ascending=[False, True, True, True, True],
            kind="stable",
        )
        return coords["region_id"].tolist()

    if region_sort == "cluster_fraction":
        # Assign each region to its dominant cluster, then rank by dominant fraction.
        cluster_order = _read_cluster_region_association_cluster_order(table)
        cluster_rank = {cluster: index for index, cluster in enumerate(cluster_order)}
        dominant_rows = table.copy()
        dominant_rows["cluster_rank"] = (
            dominant_rows["cluster"]
            .map(cluster_rank)
            .fillna(len(cluster_order))
            .astype(int)
        )
        dominant_rows["value"] = pd.to_numeric(
            dominant_rows["value"], errors="coerce"
        ).fillna(0.0)
        dominant_rows = dominant_rows.sort_values(
            ["region_id", "value", "cluster_rank"],
            ascending=[True, False, True],
            kind="stable",
        )
        dominant_rows = dominant_rows.drop_duplicates(
            subset=["region_id"], keep="first"
        )
        dominant_rows = dominant_rows.rename(
            columns={"cluster": "dominant_cluster", "value": "dominant_sort_value"}
        ).loc[
            :, ["region_id", "dominant_cluster", "dominant_sort_value", "cluster_rank"]
        ]

        coords = coords.merge(dominant_rows, on="region_id", how="left")
        coords["cluster_rank"] = (
            pd.to_numeric(coords["cluster_rank"], errors="coerce")
            .fillna(len(cluster_order))
            .astype(int)
        )
        coords["dominant_sort_value"] = pd.to_numeric(
            coords["dominant_sort_value"], errors="coerce"
        ).fillna(0.0)
        coords = coords.sort_values(
            [
                "cluster_rank",
                "dominant_sort_value",
                "chrom_sort",
                "start",
                "end",
                "region_id",
            ],
            ascending=[True, False, True, True, True, True],
            kind="stable",
        )
        return coords["region_id"].tolist()

    raise ValueError(
        "region_sort must be 'cluster_fraction', 'input', 'genomic', or 'association_strength'."
    )


def prepare_read_cluster_region_association_data(
    association_table: pd.DataFrame,
    *,
    value_mode: str = "fraction",
    top_n_regions_per_cluster: int | None = 5,
    region_sort: str = "cluster_fraction",
    association_strength_aggregate: str = "max",
) -> dict[str, pd.DataFrame | dict[str, object]]:
    if not isinstance(association_table, pd.DataFrame):
        raise TypeError("association_table must be a pandas DataFrame.")
    if value_mode not in _READ_CLUSTER_ASSOCIATION_VALUE_MODES:
        raise ValueError("value_mode must be 'fraction' or 'log2_enrichment'.")
    if top_n_regions_per_cluster is not None and top_n_regions_per_cluster < 0:
        raise ValueError("top_n_regions_per_cluster must be non-negative.")
    if region_sort not in {
        "cluster_fraction",
        "input",
        "genomic",
        "association_strength",
    }:
        raise ValueError(
            "region_sort must be 'cluster_fraction', 'input', 'genomic', or 'association_strength'."
        )
    if association_strength_aggregate not in {"max", "mean"}:
        raise ValueError("association_strength_aggregate must be 'max' or 'mean'.")

    if "cluster" in association_table.columns:
        normalized = _normalize_long_form_read_cluster_region_association(
            association_table,
            value_mode=value_mode,
        )
    else:
        normalized = _normalize_legacy_read_cluster_region_association(
            association_table,
            value_mode=value_mode,
        )

    normalized = normalized.copy()
    if normalized.duplicated(["region_id", "cluster"]).any():
        raise ValueError(
            "association_table contains duplicate region and cluster rows."
        )

    region_order = _read_cluster_region_association_region_order(
        normalized,
        region_sort=region_sort,
        strength_aggregate=association_strength_aggregate,
    )
    cluster_order = _read_cluster_region_association_cluster_order(normalized)
    normalized["region_id"] = pd.Categorical(
        normalized["region_id"], categories=region_order, ordered=True
    )
    normalized["cluster"] = pd.Categorical(
        normalized["cluster"], categories=cluster_order, ordered=True
    )
    normalized = normalized.sort_values(
        ["region_id", "cluster"], kind="stable"
    ).reset_index(drop=True)
    normalized["region_id"] = normalized["region_id"].astype(str)
    normalized["cluster"] = normalized["cluster"].astype(object)

    matrix_table = (
        normalized.loc[:, ["region_id", "cluster", "value"]]
        .pivot(index="region_id", columns="cluster", values="value")
        .reindex(index=region_order, columns=cluster_order)
        .reset_index()
    )
    matrix_table.columns.name = None
    region_axis_table = _region_coordinate_table_for_association(normalized).copy()
    region_axis_table = (
        region_axis_table.set_index("region_id").reindex(region_order).reset_index()
    )
    region_axis_table["axis_index"] = np.arange(len(region_axis_table))
    chromosome_blocks: list[dict[str, object]] = []
    if not region_axis_table.empty:
        grouped = region_axis_table.groupby("chrom", sort=False)["axis_index"]
        for chrom, values in grouped:
            chromosome_blocks.append(
                {
                    "chrom": str(chrom),
                    "start_index": int(values.min()),
                    "end_index": int(values.max()),
                    "n_regions": int(values.max() - values.min() + 1),
                }
            )

    top_columns = [
        column
        for column in [
            "region_id",
            "chrom",
            "chromosome",
            "start",
            "end",
            "strand",
            "cluster",
            "count",
            "fraction",
            "log2_enrichment",
            "total_reads",
            "value",
        ]
        if column in normalized.columns
    ]
    if top_n_regions_per_cluster is None or top_n_regions_per_cluster == 0:
        top_regions_table = normalized.loc[:, top_columns].head(0).copy()
        top_regions_table["rank"] = pd.Series(dtype="int64")
    else:
        top_regions_table = (
            normalized.loc[:, top_columns]
            .sort_values(
                ["cluster", "value", "region_id"],
                ascending=[True, False, True],
                kind="stable",
            )
            .groupby("cluster", as_index=False, sort=False, group_keys=False)
            .head(top_n_regions_per_cluster)
            .copy()
        )
        top_regions_table["rank"] = (
            top_regions_table.groupby("cluster", sort=False).cumcount() + 1
        )
        top_regions_table = top_regions_table.sort_values(
            ["cluster", "rank"], kind="stable"
        ).reset_index(drop=True)

    metadata = {
        "plot_family": "read_cluster_region_association_heatmap",
        "value_mode": value_mode,
        "source_format": normalized["source_format"].iloc[0]
        if not normalized.empty
        else "long_form",
        "region_order": [str(value) for value in region_order],
        "cluster_order": cluster_order,
        "top_n_regions_per_cluster": top_n_regions_per_cluster,
        "has_top_regions_table": not top_regions_table.empty,
        "region_sort": region_sort,
        "association_strength_aggregate": association_strength_aggregate,
        "chromosome_blocks": chromosome_blocks,
    }
    return {
        "association_table": normalized.reset_index(drop=True),
        "matrix_table": matrix_table,
        "region_axis_table": region_axis_table,
        "top_regions_table": top_regions_table,
        "metadata": metadata,
    }
