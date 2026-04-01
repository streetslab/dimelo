from __future__ import annotations

from typing import Iterable

import pandas as pd

from . import global_analysis
from .models import ContrastSpec, RegionDiscoveryResult
from .region_contrasts import (
    _adjust_p_values_bh,
    _beta_binomial_two_sided_p_value,
    _estimate_beta_binomial_prior,
)

_WINDOW_KEY_COLUMNS = ["window_id", "chromosome", "start", "end", "strand"]


def _validate_motifs(motifs: Iterable[str]) -> list[str]:
    motif_list = list(motifs)
    if len(motif_list) != 1:
        raise ValueError("scan_genome currently supports exactly one motif.")
    return motif_list


def _safe_fraction(modified_count: pd.Series, valid_count: pd.Series) -> pd.Series:
    return modified_count.div(valid_count.where(valid_count != 0), fill_value=0).fillna(0.0)


def _aggregate_window_counts(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame(columns=_WINDOW_KEY_COLUMNS + ["modified_count", "valid_count", "window_fraction"])

    aggregated = (
        summary.groupby(_WINDOW_KEY_COLUMNS, as_index=False, sort=False)
        .agg(
            modified_count=("modified_count", "sum"),
            valid_count=("valid_count", "sum"),
        )
        .copy()
    )
    aggregated["window_fraction"] = _safe_fraction(
        aggregated["modified_count"], aggregated["valid_count"]
    )
    return aggregated


def _aggregate_condition_counts(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame(
            columns=_WINDOW_KEY_COLUMNS + ["condition", "modified_count", "valid_count", "window_fraction"]
        )

    aggregated = (
        summary.groupby(_WINDOW_KEY_COLUMNS + ["condition"], as_index=False, sort=False)
        .agg(
            modified_count=("modified_count", "sum"),
            valid_count=("valid_count", "sum"),
        )
        .copy()
    )
    aggregated["window_fraction"] = _safe_fraction(
        aggregated["modified_count"], aggregated["valid_count"]
    )
    return aggregated


def _condition_spread_scores(condition_counts: pd.DataFrame) -> pd.DataFrame:
    if condition_counts.empty:
        return pd.DataFrame(columns=_WINDOW_KEY_COLUMNS + ["score_value"])

    spread = (
        condition_counts.groupby(_WINDOW_KEY_COLUMNS, as_index=False, sort=False)
        .agg(
            score_value=("window_fraction", lambda values: float(values.max() - values.min())),
        )
        .copy()
    )
    return spread


def _side_counts(
    condition_counts: pd.DataFrame,
    *,
    conditions: Iterable[str],
    side_name: str,
) -> pd.DataFrame:
    side_conditions = list(conditions)
    if not side_conditions:
        return pd.DataFrame(columns=_WINDOW_KEY_COLUMNS + [f"{side_name}_modified_count", f"{side_name}_valid_count"])

    side = condition_counts.loc[condition_counts["condition"].isin(side_conditions)].copy()
    if side.empty:
        return pd.DataFrame(columns=_WINDOW_KEY_COLUMNS + [f"{side_name}_modified_count", f"{side_name}_valid_count"])

    grouped = (
        side.groupby(_WINDOW_KEY_COLUMNS, as_index=False, sort=False)
        .agg(
            **{
                f"{side_name}_modified_count": ("modified_count", "sum"),
                f"{side_name}_valid_count": ("valid_count", "sum"),
            }
        )
        .copy()
    )
    return grouped


def _score_with_contrast(
    window_totals: pd.DataFrame,
    condition_counts: pd.DataFrame,
    *,
    contrast: ContrastSpec | None,
    score: str,
) -> pd.DataFrame:
    if score == "beta_binomial":
        available_conditions = sorted(condition_counts["condition"].dropna().unique().tolist())
        if contrast is not None and contrast.numerator and contrast.denominator:
            numerator_conditions = list(contrast.numerator)
            denominator_conditions = list(contrast.denominator)
        elif len(available_conditions) >= 2:
            numerator_conditions = [available_conditions[0]]
            denominator_conditions = [available_conditions[1]]
        else:
            numerator_conditions = []
            denominator_conditions = []

        if not numerator_conditions and not denominator_conditions:
            scored = window_totals.copy()
            scored["score_value"] = 0.0
            scored["p_value"] = 1.0
            scored["adjusted_p_value"] = 1.0
            return scored

        numerator_counts = _side_counts(
            condition_counts, conditions=numerator_conditions, side_name="numerator"
        )
        denominator_counts = _side_counts(
            condition_counts, conditions=denominator_conditions, side_name="denominator"
        )

        scored = window_totals.merge(numerator_counts, on=_WINDOW_KEY_COLUMNS, how="left")
        scored = scored.merge(denominator_counts, on=_WINDOW_KEY_COLUMNS, how="left")

        for column in [
            "numerator_modified_count",
            "numerator_valid_count",
            "denominator_modified_count",
            "denominator_valid_count",
        ]:
            scored[column] = scored[column].fillna(0).astype(int)

        scored["score_value"] = _safe_fraction(
            scored["numerator_modified_count"], scored["numerator_valid_count"]
        ) - _safe_fraction(scored["denominator_modified_count"], scored["denominator_valid_count"])

        scored["p_value"] = [
            _beta_binomial_two_sided_p_value(
                int(numerator_modified_count),
                int(numerator_valid_count),
                *_estimate_beta_binomial_prior(
                    int(denominator_modified_count),
                    int(denominator_valid_count),
                ),
            )
            for numerator_modified_count, numerator_valid_count, denominator_modified_count, denominator_valid_count in zip(
                scored["numerator_modified_count"],
                scored["numerator_valid_count"],
                scored["denominator_modified_count"],
                scored["denominator_valid_count"],
            )
        ]
        scored["adjusted_p_value"] = _adjust_p_values_bh(scored["p_value"])
        return scored

    if contrast is None or not contrast.numerator or not contrast.denominator:
        spread = _condition_spread_scores(condition_counts)
        scored = window_totals.merge(spread, on=_WINDOW_KEY_COLUMNS, how="left")
        scored["score_value"] = scored["score_value"].fillna(0.0)
        scored["p_value"] = pd.NA
        scored["adjusted_p_value"] = pd.NA
        return scored

    numerator_counts = _side_counts(
        condition_counts, conditions=contrast.numerator, side_name="numerator"
    )
    denominator_counts = _side_counts(
        condition_counts, conditions=contrast.denominator, side_name="denominator"
    )

    scored = window_totals.merge(numerator_counts, on=_WINDOW_KEY_COLUMNS, how="left")
    scored = scored.merge(denominator_counts, on=_WINDOW_KEY_COLUMNS, how="left")

    for column in [
        "numerator_modified_count",
        "numerator_valid_count",
        "denominator_modified_count",
        "denominator_valid_count",
    ]:
        scored[column] = scored[column].fillna(0).astype(int)

    scored["score_value"] = _safe_fraction(
        scored["numerator_modified_count"], scored["numerator_valid_count"]
    ) - _safe_fraction(scored["denominator_modified_count"], scored["denominator_valid_count"])

    if score == "beta_binomial":
        scored["p_value"] = [
            _beta_binomial_two_sided_p_value(
                int(numerator_modified_count),
                int(numerator_valid_count),
                *_estimate_beta_binomial_prior(
                    int(denominator_modified_count),
                    int(denominator_valid_count),
                ),
            )
            for numerator_modified_count, numerator_valid_count, denominator_modified_count, denominator_valid_count in zip(
                scored["numerator_modified_count"],
                scored["numerator_valid_count"],
                scored["denominator_modified_count"],
                scored["denominator_valid_count"],
            )
        ]
    else:
        scored["p_value"] = pd.NA
        scored["adjusted_p_value"] = pd.NA

    return scored


def _rank_windows(scored: pd.DataFrame, *, score: str) -> pd.DataFrame:
    if scored.empty:
        return scored.assign(rank=pd.Series(dtype="float64"))

    if score == "beta_binomial":
        ranked = scored.sort_values(
            by=["adjusted_p_value", "p_value", "score_value", "chromosome", "start", "end"],
            ascending=[True, True, False, True, True, True],
            kind="mergesort",
            na_position="last",
        ).reset_index(drop=True)
    else:
        ranked = scored.sort_values(
            by=["score_value", "chromosome", "start", "end"],
            ascending=[False, True, True, True],
            kind="mergesort",
            na_position="last",
        ).reset_index(drop=True)

    ranked["rank"] = range(1, len(ranked) + 1)
    return ranked


def _apply_min_coverage(
    scored: pd.DataFrame,
    *,
    min_coverage: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if scored.empty:
        return scored.copy(), scored.copy()

    filtered = scored.copy()
    covered = filtered["valid_count"] >= min_coverage
    filtered.loc[~covered, ["score_value", "p_value", "adjusted_p_value", "rank"]] = pd.NA
    hits = filtered.loc[covered].copy()
    return filtered, hits


def _merge_adjacent_hits(hits: pd.DataFrame, *, merge_distance: int) -> pd.DataFrame:
    if hits.empty or len(hits) == 1:
        return hits.copy()

    ordered = hits.sort_values(
        by=["chromosome", "start", "end", "rank"],
        ascending=[True, True, True, True],
        kind="mergesort",
    ).reset_index(drop=True)

    merged_rows: list[dict[str, object]] = []
    current = ordered.iloc[0].to_dict()
    current["merged_window_count"] = 1

    def _as_float(value: object) -> float | None:
        return None if pd.isna(value) else float(value)

    for _, row in ordered.iloc[1:].iterrows():
        same_chromosome = row["chromosome"] == current["chromosome"]
        within_distance = row["start"] <= int(current["end"]) + merge_distance
        if same_chromosome and within_distance:
            current["end"] = max(int(current["end"]), int(row["end"]))
            current["modified_count"] = int(current["modified_count"]) + int(row["modified_count"])
            current["valid_count"] = int(current["valid_count"]) + int(row["valid_count"])
            current_score = _as_float(current.get("score_value"))
            row_score = _as_float(row.get("score_value"))
            if current_score is None:
                current["score_value"] = row_score
            elif row_score is not None:
                current["score_value"] = max(current_score, row_score)
            current_p_value = _as_float(current.get("p_value"))
            row_p_value = _as_float(row.get("p_value"))
            if current_p_value is None:
                current["p_value"] = row_p_value
            elif row_p_value is not None:
                current["p_value"] = min(current_p_value, row_p_value)

            current_adjusted = _as_float(current.get("adjusted_p_value"))
            row_adjusted = _as_float(row.get("adjusted_p_value"))
            if current_adjusted is None:
                current["adjusted_p_value"] = row_adjusted
            elif row_adjusted is not None:
                current["adjusted_p_value"] = min(current_adjusted, row_adjusted)
            current["rank"] = min(int(current["rank"]), int(row["rank"]))
            current["merged_window_count"] += 1
            continue

        merged_rows.append(current)
        current = row.to_dict()
        current["merged_window_count"] = 1

    merged_rows.append(current)
    merged = pd.DataFrame(merged_rows)
    merged["rank"] = range(1, len(merged) + 1)
    return merged


def scan_genome(
    *,
    samples,
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    include_contigs: Iterable[str] | None = None,
    exclude_contigs: Iterable[str] | None = None,
    min_coverage: int = 0,
    merge_hits: bool = False,
    merge_distance: int = 0,
    score: str = "effect_size_only",
    contrast: ContrastSpec | None = None,
    quiet: bool = True,
    cores: int | None = None,
) -> RegionDiscoveryResult:
    motif_list = _validate_motifs(motifs)
    if score not in {"effect_size_only", "beta_binomial"}:
        raise ValueError("scan_genome requires score in {'effect_size_only', 'beta_binomial'}.")

    window_summary = global_analysis.build_window_summary(
        samples=samples,
        motifs=motif_list,
        genome_sizes=genome_sizes,
        window_size=window_size,
        step_size=step_size,
        include_contigs=include_contigs,
        exclude_contigs=exclude_contigs,
        quiet=quiet,
        cores=cores,
    )

    window_totals = _aggregate_window_counts(window_summary)
    condition_counts = _aggregate_condition_counts(window_summary)
    scored = _score_with_contrast(
        window_totals=window_totals,
        condition_counts=condition_counts,
        contrast=contrast,
        score=score,
    )
    ranked = _rank_windows(scored, score=score)
    window_table, hits = _apply_min_coverage(ranked, min_coverage=min_coverage)

    if merge_hits:
        hits = _merge_adjacent_hits(hits, merge_distance=merge_distance)

    plot_data = {
        "window_score_table": window_table.copy(),
        "top_hits_table": hits.copy(),
    }
    metadata = {
        "analysis_unit": "de_novo_window",
        "representation": "modified_fraction",
        "signal_source": "pileup_counts",
        "score": score,
        "window_size": window_size,
        "step_size": step_size,
        "min_coverage": min_coverage,
        "merge_hits": merge_hits,
        "merge_distance": merge_distance,
        "motifs": motif_list,
        "include_contigs": list(include_contigs) if include_contigs is not None else None,
        "exclude_contigs": list(exclude_contigs) if exclude_contigs is not None else None,
    }
    if contrast is not None:
        metadata["contrast_mode"] = contrast.mode
        metadata["contrast_numerator"] = list(contrast.numerator or [])
        metadata["contrast_denominator"] = list(contrast.denominator or [])

    return RegionDiscoveryResult(
        hits=hits.reset_index(drop=True),
        windows=window_table.reset_index(drop=True),
        contrast=contrast,
        plot_data=plot_data,
        metadata=metadata,
        figures={},
    )
