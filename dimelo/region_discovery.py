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
_MERGE_SORT_COLUMNS = ["chromosome", "strand", "start", "end", "rank"]
_OUTPUT_SORT_COLUMNS = ["rank", "chromosome", "strand", "start", "end"]
_BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


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


def _validate_contrast_conditions(
    *,
    condition_counts: pd.DataFrame,
    contrast: ContrastSpec | None,
) -> None:
    if contrast is None:
        return

    available_conditions = set(condition_counts["condition"].dropna().tolist())
    missing: list[str] = []
    for side_name, conditions in (
        ("numerator", contrast.numerator or []),
        ("denominator", contrast.denominator or []),
    ):
        absent = sorted(set(conditions) - available_conditions)
        if absent:
            missing.append(f"{side_name}: {', '.join(absent)}")

    if missing:
        raise ValueError(
            "scan_genome contrast requested missing condition(s): "
            + "; ".join(missing)
        )


def _pairing_policy_value(pairing_policy: str | None) -> str:
    return pairing_policy or "complete_pairs_only"


def _is_paired_contrast(contrast: ContrastSpec | None) -> bool:
    return contrast is not None and contrast.mode in {"matched_pairwise", "time_course"}


def _resolve_pair_ids(samples, pairing_key: str) -> dict[str, object]:
    if not pairing_key:
        raise ValueError("scan_genome paired discovery requires an explicit pairing_key.")

    pair_ids: dict[str, object] = {}
    for sample in samples:
        metadata = sample.metadata or {}
        if pairing_key not in metadata:
            raise ValueError(
                f"scan_genome paired discovery requires sample.metadata['{pairing_key}'] for every sample."
            )
        pair_ids[sample.sample_id] = metadata[pairing_key]
    return pair_ids


def _build_paired_window_table(
    summary: pd.DataFrame,
    *,
    samples,
    pairing_key: str | None,
    required_conditions: list[str],
    pairing_policy: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    sample_to_pair = _resolve_pair_ids(samples, pairing_key or "")

    paired = summary.copy()
    paired["pair_id"] = paired["sample_id"].map(sample_to_pair)
    paired = paired.dropna(subset=["pair_id"])

    pair_conditions = (
        paired.loc[:, ["pair_id", "condition"]]
        .drop_duplicates()
        .groupby("pair_id")["condition"]
        .agg(lambda values: set(values))
    )
    required_condition_set = set(required_conditions)
    complete_pair_ids = [
        pair_id for pair_id, conditions in pair_conditions.items() if required_condition_set.issubset(conditions)
    ]
    dropped_pair_count = int(len(pair_conditions) - len(complete_pair_ids))

    if pairing_policy == "error_on_missing" and dropped_pair_count:
        raise ValueError("scan_genome paired discovery found incomplete matched units.")
    if not complete_pair_ids:
        raise ValueError("scan_genome paired discovery found no complete matched units.")

    paired = paired.loc[paired["pair_id"].isin(complete_pair_ids)].copy()
    paired = (
        paired.groupby(_WINDOW_KEY_COLUMNS + ["pair_id", "condition"], as_index=False, sort=False)
        .agg(
            modified_count=("modified_count", "sum"),
            valid_count=("valid_count", "sum"),
        )
        .copy()
    )
    paired["window_fraction"] = _safe_fraction(paired["modified_count"], paired["valid_count"])
    return paired, {
        "n_pairs_used": len(complete_pair_ids),
        "n_pairs_dropped": dropped_pair_count,
    }


def _score_matched_pairwise(
    paired_window_table: pd.DataFrame,
    *,
    contrast: ContrastSpec,
    rank_by: str = "mean_abs_delta",
) -> pd.DataFrame:
    if rank_by != "mean_abs_delta":
        raise ValueError("scan_genome matched_pairwise currently supports rank_by='mean_abs_delta'.")

    numerator = paired_window_table.loc[
        paired_window_table["condition"].isin(contrast.numerator or [])
    ].copy()
    denominator = paired_window_table.loc[
        paired_window_table["condition"].isin(contrast.denominator or [])
    ].copy()

    if numerator.empty or denominator.empty:
        scored = pd.DataFrame(columns=_WINDOW_KEY_COLUMNS)
        scored["mean_delta"] = pd.Series(dtype="float64")
        scored["mean_abs_delta"] = pd.Series(dtype="float64")
        scored["delta_sd"] = pd.Series(dtype="float64")
        scored["sign_agreement"] = pd.Series(dtype="float64")
        scored["n_pairs_used"] = pd.Series(dtype="int64")
        scored["score_value"] = pd.Series(dtype="float64")
        scored["p_value"] = pd.Series(dtype="object")
        scored["adjusted_p_value"] = pd.Series(dtype="object")
        return scored

    numerator = numerator.rename(
        columns={
            "modified_count": "numerator_modified_count",
            "valid_count": "numerator_valid_count",
            "window_fraction": "numerator_fraction",
        }
    )
    denominator = denominator.rename(
        columns={
            "modified_count": "denominator_modified_count",
            "valid_count": "denominator_valid_count",
            "window_fraction": "denominator_fraction",
        }
    )

    merged = numerator.merge(
        denominator[
            _WINDOW_KEY_COLUMNS
            + ["pair_id", "denominator_modified_count", "denominator_valid_count", "denominator_fraction"]
        ],
        on=_WINDOW_KEY_COLUMNS + ["pair_id"],
        how="inner",
    )
    if merged.empty:
        scored = pd.DataFrame(columns=_WINDOW_KEY_COLUMNS)
        scored["mean_delta"] = pd.Series(dtype="float64")
        scored["mean_abs_delta"] = pd.Series(dtype="float64")
        scored["delta_sd"] = pd.Series(dtype="float64")
        scored["sign_agreement"] = pd.Series(dtype="float64")
        scored["n_pairs_used"] = pd.Series(dtype="int64")
        scored["score_value"] = pd.Series(dtype="float64")
        scored["p_value"] = pd.Series(dtype="object")
        scored["adjusted_p_value"] = pd.Series(dtype="object")
        return scored

    merged["delta"] = merged["numerator_fraction"] - merged["denominator_fraction"]

    scored = (
        merged.groupby(_WINDOW_KEY_COLUMNS, as_index=False, sort=False)
        .agg(
            mean_delta=("delta", "mean"),
            mean_abs_delta=("delta", lambda values: float(values.abs().mean())),
            delta_sd=("delta", lambda values: float(values.std(ddof=0))),
            sign_agreement=(
                "delta",
                lambda values: float(
                    max((values.gt(0)).mean(), (values.lt(0)).mean()) if len(values) else 0.0
                ),
            ),
            n_pairs_used=("pair_id", "nunique"),
        )
        .copy()
    )
    scored["score_value"] = scored["mean_abs_delta"]
    scored["p_value"] = pd.NA
    scored["adjusted_p_value"] = pd.NA
    return scored


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
        ).sub(
            _safe_fraction(scored["denominator_modified_count"], scored["denominator_valid_count"])
        ).abs()

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
    ).sub(
        _safe_fraction(scored["denominator_modified_count"], scored["denominator_valid_count"])
    ).abs()

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


def _sort_hits_for_merge(hits: pd.DataFrame) -> pd.DataFrame:
    sort_columns = [column for column in _MERGE_SORT_COLUMNS if column in hits.columns]
    if not sort_columns:
        return hits.copy()

    ordered = hits.sort_values(
        by=sort_columns,
        ascending=[True] * len(sort_columns),
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)
    return ordered


def _sort_hits_for_output(hits: pd.DataFrame) -> pd.DataFrame:
    sort_columns = [column for column in _OUTPUT_SORT_COLUMNS if column in hits.columns]
    if not sort_columns:
        return hits.copy()

    ascending = [True] * len(sort_columns)
    ordered = hits.sort_values(
        by=sort_columns,
        ascending=ascending,
        kind="mergesort",
        na_position="last",
    ).reset_index(drop=True)
    return ordered


def _merge_value(values: pd.Series, *, prefer: str) -> object:
    non_null = values.dropna()
    if non_null.empty:
        return pd.NA
    if prefer == "max":
        return non_null.max()
    if prefer == "min":
        return non_null.min()
    raise ValueError(f"Unsupported merge preference: {prefer}")


def _build_merged_hit(group: pd.DataFrame) -> dict[str, object]:
    merged = group.iloc[0].to_dict()

    if {"chromosome", "start", "end"}.issubset(group.columns):
        merged["chromosome"] = group.iloc[0]["chromosome"]
        merged["start"] = int(group["start"].min())
        merged["end"] = int(group["end"].max())

    if "modified_count" in group.columns:
        merged["modified_count"] = int(group["modified_count"].fillna(0).sum())
    if "valid_count" in group.columns:
        merged["valid_count"] = int(group["valid_count"].fillna(0).sum())
    if {"modified_count", "valid_count"}.issubset(merged):
        valid_count = int(merged.get("valid_count", 0))
        modified_count = int(merged.get("modified_count", 0))
        merged["window_fraction"] = 0.0 if valid_count == 0 else modified_count / valid_count

    if "score_value" in group.columns:
        merged["score_value"] = _merge_value(group["score_value"], prefer="max")
    if "p_value" in group.columns:
        merged["p_value"] = _merge_value(group["p_value"], prefer="min")
    if "adjusted_p_value" in group.columns:
        merged["adjusted_p_value"] = _merge_value(group["adjusted_p_value"], prefer="min")
    if "rank" in group.columns:
        rank_values = group["rank"].dropna()
        if not rank_values.empty:
            merged["rank"] = int(rank_values.min())

    contrast_count_fields = [
        "numerator_modified_count",
        "numerator_valid_count",
        "denominator_modified_count",
        "denominator_valid_count",
    ]
    contrast_counts_present = any(
        field in group.columns and group[field].notna().any() for field in contrast_count_fields
    )
    if contrast_counts_present:
        for field in contrast_count_fields:
            if field in group.columns and group[field].notna().any():
                merged[field] = int(group[field].fillna(0).sum())

    merged_window_count = 0
    if "merged_window_count" in group.columns:
        merged_window_count = int(group["merged_window_count"].fillna(1).sum())
    else:
        merged_window_count = len(group)
    merged["merged_window_count"] = merged_window_count

    if contrast_counts_present and all(
        field in merged and not pd.isna(merged[field]) for field in contrast_count_fields
    ):
        numerator_fraction = _safe_fraction(
            pd.Series([merged["numerator_modified_count"]], dtype="float64"),
            pd.Series([merged["numerator_valid_count"]], dtype="float64"),
        ).iloc[0]
        denominator_fraction = _safe_fraction(
            pd.Series([merged["denominator_modified_count"]], dtype="float64"),
            pd.Series([merged["denominator_valid_count"]], dtype="float64"),
        ).iloc[0]
        merged["score_value"] = abs(float(numerator_fraction) - float(denominator_fraction))

    if {"chromosome", "start", "end"}.issubset(merged):
        merged["window_id"] = f"{merged['chromosome']}:{int(merged['start'])}-{int(merged['end'])}"

    return merged


def merge_adjacent_hits(hits: pd.DataFrame, merge_distance: int) -> pd.DataFrame:
    if hits.empty:
        merged = hits.copy()
        if "merged_window_count" not in merged.columns:
            merged["merged_window_count"] = pd.Series(dtype="int64")
        return merged

    ordered = _sort_hits_for_merge(hits)
    if len(ordered) == 1:
        merged = ordered.copy()
        merged["merged_window_count"] = (
            merged["merged_window_count"]
            if "merged_window_count" in merged.columns
            else 1
        )
        if {"chromosome", "start", "end"}.issubset(merged.columns):
            merged["window_id"] = merged.apply(
                lambda row: f"{row['chromosome']}:{int(row['start'])}-{int(row['end'])}",
                axis=1,
            )
        if {"modified_count", "valid_count"}.issubset(merged.columns):
            merged["window_fraction"] = _safe_fraction(
                merged["modified_count"], merged["valid_count"]
            )
        return merged

    merged_rows: list[dict[str, object]] = []
    current_group: list[pd.Series] = [ordered.iloc[0]]

    def _can_merge(previous: pd.Series, row: pd.Series) -> bool:
        same_chromosome = row.get("chromosome") == previous.get("chromosome")
        same_strand = row.get("strand") == previous.get("strand")
        within_distance = int(row["start"]) <= int(previous["end"]) + merge_distance
        return same_chromosome and same_strand and within_distance

    for _, row in ordered.iloc[1:].iterrows():
        if _can_merge(current_group[-1], row):
            current_group.append(row)
            continue
        merged_rows.append(_build_merged_hit(pd.DataFrame(current_group)))
        current_group = [row]

    merged_rows.append(_build_merged_hit(pd.DataFrame(current_group)))
    merged = pd.DataFrame(merged_rows)
    merged = _sort_hits_for_output(merged)

    if "window_id" not in merged.columns and {"chromosome", "start", "end"}.issubset(merged.columns):
        merged["window_id"] = merged.apply(
            lambda row: f"{row['chromosome']}:{int(row['start'])}-{int(row['end'])}",
            axis=1,
        )

    if "merged_window_count" not in merged.columns:
        merged["merged_window_count"] = 1

    if {"modified_count", "valid_count"}.issubset(merged.columns) and "window_fraction" not in merged.columns:
        merged["window_fraction"] = _safe_fraction(merged["modified_count"], merged["valid_count"])

    ordered_columns = list(hits.columns)
    if "rank" in merged.columns and "rank" not in ordered_columns:
        ordered_columns.append("rank")
    for column in ["window_id", "window_fraction", "merged_window_count"]:
        if column in merged.columns and column not in ordered_columns:
            ordered_columns.append(column)
    merged = merged.loc[:, [column for column in ordered_columns if column in merged.columns]]
    return merged.reset_index(drop=True)


def hits_to_bed(hits: pd.DataFrame) -> pd.DataFrame:
    if hits.empty:
        return pd.DataFrame(columns=_BED_COLUMNS)

    ordered = _sort_hits_for_output(hits)

    chrom_column = "chromosome" if "chromosome" in ordered.columns else "chrom"
    chrom = ordered[chrom_column]
    if "window_id" in ordered.columns:
        name = ordered["window_id"]
    else:
        name = ordered.apply(
            lambda row: f"{row.get('chromosome', row.get('chrom'))}:{int(row['start'])}-{int(row['end'])}",
            axis=1,
        )
    if "score_value" in ordered.columns:
        score = (
            ordered["score_value"]
            .fillna(0.0)
            .mul(1000)
            .round()
            .clip(lower=0, upper=1000)
            .astype(int)
        )
    else:
        score = pd.Series([0] * len(ordered), index=ordered.index, dtype="int64")
    strand = ordered["strand"] if "strand" in ordered.columns else "."
    strand = strand.where(strand.isin({"+", "-"}), ".") if isinstance(strand, pd.Series) else strand

    bed = pd.DataFrame(
        {
            "chrom": chrom,
            "start": ordered["start"],
            "end": ordered["end"],
            "name": name,
            "score": score,
            "strand": strand,
        }
    )
    return bed.loc[:, _BED_COLUMNS].reset_index(drop=True)


def _merge_adjacent_hits(hits: pd.DataFrame, *, merge_distance: int) -> pd.DataFrame:
    return merge_adjacent_hits(hits, merge_distance=merge_distance)


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
    pairing_policy: str | None = None,
    rank_by: str | None = None,
    quiet: bool = True,
    cores: int | None = None,
) -> RegionDiscoveryResult:
    motif_list = _validate_motifs(motifs)
    if score not in {"effect_size_only", "beta_binomial"}:
        raise ValueError("scan_genome requires score in {'effect_size_only', 'beta_binomial'}.")
    if score == "beta_binomial" and contrast is None:
        raise ValueError("scan_genome score='beta_binomial' requires an explicit contrast.")
    if contrast is not None and contrast.mode not in {
        "pairwise",
        "group_vs_group",
        "matched_pairwise",
        "time_course",
    }:
        raise ValueError(
            "scan_genome currently supports only pairwise/group_vs_group and paired contrast modes."
        )

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

    pairing_metadata: dict[str, object] = {}
    active_pairing_policy = _pairing_policy_value(pairing_policy)
    active_rank_by = rank_by
    if _is_paired_contrast(contrast):
        if contrast.mode == "matched_pairwise":
            numerator_conditions = list(contrast.numerator or [])
            denominator_conditions = list(contrast.denominator or [])
            if len(numerator_conditions) != 1 or len(denominator_conditions) != 1:
                raise ValueError(
                    "scan_genome matched_pairwise requires exactly one numerator and one denominator condition."
                )
        required_conditions = list(contrast.time_order or []) if contrast.mode == "time_course" else list(
            dict.fromkeys((contrast.numerator or []) + (contrast.denominator or []))
        )
        window_summary, pairing_metadata = _build_paired_window_table(
            window_summary,
            samples=samples,
            pairing_key=contrast.pairing_key,
            required_conditions=required_conditions,
            pairing_policy=active_pairing_policy,
        )
        if contrast.mode == "matched_pairwise":
            if score != "effect_size_only":
                raise ValueError("scan_genome matched_pairwise currently supports score='effect_size_only'.")
            if merge_hits:
                raise ValueError("scan_genome matched_pairwise does not support merge_hits=True.")
            active_rank_by = active_rank_by or "mean_abs_delta"
            window_totals = _aggregate_window_counts(window_summary)
            scored = _score_matched_pairwise(
                window_summary,
                contrast=contrast,
                rank_by=active_rank_by,
            )
            ranked = _rank_windows(scored, score=score)
            scored_columns = [
                column
                for column in ranked.columns
                if column not in _WINDOW_KEY_COLUMNS
                and column not in {"modified_count", "valid_count", "window_fraction"}
            ]
            window_table = window_totals.merge(
                ranked.loc[:, _WINDOW_KEY_COLUMNS + scored_columns],
                on=_WINDOW_KEY_COLUMNS,
                how="left",
                sort=False,
            )
            covered_mask = window_table["valid_count"] >= min_coverage
            paired_score_columns = [
                column
                for column in [
                    "mean_delta",
                    "mean_abs_delta",
                    "delta_sd",
                    "sign_agreement",
                    "n_pairs_used",
                    "score_value",
                    "p_value",
                    "adjusted_p_value",
                    "rank",
                ]
                if column in window_table.columns
            ]
            window_table.loc[~covered_mask, paired_score_columns] = pd.NA
            hits = _sort_hits_for_output(window_table.loc[covered_mask].copy())
            if not hits.empty:
                hits["rank"] = range(1, len(hits) + 1)

            plot_data = {
                "window_score_table": window_table.copy(),
                "top_hits_table": hits.copy(),
            }
            metadata = {
                "analysis_unit": "ensemble_region",
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
                "contrast_mode": contrast.mode,
                "contrast_numerator": list(contrast.numerator or []),
                "contrast_denominator": list(contrast.denominator or []),
                "pairing_policy": active_pairing_policy,
                "n_pairs_used": pairing_metadata["n_pairs_used"],
                "n_pairs_dropped": pairing_metadata["n_pairs_dropped"],
                "paired_mode": contrast.mode,
                "rank_by": active_rank_by,
            }
            return RegionDiscoveryResult(
                hits=hits.reset_index(drop=True),
                windows=window_table.reset_index(drop=True),
                contrast=contrast,
                plot_data=plot_data,
                metadata=metadata,
                figures={},
            )

    window_totals = _aggregate_window_counts(window_summary)
    condition_counts = _aggregate_condition_counts(window_summary)
    _validate_contrast_conditions(condition_counts=condition_counts, contrast=contrast)

    covered_mask = window_totals["valid_count"] >= min_coverage
    covered_window_totals = window_totals.loc[covered_mask].copy()
    covered_keys = covered_window_totals.loc[:, _WINDOW_KEY_COLUMNS]
    covered_condition_counts = condition_counts.merge(
        covered_keys,
        on=_WINDOW_KEY_COLUMNS,
        how="inner",
    )

    scored = _score_with_contrast(
        window_totals=covered_window_totals,
        condition_counts=covered_condition_counts,
        contrast=contrast,
        score=score,
    )
    ranked = _rank_windows(scored, score=score)
    scored_columns = [
        column
        for column in ranked.columns
        if column not in _WINDOW_KEY_COLUMNS
        and column not in {"modified_count", "valid_count", "window_fraction"}
    ]
    window_table = window_totals.merge(
        ranked.loc[:, _WINDOW_KEY_COLUMNS + scored_columns],
        on=_WINDOW_KEY_COLUMNS,
        how="left",
        sort=False,
    )
    hits = ranked.copy()

    if merge_hits:
        hits = merge_adjacent_hits(hits, merge_distance=merge_distance)

    plot_data = {
        "window_score_table": window_table.copy(),
        "top_hits_table": hits.copy(),
    }
    metadata = {
        "analysis_unit": "ensemble_region",
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
    if pairing_metadata:
        metadata["pairing_policy"] = active_pairing_policy
        metadata["n_pairs_used"] = pairing_metadata["n_pairs_used"]
        metadata["n_pairs_dropped"] = pairing_metadata["n_pairs_dropped"]
        metadata["paired_mode"] = contrast.mode if contrast is not None else None
        metadata["rank_by"] = active_rank_by

    return RegionDiscoveryResult(
        hits=hits.reset_index(drop=True),
        windows=window_table.reset_index(drop=True),
        contrast=contrast,
        plot_data=plot_data,
        metadata=metadata,
        figures={},
    )
