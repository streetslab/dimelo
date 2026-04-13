from __future__ import annotations

from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

from dimelo.models import ContrastSpec, SharedClusterContrastResult, SharedClusterResult

_SUPPORTED_SHARED_CLUSTER_MODES = {
    "pairwise",
    "matched_pairwise",
    "group_vs_group",
    "time_course",
}

_SUPPORTED_POOLED_SCREEN_TESTS = {"chi_squared", "g_test"}


def _require_supported_shared_cluster_mode(contrast: ContrastSpec) -> None:
    if contrast.mode not in _SUPPORTED_SHARED_CLUSTER_MODES:
        raise NotImplementedError(
            f"Shared cluster tests are not implemented for contrast mode '{contrast.mode}'."
        )


def _require_requested_conditions_present(
    sample_table: pd.DataFrame,
    *,
    side_specs: dict[str, list[str]],
) -> None:
    available_conditions = set(sample_table["condition"].dropna().unique())
    for side, conditions in side_specs.items():
        missing = sorted(set(conditions) - available_conditions)
        if missing:
            raise ValueError(
                f"Missing {side} evidence for requested condition(s): {', '.join(missing)}."
            )


def _require_time_order_present(sample_table: pd.DataFrame, time_order: list[str]) -> None:
    available_conditions = set(sample_table["condition"].dropna().unique())
    missing = [condition for condition in time_order if condition not in available_conditions]
    if missing:
        raise ValueError(
            "Shared cluster time_course requested missing time_order condition(s): "
            + ", ".join(missing)
        )


def _sample_fraction_table(result: SharedClusterResult) -> pd.DataFrame:
    required = {"sample_id", "condition", "cluster", "fraction"}
    missing = required - set(result.cluster_distribution.columns)
    if missing:
        raise ValueError(
            "Shared cluster tests require cluster_distribution columns: "
            + ", ".join(sorted(missing))
            + "."
        )
    return result.cluster_distribution.loc[
        :, ["sample_id", "condition", "cluster", "fraction"]
    ].copy()


def _condition_count_table(result: SharedClusterResult) -> pd.DataFrame:
    required = {"condition", "cluster", "count", "fraction"}
    missing = required - set(result.condition_distribution.columns)
    if missing:
        raise ValueError(
            "Shared cluster tests require condition_distribution columns: "
            + ", ".join(sorted(missing))
            + "."
        )
    return result.condition_distribution.loc[
        :, ["condition", "cluster", "count", "fraction"]
    ].copy()


def _cluster_order(result: SharedClusterResult) -> list[str]:
    cluster_labels = list(result.model.cluster_labels)
    if cluster_labels:
        return cluster_labels
    return (
        result.cluster_distribution["cluster"].drop_duplicates().astype(str).sort_values().tolist()
    )


def _composition_effect_size(summary_table: pd.DataFrame) -> float:
    return float(summary_table["delta_fraction"].abs().sum() / 2.0)


def _composition_effect_size_from_vectors(observed: np.ndarray, reference: np.ndarray) -> float:
    return float(np.abs(observed - reference).sum() / 2.0)


def _trend_statistic(timepoint_fractions: np.ndarray) -> float:
    if timepoint_fractions.shape[0] < 2:
        return 0.0

    positions = np.arange(timepoint_fractions.shape[0], dtype=float)
    slopes = [
        np.polyfit(positions, timepoint_fractions[:, cluster_index], 1)[0]
        for cluster_index in range(timepoint_fractions.shape[1])
    ]
    return float(np.abs(slopes).sum() / 2.0)


def _adjust_p_values_bh(p_values: pd.Series) -> pd.Series:
    if p_values.empty:
        return pd.Series(dtype=float, index=p_values.index)

    numeric = pd.to_numeric(p_values, errors="coerce")
    valid_mask = numeric.notna()
    if not valid_mask.any():
        return pd.Series(np.nan, index=p_values.index, dtype=float)

    valid_values = numeric[valid_mask].to_numpy(dtype=float)
    order = np.argsort(valid_values, kind="mergesort")
    ranked = valid_values[order]
    adjusted = np.empty_like(ranked)
    running_min = 1.0
    total = len(ranked)

    for index in range(total - 1, -1, -1):
        rank = index + 1
        candidate = min(1.0, ranked[index] * total / rank)
        running_min = min(running_min, candidate)
        adjusted[index] = running_min

    restored = np.full(len(p_values), np.nan, dtype=float)
    restored_indices = np.flatnonzero(valid_mask.to_numpy())
    restored[restored_indices[order]] = adjusted
    return pd.Series(restored, index=p_values.index, dtype=float)


def _contrast_id(contrast: ContrastSpec) -> str:
    numerator = "+".join(contrast.numerator or [])
    denominator = "+".join(contrast.denominator or [])
    if contrast.mode == "time_course":
        time_order = "+".join(contrast.time_order or [])
        return f"time_course:{time_order}"
    return f"{numerator}_vs_{denominator}"


def _sample_cluster_matrix(result: SharedClusterResult) -> pd.DataFrame:
    sample_table = _sample_fraction_table(result)
    cluster_order = _cluster_order(result)
    sample_info = sample_table.loc[:, ["sample_id", "condition"]].drop_duplicates()
    matrix = (
        sample_table.pivot_table(
            index="sample_id",
            columns="cluster",
            values="fraction",
            aggfunc="first",
            fill_value=0.0,
        )
        .reindex(columns=cluster_order, fill_value=0.0)
        .reset_index()
    )
    return sample_info.merge(matrix, on="sample_id", how="left")


def _add_pairing_metadata(
    sample_matrix: pd.DataFrame,
    *,
    assignments: pd.DataFrame,
    pairing_key: str,
) -> pd.DataFrame:
    if pairing_key not in assignments.columns:
        raise ValueError(
            f"Shared cluster matched_pairwise requires assignments column {pairing_key!r}."
        )

    pairing = assignments.loc[:, ["sample_id", pairing_key]].drop_duplicates()
    pair_counts = pairing.groupby("sample_id", dropna=False)[pairing_key].nunique(dropna=False)
    if (pair_counts > 1).any():
        raise ValueError(
            "Shared cluster matched_pairwise requires each sample_id to map to exactly one pairing key."
        )
    if pairing[pairing_key].isna().any():
        raise ValueError(
            f"Shared cluster matched_pairwise requires non-null values in column {pairing_key!r}."
        )
    return sample_matrix.merge(pairing, on="sample_id", how="left")


def _prepare_unpaired_group_table(
    result: SharedClusterResult,
    *,
    contrast: ContrastSpec,
) -> tuple[pd.DataFrame, list[str]]:
    sample_matrix = _sample_cluster_matrix(result)
    cluster_order = _cluster_order(result)
    side_specs = {
        "numerator": contrast.numerator or [],
        "denominator": contrast.denominator or [],
    }
    _require_requested_conditions_present(sample_matrix, side_specs=side_specs)
    requested_conditions = set(side_specs["numerator"]) | set(side_specs["denominator"])
    filtered = sample_matrix.loc[sample_matrix["condition"].isin(requested_conditions)].copy()
    filtered["contrast_side"] = np.where(
        filtered["condition"].isin(side_specs["numerator"]),
        "numerator",
        "denominator",
    )
    return filtered, cluster_order


def _prepare_paired_group_table(
    result: SharedClusterResult,
    *,
    contrast: ContrastSpec,
) -> tuple[pd.DataFrame, list[str]]:
    if len(contrast.numerator or []) != 1 or len(contrast.denominator or []) != 1:
        raise ValueError(
            "Shared cluster matched_pairwise currently requires exactly one numerator "
            "and one denominator condition."
        )

    sample_matrix = _sample_cluster_matrix(result)
    sample_matrix = _add_pairing_metadata(
        sample_matrix,
        assignments=result.assignments,
        pairing_key=contrast.pairing_key or "",
    )
    cluster_order = _cluster_order(result)
    side_specs = {
        "numerator": contrast.numerator or [],
        "denominator": contrast.denominator or [],
    }
    _require_requested_conditions_present(sample_matrix, side_specs=side_specs)
    requested_conditions = set(side_specs["numerator"]) | set(side_specs["denominator"])
    filtered = sample_matrix.loc[sample_matrix["condition"].isin(requested_conditions)].copy()
    filtered["contrast_side"] = np.where(
        filtered["condition"].isin(side_specs["numerator"]),
        "numerator",
        "denominator",
    )

    grouped = (
        filtered.groupby([contrast.pairing_key, "contrast_side"], as_index=False, sort=False)[
            cluster_order
        ]
        .mean()
    )
    side_sets = (
        grouped.loc[:, [contrast.pairing_key, "contrast_side"]]
        .drop_duplicates()
        .groupby(contrast.pairing_key)["contrast_side"]
        .agg(lambda values: set(values))
    )
    complete_pairs = [
        pair_id for pair_id, sides in side_sets.items() if {"numerator", "denominator"} <= sides
    ]
    if not complete_pairs:
        raise ValueError("Shared cluster matched_pairwise found no complete matched units.")

    return grouped.loc[grouped[contrast.pairing_key].isin(complete_pairs)].copy(), cluster_order


def _permutation_p_value(observed: float, permuted: np.ndarray) -> float:
    return float((1 + np.count_nonzero(permuted >= observed - 1e-12)) / (len(permuted) + 1))


def _run_unpaired_permutations(
    values: np.ndarray,
    *,
    n_numerator: int,
    n_permutations: int,
    random_state: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(random_state)
    cluster_stats = np.zeros((n_permutations, values.shape[1]), dtype=float)
    omnibus_stats = np.zeros(n_permutations, dtype=float)

    for permutation_index in range(n_permutations):
        permuted_indices = rng.permutation(values.shape[0])
        numerator = values[permuted_indices[:n_numerator]]
        denominator = values[permuted_indices[n_numerator:]]
        delta = numerator.mean(axis=0) - denominator.mean(axis=0)
        cluster_stats[permutation_index] = np.abs(delta)
        omnibus_stats[permutation_index] = np.abs(delta).sum() / 2.0

    return cluster_stats, omnibus_stats


def _run_paired_permutations(
    pair_deltas: np.ndarray,
    *,
    n_permutations: int,
    random_state: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(random_state)
    cluster_stats = np.zeros((n_permutations, pair_deltas.shape[1]), dtype=float)
    omnibus_stats = np.zeros(n_permutations, dtype=float)

    for permutation_index in range(n_permutations):
        signs = rng.choice(np.array([-1.0, 1.0]), size=pair_deltas.shape[0])[:, None]
        delta = (pair_deltas * signs).mean(axis=0)
        cluster_stats[permutation_index] = np.abs(delta)
        omnibus_stats[permutation_index] = np.abs(delta).sum() / 2.0

    return cluster_stats, omnibus_stats


def _score_unpaired_group(
    group_table: pd.DataFrame,
    *,
    cluster_order: list[str],
    n_permutations: int,
    random_state: int | None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    numerator = group_table.loc[group_table["contrast_side"] == "numerator", cluster_order]
    denominator = group_table.loc[group_table["contrast_side"] == "denominator", cluster_order]
    if numerator.empty or denominator.empty:
        raise ValueError("Shared cluster tests require evidence for both contrast sides.")

    observed_fraction = numerator.mean(axis=0)
    reference_fraction = denominator.mean(axis=0)
    delta = observed_fraction - reference_fraction
    details = pd.DataFrame(
        {
            "cluster": cluster_order,
            "fraction": observed_fraction.to_numpy(dtype=float),
            "reference_fraction": reference_fraction.to_numpy(dtype=float),
            "delta_fraction": delta.to_numpy(dtype=float),
        }
    )

    perm_cluster, perm_omnibus = _run_unpaired_permutations(
        group_table.loc[:, cluster_order].to_numpy(dtype=float),
        n_numerator=len(numerator),
        n_permutations=n_permutations,
        random_state=random_state,
    )
    return details, {
        "permuted_cluster_stats": perm_cluster,
        "permuted_omnibus_stats": perm_omnibus,
        "paired": False,
        "n_numerator": int(len(numerator)),
        "n_denominator": int(len(denominator)),
    }


def _score_paired_group(
    pair_table: pd.DataFrame,
    *,
    cluster_order: list[str],
    contrast: ContrastSpec,
    n_permutations: int,
    random_state: int | None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    pair_key = contrast.pairing_key or ""
    numerator = (
        pair_table.loc[pair_table["contrast_side"] == "numerator"]
        .set_index(pair_key)
        .reindex(columns=cluster_order)
        .sort_index()
    )
    denominator = (
        pair_table.loc[pair_table["contrast_side"] == "denominator"]
        .set_index(pair_key)
        .reindex(columns=cluster_order)
        .sort_index()
    )
    common_pairs = numerator.index.intersection(denominator.index)
    if common_pairs.empty:
        raise ValueError("Shared cluster matched_pairwise found no complete matched units.")

    numerator = numerator.loc[common_pairs]
    denominator = denominator.loc[common_pairs]
    observed_fraction = numerator.mean(axis=0)
    reference_fraction = denominator.mean(axis=0)
    delta = observed_fraction - reference_fraction
    details = pd.DataFrame(
        {
            "cluster": cluster_order,
            "fraction": observed_fraction.to_numpy(dtype=float),
            "reference_fraction": reference_fraction.to_numpy(dtype=float),
            "delta_fraction": delta.to_numpy(dtype=float),
        }
    )

    pair_deltas = numerator.to_numpy(dtype=float) - denominator.to_numpy(dtype=float)
    perm_cluster, perm_omnibus = _run_paired_permutations(
        pair_deltas,
        n_permutations=n_permutations,
        random_state=random_state,
    )
    return details, {
        "permuted_cluster_stats": perm_cluster,
        "permuted_omnibus_stats": perm_omnibus,
        "paired": True,
        "n_pairs_used": int(len(common_pairs)),
    }


def _finalize_details(
    details: pd.DataFrame,
    *,
    permuted_cluster_stats: np.ndarray,
) -> pd.DataFrame:
    finalized = details.copy()
    finalized["log2_fc"] = np.log2(
        (finalized["fraction"] + 1e-9) / (finalized["reference_fraction"] + 1e-9)
    )
    finalized["effect_size"] = finalized["delta_fraction"].abs()
    finalized["p_value"] = [
        _permutation_p_value(
            float(observed),
            permuted_cluster_stats[:, cluster_index],
        )
        for cluster_index, observed in enumerate(finalized["effect_size"])
    ]
    finalized["adjusted_p_value"] = _adjust_p_values_bh(finalized["p_value"])
    finalized = finalized.sort_values(
        ["adjusted_p_value", "p_value", "effect_size", "cluster"],
        ascending=[True, True, False, True],
        kind="stable",
    ).reset_index(drop=True)
    finalized["rank"] = np.arange(1, len(finalized) + 1, dtype=int)
    return finalized


def _build_summary(
    details: pd.DataFrame,
    *,
    contrast: ContrastSpec,
    omnibus_p_value: float,
    stats_metadata: dict[str, Any],
    test: str,
    trend_p_value: float | None = None,
) -> pd.DataFrame:
    top_cluster_row = details.sort_values(
        ["effect_size", "p_value", "cluster"],
        ascending=[False, True, True],
        kind="stable",
    ).iloc[0]
    summary_row: dict[str, Any] = {
        "contrast_id": _contrast_id(contrast),
        "contrast_mode": contrast.mode,
        "test": test,
        "numerator_conditions": "+".join(contrast.numerator or []),
        "denominator_conditions": "+".join(contrast.denominator or []),
        "composition_effect_size": _composition_effect_size(details),
        "omnibus_p_value": float(omnibus_p_value),
        "top_cluster": top_cluster_row["cluster"],
        "top_cluster_delta_fraction": float(top_cluster_row["delta_fraction"]),
        "effect_size_metric": "total_variation_distance",
        "paired": bool(stats_metadata["paired"]),
    }
    if trend_p_value is not None:
        summary_row["trend_p_value"] = float(trend_p_value)

    for key, value in stats_metadata.items():
        if key.startswith("n_"):
            summary_row[key] = value
    return pd.DataFrame([summary_row])


def _build_pooled_omnibus_p_value(
    *,
    result: SharedClusterResult,
    contrast: ContrastSpec,
    test: str,
) -> float:
    if test not in _SUPPORTED_POOLED_SCREEN_TESTS:
        raise ValueError(
            "shared_cluster_tests currently requires test='permutation', "
            "'chi_squared', or 'g_test'."
        )

    condition_table = _condition_count_table(result)
    side_specs = {
        "numerator": contrast.numerator or [],
        "denominator": contrast.denominator or [],
    }
    _require_requested_conditions_present(condition_table, side_specs=side_specs)
    requested_conditions = side_specs["numerator"] + side_specs["denominator"]
    filtered = condition_table.loc[condition_table["condition"].isin(requested_conditions)].copy()
    contingency = (
        filtered.pivot_table(
            index="condition",
            columns="cluster",
            values="count",
            aggfunc="sum",
            fill_value=0,
        )
        .reindex(index=requested_conditions, columns=_cluster_order(result), fill_value=0)
        .astype(float)
    )
    if contingency.shape[0] < 2 or contingency.shape[1] < 2:
        raise ValueError("shared_cluster_tests pooled screening requires at least 2x2 counts.")

    if test == "chi_squared":
        _, omnibus_p_value, _, _ = stats.chi2_contingency(
            contingency.to_numpy(dtype=float),
            correction=False,
        )
    else:
        _, omnibus_p_value, _, _ = stats.chi2_contingency(
            contingency.to_numpy(dtype=float),
            correction=False,
            lambda_="log-likelihood",
        )
    return float(omnibus_p_value)


def _score_time_course(
    *,
    result: SharedClusterResult,
    contrast: ContrastSpec,
    n_permutations: int,
    random_state: int | None,
    include_pairwise: bool,
) -> SharedClusterContrastResult:
    time_order = list(dict.fromkeys(contrast.time_order or []))
    if not time_order:
        raise ValueError("Shared cluster time_course requires contrast.time_order.")

    sample_matrix = _sample_cluster_matrix(result)
    _require_time_order_present(sample_matrix, time_order)
    cluster_order = _cluster_order(result)

    filtered = sample_matrix.loc[sample_matrix["condition"].isin(time_order)].copy()
    filtered["condition"] = pd.Categorical(filtered["condition"], categories=time_order, ordered=True)
    filtered = filtered.sort_values(["condition", "sample_id"], kind="stable").reset_index(drop=True)

    grouped = (
        filtered.groupby("condition", sort=False, observed=False)[cluster_order]
        .mean()
        .reindex(time_order)
    )
    counts = (
        filtered.groupby("condition", sort=False, observed=False)["sample_id"]
        .nunique()
        .reindex(time_order)
    )
    time_course_table = grouped.reset_index(names="timepoint")
    time_course_table["n_samples"] = counts.to_numpy(dtype=int)

    observed_matrix = grouped.to_numpy(dtype=float)
    observed_first = observed_matrix[0]
    observed_last = observed_matrix[-1]
    observed_delta = observed_last - observed_first
    observed_omnibus = _composition_effect_size_from_vectors(observed_last, observed_first)
    observed_trend = _trend_statistic(observed_matrix)

    rng = np.random.default_rng(random_state)
    permuted_cluster_stats = np.zeros((n_permutations, len(cluster_order)), dtype=float)
    permuted_omnibus_stats = np.zeros(n_permutations, dtype=float)
    permuted_trend_stats = np.zeros(n_permutations, dtype=float)
    condition_values = filtered["condition"].to_numpy()
    value_matrix = filtered.loc[:, cluster_order].to_numpy(dtype=float)

    for permutation_index in range(n_permutations):
        permuted_conditions = rng.permutation(condition_values)
        permuted = pd.DataFrame(value_matrix, columns=cluster_order)
        permuted["condition"] = permuted_conditions
        permuted_grouped = (
            permuted.groupby("condition", sort=False, observed=False)[cluster_order]
            .mean()
            .reindex(time_order)
        )
        permuted_matrix = permuted_grouped.to_numpy(dtype=float)
        permuted_delta = permuted_matrix[-1] - permuted_matrix[0]
        permuted_cluster_stats[permutation_index] = np.abs(permuted_delta)
        permuted_omnibus_stats[permutation_index] = _composition_effect_size_from_vectors(
            permuted_matrix[-1],
            permuted_matrix[0],
        )
        permuted_trend_stats[permutation_index] = _trend_statistic(permuted_matrix)

    details = pd.DataFrame(
        {
            "cluster": cluster_order,
            "fraction": observed_last,
            "reference_fraction": observed_first,
            "delta_fraction": observed_delta,
        }
    )
    details = _finalize_details(details, permuted_cluster_stats=permuted_cluster_stats)

    omnibus_p_value = _permutation_p_value(observed_omnibus, permuted_omnibus_stats)
    trend_p_value = _permutation_p_value(observed_trend, permuted_trend_stats)
    summary = _build_summary(
        details,
        contrast=contrast,
        omnibus_p_value=omnibus_p_value,
        stats_metadata={
            "paired": False,
            "n_timepoints": int(len(time_order)),
            "n_samples": int(filtered["sample_id"].nunique()),
        },
        test="permutation",
        trend_p_value=trend_p_value,
    )
    summary["test"] = "permutation"
    summary["trend_p_value"] = float(trend_p_value)

    pairwise = None
    plot_data: dict[str, pd.DataFrame | dict[str, Any]] = {
        "summary_table": summary.copy(),
        "cluster_effect_table": details.copy(),
        "time_course_table": time_course_table.copy(),
    }

    if include_pairwise:
        pairwise_rows = []
        for left_timepoint, right_timepoint in zip(time_order[:-1], time_order[1:]):
            pair_table = filtered.loc[
                filtered["condition"].isin([left_timepoint, right_timepoint])
            ].copy()
            pair_table["contrast_side"] = np.where(
                pair_table["condition"] == left_timepoint,
                "numerator",
                "denominator",
            )
            pair_details, stats_metadata = _score_unpaired_group(
                pair_table,
                cluster_order=cluster_order,
                n_permutations=n_permutations,
                random_state=random_state,
            )
            pair_details = _finalize_details(
                pair_details,
                permuted_cluster_stats=stats_metadata["permuted_cluster_stats"],
            )
            pair_omnibus = _permutation_p_value(
                _composition_effect_size(pair_details),
                stats_metadata["permuted_omnibus_stats"],
            )
            pair_summary = _build_summary(
                pair_details,
                contrast=ContrastSpec(
                    mode="pairwise",
                    numerator=[left_timepoint],
                    denominator=[right_timepoint],
                ),
                omnibus_p_value=pair_omnibus,
                stats_metadata=stats_metadata,
                test="permutation",
            )
            row = pair_summary.iloc[0].to_dict()
            row["left_timepoint"] = left_timepoint
            row["right_timepoint"] = right_timepoint
            pairwise_rows.append(row)

        pairwise = pd.DataFrame(pairwise_rows)
        plot_data["pairwise_table"] = pairwise.copy()

    metadata: dict[str, Any] = {
        "contrast_mode": contrast.mode,
        "contrast_id": _contrast_id(contrast),
        "paired": False,
        "pairing_key": None,
        "test": "permutation",
        "multiple_testing": "fdr_bh",
        "n_permutations": int(n_permutations),
        "random_state": random_state,
        "inference_level": "replicate_aware",
        "n_timepoints": int(len(time_order)),
        "n_samples": int(filtered["sample_id"].nunique()),
    }
    return SharedClusterContrastResult(
        summary=summary,
        details=details,
        pairwise=pairwise,
        plot_data=plot_data,
        metadata=metadata,
    )


def shared_cluster_tests(
    *,
    result: SharedClusterResult,
    contrast: ContrastSpec,
    test: str = "permutation",
    multiple_testing: str = "fdr_bh",
    n_permutations: int = 1000,
    random_state: int | None = 42,
    include_pairwise: bool = False,
) -> SharedClusterContrastResult:
    _require_supported_shared_cluster_mode(contrast)
    if multiple_testing != "fdr_bh":
        raise ValueError("shared_cluster_tests currently requires multiple_testing='fdr_bh'.")
    if n_permutations <= 0:
        raise ValueError("shared_cluster_tests requires n_permutations > 0.")

    if contrast.mode == "time_course":
        if test != "permutation":
            raise ValueError("shared_cluster_tests time_course currently requires test='permutation'.")
        return _score_time_course(
            result=result,
            contrast=contrast,
            n_permutations=n_permutations,
            random_state=random_state,
            include_pairwise=include_pairwise,
        )

    if test == "permutation":
        if contrast.mode == "matched_pairwise":
            pair_table, cluster_order = _prepare_paired_group_table(result, contrast=contrast)
            details, stats_metadata = _score_paired_group(
                pair_table,
                cluster_order=cluster_order,
                contrast=contrast,
                n_permutations=n_permutations,
                random_state=random_state,
            )
            sample_plot_table = pair_table.copy()
        else:
            group_table, cluster_order = _prepare_unpaired_group_table(result, contrast=contrast)
            details, stats_metadata = _score_unpaired_group(
                group_table,
                cluster_order=cluster_order,
                n_permutations=n_permutations,
                random_state=random_state,
            )
            sample_plot_table = group_table.copy()

        details = _finalize_details(
            details,
            permuted_cluster_stats=stats_metadata["permuted_cluster_stats"],
        )
        omnibus_p_value = _permutation_p_value(
            _composition_effect_size(details),
            stats_metadata["permuted_omnibus_stats"],
        )
        summary = _build_summary(
            details,
            contrast=contrast,
            omnibus_p_value=omnibus_p_value,
            stats_metadata=stats_metadata,
            test=test,
        )
        metadata: dict[str, Any] = {
            "contrast_mode": contrast.mode,
            "contrast_id": _contrast_id(contrast),
            "paired": bool(stats_metadata["paired"]),
            "pairing_key": contrast.pairing_key if contrast.mode == "matched_pairwise" else None,
            "test": test,
            "multiple_testing": multiple_testing,
            "n_permutations": int(n_permutations),
            "random_state": random_state,
            "inference_level": "replicate_aware",
        }
        metadata.update({key: value for key, value in stats_metadata.items() if key.startswith("n_")})
        plot_data = {
            "summary_table": summary.copy(),
            "cluster_effect_table": details.copy(),
            "sample_fraction_table": sample_plot_table,
        }
        return SharedClusterContrastResult(
            summary=summary,
            details=details,
            plot_data=plot_data,
            metadata=metadata,
        )

    if test not in _SUPPORTED_POOLED_SCREEN_TESTS:
        raise ValueError(
            "shared_cluster_tests currently requires test='permutation', "
            "'chi_squared', or 'g_test'."
        )

    group_table, cluster_order = _prepare_unpaired_group_table(result, contrast=contrast)
    details, stats_metadata = _score_unpaired_group(
        group_table,
        cluster_order=cluster_order,
        n_permutations=n_permutations,
        random_state=random_state,
    )
    details = _finalize_details(
        details,
        permuted_cluster_stats=stats_metadata["permuted_cluster_stats"],
    )
    omnibus_p_value = _build_pooled_omnibus_p_value(result=result, contrast=contrast, test=test)
    summary = _build_summary(
        details,
        contrast=contrast,
        omnibus_p_value=omnibus_p_value,
        stats_metadata=stats_metadata,
        test=test,
    )
    metadata: dict[str, Any] = {
        "contrast_mode": contrast.mode,
        "contrast_id": _contrast_id(contrast),
        "paired": bool(stats_metadata["paired"]),
        "pairing_key": None,
        "test": test,
        "multiple_testing": multiple_testing,
        "n_permutations": int(n_permutations),
        "random_state": random_state,
        "inference_level": "pooled_screen",
    }
    metadata.update({key: value for key, value in stats_metadata.items() if key.startswith("n_")})
    plot_data = {
        "summary_table": summary.copy(),
        "cluster_effect_table": details.copy(),
        "sample_fraction_table": group_table.copy(),
    }
    return SharedClusterContrastResult(
        summary=summary,
        details=details,
        plot_data=plot_data,
        metadata=metadata,
    )
