from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

_CLUSTER_DISTRIBUTION_COLUMNS = [
    "sample_id",
    "condition",
    "cluster",
    "count",
    "fraction",
]

_CONDITION_DISTRIBUTION_COLUMNS = [
    "condition",
    "cluster",
    "count",
    "fraction",
    "replicate_n",
]


def _require_columns(frame: pd.DataFrame, required: Iterable[str], name: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {', '.join(missing)}")


def build_cluster_distribution(assignments: pd.DataFrame) -> pd.DataFrame:
    _require_columns(assignments, ("sample_id", "condition", "cluster"), "assignments")
    if assignments.empty:
        return pd.DataFrame(columns=_CLUSTER_DISTRIBUTION_COLUMNS)

    grouped = (
        assignments.groupby(["sample_id", "condition", "cluster"], sort=True)
        .size()
        .reset_index(name="count")
    )
    totals = grouped.groupby("sample_id")["count"].transform("sum")
    grouped["fraction"] = grouped["count"] / totals
    return grouped[_CLUSTER_DISTRIBUTION_COLUMNS]


def build_condition_distribution(cluster_distribution: pd.DataFrame) -> pd.DataFrame:
    _require_columns(
        cluster_distribution,
        ("sample_id", "condition", "cluster", "count"),
        "cluster_distribution",
    )
    if cluster_distribution.empty:
        return pd.DataFrame(columns=_CONDITION_DISTRIBUTION_COLUMNS)

    replicate_counts = (
        cluster_distribution.groupby("condition", sort=True)["sample_id"]
        .nunique()
        .rename("replicate_n")
        .reset_index()
    )
    grouped = (
        cluster_distribution.groupby(["condition", "cluster"], sort=True)
        .agg(count=("count", "sum"))
        .reset_index()
    )
    grouped = grouped.merge(replicate_counts, on="condition", how="left")
    totals = grouped.groupby("condition")["count"].transform("sum")
    grouped["fraction"] = grouped["count"] / totals
    return grouped[_CONDITION_DISTRIBUTION_COLUMNS]


def build_distribution_change(
    condition_distribution: pd.DataFrame,
    *,
    reference_condition: str,
    pseudo_count: float = 1e-9,
) -> pd.DataFrame:
    _require_columns(
        condition_distribution,
        ("condition", "cluster", "fraction"),
        "condition_distribution",
    )
    if condition_distribution.empty:
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

    reference_rows = condition_distribution[
        condition_distribution["condition"] == reference_condition
    ]
    if reference_rows.empty:
        raise ValueError(
            f"reference_condition {reference_condition!r} not present in condition_distribution"
        )

    reference_fractions = reference_rows.loc[:, ["cluster", "fraction"]].rename(
        columns={"fraction": "reference_fraction"}
    )
    merged = condition_distribution.merge(reference_fractions, on="cluster", how="left")
    merged = merged[merged["condition"] != reference_condition].copy()
    if "count" not in merged.columns:
        merged["count"] = pd.NA
    if "replicate_n" not in merged.columns:
        merged["replicate_n"] = pd.NA
    merged["reference_fraction"] = merged["reference_fraction"].fillna(0.0)
    merged["delta_fraction"] = merged["fraction"] - merged["reference_fraction"]
    merged["log2_fc"] = np.log2(
        (merged["fraction"] + pseudo_count) / (merged["reference_fraction"] + pseudo_count)
    )
    return merged[
        [
            "condition",
            "cluster",
            "count",
            "fraction",
            "replicate_n",
            "reference_fraction",
            "delta_fraction",
            "log2_fc",
        ]
    ].sort_values(["condition", "cluster"], kind="stable").reset_index(drop=True)
