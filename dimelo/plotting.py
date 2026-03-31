from __future__ import annotations

import pandas as pd

from dimelo.distribution import _require_columns


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
