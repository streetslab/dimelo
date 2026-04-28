import pandas as pd

from dimelo.distribution import (
    build_cluster_distribution,
    build_condition_distribution,
    build_distribution_change,
)
from dimelo.plotting import (
    prepare_cluster_distribution_bar_data,
    prepare_cluster_distribution_heatmap_data,
)


def test_build_cluster_distribution_counts_and_fractions():
    assignments = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s1", "s2"],
            "condition": ["NS", "NS", "NS", "15min"],
            "cluster": ["C0", "C0", "C1", "C1"],
        }
    )

    result = build_cluster_distribution(assignments)

    assert list(result.columns) == [
        "sample_id",
        "condition",
        "cluster",
        "count",
        "fraction",
    ]
    s1 = (
        result[result["sample_id"] == "s1"]
        .sort_values("cluster")
        .reset_index(drop=True)
    )
    assert list(s1["count"]) == [2, 1]
    assert list(s1["fraction"]) == [2 / 3, 1 / 3]


def test_build_condition_distribution_aggregates_sample_distributions():
    cluster_distribution = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s2"],
            "condition": ["NS", "NS", "NS", "NS"],
            "cluster": ["C0", "C1", "C0", "C1"],
            "count": [2, 1, 1, 2],
            "fraction": [2 / 3, 1 / 3, 1 / 3, 2 / 3],
        }
    )

    result = build_condition_distribution(cluster_distribution)

    assert list(result.columns) == [
        "condition",
        "cluster",
        "count",
        "fraction",
        "replicate_n",
    ]
    assert set(result["condition"]) == {"NS"}
    assert list(result.sort_values("cluster")["count"]) == [3, 3]
    assert list(result.sort_values("cluster")["fraction"]) == [0.5, 0.5]
    assert set(result["replicate_n"]) == {2}


def test_build_condition_distribution_tracks_condition_level_replicate_count():
    cluster_distribution = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2"],
            "condition": ["NS", "NS", "NS"],
            "cluster": ["C0", "C1", "C0"],
            "count": [2, 1, 1],
            "fraction": [2 / 3, 1 / 3, 1.0],
        }
    )

    result = build_condition_distribution(cluster_distribution)

    assert set(result["replicate_n"]) == {2}


def test_build_distribution_change_against_reference_condition():
    condition_distribution = pd.DataFrame(
        {
            "condition": ["NS", "NS", "15min", "15min"],
            "cluster": ["C0", "C1", "C0", "C1"],
            "count": [3, 1, 1, 3],
            "fraction": [0.75, 0.25, 0.25, 0.75],
            "replicate_n": [2, 2, 2, 2],
        }
    )

    result = build_distribution_change(
        condition_distribution,
        reference_condition="NS",
    )

    row = result[(result["condition"] == "15min") & (result["cluster"] == "C1")].iloc[0]
    assert row["reference_fraction"] == 0.25
    assert row["delta_fraction"] == 0.5
    assert row["log2_fc"] > 0


def test_build_distribution_change_accepts_minimal_condition_distribution_shape():
    condition_distribution = pd.DataFrame(
        {
            "condition": ["NS", "NS", "15min", "15min"],
            "cluster": ["C0", "C1", "C0", "C1"],
            "fraction": [0.75, 0.25, 0.25, 0.75],
        }
    )

    result = build_distribution_change(
        condition_distribution,
        reference_condition="NS",
    )

    assert list(result.columns) == [
        "condition",
        "cluster",
        "count",
        "fraction",
        "replicate_n",
        "reference_fraction",
        "delta_fraction",
        "log2_fc",
    ]
    assert result["count"].isna().all()
    assert result["replicate_n"].isna().all()


def test_prepare_cluster_distribution_bar_data_returns_sorted_frame():
    cluster_distribution = pd.DataFrame(
        {
            "sample_id": ["s2", "s1"],
            "condition": ["15min", "NS"],
            "cluster": ["C1", "C0"],
            "count": [1, 2],
            "fraction": [1 / 3, 2 / 3],
        }
    )

    result = prepare_cluster_distribution_bar_data(cluster_distribution)

    assert list(result.columns) == [
        "sample_id",
        "condition",
        "cluster",
        "count",
        "fraction",
    ]
    assert list(result["sample_id"]) == ["s1", "s2"]
    assert list(result["cluster"]) == ["C0", "C1"]


def test_prepare_cluster_distribution_heatmap_data_pivots_conditions_and_clusters():
    condition_distribution = pd.DataFrame(
        {
            "condition": ["NS", "NS", "15min", "15min"],
            "cluster": ["C0", "C1", "C0", "C1"],
            "count": [3, 1, 1, 3],
            "fraction": [0.75, 0.25, 0.25, 0.75],
            "replicate_n": [2, 2, 2, 2],
        }
    )

    result = prepare_cluster_distribution_heatmap_data(condition_distribution)

    assert list(result.columns) == ["condition", "C0", "C1"]
    ns_row = result[result["condition"] == "NS"].iloc[0]
    assert ns_row["C0"] == 0.75
    assert ns_row["C1"] == 0.25
