import math
from types import SimpleNamespace

import pandas as pd
import pytest

from dimelo import region_contrasts
from dimelo.models import ContrastSpec, SampleSpec


def test_validate_supported_v1_combination():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="ensemble_region",
        representation="modified_fraction",
        signal_source="pileup_counts",
        test="effect_size_only",
    )


def test_validate_accepts_beta_binomial_for_supported_v1_combination():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="ensemble_region",
        representation="modified_fraction",
        signal_source="pileup_counts",
        test="beta_binomial",
    )


def test_validate_region_contrast_request_accepts_cluster_occupancy_fraction_mode():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="effect_size_only",
    )


def test_validate_region_contrast_request_rejects_beta_binomial_for_cluster_occupancy():
    with pytest.raises(ValueError, match="cluster_occupancy"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="cluster_occupancy",
            representation="cluster_fraction",
            signal_source="cluster_occupancy",
            test="beta_binomial",
        )


def test_validate_region_contrast_request_accepts_single_read_mod_fraction():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="single_read",
        representation="read_mod_fraction",
        signal_source="extract_reads",
        test="sample_distribution_shift",
    )


def test_validate_region_contrast_request_accepts_single_read_window_features():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="feature_summary_shift",
    )


def test_validate_region_contrast_request_rejects_single_read_wrong_signal_source():
    with pytest.raises(ValueError, match="extract_reads"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="pileup_counts",
            test="sample_distribution_shift",
        )


def test_validate_region_contrast_request_rejects_single_read_unknown_representation():
    with pytest.raises(ValueError, match="read_mod_fraction"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_shape",
            signal_source="extract_reads",
            test="sample_distribution_shift",
        )


def test_validate_region_contrast_request_rejects_single_read_unknown_test():
    with pytest.raises(ValueError, match="sample_distribution_shift"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="extract_reads",
            test="beta_binomial",
        )


def test_validate_rejects_unsupported_single_read_beta_binomial():
    with pytest.raises(ValueError, match="extract_reads"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="pileup_counts",
            test="beta_binomial",
        )


def _mock_region_summaries():
    return pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "cluster": "C1",
                "count": 6,
                "fraction": 0.6,
            },
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "cluster": "C2",
                "count": 4,
                "fraction": 0.4,
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "15min",
                "cluster": "C1",
                "count": 1,
                "fraction": 0.2,
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "15min",
                "cluster": "C2",
                "count": 4,
                "fraction": 0.8,
            },
            {
                "region_id": "reg2",
                "sample_id": "s1",
                "condition": "NS",
                "cluster": "C3",
                "count": 5,
                "fraction": 1.0,
            },
        ]
    )


def _mock_cluster_occupancy_evidence():
    region_summaries = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "ns1",
                "condition": "NS",
                "cluster": "C1",
                "count": 2,
            },
            {
                "region_id": "reg1",
                "sample_id": "ns1",
                "condition": "NS",
                "cluster": "C2",
                "count": 8,
            },
            {
                "region_id": "reg1",
                "sample_id": "ns2",
                "condition": "NS",
                "cluster": "C1",
                "count": 3,
            },
            {
                "region_id": "reg1",
                "sample_id": "ns2",
                "condition": "NS",
                "cluster": "C2",
                "count": 7,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx1",
                "condition": "15min",
                "cluster": "C1",
                "count": 8,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx1",
                "condition": "15min",
                "cluster": "C2",
                "count": 2,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx2",
                "condition": "15min",
                "cluster": "C1",
                "count": 7,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx2",
                "condition": "15min",
                "cluster": "C2",
                "count": 3,
            },
            {
                "region_id": "reg2",
                "sample_id": "ns1",
                "condition": "NS",
                "cluster": "C1",
                "count": 0,
            },
            {
                "region_id": "reg2",
                "sample_id": "ns1",
                "condition": "NS",
                "cluster": "C2",
                "count": 0,
            },
            {
                "region_id": "reg2",
                "sample_id": "ns2",
                "condition": "NS",
                "cluster": "C1",
                "count": 0,
            },
            {
                "region_id": "reg2",
                "sample_id": "ns2",
                "condition": "NS",
                "cluster": "C2",
                "count": 0,
            },
            {
                "region_id": "reg2",
                "sample_id": "tx1",
                "condition": "15min",
                "cluster": "C1",
                "count": 10,
            },
            {
                "region_id": "reg2",
                "sample_id": "tx1",
                "condition": "15min",
                "cluster": "C2",
                "count": 0,
            },
            {
                "region_id": "reg2",
                "sample_id": "tx2",
                "condition": "15min",
                "cluster": "C1",
                "count": 0,
            },
            {
                "region_id": "reg2",
                "sample_id": "tx2",
                "condition": "15min",
                "cluster": "C2",
                "count": 0,
            },
        ]
    )
    return region_contrasts.build_cluster_occupancy_evidence_table(
        region_summaries=region_summaries,
    )


def _mock_group_vs_group_cluster_fraction_evidence():
    return pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "na1",
                "condition": "treated_a",
                "cluster": "C1",
                "fraction": 0.2,
            },
            {
                "region_id": "reg1",
                "sample_id": "na2",
                "condition": "treated_a",
                "cluster": "C1",
                "fraction": 0.4,
            },
            {
                "region_id": "reg1",
                "sample_id": "nb1",
                "condition": "treated_b",
                "cluster": "C1",
                "fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "sample_id": "nb2",
                "condition": "treated_b",
                "cluster": "C1",
                "fraction": 1.0,
            },
            {
                "region_id": "reg1",
                "sample_id": "ca1",
                "condition": "control_a",
                "cluster": "C1",
                "fraction": 0.0,
            },
            {
                "region_id": "reg1",
                "sample_id": "ca2",
                "condition": "control_a",
                "cluster": "C1",
                "fraction": 0.2,
            },
            {
                "region_id": "reg1",
                "sample_id": "cb1",
                "condition": "control_b",
                "cluster": "C1",
                "fraction": 0.4,
            },
            {
                "region_id": "reg1",
                "sample_id": "cb2",
                "condition": "control_b",
                "cluster": "C1",
                "fraction": 0.6,
            },
        ]
    )


def test_build_cluster_occupancy_evidence_table_summarizes_region_sample_clusters():
    evidence = region_contrasts.build_cluster_occupancy_evidence_table(
        region_summaries=_mock_region_summaries(),
    )

    assert {
        "region_id",
        "sample_id",
        "condition",
        "cluster",
        "fraction",
        "dominant_cluster",
        "cluster_entropy",
    } <= set(evidence.columns)


def test_build_cluster_occupancy_evidence_table_computes_fraction_dominant_and_entropy():
    evidence = _mock_cluster_occupancy_evidence()

    reg1_ns1_c1 = evidence[
        (evidence["region_id"] == "reg1")
        & (evidence["sample_id"] == "ns1")
        & (evidence["cluster"] == "C1")
    ].iloc[0]
    assert reg1_ns1_c1["fraction"] == pytest.approx(0.2)
    assert reg1_ns1_c1["dominant_cluster"] == "C2"
    assert reg1_ns1_c1["cluster_entropy"] == pytest.approx(
        -(0.2 * math.log2(0.2) + 0.8 * math.log2(0.8))
    )


def test_build_cluster_occupancy_evidence_table_zero_count_groups_stay_zero_entropy():
    evidence = _mock_cluster_occupancy_evidence()

    reg2_ns1 = evidence[
        (evidence["region_id"] == "reg2") & (evidence["sample_id"] == "ns1")
    ]
    assert list(reg2_ns1["fraction"]) == [0.0, 0.0]
    assert reg2_ns1["dominant_cluster"].isna().all()
    assert reg2_ns1["cluster_entropy"].tolist() == [0.0, 0.0]


@pytest.mark.parametrize("count", [-1, float("nan"), float("inf"), float("-inf")])
def test_build_cluster_occupancy_evidence_table_rejects_invalid_counts(count):
    with pytest.raises(ValueError, match="count values must be finite and >= 0"):
        region_contrasts.build_cluster_occupancy_evidence_table(
            region_summaries=pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1",
                        "condition": "NS",
                        "cluster": "C1",
                        "count": count,
                    }
                ]
            ),
        )


def test_build_single_read_mod_fraction_evidence_table():
    extract_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "read_id": "r1",
                "modified_count": 2,
                "valid_count": 4,
            },
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "read_id": "r2",
                "modified_count": 1,
                "valid_count": 4,
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "treated",
                "read_id": "r3",
                "modified_count": 4,
                "valid_count": 4,
            },
        ]
    )

    result = region_contrasts.build_single_read_mod_fraction_evidence_table(
        extract_table=extract_table
    )

    assert result.to_dict("records") == [
        {
            "region_id": "reg1",
            "sample_id": "s1",
            "condition": "NS",
            "read_id": "r1",
            "modified_count": 2,
            "valid_count": 4,
            "read_mod_fraction": pytest.approx(0.5),
        },
        {
            "region_id": "reg1",
            "sample_id": "s1",
            "condition": "NS",
            "read_id": "r2",
            "modified_count": 1,
            "valid_count": 4,
            "read_mod_fraction": pytest.approx(0.25),
        },
        {
            "region_id": "reg1",
            "sample_id": "s2",
            "condition": "treated",
            "read_id": "r3",
            "modified_count": 4,
            "valid_count": 4,
            "read_mod_fraction": pytest.approx(1.0),
        },
    ]


def test_build_single_read_mod_fraction_evidence_table_rejects_missing_columns():
    with pytest.raises(ValueError, match="modified_count"):
        region_contrasts.build_single_read_mod_fraction_evidence_table(
            extract_table=pd.DataFrame([{"region_id": "reg1", "sample_id": "s1"}])
        )


@pytest.mark.parametrize(
    ("modified_count", "valid_count"),
    [
        (-1, 4),
        (1, -1),
        (1.5, 4),
        (1, 4.5),
        (float("nan"), 4),
        (1, float("nan")),
        (float("inf"), 4),
        (1, float("inf")),
        (5, 4),
    ],
)
def test_build_single_read_mod_fraction_evidence_table_rejects_invalid_counts(
    modified_count,
    valid_count,
):
    with pytest.raises(ValueError, match="modified_count and valid_count"):
        region_contrasts.build_single_read_mod_fraction_evidence_table(
            extract_table=pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1",
                        "condition": "NS",
                        "read_id": "r1",
                        "modified_count": modified_count,
                        "valid_count": valid_count,
                    }
                ]
            )
        )


def test_score_regions_single_read_mod_fraction_effect_size_only():
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "read_id": "r1",
                "modified_count": 1,
                "valid_count": 4,
                "read_mod_fraction": 0.25,
            },
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "read_id": "r2",
                "modified_count": 2,
                "valid_count": 4,
                "read_mod_fraction": 0.50,
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "treated",
                "read_id": "r3",
                "modified_count": 4,
                "valid_count": 4,
                "read_mod_fraction": 1.00,
            },
            {
                "region_id": "reg1",
                "sample_id": "s3",
                "condition": "treated",
                "read_id": "r4",
                "modified_count": 3,
                "valid_count": 4,
                "read_mod_fraction": 0.75,
            },
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(
            mode="group_vs_group", numerator=["treated"], denominator=["NS"]
        ),
        analysis_unit="single_read",
        representation="read_mod_fraction",
        signal_source="extract_reads",
        test="effect_size_only",
        read_table=evidence,
    )

    row = result.summary.iloc[0]
    assert row["region_id"] == "reg1"
    assert row["sample_summary_numerator_mean"] == pytest.approx(0.875)
    assert row["sample_summary_denominator_mean"] == pytest.approx(0.375)
    assert row["delta_summary_mean"] == pytest.approx(0.5)


def test_score_regions_single_read_mod_fraction_supports_matched_pairwise():
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1_before",
                "condition": "before",
                "read_id": "r1",
                "modified_count": 0,
                "valid_count": 4,
                "read_mod_fraction": 0.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s1_after",
                "condition": "after",
                "read_id": "r2",
                "modified_count": 4,
                "valid_count": 4,
                "read_mod_fraction": 1.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_before",
                "condition": "before",
                "read_id": "r3",
                "modified_count": 4,
                "valid_count": 4,
                "read_mod_fraction": 1.0,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_after",
                "condition": "after",
                "read_id": "r4",
                "modified_count": 4,
                "valid_count": 4,
                "read_mod_fraction": 1.0,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s1_before",
                "condition": "before",
                "read_id": "r5",
                "modified_count": 0,
                "valid_count": 4,
                "read_mod_fraction": 0.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s1_before",
                "condition": "before",
                "read_id": "r6",
                "modified_count": 0,
                "valid_count": 4,
                "read_mod_fraction": 0.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_after",
                "condition": "after",
                "read_id": "r7",
                "modified_count": 4,
                "valid_count": 4,
                "read_mod_fraction": 1.0,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_after",
                "condition": "after",
                "read_id": "r8",
                "modified_count": 3,
                "valid_count": 4,
                "read_mod_fraction": 0.75,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s3_after",
                "condition": "after",
                "read_id": "r5",
                "modified_count": 4,
                "valid_count": 4,
                "read_mod_fraction": 1.0,
                "pair_id": "p3",
            },
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["after"],
            denominator=["before"],
            pairing_key="pair_id",
        ),
        analysis_unit="single_read",
        representation="read_mod_fraction",
        signal_source="extract_reads",
        test="sample_distribution_shift",
        read_table=evidence,
    )

    row = result.summary.iloc[0]
    assert row["sample_summary_numerator_mean"] == pytest.approx(0.9583333333333333)
    assert row["sample_summary_denominator_mean"] == pytest.approx(0.5)
    assert row["delta_summary_mean"] == pytest.approx(0.45833333333333326)
    assert row["numerator_replicate_n"] == 2
    assert row["denominator_replicate_n"] == 2


def test_build_single_read_feature_evidence_table_accepts_user_features():
    feature_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "read_id": "r1",
                "f0": 0.1,
                "f1": 0.2,
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "treated",
                "read_id": "r2",
                "f0": 0.8,
                "f1": 0.9,
            },
        ]
    )

    result = region_contrasts.build_single_read_feature_evidence_table(
        feature_table=feature_table
    )

    assert result.to_dict("records") == feature_table.to_dict("records")


def test_build_single_read_feature_evidence_table_rejects_missing_read_id():
    with pytest.raises(ValueError, match="read_id"):
        region_contrasts.build_single_read_feature_evidence_table(
            feature_table=pd.DataFrame(
                [{"region_id": "reg1", "sample_id": "s1", "condition": "NS", "f0": 0.1}]
            )
        )


def test_build_single_read_feature_evidence_table_rejects_non_numeric_features():
    with pytest.raises(ValueError, match="numeric"):
        region_contrasts.build_single_read_feature_evidence_table(
            feature_table=pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1",
                        "condition": "NS",
                        "read_id": "r1",
                        "f0": "high",
                    }
                ]
            )
        )


def test_score_regions_single_read_window_features_effect_size_only():
    feature_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "read_id": "r1",
                "f0": 0.1,
                "f1": 0.2,
            },
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "read_id": "r2",
                "f0": 0.2,
                "f1": 0.3,
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "treated",
                "read_id": "r3",
                "f0": 0.8,
                "f1": 0.9,
            },
            {
                "region_id": "reg1",
                "sample_id": "s3",
                "condition": "treated",
                "read_id": "r4",
                "f0": 0.7,
                "f1": 0.8,
            },
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(
            mode="group_vs_group", numerator=["treated"], denominator=["NS"]
        ),
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="effect_size_only",
        feature_table=feature_table,
    )

    row = result.summary.iloc[0]
    assert row["region_id"] == "reg1"
    assert row["f0_delta_mean"] == pytest.approx(0.6)
    assert row["f1_delta_mean"] == pytest.approx(0.6)


def test_score_regions_single_read_window_features_supports_matched_pairwise():
    feature_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1_before",
                "condition": "before",
                "read_id": "r1",
                "f0": 0.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s1_after",
                "condition": "after",
                "read_id": "r2",
                "f0": 1.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_before",
                "condition": "before",
                "read_id": "r3",
                "f0": 1.0,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_after",
                "condition": "after",
                "read_id": "r4",
                "f0": 0.6,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s1_before",
                "condition": "before",
                "read_id": "r5",
                "f0": 0.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s1_before",
                "condition": "before",
                "read_id": "r6",
                "f0": 0.0,
                "pair_id": "p1",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_after",
                "condition": "after",
                "read_id": "r7",
                "f0": 0.6,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2_after",
                "condition": "after",
                "read_id": "r8",
                "f0": 0.6,
                "pair_id": "p2",
            },
            {
                "region_id": "reg1",
                "sample_id": "s3_after",
                "condition": "after",
                "read_id": "r9",
                "f0": 0.5,
                "pair_id": "p3",
            },
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["after"],
            denominator=["before"],
            pairing_key="pair_id",
        ),
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="feature_summary_shift",
        feature_table=feature_table,
    )

    row = result.summary.iloc[0]
    assert row["f0_numerator_mean"] == pytest.approx(0.8)
    assert row["f0_denominator_mean"] == pytest.approx(0.5)
    assert row["f0_delta_mean"] == pytest.approx(0.3)


@pytest.mark.parametrize(
    "representation,table_factory",
    [
        ("read_mod_fraction", "read_table"),
        ("read_window_features", "feature_table"),
    ],
)
def test_score_regions_single_read_matched_pairwise_requires_pairing_key_column(
    representation,
    table_factory,
):
    contrast = ContrastSpec(
        mode="matched_pairwise",
        numerator=["after"],
        denominator=["before"],
        pairing_key="pair_id",
    )
    if table_factory == "read_table":
        evidence_kwargs = {
            "read_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "modified_count": 1,
                        "valid_count": 4,
                        "read_mod_fraction": 0.25,
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r2",
                        "modified_count": 4,
                        "valid_count": 4,
                        "read_mod_fraction": 1.0,
                    },
                ]
            )
        }
        signal_source = "extract_reads"
        test_name = "sample_distribution_shift"
    else:
        evidence_kwargs = {
            "feature_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "f0": 0.1,
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r2",
                        "f0": 0.9,
                    },
                ]
            )
        }
        signal_source = "extract_features"
        test_name = "feature_summary_shift"

    with pytest.raises(ValueError, match="pair_id"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=contrast,
            analysis_unit="single_read",
            representation=representation,
            signal_source=signal_source,
            test=test_name,
            **evidence_kwargs,
        )


@pytest.mark.parametrize(
    "representation,table_factory,signal_source,test_name",
    [
        (
            "read_mod_fraction",
            "read_table",
            "extract_reads",
            "sample_distribution_shift",
        ),
        (
            "read_window_features",
            "feature_table",
            "extract_features",
            "feature_summary_shift",
        ),
    ],
)
def test_score_regions_single_read_matched_pairwise_requires_non_null_pairing_keys(
    representation,
    table_factory,
    signal_source,
    test_name,
):
    contrast = ContrastSpec(
        mode="matched_pairwise",
        numerator=["after"],
        denominator=["before"],
        pairing_key="pair_id",
    )
    if table_factory == "read_table":
        evidence_kwargs = {
            "read_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "modified_count": 1,
                        "valid_count": 4,
                        "read_mod_fraction": 0.25,
                        "pair_id": None,
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r2",
                        "modified_count": 4,
                        "valid_count": 4,
                        "read_mod_fraction": 1.0,
                        "pair_id": "p1",
                    },
                ]
            )
        }
    else:
        evidence_kwargs = {
            "feature_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "f0": 0.1,
                        "pair_id": None,
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r2",
                        "f0": 0.9,
                        "pair_id": "p1",
                    },
                ]
            )
        }

    with pytest.raises(ValueError, match="non-null"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=contrast,
            analysis_unit="single_read",
            representation=representation,
            signal_source=signal_source,
            test=test_name,
            **evidence_kwargs,
        )


@pytest.mark.parametrize(
    "representation,table_factory,signal_source,test_name",
    [
        (
            "read_mod_fraction",
            "read_table",
            "extract_reads",
            "sample_distribution_shift",
        ),
        (
            "read_window_features",
            "feature_table",
            "extract_features",
            "feature_summary_shift",
        ),
    ],
)
def test_score_regions_single_read_matched_pairwise_rejects_multi_condition_sides(
    representation,
    table_factory,
    signal_source,
    test_name,
):
    contrast = ContrastSpec(
        mode="matched_pairwise",
        numerator=["after", "later"],
        denominator=["before"],
        pairing_key="pair_id",
    )
    if table_factory == "read_table":
        evidence_kwargs = {
            "read_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "modified_count": 1,
                        "valid_count": 4,
                        "read_mod_fraction": 0.25,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r2",
                        "modified_count": 4,
                        "valid_count": 4,
                        "read_mod_fraction": 1.0,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_later",
                        "condition": "later",
                        "read_id": "r3",
                        "modified_count": 3,
                        "valid_count": 4,
                        "read_mod_fraction": 0.75,
                        "pair_id": "p1",
                    },
                ]
            )
        }
    else:
        evidence_kwargs = {
            "feature_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "f0": 0.1,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r2",
                        "f0": 0.9,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_later",
                        "condition": "later",
                        "read_id": "r3",
                        "f0": 0.7,
                        "pair_id": "p1",
                    },
                ]
            )
        }

    with pytest.raises(ValueError, match="exactly one numerator and one denominator"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=contrast,
            analysis_unit="single_read",
            representation=representation,
            signal_source=signal_source,
            test=test_name,
            **evidence_kwargs,
        )


@pytest.mark.parametrize(
    "representation,table_factory,signal_source,test_name",
    [
        (
            "read_mod_fraction",
            "read_table",
            "extract_reads",
            "sample_distribution_shift",
        ),
        (
            "read_window_features",
            "feature_table",
            "extract_features",
            "feature_summary_shift",
        ),
    ],
)
def test_score_regions_single_read_matched_pairwise_rejects_multiple_samples_on_one_side(
    representation,
    table_factory,
    signal_source,
    test_name,
):
    contrast = ContrastSpec(
        mode="matched_pairwise",
        numerator=["after"],
        denominator=["before"],
        pairing_key="pair_id",
    )
    if table_factory == "read_table":
        evidence_kwargs = {
            "read_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "modified_count": 1,
                        "valid_count": 4,
                        "read_mod_fraction": 0.25,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after_a",
                        "condition": "after",
                        "read_id": "r2",
                        "modified_count": 4,
                        "valid_count": 4,
                        "read_mod_fraction": 1.0,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after_b",
                        "condition": "after",
                        "read_id": "r3",
                        "modified_count": 3,
                        "valid_count": 4,
                        "read_mod_fraction": 0.75,
                        "pair_id": "p1",
                    },
                ]
            )
        }
    else:
        evidence_kwargs = {
            "feature_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_before",
                        "condition": "before",
                        "read_id": "r1",
                        "f0": 0.1,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after_a",
                        "condition": "after",
                        "read_id": "r2",
                        "f0": 0.9,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after_b",
                        "condition": "after",
                        "read_id": "r3",
                        "f0": 0.7,
                        "pair_id": "p1",
                    },
                ]
            )
        }

    with pytest.raises(
        ValueError, match="exactly one sample per region, pair, and condition"
    ):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=contrast,
            analysis_unit="single_read",
            representation=representation,
            signal_source=signal_source,
            test=test_name,
            **evidence_kwargs,
        )


@pytest.mark.parametrize(
    "representation,table_factory,signal_source,test_name",
    [
        (
            "read_mod_fraction",
            "read_table",
            "extract_reads",
            "sample_distribution_shift",
        ),
        (
            "read_window_features",
            "feature_table",
            "extract_features",
            "feature_summary_shift",
        ),
    ],
)
def test_score_regions_single_read_matched_pairwise_rejects_sample_ids_mapped_to_multiple_pairs(
    representation,
    table_factory,
    signal_source,
    test_name,
):
    contrast = ContrastSpec(
        mode="matched_pairwise",
        numerator=["after"],
        denominator=["before"],
        pairing_key="pair_id",
    )
    if table_factory == "read_table":
        evidence_kwargs = {
            "read_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "shared_sample",
                        "condition": "before",
                        "read_id": "r1",
                        "modified_count": 1,
                        "valid_count": 4,
                        "read_mod_fraction": 0.25,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "shared_sample",
                        "condition": "before",
                        "read_id": "r2",
                        "modified_count": 2,
                        "valid_count": 4,
                        "read_mod_fraction": 0.5,
                        "pair_id": "p2",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r3",
                        "modified_count": 4,
                        "valid_count": 4,
                        "read_mod_fraction": 1.0,
                        "pair_id": "p1",
                    },
                ]
            )
        }
    else:
        evidence_kwargs = {
            "feature_table": pd.DataFrame(
                [
                    {
                        "region_id": "reg1",
                        "sample_id": "shared_sample",
                        "condition": "before",
                        "read_id": "r1",
                        "f0": 0.1,
                        "pair_id": "p1",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "shared_sample",
                        "condition": "before",
                        "read_id": "r2",
                        "f0": 0.2,
                        "pair_id": "p2",
                    },
                    {
                        "region_id": "reg1",
                        "sample_id": "s1_after",
                        "condition": "after",
                        "read_id": "r3",
                        "f0": 0.9,
                        "pair_id": "p1",
                    },
                ]
            )
        }

    with pytest.raises(ValueError, match="sample_id to map to exactly one pairing key"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=contrast,
            analysis_unit="single_read",
            representation=representation,
            signal_source=signal_source,
            test=test_name,
            **evidence_kwargs,
        )


def test_score_regions_single_read_window_features_uses_builtin_loader(monkeypatch):
    extracted_by_path = {
        "ns.h5": SimpleNamespace(
            metadata=[
                {
                    "read_name": "r1",
                    "chromosome": "chr1",
                    "region_start": 0,
                    "region_end": 10,
                    "region_strand": "+",
                }
            ]
        ),
        "treated.h5": SimpleNamespace(
            metadata=[
                {
                    "read_name": "r2",
                    "chromosome": "chr1",
                    "region_start": 0,
                    "region_end": 10,
                    "region_strand": "+",
                }
            ]
        ),
    }

    feature_rows_by_path = {
        "ns.h5": ([[0.1]], ["f0"]),
        "treated.h5": ([[0.9]], ["f0"]),
    }

    monkeypatch.setattr(
        region_contrasts,
        "cluster",
        SimpleNamespace(
            extract_read_windows=lambda **kwargs: extracted_by_path[
                kwargs["hdf5_file"]
            ],
            read_window_feature_matrix=lambda extracted: feature_rows_by_path[
                "ns.h5" if extracted is extracted_by_path["ns.h5"] else "treated.h5"
            ],
        ),
        raising=False,
    )

    result = region_contrasts.score_regions(
        samples=[
            SampleSpec(sample_id="s1", condition="NS", extract_h5="ns.h5"),
            SampleSpec(sample_id="s2", condition="treated", extract_h5="treated.h5"),
        ],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=ContrastSpec(
            mode="pairwise", numerator=["treated"], denominator=["NS"]
        ),
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="feature_summary_shift",
    )

    assert result.summary.iloc[0]["f0_delta_mean"] == pytest.approx(0.8)


def test_build_region_evidence_table_from_pileup_counts(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="s2.h5",
            replicate=2,
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]

    monkeypatch.setattr(
        region_contrasts.utils,
        "regions_dict_from_input",
        lambda regions, window_size=None: {
            "chr1": [(0, 10, "+")],
            "chr2": [(20, 30, "-")],
        },
    )

    counts_by_file = {
        "s1.bed.gz": [(2, 10), (0, 0)],
        "s2.bed.gz": [(7, 10), (8, 10)],
    }

    def fake_regions_to_list(
        function_handle,
        regions,
        window_size,
        quiet,
        cores,
        split_large_regions=False,
    ):
        return counts_by_file[function_handle.keywords["bedmethyl_file"]]

    monkeypatch.setattr(
        region_contrasts.load_processed,
        "regions_to_list",
        fake_regions_to_list,
    )

    evidence = region_contrasts.build_region_evidence_table(
        samples=samples,
        regions="matched.bed",
        motifs=["A,0"],
    )

    expected = pd.DataFrame(
        [
            {
                "region_id": "chr1:0-10,+",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "chr2:20-30,-",
                "chromosome": "chr2",
                "start": 20,
                "end": 30,
                "strand": "-",
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "modified_count": 0,
                "valid_count": 0,
                "mod_fraction": 0.0,
            },
            {
                "region_id": "chr1:0-10,+",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "s2",
                "condition": "15min",
                "replicate": 2,
                "modified_count": 7,
                "valid_count": 10,
                "mod_fraction": 0.7,
            },
            {
                "region_id": "chr2:20-30,-",
                "chromosome": "chr2",
                "start": 20,
                "end": 30,
                "strand": "-",
                "sample_id": "s2",
                "condition": "15min",
                "replicate": 2,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
        ]
    )

    pd.testing.assert_frame_equal(evidence.reset_index(drop=True), expected)


def test_score_regions_effect_size_only_pairwise_ranks_largest_delta_first(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 6,
                "valid_count": 10,
                "mod_fraction": 0.6,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
    )

    expected_summary = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "fraction": 0.8,
                "reference_fraction": 0.2,
                "delta_fraction": 0.6,
                "rank": 1,
            },
            {
                "region_id": "reg2",
                "fraction": 0.6,
                "reference_fraction": 0.4,
                "delta_fraction": 0.2,
                "rank": 2,
            },
        ]
    )

    actual_summary = result.summary[
        [
            "region_id",
            "fraction",
            "reference_fraction",
            "delta_fraction",
            "log2_fc",
            "rank",
        ]
    ].reset_index(drop=True)

    pd.testing.assert_frame_equal(
        actual_summary.drop(columns=["log2_fc"]),
        expected_summary,
    )
    assert actual_summary.loc[0, "log2_fc"] == pytest.approx(1.99999459)
    assert actual_summary.loc[1, "log2_fc"] == pytest.approx(0.58496124)
    assert list(result.regions["region_id"]) == ["reg1", "reg2"]
    assert set(result.plot_data) == {"region_effect_sizes"}
    assert result.metadata["contrast_mode"] == "pairwise"


def test_score_regions_defaults_to_effect_size_only(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
    )

    assert result.summary.loc[0, "region_id"] == "reg1"
    assert result.summary.loc[0, "delta_fraction"] == pytest.approx(0.6)
    assert result.metadata["test"] == "effect_size_only"


def test_score_regions_effect_size_only_group_vs_group_pools_conditions(monkeypatch):
    contrast = ContrastSpec(
        mode="group_vs_group",
        numerator=["treated_a", "treated_b"],
        denominator=["control_a", "control_b"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated_a",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n2",
                "condition": "treated_b",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control_a",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d2",
                "condition": "control_b",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
    )

    summary_row = result.summary.iloc[0]
    assert summary_row["region_id"] == "reg1"
    assert summary_row["fraction"] == pytest.approx(0.6)
    assert summary_row["reference_fraction"] == pytest.approx(0.3)
    assert summary_row["delta_fraction"] == pytest.approx(0.3)
    assert summary_row["numerator_replicate_n"] == 2
    assert summary_row["denominator_replicate_n"] == 2


def test_score_regions_beta_binomial_adds_p_values(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 9,
                "valid_count": 10,
                "mod_fraction": 0.9,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 5,
                "valid_count": 10,
                "mod_fraction": 0.5,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="beta_binomial",
    )

    for table in (result.regions, result.summary):
        assert "p_value" in table.columns
        assert "adjusted_p_value" in table.columns
        assert ((table["p_value"] >= 0) & (table["p_value"] <= 1)).all()
        assert (
            (table["adjusted_p_value"] >= 0) & (table["adjusted_p_value"] <= 1)
        ).all()


def test_score_regions_beta_binomial_rejects_unsupported_multiple_testing(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 9,
                "valid_count": 10,
                "mod_fraction": 0.9,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    with pytest.raises(ValueError, match="fdr_bh"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
            test="beta_binomial",
            multiple_testing="bonferroni",
        )


def test_score_regions_beta_binomial_uses_row_specific_denominator_counts(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 5,
                "valid_count": 10,
                "mod_fraction": 0.5,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 1,
                "valid_count": 10,
                "mod_fraction": 0.1,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 5,
                "valid_count": 10,
                "mod_fraction": 0.5,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="beta_binomial",
    )

    row1 = result.summary.loc[result.summary["region_id"] == "reg1", "p_value"].iloc[0]
    row2 = result.summary.loc[result.summary["region_id"] == "reg2", "p_value"].iloc[0]
    assert row1 != row2


def test_score_regions_beta_binomial_ranks_by_adjusted_p_value(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "high_effect_low_significance",
                "chromosome": "chr1",
                "start": 0,
                "end": 1,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 1,
                "valid_count": 1,
                "mod_fraction": 1.0,
            },
            {
                "region_id": "high_effect_low_significance",
                "chromosome": "chr1",
                "start": 0,
                "end": 1,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 0,
                "valid_count": 1,
                "mod_fraction": 0.0,
            },
            {
                "region_id": "lower_effect_higher_significance",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 18,
                "valid_count": 20,
                "mod_fraction": 0.9,
            },
            {
                "region_id": "lower_effect_higher_significance",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 10,
                "valid_count": 20,
                "mod_fraction": 0.5,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="beta_binomial",
    )

    assert result.summary.iloc[0]["region_id"] == "lower_effect_higher_significance"
    assert (
        result.summary.iloc[0]["adjusted_p_value"]
        <= result.summary.iloc[1]["adjusted_p_value"]
    )
    assert (
        result.summary.iloc[0]["delta_fraction"]
        < result.summary.iloc[1]["delta_fraction"]
    )
    assert result.summary.iloc[0]["rank"] == 1


def test_score_regions_beta_binomial_uses_p_value_tiebreak_when_adjusted_values_match(
    monkeypatch,
):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 5,
                "valid_count": 10,
                "mod_fraction": 0.5,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    original_add = region_contrasts._add_beta_binomial_scores

    def fake_add_beta_binomial_scores(regions_table, *, multiple_testing):
        scored = original_add(regions_table, multiple_testing=multiple_testing)
        scored.loc[scored["region_id"] == "reg1", "p_value"] = 0.02
        scored.loc[scored["region_id"] == "reg2", "p_value"] = 0.01
        scored["adjusted_p_value"] = 0.05
        return scored

    monkeypatch.setattr(
        region_contrasts,
        "_add_beta_binomial_scores",
        fake_add_beta_binomial_scores,
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="beta_binomial",
    )

    assert list(result.summary["region_id"]) == ["reg2", "reg1"]
    assert list(result.summary["rank"]) == [1, 2]


def test_beta_binomial_two_sided_p_value_rejects_invalid_counts():
    with pytest.raises(ValueError, match="modified_count"):
        region_contrasts._beta_binomial_two_sided_p_value(6, 5, 2.0, 2.0)


@pytest.mark.parametrize(
    ("denominator_modified_count", "denominator_valid_count"),
    [
        (-1, 5),
        (3, -1),
        (6, 5),
    ],
)
def test_estimate_beta_binomial_prior_rejects_invalid_denominator_counts(
    denominator_modified_count, denominator_valid_count
):
    with pytest.raises(ValueError):
        region_contrasts._estimate_beta_binomial_prior(
            denominator_modified_count,
            denominator_valid_count,
        )


@pytest.mark.parametrize(
    ("modified_count", "valid_count"),
    [
        (-1, 5),
        (1, -1),
    ],
)
def test_beta_binomial_two_sided_p_value_rejects_negative_counts(
    modified_count, valid_count
):
    with pytest.raises(ValueError):
        region_contrasts._beta_binomial_two_sided_p_value(
            modified_count,
            valid_count,
            2.0,
            2.0,
        )


def test_score_regions_rejects_missing_denominator_condition(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["cntrol"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    with pytest.raises(ValueError, match="denominator.*cntrol"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
        )


def test_score_regions_rejects_missing_numerator_condition(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["trated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    with pytest.raises(ValueError, match="numerator.*trated"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
        )


def test_score_regions_modified_count_representation_ranks_by_delta_count(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 11,
                "valid_count": 20,
                "mod_fraction": 0.55,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 1,
                "valid_count": 10,
                "mod_fraction": 0.1,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 7,
                "valid_count": 10,
                "mod_fraction": 0.7,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 0,
                "valid_count": 5,
                "mod_fraction": 0.0,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        representation="modified_count",
    )

    assert list(result.summary["region_id"]) == ["reg1", "reg2"]
    assert list(result.summary["rank"]) == [1, 2]
    assert result.summary.loc[0, "count"] == 11
    assert result.summary.loc[0, "reference_count"] == 1
    assert result.summary.loc[0, "delta_count"] == 10
    assert result.summary.loc[1, "delta_count"] == 7
    assert result.summary.loc[0, "log2_fc_count"] == pytest.approx(
        math.log2((11 + 1e-6) / (1 + 1e-6))
    )
    assert result.metadata["representation"] == "modified_count"


def test_score_regions_effect_size_only_rejects_unsupported_contrast_mode(monkeypatch):
    contrast = ContrastSpec(
        mode="matched_pairwise",
        numerator=["treated"],
        denominator=["control"],
        pairing_key="donor_id",
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: pd.DataFrame(),
    )

    with pytest.raises(NotImplementedError, match="matched_pairwise"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
            test="effect_size_only",
        )


def test_score_regions_cluster_fraction_effect_size_only_ranks_largest_fraction_shift_first():
    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(mode="pairwise", numerator=["15min"], denominator=["NS"]),
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="effect_size_only",
        occupancy_table=_mock_cluster_occupancy_evidence(),
    )

    first_row = result.regions.iloc[0]
    assert first_row["region_id"] == "reg1"
    assert first_row["cluster"] == "C1"
    assert first_row["fraction"] == pytest.approx(0.75)
    assert first_row["reference_fraction"] == pytest.approx(0.25)
    assert first_row["delta_fraction"] == pytest.approx(0.5)
    assert first_row["numerator_replicate_n"] == 2
    assert first_row["denominator_replicate_n"] == 2
    assert result.summary.iloc[0]["delta_fraction"] == pytest.approx(0.5)


def test_score_regions_cluster_fraction_fraction_test_adds_p_values():
    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(mode="pairwise", numerator=["15min"], denominator=["NS"]),
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="fraction_test",
        occupancy_table=_mock_cluster_occupancy_evidence(),
    )

    assert {"p_value", "adjusted_p_value"} <= set(result.regions.columns)
    assert ((result.regions["p_value"] >= 0) & (result.regions["p_value"] <= 1)).all()
    assert (
        (result.regions["adjusted_p_value"] >= 0)
        & (result.regions["adjusted_p_value"] <= 1)
    ).all()


def test_score_regions_cluster_occupancy_dominant_cluster_returns_descriptive_summary_only():
    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(mode="pairwise", numerator=["15min"], denominator=["NS"]),
        analysis_unit="cluster_occupancy",
        representation="dominant_cluster",
        signal_source="cluster_occupancy",
        test="effect_size_only",
        occupancy_table=_mock_cluster_occupancy_evidence(),
    )

    reg1 = result.summary.loc[result.summary["region_id"] == "reg1"].iloc[0]
    reg2 = result.summary.loc[result.summary["region_id"] == "reg2"].iloc[0]
    assert reg1["dominant_cluster"] == "C1"
    assert reg1["reference_dominant_cluster"] == "C2"
    assert reg2["dominant_cluster"] == "C1"
    assert pd.isna(reg2["reference_dominant_cluster"])
    assert reg2["dominant_cluster_changed"]
    assert "p_value" not in result.regions.columns


def test_score_regions_cluster_occupancy_dominant_cluster_rejects_missing_requested_conditions():
    occupancy_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "tx1",
                "condition": "15min",
                "dominant_cluster": "C1",
            }
        ]
    )

    with pytest.raises(ValueError, match="Missing denominator evidence"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=ContrastSpec(
                mode="pairwise", numerator=["15min"], denominator=["NS"]
            ),
            analysis_unit="cluster_occupancy",
            representation="dominant_cluster",
            signal_source="cluster_occupancy",
            test="effect_size_only",
            occupancy_table=occupancy_table,
        )


def test_score_regions_cluster_occupancy_dominant_cluster_rejects_region_specific_missing_side():
    occupancy_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "tx1",
                "condition": "15min",
                "dominant_cluster": "C1",
            },
            {
                "region_id": "reg2",
                "sample_id": "ns1",
                "condition": "NS",
                "dominant_cluster": "C2",
            },
        ]
    )

    with pytest.raises(ValueError, match="Missing dominant_cluster evidence"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=ContrastSpec(
                mode="pairwise", numerator=["15min"], denominator=["NS"]
            ),
            analysis_unit="cluster_occupancy",
            representation="dominant_cluster",
            signal_source="cluster_occupancy",
            test="effect_size_only",
            occupancy_table=occupancy_table,
        )


def test_score_regions_cluster_occupancy_cluster_entropy_returns_descriptive_summary_only():
    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(mode="pairwise", numerator=["15min"], denominator=["NS"]),
        analysis_unit="cluster_occupancy",
        representation="cluster_entropy",
        signal_source="cluster_occupancy",
        test="effect_size_only",
        occupancy_table=_mock_cluster_occupancy_evidence(),
    )

    reg1 = result.summary.loc[result.summary["region_id"] == "reg1"].iloc[0]
    reg2 = result.summary.loc[result.summary["region_id"] == "reg2"].iloc[0]
    expected_entropy = (
        -(
            0.8 * math.log2(0.8)
            + 0.2 * math.log2(0.2)
            + 0.7 * math.log2(0.7)
            + 0.3 * math.log2(0.3)
        )
        / 2.0
    )
    assert reg1["cluster_entropy"] == pytest.approx(expected_entropy)
    assert reg1["reference_cluster_entropy"] == pytest.approx(expected_entropy)
    assert reg1["delta_cluster_entropy"] == pytest.approx(0.0)
    assert reg2["cluster_entropy"] == pytest.approx(0.0)
    assert reg2["reference_cluster_entropy"] == pytest.approx(0.0)
    assert "p_value" not in result.regions.columns


def test_score_regions_cluster_occupancy_rejects_matched_pairwise():
    with pytest.raises(NotImplementedError, match="matched_pairwise"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=ContrastSpec(
                mode="matched_pairwise",
                numerator=["15min"],
                denominator=["NS"],
                pairing_key="donor_id",
            ),
            analysis_unit="cluster_occupancy",
            representation="cluster_fraction",
            signal_source="cluster_occupancy",
            test="effect_size_only",
            occupancy_table=_mock_cluster_occupancy_evidence(),
        )


@pytest.mark.parametrize(
    "representation",
    ["cluster_fraction", "dominant_cluster", "cluster_entropy"],
)
def test_score_regions_cluster_occupancy_rejects_unsupported_contrast_modes_for_all_representations(
    representation,
):
    with pytest.raises(NotImplementedError, match="matched_pairwise"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=ContrastSpec(
                mode="matched_pairwise",
                numerator=["15min"],
                denominator=["NS"],
                pairing_key="donor_id",
            ),
            analysis_unit="cluster_occupancy",
            representation=representation,
            signal_source="cluster_occupancy",
            test="effect_size_only",
            occupancy_table=_mock_cluster_occupancy_evidence(),
        )


def test_score_regions_cluster_occupancy_rejects_fraction_test_for_dominant_cluster():
    with pytest.raises(ValueError, match="cluster_fraction"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=ContrastSpec(
                mode="pairwise", numerator=["15min"], denominator=["NS"]
            ),
            analysis_unit="cluster_occupancy",
            representation="dominant_cluster",
            signal_source="cluster_occupancy",
            test="fraction_test",
            occupancy_table=_mock_cluster_occupancy_evidence(),
        )


def test_score_regions_cluster_occupancy_group_vs_group_aggregates_condition_groups():
    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(
            mode="group_vs_group",
            numerator=["treated_a", "treated_b"],
            denominator=["control_a", "control_b"],
        ),
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="effect_size_only",
        occupancy_table=_mock_group_vs_group_cluster_fraction_evidence(),
    )

    row = result.regions.iloc[0]
    assert row["fraction"] == pytest.approx(0.6)
    assert row["reference_fraction"] == pytest.approx(0.3)
    assert row["delta_fraction"] == pytest.approx(0.3)
    assert row["numerator_replicate_n"] == 4
    assert row["denominator_replicate_n"] == 4
    assert result.summary.iloc[0]["delta_fraction"] == pytest.approx(0.3)


def test_score_regions_cluster_occupancy_group_vs_group_uses_sample_weighted_dominant_cluster():
    occupancy_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "na1",
                "condition": "treated_a",
                "dominant_cluster": "C1",
            },
            {
                "region_id": "reg1",
                "sample_id": "nb1",
                "condition": "treated_b",
                "dominant_cluster": "C2",
            },
            {
                "region_id": "reg1",
                "sample_id": "nb2",
                "condition": "treated_b",
                "dominant_cluster": "C2",
            },
            {
                "region_id": "reg1",
                "sample_id": "ca1",
                "condition": "control_a",
                "dominant_cluster": "C1",
            },
            {
                "region_id": "reg1",
                "sample_id": "ca2",
                "condition": "control_a",
                "dominant_cluster": "C1",
            },
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(
            mode="group_vs_group",
            numerator=["treated_a", "treated_b"],
            denominator=["control_a"],
        ),
        analysis_unit="cluster_occupancy",
        representation="dominant_cluster",
        signal_source="cluster_occupancy",
        test="effect_size_only",
        occupancy_table=occupancy_table,
    )

    row = result.regions.iloc[0]
    assert row["dominant_cluster"] == "C2"
    assert row["reference_dominant_cluster"] == "C1"
    assert row["dominant_cluster_changed"]
    assert row["numerator_replicate_n"] == 3
    assert row["denominator_replicate_n"] == 2


def test_score_regions_cluster_fraction_treats_missing_sample_cluster_rows_as_zero():
    occupancy_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "ns1",
                "condition": "NS",
                "cluster": "C1",
                "fraction": 1.0,
            },
            {
                "region_id": "reg1",
                "sample_id": "ns2",
                "condition": "NS",
                "cluster": "C1",
                "fraction": 1.0,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx1",
                "condition": "15min",
                "cluster": "C1",
                "fraction": 1.0,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx1",
                "condition": "15min",
                "cluster": "C2",
                "fraction": 0.0,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx2",
                "condition": "15min",
                "cluster": "C1",
                "fraction": 0.0,
            },
            {
                "region_id": "reg1",
                "sample_id": "tx2",
                "condition": "15min",
                "cluster": "C2",
                "fraction": 1.0,
            },
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(mode="pairwise", numerator=["15min"], denominator=["NS"]),
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="fraction_test",
        occupancy_table=occupancy_table,
    )

    c2_row = result.regions.loc[result.regions["cluster"] == "C2"].iloc[0]
    assert c2_row["fraction"] == pytest.approx(0.5)
    assert c2_row["reference_fraction"] == pytest.approx(0.0)
    assert c2_row["numerator_replicate_n"] == 2
    assert c2_row["denominator_replicate_n"] == 2
    assert c2_row["numerator_sample_values"] == (0.0, 1.0)
    assert c2_row["denominator_sample_values"] == (0.0, 0.0)
    assert c2_row["p_value"] == pytest.approx(0.5)


@pytest.mark.parametrize(
    ("representation", "column", "value", "message"),
    [
        (
            "cluster_fraction",
            "fraction",
            1.2,
            "fraction values must be finite and between 0 and 1",
        ),
        (
            "cluster_fraction",
            "fraction",
            float("nan"),
            "fraction values must be finite and between 0 and 1",
        ),
        (
            "cluster_entropy",
            "cluster_entropy",
            float("inf"),
            "cluster_entropy values must be finite and >= 0",
        ),
        (
            "cluster_entropy",
            "cluster_entropy",
            -0.1,
            "cluster_entropy values must be finite and >= 0",
        ),
    ],
)
def test_score_regions_cluster_occupancy_rejects_invalid_numeric_values(
    representation,
    column,
    value,
    message,
):
    occupancy_table = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                column: value,
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "15min",
                column: 0.5 if column == "fraction" else 0.5,
            },
        ]
    )
    if representation == "cluster_fraction":
        occupancy_table["cluster"] = ["C1", "C1"]

    with pytest.raises(ValueError, match=message):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=ContrastSpec(
                mode="pairwise", numerator=["15min"], denominator=["NS"]
            ),
            analysis_unit="cluster_occupancy",
            representation=representation,
            signal_source="cluster_occupancy",
            test="effect_size_only",
            occupancy_table=occupancy_table,
        )


def test_score_regions_cluster_occupancy_rejects_missing_occupancy_columns():
    with pytest.raises(ValueError, match="occupancy_table requires columns"):
        region_contrasts.score_regions(
            samples=[],
            regions=None,
            motifs=[],
            contrast=ContrastSpec(
                mode="pairwise", numerator=["15min"], denominator=["NS"]
            ),
            analysis_unit="cluster_occupancy",
            representation="cluster_fraction",
            signal_source="cluster_occupancy",
            test="effect_size_only",
            occupancy_table=pd.DataFrame({"region_id": ["reg1"]}),
        )
