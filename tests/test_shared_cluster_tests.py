import pandas as pd
import pytest

from dimelo.models import ContrastSpec, SharedClusterModel, SharedClusterResult


def _make_shared_cluster_test_result() -> SharedClusterResult:
    return SharedClusterResult(
        model=SharedClusterModel(
            mode="read_global",
            motifs=["A,0"],
            feature_names=["f0", "f1"],
            preprocessing={"signal_normalization": "none"},
            estimator=object(),
            cluster_labels=["C0", "C1"],
            fit_metadata={"clusterer": "minibatch_kmeans", "n_clusters": 2},
        ),
        assignments=pd.DataFrame(
            {
                "sample_id": ["s1", "s2", "s3", "s4"],
                "condition": ["NS", "NS", "treated", "treated"],
                "subject_id": ["p1", "p2", "p1", "p2"],
                "cluster": ["C0", "C0", "C1", "C1"],
            }
        ),
        cluster_distribution=pd.DataFrame(
            [
                {"sample_id": "s1", "condition": "NS", "cluster": "C0", "count": 80, "fraction": 0.80},
                {"sample_id": "s1", "condition": "NS", "cluster": "C1", "count": 20, "fraction": 0.20},
                {"sample_id": "s2", "condition": "NS", "cluster": "C0", "count": 75, "fraction": 0.75},
                {"sample_id": "s2", "condition": "NS", "cluster": "C1", "count": 25, "fraction": 0.25},
                {"sample_id": "s3", "condition": "treated", "cluster": "C0", "count": 30, "fraction": 0.30},
                {"sample_id": "s3", "condition": "treated", "cluster": "C1", "count": 70, "fraction": 0.70},
                {"sample_id": "s4", "condition": "treated", "cluster": "C0", "count": 25, "fraction": 0.25},
                {"sample_id": "s4", "condition": "treated", "cluster": "C1", "count": 75, "fraction": 0.75},
            ]
        ),
        condition_distribution=pd.DataFrame(
            [
                {"condition": "NS", "cluster": "C0", "count": 155, "fraction": 0.775, "replicate_n": 2},
                {"condition": "NS", "cluster": "C1", "count": 45, "fraction": 0.225, "replicate_n": 2},
                {"condition": "treated", "cluster": "C0", "count": 55, "fraction": 0.275, "replicate_n": 2},
                {"condition": "treated", "cluster": "C1", "count": 145, "fraction": 0.725, "replicate_n": 2},
            ]
        ),
        distribution_change=None,
        cluster_profiles=pd.DataFrame(columns=["cluster", "count", "f0", "f1"]),
        region_summaries=None,
        plot_data={},
        figures={},
        metadata={},
    )


def _make_shared_cluster_time_course_result() -> SharedClusterResult:
    result = _make_shared_cluster_test_result()
    result.assignments = pd.DataFrame(
        {
            "sample_id": ["t0_a", "t1_a", "t2_a"],
            "condition": ["t0", "t1", "t2"],
            "subject_id": ["p1", "p1", "p1"],
            "cluster": ["C0", "C0", "C1"],
        }
    )
    result.cluster_distribution = pd.DataFrame(
        [
            {"sample_id": "t0_a", "condition": "t0", "cluster": "C0", "count": 80, "fraction": 0.80},
            {"sample_id": "t0_a", "condition": "t0", "cluster": "C1", "count": 20, "fraction": 0.20},
            {"sample_id": "t1_a", "condition": "t1", "cluster": "C0", "count": 55, "fraction": 0.55},
            {"sample_id": "t1_a", "condition": "t1", "cluster": "C1", "count": 45, "fraction": 0.45},
            {"sample_id": "t2_a", "condition": "t2", "cluster": "C0", "count": 25, "fraction": 0.25},
            {"sample_id": "t2_a", "condition": "t2", "cluster": "C1", "count": 75, "fraction": 0.75},
        ]
    )
    result.condition_distribution = pd.DataFrame(
        [
            {"condition": "t0", "cluster": "C0", "count": 80, "fraction": 0.80, "replicate_n": 1},
            {"condition": "t0", "cluster": "C1", "count": 20, "fraction": 0.20, "replicate_n": 1},
            {"condition": "t1", "cluster": "C0", "count": 55, "fraction": 0.55, "replicate_n": 1},
            {"condition": "t1", "cluster": "C1", "count": 45, "fraction": 0.45, "replicate_n": 1},
            {"condition": "t2", "cluster": "C0", "count": 25, "fraction": 0.25, "replicate_n": 1},
            {"condition": "t2", "cluster": "C1", "count": 75, "fraction": 0.75, "replicate_n": 1},
        ]
    )
    return result


def test_shared_cluster_tests_module_exports_entry_point():
    from dimelo import shared_cluster_tests

    assert hasattr(shared_cluster_tests, "shared_cluster_tests")


def test_shared_cluster_tests_rejects_unsupported_contrast_mode():
    from dimelo import shared_cluster_tests

    result = _make_shared_cluster_test_result()

    with pytest.raises(NotImplementedError, match="background_adjusted"):
        shared_cluster_tests.shared_cluster_tests(
            result=result,
            contrast=ContrastSpec(
                mode="background_adjusted",
                numerator=["treated"],
                denominator=["NS"],
                background=["bg"],
            ),
        )


def test_shared_cluster_tests_pairwise_returns_summary_details_and_plot_data():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
        test="permutation",
        n_permutations=50,
        random_state=7,
    )

    assert set(result.summary.columns) >= {
        "contrast_id",
        "composition_effect_size",
        "omnibus_p_value",
        "top_cluster",
    }
    assert set(result.details.columns) >= {
        "cluster",
        "fraction",
        "reference_fraction",
        "delta_fraction",
        "p_value",
        "adjusted_p_value",
    }
    assert set(result.plot_data) >= {"summary_table", "cluster_effect_table"}


def test_shared_cluster_tests_matched_pairwise_uses_contrast_pairing_key():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["treated"],
            denominator=["NS"],
            pairing_key="subject_id",
        ),
        test="permutation",
        n_permutations=50,
        random_state=7,
    )

    assert result.metadata["paired"] is True
    assert result.metadata["pairing_key"] == "subject_id"


def test_shared_cluster_tests_supports_chi_squared_screen():
    from dimelo import shared_cluster_tests

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
        test="chi_squared",
    )

    assert result.metadata["inference_level"] == "pooled_screen"
    assert result.summary.loc[0, "test"] == "chi_squared"
    assert 0.0 <= result.summary.loc[0, "omnibus_p_value"] <= 1.0


def test_shared_cluster_tests_supports_g_test_screen():
    from dimelo import shared_cluster_tests

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
        test="g_test",
    )

    assert result.metadata["inference_level"] == "pooled_screen"
    assert result.summary.loc[0, "test"] == "g_test"
    assert 0.0 <= result.summary.loc[0, "omnibus_p_value"] <= 1.0


def test_shared_cluster_tests_time_course_returns_omnibus_and_trend_outputs():
    from dimelo import shared_cluster_tests

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_time_course_result(),
        contrast=ContrastSpec(mode="time_course", time_order=["t0", "t1", "t2"]),
        test="permutation",
        n_permutations=50,
        random_state=7,
    )

    assert {"omnibus_p_value", "trend_p_value", "composition_effect_size"} <= set(
        result.summary.columns
    )
    assert "time_course_table" in result.plot_data


def test_shared_cluster_tests_time_course_optional_pairwise_follow_up():
    from dimelo import shared_cluster_tests

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_time_course_result(),
        contrast=ContrastSpec(mode="time_course", time_order=["t0", "t1", "t2"]),
        test="permutation",
        include_pairwise=True,
        n_permutations=20,
        random_state=7,
    )

    assert result.pairwise is not None
