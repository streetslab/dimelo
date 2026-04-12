import pandas as pd
import pytest

from dimelo.models import ContrastSpec, SharedClusterModel, SharedClusterResult


def _make_shared_cluster_test_result() -> SharedClusterResult:
    model = SharedClusterModel(
        mode="shared",
        motifs=["A,0"],
        feature_names=["A,0_mod_fraction"],
        preprocessing={"scale": "standard"},
        estimator=object(),
        cluster_labels=["cluster-1"],
        fit_metadata={"random_state": 7},
    )
    return SharedClusterResult(
        model=model,
        assignments=pd.DataFrame({"cluster": ["cluster-1"]}),
        cluster_distribution=pd.DataFrame({"cluster": ["cluster-1"]}),
        condition_distribution=pd.DataFrame({"condition": ["treated"]}),
        distribution_change=None,
        cluster_profiles=pd.DataFrame({"profile": [1.0]}),
        region_summaries=None,
        plot_data={"cluster_distribution_bar": {"kind": "bar"}},
        figures={},
        metadata={"notes": "ok"},
    )


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
