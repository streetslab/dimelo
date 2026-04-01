from pathlib import Path

import pandas as pd
import pytest

from dimelo import models
from dimelo.models import (
    BatchJob,
    CohortSpec,
    ContrastSpec,
    GlobalAnalysisResult,
    DatasetArtifact,
    SampleSpec,
    RegionContrastResult,
    SharedClusterModel,
    SharedClusterResult,
)


def test_dimelo_package_exports_models():
    assert models.SampleSpec is SampleSpec


def test_sample_spec_fields():
    sample = SampleSpec(
        sample_id="sample-1",
        condition="treated",
        extract_h5=Path("sample-1.h5"),
        regions_bed=Path("sample-1.bed"),
        replicate=2,
        metadata={"batch": "A"},
    )

    assert sample.sample_id == "sample-1"
    assert sample.condition == "treated"
    assert sample.extract_h5 == Path("sample-1.h5")
    assert sample.regions_bed == Path("sample-1.bed")
    assert sample.replicate == 2
    assert sample.metadata == {"batch": "A"}


def test_dataset_artifact_stores_metadata():
    artifact = DatasetArtifact(
        sample_id="sample-1",
        artifact_type="extract",
        path=Path("sample-1.h5"),
        format="hdf5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"source": "parse_bam"},
    )

    assert artifact.metadata == {"source": "parse_bam"}
    assert artifact.params == {"window_size": 200}
    assert artifact.provenance == {"pipeline": "parse_bam"}


def test_dataset_artifact_rejects_none_metadata():
    with pytest.raises(ValueError, match="metadata"):
        DatasetArtifact(
            sample_id="sample-1",
            artifact_type="extract",
            path=Path("sample-1.h5"),
            format="hdf5",
            params={"window_size": 200},
            provenance={"pipeline": "parse_bam"},
            metadata=None,
        )


def test_dataset_artifact_rejects_none_provenance():
    with pytest.raises(ValueError):
        DatasetArtifact(
            sample_id="sample-1",
            artifact_type="extract",
            path=Path("sample-1.h5"),
            format="hdf5",
            params={"window_size": 200},
            provenance=None,
        )


def test_dataset_artifact_rejects_none_params():
    with pytest.raises(ValueError, match="params"):
        DatasetArtifact(
            sample_id="sample-1",
            artifact_type="extract",
            path=Path("sample-1.h5"),
            format="hdf5",
            params=None,
            provenance={"pipeline": "parse_bam"},
        )


def test_contrast_spec_accepts_pairwise_mode():
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["15min"],
        denominator=["NS"],
        reference_condition="NS",
    )

    assert contrast.mode == "pairwise"
    assert contrast.numerator == ["15min"]
    assert contrast.denominator == ["NS"]
    assert contrast.reference_condition == "NS"


def test_contrast_spec_accepts_group_vs_group_mode():
    contrast = ContrastSpec(
        mode="group_vs_group",
        numerator=["15min", "30min"],
        denominator=["NS", "0min"],
    )

    assert contrast.mode == "group_vs_group"
    assert contrast.numerator == ["15min", "30min"]
    assert contrast.denominator == ["NS", "0min"]


def test_contrast_spec_accepts_background_adjusted_mode():
    contrast = ContrastSpec(
        mode="background_adjusted",
        numerator=["15min"],
        denominator=["NS"],
        background=["bg"],
    )

    assert contrast.mode == "background_adjusted"
    assert contrast.background == ["bg"]


def test_contrast_spec_accepts_time_course_mode():
    contrast = ContrastSpec(
        mode="time_course",
        time_order=["NS", "15min", "30min"],
    )

    assert contrast.mode == "time_course"
    assert contrast.time_order == ["NS", "15min", "30min"]


def test_contrast_spec_rejects_missing_groups_for_pairwise():
    with pytest.raises(ValueError, match="numerator and denominator"):
        ContrastSpec(mode="pairwise")


def test_contrast_spec_rejects_missing_groups_for_group_vs_group():
    with pytest.raises(ValueError, match="numerator and denominator"):
        ContrastSpec(mode="group_vs_group", numerator=["15min"])


def test_contrast_spec_rejects_missing_pairing_key_for_matched_pairwise():
    with pytest.raises(ValueError, match="pairing_key"):
        ContrastSpec(
            mode="matched_pairwise",
            numerator=["15min"],
            denominator=["NS"],
        )


def test_contrast_spec_rejects_missing_groups_for_matched_pairwise():
    with pytest.raises(ValueError, match="matched_pairwise mode requires"):
        ContrastSpec(
            mode="matched_pairwise",
            pairing_key="pair-1",
        )


def test_contrast_spec_rejects_missing_background_for_background_adjusted():
    with pytest.raises(ValueError, match="background"):
        ContrastSpec(
            mode="background_adjusted",
            numerator=["15min"],
            denominator=["NS"],
        )


def test_contrast_spec_rejects_missing_groups_for_background_adjusted():
    with pytest.raises(ValueError, match="numerator and denominator"):
        ContrastSpec(
            mode="background_adjusted",
            background=["bg"],
        )


def test_contrast_spec_rejects_missing_time_order_for_time_course():
    with pytest.raises(ValueError, match="time_order"):
        ContrastSpec(mode="time_course")


def test_shared_cluster_result_supports_plot_data():
    model = SharedClusterModel(
        mode="shared",
        motifs=["A,0"],
        feature_names=["A,0_mod_fraction"],
        preprocessing={"scale": "standard"},
        estimator=object(),
        cluster_labels=["cluster-1"],
        fit_metadata={"random_state": 7},
    )
    result = SharedClusterResult(
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

    assert result.plot_data["cluster_distribution_bar"] == {"kind": "bar"}
    assert result.model is model


def test_shared_cluster_result_rejects_none_model():
    with pytest.raises(ValueError, match="model"):
        SharedClusterResult(
            model=None,
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


@pytest.mark.parametrize(
    "field_name",
    [
        "assignments",
        "cluster_distribution",
        "condition_distribution",
        "cluster_profiles",
        "plot_data",
    ],
)
def test_shared_cluster_result_rejects_none_core_outputs(field_name):
    model = SharedClusterModel(
        mode="shared",
        motifs=["A,0"],
        feature_names=["A,0_mod_fraction"],
        preprocessing={"scale": "standard"},
        estimator=object(),
        cluster_labels=["cluster-1"],
        fit_metadata={"random_state": 7},
    )
    kwargs = {
        "model": model,
        "assignments": pd.DataFrame({"cluster": ["cluster-1"]}),
        "cluster_distribution": pd.DataFrame({"cluster": ["cluster-1"]}),
        "condition_distribution": pd.DataFrame({"condition": ["treated"]}),
        "distribution_change": None,
        "cluster_profiles": pd.DataFrame({"profile": [1.0]}),
        "region_summaries": None,
        "plot_data": {"cluster_distribution_bar": {"kind": "bar"}},
        "figures": {},
        "metadata": {"notes": "ok"},
    }
    kwargs[field_name] = None

    with pytest.raises(ValueError):
        SharedClusterResult(**kwargs)


def test_global_analysis_result_supports_summary_windows_and_normalization():
    result = GlobalAnalysisResult(
        summary=pd.DataFrame({"sample_id": ["s1"]}),
        windows=pd.DataFrame({"window_id": ["chr1:0-1000"]}),
        normalization_factors=pd.DataFrame(
            {"sample_id": ["s1"], "global_offset": [0.1]}
        ),
        plot_data={"global_fraction_bar": pd.DataFrame({"sample_id": ["s1"]})},
        metadata={"normalization_mode": "per_sample_global"},
        figures={},
    )

    assert list(result.summary["sample_id"]) == ["s1"]
    assert list(result.windows["window_id"]) == ["chr1:0-1000"]
    assert list(result.normalization_factors["sample_id"]) == ["s1"]


def test_global_analysis_result_rejects_none_core_outputs():
    with pytest.raises(
        ValueError, match="summary, normalization_factors, plot_data"
    ):
        GlobalAnalysisResult(
            summary=None,
            windows=None,
            normalization_factors=None,
            plot_data=None,
            metadata={},
            figures={},
        )


def test_region_contrast_result_rejects_none_core_tables():
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["15min"],
        denominator=["NS"],
        reference_condition="NS",
    )

    with pytest.raises(ValueError, match="regions, summary, plot_data"):
        RegionContrastResult(
            regions=None,
            summary=None,
            contrast=contrast,
            plot_data=None,
        )


def test_region_contrast_result_rejects_none_contrast():
    with pytest.raises(ValueError, match="contrast"):
        RegionContrastResult(
            regions=pd.DataFrame({"region": ["r1"]}),
            summary=pd.DataFrame({"metric": [1.0]}),
            contrast=None,
            plot_data={"volcano": pd.DataFrame({"x": [1.0]})},
        )


def test_cohort_spec_stores_workflow_and_params():
    cohort = CohortSpec(
        cohort_id="cohort-1",
        sample_ids=["sample-1", "sample-2"],
        workflow="shared_cluster_distribution",
        params={"window_size": 200},
        metadata={"owner": "team-a"},
    )

    assert cohort.workflow == "shared_cluster_distribution"
    assert cohort.params == {"window_size": 200}
    assert cohort.sample_ids == ["sample-1", "sample-2"]
    assert cohort.metadata == {"owner": "team-a"}


def test_batch_job_groups_cohorts():
    cohort = CohortSpec(
        cohort_id="cohort-1",
        sample_ids=["sample-1"],
        workflow="shared_cluster_distribution",
        params={},
    )
    job = BatchJob(
        job_id="job-1",
        workflow="shared_cluster_distribution",
        cohorts=[cohort],
    )

    assert job.job_id == "job-1"
    assert job.workflow == "shared_cluster_distribution"
    assert job.cohorts == [cohort]
    assert job.artifact_policy == "prefer_cached"
