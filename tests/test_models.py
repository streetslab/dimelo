from pathlib import Path

import pandas as pd

from dimelo import models
from dimelo.models import (
    BatchJob,
    CohortSpec,
    DatasetArtifact,
    SampleSpec,
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
        metadata={"source": "parse_bam"},
    )

    assert artifact.metadata == {"source": "parse_bam"}
    assert artifact.params == {"window_size": 200}


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
