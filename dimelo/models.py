from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class SampleSpec:
    sample_id: str
    condition: str
    extract_h5: str | Path
    regions_bed: str | Path | None = None
    replicate: int | None = None
    metadata: dict[str, Any] | None = None


@dataclass
class DatasetArtifact:
    sample_id: str
    artifact_type: str
    path: str | Path
    format: str
    params: dict[str, Any]
    provenance: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class SharedClusterModel:
    mode: str
    motifs: list[str]
    feature_names: list[str]
    preprocessing: dict[str, Any]
    estimator: Any
    cluster_labels: list[str]
    fit_metadata: dict[str, Any]


@dataclass
class SharedClusterResult:
    model: SharedClusterModel
    assignments: pd.DataFrame | None
    cluster_distribution: pd.DataFrame | None
    condition_distribution: pd.DataFrame | None
    distribution_change: pd.DataFrame | None
    cluster_profiles: pd.DataFrame | None
    region_summaries: pd.DataFrame | None
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    figures: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class CohortSpec:
    cohort_id: str
    sample_ids: list[str]
    workflow: str
    params: dict[str, Any]
    metadata: dict[str, Any] | None = None


@dataclass
class BatchJob:
    job_id: str
    workflow: str
    cohorts: list[CohortSpec]
    artifact_policy: str = "prefer_cached"
    metadata: dict[str, Any] | None = None
