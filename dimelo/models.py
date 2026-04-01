from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


_SELECTED_REGION_CHROM_COLUMNS = {"chrom", "chromosome"}
_SELECTED_REGION_REQUIRED_COLUMNS = {"start", "end"}


def _validate_selected_regions_dataframe(
    selected_regions: pd.DataFrame,
    *,
    owner: str,
) -> None:
    if not isinstance(selected_regions, pd.DataFrame):
        raise TypeError(f"{owner}.selected_regions must be a pandas DataFrame")

    has_chrom_column = bool(_SELECTED_REGION_CHROM_COLUMNS & set(selected_regions.columns))
    has_required_columns = _SELECTED_REGION_REQUIRED_COLUMNS.issubset(
        selected_regions.columns
    )
    if not has_chrom_column or not has_required_columns:
        raise ValueError(
            f"{owner}.selected_regions must include 'start', 'end', and either "
            "'chrom' or 'chromosome'"
        )


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
    provenance: dict[str, Any]
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.params is None:
            raise ValueError("DatasetArtifact.params cannot be None")
        if self.provenance is None:
            raise ValueError("DatasetArtifact.provenance cannot be None")
        if self.metadata is None:
            raise ValueError("DatasetArtifact.metadata cannot be None")


@dataclass
class ContrastSpec:
    mode: str
    numerator: list[str] | None = None
    denominator: list[str] | None = None
    background: list[str] | None = None
    time_order: list[str] | None = None
    pairing_key: str | None = None
    reference_condition: str | None = None
    metadata: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        allowed_modes = {
            "single_dataset",
            "pairwise",
            "matched_pairwise",
            "group_vs_group",
            "background_adjusted",
            "time_course",
        }
        if self.mode not in allowed_modes:
            raise ValueError(f"Unsupported contrast mode: {self.mode}")
        if self.mode in {"pairwise", "group_vs_group"} and (
            not self.numerator or not self.denominator
        ):
            raise ValueError(
                "ContrastSpec pairwise/group_vs_group modes require numerator and denominator."
            )
        if self.mode == "matched_pairwise" and (
            not self.numerator or not self.denominator or not self.pairing_key
        ):
            raise ValueError(
                "ContrastSpec matched_pairwise mode requires numerator, "
                "denominator, and pairing_key."
            )
        if self.mode == "background_adjusted":
            if not self.numerator or not self.denominator:
                raise ValueError(
                    "ContrastSpec background_adjusted mode requires numerator and denominator."
                )
            if not self.background:
                raise ValueError(
                    "ContrastSpec background_adjusted mode requires background."
                )
        if self.mode == "time_course" and not self.time_order:
            raise ValueError("ContrastSpec time_course mode requires time_order.")


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
    assignments: pd.DataFrame
    cluster_distribution: pd.DataFrame
    condition_distribution: pd.DataFrame
    distribution_change: pd.DataFrame | None
    cluster_profiles: pd.DataFrame
    region_summaries: pd.DataFrame | None
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    figures: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "model": self.model,
            "assignments": self.assignments,
            "cluster_distribution": self.cluster_distribution,
            "condition_distribution": self.condition_distribution,
            "cluster_profiles": self.cluster_profiles,
            "plot_data": self.plot_data,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "SharedClusterResult requires non-None values for: "
                f"{', '.join(missing)}"
            )


@dataclass
class GlobalAnalysisResult:
    summary: pd.DataFrame
    windows: pd.DataFrame | None
    normalization_factors: pd.DataFrame
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "summary": self.summary,
            "normalization_factors": self.normalization_factors,
            "plot_data": self.plot_data,
            "metadata": self.metadata,
            "figures": self.figures,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "GlobalAnalysisResult requires non-None values for: "
                f"{', '.join(missing)}"
            )


@dataclass
class RegionContrastResult:
    regions: pd.DataFrame
    summary: pd.DataFrame
    contrast: ContrastSpec
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "regions": self.regions,
            "summary": self.summary,
            "contrast": self.contrast,
            "plot_data": self.plot_data,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "RegionContrastResult requires non-None values for: "
                f"{', '.join(missing)}"
            )


@dataclass
class RegionDiscoveryResult:
    hits: pd.DataFrame
    windows: pd.DataFrame
    contrast: ContrastSpec | None
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "hits": self.hits,
            "windows": self.windows,
            "plot_data": self.plot_data,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "RegionDiscoveryResult requires non-None values for: "
                f"{', '.join(missing)}"
            )


@dataclass
class RegionDiscoveryClusterContrastResult:
    discovery: RegionDiscoveryResult
    clustering: SharedClusterResult
    contrasts: RegionContrastResult
    selected_regions: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "discovery": self.discovery,
            "clustering": self.clustering,
            "contrasts": self.contrasts,
            "selected_regions": self.selected_regions,
            "metadata": self.metadata,
            "figures": self.figures,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "RegionDiscoveryClusterContrastResult requires non-None values for: "
                f"{', '.join(missing)}"
            )
        if not isinstance(self.discovery, RegionDiscoveryResult):
            raise TypeError(
                "RegionDiscoveryClusterContrastResult.discovery must be a "
                "RegionDiscoveryResult"
            )
        if not isinstance(self.clustering, SharedClusterResult):
            raise TypeError(
                "RegionDiscoveryClusterContrastResult.clustering must be a "
                "SharedClusterResult"
            )
        if not isinstance(self.contrasts, RegionContrastResult):
            raise TypeError(
                "RegionDiscoveryClusterContrastResult.contrasts must be a "
                "RegionContrastResult"
            )
        _validate_selected_regions_dataframe(
            self.selected_regions,
            owner="RegionDiscoveryClusterContrastResult",
        )


@dataclass
class RegionDiscoveryClusterResult:
    discovery: RegionDiscoveryResult
    clustering: SharedClusterResult
    selected_regions: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "discovery": self.discovery,
            "clustering": self.clustering,
            "selected_regions": self.selected_regions,
            "metadata": self.metadata,
            "figures": self.figures,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "RegionDiscoveryClusterResult requires non-None values for: "
                f"{', '.join(missing)}"
            )
        if not isinstance(self.discovery, RegionDiscoveryResult):
            raise TypeError(
                "RegionDiscoveryClusterResult.discovery must be a "
                "RegionDiscoveryResult"
            )
        if not isinstance(self.clustering, SharedClusterResult):
            raise TypeError(
                "RegionDiscoveryClusterResult.clustering must be a "
                "SharedClusterResult"
            )


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
