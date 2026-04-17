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
class SharedClusterContrastResult:
    summary: pd.DataFrame
    details: pd.DataFrame
    pairwise: pd.DataFrame | None = None
    plot_data: dict[str, pd.DataFrame | dict[str, Any]] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "summary": self.summary,
            "details": self.details,
            "plot_data": self.plot_data,
            "figures": self.figures,
            "metadata": self.metadata,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "SharedClusterContrastResult requires non-None values for: "
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
class ChipAtlasEnrichmentResult:
    request_id: str
    status: str
    results: pd.DataFrame | None
    query: dict[str, Any] = field(default_factory=dict)
    status_history: list[dict[str, Any]] = field(default_factory=list)
    submit_url: str | None = None
    status_url: str | None = None
    result_url: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.request_id:
            raise ValueError("ChipAtlasEnrichmentResult.request_id cannot be empty.")
        if not self.status:
            raise ValueError("ChipAtlasEnrichmentResult.status cannot be empty.")
        if self.query is None:
            raise ValueError("ChipAtlasEnrichmentResult.query cannot be None.")
        if self.status_history is None:
            raise ValueError("ChipAtlasEnrichmentResult.status_history cannot be None.")
        if self.metadata is None:
            raise ValueError("ChipAtlasEnrichmentResult.metadata cannot be None.")


@dataclass
class UniBindJobResult:
    job_id: str
    status: str
    job_url: str
    endpoint_url: str
    results_url: str | None = None
    download_urls: list[str] = field(default_factory=list)
    query: dict[str, Any] = field(default_factory=dict)
    status_history: list[dict[str, Any]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.job_id:
            raise ValueError("UniBindJobResult.job_id cannot be empty.")
        if not self.status:
            raise ValueError("UniBindJobResult.status cannot be empty.")
        if not self.job_url:
            raise ValueError("UniBindJobResult.job_url cannot be empty.")
        if not self.endpoint_url:
            raise ValueError("UniBindJobResult.endpoint_url cannot be empty.")
        if self.download_urls is None:
            raise ValueError("UniBindJobResult.download_urls cannot be None.")
        if self.query is None:
            raise ValueError("UniBindJobResult.query cannot be None.")
        if self.status_history is None:
            raise ValueError("UniBindJobResult.status_history cannot be None.")
        if self.metadata is None:
            raise ValueError("UniBindJobResult.metadata cannot be None.")


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


@dataclass
class ModkitDMRPairResult:
    output_path: Path
    segment_path: Path | None
    command: list[str]
    sites: pd.DataFrame
    segments: pd.DataFrame | None
    high_confidence_sites: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.output_path is None:
            raise ValueError("ModkitDMRPairResult.output_path cannot be None.")
        if self.command is None:
            raise ValueError("ModkitDMRPairResult.command cannot be None.")
        if self.sites is None:
            raise ValueError("ModkitDMRPairResult.sites cannot be None.")
        if self.high_confidence_sites is None:
            raise ValueError("ModkitDMRPairResult.high_confidence_sites cannot be None.")
        if self.metadata is None:
            raise ValueError("ModkitDMRPairResult.metadata cannot be None.")


@dataclass
class ModkitDMRMultiResult:
    out_dir: Path
    command: list[str]
    pair_files: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.out_dir is None:
            raise ValueError("ModkitDMRMultiResult.out_dir cannot be None.")
        if self.command is None:
            raise ValueError("ModkitDMRMultiResult.command cannot be None.")
        if self.pair_files is None:
            raise ValueError("ModkitDMRMultiResult.pair_files cannot be None.")
        if self.metadata is None:
            raise ValueError("ModkitDMRMultiResult.metadata cannot be None.")
