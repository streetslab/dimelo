from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from . import (
    chip_atlas,
    cluster,
    distribution,
    dmr,
    plotting,
    region_analysis,
    region_contrasts,
    region_discovery,
    regulatory_enrichment,
)
from .artifacts import resolve_artifact
from .models import (
    ChipAtlasEnrichmentResult,
    DatasetArtifact,
    ModkitDMRMultiResult,
    ModkitDMRPairResult,
    RegionDiscoveryClusterContrastResult,
    RegionDiscoveryClusterResult,
    SampleSpec,
    SharedClusterModel,
    SharedClusterResult,
    UniBindJobResult,
)

_LOGGER = logging.getLogger(__name__)

_SUPPORTED_SIGNAL_NORMALIZATION = {"none", "per_sample_global", "control_regions"}
_SUPPORTED_FEATURE_SCALING = {"none", "robust_zscore"}
_SUPPORTED_CLUSTER_BASIS = {"shape_only", "shape_plus_level", "level_only"}
_PACKAGE_VERSION = "1.0.0"
_PREDICTION_CHUNK_SIZE = 100_000
_DISCOVERY_SELECTION_DEFAULT_TOP_N = 250
_LEVEL_FEATURES = {
    "global_mean",
    "global_var",
    "global_median",
    "q25",
    "q75",
    "iqr",
    "global_mod_fraction",
}
RegionSpec = str | Path | list[str | Path] | None


def _source_fingerprint(path: Path) -> dict[str, Any]:
    if path.exists():
        stat = path.stat()
        return {
            "path": str(path),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        }
    return {"path": str(path), "missing": True}


def _serialize_region_spec(region_spec: Any) -> Any:
    if region_spec is None:
        return None
    if isinstance(region_spec, (str, Path)):
        return str(region_spec)
    if isinstance(region_spec, list):
        return [str(item) if isinstance(item, Path) else item for item in region_spec]
    return region_spec


def _region_spec_content_key(region_spec: Any) -> Any:
    """Normalize a region spec to a comparable content key.

    Review fix #4: control-region vs matched-region equality must compare
    resolved region *content*, not object identity/type. A path is resolved to
    its absolute filesystem location; a list of region strings is normalized to a
    sorted tuple; anything else falls back to its serialized form. This lets a
    ``str``/``Path`` control_regions be compared meaningfully against a
    ``list``/``str``/``Path`` matched_regions (previously ``list == Path`` was
    always False, so the "must be separate" guard never fired).
    """
    if region_spec is None:
        return None
    if isinstance(region_spec, (str, Path)):
        try:
            return ("path", str(Path(region_spec).expanduser().resolve()))
        except OSError:
            return ("path", str(region_spec))
    if isinstance(region_spec, (list, tuple)):
        return (
            "list",
            tuple(sorted(_region_spec_content_key(item) for item in region_spec)),
        )
    return ("other", str(_serialize_region_spec(region_spec)))


def _coerce_artifacts(sample: SampleSpec) -> list[DatasetArtifact]:
    if not sample.metadata:
        return []
    artifacts = sample.metadata.get("artifacts", [])
    coerced: list[DatasetArtifact] = []
    for artifact in artifacts:
        if isinstance(artifact, DatasetArtifact):
            coerced.append(artifact)
        elif isinstance(artifact, dict):
            coerced.append(DatasetArtifact(**artifact))
    return coerced


def _merge_sample_metadata(
    row: dict[str, Any],
    sample_metadata: dict[str, Any] | None,
) -> dict[str, Any]:
    if not sample_metadata:
        return row
    for key, value in sample_metadata.items():
        if key not in row:
            row[key] = value
    return row


def _requested_extract_artifact(
    sample: SampleSpec,
    *,
    motifs: list[str],
    matched_regions: Any,
    signal_normalization: str,
    feature_scaling: str,
    cluster_basis: str,
    window_size: int | None = None,
) -> DatasetArtifact:
    extract_path = Path(sample.extract_h5)
    return DatasetArtifact(
        sample_id=sample.sample_id,
        artifact_type="extract",
        path=extract_path,
        format=extract_path.suffix.lstrip(".") or "h5",
        params={
            "motifs": motifs,
            "matched_regions": _serialize_region_spec(matched_regions),
            "window_size": window_size,
            "signal_normalization": signal_normalization,
            "feature_scaling": feature_scaling,
            "cluster_basis": cluster_basis,
        },
        provenance={
            "pipeline": "parse_bam",
            "source_files": [str(extract_path)],
            "source_fingerprints": [_source_fingerprint(extract_path)],
            "upstream_lineage": [],
        },
        metadata={
            "schema_version": "artifact-v1",
            "package_version": _PACKAGE_VERSION,
        },
    )


def _requested_pileup_artifact(
    sample: SampleSpec,
    *,
    motifs: list[str],
    matched_regions: Any,
    signal_normalization: str,
    feature_scaling: str,
    cluster_basis: str,
) -> DatasetArtifact:
    if not sample.metadata or "pileup_path" not in sample.metadata:
        raise ValueError(
            f"Sample {sample.sample_id!r} is missing metadata['pileup_path'] for region_anchored mode."
        )
    pileup_path = Path(sample.metadata["pileup_path"])
    return DatasetArtifact(
        sample_id=sample.sample_id,
        artifact_type="pileup",
        path=pileup_path,
        format=pileup_path.suffix.lstrip(".") or "bed.gz",
        params={
            "motifs": motifs,
            "matched_regions": _serialize_region_spec(matched_regions),
            "signal_normalization": signal_normalization,
            "feature_scaling": feature_scaling,
            "cluster_basis": cluster_basis,
        },
        provenance={
            "pipeline": "parse_bam",
            "source_files": [str(pileup_path)],
            "source_fingerprints": [_source_fingerprint(pileup_path)],
            "upstream_lineage": [],
        },
        metadata={
            "schema_version": "artifact-v1",
            "package_version": _PACKAGE_VERSION,
        },
    )


def _require_pileup_path(sample: SampleSpec) -> str | Path:
    if sample.metadata and "pileup_path" in sample.metadata:
        return sample.metadata["pileup_path"]
    raise ValueError(
        f"Sample {sample.sample_id!r} is missing metadata['pileup_path'] for region_anchored mode."
    )


def _coverage_weighted_offset(
    data_matrix: np.ndarray,
    val_matrix: np.ndarray | None,
) -> float:
    """Single offset formula shared by read_global and region_anchored modes.

    Review fix #4: both normalization paths must use the same offset definition.
    When per-position valid counts are available and non-zero, the offset is the
    coverage-weighted mean modified fraction (sum of modified / sum of valid);
    otherwise it falls back to the unweighted mean of the data matrix. The
    region_anchored path has no valid-count matrix, so it consistently uses the
    unweighted-mean fallback via ``val_matrix=None``.
    """
    data = np.asarray(data_matrix, dtype=float)
    if val_matrix is not None:
        valid = np.asarray(val_matrix, dtype=float)
        if valid.sum() > 0:
            return float(data.sum() / valid.sum())
    return float(data.mean())


def _normalize_read_windows(
    result: cluster.ReadWindowExtractionResult,
    *,
    signal_normalization: str,
    control_result: cluster.ReadWindowExtractionResult | None = None,
) -> tuple[cluster.ReadWindowExtractionResult, dict[str, float | None]]:
    if signal_normalization == "none":
        return result, {"global_offset": None}
    if signal_normalization not in {"per_sample_global", "control_regions"}:
        raise ValueError(f"Unsupported signal_normalization: {signal_normalization}")

    data_matrix = np.asarray(result.data_matrix, dtype=float)
    val_matrix = (
        None
        if result.val_matrix is None
        else np.asarray(result.val_matrix, dtype=float)
    )
    offset_source = (
        result if signal_normalization == "per_sample_global" else control_result
    )
    if offset_source is None:
        raise ValueError(
            "control_regions normalization requires control-region read windows."
        )
    offset_matrix = np.asarray(offset_source.data_matrix, dtype=float)
    offset_val_matrix = (
        None
        if offset_source.val_matrix is None
        else np.asarray(offset_source.val_matrix, dtype=float)
    )
    global_offset = _coverage_weighted_offset(offset_matrix, offset_val_matrix)

    normalized = cluster.ReadWindowExtractionResult(
        data_matrix=data_matrix - global_offset,
        val_matrix=val_matrix,
        metadata=list(result.metadata),
        datasets=list(result.datasets),
        regions_dict=result.regions_dict,
    )
    return normalized, {"global_offset": global_offset}


def _select_feature_columns(
    feature_matrix: np.ndarray,
    feature_names: list[str],
    *,
    cluster_basis: str,
) -> tuple[np.ndarray, list[str]]:
    if cluster_basis not in _SUPPORTED_CLUSTER_BASIS:
        raise ValueError(f"Unsupported cluster_basis: {cluster_basis}")
    if cluster_basis == "shape_plus_level":
        return feature_matrix, feature_names

    if cluster_basis == "shape_only":
        keep = [name not in _LEVEL_FEATURES for name in feature_names]
    else:
        keep = [name in _LEVEL_FEATURES for name in feature_names]

    selected_names = [
        name for name, keep_name in zip(feature_names, keep, strict=False) if keep_name
    ]
    if not selected_names:
        raise ValueError(f"No features available for cluster_basis={cluster_basis!r}.")
    selected_matrix = feature_matrix[:, keep]
    return selected_matrix, selected_names


def _scale_features(
    training_matrix: np.ndarray,
    full_matrix: np.ndarray,
    *,
    feature_scaling: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    if feature_scaling == "none":
        return training_matrix, full_matrix, {"scaler": None}
    if feature_scaling != "robust_zscore":
        raise ValueError(f"Unsupported feature_scaling: {feature_scaling}")

    from sklearn.preprocessing import RobustScaler

    scaler = RobustScaler()
    training_scaled = scaler.fit_transform(training_matrix)
    full_scaled = scaler.transform(full_matrix)
    return training_scaled, full_scaled, {"scaler": "RobustScaler"}


def _cluster_label_strings(labels: np.ndarray) -> list[str]:
    values = []
    for label in np.asarray(labels):
        label_int = int(label)
        values.append("noise" if label_int < 0 else f"C{label_int}")
    return values


def _cluster_profiles(
    feature_matrix: np.ndarray,
    feature_names: list[str],
    cluster_labels: np.ndarray,
) -> pd.DataFrame:
    if feature_matrix.shape[0] == 0:
        return pd.DataFrame(columns=["cluster", "count", *feature_names])

    labels = np.asarray(cluster_labels, dtype=str)
    unique_labels = np.unique(labels)
    codes = pd.Categorical(labels, categories=unique_labels, ordered=True).codes
    counts = np.bincount(codes, minlength=len(unique_labels)).astype(
        np.int64, copy=False
    )

    sums = np.zeros((len(unique_labels), feature_matrix.shape[1]), dtype=np.float64)
    np.add.at(sums, codes, np.asarray(feature_matrix, dtype=np.float64))
    means = np.divide(
        sums,
        counts[:, None],
        out=np.zeros_like(sums),
        where=counts[:, None] > 0,
    )

    profiles = pd.DataFrame(means, columns=feature_names)
    profiles.insert(0, "count", counts)
    profiles.insert(0, "cluster", unique_labels)
    return profiles


def _sample_training_rows(
    feature_matrix: np.ndarray,
    *,
    training_sample_per_dataset: int,
    random_state: int,
) -> tuple[np.ndarray, np.ndarray]:
    if training_sample_per_dataset <= 0:
        raise ValueError("training_sample_per_dataset must be positive.")
    if feature_matrix.shape[0] <= training_sample_per_dataset:
        indices = np.arange(feature_matrix.shape[0], dtype=int)
        return feature_matrix, indices
    sampled_matrix, _, indices = cluster.sample_rows(
        feature_matrix,
        n=training_sample_per_dataset,
        random_state=random_state,
    )
    return sampled_matrix, indices


def _predict_in_chunks(model: Any, feature_matrix: np.ndarray) -> np.ndarray:
    labels: list[np.ndarray] = []
    for start in range(0, feature_matrix.shape[0], _PREDICTION_CHUNK_SIZE):
        stop = start + _PREDICTION_CHUNK_SIZE
        labels.append(np.asarray(model.predict(feature_matrix[start:stop])))
    return np.concatenate(labels, axis=0)


def _build_region_summary(assignments: pd.DataFrame) -> pd.DataFrame:
    summary = (
        assignments.groupby(
            ["region_id", "sample_id", "condition", "cluster"],
            sort=True,
        )
        .size()
        .reset_index(name="count")
    )
    totals = summary.groupby(["region_id", "sample_id", "condition"])[
        "count"
    ].transform("sum")
    summary["fraction"] = summary["count"] / totals
    return summary


def _assignments_have_region_coordinates(assignments: pd.DataFrame) -> bool:
    required = ("chromosome", "region_start", "region_end")
    if not set(required).issubset(assignments.columns):
        return False
    return assignments.loc[:, list(required)].notna().all().all()


def _region_id_from_coordinates(row: pd.Series) -> str:
    chrom = row.get("chromosome", row.get("chrom"))
    start = row.get("region_start", row.get("start"))
    end = row.get("region_end", row.get("end"))
    region_id = f"{chrom}:{int(start)}-{int(end)}"
    strand = row.get("region_strand", row.get("strand"))
    if pd.notna(strand):
        strand_value = str(strand)
        if strand_value in {"+", "-", "."}:
            return f"{region_id},{strand_value}"
    return region_id


def _build_read_global_region_summary(assignments: pd.DataFrame) -> pd.DataFrame | None:
    if not _assignments_have_region_coordinates(assignments):
        return None

    summarizer = getattr(cluster, "summarize_read_cluster_region_associations", None)
    if callable(summarizer):
        association_frames: list[pd.DataFrame] = []
        include_strand = "region_strand" in assignments.columns
        grouping_columns = ["sample_id", "condition"]
        grouped = assignments.groupby(grouping_columns, sort=False, dropna=False)
        for (sample_id, condition), sample_assignments in grouped:
            try:
                sample_summary = summarizer(
                    metadata=sample_assignments.to_dict("records"),
                    labels=sample_assignments["cluster"].to_numpy(),
                    include_strand=include_strand,
                )
            except (ValueError, KeyError, TypeError) as error:
                # Review fix #5: only recoverable, data-shape errors justify the
                # silent downgrade to the coordinate fallback below. Log the
                # reason so the downgrade is not invisible, and let unexpected
                # errors (bugs, MemoryError, KeyboardInterrupt, etc.) propagate
                # instead of being swallowed.
                _LOGGER.warning(
                    "summarize_read_cluster_region_associations failed for "
                    "sample_id=%r condition=%r (%s: %s); falling back to "
                    "coordinate-based region summary.",
                    sample_id,
                    condition,
                    type(error).__name__,
                    error,
                )
                association_frames = []
                break
            if not isinstance(sample_summary, pd.DataFrame) or sample_summary.empty:
                continue
            normalized_summary = sample_summary.copy()
            if "region_id" not in normalized_summary.columns:
                coordinate_columns = {"chrom", "start", "end"}
                if coordinate_columns.issubset(normalized_summary.columns):
                    normalized_summary["region_id"] = normalized_summary.apply(
                        _region_id_from_coordinates,
                        axis=1,
                    )
            normalized_summary["sample_id"] = sample_id
            normalized_summary["condition"] = condition
            association_frames.append(normalized_summary)

        if association_frames:
            region_summaries = pd.concat(association_frames, ignore_index=True)
            required_columns = {
                "region_id",
                "sample_id",
                "condition",
                "cluster",
                "count",
                "fraction",
            }
            if required_columns.issubset(region_summaries.columns):
                return region_summaries

    summary_source = assignments.copy()
    summary_source["region_id"] = summary_source.apply(
        _region_id_from_coordinates, axis=1
    )
    return _build_region_summary(summary_source)


def _select_discovery_hits(
    hits: pd.DataFrame,
    *,
    selection_mode: str,
    top_n: int | None,
) -> pd.DataFrame:
    resolved_top_n = _resolve_discovery_selection_top_n(
        selection_mode=selection_mode,
        top_n=top_n,
    )
    if selection_mode == "all":
        return hits.copy()
    return hits.head(resolved_top_n).copy()


def _resolve_discovery_selection_top_n(
    *,
    selection_mode: str,
    top_n: int | None,
) -> int | None:
    if selection_mode == "all":
        return None
    if selection_mode != "top_n":
        raise ValueError(
            f"Unsupported selection mode: {selection_mode!r}. Supported modes are 'top_n' and 'all'."
        )

    resolved_top_n = _DISCOVERY_SELECTION_DEFAULT_TOP_N if top_n is None else int(top_n)
    if resolved_top_n < 0:
        raise ValueError("selection.top_n must be non-negative.")
    return resolved_top_n


def _validate_shared_cluster_distribution_config(
    *,
    mode: str,
    clusterer: str,
    signal_normalization: str,
    feature_scaling: str,
    cluster_basis: str,
    sample_count: int,
) -> None:
    if mode not in {"read_global", "region_anchored"}:
        raise NotImplementedError(
            "The first workflow slice implements mode='read_global' and "
            "mode='region_anchored' only."
        )
    if sample_count < 2:
        raise ValueError("shared_cluster_distribution requires at least two datasets.")
    if clusterer != "minibatch_kmeans":
        raise NotImplementedError(
            "The first workflow slice supports clusterer='minibatch_kmeans' only."
        )
    if signal_normalization not in _SUPPORTED_SIGNAL_NORMALIZATION:
        raise ValueError(f"Unsupported signal_normalization: {signal_normalization}")
    if feature_scaling not in _SUPPORTED_FEATURE_SCALING:
        raise ValueError(f"Unsupported feature_scaling: {feature_scaling}")
    if cluster_basis not in _SUPPORTED_CLUSTER_BASIS:
        raise ValueError(f"Unsupported cluster_basis: {cluster_basis}")


def _selected_regions_to_region_spec(selected_regions: pd.DataFrame) -> list[str]:
    if selected_regions.empty:
        return []

    chrom_column = "chrom" if "chrom" in selected_regions.columns else "chromosome"
    region_spec: list[str] = []
    for row in selected_regions.itertuples(index=False):
        chrom = getattr(row, chrom_column)
        start = int(row.start)
        end = int(row.end)
        strand = getattr(row, "strand", ".")
        strand_value = strand if strand in {"+", "-", "."} else "."
        region_spec.append(f"{chrom}:{start}-{end},{strand_value}")
    return region_spec


def _normalize_region_id_value(
    row: pd.Series,
    *,
    default_region_ids: dict[str, str],
) -> Any:
    region_id = row.get("region_id")
    if pd.isna(region_id):
        return region_id

    region_id_str = str(region_id)
    if region_id_str in default_region_ids:
        return default_region_ids[region_id_str]

    if "," not in region_id_str and region_id_str.count(":") >= 2:
        region_core, strand = region_id_str.rsplit(":", 1)
        strand_value = strand if strand in {"+", "-", "."} else "."
        if region_core in default_region_ids:
            return default_region_ids[region_core]
        return f"{region_core},{strand_value}"

    chrom = row.get("chromosome", row.get("chrom"))
    start = row.get("start")
    end = row.get("end")
    if pd.notna(chrom) and pd.notna(start) and pd.notna(end):
        strand = row.get("strand", ".")
        strand_value = strand if strand in {"+", "-", "."} else "."
        return f"{chrom}:{int(start)}-{int(end)},{strand_value}"

    if "," in region_id_str:
        region_core, strand = region_id_str.rsplit(",", 1)
        strand_value = strand if strand in {"+", "-", "."} else "."
        return f"{region_core},{strand_value}"

    strand = row.get("strand")
    if pd.notna(strand):
        strand_value = str(strand)
        if strand_value in {"+", "-", "."}:
            return f"{region_id_str},{strand_value}"

    return region_id_str


def _normalize_region_id_frame(
    region_frame: pd.DataFrame | None,
    *,
    default_region_ids: dict[str, str],
) -> pd.DataFrame | None:
    if region_frame is None or "region_id" not in region_frame.columns:
        return region_frame

    normalized_frame = region_frame.copy()
    normalized_frame["region_id"] = normalized_frame.apply(
        _normalize_region_id_value,
        axis=1,
        default_region_ids=default_region_ids,
    )
    return normalized_frame


def _normalize_cluster_region_ids(
    clustering_result: SharedClusterResult,
    *,
    default_region_spec: list[str],
) -> SharedClusterResult:
    default_region_ids = {
        region_id.rsplit(",", 1)[0]: region_id for region_id in default_region_spec
    }

    normalized_assignments = _normalize_region_id_frame(
        clustering_result.assignments,
        default_region_ids=default_region_ids,
    )
    normalized_region_summaries = _normalize_region_id_frame(
        clustering_result.region_summaries,
        default_region_ids=default_region_ids,
    )

    return SharedClusterResult(
        model=clustering_result.model,
        assignments=normalized_assignments,
        cluster_distribution=clustering_result.cluster_distribution,
        condition_distribution=clustering_result.condition_distribution,
        distribution_change=clustering_result.distribution_change,
        cluster_profiles=clustering_result.cluster_profiles,
        region_summaries=normalized_region_summaries,
        plot_data=clustering_result.plot_data,
        figures=clustering_result.figures,
        metadata=clustering_result.metadata,
    )


def _require_region_summary_table_for_chip_atlas(
    cluster_result: SharedClusterResult,
) -> pd.DataFrame:
    if cluster_result.region_summaries is not None:
        region_summaries = cluster_result.region_summaries.copy()
    else:
        assignments = cluster_result.assignments.copy()
        if "region_id" not in assignments.columns:
            if _assignments_have_region_coordinates(assignments):
                assignments["region_id"] = assignments.apply(
                    _region_id_from_coordinates, axis=1
                )
            else:
                raise ValueError(
                    "SharedClusterResult does not include region_summaries, and assignments "
                    "do not contain region coordinates to derive region ids."
                )
        region_summaries = _build_region_summary(assignments)

    required_columns = {"cluster", "region_id"}
    if not required_columns.issubset(region_summaries.columns):
        missing = required_columns - set(region_summaries.columns)
        raise ValueError(
            "Cluster region summary table is missing required columns for ChIP-Atlas "
            f"workflow: {sorted(missing)}"
        )
    return region_summaries


def _rank_cluster_regions_for_chip_atlas(
    *,
    region_summaries: pd.DataFrame,
    cluster_label: str,
    min_fraction: float | None,
    top_n_regions: int | None,
) -> list[str]:
    cluster_rows = region_summaries.loc[
        region_summaries["cluster"] == cluster_label
    ].copy()
    if cluster_rows.empty:
        return []

    if "fraction" in cluster_rows.columns:
        ranked = (
            cluster_rows.groupby("region_id", dropna=False)["fraction"]
            .mean()
            .sort_values(ascending=False)
            .reset_index(name="score")
        )
        if min_fraction is not None:
            ranked = ranked.loc[ranked["score"] >= float(min_fraction)]
    elif "count" in cluster_rows.columns:
        ranked = (
            cluster_rows.groupby("region_id", dropna=False)["count"]
            .sum()
            .sort_values(ascending=False)
            .reset_index(name="score")
        )
    else:
        ranked = (
            cluster_rows.groupby("region_id", dropna=False)
            .size()
            .sort_values(ascending=False)
            .reset_index(name="score")
        )

    if top_n_regions is not None:
        ranked = ranked.head(int(top_n_regions))
    return [
        str(region_id)
        for region_id in ranked["region_id"].tolist()
        if pd.notna(region_id)
    ]


def resolve_regulatory_enrichment_spec(
    *,
    providers: Iterable[str] = ("screen", "unibind"),
    species: str | None = None,
    reference_genome: str | None = None,
    target_genome: str | None = None,
    crossmap_chain_file: str | Path | None = None,
    crossmap_chain_url: str | None = None,
    crossmap_chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
    strict_provider_support: bool = False,
) -> regulatory_enrichment.RegulatoryEnrichmentSpec:
    return regulatory_enrichment.RegulatoryEnrichmentSpec(
        providers=tuple(providers),
        species=species,
        reference_genome=reference_genome,
        target_genome=target_genome,
        crossmap_chain_file=crossmap_chain_file,
        crossmap_chain_url=crossmap_chain_url,
        crossmap_chain_cache_dir=crossmap_chain_cache_dir,
        crossmap_executable=crossmap_executable,
        strict_provider_support=strict_provider_support,
    )


def resolve_unibind_track_paths(
    *,
    track_paths: Iterable[str | Path] | None = None,
    trackhub_url: str | None = None,
    collection: str = "robust",
    assembly: str | None = None,
    search_terms: Iterable[str] | None = None,
    cache_dir: str | Path = "cache/unibind_tracks",
    max_tracks: int | None = None,
    allow_cached: bool = True,
    timeout_seconds: float = 60.0,
    convert_bigbed_to_bed: bool = False,
    regulatory_spec: regulatory_enrichment.RegulatoryEnrichmentSpec | None = None,
) -> list[Path]:
    resolved_assembly = assembly
    if resolved_assembly is None and regulatory_spec is not None:
        resolved_assembly = regulatory_spec.target_genome
    return regulatory_enrichment.resolve_unibind_track_paths(
        track_paths=track_paths,
        trackhub_url=trackhub_url,
        collection=collection,
        assembly=resolved_assembly,
        search_terms=search_terms,
        cache_dir=cache_dir,
        max_tracks=max_tracks,
        allow_cached=allow_cached,
        timeout_seconds=timeout_seconds,
        convert_bigbed_to_bed=convert_bigbed_to_bed,
    )


def unibind_tfbs_extraction_workflow(
    *,
    regions: Any,
    species: str | None = None,
    collection: str = "robust",
    tf_list: str | Path | Iterable[str] | None = None,
    experiment_ids: str | Path | Iterable[str] | None = None,
    name: str | None = None,
    email: str | None = None,
    endpoint_url: str = regulatory_enrichment.DEFAULT_UNIBIND_TFBS_EXTRACTION_URL,
    submit_timeout_seconds: float = 120.0,
    wait: bool = True,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 1200.0,
    download_outputs: bool = False,
    output_dir: str | Path = "cache/unibind_jobs",
    allow_cached_downloads: bool = True,
    download_timeout_seconds: float = 120.0,
    regulatory_spec: regulatory_enrichment.RegulatoryEnrichmentSpec | None = None,
) -> UniBindJobResult:
    resolved_species = species
    if resolved_species is None and regulatory_spec is not None:
        resolved_species = regulatory_spec.species
    return regulatory_enrichment.run_unibind_tfbs_extraction(
        regions=regions,
        species=resolved_species or "homo_sapiens",
        collection=collection,
        tf_list=tf_list,
        experiment_ids=experiment_ids,
        name=name,
        email=email,
        endpoint_url=endpoint_url,
        submit_timeout_seconds=submit_timeout_seconds,
        wait=wait,
        poll_interval_seconds=poll_interval_seconds,
        timeout_seconds=timeout_seconds,
        download_outputs=download_outputs,
        output_dir=output_dir,
        allow_cached_downloads=allow_cached_downloads,
        download_timeout_seconds=download_timeout_seconds,
    )


def unibind_enrichment_workflow(
    *,
    regions: Any,
    analysis_type: str = "oneSetBg",
    background_regions: Any | None = None,
    comparison_regions: Any | None = None,
    species: str | None = None,
    collection: str = "robust",
    name: str | None = None,
    email: str | None = None,
    endpoint_url: str = regulatory_enrichment.DEFAULT_UNIBIND_ENRICHMENT_URL,
    submit_timeout_seconds: float = 120.0,
    wait: bool = True,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 1200.0,
    download_outputs: bool = False,
    output_dir: str | Path = "cache/unibind_jobs",
    allow_cached_downloads: bool = True,
    download_timeout_seconds: float = 120.0,
    regulatory_spec: regulatory_enrichment.RegulatoryEnrichmentSpec | None = None,
) -> UniBindJobResult:
    resolved_species = species
    if resolved_species is None and regulatory_spec is not None:
        resolved_species = regulatory_spec.species
    return regulatory_enrichment.run_unibind_enrichment(
        regions=regions,
        analysis_type=analysis_type,
        background_regions=background_regions,
        comparison_regions=comparison_regions,
        species=resolved_species or "homo_sapiens",
        collection=collection,
        name=name,
        email=email,
        endpoint_url=endpoint_url,
        submit_timeout_seconds=submit_timeout_seconds,
        wait=wait,
        poll_interval_seconds=poll_interval_seconds,
        timeout_seconds=timeout_seconds,
        download_outputs=download_outputs,
        output_dir=output_dir,
        allow_cached_downloads=allow_cached_downloads,
        download_timeout_seconds=download_timeout_seconds,
    )


def chip_atlas_search_peak_datasets_workflow(
    *,
    antigen: str,
    genome: str = "hg38",
    cell_type: str | None = None,
    antigen_class: str | None = None,
    cell_type_class: str | None = None,
    threshold: str = "05",
    match_mode: str = "contains",
    max_results: int | None = None,
    experiment_list_url: str = chip_atlas.DEFAULT_EXPERIMENT_LIST_URL,
    cache_dir: str | Path = chip_atlas.DEFAULT_METADATA_CACHE_DIR,
    allow_cached_metadata: bool = True,
    timeout_seconds: float = 120.0,
) -> pd.DataFrame:
    return chip_atlas.search_peak_datasets(
        antigen=antigen,
        genome=genome,
        cell_type=cell_type,
        antigen_class=antigen_class,
        cell_type_class=cell_type_class,
        threshold=threshold,
        match_mode=match_mode,
        max_results=max_results,
        experiment_list_url=experiment_list_url,
        cache_dir=cache_dir,
        allow_cached_metadata=allow_cached_metadata,
        timeout_seconds=timeout_seconds,
    )


def chip_atlas_download_peak_datasets_workflow(
    *,
    datasets: pd.DataFrame,
    dataset_ids: Iterable[str] | None = None,
    output_dir: str | Path = "cache/chip_atlas/peak_sets",
    include_complete_sorted: bool = True,
    include_top_n: int | None = 3000,
    include_bottom_n: int | None = 3000,
    stratify: str | int | None = None,
    allow_cached: bool = True,
    timeout_seconds: float = 180.0,
    crossmap_target_genome: str | None = None,
    crossmap_chain_file: str | Path | None = None,
    crossmap_chain_url: str | None = None,
    crossmap_chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
) -> pd.DataFrame:
    return chip_atlas.download_peak_datasets(
        datasets=datasets,
        dataset_ids=dataset_ids,
        output_dir=output_dir,
        include_complete_sorted=include_complete_sorted,
        include_top_n=include_top_n,
        include_bottom_n=include_bottom_n,
        stratify=stratify,
        allow_cached=allow_cached,
        timeout_seconds=timeout_seconds,
        crossmap_target_genome=crossmap_target_genome,
        crossmap_chain_file=crossmap_chain_file,
        crossmap_chain_url=crossmap_chain_url,
        crossmap_chain_cache_dir=crossmap_chain_cache_dir,
        crossmap_executable=crossmap_executable,
    )


def chip_atlas_enrichment_workflow(
    *,
    regions: pd.DataFrame | str | Path | list[str],
    genome: str = "hg38",
    regions_genome: str | None = None,
    antigen_class: str | None = "TFs and others",
    antigen: str | None = None,
    cell_type_class: str | None = "No description",
    cell_type: str | None = None,
    distance: int | None = None,
    threshold: str | None = "100",
    crossmap_chain_file: str | Path | None = None,
    crossmap_chain_url: str | None = None,
    crossmap_chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
    params: dict[str, Any] | None = None,
    wait: bool = True,
    fetch_results: bool = True,
    raise_on_failure: bool = True,
    submit_url: str = chip_atlas.DEFAULT_SUBMIT_URL,
    status_url: str = chip_atlas.DEFAULT_STATUS_URL,
    result_url: str = chip_atlas.DEFAULT_RESULT_URL,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 600.0,
) -> ChipAtlasEnrichmentResult:
    return chip_atlas.run_enrichment(
        regions=regions,
        genome=genome,
        regions_genome=regions_genome,
        antigen_class=antigen_class,
        antigen=antigen,
        cell_type_class=cell_type_class,
        cell_type=cell_type,
        distance=distance,
        threshold=threshold,
        crossmap_chain_file=crossmap_chain_file,
        crossmap_chain_url=crossmap_chain_url,
        crossmap_chain_cache_dir=crossmap_chain_cache_dir,
        crossmap_executable=crossmap_executable,
        params=params,
        wait=wait,
        fetch_results=fetch_results,
        raise_on_failure=raise_on_failure,
        submit_url=submit_url,
        status_url=status_url,
        result_url=result_url,
        poll_interval_seconds=poll_interval_seconds,
        timeout_seconds=timeout_seconds,
    )


def chip_atlas_cluster_enrichment_workflow(
    *,
    cluster_result: SharedClusterResult,
    genome: str = "hg38",
    regions_genome: str | None = None,
    clusters: Iterable[str] | None = None,
    min_fraction: float | None = None,
    top_n_regions: int | None = None,
    mode: str = "per_cluster",
    antigen_class: str | None = "TFs and others",
    antigen: str | None = None,
    cell_type_class: str | None = "No description",
    cell_type: str | None = None,
    distance: int | None = None,
    threshold: str | None = "100",
    crossmap_chain_file: str | Path | None = None,
    crossmap_chain_url: str | None = None,
    crossmap_chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
    params: dict[str, Any] | None = None,
    wait: bool = True,
    fetch_results: bool = True,
    raise_on_failure: bool = True,
    submit_url: str = chip_atlas.DEFAULT_SUBMIT_URL,
    status_url: str = chip_atlas.DEFAULT_STATUS_URL,
    result_url: str = chip_atlas.DEFAULT_RESULT_URL,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 600.0,
) -> dict[str, ChipAtlasEnrichmentResult]:
    if mode not in {"per_cluster", "combined"}:
        raise ValueError("mode must be 'per_cluster' or 'combined'.")

    region_summaries = _require_region_summary_table_for_chip_atlas(cluster_result)
    available_clusters = [
        str(value) for value in pd.unique(region_summaries["cluster"])
    ]
    selected_clusters = (
        [str(cluster_label) for cluster_label in clusters]
        if clusters is not None
        else available_clusters
    )
    unknown = sorted(set(selected_clusters) - set(available_clusters))
    if unknown:
        raise ValueError(
            "Requested clusters were not present in cluster_result.region_summaries: "
            f"{unknown}"
        )

    per_cluster_region_ids: dict[str, list[str]] = {}
    for cluster_label in selected_clusters:
        region_ids = _rank_cluster_regions_for_chip_atlas(
            region_summaries=region_summaries,
            cluster_label=cluster_label,
            min_fraction=min_fraction,
            top_n_regions=top_n_regions,
        )
        if region_ids:
            per_cluster_region_ids[cluster_label] = region_ids

    if not per_cluster_region_ids:
        raise ValueError(
            "No regions met the requested filters for ChIP-Atlas enrichment."
        )

    queries: dict[str, list[str]]
    if mode == "combined":
        union_region_ids = sorted(
            {region_id for ids in per_cluster_region_ids.values() for region_id in ids}
        )
        queries = {"combined": union_region_ids}
    else:
        queries = per_cluster_region_ids

    results: dict[str, ChipAtlasEnrichmentResult] = {}
    for query_key, region_ids in queries.items():
        bed_frame = chip_atlas.region_ids_to_bed_dataframe(region_ids)
        enrichment = chip_atlas_enrichment_workflow(
            regions=bed_frame,
            genome=genome,
            regions_genome=regions_genome,
            antigen_class=antigen_class,
            antigen=antigen,
            cell_type_class=cell_type_class,
            cell_type=cell_type,
            distance=distance,
            threshold=threshold,
            crossmap_chain_file=crossmap_chain_file,
            crossmap_chain_url=crossmap_chain_url,
            crossmap_chain_cache_dir=crossmap_chain_cache_dir,
            crossmap_executable=crossmap_executable,
            params=params,
            wait=wait,
            fetch_results=fetch_results,
            raise_on_failure=raise_on_failure,
            submit_url=submit_url,
            status_url=status_url,
            result_url=result_url,
            poll_interval_seconds=poll_interval_seconds,
            timeout_seconds=timeout_seconds,
        )
        enrichment.metadata["workflow"] = {
            "query_key": query_key,
            "cluster_mode": mode,
            "source_clusters": selected_clusters,
            "region_count": len(region_ids),
            "min_fraction": min_fraction,
            "top_n_regions": top_n_regions,
        }
        results[query_key] = enrichment
    return results


def modkit_dmr_pair_workflow(
    *,
    control_bed_methyl: str | Path,
    experiment_bed_methyl: str | Path,
    ref_genome: str | Path,
    out_path: str | Path,
    regions_bed: str | Path | None = None,
    segment_path: str | Path | None = None,
    bases: Iterable[str] = ("A",),
    assign_codes: Iterable[str] | None = None,
    min_valid_coverage: int = 0,
    dmr_prior: float | None = None,
    diff_stay: float | None = None,
    significance_factor: float | None = None,
    decay_distance: int | None = None,
    max_gap_size: int | None = None,
    log_transition_decay: bool = False,
    fine_grained: bool = False,
    prior_alpha: float | None = None,
    prior_beta: float | None = None,
    delta: float | None = None,
    n_sample_records: int | None = None,
    max_coverages: tuple[int, int] | None = None,
    cap_coverages: bool = False,
    missing: str | None = None,
    threads: int | None = None,
    io_threads: int | None = None,
    batch_size: int | None = None,
    interval_size: int | None = None,
    header: bool = True,
    force: bool = True,
    suppress_progress: bool = True,
    log_filepath: str | Path | None = None,
    modkit_executable: str | Path | None = None,
    pvalue_max: float = 0.01,
    abs_effect_size_min: float = 0.1,
    min_total_coverage: int | None = None,
) -> ModkitDMRPairResult:
    return dmr.run_dmr_pair(
        control_bed_methyl=control_bed_methyl,
        experiment_bed_methyl=experiment_bed_methyl,
        ref_genome=ref_genome,
        out_path=out_path,
        regions_bed=regions_bed,
        segment_path=segment_path,
        bases=list(bases),
        assign_codes=None if assign_codes is None else list(assign_codes),
        min_valid_coverage=min_valid_coverage,
        dmr_prior=dmr_prior,
        diff_stay=diff_stay,
        significance_factor=significance_factor,
        decay_distance=decay_distance,
        max_gap_size=max_gap_size,
        log_transition_decay=log_transition_decay,
        fine_grained=fine_grained,
        prior_alpha=prior_alpha,
        prior_beta=prior_beta,
        delta=delta,
        n_sample_records=n_sample_records,
        max_coverages=max_coverages,
        cap_coverages=cap_coverages,
        missing=missing,
        threads=threads,
        io_threads=io_threads,
        batch_size=batch_size,
        interval_size=interval_size,
        header=header,
        force=force,
        suppress_progress=suppress_progress,
        log_filepath=log_filepath,
        modkit_executable=modkit_executable,
        pvalue_max=pvalue_max,
        abs_effect_size_min=abs_effect_size_min,
        min_total_coverage=min_total_coverage,
    )


def modkit_dmr_multi_workflow(
    *,
    samples: dict[str, str | Path] | Iterable[tuple[str, str | Path]],
    regions_bed: str | Path,
    ref_genome: str | Path,
    out_dir: str | Path,
    bases: Iterable[str] = ("A",),
    assign_codes: Iterable[str] | None = None,
    min_valid_coverage: int = 0,
    missing: str | None = None,
    threads: int | None = None,
    io_threads: int | None = None,
    prefix: str | None = None,
    header: bool = True,
    force: bool = True,
    suppress_progress: bool = True,
    log_filepath: str | Path | None = None,
    modkit_executable: str | Path | None = None,
) -> ModkitDMRMultiResult:
    return dmr.run_dmr_multi(
        samples=samples,
        regions_bed=regions_bed,
        ref_genome=ref_genome,
        out_dir=out_dir,
        bases=list(bases),
        assign_codes=None if assign_codes is None else list(assign_codes),
        min_valid_coverage=min_valid_coverage,
        missing=missing,
        threads=threads,
        io_threads=io_threads,
        prefix=prefix,
        header=header,
        force=force,
        suppress_progress=suppress_progress,
        log_filepath=log_filepath,
        modkit_executable=modkit_executable,
    )


def modkit_dmr_multi_from_samples_workflow(
    *,
    samples: Iterable[SampleSpec],
    regions_bed: str | Path,
    ref_genome: str | Path,
    out_dir: str | Path,
    bases: Iterable[str] = ("A",),
    assign_codes: Iterable[str] | None = None,
    min_valid_coverage: int = 0,
    missing: str | None = None,
    threads: int | None = None,
    io_threads: int | None = None,
    prefix: str | None = None,
    header: bool = True,
    force: bool = True,
    suppress_progress: bool = True,
    log_filepath: str | Path | None = None,
    modkit_executable: str | Path | None = None,
) -> ModkitDMRMultiResult:
    sample_list = list(samples)
    if len(sample_list) < 2:
        raise ValueError(
            "modkit_dmr_multi_from_samples_workflow requires at least two samples."
        )

    dmr_samples: dict[str, str | Path] = {}
    for sample in sample_list:
        metadata = sample.metadata or {}
        pileup_path = metadata.get("pileup_path")
        if pileup_path is None:
            raise ValueError(
                f"Sample {sample.sample_id!r} is missing metadata['pileup_path'] "
                "required for modkit DMR multi workflow."
            )
        dmr_samples[sample.sample_id] = pileup_path

    return modkit_dmr_multi_workflow(
        samples=dmr_samples,
        regions_bed=regions_bed,
        ref_genome=ref_genome,
        out_dir=out_dir,
        bases=list(bases),
        assign_codes=None if assign_codes is None else list(assign_codes),
        min_valid_coverage=min_valid_coverage,
        missing=missing,
        threads=threads,
        io_threads=io_threads,
        prefix=prefix,
        header=header,
        force=force,
        suppress_progress=suppress_progress,
        log_filepath=log_filepath,
        modkit_executable=modkit_executable,
    )


def discovery_cluster_workflow(
    *,
    samples: Iterable[SampleSpec],
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    discovery: dict[str, Any],
    clustering: dict[str, Any],
    selection: dict[str, Any] | None = None,
    quiet: bool = True,
) -> RegionDiscoveryClusterResult:
    sample_list = list(samples)
    motif_list = list(motifs)
    selection_config = dict(selection or {})
    selection_mode = str(selection_config.get("mode", "top_n"))
    selection_top_n = selection_config.get("top_n")
    resolved_top_n = _resolve_discovery_selection_top_n(
        selection_mode=selection_mode,
        top_n=selection_top_n,
    )
    _validate_shared_cluster_distribution_config(
        mode=str(clustering.get("mode", "read_global")),
        clusterer=str(clustering.get("clusterer", "minibatch_kmeans")),
        signal_normalization=str(clustering.get("signal_normalization", "none")),
        feature_scaling=str(clustering.get("feature_scaling", "robust_zscore")),
        cluster_basis=str(clustering.get("cluster_basis", "shape_plus_level")),
        sample_count=len(sample_list),
    )

    discovery_result = region_discovery.scan_genome(
        samples=sample_list,
        motifs=motif_list,
        genome_sizes=genome_sizes,
        **discovery,
    )
    selected_hits = _select_discovery_hits(
        discovery_result.hits,
        selection_mode=selection_mode,
        top_n=selection_top_n,
    )
    if selected_hits.empty:
        raise ValueError("No discovery hits remained after selection.")

    selected_regions = region_discovery.hits_to_bed(selected_hits)
    matched_regions = _selected_regions_to_region_spec(selected_regions)
    clustering_result = shared_cluster_distribution(
        samples=sample_list,
        motifs=motif_list,
        matched_regions=matched_regions,
        quiet=quiet,
        **clustering,
    )
    return RegionDiscoveryClusterResult(
        discovery=discovery_result,
        clustering=clustering_result,
        selected_regions=selected_regions,
        metadata={
            "selection": {
                "mode": selection_mode,
                "top_n": resolved_top_n,
                "discovered_hit_count": int(len(discovery_result.hits)),
                "selected_hit_count": int(len(selected_hits)),
            }
        },
    )


def discovery_cluster_contrast_workflow(
    *,
    samples: Iterable[SampleSpec],
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    discovery: dict[str, Any],
    clustering: dict[str, Any],
    contrasts: dict[str, Any],
    selection: dict[str, Any] | None = None,
) -> RegionDiscoveryClusterContrastResult:
    sample_list = list(samples)
    motif_list = list(motifs)
    contrast_config = dict(contrasts)
    contrast_spec = contrast_config.get("contrast")
    if contrast_spec is None:
        raise ValueError(
            "discovery_cluster_contrast_workflow requires contrasts['contrast']."
        )

    region_contrasts.validate_region_contrast_request(
        analysis_unit=str(contrast_config.get("analysis_unit", "ensemble_region")),
        representation=str(contrast_config.get("representation", "modified_fraction")),
        signal_source=str(contrast_config.get("signal_source", "pileup_counts")),
        test=str(contrast_config.get("test", "effect_size_only")),
    )

    discovery_cluster_result = discovery_cluster_workflow(
        samples=sample_list,
        motifs=motif_list,
        genome_sizes=genome_sizes,
        discovery=discovery,
        clustering=clustering,
        selection=selection,
    )
    selected_region_spec = _selected_regions_to_region_spec(
        discovery_cluster_result.selected_regions
    )
    contrast_regions = contrast_config.pop("regions", None)
    contrast_scope = "selected"
    if contrast_regions is None:
        contrast_regions = selected_region_spec
    else:
        contrast_scope = "custom"

    contrast_result = region_contrasts.score_regions(
        samples=sample_list,
        regions=contrast_regions,
        motifs=motif_list,
        **contrast_config,
    )
    normalized_clustering = _normalize_cluster_region_ids(
        discovery_cluster_result.clustering,
        default_region_spec=selected_region_spec,
    )

    return RegionDiscoveryClusterContrastResult(
        discovery=discovery_cluster_result.discovery,
        clustering=normalized_clustering,
        contrasts=contrast_result,
        selected_regions=discovery_cluster_result.selected_regions,
        metadata={
            **discovery_cluster_result.metadata,
            "contrast_scope": contrast_scope,
            "full_scan_windows": discovery_cluster_result.discovery.windows.copy(),
        },
    )


def _build_shared_cluster_result(
    *,
    mode: str,
    motifs: list[str],
    feature_blocks: list[np.ndarray],
    training_blocks: list[np.ndarray],
    metadata_rows: list[dict[str, Any]],
    feature_names: list[str],
    signal_normalization: str,
    feature_scaling: str,
    cluster_basis: str,
    clusterer: str,
    n_clusters: int,
    training_sample_per_dataset: int,
    artifact_policy: str,
    random_state: int,
    make_plots: bool,
    matched_regions: Any,
    sample_training_rows: dict[str, int],
    sample_dataset_sizes: dict[str, int],
    sample_normalization: dict[str, dict[str, float | None]],
    cache_hits: dict[str, str],
    cache_misses: list[str],
    sample_rebuild_decisions: dict[str, str],
    quiet: bool,
    region_summaries: pd.DataFrame | None = None,
) -> SharedClusterResult:
    full_matrix = np.vstack(feature_blocks)
    training_matrix = np.vstack(training_blocks)
    if training_matrix.shape[0] < n_clusters:
        raise ValueError(
            "Pooled training subset has fewer rows than n_clusters. "
            "Reduce n_clusters or increase training_sample_per_dataset."
        )

    training_scaled, full_scaled, scaler_meta = _scale_features(
        training_matrix,
        full_matrix,
        feature_scaling=feature_scaling,
    )
    clustering_result = cluster.cluster_read_windows(
        training_scaled,
        method=clusterer,
        n_clusters=n_clusters,
        random_state=random_state,
    )
    label_mapping = cluster.cluster_label_mapping(
        clustering_result.labels_raw,
        clustering_result.labels_size_ordered,
    )
    if not hasattr(clustering_result.model, "predict"):
        raise TypeError(
            "Fitted clustering model does not support prediction for full assignment."
        )
    predicted_raw = _predict_in_chunks(clustering_result.model, full_scaled)
    predicted_ordered = cluster.apply_cluster_label_mapping(
        predicted_raw, label_mapping
    )

    assignments = pd.DataFrame(metadata_rows)
    assignments["cluster"] = _cluster_label_strings(predicted_ordered)

    cluster_distribution = distribution.build_cluster_distribution(assignments)
    condition_distribution = distribution.build_condition_distribution(
        cluster_distribution
    )
    cluster_profiles = _cluster_profiles(
        full_matrix,
        feature_names,
        assignments["cluster"].to_numpy(),
    )
    plot_data = {
        "cluster_distribution_bar": plotting.prepare_cluster_distribution_bar_data(
            cluster_distribution
        ),
        "cluster_distribution_heatmap": plotting.prepare_cluster_distribution_heatmap_data(
            condition_distribution
        ),
    }
    if region_summaries is None and mode == "region_anchored":
        region_summaries = _build_region_summary(assignments)
    if region_summaries is None and mode == "read_global":
        region_summaries = _build_read_global_region_summary(assignments)

    model = SharedClusterModel(
        mode=mode,
        motifs=motifs,
        feature_names=list(feature_names),
        preprocessing={
            "signal_normalization": signal_normalization,
            "feature_scaling": feature_scaling,
            "cluster_basis": cluster_basis,
            **scaler_meta,
        },
        estimator=clustering_result.model,
        cluster_labels=sorted(assignments["cluster"].unique()),
        fit_metadata={
            "clusterer": clusterer,
            "n_clusters": n_clusters,
            "artifact_policy": artifact_policy,
            "training_sample_per_dataset": training_sample_per_dataset,
            "training_rows": int(training_scaled.shape[0]),
            "sample_training_rows": sample_training_rows,
            "metrics": clustering_result.metrics,
        },
    )
    return SharedClusterResult(
        model=model,
        assignments=assignments,
        cluster_distribution=cluster_distribution,
        condition_distribution=condition_distribution,
        distribution_change=None,
        cluster_profiles=cluster_profiles,
        region_summaries=region_summaries,
        plot_data=plot_data,
        figures={},
        metadata={
            "mode": mode,
            "artifact_policy": artifact_policy,
            "cache_hits": cache_hits,
            "cache_misses": cache_misses,
            "sample_rebuild_decisions": sample_rebuild_decisions,
            "signal_normalization": signal_normalization,
            "feature_scaling": feature_scaling,
            "cluster_basis": cluster_basis,
            "make_plots": make_plots,
            "matched_regions": matched_regions,
            "sample_dataset_sizes": sample_dataset_sizes,
            "rows_after_filtering": sample_dataset_sizes,
            "sample_normalization": sample_normalization,
            "quiet": quiet,
        },
    )


def shared_cluster_distribution(
    *,
    samples: Iterable[SampleSpec],
    mode: str,
    motifs: Iterable[str],
    matched_regions: RegionSpec = None,
    signal_normalization: str = "none",
    feature_scaling: str = "robust_zscore",
    cluster_basis: str = "shape_plus_level",
    clusterer: str = "minibatch_kmeans",
    n_clusters: int = 8,
    training_sample_per_dataset: int = 100_000,
    artifact_policy: str = "prefer_cached",
    random_state: int = 42,
    make_plots: bool = True,
    cores: int | None = None,
    quiet: bool = True,
) -> SharedClusterResult:
    """
    Fit one shared clustering model on pooled read windows and assign all reads back
    onto those fixed cluster boundaries.
    """

    sample_list = list(samples)
    motif_list = list(motifs)
    if not sample_list:
        raise ValueError("samples must contain at least one sample.")
    if not motif_list:
        raise ValueError("motifs must contain at least one motif.")
    _validate_shared_cluster_distribution_config(
        mode=mode,
        clusterer=clusterer,
        signal_normalization=signal_normalization,
        feature_scaling=feature_scaling,
        cluster_basis=cluster_basis,
        sample_count=len(sample_list),
    )

    if mode == "region_anchored":
        if matched_regions is None:
            raise ValueError(
                "mode='region_anchored' requires matched_regions; per-sample regions_bed is not used "
                "as an implicit fallback."
            )
        pileup_paths: dict[str, str | Path] = {}
        cache_hits: dict[str, str] = {}
        cache_misses: list[str] = []
        sample_rebuild_decisions: dict[str, str] = {}

        for sample in tqdm(sample_list, desc="Processing samples", disable=quiet):
            requested_artifact = _requested_pileup_artifact(
                sample,
                motifs=motif_list,
                matched_regions=matched_regions,
                signal_normalization=signal_normalization,
                feature_scaling=feature_scaling,
                cluster_basis=cluster_basis,
            )
            resolved_artifact = resolve_artifact(
                requested_artifact,
                _coerce_artifacts(sample),
                artifact_policy=artifact_policy,
            )
            if resolved_artifact is not None:
                pileup_paths[sample.sample_id] = resolved_artifact.path
                cache_hits[sample.sample_id] = str(resolved_artifact.path)
                sample_rebuild_decisions[sample.sample_id] = "cache_hit"
            else:
                pileup_paths[sample.sample_id] = _require_pileup_path(sample)
                cache_misses.append(sample.sample_id)
                sample_rebuild_decisions[sample.sample_id] = "rebuilt_from_raw"

        feature_matrix, metadata_rows = region_analysis.build_region_feature_table(
            samples=sample_list,
            motifs=motif_list,
            matched_regions=matched_regions,
            pileup_paths=pileup_paths,
            cores=cores,
            quiet=quiet,
        )
        sample_by_id = {sample.sample_id: sample for sample in sample_list}
        for row in metadata_rows:
            sample = sample_by_id[row["sample_id"]]
            _merge_sample_metadata(row, sample.metadata)
        sample_row_positions: dict[str, list[int]] = {
            sample.sample_id: [] for sample in sample_list
        }
        for index, row in enumerate(metadata_rows):
            sample_id = row["sample_id"]
            if sample_id in sample_row_positions:
                sample_row_positions[sample_id].append(index)
        sample_row_indices = {
            sample_id: np.asarray(indices, dtype=int)
            for sample_id, indices in sample_row_positions.items()
        }
        sample_dataset_sizes = {
            sample.sample_id: int(sample_row_indices[sample.sample_id].size)
            for sample in sample_list
        }
        sample_normalization = {
            sample.sample_id: {"global_offset": None} for sample in sample_list
        }
        region_matrix = np.asarray(feature_matrix, dtype=float)
        if signal_normalization in {"per_sample_global", "control_regions"}:
            for sample in sample_list:
                sample_id = sample.sample_id
                sample_indices = sample_row_indices[sample_id]
                if sample_indices.size == 0:
                    continue
                if signal_normalization == "control_regions":
                    control_regions = sample.regions_bed
                    # Review fix #4: compare resolved region content, not object
                    # identity. Previously ``list == Path`` was always False so the
                    # guard never triggered and normalization could silently divide
                    # by the matched regions themselves.
                    if control_regions is None or _region_spec_content_key(
                        control_regions
                    ) == _region_spec_content_key(matched_regions):
                        raise ValueError(
                            "signal_normalization='control_regions' in region_anchored mode "
                            "requires sample.regions_bed to provide separate control regions."
                        )
                    control_artifact = _requested_pileup_artifact(
                        sample,
                        motifs=motif_list,
                        matched_regions=control_regions,
                        signal_normalization=signal_normalization,
                        feature_scaling=feature_scaling,
                        cluster_basis=cluster_basis,
                    )
                    resolved_control_artifact = resolve_artifact(
                        control_artifact,
                        _coerce_artifacts(sample),
                        artifact_policy=artifact_policy,
                    )
                    control_path = (
                        resolved_control_artifact.path
                        if resolved_control_artifact is not None
                        else _require_pileup_path(sample)
                    )
                    control_matrix, _ = region_analysis.build_region_feature_table(
                        samples=[sample],
                        motifs=motif_list,
                        matched_regions=control_regions,
                        pileup_paths={sample.sample_id: control_path},
                        cores=cores,
                        quiet=quiet,
                    )
                    # Review fix #4: use the shared offset formula. The region
                    # feature table carries no per-position valid-count matrix,
                    # so val_matrix=None yields the unweighted-mean branch,
                    # matching read_global's fallback behavior.
                    offset = _coverage_weighted_offset(
                        np.asarray(control_matrix, dtype=float), None
                    )
                else:
                    offset = _coverage_weighted_offset(
                        region_matrix[sample_indices], None
                    )
                region_matrix[sample_indices] = region_matrix[sample_indices] - offset
                sample_normalization[sample_id] = {"global_offset": offset}

        if cluster_basis == "shape_only":
            region_matrix = region_matrix - region_matrix.mean(axis=1, keepdims=True)
            feature_names = [f"shape_{idx}" for idx in range(region_matrix.shape[1])]
        elif cluster_basis == "level_only":
            region_matrix = region_matrix.mean(axis=1, keepdims=True)
            feature_names = ["region_mean_mod_fraction"]
        else:
            feature_names = [f"pos_{idx}" for idx in range(region_matrix.shape[1])]

        feature_blocks: list[np.ndarray] = []
        training_blocks: list[np.ndarray] = []
        sample_training_rows: dict[str, int] = {}
        for sample in sample_list:
            sample_matrix = region_matrix[sample_row_indices[sample.sample_id]]
            feature_blocks.append(sample_matrix)
            training_matrix, _ = _sample_training_rows(
                sample_matrix,
                training_sample_per_dataset=training_sample_per_dataset,
                random_state=random_state + len(training_blocks),
            )
            training_blocks.append(training_matrix)
            sample_training_rows[sample.sample_id] = int(training_matrix.shape[0])

        return _build_shared_cluster_result(
            mode=mode,
            motifs=motif_list,
            feature_blocks=feature_blocks,
            training_blocks=training_blocks,
            metadata_rows=metadata_rows,
            feature_names=feature_names,
            signal_normalization=signal_normalization,
            feature_scaling=feature_scaling,
            cluster_basis=cluster_basis,
            clusterer=clusterer,
            n_clusters=n_clusters,
            training_sample_per_dataset=training_sample_per_dataset,
            artifact_policy=artifact_policy,
            random_state=random_state,
            make_plots=make_plots,
            matched_regions=matched_regions,
            sample_training_rows=sample_training_rows,
            sample_dataset_sizes=sample_dataset_sizes,
            sample_normalization=sample_normalization,
            cache_hits=cache_hits,
            cache_misses=cache_misses,
            sample_rebuild_decisions=sample_rebuild_decisions,
            quiet=quiet,
        )

    feature_blocks: list[np.ndarray] = []
    metadata_rows: list[dict[str, Any]] = []
    training_blocks: list[np.ndarray] = []
    selected_feature_names: list[str] | None = None
    sample_training_rows: dict[str, int] = {}
    sample_dataset_sizes: dict[str, int] = {}
    sample_normalization: dict[str, dict[str, float | None]] = {}
    cache_hits: dict[str, str] = {}
    cache_misses: list[str] = []
    sample_rebuild_decisions: dict[str, str] = {}

    for sample_index, sample in enumerate(sample_list):
        requested_artifact = _requested_extract_artifact(
            sample,
            motifs=motif_list,
            matched_regions=matched_regions,
            signal_normalization=signal_normalization,
            feature_scaling=feature_scaling,
            cluster_basis=cluster_basis,
        )
        resolved_artifact = resolve_artifact(
            requested_artifact,
            _coerce_artifacts(sample),
            artifact_policy=artifact_policy,
        )
        extract_path = sample.extract_h5
        if resolved_artifact is not None:
            extract_path = resolved_artifact.path
            cache_hits[sample.sample_id] = str(resolved_artifact.path)
            sample_rebuild_decisions[sample.sample_id] = "cache_hit"
        else:
            cache_misses.append(sample.sample_id)
            sample_rebuild_decisions[sample.sample_id] = "rebuilt_from_raw"

        extracted = cluster.extract_read_windows(
            hdf5_file=extract_path,
            motifs=motif_list,
            regions=matched_regions,
        )
        control_result = None
        if signal_normalization == "control_regions":
            control_regions = sample.regions_bed or matched_regions
            if control_regions is None:
                raise ValueError(
                    "signal_normalization='control_regions' requires sample.regions_bed "
                    "or matched_regions."
                )
            control_artifact = _requested_extract_artifact(
                sample,
                motifs=motif_list,
                matched_regions=control_regions,
                signal_normalization=signal_normalization,
                feature_scaling=feature_scaling,
                cluster_basis=cluster_basis,
            )
            resolved_control_artifact = resolve_artifact(
                control_artifact,
                _coerce_artifacts(sample),
                artifact_policy=artifact_policy,
            )
            control_path = (
                resolved_control_artifact.path
                if resolved_control_artifact is not None
                else extract_path
            )
            control_result = cluster.extract_read_windows(
                hdf5_file=control_path,
                motifs=motif_list,
                regions=control_regions,
            )
        extracted, normalization_meta = _normalize_read_windows(
            extracted,
            signal_normalization=signal_normalization,
            control_result=control_result,
        )
        sample_normalization[sample.sample_id] = normalization_meta

        feature_matrix, feature_names = cluster.read_window_feature_matrix(extracted)
        feature_matrix, feature_names = _select_feature_columns(
            feature_matrix,
            feature_names,
            cluster_basis=cluster_basis,
        )
        if selected_feature_names is None:
            selected_feature_names = feature_names
        elif selected_feature_names != feature_names:
            raise ValueError("Feature names must match across all samples.")

        feature_blocks.append(feature_matrix)
        sample_dataset_sizes[sample.sample_id] = int(feature_matrix.shape[0])
        training_matrix, _ = _sample_training_rows(
            feature_matrix,
            training_sample_per_dataset=training_sample_per_dataset,
            random_state=random_state + sample_index,
        )
        training_blocks.append(training_matrix)
        sample_training_rows[sample.sample_id] = int(training_matrix.shape[0])

        for metadata in extracted.metadata:
            row = {
                "sample_id": sample.sample_id,
                "condition": sample.condition,
                "replicate": sample.replicate,
            }
            _merge_sample_metadata(row, sample.metadata)
            row.update(metadata)
            metadata_rows.append(row)

    if selected_feature_names is None:
        raise ValueError("No features were generated for the requested samples.")
    return _build_shared_cluster_result(
        mode=mode,
        motifs=motif_list,
        feature_blocks=feature_blocks,
        training_blocks=training_blocks,
        metadata_rows=metadata_rows,
        feature_names=list(selected_feature_names),
        signal_normalization=signal_normalization,
        feature_scaling=feature_scaling,
        cluster_basis=cluster_basis,
        clusterer=clusterer,
        n_clusters=n_clusters,
        training_sample_per_dataset=training_sample_per_dataset,
        artifact_policy=artifact_policy,
        random_state=random_state,
        make_plots=make_plots,
        matched_regions=matched_regions,
        sample_training_rows=sample_training_rows,
        sample_dataset_sizes=sample_dataset_sizes,
        sample_normalization=sample_normalization,
        cache_hits=cache_hits,
        cache_misses=cache_misses,
        sample_rebuild_decisions=sample_rebuild_decisions,
        quiet=quiet,
    )
