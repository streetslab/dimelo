from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .artifacts import resolve_artifact
from . import (
    cluster,
    distribution,
    plotting,
    region_analysis,
    region_contrasts,
    region_discovery,
)
from .models import (
    DatasetArtifact,
    RegionDiscoveryClusterContrastResult,
    RegionDiscoveryClusterResult,
    SampleSpec,
    SharedClusterModel,
    SharedClusterResult,
)

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
    val_matrix = None if result.val_matrix is None else np.asarray(result.val_matrix, dtype=float)
    offset_source = result if signal_normalization == "per_sample_global" else control_result
    if offset_source is None:
        raise ValueError("control_regions normalization requires control-region read windows.")
    offset_matrix = np.asarray(offset_source.data_matrix, dtype=float)
    offset_val_matrix = (
        None if offset_source.val_matrix is None else np.asarray(offset_source.val_matrix, dtype=float)
    )
    if offset_val_matrix is not None and offset_val_matrix.sum() > 0:
        global_offset = float(offset_matrix.sum() / offset_val_matrix.sum())
    else:
        global_offset = float(offset_matrix.mean())

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

    selected_names = [name for name, keep_name in zip(feature_names, keep) if keep_name]
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
    feature_frame: pd.DataFrame,
    assignments: pd.DataFrame,
) -> pd.DataFrame:
    profiles = feature_frame.copy()
    profiles["cluster"] = assignments["cluster"].to_numpy()
    aggregated = profiles.groupby("cluster", sort=True).mean(numeric_only=True).reset_index()
    counts = assignments.groupby("cluster", sort=True).size().reset_index(name="count")
    return counts.merge(aggregated, on="cluster", how="left")


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
    totals = summary.groupby(["region_id", "sample_id", "condition"])["count"].transform("sum")
    summary["fraction"] = summary["count"] / totals
    return summary


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
        start = int(getattr(row, "start"))
        end = int(getattr(row, "end"))
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


def discovery_cluster_workflow(
    *,
    samples: Iterable[SampleSpec],
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    discovery: dict[str, Any],
    clustering: dict[str, Any],
    selection: dict[str, Any] | None = None,
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
        raise ValueError("discovery_cluster_contrast_workflow requires contrasts['contrast'].")

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
    selected_region_spec = _selected_regions_to_region_spec(discovery_cluster_result.selected_regions)
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
        raise TypeError("Fitted clustering model does not support prediction for full assignment.")
    predicted_raw = _predict_in_chunks(clustering_result.model, full_scaled)
    predicted_ordered = cluster.apply_cluster_label_mapping(predicted_raw, label_mapping)

    assignments = pd.DataFrame(metadata_rows)
    assignments["cluster"] = _cluster_label_strings(predicted_ordered)

    feature_frame = pd.DataFrame(full_matrix, columns=feature_names)
    cluster_distribution = distribution.build_cluster_distribution(assignments)
    condition_distribution = distribution.build_condition_distribution(cluster_distribution)
    cluster_profiles = _cluster_profiles(feature_frame, assignments)
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

        for sample in sample_list:
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
        )
        sample_ids = [row["sample_id"] for row in metadata_rows]
        sample_dataset_sizes = {
            sample_id: int(sum(1 for row_sample in sample_ids if row_sample == sample_id))
            for sample_id in {sample.sample_id for sample in sample_list}
        }
        sample_normalization = {
            sample.sample_id: {"global_offset": None} for sample in sample_list
        }
        region_matrix = np.asarray(feature_matrix, dtype=float)
        if signal_normalization in {"per_sample_global", "control_regions"}:
            for sample_id in sample_dataset_sizes:
                mask = np.array([row["sample_id"] == sample_id for row in metadata_rows], dtype=bool)
                if signal_normalization == "control_regions":
                    sample = next(
                        sample_item for sample_item in sample_list if sample_item.sample_id == sample_id
                    )
                    control_regions = sample.regions_bed
                    if control_regions is None or control_regions == matched_regions:
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
                    )
                    offset = float(np.asarray(control_matrix, dtype=float).mean())
                else:
                    offset = float(region_matrix[mask].mean())
                region_matrix[mask] = region_matrix[mask] - offset
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
            mask = np.array(
                [row["sample_id"] == sample.sample_id for row in metadata_rows],
                dtype=bool,
            )
            sample_matrix = region_matrix[mask]
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
    )
