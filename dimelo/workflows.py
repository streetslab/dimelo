from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .artifacts import resolve_artifact
from . import cluster, distribution, plotting
from .models import DatasetArtifact, SampleSpec, SharedClusterModel, SharedClusterResult

_SUPPORTED_SIGNAL_NORMALIZATION = {"none", "per_sample_global", "control_regions"}
_SUPPORTED_FEATURE_SCALING = {"none", "robust_zscore"}
_SUPPORTED_CLUSTER_BASIS = {"shape_only", "shape_plus_level", "level_only"}
_PACKAGE_VERSION = "1.0.0"
_PREDICTION_CHUNK_SIZE = 100_000
_LEVEL_FEATURES = {
    "global_mean",
    "global_var",
    "global_median",
    "q25",
    "q75",
    "iqr",
    "global_mod_fraction",
}


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


def shared_cluster_distribution(
    *,
    samples: Iterable[SampleSpec],
    mode: str,
    motifs: Iterable[str],
    matched_regions: str | None = None,
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
    if mode != "read_global":
        raise NotImplementedError("The first workflow slice implements mode='read_global' only.")
    if len(sample_list) < 2:
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

    feature_frame = pd.DataFrame(full_matrix, columns=selected_feature_names)
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

    model = SharedClusterModel(
        mode=mode,
        motifs=motif_list,
        feature_names=list(selected_feature_names),
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
        region_summaries=None,
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
