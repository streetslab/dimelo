from __future__ import annotations

import contextlib
import math
import warnings
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from . import load_processed, utils

try:  # Optional SciPy dependency for advanced peak metrics
    from scipy.signal import find_peaks
except Exception:  # pragma: no cover - SciPy is optional
    find_peaks = None

try:  # Optional dependency for density-based clustering
    import hdbscan as _hdbscan_lib

    _HAS_HDBSCAN = True
except Exception:  # pragma: no cover - HDBSCAN is optional
    _HAS_HDBSCAN = False

RegionMetadata = tuple[str, int, int, str]
# Read metadata fields pulled alongside feature matrices by default
DEFAULT_METADATA_FIELDS: tuple[str, ...] = (
    "read_name",
    "chromosome",
    "region_start",
    "region_end",
    "motif",
    "region_strand",
    "read_length",
)

DEFAULT_AUTOCORR_LAGS = (20, 75, 100, 150, 250)
DEFAULT_DENSITY_WINDOWS = (
    ("density_pm50", -50, 50),
    ("density_pm75", -75, 75),
    ("density_pm100", -100, 100),
    ("density_pm300", -300, 300),
    ("density_pm500", -500, 500),
    ("density_plus240", 140, 340),
    ("density_minus240", -340, -140),
    ("density_plus320", 220, 420),
    ("density_minus320", -420, -220),
    ("density_plus800", 700, 900),
    ("density_minus800", -900, -700),
    ("density_plus150", 50, 250),
    ("density_minus150", -250, -50),
)


@dataclass
class ReadWindowExtractionConfig:
    # Legacy-style half-window around region center (full span = 2 * window_size)
    window_size: int | None = None
    orientation_aware: bool = True
    filter_multi_region_reads: bool = False
    # Require thresholded/binary read vectors for downstream clustering/classification.
    # If raw probabilities are detected, apply this threshold automatically unless set to None.
    auto_threshold_if_raw: float | int | None = 190
    enforce_thresholded_vectors: bool = True
    warn_on_auto_threshold: bool = True


@dataclass
class ReadWindowExtractionResult:
    data_matrix: np.ndarray
    val_matrix: np.ndarray | None
    metadata: list[dict[str, Any]]
    datasets: list[str]
    regions_dict: dict | None


@dataclass
class ClusterResult:
    labels_raw: np.ndarray
    labels_size_ordered: np.ndarray
    model: Any
    metrics: dict[str, float | None]


def region_feature_matrix_from_pileup(
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = True,
    *,
    pseudo_count: float = 1e-3,
    chunk_size: int = load_processed.DEFAULT_CHUNK_SIZE,
    split_large_regions: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    regions_executor=None,
) -> tuple[np.ndarray, list[RegionMetadata]]:
    """
    Convert every region in a pileup into a vector of modification fractions suitable for clustering.
    Uses pileup loaders so we never touch raw BAMs here.

    Args:
        bedmethyl_file: path to the processed bedmethyl file produced by parse_bam.pileup
        motif: motif string (e.g. "A,0") to pull from the pileup
        regions: region specifier accepted by load_processed (bed path, list, or str)
        window_size: optional window applied to every region before loading
        single_strand: whether to only use reads on the strand specified per region
        regions_5to3prime: align all regions to the 5'->3' direction before aggregation
        pseudo_count: small constant used to avoid division by zero when computing fractions
        chunk_size: chunk size forwarded to load_processed.pileup_vectors_from_bedmethyl
        split_large_regions: mirror of regions_to_list parameter for large regions
        quiet: disable progress bars in downstream loader
        cores: number of worker processes for loaders

    Returns:
        Tuple containing:
            * 2D numpy array (n_regions, region_length) of modification fractions
            * ordered list with metadata for every region (chrom, start, end, strand)
    """

    if regions is None:
        raise ValueError("regions must be provided when building a clustering matrix.")

    regions_dict = utils.regions_dict_from_input(regions, window_size)
    region_metadata: list[RegionMetadata] = []
    region_strings: list[str] = []
    for chromosome, region_list in regions_dict.items():
        for start, end, strand in region_list:
            region_metadata.append((chromosome, start, end, strand))
            region_strings.append(f"{chromosome}:{start}-{end},{strand}")
    if len(region_metadata) == 0:
        raise ValueError("No regions were found to build clustering features.")

    loader = partial(
        _pileup_fraction_vector_from_bedmethyl,
        bedmethyl_file=bedmethyl_file,
        motif=motif,
        window_size=window_size,
        single_strand=single_strand,
        regions_5to3prime=regions_5to3prime,
        pseudo_count=pseudo_count,
        quiet=quiet,
        cores=cores,
        chunk_size=chunk_size,
    )
    per_region_vectors = load_processed.regions_to_list(
        function_handle=loader,
        regions=region_strings,
        window_size=window_size,
        quiet=quiet,
        cores=cores,
        split_large_regions=split_large_regions,
        executor=regions_executor,
    )

    if len(per_region_vectors) != len(region_metadata):
        raise RuntimeError(
            "Region metadata and pileup vectors have mismatched lengths. "
            "Double-check the region inputs."
        )

    feature_vectors = []
    expected_length: int | None = None
    for vector in tqdm(per_region_vectors, desc="Processing regions", disable=quiet):
        # Backward-compatible path if upstream loaders still return (modified, valid).
        if isinstance(vector, tuple) and len(vector) == 2:
            modified = np.asarray(vector[0], dtype=np.float32)
            valid = np.asarray(vector[1], dtype=np.float32)
            if modified.shape != valid.shape:
                raise ValueError(
                    "Modified and valid vectors must have matching shapes."
                )
            denominator = valid + (2 * pseudo_count)
            with np.errstate(divide="ignore", invalid="ignore"):
                fraction = np.divide(
                    modified + pseudo_count,
                    denominator,
                    out=np.zeros_like(modified),
                    where=denominator > 0,
                )
        else:
            fraction = np.asarray(vector, dtype=np.float32)

        if expected_length is None:
            expected_length = fraction.shape[0]
        elif fraction.shape[0] != expected_length:
            raise ValueError(
                "All regions must have the same length. "
                "Pass window_size to enforce length-matched regions."
            )
        feature_vectors.append(fraction)

    feature_matrix = np.vstack(feature_vectors).astype(np.float32, copy=False)
    return feature_matrix, region_metadata


def _pileup_fraction_vector_from_bedmethyl(
    *,
    bedmethyl_file: str | Path,
    motif: str,
    regions: str | Path | list[str | Path],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    pseudo_count: float = 1e-3,
    quiet: bool = False,
    cores: int | None = None,
    chunk_size: int = load_processed.DEFAULT_CHUNK_SIZE,
) -> np.ndarray:
    modified_vec, valid_vec = load_processed.pileup_vectors_from_bedmethyl(
        bedmethyl_file=bedmethyl_file,
        motif=motif,
        regions=regions,
        window_size=window_size,
        single_strand=single_strand,
        regions_5to3prime=regions_5to3prime,
        quiet=quiet,
        cores=cores,
        chunk_size=chunk_size,
    )
    modified = np.asarray(modified_vec, dtype=np.float32)
    valid = np.asarray(valid_vec, dtype=np.float32)
    if modified.shape != valid.shape:
        raise ValueError("Modified and valid vectors must have matching shapes.")
    denominator = valid + (2 * pseudo_count)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(
            modified + pseudo_count,
            denominator,
            out=np.zeros_like(modified),
            where=denominator > 0,
        )


def read_mod_fraction_table(
    hdf5_file: str | Path,
    motifs: Sequence[str],
    regions: str | Path | list[str | Path] | None = None,
    *,
    window_size: int | None = None,
    single_strand: bool = False,
    sort_by: Sequence[str | tuple[str, str]] | str = (
        "chromosome",
        "region_start",
        "read_start",
    ),
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
    metadata_fields: Sequence[str] | None = DEFAULT_METADATA_FIELDS,
) -> tuple[np.ndarray, list[str], list[dict], dict | None]:
    """
    Pull per-read modification fractions from an extract .h5 file.
    This is the lightest-weight read-level feature matrix: one column per motif.

    Args:
        hdf5_file: path to the .h5 output from parse_bam.extract
        motifs: motifs whose modification fractions should be returned
        regions: optional region specifier to filter reads
        window_size: optional window applied to every region before loading
        single_strand: whether to enforce strand agreement with each region
        sort_by: forwarded to load_processed.read_vectors_from_hdf5
        subset_parameters: forwarded to utils.random_sample through the loader
        span_full_window: only keep reads fully covering the requested window
        metadata_fields: metadata columns to capture per read

    Returns:
        Tuple containing:
            * numpy array of shape (n_reads, len(motifs)) with per-read fractions
            * ordered list of feature names
            * list of metadata dictionaries aligned to the rows
            * regions_dict returned by load_processed.read_vectors_from_hdf5
    """

    read_records, datasets, regions_dict = load_processed.read_vectors_from_hdf5(
        file=hdf5_file,
        motifs=list(motifs),
        regions=regions,
        window_size=window_size,
        single_strand=single_strand,
        sort_by=list(sort_by) if isinstance(sort_by, (list, tuple)) else [sort_by],
        subset_parameters=subset_parameters,
        span_full_window=span_full_window,
    )

    dataset_indices = {name: idx for idx, name in enumerate(datasets)}
    feature_names = [f"{motif}_mod_fraction" for motif in motifs]
    missing_features = [name for name in feature_names if name not in dataset_indices]
    if missing_features:
        raise ValueError(
            "Missing modification fraction columns in HDF5 file: "
            f"{', '.join(missing_features)}. "
            "Ensure calculate_mod_fractions=True when exporting reads."
        )

    feature_matrix = np.zeros((len(read_records), len(feature_names)), dtype=float)
    capture_metadata = metadata_fields is not None
    metadata: list[dict] = []

    for idx, record in enumerate(read_records):
        feature_matrix[idx, :] = [
            float(record[dataset_indices[name]]) for name in feature_names
        ]
        if capture_metadata:
            metadata.append(
                {
                    field: record[dataset_indices[field]]
                    if field in dataset_indices
                    else None
                    for field in metadata_fields
                }
            )
        else:
            metadata.append({})

    return feature_matrix, feature_names, metadata, regions_dict


def cluster_features(
    feature_matrix: np.ndarray,
    *,
    method: str = "kmeans",
    n_clusters: int = 5,
    random_state: int | None = 0,
    **kwargs,
) -> tuple[np.ndarray, object]:
    """
    Run a simple clustering algorithm on a feature matrix.

    Args:
        feature_matrix: 2D array of shape (n_samples, n_features)
        method: clustering algorithm to run, currently only 'kmeans' is supported
        n_clusters: number of clusters to request from the estimator
        random_state: forwarded to the estimator for reproducibility
        kwargs: forwarded to the sklearn estimator

    Returns:
        Tuple of (cluster_labels, fitted_estimator)
    """

    if feature_matrix.ndim != 2:
        raise ValueError("feature_matrix must be 2-dimensional.")
    if feature_matrix.shape[0] == 0:
        raise ValueError("feature_matrix must contain at least one sample.")

    method = method.lower()
    if method != "kmeans":
        raise ValueError(f"Unsupported clustering method '{method}'.")

    kmeans_cls = _get_kmeans()
    estimator = kmeans_cls(
        n_clusters=n_clusters,
        random_state=random_state,
        **kwargs,
    )
    labels = estimator.fit_predict(feature_matrix)
    return labels, estimator


def _get_kmeans():
    try:
        from sklearn.cluster import KMeans
    except (
        ImportError
    ) as exc:  # pragma: no cover - exercised in environments w/o sklearn
        raise ImportError(
            "scikit-learn is required for k-means clustering. "
            "Install it with `pip install scikit-learn`."
        ) from exc
    return KMeans


def _get_xgb_classifier():
    try:
        from xgboost import XGBClassifier
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "Install xgboost to use classifier='xgboost'. Try `pip install xgboost`."
        ) from exc
    return XGBClassifier


def _build_dataset_index(dataset_names: Sequence[str]) -> dict[str, int]:
    # Map dataset names from the HDF5 file into tuple indices for quick lookup
    return {name: idx for idx, name in enumerate(dataset_names)}


def _coerce_strand(value: Any) -> str | None:
    # Normalize strand values to '+' / '-' / None
    if isinstance(value, str):
        stripped = value.strip()
        if stripped in {"+", "-"}:
            return stripped
    return None


def _should_flip(region_strand: str | None, read_strand: str | None) -> bool:
    return (
        region_strand in {"+", "-"}
        and read_strand in {"+", "-"}
        and region_strand != read_strand
    )


def _region_key(record: tuple, idx: dict[str, int]) -> tuple[Any, Any, Any]:
    # Defines "same region" for multi-region filtering
    chrom = record[idx["chromosome"]] if "chromosome" in idx else None
    start = record[idx["region_start"]] if "region_start" in idx else None
    end = record[idx["region_end"]] if "region_end" in idx else None
    return (chrom, start, end)


def _identify_multi_region_reads(
    records: Sequence[tuple], idx: dict[str, int]
) -> set[str]:
    if "read_name" not in idx:
        return set()
    regions_by_name: defaultdict[str, set[tuple[Any, Any, Any]]] = defaultdict(set)
    for rec in records:
        name = rec[idx["read_name"]]
        regions_by_name[name].add(_region_key(rec, idx))
    return {name for name, keys in regions_by_name.items() if len(keys) > 1}


def _infer_window_span_from_records(
    records: Sequence[tuple],
    idx: dict[str, int],
) -> int:
    if "region_start" not in idx or "region_end" not in idx:
        raise ValueError(
            "Could not infer window span from records. Pass window_size explicitly."
        )
    lengths = []
    for rec in records:
        try:
            start = int(rec[idx["region_start"]])
            end = int(rec[idx["region_end"]])
        except Exception:
            continue
        length = end - start
        if length > 0:
            lengths.append(length)
    if not lengths:
        raise ValueError(
            "Could not infer a default window span from region metadata. "
            "Pass window_size explicitly."
        )
    return int(min(lengths))


def _infer_window_span_from_metadata(
    metadata: Sequence[dict[str, Any]] | None,
) -> int | None:
    if metadata is None:
        return None
    lengths = []
    for row in metadata:
        try:
            start = int(row.get("region_start"))
            end = int(row.get("region_end"))
        except Exception:
            continue
        width = end - start
        if width > 0:
            lengths.append(width)
    if not lengths:
        return None
    return int(min(lengths))


def _resolve_window_size(window_size: int | None = None) -> int | None:
    if window_size is None:
        return None
    resolved = int(window_size)
    if resolved <= 0:
        raise ValueError("window_size must be a positive integer when provided.")
    return resolved


def _window_span_from_size(window_size: int) -> int:
    return int(window_size) * 2


def _centered_x_axis(length: int, span_bp: int | None) -> np.ndarray:
    if length <= 0:
        return np.array([], dtype=float)
    if length == 1:
        return np.array([0.0], dtype=float)
    span = float(span_bp if span_bp is not None else length)
    step = span / float(length)
    return (np.arange(length, dtype=float) * step) - (span / 2.0)


def _smooth_profile_vector(
    values: np.ndarray,
    *,
    smoothing: str | None = None,
    smooth_win: int = 21,
    smooth_sigma: float = 6.0,
) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if smoothing is None or arr.size < 3:
        return arr
    mode = smoothing.lower()
    if mode not in {"gaussian", "boxcar"}:
        raise ValueError("smoothing must be None, 'gaussian', or 'boxcar'.")

    win = int(max(3, smooth_win))
    if win % 2 == 0:
        win += 1
    if win >= arr.size:
        win = arr.size - 1 if arr.size % 2 == 0 else arr.size
    if win < 3:
        return arr

    if mode == "gaussian":
        radius = win // 2
        x = np.arange(-radius, radius + 1, dtype=float)
        sigma = max(float(smooth_sigma), 1e-6)
        kernel = np.exp(-0.5 * (x / sigma) ** 2)
    else:
        kernel = np.ones(win, dtype=float)
    kernel = kernel / kernel.sum()

    pad = len(kernel) // 2
    padded = np.pad(arr, (pad, pad), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def _infer_shared_window_size(
    regions: Sequence[str | Path | Sequence[str | Path]],
) -> int:
    """
    Infer one half-window size that is safe across multiple region selectors.
    Returns half of the shortest positive region width across all inputs.
    """
    widths: list[int] = []
    for region_input in regions:
        regions_dict = utils.regions_dict_from_input(region_input, window_size=None)
        for region_list in regions_dict.values():
            for start, end, _ in region_list:
                width = int(end) - int(start)
                if width > 0:
                    widths.append(width)
    if not widths:
        raise ValueError("Could not infer shared window size from provided regions.")
    shortest_width = int(min(widths))
    if shortest_width < 2:
        raise ValueError(
            "Could not infer a usable shared window_size because the shortest region "
            f"width is {shortest_width} bp."
        )
    return shortest_width // 2


def _center_crop_matrix(matrix: np.ndarray, target_width: int) -> np.ndarray:
    arr = np.asarray(matrix)
    if arr.ndim != 2:
        raise ValueError("Expected a 2D matrix for center-cropping.")
    width = arr.shape[1]
    if width == target_width:
        return arr
    if width < target_width:
        raise ValueError(
            f"Cannot center-crop from width {width} to larger target_width {target_width}."
        )
    start = (width - target_width) // 2
    end = start + target_width
    return arr[:, start:end]


def merge_read_window_results(
    results: Sequence[ReadWindowExtractionResult],
    *,
    source_labels: Sequence[str] | None = None,
    align: str = "error",
) -> ReadWindowExtractionResult:
    """
    Merge multiple ReadWindowExtractionResult objects safely.

    Args:
        results: extraction results to concatenate
        source_labels: optional labels to append to metadata as ``source_label``
        align: one of:
            - 'error': require identical window widths (strict)
            - 'center_crop': center-crop all matrices to the smallest width before merging

    Returns:
        A merged ReadWindowExtractionResult.
    """
    if len(results) == 0:
        raise ValueError(
            "results must contain at least one ReadWindowExtractionResult."
        )
    if source_labels is not None and len(source_labels) != len(results):
        raise ValueError("source_labels length must match results length.")
    if align not in {"error", "center_crop"}:
        raise ValueError("align must be one of {'error', 'center_crop'}.")

    widths = [int(r.data_matrix.shape[1]) for r in results]
    unique_widths = sorted(set(widths))
    target_width = min(widths)
    if len(unique_widths) > 1 and align == "error":
        raise ValueError(
            "Read-window widths do not match across results: "
            f"{unique_widths}. Use a shared window_size upstream "
            "or call merge_read_window_results(..., align='center_crop')."
        )

    data_blocks: list[np.ndarray] = []
    val_blocks: list[np.ndarray] = []
    all_have_val = all(r.val_matrix is not None for r in results)
    metadata_merged: list[dict[str, Any]] = []
    datasets_merged: list[str] = []

    for i, res in enumerate(results):
        data = np.asarray(res.data_matrix)
        if align == "center_crop":
            data = _center_crop_matrix(data, target_width)
        data_blocks.append(data)

        if all_have_val and res.val_matrix is not None:
            val = np.asarray(res.val_matrix)
            if align == "center_crop":
                val = _center_crop_matrix(val, target_width)
            val_blocks.append(val)

        label = source_labels[i] if source_labels is not None else None
        for row in res.metadata:
            if label is None:
                metadata_merged.append(dict(row))
            else:
                metadata_merged.append(dict(row, source_label=label))

        if not datasets_merged and res.datasets:
            datasets_merged = list(res.datasets)

    merged_data = np.vstack(data_blocks)
    merged_val = np.vstack(val_blocks) if all_have_val and len(val_blocks) > 0 else None

    return ReadWindowExtractionResult(
        data_matrix=merged_data,
        val_matrix=merged_val,
        metadata=metadata_merged,
        datasets=datasets_merged,
        regions_dict=None,
    )


def _prepare_group_labels(
    labels: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    labels_arr = np.asarray(labels)
    codes, uniques = pd.factorize(labels_arr, sort=True)
    # Handle NaN/None labels from factorize (-1 code): convert to explicit category name
    if np.any(codes < 0):
        labels_arr = labels_arr.astype(object)
        labels_arr[codes < 0] = "NA"
        codes, uniques = pd.factorize(labels_arr, sort=True)
    unique_codes = np.unique(codes)
    return labels_arr, codes, unique_codes


def _resolve_motif_slices(
    width: int,
    *,
    motif_count: int | None = None,
    motif_labels: Sequence[str] | None = None,
    window_size: int | None = None,
    view_window_size: int | None = None,
) -> tuple[int, int]:
    window_span = _window_span_from_size(window_size) if window_size else None
    view_span = _window_span_from_size(view_window_size) if view_window_size else None
    inferred_count: int | None = None
    if motif_count is not None:
        inferred_count = int(max(1, motif_count))
    elif motif_labels is not None and len(motif_labels) > 0:
        inferred_count = int(len(motif_labels))
    elif window_span and window_span > 0 and width % window_span == 0:
        inferred_count = width // window_span
    elif view_span and view_span > 0 and width % view_span == 0:
        inferred_count = width // view_span
    n_motifs = inferred_count or 1
    slice_width = width // n_motifs if n_motifs > 0 else width
    if slice_width <= 0:
        return 1, width
    return n_motifs, slice_width


def _extract_window_from_record(
    record: tuple,
    idx: dict[str, int],
    window_span: int,
    orientation_aware: bool,
    *,
    enforce_thresholded_vectors: bool,
    auto_threshold_if_raw: float | None,
    notify_auto_threshold,
) -> tuple[np.ndarray, np.ndarray | None] | None:
    required = ["read_start", "read_end", "region_start", "region_end", "mod_vector"]
    for key in required:
        if key not in idx:
            raise ValueError(f"read_vectors_from_hdf5 output missing '{key}' dataset.")

    read_start = int(record[idx["read_start"]])
    read_end = int(record[idx["read_end"]])
    region_start = int(record[idx["region_start"]])
    region_end = int(record[idx["region_end"]])
    mod_vector = np.asarray(record[idx["mod_vector"]])
    val_vector = np.asarray(record[idx["val_vector"]]) if "val_vector" in idx else None

    if read_end <= read_start or region_end <= region_start:
        return None

    center = (region_start + region_end) // 2
    half = window_span // 2
    # Grab a symmetric window around the region center
    window_start = center - half
    window_end = window_start + window_span
    # Skip reads that do not fully span the requested window
    if window_start < read_start or window_end > read_end:
        return None

    slice_start = window_start - read_start
    slice_end = slice_start + window_span
    if slice_start < 0 or slice_end > len(mod_vector):
        return None

    mod_window = mod_vector[slice_start:slice_end]
    val_window = val_vector[slice_start:slice_end] if val_vector is not None else None

    if enforce_thresholded_vectors and not _is_binary_mod_window(
        mod_window=mod_window,
        val_window=val_window,
    ):
        if auto_threshold_if_raw is None:
            raise ValueError(
                "Loaded extract vectors appear unthresholded (ML/probability values), "
                "but thresholded vectors are required. "
                "Set ReadWindowExtractionConfig.auto_threshold_if_raw to a value "
                "(for example 190) or disable enforcement explicitly."
            )
        mod_window = _threshold_mod_window(
            mod_window=mod_window,
            val_window=val_window,
            threshold=auto_threshold_if_raw,
        )
        notify_auto_threshold(auto_threshold_if_raw)

    if orientation_aware:
        region_strand = (
            _coerce_strand(record[idx["region_strand"]])
            if "region_strand" in idx
            else None
        )
        read_strand = _coerce_strand(record[idx["strand"]]) if "strand" in idx else None
        # If read/reference strands disagree, flip to align everything 5'->3'
        if _should_flip(region_strand, read_strand):
            mod_window = np.flip(mod_window)
            if val_window is not None:
                val_window = np.flip(val_window)

    return mod_window, val_window


def _resolve_raw_vector_threshold(threshold: float | int | None) -> float | None:
    if threshold is None:
        return None
    try:
        numeric = float(threshold)
    except Exception as exc:
        raise ValueError("auto_threshold_if_raw must be numeric or None.") from exc
    if numeric <= 0:
        raise ValueError("auto_threshold_if_raw must be > 0 when provided.")
    return float(utils.adjust_threshold(numeric, quiet=True))


def _is_binary_mod_window(
    mod_window: np.ndarray, val_window: np.ndarray | None
) -> bool:
    arr = np.asarray(mod_window)
    if arr.size == 0:
        return True
    if arr.dtype == np.bool_:
        return True
    if val_window is not None:
        valid = np.asarray(val_window).astype(bool)
        if valid.shape == arr.shape:
            arr = arr[valid]
            if arr.size == 0:
                return True
    unique_vals = np.unique(arr)
    return bool(np.all(np.isin(unique_vals, [0, 1])))


def _threshold_mod_window(
    *,
    mod_window: np.ndarray,
    val_window: np.ndarray | None,
    threshold: float,
) -> np.ndarray:
    arr = np.asarray(mod_window, dtype=float)
    if val_window is None:
        return (arr >= threshold).astype(float)
    valid = np.asarray(val_window).astype(bool)
    out = np.zeros(arr.shape, dtype=float)
    out[valid] = (arr[valid] >= threshold).astype(float)
    return out


def build_multimotif_read_windows(
    hdf5_file: str | Path,
    motifs: Sequence[str],
    regions: str | Path | list[str | Path] | None = None,
    *,
    window_size: int | None = None,
    orientation_aware: bool = True,
    single_strand: bool = False,
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
    require_all_motifs: bool = True,
    enforce_thresholded_vectors: bool = True,
    auto_threshold_if_raw: float | int | None = 190,
    warn_on_auto_threshold: bool = True,
) -> ReadWindowExtractionResult:
    """
    Group per-motif rows by read and return a combined window per read, concatenating motifs.

    Each requested motif contributes a centered window of length `2 * window_size`
    (or inferred full span when window_size is None); missing motifs are filled with zeros.
    Windows are concatenated in the order provided by `motifs`.

    Args:
        hdf5_file: path to extract .h5 file
        motifs: motifs to include (order matters for concatenation)
        regions: optional region filter
        window_size: half-window in bp around region center to extract per motif.
            Full extracted span is ``2 * window_size``.
            If None, full span is inferred from the shortest selected region length.
        orientation_aware: flip windows if read strand != region strand
        single_strand: passed to loader
        subset_parameters: passed to loader for subsetting
        span_full_window: passed to loader (if True, only reads spanning the region are loaded)
        require_all_motifs: if True, drop reads that are missing any requested motif
        enforce_thresholded_vectors: require thresholded/binary vectors for clustering-style workflows
        auto_threshold_if_raw: threshold applied when raw probability vectors are detected
            (default 190 interpreted in 0-255 space)
        warn_on_auto_threshold: emit a warning when auto-thresholding is applied

    Returns:
        ReadWindowExtractionResult with data_matrix shape
        ``(n_reads, len(motifs) * full_window_span_bp)``.
        and val_matrix if available, metadata per combined read, and datasets info.
    """

    # Load all per-motif rows
    read_tuples, dataset_names, regions_dict = load_processed.read_vectors_from_hdf5(
        file=hdf5_file,
        motifs=list(motifs),
        regions=regions,
        window_size=None,  # we handle windowing here
        single_strand=single_strand,
        subset_parameters=subset_parameters,
        span_full_window=span_full_window,
    )
    idx = _build_dataset_index(dataset_names)

    # Group rows by read+region key
    groups: defaultdict[tuple[Any, Any, Any, Any, Any], dict[str, tuple]] = defaultdict(
        dict
    )
    for rec in read_tuples:
        key = (
            rec[idx.get("read_name")],
            rec[idx.get("chromosome")],
            rec[idx.get("region_start")],
            rec[idx.get("region_end")],
            rec[idx.get("region_strand")] if "region_strand" in idx else None,
        )
        motif = rec[idx.get("motif")]
        groups[key][motif] = rec

    combined_windows: list[np.ndarray] = []
    combined_vals: list[np.ndarray] = []
    metadata: list[dict[str, Any]] = []
    has_val = "val_vector" in idx

    requested_window_size = _resolve_window_size(window_size=window_size)
    requested_window_span = (
        _window_span_from_size(requested_window_size)
        if requested_window_size is not None
        else None
    )
    effective_window_span = (
        requested_window_span
        if requested_window_span is not None
        else _infer_window_span_from_records(read_tuples, idx)
    )
    resolved_threshold = _resolve_raw_vector_threshold(auto_threshold_if_raw)
    auto_threshold_notified = False

    def _notify_auto_threshold(thresh_scaled: float) -> None:
        nonlocal auto_threshold_notified
        if auto_threshold_notified or not warn_on_auto_threshold:
            return
        auto_threshold_notified = True
        threshold_255 = int(round(thresh_scaled * 255))
        warnings.warn(
            "Detected unthresholded extract vectors; applying automatic thresholding "
            f"at {thresh_scaled:.6f} (~{threshold_255}/255) for clustering/classification.",
            RuntimeWarning,
            stacklevel=2,
        )

    for key, motif_map in groups.items():
        motif_windows = []
        motif_val_windows = []
        motifs_present = []
        representative_record: tuple | None = None
        for motif in motifs:
            if motif in motif_map:
                if representative_record is None:
                    representative_record = motif_map[motif]
                extracted = _extract_window_from_record(
                    motif_map[motif],
                    idx,
                    effective_window_span,
                    orientation_aware,
                    enforce_thresholded_vectors=enforce_thresholded_vectors,
                    auto_threshold_if_raw=resolved_threshold,
                    notify_auto_threshold=_notify_auto_threshold,
                )
                if extracted is None:
                    # If window cannot be extracted, treat as missing
                    motif_windows.append(np.zeros(effective_window_span, dtype=float))
                    motif_val_windows.append(
                        np.zeros(effective_window_span, dtype=float)
                        if has_val
                        else None
                    )
                else:
                    mw, vw = extracted
                    motif_windows.append(np.asarray(mw, dtype=float))
                    motif_val_windows.append(
                        np.asarray(vw, dtype=float) if vw is not None else None
                    )
                    motifs_present.append(motif)
            else:
                motif_windows.append(np.zeros(effective_window_span, dtype=float))
                motif_val_windows.append(
                    np.zeros(effective_window_span, dtype=float) if has_val else None
                )

        if require_all_motifs and len(motifs_present) < len(motifs):
            continue

        combined_windows.append(np.concatenate(motif_windows, axis=0))
        if has_val:
            combined_vals.append(
                np.concatenate(
                    [
                        vw
                        if vw is not None
                        else np.zeros(effective_window_span, dtype=float)
                        for vw in motif_val_windows
                    ],
                    axis=0,
                )
            )
        metadata.append(
            {
                "read_name": key[0],
                "chromosome": key[1],
                "read_start": (
                    int(representative_record[idx["read_start"]])
                    if representative_record is not None and "read_start" in idx
                    else None
                ),
                "read_end": (
                    int(representative_record[idx["read_end"]])
                    if representative_record is not None and "read_end" in idx
                    else None
                ),
                "read_length": (
                    int(representative_record[idx["read_length"]])
                    if representative_record is not None and "read_length" in idx
                    else None
                ),
                "region_start": key[2],
                "region_end": key[3],
                "region_strand": key[4],
                "motifs_present": motifs_present,
            }
        )

    if not combined_windows:
        raise ValueError(
            "No reads produced combined motif windows; check inputs and window_size."
        )

    data_matrix = np.vstack(combined_windows)
    val_matrix = np.vstack(combined_vals) if has_val and combined_vals else None

    return ReadWindowExtractionResult(
        data_matrix=data_matrix,
        val_matrix=val_matrix,
        metadata=metadata,
        datasets=list(dataset_names),
        regions_dict=regions_dict,
    )


def extract_read_windows(
    hdf5_file: str | Path,
    motifs: Sequence[str],
    regions: str | Path | list[str | Path] | None = None,
    *,
    config: ReadWindowExtractionConfig | None = None,
    window_size: int | None = None,
    single_strand: bool = False,
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
    quiet: bool = True,
) -> ReadWindowExtractionResult:
    """
    Extract fixed-length windows from single-read vectors, optionally flipping reads
    to align with the reference strand.

    Args:
        hdf5_file: path to the extract .h5 file
        motifs: motifs to pull from the file
        regions: optional region specifier
        config: controls windowing, orientation handling, and raw-vector threshold enforcement
        window_size: half-window override in bp (full span = 2 * window_size);
            when omitted, falls back to config.window_size.
            When omitted in both places, full span is inferred from the shortest selected region length.

    Returns:
        ReadWindowExtractionResult containing thresholded/binary mod windows by default,
        plus val matrices and metadata.
    """

    cfg = config or ReadWindowExtractionConfig()
    requested_window_size = _resolve_window_size(
        window_size=window_size if window_size is not None else cfg.window_size,
    )
    requested_window_span = (
        _window_span_from_size(requested_window_size)
        if requested_window_size is not None
        else None
    )
    effective_window_span = requested_window_span
    resolved_threshold = _resolve_raw_vector_threshold(cfg.auto_threshold_if_raw)
    auto_threshold_notified = False

    def _notify_auto_threshold(thresh_scaled: float) -> None:
        nonlocal auto_threshold_notified
        if auto_threshold_notified or not cfg.warn_on_auto_threshold:
            return
        auto_threshold_notified = True
        threshold_255 = int(round(thresh_scaled * 255))
        warnings.warn(
            "Detected unthresholded extract vectors; applying automatic thresholding "
            f"at {thresh_scaled:.6f} (~{threshold_255}/255) for clustering/classification.",
            RuntimeWarning,
            stacklevel=2,
        )

    # Load all requested reads/vectors from extract output
    read_tuples, dataset_names, regions_dict = load_processed.read_vectors_from_hdf5(
        file=hdf5_file,
        motifs=list(motifs),
        regions=regions,
        window_size=None,
        single_strand=single_strand,
        subset_parameters=subset_parameters,
        span_full_window=span_full_window,
    )
    idx = _build_dataset_index(dataset_names)
    if effective_window_span is None:
        effective_window_span = _infer_window_span_from_records(read_tuples, idx)
    effective_window_span = int(effective_window_span)

    drop_names: set[str] = set()
    if cfg.filter_multi_region_reads:
        drop_names = _identify_multi_region_reads(read_tuples, idx)

    matrices: list[np.ndarray] = []
    val_matrices: list[np.ndarray] = []
    metadata: list[dict[str, Any]] = []
    has_val = "val_vector" in idx

    for rec in read_tuples:
        read_name = rec[idx["read_name"]] if "read_name" in idx else None
        if read_name in drop_names:
            continue

        extracted = _extract_window_from_record(
            rec,
            idx,
            effective_window_span,
            cfg.orientation_aware,
            enforce_thresholded_vectors=cfg.enforce_thresholded_vectors,
            auto_threshold_if_raw=resolved_threshold,
            notify_auto_threshold=_notify_auto_threshold,
        )
        if extracted is None:
            continue

        mod_window, val_window = extracted
        matrices.append(mod_window)
        if has_val and val_window is not None:
            val_matrices.append(val_window)

        # Pull a small set of fields that are helpful for QC/relabelling
        meta_fields = [
            "read_name",
            "chromosome",
            "read_start",
            "read_end",
            "read_length",
            "region_start",
            "region_end",
            "region_strand",
            "strand",
            "motif",
        ]
        metadata.append(
            {field: rec[idx[field]] if field in idx else None for field in meta_fields}
        )

    if not matrices:
        raise ValueError("No reads produced a full window; try reducing window_size.")

    data_matrix = np.vstack([row.astype(float) for row in matrices])
    val_matrix = (
        np.vstack([row.astype(float) for row in val_matrices])
        if has_val and val_matrices
        else None
    )

    return ReadWindowExtractionResult(
        data_matrix=data_matrix,
        val_matrix=val_matrix,
        metadata=metadata,
        datasets=list(dataset_names),
        regions_dict=regions_dict,
    )


def compute_autocorrelation_feature(vec: np.ndarray, lag: int) -> float:
    # Normalized autocorrelation at a given lag (handles short/flat vectors)
    arr = np.asarray(vec, dtype=float)
    if arr.size <= abs(lag):
        return 0.0
    centered = arr - arr.mean()
    denom = np.dot(centered, centered)
    if denom == 0:
        return 0.0
    if lag >= 0:
        return float(np.dot(centered[:-lag], centered[lag:]) / denom)
    return float(np.dot(centered[-lag:], centered[:lag]) / denom)


def read_window_feature_matrix(
    result: ReadWindowExtractionResult,
    *,
    n_pca: int = 6,
    random_state: int = 42,
    autocorr_lags: Sequence[int] = DEFAULT_AUTOCORR_LAGS,
    density_windows: Sequence[tuple[str, int, int]] = DEFAULT_DENSITY_WINDOWS,
    peak_prominence: float = 0.005,
    use_peak_features: bool = False,
    require_nonzero_valid: bool = False,
    min_valid_fraction: float = 0.0,
) -> tuple[np.ndarray, list[str]]:
    """
    Convert read windows into an augmented feature matrix including PCA components,
    autocorrelation, density summaries, and basic statistics.

    Set use_peak_features=True to append peak counts/prominences (requires SciPy).
    If require_nonzero_valid is True, rows with no valid sites (or below min_valid_fraction)
    are dropped before feature computation (requires val_matrix to be present).
    """

    data_matrix = result.data_matrix.astype(float, copy=False)
    val_matrix = result.val_matrix

    if require_nonzero_valid and val_matrix is not None:
        valid_sums = val_matrix.sum(axis=1)
        mask = valid_sums > 0
        if min_valid_fraction > 0 and val_matrix.shape[1] > 0:
            mask &= (valid_sums / val_matrix.shape[1]) >= min_valid_fraction
        data_matrix = data_matrix[mask]
        val_matrix = val_matrix[mask]
    if data_matrix.size == 0:
        raise ValueError("No rows remaining after filtering for valid sites.")

    n_reads, window_size = data_matrix.shape
    center = window_size // 2

    columns: list[np.ndarray] = []
    names: list[str] = []

    n_components = min(n_pca, n_reads, window_size)
    if n_components > 0:
        from sklearn.decomposition import PCA

        pca = PCA(n_components=n_components, random_state=random_state)
        pca_vals = pca.fit_transform(data_matrix)
        columns.append(pca_vals)
        names.extend([f"pca_{i}" for i in range(n_components)])

    # Precompute cumulative sums to accelerate density windows
    cumsum = np.cumsum(data_matrix, axis=1)
    cumsum = np.pad(cumsum, ((0, 0), (1, 0)), mode="constant", constant_values=0)

    for lag in autocorr_lags:
        values = np.array(
            [compute_autocorrelation_feature(row, lag) for row in data_matrix]
        )
        columns.append(values[:, None])
        names.append(f"autocorr_{lag}")

    for label, start, end in density_windows:
        start_idx = max(0, center + start)
        end_idx = min(window_size, center + end)
        if end_idx <= start_idx:
            values = np.zeros(n_reads)
        else:
            length = end_idx - start_idx
            window_sum = cumsum[:, end_idx] - cumsum[:, start_idx]
            values = window_sum / length
        columns.append(values[:, None])
        names.append(label)

    global_mean = data_matrix.mean(axis=1)
    global_var = data_matrix.var(axis=1)
    global_median = np.median(data_matrix, axis=1)
    q25 = np.percentile(data_matrix, 25, axis=1)
    q75 = np.percentile(data_matrix, 75, axis=1)

    columns.extend(
        [
            global_mean[:, None],
            global_var[:, None],
            global_median[:, None],
            q25[:, None],
            q75[:, None],
            (q75 - q25)[:, None],
        ]
    )
    names.extend(
        [
            "global_mean",
            "global_var",
            "global_median",
            "q25",
            "q75",
            "iqr",
        ]
    )

    if val_matrix is not None:
        valid_sum = val_matrix.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.divide(
                data_matrix.sum(axis=1),
                valid_sum,
                out=np.zeros_like(valid_sum),
                where=valid_sum > 0,
            )
        columns.append(frac[:, None])
        names.append("global_mod_fraction")

    if use_peak_features and find_peaks is not None:  # pragma: no branch
        peak_counts = []
        peak_prominences = []
        for row in data_matrix:
            peaks, props = find_peaks(row, prominence=peak_prominence)
            peak_counts.append(len(peaks))
            if peaks.size > 0:
                peak_prominences.append(float(np.mean(props["prominences"])))
            else:
                peak_prominences.append(0.0)
        columns.append(np.array(peak_counts)[:, None])
        columns.append(np.array(peak_prominences)[:, None])
        names.extend(["peak_count", "peak_prominence"])

    feature_matrix = np.hstack(columns)
    return feature_matrix, names


def _safe_scores(X: np.ndarray, labels: np.ndarray) -> dict[str, float | None]:
    from sklearn.metrics import (
        calinski_harabasz_score,
        davies_bouldin_score,
        silhouette_score,
    )

    scores = {"silhouette": None, "calinski_harabasz": None, "davies_bouldin": None}
    unique_labels = np.unique(labels)
    valid = len(unique_labels) >= 2 and not (
        len(unique_labels) == 1 and unique_labels[0] == -1
    )
    if not valid:
        return scores
    with contextlib.suppress(Exception):
        scores["silhouette"] = float(silhouette_score(X, labels))
    with contextlib.suppress(Exception):
        scores["calinski_harabasz"] = float(calinski_harabasz_score(X, labels))
    with contextlib.suppress(Exception):
        scores["davies_bouldin"] = float(davies_bouldin_score(X, labels))
    return scores


def _renumber_by_size(labels: np.ndarray, noise_label: int | None = None) -> np.ndarray:
    labels = np.asarray(labels)
    counts = Counter(labels)
    order = [
        label
        for label, _ in sorted(counts.items(), key=lambda x: x[1], reverse=True)
        if noise_label is None or label != noise_label
    ]
    mapping = {old: new for new, old in enumerate(order)}
    remapped = []
    for label in labels:
        if noise_label is not None and label == noise_label:
            remapped.append(-1)
        else:
            remapped.append(mapping.get(label, -1))
    return np.array(remapped, dtype=int)


def cluster_label_mapping(
    labels_raw: np.ndarray,
    labels_size_ordered: np.ndarray,
) -> dict[int, int]:
    """
    Build a stable mapping from raw estimator labels to size-ordered labels.
    """

    mapping: dict[int, int] = {}
    for raw_label, ordered_label in zip(labels_raw, labels_size_ordered, strict=False):
        raw = int(raw_label)
        ordered = int(ordered_label)
        if raw in mapping and mapping[raw] != ordered:
            raise ValueError("Raw cluster labels map to multiple size-ordered labels.")
        mapping[raw] = ordered
    return mapping


def apply_cluster_label_mapping(
    labels: np.ndarray,
    mapping: dict[int, int],
    *,
    unknown_label: int = -1,
) -> np.ndarray:
    """
    Apply a raw-to-size-ordered label mapping to new assignments.
    """

    return np.array(
        [mapping.get(int(label), unknown_label) for label in np.asarray(labels)],
        dtype=int,
    )


def sample_rows(
    data: np.ndarray,
    labels: Sequence[Any] | None = None,
    *,
    n: int | None = None,
    frac: float | None = None,
    random_state: int | None = 42,
    stratify: bool = False,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
    """
    Downsample rows for faster clustering/classification.

    Args:
        data: 2D array of shape (n_rows, n_features)
        labels: optional parallel labels to return (e.g., cluster or sample labels)
        n: fixed number of rows to sample
        frac: fraction of rows to sample (ignored if n is provided)
        random_state: RNG seed
        stratify: if True and labels are provided, sample proportionally to label frequency

    Returns:
        sampled_data, sampled_labels (or None), indices of selected rows
    """

    rng = np.random.default_rng(random_state)
    n_rows = data.shape[0]
    if n is None and frac is None:
        raise ValueError("Provide n or frac to sample.")
    if n is None:
        n = max(1, int(n_rows * frac))
    n = min(n, n_rows)

    if stratify and labels is not None:
        labels_arr = np.asarray(labels)
        unique, counts = np.unique(labels_arr, return_counts=True)
        probs = counts / counts.sum()
        selected = []
        for lbl, p in zip(unique, probs, strict=False):
            k = max(1, int(round(p * n)))
            choices = np.flatnonzero(labels_arr == lbl)
            k = min(k, len(choices))
            selected.append(rng.choice(choices, size=k, replace=False))
        idx = np.unique(np.concatenate(selected))
        if len(idx) > n:
            idx = rng.choice(idx, size=n, replace=False)
    else:
        idx = rng.choice(n_rows, size=n, replace=False)

    sampled_data = data[idx]
    sampled_labels = np.asarray(labels)[idx] if labels is not None else None
    return sampled_data, sampled_labels, idx


def cluster_read_windows(
    feature_matrix: np.ndarray,
    *,
    method: str = "kmeans",
    n_clusters: int = 8,
    random_state: int = 42,
    auto_k: bool = False,
    k_grid: Sequence[int] | None = None,
    score: str = "silhouette",
    **kwargs,
) -> ClusterResult:
    """
    Run clustering on a read feature matrix with support for several algorithms.
    auto_k=True grid-searches K for methods that require it using the requested score.
    """

    from sklearn.cluster import (
        DBSCAN,
        OPTICS,
        AgglomerativeClustering,
        Birch,
        KMeans,
        MiniBatchKMeans,
        SpectralClustering,
    )
    from sklearn.mixture import GaussianMixture

    X = np.asarray(feature_matrix)

    def fit(
        n_components: int,
    ) -> tuple[np.ndarray, Any, dict[str, float | None], int | None]:
        algo = method.lower()
        model = None
        labels: np.ndarray
        noise_label: int | None = None

        if algo == "kmeans":
            model = KMeans(
                n_clusters=n_components,
                random_state=random_state,
                n_init=kwargs.get("n_init", 10),
            )
            labels = model.fit_predict(X)
        elif algo == "minibatch_kmeans":
            model = MiniBatchKMeans(
                n_clusters=n_components,
                random_state=random_state,
                n_init=kwargs.get("n_init", 10),
                batch_size=kwargs.get("batch_size", 1024),
            )
            labels = model.fit_predict(X)
        elif algo == "gmm":
            model = GaussianMixture(
                n_components=n_components,
                covariance_type=kwargs.get("covariance_type", "full"),
                random_state=random_state,
            )
            model.fit(X)
            labels = model.predict(X)
        elif algo == "agglomerative":
            model = AgglomerativeClustering(
                n_clusters=n_components,
                linkage=kwargs.get("linkage", "ward"),
                metric="euclidean",
            )
            labels = model.fit_predict(X)
        elif algo == "spectral":
            model = SpectralClustering(
                n_clusters=n_components,
                assign_labels=kwargs.get("assign_labels", "kmeans"),
                affinity=kwargs.get("affinity", "rbf"),
                random_state=random_state,
            )
            labels = model.fit_predict(X)
        elif algo == "birch":
            model = Birch(
                n_clusters=n_components, threshold=kwargs.get("threshold", 0.5)
            )
            labels = model.fit_predict(X)
        elif algo == "dbscan":
            model = DBSCAN(
                eps=kwargs.get("eps", 0.8),
                min_samples=kwargs.get("min_samples", 20),
            )
            labels = model.fit_predict(X)
            noise_label = -1
        elif algo == "optics":
            model = OPTICS(
                min_samples=kwargs.get("min_samples", 20),
                max_eps=kwargs.get("max_eps", np.inf),
                cluster_method="xi",
                xi=kwargs.get("xi", 0.05),
            )
            labels = model.fit_predict(X)
            noise_label = -1
        elif algo == "hdbscan":
            if not _HAS_HDBSCAN:
                raise ImportError("Install hdbscan to use method='hdbscan'.")
            model = _hdbscan_lib.HDBSCAN(
                min_cluster_size=kwargs.get("min_cluster_size", 30),
                min_samples=kwargs.get("min_samples"),
            )
            labels = model.fit_predict(X)
            noise_label = -1
        elif algo == "umap_kmeans":
            try:
                import umap
            except Exception as exc:  # pragma: no cover - optional dependency
                raise ImportError(
                    "Install umap-learn to use method='umap_kmeans'."
                ) from exc
            reducer = umap.UMAP(
                n_neighbors=kwargs.get("umap_n_neighbors", 15),
                min_dist=kwargs.get("umap_min_dist", 0.1),
                n_components=kwargs.get("umap_n_components", 2),
                metric=kwargs.get("umap_metric", "euclidean"),
                random_state=random_state,
            )
            embedded = reducer.fit_transform(X)
            kmeans_cls = _get_kmeans()
            model = kmeans_cls(
                n_clusters=n_components,
                random_state=random_state,
                n_init=kwargs.get("n_init", 10),
            )
            labels = model.fit_predict(embedded)
        else:
            raise ValueError(f"Unknown clustering method '{method}'.")

        metrics = _safe_scores(X, labels)
        return labels, model, metrics, noise_label

    if auto_k:
        grid = list(k_grid or range(2, 16))
        best: ClusterResult | None = None
        best_score: float | None = None

        def score_value(metric_dict: dict[str, float | None]) -> float | None:
            val = metric_dict.get(score)
            if val is None:
                return None
            if score == "davies_bouldin":
                return -val
            return val

        for k in grid:
            labels, model, metrics, noise_label = fit(k)
            val = score_value(metrics)
            if val is None:
                continue
            if best is None or val > best_score:
                best = ClusterResult(
                    labels_raw=labels,
                    labels_size_ordered=_renumber_by_size(labels, noise_label),
                    model=model,
                    metrics=metrics,
                )
                best_score = val
        if best is None:
            labels, model, metrics, noise_label = fit(grid[-1])
            best = ClusterResult(
                labels_raw=labels,
                labels_size_ordered=_renumber_by_size(labels, noise_label),
                model=model,
                metrics=metrics,
            )
        return best

    labels, model, metrics, noise_label = fit(n_clusters)
    return ClusterResult(
        labels_raw=labels,
        labels_size_ordered=_renumber_by_size(labels, noise_label),
        model=model,
        metrics=metrics,
    )


def plot_cluster_profiles(
    data_matrix: np.ndarray,
    labels: np.ndarray,
    *,
    val_matrix: np.ndarray | None = None,
    include_clusters: Sequence[int] | None = None,
    view_window_size: int | None = None,
    window_size: int | None = None,
    motif_index: int = 0,
    motif_count: int | None = None,
    metadata: Sequence[dict[str, Any]] | None = None,
    motif_labels: Sequence[str] | None = None,
    motif_colors: Sequence[str] | None = None,
    plot_all_motifs: bool = False,
    motif_profile_mode: str = "single_axis",
    color_points_by: str = "cluster",
    show_valid_profile: bool = False,
    signal_label: str = "Modification signal",
    valid_label: str = "Valid-site fraction",
    smoothing: str | None = None,
    smooth_win: int = 21,
    smooth_sigma: float = 6.0,
    show_unsmoothed_overlay: bool = False,
    cmap_name: str = "viridis",
    invert_y: bool = True,
    point_size: float = 1.0,
    point_alpha: float = 0.01,
):
    """
    Scatter plot of per-read modification calls with side profile summaries per label.
    Supports concatenated multi-motif windows and optional overlays per motif.
    """

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    X_full = np.asarray(data_matrix)
    labs = np.asarray(labels)
    if X_full.shape[0] != labs.shape[0]:
        raise ValueError("data_matrix and labels must have the same number of rows.")
    if include_clusters is not None:
        mask = np.isin(labs, list(include_clusters))
        if not np.any(mask):
            raise ValueError("No reads remain after filtering include_clusters.")
        X_full = X_full[mask]
        labs = labs[mask]
        if val_matrix is not None:
            val_matrix = np.asarray(val_matrix)[mask]
        if metadata is not None:
            metadata = [row for row, keep in zip(metadata, mask, strict=False) if keep]

    V_full = np.asarray(val_matrix) if val_matrix is not None else None
    if V_full is not None and V_full.shape != X_full.shape:
        raise ValueError("val_matrix must match data_matrix shape.")

    if color_points_by not in {"cluster", "motif"}:
        raise ValueError("color_points_by must be 'cluster' or 'motif'.")
    if motif_profile_mode not in {"single_axis", "separate_axes"}:
        raise ValueError("motif_profile_mode must be 'single_axis' or 'separate_axes'.")

    window_len = X_full.shape[1]
    n_motifs, slice_width = _resolve_motif_slices(
        window_len,
        motif_count=motif_count,
        motif_labels=motif_labels,
        window_size=window_size,
        view_window_size=view_window_size,
    )
    if motif_index < 0 or motif_index >= n_motifs:
        raise ValueError(
            f"motif_index {motif_index} out of range for {n_motifs} motif slices."
        )
    motif_ids = (
        list(range(n_motifs)) if plot_all_motifs and n_motifs > 1 else [motif_index]
    )

    if motif_labels is None:
        motif_labels = [f"motif_{i}" for i in range(n_motifs)]
    else:
        motif_labels = list(motif_labels)
        if len(motif_labels) < n_motifs:
            motif_labels.extend(
                [f"motif_{i}" for i in range(len(motif_labels), n_motifs)]
            )
    if motif_colors is None:
        motif_cmap = plt.get_cmap("tab10")
        motif_colors = [motif_cmap(i % 10) for i in range(n_motifs)]
    else:
        motif_colors = list(motif_colors)
        if len(motif_colors) < n_motifs:
            motif_cmap = plt.get_cmap("tab10")
            motif_colors.extend(
                motif_cmap(i % 10) for i in range(len(motif_colors), n_motifs)
            )

    selected_start = motif_index * slice_width
    selected_end = min(window_len, selected_start + slice_width)
    primary_matrix = X_full[:, selected_start:selected_end]
    inferred_window_span = (
        _window_span_from_size(window_size)
        if window_size is not None
        else primary_matrix.shape[1]
    )
    view_span = (
        _window_span_from_size(view_window_size)
        if view_window_size is not None
        else None
    )
    x_axis = _centered_x_axis(
        primary_matrix.shape[1], view_span or inferred_window_span
    )

    mod_fraction = primary_matrix.mean(axis=1)
    labs_arr, lab_codes, unique_codes = _prepare_group_labels(labs)
    order = np.lexsort((-mod_fraction, lab_codes))
    X_full = X_full[order]
    labs_arr = labs_arr[order]
    lab_codes = lab_codes[order]
    V_full = V_full[order] if V_full is not None else None

    unique_labels = np.array(
        [np.unique(labs_arr[lab_codes == code])[0] for code in unique_codes]
    )
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=int(unique_codes.min()), vmax=int(unique_codes.max() or 1))

    fig = plt.figure(figsize=(12, max(4, len(unique_codes))))
    gs = gridspec.GridSpec(nrows=len(unique_codes), ncols=2, width_ratios=[3, 1])
    ax_left = fig.add_subplot(gs[:, 0])
    if color_points_by == "cluster":
        rows, cols = np.nonzero(primary_matrix[order])
        colors = cmap(norm(lab_codes[rows]))
        ax_left.scatter(
            x_axis[cols],
            rows,
            s=point_size,
            alpha=point_alpha,
            c=colors,
            rasterized=True,
        )
    else:
        for motif_id in motif_ids:
            s = motif_id * slice_width
            e = min(window_len, s + slice_width)
            motif_matrix = X_full[:, s:e]
            rows, cols = np.nonzero(motif_matrix)
            ax_left.scatter(
                x_axis[cols],
                rows,
                s=point_size,
                alpha=point_alpha,
                c=motif_colors[motif_id],
                label=str(motif_labels[motif_id]),
                rasterized=True,
            )
    ax_left.set_xlabel("Distance from region center (bp)")
    ax_left.set_ylabel("Sorted read index")
    if invert_y:
        ax_left.invert_yaxis()

    change_points = np.flatnonzero(np.diff(lab_codes)) + 1
    for cp in change_points:
        ax_left.axhline(cp, color="0.2", linestyle="--", linewidth=0.3)

    for i, code in enumerate(unique_codes):
        ax = fig.add_subplot(gs[i, 1])
        row_mask = lab_codes == code
        panel_max = 0.0
        if motif_profile_mode == "separate_axes" and len(motif_ids) > 1:
            axes_for_motifs = [ax]
            motif_maxima = [0.0 for _ in motif_ids]
            for motif_pos in range(1, len(motif_ids)):
                twin = ax.twinx()
                twin.spines["right"].set_position(("outward", 35 * motif_pos))
                axes_for_motifs.append(twin)
            for motif_pos, motif_id in enumerate(motif_ids):
                s = motif_id * slice_width
                e = min(window_len, s + slice_width)
                motif_matrix = X_full[:, s:e]
                mean_profile_raw = motif_matrix[row_mask].mean(axis=0)
                mean_profile = _smooth_profile_vector(
                    mean_profile_raw,
                    smoothing=smoothing,
                    smooth_win=smooth_win,
                    smooth_sigma=smooth_sigma,
                )
                motif_maxima[motif_pos] = max(
                    motif_maxima[motif_pos], float(np.nanmax(mean_profile_raw))
                )
                panel_max = max(panel_max, motif_maxima[motif_pos])
                motif_ax = axes_for_motifs[motif_pos]
                if show_unsmoothed_overlay and smoothing is not None:
                    motif_ax.plot(
                        x_axis,
                        mean_profile_raw,
                        color=motif_colors[motif_id],
                        linewidth=1.0,
                        alpha=0.25,
                    )
                motif_ax.plot(
                    x_axis,
                    mean_profile,
                    color=motif_colors[motif_id],
                    linewidth=1.5,
                )
                motif_ax.set_ylabel(
                    str(motif_labels[motif_id]), color=motif_colors[motif_id]
                )
                motif_ax.tick_params(
                    axis="y", colors=motif_colors[motif_id], labelsize=8
                )
                if show_valid_profile and V_full is not None:
                    mean_val_raw = V_full[row_mask, s:e].mean(axis=0)
                    mean_val = _smooth_profile_vector(
                        mean_val_raw,
                        smoothing=smoothing,
                        smooth_win=smooth_win,
                        smooth_sigma=smooth_sigma,
                    )
                    motif_ax.plot(
                        x_axis,
                        mean_val,
                        color=motif_colors[motif_id],
                        linestyle="--",
                        linewidth=1.0,
                        alpha=0.85,
                    )
            for motif_pos, motif_ax in enumerate(axes_for_motifs):
                motif_ax.set_ylim(0, max(motif_maxima[motif_pos], 0.05) * 1.05)
        else:
            for motif_id in motif_ids:
                s = motif_id * slice_width
                e = min(window_len, s + slice_width)
                motif_matrix = X_full[:, s:e]
                mean_profile_raw = motif_matrix[row_mask].mean(axis=0)
                mean_profile = _smooth_profile_vector(
                    mean_profile_raw,
                    smoothing=smoothing,
                    smooth_win=smooth_win,
                    smooth_sigma=smooth_sigma,
                )
                panel_max = max(panel_max, float(np.nanmax(mean_profile_raw)))
                if show_unsmoothed_overlay and smoothing is not None:
                    ax.plot(
                        x_axis,
                        mean_profile_raw,
                        color=motif_colors[motif_id]
                        if len(motif_ids) > 1
                        else cmap(norm(code)),
                        linewidth=1.0,
                        alpha=0.25,
                    )
                ax.plot(
                    x_axis,
                    mean_profile,
                    color=motif_colors[motif_id]
                    if len(motif_ids) > 1
                    else cmap(norm(code)),
                    linewidth=1.5,
                )
                if show_valid_profile and V_full is not None:
                    mean_val_raw = V_full[row_mask, s:e].mean(axis=0)
                    mean_val = _smooth_profile_vector(
                        mean_val_raw,
                        smoothing=smoothing,
                        smooth_win=smooth_win,
                        smooth_sigma=smooth_sigma,
                    )
                    ax.plot(
                        x_axis,
                        mean_val,
                        color=motif_colors[motif_id] if len(motif_ids) > 1 else "0.35",
                        linestyle="--",
                        linewidth=1.0,
                        alpha=0.85,
                    )
        ax.set_title(f"{unique_labels[i]} (n={row_mask.sum()})")
        ax.set_xlim(x_axis[0], x_axis[-1])
        y_max = max(panel_max, 0.05)
        ax.set_ylim(0, y_max * 1.05)

    legend_handles: list[Line2D] = []
    if len(motif_ids) > 1 or color_points_by == "motif":
        legend_handles.extend(
            [
                Line2D(
                    [0],
                    [0],
                    color=motif_colors[m],
                    linewidth=2,
                    label=str(motif_labels[m]),
                )
                for m in motif_ids
            ]
        )
    if show_valid_profile and V_full is not None:
        legend_handles.extend(
            [
                Line2D(
                    [0],
                    [0],
                    color="black",
                    linewidth=2,
                    linestyle="-",
                    label=signal_label,
                ),
                Line2D(
                    [0],
                    [0],
                    color="black",
                    linewidth=1,
                    linestyle="--",
                    label=valid_label,
                ),
            ]
        )
    if legend_handles:
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            ncol=min(4, len(legend_handles)),
            frameon=False,
            bbox_to_anchor=(0.5, 1.01),
        )

    fig.tight_layout()
    return fig


def classify_read_features_binary(
    feature_matrix: np.ndarray,
    sample_labels: Sequence[str | int],
    *,
    classifier: str = "xgboost",
    test_size: float = 0.2,
    random_state: int = 42,
    **kwargs,
) -> dict[str, Any]:
    """
    Train/test a binary classifier on read-level features to quantify separability between two samples.

    Args:
        feature_matrix: 2D array (n_reads, n_features)
        sample_labels: iterable of length n_reads with exactly two distinct values (sample IDs)
        classifier: 'xgboost' (default) or one of:
            'logreg' (sklearn LogisticRegression),
            'sgd' (sklearn SGDClassifier),
            'random_forest' (sklearn RandomForestClassifier),
            'svc' (sklearn SVC with probability=True),
            'knn' (sklearn KNeighborsClassifier)
        test_size: fraction held out for testing
        random_state: RNG seed for reproducibility
        kwargs: forwarded to the classifier constructor

    Returns:
        dict with keys:
            model: fitted estimator
            metrics: accuracy, roc_auc, confusion_matrix
            predictions: DataFrame with columns ['true_label','pred_label','proba','split','sample_label']
            label_mapping: mapping from class name to encoded int
    """

    from sklearn.metrics import accuracy_score, confusion_matrix, roc_auc_score
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import LabelEncoder

    X = np.asarray(feature_matrix)
    y_raw = np.asarray(sample_labels)
    if X.shape[0] != y_raw.shape[0]:
        raise ValueError("feature_matrix and sample_labels must have matching length.")
    row_index = np.arange(X.shape[0], dtype=int)

    le = LabelEncoder()
    y = le.fit_transform(y_raw)
    classes = le.classes_
    if len(classes) != 2:
        raise ValueError("Binary classification requires exactly two sample labels.")

    X_train, X_test, y_train, y_test, lbl_train, lbl_test, idx_train, idx_test = (
        train_test_split(
            X,
            y,
            y_raw,
            row_index,
            test_size=test_size,
            random_state=random_state,
            stratify=y,
        )
    )

    clf_name = classifier.lower()
    if clf_name == "xgboost":
        XGB = _get_xgb_classifier()
        # Sensible defaults; kwargs can override
        model = XGB(
            n_estimators=kwargs.get("n_estimators", 200),
            max_depth=kwargs.get("max_depth", 4),
            learning_rate=kwargs.get("learning_rate", 0.1),
            subsample=kwargs.get("subsample", 0.8),
            colsample_bytree=kwargs.get("colsample_bytree", 0.8),
            random_state=random_state,
            n_jobs=kwargs.get("n_jobs", 4),
        )
    elif clf_name == "logreg":
        from sklearn.linear_model import LogisticRegression

        model = LogisticRegression(
            max_iter=kwargs.get("max_iter", 200),
            solver=kwargs.get("solver", "lbfgs"),
            class_weight=kwargs.get("class_weight", "balanced"),
            random_state=random_state,
        )
    elif clf_name == "sgd":
        from sklearn.linear_model import SGDClassifier

        model = SGDClassifier(
            loss=kwargs.get("loss", "log_loss"),
            max_iter=kwargs.get("max_iter", 1000),
            alpha=kwargs.get("alpha", 1e-4),
            random_state=random_state,
            class_weight=kwargs.get("class_weight", "balanced"),
        )
    elif clf_name == "random_forest":
        from sklearn.ensemble import RandomForestClassifier

        model = RandomForestClassifier(
            n_estimators=kwargs.get("n_estimators", 200),
            max_depth=kwargs.get("max_depth"),
            random_state=random_state,
            n_jobs=kwargs.get("n_jobs", -1),
            class_weight=kwargs.get("class_weight", "balanced_subsample"),
        )
    elif clf_name == "svc":
        from sklearn.svm import SVC

        model = SVC(
            C=kwargs.get("C", 1.0),
            kernel=kwargs.get("kernel", "rbf"),
            probability=True,
            random_state=random_state,
            class_weight=kwargs.get("class_weight", "balanced"),
        )
    elif clf_name == "knn":
        from sklearn.neighbors import KNeighborsClassifier

        model = KNeighborsClassifier(
            n_neighbors=kwargs.get("n_neighbors", 5),
            weights=kwargs.get("weights", "distance"),
            n_jobs=kwargs.get("n_jobs"),
        )
    else:
        raise ValueError(f"Unknown classifier '{classifier}'.")

    model.fit(X_train, y_train)
    proba = model.predict_proba(X_test)[:, 1]
    y_pred = (proba >= 0.5).astype(int)

    proba_train = model.predict_proba(X_train)[:, 1]
    y_pred_train = (proba_train >= 0.5).astype(int)

    metrics = {
        "test": {
            "accuracy": float(accuracy_score(y_test, y_pred)),
            "roc_auc": float(roc_auc_score(y_test, proba)),
            "confusion_matrix": confusion_matrix(y_test, y_pred).tolist(),
        },
        "train": {
            "accuracy": float(accuracy_score(y_train, y_pred_train)),
            "roc_auc": float(roc_auc_score(y_train, proba_train)),
            "confusion_matrix": confusion_matrix(y_train, y_pred_train).tolist(),
        },
    }

    preds_df = pd.DataFrame(
        {
            "true_label": le.inverse_transform(y_test),
            "pred_label": le.inverse_transform(y_pred),
            "proba": proba,
            "split": "test",
            "sample_label": lbl_test,
            "row_index": idx_test,
        }
    )

    # Also record train predictions for diagnostics
    train_df = pd.DataFrame(
        {
            "true_label": le.inverse_transform(y_train),
            "pred_label": le.inverse_transform(y_pred_train),
            "proba": proba_train,
            "split": "train",
            "sample_label": lbl_train,
            "row_index": idx_train,
        }
    )

    all_preds = pd.concat([train_df, preds_df], ignore_index=True)

    return {
        "model": model,
        "metrics": metrics,
        "predictions": all_preds,
        "label_mapping": {
            cls: int(code)
            for cls, code in zip(classes, le.transform(classes), strict=False)
        },
    }


def plot_confusion_matrices(predictions: pd.DataFrame, *, cmap: str = "Blues"):
    """
    Plot train/test confusion matrices from the predictions DataFrame returned by classify_read_features_binary.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.metrics import confusion_matrix

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for ax, split in zip(axes, ["train", "test"], strict=False):
        df_split = predictions[predictions["split"] == split]
        labels_true = df_split["true_label"]
        labels_pred = df_split["pred_label"]
        uniq = sorted(labels_true.unique())
        cm = confusion_matrix(labels_true, labels_pred, labels=uniq)
        sns.heatmap(
            cm,
            annot=True,
            fmt="d",
            cmap=cmap,
            xticklabels=uniq,
            yticklabels=uniq,
            ax=ax,
        )
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_title(f"{split.capitalize()} confusion")
    fig.tight_layout()
    return fig


def _resolve_prediction_row_selection(
    matrix: np.ndarray,
    preds: pd.DataFrame,
    *,
    owner: str,
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Align a feature/signal matrix to a predictions split.

    Supports two calling conventions:
    1) matrix is the full dataset and ``preds.row_index`` maps split rows.
    2) matrix is already split-aligned and has the same row count as ``preds``.
    """
    matrix_full = np.asarray(matrix)
    if "row_index" not in preds.columns:
        if len(preds) != matrix_full.shape[0]:
            raise ValueError(
                f"{owner}: predictions and matrix must have matching rows for the chosen split "
                "or predictions must include a row_index column."
            )
        return matrix_full, None

    row_idx = preds["row_index"].to_numpy(dtype=int)
    if row_idx.size == 0:
        return matrix_full[:0], row_idx

    in_bounds = (row_idx >= 0) & (row_idx < matrix_full.shape[0])
    if in_bounds.all():
        return matrix_full[row_idx], row_idx

    # If matrix is already split-aligned, keep it as-is and ignore out-of-range
    # row indices (which refer to original global row IDs).
    if len(preds) == matrix_full.shape[0]:
        return matrix_full, None

    raise ValueError(
        f"{owner}: row_index contains out-of-bounds values for the provided matrix "
        f"(max row_index={int(row_idx.max())}, matrix_rows={matrix_full.shape[0]})."
    )


def plot_classification_profiles(
    data_matrix: np.ndarray,
    predictions: pd.DataFrame,
    *,
    val_matrix: np.ndarray | None = None,
    split: str = "test",
    group_by: str = "pred_label",
    positive_label: str | None = None,
    motif_index: int = 0,
    motif_count: int | None = None,
    metadata: Sequence[dict[str, Any]] | None = None,
    motif_labels: Sequence[str] | None = None,
    motif_colors: Sequence[str] | None = None,
    plot_all_motifs: bool = False,
    motif_profile_mode: str = "single_axis",
    color_points_by: str = "cluster",
    view_window_size: int | None = None,
    window_size: int | None = None,
    show_valid_profile: bool = False,
    smoothing: str | None = None,
    smooth_win: int = 21,
    smooth_sigma: float = 6.0,
    show_unsmoothed_overlay: bool = False,
    cmap_name: str = "viridis",
    invert_y: bool = True,
    point_size: float = 1.0,
    point_alpha: float = 0.01,
):
    """
    Visualize per-read modification calls colored by prediction-derived groups.
    For binary models, ``group_by='confusion'`` renders TP/TN/FP/FN groups.

    Args:
        data_matrix: per-read windows (n_reads, window_len)
        predictions: DataFrame from classify_read_features_binary
        val_matrix: optional valid-site matrix to overlay
        split: which split to visualize ('train' or 'test')
    """

    preds = predictions[predictions["split"] == split].reset_index(drop=True)
    X_full = np.asarray(data_matrix)
    X, row_idx = _resolve_prediction_row_selection(
        X_full,
        preds,
        owner="plot_classification_profiles",
    )
    V_full = np.asarray(val_matrix) if val_matrix is not None else None
    M_full = list(metadata) if metadata is not None else None

    if V_full is not None:
        row_idx_max = (
            int(row_idx.max()) if row_idx is not None and row_idx.size > 0 else -1
        )
        if row_idx is not None and row_idx.size > 0 and row_idx_max < V_full.shape[0]:
            V = V_full[row_idx]
        elif V_full.shape[0] == len(preds):
            V = V_full
        else:
            raise ValueError(
                "plot_classification_profiles: val_matrix could not be aligned to predictions. "
                "Pass either full val_matrix (with row_index) or split-aligned val_matrix."
            )
    else:
        V = None

    if M_full is not None:
        row_idx_max = (
            int(row_idx.max()) if row_idx is not None and row_idx.size > 0 else -1
        )
        if row_idx is not None and row_idx.size > 0 and row_idx_max < len(M_full):
            M = [M_full[i] for i in row_idx]
        elif len(M_full) == len(preds):
            M = M_full
        else:
            raise ValueError(
                "plot_classification_profiles: metadata could not be aligned to predictions. "
                "Pass either full metadata (with row_index) or split-aligned metadata."
            )
    else:
        M = None

    if group_by not in {"pred_label", "true_label", "confusion"}:
        raise ValueError(
            "group_by must be one of {'pred_label', 'true_label', 'confusion'}."
        )

    if group_by == "confusion":
        class_values = sorted(pd.unique(preds["true_label"]))
        if len(class_values) != 2:
            raise ValueError("group_by='confusion' requires exactly two classes.")
        positive = positive_label if positive_label is not None else class_values[-1]
        true_vals = preds["true_label"].to_numpy()
        pred_vals = preds["pred_label"].to_numpy()
        confusion_labels = np.empty(len(preds), dtype=object)
        confusion_labels[(true_vals != positive) & (pred_vals != positive)] = "TN"
        confusion_labels[(true_vals != positive) & (pred_vals == positive)] = "FP"
        confusion_labels[(true_vals == positive) & (pred_vals != positive)] = "FN"
        confusion_labels[(true_vals == positive) & (pred_vals == positive)] = "TP"
        labels = confusion_labels
    else:
        labels = preds[group_by].to_numpy()

    return plot_cluster_profiles(
        X,
        labels,
        val_matrix=V,
        view_window_size=view_window_size,
        window_size=window_size,
        motif_index=motif_index,
        motif_count=motif_count,
        metadata=M,
        motif_labels=motif_labels,
        motif_colors=motif_colors,
        plot_all_motifs=plot_all_motifs,
        motif_profile_mode=motif_profile_mode,
        color_points_by=color_points_by,
        show_valid_profile=show_valid_profile,
        smoothing=smoothing,
        smooth_win=smooth_win,
        smooth_sigma=smooth_sigma,
        show_unsmoothed_overlay=show_unsmoothed_overlay,
        cmap_name=cmap_name,
        invert_y=invert_y,
        point_size=point_size,
        point_alpha=point_alpha,
    )


def plot_region_classification_profiles(
    pileup_matrix: np.ndarray,
    predictions: pd.DataFrame,
    *,
    split: str = "test",
    group_by: str = "pred_label",
    positive_label: str | None = None,
    window_size: int | None = None,
    motif_index: int | None = 0,
    motif_count: int | None = None,
    plot_all_motifs: bool = False,
    motif_labels: Sequence[str] | None = None,
    motif_colors: Sequence[str] | None = None,
    motif_profile_mode: str = "single_axis",
    color_points_by: str = "cluster",
    point_size: float = 1.0,
    point_alpha: float = 0.01,
    smoothing: str | None = None,
    smooth_win: int = 21,
    smooth_sigma: float = 6.0,
    show_unsmoothed_overlay: bool = False,
    cmap_name: str = "viridis",
):
    """
    Region-level counterpart to plot_classification_profiles.
    Uses classifier predictions to group regions by predicted label, true label,
    or confusion-class (TP/TN/FP/FN for binary tasks) and renders cluster-style profiles.
    """

    preds = predictions[predictions["split"] == split].reset_index(drop=True)
    X_full = np.asarray(pileup_matrix)
    X, _ = _resolve_prediction_row_selection(
        X_full,
        preds,
        owner="plot_region_classification_profiles",
    )

    if group_by not in {"pred_label", "true_label", "confusion"}:
        raise ValueError(
            "group_by must be one of {'pred_label', 'true_label', 'confusion'}."
        )

    if group_by == "confusion":
        class_values = sorted(pd.unique(preds["true_label"]))
        if len(class_values) != 2:
            raise ValueError("group_by='confusion' requires exactly two classes.")
        positive = positive_label if positive_label is not None else class_values[-1]
        true_vals = preds["true_label"].to_numpy()
        pred_vals = preds["pred_label"].to_numpy()
        labels = np.empty(len(preds), dtype=object)
        labels[(true_vals != positive) & (pred_vals != positive)] = "TN"
        labels[(true_vals != positive) & (pred_vals == positive)] = "FP"
        labels[(true_vals == positive) & (pred_vals != positive)] = "FN"
        labels[(true_vals == positive) & (pred_vals == positive)] = "TP"
    else:
        labels = preds[group_by].to_numpy()

    return plot_region_cluster_profiles(
        X,
        labels,
        window_size=window_size,
        motif_index=motif_index,
        motif_count=motif_count,
        plot_all_motifs=plot_all_motifs,
        motif_labels=motif_labels,
        motif_colors=motif_colors,
        motif_profile_mode=motif_profile_mode,
        color_points_by=color_points_by,
        point_size=point_size,
        point_alpha=point_alpha,
        smoothing=smoothing,
        smooth_win=smooth_win,
        smooth_sigma=smooth_sigma,
        show_unsmoothed_overlay=show_unsmoothed_overlay,
        cmap_name=cmap_name,
    )


def plot_multisite_read_raster(
    read_windows: ReadWindowExtractionResult,
    *,
    n_windows: int = 2,
    min_separation_bp: int = 5000,
    distance_mode: str = "center",  # or "bounds"
    window_size: int | None = None,
    motif_index: int = 0,
    motif_count: int | None = None,
    plot_all_motifs: bool = False,
    motif_labels: Sequence[str] | None = None,
    window_centers_bp: Sequence[int] | None = None,
    window_widths_bp: int | Sequence[int] | None = None,
    window_match_tolerance_bp: int | None = None,
    symmetric_side_windows: int | None = None,
    symmetric_max_offset_bp: int | None = None,
    enforce_full_window_span: bool = True,
    # Legacy aliases
    window_center_offsets: Sequence[int] | None = None,
    center_tolerance_bp: int | None = None,
    min_read_length_bp: int | None = None,
    smoothing: str | None = "gaussian",  # None, "boxcar", "gaussian"
    smooth_win: int = 21,
    smooth_sigma_bp: float = 6.0,
    max_rows: int | None = 500,
    downsample_method: str = "bin_mean",  # "bin_mean" | "uniform"
    cmap: str = "magma",
    vmin: float | None = None,
    vmax: float | None = None,
    render_mode: str = "scatter",  # "scatter" | "heatmap"
    render_values: str = "auto",  # "auto" | "raw" | "smoothed"
    axis_orientation: str | None = None,  # "position_x" | "position_y"
    scatter_size: float = 4.0,
    scatter_alpha: float = 0.7,
    scatter_color_values: str = "ml_score",  # "auto" | "ml_score" | "raw" | "smoothed"
    ml_score_thresholds: Sequence[float] | None = None,
    motif_colors: Sequence[str] | None = None,
    x_axis_mode: str | None = None,  # "relative_to_primary" | "centered"
    sort_by: str = "mod_fraction",  # "mod_fraction" | "cluster" | "window_center" | "read_name" | "region_start" | "read_length" | "none"
    sort_window_index: int | None = None,
    sort_descending: bool | None = None,
    sort_labels: Sequence[int] | None = None,
    beds: Sequence[str | Path] | None = None,
    rotate: bool = True,
):
    """
    Plot multi-window read rasters, optionally across multiple motifs.

    Args:
        read_windows: ReadWindowExtractionResult from extract_read_windows / build_multimotif_read_windows
        n_windows: number of windows per read to plot when explicit centers are not provided
        min_separation_bp: minimum separation between site centers (distance_mode="center") or read starts (bounds)
        distance_mode: "center" uses region_start/end center; "bounds" uses read_start/read_end
        window_size: half-window in bp; if provided, per-motif slice width is ``2 * window_size``.
        motif_index: reference motif slice index for single-motif plotting and ordering
        motif_count: total motif slices concatenated in data_matrix; inferred when possible
        plot_all_motifs: when True, render all motif slices in a motif x window panel grid
        motif_labels: optional motif names aligned to motif slices
        window_centers_bp: explicit desired window centers in primary-region coordinates
            (for example ``[-2500, 0, 2500]``). When provided, this replaces ``n_windows``.
        window_widths_bp: per-window display widths in bp (common width when scalar). Defaults
            to the extracted span for each panel.
        window_match_tolerance_bp: matching tolerance for explicit centers. Defaults to half
            the extracted span.
        symmetric_side_windows: when set with ``symmetric_max_offset_bp``, generate centers as
            evenly spaced offsets from ``-symmetric_max_offset_bp`` to ``+symmetric_max_offset_bp``
            with ``2 * symmetric_side_windows + 1`` windows total.
        symmetric_max_offset_bp: max absolute offset used for symmetric center generation.
        enforce_full_window_span: when True (default), drop read-sets that cannot span the
            full displayed windows based on read_start/read_end (or read_length when available).
        window_center_offsets: legacy alias for ``window_centers_bp``.
        center_tolerance_bp: legacy alias for ``window_match_tolerance_bp``.
        min_read_length_bp: optional minimum read length filter. Must be at least the required
            span implied by requested centers/separation plus plotted window width.
        smoothing: None, "boxcar", or "gaussian"
        smooth_win: smoothing window length
        smooth_sigma_bp: sigma for gaussian smoothing
        max_rows: cap rows; when exceeded, downsample according to ``downsample_method``.
            Pass None to disable downsampling.
        downsample_method: "bin_mean" averages adjacent sorted read rows; "uniform" takes
            evenly spaced sorted rows.
        cmap: matplotlib colormap
        vmin/vmax: color scale limits; None auto-scales
        render_mode: "scatter" or "heatmap"
        render_values: value source used for plotting.
            - "auto": scatter uses smoothed values when smoothing is enabled; heatmap uses raw per-read values
            - "raw": always plot raw per-read values
            - "smoothed": always plot smoothed values
        axis_orientation: controls which axis encodes position.
            - "position_x": x=position, y=reads
            - "position_y": x=reads, y=position
            When None, preserves legacy behavior for heatmaps based on ``rotate`` and uses ``position_x`` for scatter.
        scatter_size: marker size for scatter mode
        scatter_alpha: marker alpha for scatter mode
        scatter_color_values: scatter color source.
            - "ml_score": normalized ML score from mod_vector values (clipped to 0..1)
            - "raw": raw per-read values
            - "smoothed": smoothed per-read values
            - "auto": ml_score (default behavior)
        ml_score_thresholds: optional per-motif ML score thresholds for scatter mode. When provided,
            only points with ML score >= threshold are plotted in motif-specific colors.
        motif_colors: optional colors aligned to motifs for thresholded scatter display
        x_axis_mode: coordinate mode. ``relative_to_primary`` preserves primary-region offsets;
            ``centered`` keeps each panel centered at zero for legacy compatibility.
        sort_by: ordering strategy for paired reads
        sort_window_index: optional window index used by window-aware sorting modes
        sort_descending: descending/ascending order toggle for supported sort modes.
            If None, defaults are mode-aware (descending for mod_fraction, ascending otherwise).
        sort_labels: cluster labels to sort by when sort_by="cluster"
        beds: optional list of BED/GTF/GFF paths; only reads intersecting all beds are retained
        rotate: subplot layout orientation (legacy compatibility); does not force data-axis orientation
    """

    import math

    import matplotlib.pyplot as plt

    if render_mode not in {"scatter", "heatmap"}:
        raise ValueError("render_mode must be 'scatter' or 'heatmap'.")
    if render_values not in {"auto", "raw", "smoothed"}:
        raise ValueError("render_values must be 'auto', 'raw', or 'smoothed'.")
    if axis_orientation is not None and axis_orientation not in {
        "position_x",
        "position_y",
    }:
        raise ValueError(
            "axis_orientation must be None, 'position_x', or 'position_y'."
        )
    if scatter_color_values not in {"auto", "ml_score", "raw", "smoothed"}:
        raise ValueError(
            "scatter_color_values must be 'auto', 'ml_score', 'raw', or 'smoothed'."
        )
    if x_axis_mode is not None and x_axis_mode not in {
        "relative_to_primary",
        "centered",
    }:
        raise ValueError(
            "x_axis_mode must be None, 'relative_to_primary', or 'centered'."
        )
    effective_x_axis_mode = x_axis_mode or "relative_to_primary"
    if distance_mode not in {"center", "bounds"}:
        raise ValueError("distance_mode must be 'center' or 'bounds'.")
    if n_windows < 1:
        raise ValueError("n_windows must be >= 1.")
    if min_separation_bp < 0:
        raise ValueError("min_separation_bp must be >= 0.")
    if downsample_method not in {"bin_mean", "uniform"}:
        raise ValueError("downsample_method must be 'bin_mean' or 'uniform'.")
    if center_tolerance_bp is not None and center_tolerance_bp < 0:
        raise ValueError("center_tolerance_bp must be >= 0 when provided.")
    if window_match_tolerance_bp is not None and window_match_tolerance_bp < 0:
        raise ValueError("window_match_tolerance_bp must be >= 0 when provided.")
    if max_rows is not None and max_rows <= 0:
        raise ValueError("max_rows must be > 0 when provided.")
    if sort_window_index is not None and sort_window_index < 0:
        raise ValueError("sort_window_index must be >= 0 when provided.")

    def _kernel(kind: str, win: int, sigma: float):
        if win % 2 == 0:
            win += 1
        if kind == "gaussian":
            r = win // 2
            x = np.arange(-r, r + 1)
            k = np.exp(-0.5 * (x / sigma) ** 2)
            return k / k.sum()
        else:
            k = np.ones(win, dtype=float)
            return k / k.sum()

    def _smooth(M: np.ndarray) -> np.ndarray:
        if smoothing is None:
            return M
        k = _kernel(
            "gaussian" if smoothing == "gaussian" else "boxcar",
            smooth_win,
            smooth_sigma_bp,
        )
        r = len(k) // 2
        out = np.empty_like(M, dtype=float)
        for i in range(M.shape[0]):
            row = M[i]
            pad_left = row[1 : r + 1][::-1] if r > 0 else row[:0]
            pad_right = row[-r - 1 : -1][::-1] if r > 0 else row[:0]
            padded = np.concatenate([pad_left, row, pad_right])
            out[i] = np.convolve(padded, k, mode="same")[r : r + row.size]
        return out

    # Optional bed filtering
    bed_filters = []
    if beds:
        try:
            import pyranges as pr
        except Exception as exc:  # pragma: no cover - optional dependency
            raise ImportError("pyranges is required for bed filtering.") from exc
        for bed in beds:
            bed_filters.append(pr.read_bed(bed))

    meta = list(read_windows.metadata)
    filtered_original_indices = np.arange(len(meta), dtype=int)
    X_full = np.asarray(read_windows.data_matrix)
    if X_full.ndim != 2:
        raise ValueError("read_windows.data_matrix must be a 2D array.")
    width = X_full.shape[1]
    n_motifs, slice_width = _resolve_motif_slices(
        width,
        motif_count=motif_count,
        motif_labels=motif_labels,
        window_size=window_size,
    )
    if n_motifs <= 0:
        n_motifs = 1
    if slice_width <= 0:
        raise ValueError(
            "Unable to infer motif slice width from read_windows.data_matrix."
        )
    if width % slice_width != 0:
        raise ValueError(
            "read_windows.data_matrix width is not divisible into consistent motif slices. "
            "Pass an explicit motif_count or compatible window_size."
        )
    if motif_index < 0 or motif_index >= n_motifs:
        raise ValueError(
            f"motif_index {motif_index} out of range for {n_motifs} motif slices."
        )
    motif_ids = (
        list(range(n_motifs)) if plot_all_motifs and n_motifs > 1 else [motif_index]
    )
    if motif_labels is None:
        motif_labels_resolved = [f"motif_{i}" for i in range(n_motifs)]
    else:
        motif_labels_resolved = list(motif_labels)
        if len(motif_labels_resolved) < n_motifs:
            motif_labels_resolved.extend(
                [f"motif_{i}" for i in range(len(motif_labels_resolved), n_motifs)]
            )
    X_by_motif = [X_full[:, m * slice_width : (m + 1) * slice_width] for m in motif_ids]
    ML_by_motif = [np.clip(matrix.astype(float), 0.0, 1.0) for matrix in X_by_motif]

    if (
        window_center_offsets is not None
        and window_centers_bp is not None
        and list(window_center_offsets) != list(window_centers_bp)
    ):
        raise ValueError(
            "Pass either window_centers_bp or window_center_offsets (legacy alias), not both with different values."
        )
    resolved_centers = (
        window_centers_bp if window_centers_bp is not None else window_center_offsets
    )
    if (
        center_tolerance_bp is not None
        and window_match_tolerance_bp is not None
        and int(center_tolerance_bp) != int(window_match_tolerance_bp)
    ):
        raise ValueError(
            "Pass either window_match_tolerance_bp or center_tolerance_bp (legacy alias), not both with different values."
        )
    match_tolerance_bp = (
        int(window_match_tolerance_bp)
        if window_match_tolerance_bp is not None
        else (int(center_tolerance_bp) if center_tolerance_bp is not None else None)
    )

    if symmetric_side_windows is not None or symmetric_max_offset_bp is not None:
        if resolved_centers is not None:
            raise ValueError(
                "Do not combine symmetric center spec with explicit window_centers_bp."
            )
        if symmetric_side_windows is None or symmetric_max_offset_bp is None:
            raise ValueError(
                "symmetric_side_windows and symmetric_max_offset_bp must be provided together."
            )
        if int(symmetric_side_windows) < 0:
            raise ValueError("symmetric_side_windows must be >= 0.")
        if int(symmetric_max_offset_bp) < 0:
            raise ValueError("symmetric_max_offset_bp must be >= 0.")
        n_side = int(symmetric_side_windows)
        max_offset = float(symmetric_max_offset_bp)
        if n_side == 0:
            resolved_centers = [0]
        else:
            resolved_centers = (
                np.linspace(-max_offset, max_offset, (2 * n_side) + 1)
                .round()
                .astype(int)
                .tolist()
            )

    if resolved_centers is None:
        n_windows_effective = int(n_windows)
        target_offsets = None
    else:
        offsets = [int(offset) for offset in resolved_centers]
        if len(offsets) == 0:
            raise ValueError(
                "window_centers_bp must contain at least one center when provided."
            )
        n_windows_effective = len(offsets)
        target_offsets = np.asarray(offsets, dtype=float)

    if target_offsets is None:
        required_read_length_bp = int(
            max(
                slice_width,
                ((n_windows_effective - 1) * min_separation_bp) + slice_width,
            )
        )
    else:
        required_read_length_bp = int(
            max(
                slice_width,
                np.ceil(
                    float(target_offsets.max() - target_offsets.min())
                    + float(slice_width)
                ),
            )
        )

    if min_read_length_bp is not None and min_read_length_bp < required_read_length_bp:
        raise ValueError(
            "min_read_length_bp must be >= required span implied by requested windows "
            f"({required_read_length_bp} bp)."
        )
    effective_min_read_length = (
        required_read_length_bp
        if target_offsets is not None and min_read_length_bp is None
        else min_read_length_bp
    )

    if window_widths_bp is None:
        default_width_bp = int(
            slice_width if window_size is None else max(1, 2 * int(window_size))
        )
        panel_widths_bp = [default_width_bp for _ in range(n_windows_effective)]
    elif isinstance(window_widths_bp, (int, np.integer)):
        panel_widths_bp = [int(window_widths_bp) for _ in range(n_windows_effective)]
    else:
        panel_widths_bp = [int(width) for width in window_widths_bp]
        if len(panel_widths_bp) != n_windows_effective:
            raise ValueError(
                "window_widths_bp must be a scalar or a sequence with one entry per window center."
            )
    if any(width <= 0 for width in panel_widths_bp):
        raise ValueError("window_widths_bp values must be > 0.")
    if any(width > slice_width for width in panel_widths_bp):
        raise ValueError(
            "window_widths_bp values cannot exceed extracted window span; increase extract/parse window_size."
        )
    panel_half_widths = [0.5 * float(width) for width in panel_widths_bp]

    def _apply_index_filter(index_vector: np.ndarray) -> None:
        nonlocal X_by_motif, ML_by_motif, meta, filtered_original_indices
        X_by_motif = [matrix[index_vector] for matrix in X_by_motif]
        ML_by_motif = [matrix[index_vector] for matrix in ML_by_motif]
        meta = [meta[i] for i in index_vector]
        filtered_original_indices = filtered_original_indices[index_vector]

    # Filter by beds if provided
    keep_idx = np.arange(len(meta))
    if bed_filters:
        interval_rows = []
        for idx, m in enumerate(meta):
            try:
                chrom = m.get("chromosome")
                start_bp = int(m.get("region_start", -1))
                end_bp = int(m.get("region_end", -1))
            except Exception:
                continue
            if chrom is None or start_bp < 0 or end_bp <= start_bp:
                continue
            interval_rows.append(
                {
                    "read_idx": idx,
                    "Chromosome": chrom,
                    "Start": start_bp,
                    "End": end_bp,
                }
            )

        if interval_rows:
            intervals = pd.DataFrame(interval_rows)
            candidate_indices = intervals["read_idx"].to_numpy(dtype=int)
            for bed_filter in bed_filters:
                if candidate_indices.size == 0:
                    break
                candidate_df = intervals.loc[
                    intervals["read_idx"].isin(candidate_indices),
                    ["read_idx", "Chromosome", "Start", "End"],
                ]
                joined = pr.PyRanges(candidate_df).join(bed_filter)
                if joined.df.empty:
                    candidate_indices = np.array([], dtype=int)
                    break
                candidate_indices = (
                    joined.df["read_idx"].drop_duplicates().to_numpy(dtype=int)
                )

            keep_idx = np.sort(candidate_indices)
        else:
            keep_idx = np.array([], dtype=int)
        _apply_index_filter(keep_idx)

    if effective_min_read_length is not None:
        keep_by_length: list[int] = []
        for idx, m in enumerate(meta):
            read_length_val = m.get("read_length")
            parsed_length: int | None = None
            try:
                if read_length_val is not None:
                    parsed_length = int(read_length_val)
            except Exception:
                parsed_length = None
            if (parsed_length is None or parsed_length <= 0) and {
                "read_start",
                "read_end",
            }.issubset(m.keys()):
                try:
                    parsed_length = int(m["read_end"]) - int(m["read_start"])
                except Exception:
                    parsed_length = None
            if parsed_length is not None and parsed_length >= int(
                effective_min_read_length
            ):
                keep_by_length.append(idx)
        _apply_index_filter(np.asarray(keep_by_length, dtype=int))

    span_check_warning_emitted = False
    dropped_for_window_span = 0

    def _safe_int(value: Any) -> int | None:
        try:
            if value is None:
                return None
            return int(value)
        except Exception:
            return None

    def _read_bounds(meta_row: dict[str, Any]) -> tuple[int, int] | None:
        start = _safe_int(meta_row.get("read_start"))
        end = _safe_int(meta_row.get("read_end"))
        length = _safe_int(meta_row.get("read_length"))
        if start is not None and end is not None and end > start:
            return (start, end)
        if start is not None and length is not None and length > 0:
            return (start, start + length)
        if end is not None and length is not None and length > 0:
            return (end - length, end)
        return None

    def _read_length(meta_row: dict[str, Any]) -> int | None:
        length = _safe_int(meta_row.get("read_length"))
        if length is not None and length > 0:
            return length
        bounds = _read_bounds(meta_row)
        if bounds is None:
            return None
        return bounds[1] - bounds[0]

    def _window_is_covered(
        meta_row: dict[str, Any], center_bp: float, half_width_bp: float
    ) -> bool | None:
        bounds = _read_bounds(meta_row)
        if bounds is None:
            return None
        left = float(center_bp) - float(half_width_bp)
        right = float(center_bp) + float(half_width_bp)
        return left >= float(bounds[0]) and right <= float(bounds[1])

    def _selection_spans_all_windows(
        selected: Sequence[tuple[int, float, int]],
    ) -> bool:
        nonlocal span_check_warning_emitted
        if not enforce_full_window_span:
            return True

        unresolved_present = False
        left_edges: list[float] = []
        right_edges: list[float] = []
        candidate_lengths: list[int] = []
        for row_idx, center_bp, window_pos in selected:
            half_width_bp = float(panel_half_widths[window_pos])
            covered = _window_is_covered(meta[row_idx], center_bp, half_width_bp)
            if covered is False:
                return False
            if covered is None:
                unresolved_present = True
            left_edges.append(float(center_bp) - half_width_bp)
            right_edges.append(float(center_bp) + half_width_bp)
            length_value = _read_length(meta[row_idx])
            if length_value is not None:
                candidate_lengths.append(int(length_value))

        if unresolved_present:
            if candidate_lengths:
                required_span = max(right_edges) - min(left_edges)
                if max(candidate_lengths) < required_span:
                    return False
            elif not span_check_warning_emitted:
                warnings.warn(
                    "Unable to fully verify window-span coverage for some reads because "
                    "read_start/read_end/read_length metadata are missing. "
                    "Re-extract windows with current APIs for strict filtering.",
                    RuntimeWarning,
                    stacklevel=2,
                )
                span_check_warning_emitted = True
        return True

    # Group by read key
    groups = defaultdict(list)
    for i, m in zip(range(len(meta)), meta, strict=False):
        key = (
            m.get("read_name"),
            m.get("chromosome"),
            m.get("region_strand"),
        )
        if distance_mode == "center":
            center = (int(m.get("region_start", 0)) + int(m.get("region_end", 0))) // 2
        else:
            center = int(m.get("read_start", 0))
        groups[key].append((i, center))

    pair_indices_per_window = [[] for _ in range(n_windows_effective)]
    panel_centers = [[] for _ in range(n_windows_effective)]
    for _, items in groups.items():
        if len(items) < n_windows_effective:
            continue
        items = sorted(items, key=lambda x: x[1])
        centers_arr = np.asarray([entry[1] for entry in items], dtype=float)
        idx_arr = np.asarray([entry[0] for entry in items], dtype=int)

        if target_offsets is None:
            # sliding window of n_windows with minimum separation
            for j in range(len(items) - n_windows_effective + 1):
                centers = [items[j + k][1] for k in range(n_windows_effective)]
                if any(
                    (centers[k + 1] - centers[k]) < min_separation_bp
                    for k in range(n_windows_effective - 1)
                ):
                    continue
                idxs = [items[j + k][0] for k in range(n_windows_effective)]
                candidate = [
                    (int(idxs[k]), float(centers[k]), int(k))
                    for k in range(n_windows_effective)
                ]
                if not _selection_spans_all_windows(candidate):
                    dropped_for_window_span += 1
                    continue
                for k, idx in enumerate(idxs):
                    pair_indices_per_window[k].append(int(idx))
                    panel_centers[k].append(float(centers[k]))
            continue

        tolerance = float(
            match_tolerance_bp
            if match_tolerance_bp is not None
            else max(1, slice_width // 2)
        )
        best_candidate: tuple[float, list[int]] | None = None
        for anchor_pos in range(len(items)):
            anchor_center = centers_arr[anchor_pos]
            targets = anchor_center + target_offsets
            used_positions: set[int] = set()
            selected_positions: list[int] = []
            total_error = 0.0
            for target in targets:
                distances = np.abs(centers_arr - target)
                order = np.argsort(distances)
                chosen = None
                for candidate_pos in order:
                    pos = int(candidate_pos)
                    if pos in used_positions:
                        continue
                    if float(distances[pos]) <= tolerance:
                        candidate_idx = int(idx_arr[pos])
                        window_half = float(panel_half_widths[len(selected_positions)])
                        covered = _window_is_covered(
                            meta[candidate_idx], float(centers_arr[pos]), window_half
                        )
                        if covered is False:
                            continue
                        chosen = pos
                        break
                if chosen is None:
                    selected_positions = []
                    break
                used_positions.add(chosen)
                selected_positions.append(chosen)
                total_error += float(distances[chosen])

            if not selected_positions:
                continue
            if best_candidate is None or total_error < best_candidate[0]:
                best_candidate = (total_error, selected_positions)

        if best_candidate is None:
            continue
        _, selected_positions = best_candidate
        candidate = [
            (int(idx_arr[pos]), float(centers_arr[pos]), int(k))
            for k, pos in enumerate(selected_positions)
        ]
        if not _selection_spans_all_windows(candidate):
            dropped_for_window_span += 1
            continue
        for k, pos in enumerate(selected_positions):
            pair_indices_per_window[k].append(int(idx_arr[pos]))
            panel_centers[k].append(float(centers_arr[pos]))

    if not all(len(p) > 0 for p in pair_indices_per_window):
        if enforce_full_window_span and dropped_for_window_span > 0:
            raise ValueError(
                "No read sets found after enforcing full-window span coverage; "
                "increase read length, reduce window widths/offsets, or relax constraints."
            )
        raise ValueError("No read sets found meeting separation criteria.")

    window_indices = [
        np.asarray(values, dtype=int) for values in pair_indices_per_window
    ]
    panel_centers = [
        np.asarray(center_values, dtype=float) for center_values in panel_centers
    ]
    panels_by_motif = []
    panels_smooth_by_motif = []
    panels_ml_by_motif = []
    for motif_matrix in X_by_motif:
        motif_panels = [motif_matrix[idx] for idx in window_indices]
        panels_by_motif.append(motif_panels)
        panels_smooth_by_motif.append([_smooth(panel) for panel in motif_panels])
    for motif_ml_matrix in ML_by_motif:
        panels_ml_by_motif.append([motif_ml_matrix[idx] for idx in window_indices])

    if render_values == "raw":
        panels_render_by_motif = panels_by_motif
        render_value_mode = "raw"
    elif render_values == "smoothed":
        panels_render_by_motif = panels_smooth_by_motif
        render_value_mode = "smoothed"
    else:
        if render_mode == "heatmap" or smoothing is None:
            panels_render_by_motif = panels_by_motif
            render_value_mode = "raw"
        else:
            panels_render_by_motif = panels_smooth_by_motif
            render_value_mode = "smoothed"

    # Sorting uses the reference motif to keep all motif panels aligned.
    try:
        ref_motif_pos = motif_ids.index(motif_index)
    except ValueError:
        ref_motif_pos = 0
    reference_panels = panels_render_by_motif[ref_motif_pos]

    if sort_window_index is None:
        sort_window_pos = 0
    else:
        sort_window_pos = int(sort_window_index)
        if sort_window_pos >= n_windows_effective:
            raise ValueError(
                f"sort_window_index {sort_window_index} out of range for {n_windows_effective} windows."
            )

    if sort_descending is None:
        effective_sort_descending = sort_by == "mod_fraction"
    else:
        effective_sort_descending = bool(sort_descending)

    # Sorting
    if sort_by == "mod_fraction":
        if sort_window_index is not None:
            score = reference_panels[sort_window_pos].mean(axis=1)
            order = np.argsort(score, kind="stable")
            if effective_sort_descending:
                order = order[::-1]
        else:
            means = [panel.mean(axis=1) for panel in reference_panels]
            if effective_sort_descending:
                order = np.lexsort(tuple(-m for m in reversed(means)))
            else:
                order = np.lexsort(tuple(m for m in reversed(means)))
    elif sort_by == "cluster":
        if sort_labels is None:
            raise ValueError("sort_labels must be provided when sort_by='cluster'.")
        labels_array = np.asarray(sort_labels)
        if labels_array.shape[0] == len(read_windows.metadata):
            labels_filtered = labels_array[filtered_original_indices]
        elif labels_array.shape[0] == len(meta):
            labels_filtered = labels_array
        else:
            raise ValueError(
                "sort_labels must match either the original metadata length "
                "or the filtered metadata length."
            )
        pair_labels = labels_filtered[window_indices[0]]
        order = np.argsort(pair_labels, kind="stable")
        if effective_sort_descending:
            order = order[::-1]
    elif sort_by == "window_center":
        center_values = panel_centers[sort_window_pos]
        order = np.argsort(center_values, kind="stable")
        if effective_sort_descending:
            order = order[::-1]
    elif sort_by == "read_name":
        values = np.asarray(
            [str(meta[idx].get("read_name", "")) for idx in window_indices[0]]
        )
        order = np.argsort(values, kind="stable")
        if effective_sort_descending:
            order = order[::-1]
    elif sort_by == "region_start":
        values = np.asarray(
            [int(meta[idx].get("region_start", 0) or 0) for idx in window_indices[0]],
            dtype=float,
        )
        order = np.argsort(values, kind="stable")
        if effective_sort_descending:
            order = order[::-1]
    elif sort_by == "read_length":

        def _resolved_read_length(read_meta: dict) -> int:
            length_value = read_meta.get("read_length")
            if length_value is not None:
                try:
                    return int(length_value)
                except Exception:
                    pass
            try:
                return int(read_meta.get("read_end", 0)) - int(
                    read_meta.get("read_start", 0)
                )
            except Exception:
                return 0

        values = np.asarray(
            [_resolved_read_length(meta[idx]) for idx in window_indices[0]], dtype=float
        )
        order = np.argsort(values, kind="stable")
        if effective_sort_descending:
            order = order[::-1]
    elif sort_by == "none":
        order = np.arange(reference_panels[0].shape[0], dtype=int)
    else:
        raise ValueError(
            "sort_by must be one of 'mod_fraction', 'cluster', 'window_center', "
            "'read_name', 'region_start', 'read_length', or 'none'."
        )

    panels_sorted_by_motif = [
        [panel[order] for panel in motif_panels]
        for motif_panels in panels_render_by_motif
    ]
    panels_raw_sorted_by_motif = [
        [panel[order] for panel in motif_panels] for motif_panels in panels_by_motif
    ]
    panels_smooth_sorted_by_motif = [
        [panel[order] for panel in motif_panels]
        for motif_panels in panels_smooth_by_motif
    ]
    panels_ml_sorted_by_motif: list[list[np.ndarray]] = [
        [panel[order] for panel in motif_panels] for motif_panels in panels_ml_by_motif
    ]
    expected_rows = panels_raw_sorted_by_motif[0][0].shape[0]
    for motif_panels in panels_raw_sorted_by_motif:
        for panel in motif_panels:
            if panel.shape[0] != expected_rows:
                raise RuntimeError(
                    "Inconsistent read ordering across multisite panels; panel row counts diverged."
                )
    centers_sorted = [center_values[order] for center_values in panel_centers]
    primary_centers = centers_sorted[0]
    if effective_x_axis_mode == "centered":
        center_offsets = [
            np.zeros_like(center_values, dtype=float)
            for center_values in centers_sorted
        ]
    else:
        center_offsets = [
            center_values - primary_centers for center_values in centers_sorted
        ]
    if render_mode == "heatmap" and effective_x_axis_mode == "relative_to_primary":
        for offsets in center_offsets:
            if offsets.size > 1 and not np.allclose(
                offsets, offsets[0], equal_nan=True
            ):
                raise ValueError(
                    "Heatmap rendering with x_axis_mode='relative_to_primary' requires constant "
                    "per-read offsets within each panel because imshow uses one rectangular extent. "
                    "Use render_mode='scatter', x_axis_mode='centered', or fixed/explicit centers "
                    "producing constant offsets."
                )

    # Downsample rows if too tall
    def bin_rows(M, k):
        n_bins = math.ceil(M.shape[0] / k)
        out = np.zeros((n_bins, M.shape[1]))
        for i in range(n_bins):
            s, e = i * k, min((i + 1) * k, M.shape[0])
            out[i] = M[s:e].mean(axis=0)
        return out

    def bin_vector(vector: np.ndarray, k: int) -> np.ndarray:
        n_bins = math.ceil(vector.shape[0] / k)
        out = np.zeros(n_bins, dtype=float)
        for i in range(n_bins):
            s, e = i * k, min((i + 1) * k, vector.shape[0])
            out[i] = vector[s:e].mean()
        return out

    P = panels_sorted_by_motif[0][0].shape[0]
    rows_before_downsample = int(P)
    downsampled = False
    downsample_factor = 1
    if max_rows is not None and max_rows < P:
        step = math.ceil(P / max_rows)
        downsample_factor = step
        if downsample_method == "bin_mean":
            panels_sorted_by_motif = [
                [bin_rows(panel, step) for panel in motif_panels]
                for motif_panels in panels_sorted_by_motif
            ]
            panels_raw_sorted_by_motif = [
                [bin_rows(panel, step) for panel in motif_panels]
                for motif_panels in panels_raw_sorted_by_motif
            ]
            panels_smooth_sorted_by_motif = [
                [bin_rows(panel, step) for panel in motif_panels]
                for motif_panels in panels_smooth_sorted_by_motif
            ]
            panels_ml_sorted_by_motif = [
                [bin_rows(panel, step) for panel in motif_panels]
                for motif_panels in panels_ml_sorted_by_motif
            ]
            center_offsets = [bin_vector(v, step) for v in center_offsets]
        else:
            keep = np.unique(np.round(np.linspace(0, P - 1, max_rows)).astype(int))
            panels_sorted_by_motif = [
                [panel[keep] for panel in motif_panels]
                for motif_panels in panels_sorted_by_motif
            ]
            panels_raw_sorted_by_motif = [
                [panel[keep] for panel in motif_panels]
                for motif_panels in panels_raw_sorted_by_motif
            ]
            panels_smooth_sorted_by_motif = [
                [panel[keep] for panel in motif_panels]
                for motif_panels in panels_smooth_sorted_by_motif
            ]
            panels_ml_sorted_by_motif = [
                [panel[keep] for panel in motif_panels]
                for motif_panels in panels_ml_sorted_by_motif
            ]
            center_offsets = [vector[keep] for vector in center_offsets]
        P = panels_sorted_by_motif[0][0].shape[0]
        downsampled = True

    if scatter_color_values == "auto":
        effective_scatter_color_values = "ml_score"
    else:
        effective_scatter_color_values = scatter_color_values

    threshold_map: dict[int, float] | None = None
    if render_mode == "scatter" and ml_score_thresholds is not None:
        thresholds_list = [float(value) for value in ml_score_thresholds]
        if len(thresholds_list) == len(motif_ids):
            threshold_map = {
                motif_id: float(thresholds_list[idx])
                for idx, motif_id in enumerate(motif_ids)
            }
        elif len(thresholds_list) == n_motifs:
            threshold_map = {
                idx: float(value) for idx, value in enumerate(thresholds_list)
            }
        elif len(thresholds_list) == 1:
            threshold_map = {
                motif_id: float(thresholds_list[0]) for motif_id in motif_ids
            }
        else:
            raise ValueError(
                "ml_score_thresholds must have length 1, motif_count, or motifs_plotted length."
            )

    motif_color_map: dict[int, str] = {}
    if threshold_map is not None:
        if motif_colors is not None:
            colors_list = [str(value) for value in motif_colors]
            if len(colors_list) == len(motif_ids):
                motif_color_map = {
                    motif_id: colors_list[idx] for idx, motif_id in enumerate(motif_ids)
                }
            elif len(colors_list) == n_motifs:
                motif_color_map = {
                    motif_id: colors_list[motif_id] for motif_id in motif_ids
                }
            elif len(colors_list) == 1:
                motif_color_map = {motif_id: colors_list[0] for motif_id in motif_ids}
            else:
                raise ValueError(
                    "motif_colors must have length 1, motif_count, or motifs_plotted length."
                )
        else:
            tab10 = plt.get_cmap("tab10")
            motif_color_map = {
                motif_id: tab10(idx % 10) for idx, motif_id in enumerate(motif_ids)
            }

    def _default_ml_cmap_name(motif_label: str, motif_pos: int) -> str:
        normalized = motif_label.strip().upper()
        if normalized == "A,0":
            return "Blues"
        if normalized == "CG,0":
            return "Oranges"
        sequential_cycle = (
            "Greens",
            "Purples",
            "Reds",
            "Greys",
            "YlGnBu",
            "BuPu",
            "PuRd",
            "YlOrBr",
            "BuGn",
        )
        return sequential_cycle[motif_pos % len(sequential_cycle)]

    ml_cmap_names = {
        motif_id: _default_ml_cmap_name(str(motif_labels_resolved[motif_id]), motif_pos)
        for motif_pos, motif_id in enumerate(motif_ids)
    }

    if vmin is None or vmax is None:
        if (
            render_mode == "scatter"
            and threshold_map is not None
            or render_mode == "scatter"
            and effective_scatter_color_values == "ml_score"
        ):
            vmin = 0.0 if vmin is None else vmin
            vmax = 1.0 if vmax is None else vmax
        else:
            vmax_auto = max(
                panel.max()
                for motif_panels in panels_sorted_by_motif
                for panel in motif_panels
            )
            vmin = 0.0 if vmin is None else vmin
            vmax = vmax_auto if vmax is None else vmax

    x_positions = np.arange(slice_width, dtype=float) - (slice_width // 2)

    if axis_orientation is None:
        if render_mode == "heatmap":
            effective_axis_orientation = "position_y" if rotate else "position_x"
        else:
            effective_axis_orientation = "position_x"
    else:
        effective_axis_orientation = axis_orientation

    share_x_axes = True
    if effective_axis_orientation == "position_x" and n_windows_effective > 1:
        # When each window has a different positional frame (relative-to-primary),
        # sharing x limits forces misleading axis alignment across panels.
        share_x_axes = False

    n_plot_motifs = len(motif_ids)
    if n_plot_motifs == 1:
        if rotate:
            fig, axes = plt.subplots(
                n_windows_effective,
                1,
                figsize=(max(6, slice_width * 0.01), 2.7 + n_windows_effective * 1.5),
                sharex=share_x_axes,
            )
            axes = np.atleast_1d(axes).reshape(n_windows_effective, 1)
        else:
            fig, axes = plt.subplots(
                1,
                n_windows_effective,
                figsize=(10, max(3, P * 0.06)),
                sharey=True,
            )
            axes = np.atleast_1d(axes).reshape(1, n_windows_effective)
    else:
        if rotate:
            n_rows, n_cols = n_windows_effective, n_plot_motifs
            fig, axes = plt.subplots(
                n_rows,
                n_cols,
                figsize=(max(8, 4.2 * n_cols), max(3.5, 2.3 * n_rows)),
                sharex=share_x_axes,
                sharey=True,
            )
        else:
            n_rows, n_cols = n_plot_motifs, n_windows_effective
            fig, axes = plt.subplots(
                n_rows,
                n_cols,
                figsize=(max(8, 3.8 * n_cols), max(3.5, 2.2 * n_rows)),
                sharey=True,
            )
        axes = np.atleast_2d(axes)

    def _axis_for(motif_pos: int, window_pos: int):
        if len(motif_ids) == 1:
            return axes[window_pos, 0] if rotate else axes[0, window_pos]
        if rotate:
            return axes[window_pos, motif_pos]
        return axes[motif_pos, window_pos]

    first_mappable = None
    threshold_legend_entries: dict[str, tuple[str, float]] = {}
    for motif_pos, motif_id in enumerate(motif_ids):
        motif_panels = panels_sorted_by_motif[motif_pos]
        motif_raw_panels = panels_raw_sorted_by_motif[motif_pos]
        motif_smooth_panels = panels_smooth_sorted_by_motif[motif_pos]
        motif_ml_panels = panels_ml_sorted_by_motif[motif_pos]
        motif_label = str(motif_labels_resolved[motif_id])
        for window_pos in range(n_windows_effective):
            ax = _axis_for(motif_pos, window_pos)
            panel = motif_panels[window_pos]
            panel_raw = motif_raw_panels[window_pos]
            panel_smooth = motif_smooth_panels[window_pos]
            panel_ml = motif_ml_panels[window_pos]
            row_offsets = center_offsets[window_pos]
            axis_center = float(np.nanmedian(row_offsets)) if row_offsets.size else 0.0
            window_half_width = panel_half_widths[window_pos]
            window_lo = axis_center - window_half_width
            window_hi = axis_center + window_half_width

            if render_mode == "scatter":
                rows, cols = np.nonzero(panel_raw > 0)
                if (
                    effective_scatter_color_values == "ml_score"
                    and panel_ml is not None
                ):
                    color_panel = panel_ml
                elif effective_scatter_color_values == "smoothed":
                    color_panel = panel_smooth
                else:
                    color_panel = panel_raw
                values = color_panel[rows, cols]
                if threshold_map is not None:
                    threshold_value = float(threshold_map[motif_id])
                    keep = panel_ml[rows, cols] >= threshold_value
                    rows = rows[keep]
                    cols = cols[keep]
                    values = panel_ml[rows, cols]
                    threshold_legend_entries[motif_label] = (
                        str(motif_color_map[motif_id]),
                        threshold_value,
                    )
                if rows.size > 0:
                    x_scatter = x_positions[cols] + row_offsets[rows]
                    in_window = (x_scatter >= window_lo) & (x_scatter <= window_hi)
                    rows = rows[in_window]
                    cols = cols[in_window]
                    values = values[in_window]
                    x_scatter = x_scatter[in_window]
                if rows.size > 0:
                    if effective_axis_orientation == "position_x":
                        if threshold_map is None:
                            scatter_cmap = (
                                ml_cmap_names[motif_id]
                                if effective_scatter_color_values == "ml_score"
                                else cmap
                            )
                            first_mappable = ax.scatter(
                                x_scatter,
                                rows,
                                c=values,
                                s=scatter_size,
                                alpha=scatter_alpha,
                                cmap=scatter_cmap,
                                vmin=vmin,
                                vmax=vmax,
                                linewidths=0,
                                rasterized=True,
                            )
                        else:
                            ax.scatter(
                                x_scatter,
                                rows,
                                c=motif_color_map[motif_id],
                                s=scatter_size,
                                alpha=scatter_alpha,
                                linewidths=0,
                                rasterized=True,
                            )
                        x_min = float(np.nanmin(x_scatter))
                        x_max = float(np.nanmax(x_scatter))
                    else:
                        if threshold_map is None:
                            scatter_cmap = (
                                ml_cmap_names[motif_id]
                                if effective_scatter_color_values == "ml_score"
                                else cmap
                            )
                            first_mappable = ax.scatter(
                                rows,
                                x_scatter,
                                c=values,
                                s=scatter_size,
                                alpha=scatter_alpha,
                                cmap=scatter_cmap,
                                vmin=vmin,
                                vmax=vmax,
                                linewidths=0,
                                rasterized=True,
                            )
                        else:
                            ax.scatter(
                                rows,
                                x_scatter,
                                c=motif_color_map[motif_id],
                                s=scatter_size,
                                alpha=scatter_alpha,
                                linewidths=0,
                                rasterized=True,
                            )
                        y_min = float(np.nanmin(x_scatter))
                        y_max = float(np.nanmax(x_scatter))
                else:
                    if effective_axis_orientation == "position_x":
                        x_min = float(window_lo)
                        x_max = float(window_hi)
                    else:
                        y_min = float(window_lo)
                        y_max = float(window_hi)
                if effective_axis_orientation == "position_x":
                    x_min = min(x_min, float(window_lo))
                    x_max = max(x_max, float(window_hi))
                    x_pad = max(1.0, 0.02 * (window_hi - window_lo + 1.0))
                    ax.axvline(axis_center, linestyle="--", linewidth=1, color="0.7")
                    ax.set_ylim(-1, panel.shape[0] + 1)
                    ax.set_xlim(window_lo - x_pad, window_hi + x_pad)
                else:
                    y_min = min(y_min, float(window_lo))
                    y_max = max(y_max, float(window_hi))
                    y_pad = max(1.0, 0.02 * (window_hi - window_lo + 1.0))
                    ax.axhline(axis_center, linestyle="--", linewidth=1, color="0.7")
                    ax.set_xlim(-1, panel.shape[0] + 1)
                    ax.set_ylim(window_lo - y_pad, window_hi + y_pad)
            else:
                x_lo = float(x_positions[0] + axis_center)
                x_hi = float(x_positions[-1] + axis_center)
                if effective_axis_orientation == "position_y":
                    first_mappable = ax.imshow(
                        panel.T,
                        aspect="auto",
                        origin="lower",
                        extent=[-0.5, panel.shape[0] - 0.5, x_lo, x_hi],
                        vmin=vmin,
                        vmax=vmax,
                        cmap=cmap,
                    )
                    ax.axhline(axis_center, linestyle="--", linewidth=1, color="0.7")
                    y_pad = max(1.0, 0.02 * (window_hi - window_lo + 1.0))
                    ax.set_ylim(window_lo - y_pad, window_hi + y_pad)
                else:
                    first_mappable = ax.imshow(
                        panel,
                        aspect="auto",
                        origin="upper",
                        extent=[x_lo, x_hi, panel.shape[0] - 0.5, -0.5],
                        vmin=vmin,
                        vmax=vmax,
                        cmap=cmap,
                    )
                    ax.axvline(axis_center, linestyle="--", linewidth=1, color="0.7")
                    x_pad = max(1.0, 0.02 * (window_hi - window_lo + 1.0))
                    ax.set_xlim(window_lo - x_pad, window_hi + x_pad)

            if render_mode == "heatmap":
                heatmap_title = (
                    f"{motif_label} | Window {window_pos + 1}"
                    if len(motif_ids) > 1
                    else f"Window {window_pos + 1}"
                )
                ax.set_title(heatmap_title, fontsize=9)
            elif len(motif_ids) > 1:
                if rotate and window_pos == 0:
                    ax.set_title(motif_label)
                if not rotate and motif_pos == 0:
                    ax.set_title(f"Window {window_pos + 1}")

            if effective_axis_orientation == "position_x":
                if rotate and motif_pos == 0:
                    ax.set_ylabel(f"Window {window_pos + 1}\nReads")
                if (not rotate) and len(motif_ids) > 1 and window_pos == 0:
                    ax.set_ylabel(f"{motif_label}\nReads")
            else:
                if rotate and motif_pos == 0:
                    ax.set_ylabel(f"Window {window_pos + 1}\nPosition")
                if (not rotate) and len(motif_ids) > 1 and window_pos == 0:
                    ax.set_ylabel(f"{motif_label}\nPosition")

    if rotate:
        bottom_axes = [
            _axis_for(motif_pos, n_windows_effective - 1)
            for motif_pos in range(len(motif_ids))
        ]
        position_label = "Position relative to primary center (bp)"
        xlabel = (
            position_label
            if effective_axis_orientation == "position_x"
            else "Reads (sorted)"
        )
        for ax in bottom_axes:
            ax.set_xlabel(xlabel)
    else:
        for motif_pos in range(len(motif_ids)):
            if effective_axis_orientation == "position_x":
                _axis_for(motif_pos, 0).set_ylabel(
                    "Reads (sorted)"
                    if len(motif_ids) == 1
                    else _axis_for(motif_pos, 0).get_ylabel()
                )
            else:
                _axis_for(motif_pos, 0).set_ylabel(
                    "Position (bp)"
                    if len(motif_ids) == 1
                    else _axis_for(motif_pos, 0).get_ylabel()
                )
            for window_pos in range(n_windows_effective):
                if effective_axis_orientation == "position_x":
                    _axis_for(motif_pos, window_pos).set_xlabel(
                        "Position relative to primary center (bp)"
                    )
                else:
                    _axis_for(motif_pos, window_pos).set_xlabel("Reads (sorted)")

    if threshold_map is None:
        if first_mappable is None:
            first_axis = _axis_for(0, 0)
            fallback_cmap = (
                ml_cmap_names[motif_ids[0]]
                if effective_scatter_color_values == "ml_score"
                else cmap
            )
            first_mappable = first_axis.scatter(
                [], [], c=[], cmap=fallback_cmap, vmin=vmin, vmax=vmax
            )
        cbar = fig.colorbar(
            first_mappable, ax=np.ravel(axes).tolist(), shrink=0.6, pad=0.02
        )
        sigma_txt = f", σ={smooth_sigma_bp}" if smoothing == "gaussian" else ""
        if render_mode == "scatter" and effective_scatter_color_values == "ml_score":
            cbar.set_label(
                f"normalized ML score (0-1)\\n{'downsampled' if downsampled else 'full'}"
            )
        elif render_mode == "scatter":
            scale_label = "fraction modified signal"
            smoothing_label = (
                f"{smoothing}, win={smooth_win}{sigma_txt}"
                if render_value_mode == "smoothed"
                else "none"
            )
            cbar.set_label(
                f"{scale_label} (values={render_value_mode}, smoother={smoothing_label})\\n"
                f"{'downsampled' if downsampled else 'full'}"
            )
        else:
            scale_label = "window signal heatmap"
            smoothing_label = (
                f"{smoothing}, win={smooth_win}{sigma_txt}"
                if render_value_mode == "smoothed"
                else "none"
            )
            cbar.set_label(
                f"{scale_label} (per-read smoothing along position axis; "
                f"values={render_value_mode}, smoother={smoothing_label})\\n"
                f"{'downsampled' if downsampled else 'full'}"
            )
    else:
        from matplotlib.lines import Line2D

        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                color="none",
                markerfacecolor=color,
                markeredgecolor=color,
                markersize=6,
                label=f"{label} (ML >= {threshold:.2f})",
                alpha=scatter_alpha,
            )
            for label, (color, threshold) in threshold_legend_entries.items()
        ]
        fig.legend(
            handles=handles,
            title="Motif thresholds",
            loc="upper right",
            bbox_to_anchor=(0.98, 0.98),
            frameon=False,
        )

    return fig, {
        "pairs": panels_sorted_by_motif[0][0].shape[0],
        "downsampled": downsampled,
        "downsample_method": (downsample_method if downsampled else "none"),
        "downsample_factor": int(downsample_factor),
        "rows_before_downsample": rows_before_downsample,
        "rows_after_downsample": int(panels_sorted_by_motif[0][0].shape[0]),
        "vmin": vmin,
        "vmax": vmax,
        "render_mode": render_mode,
        "render_values": render_value_mode,
        "scatter_color_values": effective_scatter_color_values,
        "axis_orientation": effective_axis_orientation,
        "x_axis_mode": effective_x_axis_mode,
        "sort_by": sort_by,
        "sort_window_index": sort_window_pos,
        "sort_descending": effective_sort_descending,
        "motifs_plotted": [motif_labels_resolved[m] for m in motif_ids],
        "n_windows": n_windows_effective,
        "window_centers_bp": target_offsets.tolist()
        if target_offsets is not None
        else None,
        "window_center_offsets": target_offsets.tolist()
        if target_offsets is not None
        else None,
        "window_widths_bp": [int(value) for value in panel_widths_bp],
        "required_read_length_bp": required_read_length_bp,
        "min_read_length_bp": effective_min_read_length,
        "window_match_tolerance_bp": match_tolerance_bp,
        "center_tolerance_bp": match_tolerance_bp,
        "enforce_full_window_span": bool(enforce_full_window_span),
        "dropped_for_window_span": int(dropped_for_window_span),
        "read_order_consistent": True,
        "ml_score_thresholds": (
            None
            if threshold_map is None
            else {
                motif_labels_resolved[key]: value
                for key, value in threshold_map.items()
                if key in motif_ids
            }
        ),
        "ml_score_cmaps": (
            {
                motif_labels_resolved[motif_id]: ml_cmap_names[motif_id]
                for motif_id in motif_ids
            }
            if threshold_map is None and effective_scatter_color_values == "ml_score"
            else None
        ),
    }


def plot_two_site_read_raster(
    read_windows: ReadWindowExtractionResult,
    *,
    second_site_offset_bp: int,
    window_width_bp: int | None = None,
    window_match_tolerance_bp: int | None = None,
    **kwargs,
):
    """
    Convenience wrapper for two-site raster plots in primary-region coordinates.

    Args:
        read_windows: ReadWindowExtractionResult from extract_read_windows / build_multimotif_read_windows
        second_site_offset_bp: secondary site center (bp) relative to the primary site at 0
        window_width_bp: common width (bp) applied to both windows when provided
        window_match_tolerance_bp: tolerance used when matching explicit centers
        **kwargs: forwarded to plot_multisite_read_raster
    """

    if "window_centers_bp" in kwargs or "window_center_offsets" in kwargs:
        raise ValueError(
            "plot_two_site_read_raster manages centers internally; do not pass window_centers_bp."
        )
    if "n_windows" in kwargs:
        raise ValueError("plot_two_site_read_raster always uses exactly two windows.")

    forwarded: dict[str, Any] = dict(kwargs)
    forwarded["window_centers_bp"] = [0, int(second_site_offset_bp)]
    if window_width_bp is not None:
        forwarded["window_widths_bp"] = [int(window_width_bp), int(window_width_bp)]
    if window_match_tolerance_bp is not None:
        forwarded["window_match_tolerance_bp"] = int(window_match_tolerance_bp)
    return plot_multisite_read_raster(read_windows, n_windows=2, **forwarded)


def plot_region_cluster_profiles(
    pileup_matrix: np.ndarray,
    labels: np.ndarray,
    *,
    window_size: int | None = None,
    motif_index: int | None = None,
    motif_count: int | None = None,
    plot_all_motifs: bool = False,
    motif_labels: Sequence[str] | None = None,
    motif_colors: Sequence[str] | None = None,
    motif_profile_mode: str = "single_axis",
    color_points_by: str = "cluster",
    point_size: float = 1.0,
    point_alpha: float = 0.01,
    smoothing: str | None = None,
    smooth_win: int = 21,
    smooth_sigma: float = 6.0,
    show_unsmoothed_overlay: bool = False,
    invert_y: bool = True,
    cmap_name: str = "viridis",
):
    """
    Visualize region-level clustering with a sorted matrix view and side profiles.
    Supports concatenated multi-motif matrices with shared or separate profile axes.

    Args:
        pileup_matrix: 2D array (n_regions, region_length) of modification fractions
        labels: cluster assignments per region
        window_size: optional half-window in bp for axis scaling and motif-slice inference
            (full span = ``2 * window_size``).
        motif_index: if pileup_matrix is a concatenation of multiple motifs, select which motif slice to plot
        motif_count: total number of motifs concatenated; if None, inferred when possible
        plot_all_motifs: if True and motifs are concatenated, plot all motif slices per cluster
        cmap_name: matplotlib colormap name for clusters
    """

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D

    X_full = np.asarray(pileup_matrix)
    labs = np.asarray(labels)
    if X_full.shape[0] != labs.shape[0]:
        raise ValueError("pileup_matrix rows and labels must have the same length.")
    if motif_profile_mode not in {"single_axis", "separate_axes"}:
        raise ValueError("motif_profile_mode must be 'single_axis' or 'separate_axes'.")
    if color_points_by not in {"cluster", "motif"}:
        raise ValueError("color_points_by must be 'cluster' or 'motif'.")

    region_len = X_full.shape[1]
    motif_idx = motif_index if motif_index is not None else 0
    n_motifs, slice_width = _resolve_motif_slices(
        region_len,
        motif_count=motif_count,
        motif_labels=motif_labels,
        window_size=window_size,
    )
    if motif_idx < 0 or motif_idx >= n_motifs:
        raise ValueError(
            f"motif_index {motif_idx} out of range for {n_motifs} motif slices."
        )
    motif_ids = (
        list(range(n_motifs)) if plot_all_motifs and n_motifs > 1 else [motif_idx]
    )

    if motif_labels is None:
        motif_labels = [f"motif_{i}" for i in range(n_motifs)]
    else:
        motif_labels = list(motif_labels)
        if len(motif_labels) < n_motifs:
            motif_labels.extend(
                [f"motif_{i}" for i in range(len(motif_labels), n_motifs)]
            )
    if motif_colors is None:
        motif_cmap = plt.get_cmap("tab10")
        motif_colors = [motif_cmap(i % 10) for i in range(n_motifs)]
    else:
        motif_colors = list(motif_colors)
        if len(motif_colors) < n_motifs:
            motif_cmap = plt.get_cmap("tab10")
            motif_colors.extend(
                motif_cmap(i % 10) for i in range(len(motif_colors), n_motifs)
            )

    primary_start = motif_idx * slice_width
    primary_end = min(region_len, primary_start + slice_width)
    primary_matrix = X_full[:, primary_start:primary_end]
    x_axis = _centered_x_axis(
        primary_matrix.shape[1],
        _window_span_from_size(window_size)
        if window_size is not None
        else primary_matrix.shape[1],
    )

    labs_arr, lab_codes, unique_codes = _prepare_group_labels(labs)
    order = np.lexsort((-primary_matrix.mean(axis=1), lab_codes))
    X_sorted_full = X_full[order]
    lab_codes_sorted = lab_codes[order]

    unique_labels = np.array(
        [np.unique(labs_arr[lab_codes == code])[0] for code in unique_codes]
    )
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=int(unique_codes.min()), vmax=int(unique_codes.max() or 1))

    fig = plt.figure(figsize=(12, max(4, len(unique_codes))))
    gs = gridspec.GridSpec(nrows=len(unique_codes), ncols=2, width_ratios=[3, 1])
    ax_left = fig.add_subplot(gs[:, 0])

    # Left panel: scatter only (parallel to read-clustering style)
    for motif_id in motif_ids:
        s = motif_id * slice_width
        e = min(region_len, s + slice_width)
        matrix_view = X_sorted_full[:, s:e]
        rows, cols = np.nonzero(matrix_view)
        if color_points_by == "cluster":
            colors = cmap(norm(lab_codes_sorted[rows]))
        else:
            colors = motif_colors[motif_id]
        ax_left.scatter(
            x_axis[cols],
            rows,
            s=point_size,
            alpha=point_alpha,
            c=colors,
            rasterized=True,
        )

    ax_left.set_xlabel("Distance from region center (bp)")
    ax_left.set_ylabel("Region (sorted)")
    if invert_y:
        ax_left.invert_yaxis()
    change_points = np.flatnonzero(np.diff(lab_codes_sorted)) + 1
    for cp in change_points:
        ax_left.axhline(cp, color="0.2", linestyle="--", linewidth=0.35)

    for i, code in enumerate(unique_codes):
        ax = fig.add_subplot(gs[i, 1])
        row_mask = lab_codes_sorted == code
        panel_max = 0.0
        if motif_profile_mode == "separate_axes" and len(motif_ids) > 1:
            axes_for_motifs = [ax]
            motif_maxima = [0.0 for _ in motif_ids]
            for motif_pos in range(1, len(motif_ids)):
                twin = ax.twinx()
                twin.spines["right"].set_position(("outward", 35 * motif_pos))
                axes_for_motifs.append(twin)
            for motif_pos, motif_id in enumerate(motif_ids):
                s = motif_id * slice_width
                e = min(region_len, s + slice_width)
                mean_profile_raw = X_sorted_full[row_mask, s:e].mean(axis=0)
                mean_profile = _smooth_profile_vector(
                    mean_profile_raw,
                    smoothing=smoothing,
                    smooth_win=smooth_win,
                    smooth_sigma=smooth_sigma,
                )
                motif_maxima[motif_pos] = max(
                    motif_maxima[motif_pos], float(np.nanmax(mean_profile_raw))
                )
                panel_max = max(panel_max, motif_maxima[motif_pos])
                motif_ax = axes_for_motifs[motif_pos]
                if show_unsmoothed_overlay and smoothing is not None:
                    motif_ax.plot(
                        x_axis,
                        mean_profile_raw,
                        color=motif_colors[motif_id],
                        linewidth=1.0,
                        alpha=0.25,
                    )
                motif_ax.plot(
                    x_axis,
                    mean_profile,
                    color=motif_colors[motif_id],
                    linewidth=1.5,
                )
                motif_ax.set_ylabel(
                    str(motif_labels[motif_id]), color=motif_colors[motif_id]
                )
                motif_ax.tick_params(
                    axis="y", colors=motif_colors[motif_id], labelsize=8
                )
            for motif_pos, motif_ax in enumerate(axes_for_motifs):
                motif_ax.set_ylim(0, max(motif_maxima[motif_pos], 0.05) * 1.05)
        else:
            for motif_id in motif_ids:
                s = motif_id * slice_width
                e = min(region_len, s + slice_width)
                mean_profile_raw = X_sorted_full[row_mask, s:e].mean(axis=0)
                mean_profile = _smooth_profile_vector(
                    mean_profile_raw,
                    smoothing=smoothing,
                    smooth_win=smooth_win,
                    smooth_sigma=smooth_sigma,
                )
                panel_max = max(panel_max, float(np.nanmax(mean_profile_raw)))
                if show_unsmoothed_overlay and smoothing is not None:
                    ax.plot(
                        x_axis,
                        mean_profile_raw,
                        color=motif_colors[motif_id]
                        if len(motif_ids) > 1
                        else cmap(norm(code)),
                        linewidth=1.0,
                        alpha=0.25,
                    )
                ax.plot(
                    x_axis,
                    mean_profile,
                    color=motif_colors[motif_id]
                    if len(motif_ids) > 1
                    else cmap(norm(code)),
                    linewidth=1.5,
                )
        ax.set_title(f"{unique_labels[i]} (n={row_mask.sum()})")
        ax.set_xlim(x_axis[0], x_axis[-1])
        y_max = max(panel_max, 0.05)
        ax.set_ylim(0, y_max * 1.05)

    if len(motif_ids) > 1 or color_points_by == "motif":
        handles = [
            Line2D(
                [0], [0], color=motif_colors[m], linewidth=2, label=str(motif_labels[m])
            )
            for m in motif_ids
        ]
        fig.legend(
            handles=handles,
            loc="upper center",
            ncol=min(4, len(handles)),
            frameon=False,
            bbox_to_anchor=(0.5, 1.01),
        )

    fig.tight_layout()
    return fig


def export_region_clusters_to_bed(
    region_metadata: Sequence[RegionMetadata],
    labels: Sequence[int],
    output_path: str | Path,
    *,
    name_prefix: str = "cluster",
    include_header: bool = False,
):
    """
    Export region cluster assignments as a BED-like file.

    Args:
        region_metadata: iterable of (chrom, start, end, strand) tuples, typically from region_feature_matrix_from_pileup
        labels: cluster labels aligned to region_metadata
        output_path: where to write the BED
        name_prefix: prefix for the BED name field; final name becomes f\"{prefix}_{label}\"
        include_header: if True, writes a comment header line
    """

    if len(region_metadata) != len(labels):
        raise ValueError("region_metadata and labels must have the same length.")

    output_path = Path(output_path)
    lines = []
    for (chrom, start, end, strand), label in zip(
        region_metadata, labels, strict=False
    ):
        name = f"{name_prefix}_{label}"
        score = str(label)
        strand_field = strand if strand in {"+", "-"} else "."
        lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand_field}\n")

    with output_path.open("w") as fh:
        if include_header:
            fh.write(f"# chrom\tstart\tend\tname({name_prefix}_id)\tscore\tstrand\n")
        fh.writelines(lines)

    return output_path


def plot_cluster_karyotype(
    region_bed: str | Path,
    chrom_sizes: str | Path,
    *,
    cmap_name: str = "tab10",
    linewidth: float = 4.0,
    figsize_per_chrom: float = 0.7,
    min_visible_bp: int = 50_000,
    min_visible_fraction: float = 0.002,
    chromosome_order: str | Sequence[str] = "length_desc",
    invert_position_axis: bool = True,
    detect_haplotype_backbone_shading: bool = True,
    show_position_axis: bool = False,
    coordinate_label_stagger: bool = True,
):
    """
    Plot cluster-labeled regions along chromosomes (ideogram-style), colored by cluster.

    Args:
        region_bed: BED file, typically from export_region_clusters_to_bed, with Name carrying cluster id
        chrom_sizes: 2-column TSV (chrom, length) e.g. ref.fasta.fai
        cmap_name: matplotlib colormap to map cluster ids to colors
        linewidth: thickness of region segments
        figsize_per_chrom: vertical size per chromosome for sizing the figure
        min_visible_bp: minimum displayed region length in bp for visibility
        min_visible_fraction: minimum displayed region length as a fraction of chromosome length
        chromosome_order: "length_desc", "length_asc", "natural", or explicit ordered chromosome list
        invert_position_axis: if True, chromosome position increases downward
        detect_haplotype_backbone_shading: if True, apply subtle maternal/paternal-like
            backbone shading when paired haplotype naming is detected
        show_position_axis: if True, show y-axis ticks/label; default hides them for cleaner ideogram view
        coordinate_label_stagger: if True, alternate end-label vertical offsets to reduce overlap
    """
    import re

    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    clusters_df = pd.read_csv(
        region_bed,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["Chromosome", "Start", "End", "Name"],
    )
    if clusters_df.empty:
        raise ValueError(f"No records found in region_bed: {region_bed}")
    chrom_df = pd.read_csv(
        chrom_sizes,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["Chromosome", "Length"],
    )
    if chrom_df.empty:
        raise ValueError(f"No chromosome sizes found in: {chrom_sizes}")
    df = clusters_df.merge(chrom_df, on="Chromosome", how="inner")
    if df.empty:
        raise ValueError("No overlapping chromosomes between BED and chrom_sizes.")

    df["Start"] = pd.to_numeric(df["Start"], errors="coerce")
    df["End"] = pd.to_numeric(df["End"], errors="coerce")
    df["Length"] = pd.to_numeric(df["Length"], errors="coerce")
    df = df.dropna(subset=["Start", "End", "Length"]).copy()
    if df.empty:
        raise ValueError("No valid numeric coordinates after parsing BED/chrom sizes.")

    def _natural_chrom_key(chrom: str):
        text = str(chrom)
        match = re.match(r"^chr(\d+)$", text, flags=re.IGNORECASE)
        if match:
            return (0, int(match.group(1)), text)
        if text.lower() in {"chrx", "x"}:
            return (1, 23, text)
        if text.lower() in {"chry", "y"}:
            return (1, 24, text)
        if text.lower() in {"chrm", "mt", "m"}:
            return (1, 25, text)
        return (2, text.lower(), text)

    def _split_haplotype(chrom: str) -> tuple[str, str | None]:
        name = str(chrom)
        patterns = (
            r"^(?P<base>.+?)[._-](?P<hap>mat|pat|maternal|paternal)$",
            r"^(?P<base>.+?)[._-](?P<hap>hap1|hap2)$",
        )
        for pattern in patterns:
            m = re.match(pattern, name, flags=re.IGNORECASE)
            if m:
                return (m.group("base"), m.group("hap").lower())
        return (name, None)

    parsed_cluster = pd.to_numeric(
        df["Name"].astype(str).str.extract(r"(\d+)$", expand=False),
        errors="coerce",
    )
    if parsed_cluster.isna().all():
        df["cluster_id"] = pd.factorize(df["Name"].astype(str))[0]
    else:
        fallback = pd.factorize(df["Name"].astype(str))[0]
        df["cluster_id"] = parsed_cluster.fillna(
            pd.Series(fallback, index=df.index)
        ).astype(int)

    chrom_sizes_series = df.groupby("Chromosome")["Length"].max()
    if isinstance(chromosome_order, (list, tuple, pd.Index, np.ndarray)):
        requested = [str(chrom) for chrom in chromosome_order]
        observed = chrom_sizes_series.index.astype(str).tolist()
        chrom_order = [chrom for chrom in requested if chrom in observed]
        chrom_order.extend([chrom for chrom in observed if chrom not in chrom_order])
    else:
        order_mode = str(chromosome_order).lower()
        if order_mode == "length_asc":
            chrom_order = (
                chrom_sizes_series.sort_values(ascending=True)
                .index.astype(str)
                .tolist()
            )
        elif order_mode == "natural":
            chrom_order = sorted(
                chrom_sizes_series.index.astype(str).tolist(), key=_natural_chrom_key
            )
        else:
            chrom_order = (
                chrom_sizes_series.sort_values(ascending=False)
                .index.astype(str)
                .tolist()
            )

    haplotype_pairs: dict[str, set[str]] = {}
    if detect_haplotype_backbone_shading:
        for chrom in chrom_order:
            base, hap = _split_haplotype(chrom)
            if hap is None:
                continue
            haplotype_pairs.setdefault(base, set()).add(hap)

    def _backbone_color(chrom: str) -> str:
        if not detect_haplotype_backbone_shading:
            return "0.75"
        base, hap = _split_haplotype(chrom)
        if hap is None:
            return "0.75"
        pair = haplotype_pairs.get(base, set())
        # Only apply shading when we detect at least a pair for the same base.
        if len(pair) < 2:
            return "0.75"
        if hap in {"maternal", "mat", "hap1"}:
            return "0.66"
        if hap in {"paternal", "pat", "hap2"}:
            return "0.84"
        return "0.75"

    cmap = plt.get_cmap(cmap_name)
    unique_clusters = np.sort(df["cluster_id"].unique())
    if hasattr(cmap, "colors") and len(cmap.colors) > 0:
        color_cycle = list(cmap.colors)
        cluster_to_color = {
            int(cid): color_cycle[i % len(color_cycle)]
            for i, cid in enumerate(unique_clusters)
        }
    else:
        norm = Normalize(vmin=0, vmax=max(1, len(unique_clusters) - 1))
        cluster_to_color = {
            int(cid): cmap(norm(i)) for i, cid in enumerate(unique_clusters)
        }

    fig_height = max(3.0, figsize_per_chrom * len(chrom_order))
    fig, ax = plt.subplots(figsize=(max(8, len(chrom_order) * 0.7), fig_height + 1.5))
    max_length = float(df["Length"].max())
    x_positions = np.arange(len(chrom_order), dtype=float)
    cluster_counts = df["cluster_id"].value_counts().to_dict()
    # Draw dense clusters first and rarer clusters on top for visibility.
    draw_order = sorted(
        unique_clusters,
        key=lambda cid: int(cluster_counts.get(int(cid), 0)),
        reverse=True,
    )

    for xi, chrom in zip(x_positions, chrom_order, strict=False):
        chrom_len = float(df.loc[df["Chromosome"] == chrom, "Length"].max())
        # Backbone for full chromosome (proportional in bp units)
        ax.plot(
            [xi, xi],
            [0, chrom_len],
            color=_backbone_color(chrom),
            lw=linewidth + 1.5,
            alpha=0.7,
            solid_capstyle="butt",
            zorder=1,
        )
        # Cluster-assigned regions
        sub = df[df["Chromosome"] == chrom].copy()
        for z_rank, cid in enumerate(draw_order, start=2):
            sub_cluster = sub[sub["cluster_id"] == int(cid)]
            if sub_cluster.empty:
                continue
            for _, row in sub_cluster.iterrows():
                start_bp = float(row["Start"])
                end_bp = float(row["End"])
                region_len = max(0.0, end_bp - start_bp)
                min_len = max(
                    float(min_visible_bp), chrom_len * float(min_visible_fraction)
                )
                if region_len < min_len:
                    center = 0.5 * (start_bp + end_bp)
                    start_bp = max(0.0, center - (0.5 * min_len))
                    end_bp = min(chrom_len, center + (0.5 * min_len))
                ax.plot(
                    [xi, xi],
                    [start_bp, end_bp],
                    color=cluster_to_color[int(row.cluster_id)],
                    lw=linewidth + 2.0,
                    alpha=0.98,
                    solid_capstyle="butt",
                    zorder=z_rank,
                )

        # Per-chromosome coordinate labels centered on each ideogram.
        # Use per-chromosome padding (not global max length) so labels stay close.
        x_label = xi
        y_pad = max(1.0, 0.0015 * chrom_len)
        if coordinate_label_stagger:
            bottom_extra = 0.0 if (int(xi) % 2 == 0) else (0.6 * y_pad)
        else:
            bottom_extra = 0.0
        if invert_position_axis:
            top_label_y = -y_pad
            bottom_label_y = chrom_len + y_pad + bottom_extra
            top_label_text = "0"
            bottom_label_text = f"{int(chrom_len):,}"
        else:
            top_label_y = chrom_len + y_pad + bottom_extra
            bottom_label_y = -y_pad
            top_label_text = f"{int(chrom_len):,}"
            bottom_label_text = "0"
        ax.text(
            x_label,
            top_label_y,
            top_label_text,
            ha="center",
            va="bottom",
            fontsize=8,
            color="0.35",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.2},
            zorder=4,
        )
        ax.text(
            x_label,
            bottom_label_y,
            bottom_label_text,
            ha="center",
            va="top",
            fontsize=8,
            color="0.35",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.2},
            zorder=4,
        )

    ax.set_xlim(-0.7, len(chrom_order) - 0.3)
    if invert_position_axis:
        ax.set_ylim(1.08 * max_length, -0.1 * max_length)
    else:
        ax.set_ylim(-0.1 * max_length, 1.08 * max_length)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(chrom_order, rotation=45, ha="right")
    if show_position_axis:
        ax.set_ylabel("Position (bp)")
    else:
        ax.set_ylabel("")
        ax.yaxis.set_visible(False)
        ax.spines["left"].set_visible(False)
    ax.set_title("Clustered Regions by Chromosome (Vertical, Proportional Length)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=cluster_to_color[int(cid)],
            linewidth=3,
            label=f"C{int(cid)}",
        )
        for cid in unique_clusters
    ]
    ax.legend(
        handles=legend_handles,
        title="Cluster",
        loc="upper right",
        frameon=False,
        ncol=min(4, len(legend_handles)),
    )
    plt.tight_layout()
    return fig


def summarize_region_cluster_annotations(
    region_metadata: Sequence[RegionMetadata],
    labels: Sequence[int],
    annotation_file: str | Path,
    *,
    feature_field: str | None = None,
    name_prefix: str = "cluster",
) -> pd.DataFrame:
    """
    Intersect clustered regions with genome annotations and return per-cluster tallies.

    Args:
        region_metadata: iterable of (chrom, start, end, strand) tuples
        labels: cluster labels aligned to region_metadata
        annotation_file: path to a BED/GTF/GFF file with annotations
        feature_field: column to use as the feature/category. For GTF/GFF this defaults
            to "Feature" if present. For BED, this defaults to the Name/score column if
            available, else "annotation".
        name_prefix: expected prefix for cluster names in the BED export (used to strip
            when reading Name)

    Returns:
        pandas DataFrame with counts per (cluster, feature) and a column "total_cluster"
        plus a fraction-normalized view ("fraction").

    Note:
        Requires pyranges (pip install pyranges). This keeps everything in Python without
        shelling out to bedtools.
    """

    try:
        import pyranges as pr
    except Exception as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "summarize_region_cluster_annotations requires pyranges. "
            "Install it with `pip install pyranges`."
        ) from exc

    if len(region_metadata) != len(labels):
        raise ValueError("region_metadata and labels must have the same length.")

    # Build PyRanges for clusters
    cluster_df = pd.DataFrame(
        {
            "Chromosome": [c for c, _, _, _ in region_metadata],
            "Start": [s for _, s, _, _ in region_metadata],
            "End": [e for _, _, e, _ in region_metadata],
            "Strand": [st if st in {"+", "-"} else "." for *_, st in region_metadata],
            "Name": [f"{name_prefix}_{lab}" for lab in labels],
            "Score": labels,
        }
    )
    clusters = pr.PyRanges(cluster_df)

    # Read annotations
    annotation_file = Path(annotation_file)
    if annotation_file.suffix.lower() in {".gtf", ".gff", ".gff3"}:
        annot = pr.read_gff(annotation_file)
    else:
        annot = pr.read_bed(annotation_file)

    # Pick feature column
    annot_df = annot.df
    if feature_field and feature_field in annot_df.columns:
        feature_col = feature_field
    elif "Feature" in annot_df.columns:
        feature_col = "Feature"
    elif "Name" in annot_df.columns:
        feature_col = "Name"
    else:
        feature_col = None
        annot_df["feature_tmp"] = "annotation"
        feature_col = "feature_tmp"
        annot = pr.PyRanges(annot_df)

    # Intersect and count
    overlaps = clusters.join(annot)
    if overlaps.df.empty:
        return pd.DataFrame()

    df = overlaps.df.copy()
    df["cluster"] = (
        df["Name"]
        .astype(str)
        .str.replace(f"{name_prefix}_", "", regex=False)
        .astype(int)
    )
    df["feature"] = df[feature_col].astype(str)

    counts = df.groupby(["cluster", "feature"]).size().reset_index(name="count")
    totals = counts.groupby("cluster")["count"].sum().reset_index(name="total_cluster")
    merged = counts.merge(totals, on="cluster", how="left")
    merged["fraction"] = merged["count"] / merged["total_cluster"]
    return merged


def _adjust_p_values_bh(p_values: pd.Series) -> pd.Series:
    if p_values.empty:
        return pd.Series(dtype=float, index=p_values.index)

    numeric = pd.to_numeric(p_values, errors="coerce")
    valid_mask = numeric.notna()
    if not valid_mask.any():
        return pd.Series(np.nan, index=p_values.index, dtype=float)

    valid_values = numeric[valid_mask].to_numpy(dtype=float)
    order = np.argsort(valid_values, kind="mergesort")
    ranked = valid_values[order]
    adjusted = np.empty_like(ranked)
    running_min = 1.0
    total = len(ranked)

    for index in range(total - 1, -1, -1):
        rank = index + 1
        candidate = min(1.0, ranked[index] * total / rank)
        running_min = min(running_min, candidate)
        adjusted[index] = running_min

    restored = np.full(len(p_values), np.nan, dtype=float)
    restored_indices = np.flatnonzero(valid_mask.to_numpy())
    restored[restored_indices[order]] = adjusted
    return pd.Series(restored, index=p_values.index, dtype=float)


def _binomial_greater_tail(k: int, n: int, p: float) -> float:
    if n <= 0:
        return 1.0
    if k <= 0:
        return 1.0
    if p <= 0.0:
        return 0.0 if k > 0 else 1.0
    if p >= 1.0:
        return 1.0
    if k > n:
        raise ValueError("k cannot exceed n for a binomial tail probability.")

    log_p = math.log(p)
    log_q = math.log1p(-p)
    tail = 0.0
    for i in range(k, n + 1):
        log_prob = (
            math.lgamma(n + 1)
            - math.lgamma(i + 1)
            - math.lgamma(n - i + 1)
            + i * log_p
            + (n - i) * log_q
        )
        tail += math.exp(log_prob)
    return float(min(max(tail, 0.0), 1.0))


def summarize_read_cluster_region_associations(
    metadata: Sequence[dict[str, Any]],
    labels: Sequence[int | str],
    *,
    include_strand: bool = True,
    min_reads_per_region: int = 1,
    prior_count: float = 0.5,
) -> pd.DataFrame:
    """
    Summarize read-cluster associations per region with enrichment statistics.

    Returns a long-form DataFrame with one row per (region, cluster).
    """

    if len(metadata) != len(labels):
        raise ValueError("metadata and labels must have the same length.")
    if min_reads_per_region < 1:
        raise ValueError("min_reads_per_region must be at least 1.")
    if prior_count < 0:
        raise ValueError("prior_count must be non-negative.")

    rows: list[dict[str, Any]] = []
    for meta, cluster_label in zip(metadata, labels, strict=False):
        chrom = meta.get("chromosome", meta.get("chrom"))
        start = meta.get("region_start", meta.get("start"))
        end = meta.get("region_end", meta.get("end"))
        strand_value = (
            meta.get("region_strand", meta.get("strand")) if include_strand else "."
        )
        strand = strand_value if strand_value in {"+", "-"} else "."
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "cluster": cluster_label,
            }
        )

    columns = [
        "chrom",
        "start",
        "end",
        "strand",
        "cluster",
        "count",
        "total_reads",
        "fraction",
        "global_count",
        "global_total",
        "global_fraction",
        "log2_enrichment",
        "p_value",
        "q_value",
    ]
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "cluster"])
    if df.empty:
        return pd.DataFrame(columns=columns)

    region_cols = ["chrom", "start", "end", "strand"]
    cluster_order = list(pd.unique(df["cluster"]))
    region_totals = (
        df.groupby(region_cols, sort=True).size().reset_index(name="total_reads")
    )
    kept_regions = region_totals.loc[
        region_totals["total_reads"] >= min_reads_per_region, region_cols
    ]
    if kept_regions.empty:
        return pd.DataFrame(columns=columns)

    region_grid = (
        kept_regions.assign(_key=1)
        .merge(
            pd.DataFrame({"cluster": cluster_order, "_key": 1}),
            on="_key",
            how="outer",
        )
        .drop(columns="_key")
    )

    observed_counts = (
        df.groupby(region_cols + ["cluster"], sort=True)
        .size()
        .reset_index(name="count")
    )
    global_counts = (
        df.groupby("cluster", sort=True)
        .size()
        .reindex(cluster_order, fill_value=0)
        .reset_index(name="global_count")
    )
    global_total = int(len(df))
    global_counts["global_total"] = global_total
    global_counts["global_fraction"] = (
        global_counts["global_count"] / global_total if global_total > 0 else 0.0
    )

    summary = region_grid.merge(
        observed_counts, on=region_cols + ["cluster"], how="left"
    )
    summary["count"] = summary["count"].fillna(0).astype(int)
    summary = summary.merge(region_totals, on=region_cols, how="left")
    summary = summary.merge(global_counts, on="cluster", how="left")
    summary["fraction"] = summary["count"] / summary["total_reads"]

    cluster_count = max(len(cluster_order), 1)
    if prior_count == 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            summary["log2_enrichment"] = np.log2(
                (summary["count"] / summary["total_reads"])
                / (summary["global_count"] / summary["global_total"])
            )
    else:
        region_shrunk = (summary["count"] + prior_count) / (
            summary["total_reads"] + (prior_count * cluster_count)
        )
        global_shrunk = (summary["global_count"] + prior_count) / (
            summary["global_total"] + (prior_count * cluster_count)
        )
        summary["log2_enrichment"] = np.log2(region_shrunk / global_shrunk)

    summary["p_value"] = [
        _binomial_greater_tail(int(count), int(total), float(prob))
        for count, total, prob in zip(
            summary["count"],
            summary["total_reads"],
            summary["global_fraction"].fillna(0.0),
            strict=False,
        )
    ]
    summary["q_value"] = _adjust_p_values_bh(summary["p_value"])
    summary = summary.sort_values(region_cols + ["cluster"]).reset_index(drop=True)
    return summary[columns]


def summarize_read_clusters_by_region(
    metadata: Sequence[dict[str, Any]],
    labels: Sequence[int],
    *,
    include_strand: bool = True,
) -> pd.DataFrame:
    """
    Summarize read-level cluster composition per locus to assess homogeneity.

    Args:
        metadata: list of metadata dicts returned by extract_read_windows
        labels: cluster labels aligned to metadata entries
        include_strand: if True, region keys include strand; set False to pool strands

    Returns:
        DataFrame with one row per region key containing:
            - counts per cluster
            - total_reads
            - dominant_cluster
            - dominant_fraction
            - entropy (natural log base)
    """

    if len(metadata) != len(labels):
        raise ValueError("metadata and labels must have the same length.")
    summary = summarize_read_cluster_region_associations(
        metadata,
        labels,
        include_strand=include_strand,
    )
    if summary.empty:
        return pd.DataFrame(
            columns=[
                "chrom",
                "start",
                "end",
                "strand",
                "total_reads",
                "dominant_cluster",
                "dominant_fraction",
                "entropy",
            ]
        )

    pivot = summary.pivot_table(
        index=["chrom", "start", "end", "strand"],
        columns="cluster",
        values="count",
        fill_value=0,
        aggfunc="sum",
    )
    pivot["total_reads"] = pivot.sum(axis=1)

    # Homogeneity metrics per region
    def _dominance_and_entropy(row: pd.Series) -> tuple[int, float, float]:
        cluster_counts = row[row.index != "total_reads"].astype(float)
        total = cluster_counts.sum()
        if total == 0:
            return -1, 0.0, 0.0
        dominant_cluster = int(cluster_counts.idxmax())
        dominant_fraction = float(cluster_counts.max() / total)
        probs = cluster_counts / total
        entropy = float(-(probs * np.log(probs + 1e-12)).sum())
        return dominant_cluster, dominant_fraction, entropy

    dom, dom_frac, ent = zip(*pivot.apply(_dominance_and_entropy, axis=1), strict=False)
    pivot["dominant_cluster"] = dom
    pivot["dominant_fraction"] = dom_frac
    pivot["entropy"] = ent

    pivot = pivot.reset_index()
    return pivot


__all__ = [
    "cluster_features",
    "read_mod_fraction_table",
    "merge_read_window_results",
    "extract_read_windows",
    "build_multimotif_read_windows",
    "read_window_feature_matrix",
    "cluster_read_windows",
    "classify_read_features_binary",
    "plot_cluster_profiles",
    "region_feature_matrix_from_pileup",
    "plot_region_cluster_profiles",
    "export_region_clusters_to_bed",
    "summarize_region_cluster_annotations",
    "summarize_read_cluster_region_associations",
    "summarize_read_clusters_by_region",
    "plot_cluster_karyotype",
    "plot_classification_profiles",
    "plot_region_classification_profiles",
    "plot_confusion_matrices",
    "sample_rows",
    "cluster_label_mapping",
    "apply_cluster_label_mapping",
    "plot_two_site_read_raster",
    "plot_multisite_read_raster",
]
