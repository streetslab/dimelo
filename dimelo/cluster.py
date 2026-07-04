from __future__ import annotations

import contextlib
import math
import tempfile
import warnings
from collections import Counter, defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from functools import partial
from itertools import combinations
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

DEFAULT_RICH_FEATURE_SPEC: dict[str, Any] = {
    "motif_mode": "per_motif",
    "motif_count": 1,
    "motif_labels": None,
    "pooled": True,
    "pca": {"enabled": True, "n_components": 6, "scope": "pooled"},
    "summary": {"enabled": True},
    "densities": {"enabled": True, "windows": DEFAULT_DENSITY_WINDOWS},
    "asymmetry": {"enabled": True, "spans": (100, 300, 500)},
    "center_edge": {"enabled": True, "center_bp": 100, "edge_bp": 200},
    "autocorr": {"enabled": True, "lags": DEFAULT_AUTOCORR_LAGS},
    "fft": {"enabled": True, "periods_bp": (10, 50, 100, 150, 200)},
    "peaks": {"enabled": False, "prominence": 0.005},
    "cross_motif": {"enabled": True},
}


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


@dataclass(frozen=True)
class _RasterSiteSelection:
    mode: str
    n_windows: int
    min_distance_bp: int
    max_distance_bp: int | None
    selection_multiplicity: str
    choose: str
    selection_seed: int | None
    anchor: dict[str, Any]
    strand_relation: str
    exclude: dict[str, Any] | None
    orientation: str
    window_offsets_bp: tuple[int, ...] | None
    primary_window_index: int


@dataclass(frozen=True)
class _RasterSiteCandidate:
    row_idx: int
    center_bp: float
    chrom: str | None
    start_bp: int | None
    end_bp: int | None
    strand: str | None
    read_key: tuple[Any, Any]


@dataclass
class _RasterSiteSelectionResult:
    window_indices: list[np.ndarray]
    panel_centers: list[np.ndarray]
    stats: dict[str, Any]


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


def _resolve_raster_site_selection(
    site_selection: dict[str, Any] | None,
    *,
    selection_mode: str,
    window_offsets_bp: Sequence[int] | None,
    n_windows: int,
    min_separation_bp: int,
    primary_window_index: int,
    panel_widths_bp: Sequence[int],
) -> _RasterSiteSelection:
    spec = dict(site_selection or {})
    mode_aliases = {
        "cooccurring_regions": "cooccurring",
        "cooccurring": "cooccurring",
        "fixed_offsets": "fixed_offsets",
        "anchor_plus_neighbors": "anchor_plus_neighbors",
    }
    mode_raw = spec.get("mode", selection_mode)
    if mode_raw not in mode_aliases:
        raise ValueError(
            "site_selection mode must be 'cooccurring', 'fixed_offsets', or 'anchor_plus_neighbors'."
        )
    mode = mode_aliases[str(mode_raw)]

    offsets_value = spec.get("window_offsets_bp", window_offsets_bp)
    offsets = None
    if offsets_value is not None:
        offsets = tuple(int(value) for value in offsets_value)
        if not offsets:
            raise ValueError("window_offsets_bp must contain at least one offset.")

    if mode == "fixed_offsets":
        if offsets is None:
            raise ValueError(
                "window_offsets_bp is required for fixed_offsets site selection."
            )
        n_windows_effective = len(offsets)
    else:
        if offsets is not None:
            raise ValueError(
                "window_offsets_bp is only valid for fixed_offsets site selection."
            )
        n_windows_effective = int(spec.get("n_windows", n_windows))
        if n_windows_effective < 1:
            raise ValueError("site_selection n_windows must be >= 1.")

    primary_index = int(spec.get("primary_window_index", primary_window_index))
    if primary_index < 0 or primary_index >= n_windows_effective:
        raise ValueError(
            f"primary_window_index {primary_index} out of range for {n_windows_effective} windows."
        )

    default_min_distance: int | str = (
        "window_width" if site_selection is not None else min_separation_bp
    )
    min_distance_value = spec.get("min_distance_bp", default_min_distance)
    if min_distance_value == "window_width":
        min_distance = int(max(panel_widths_bp))
    else:
        min_distance = int(min_distance_value)
    if min_distance < 0:
        raise ValueError("site_selection min_distance_bp must be >= 0.")

    max_distance_value = spec.get("max_distance_bp")
    max_distance = None if max_distance_value is None else int(max_distance_value)
    if max_distance is not None and max_distance < min_distance:
        raise ValueError("site_selection max_distance_bp must be >= min_distance_bp.")

    multiplicity = str(spec.get("selection_multiplicity", "one_per_read"))
    if multiplicity == "all_valid_sets":
        raise NotImplementedError(
            "site_selection selection_multiplicity='all_valid_sets' is not implemented yet; use 'one_per_read'."
        )
    if multiplicity != "one_per_read":
        raise ValueError(
            "site_selection selection_multiplicity must be 'one_per_read' or 'all_valid_sets'."
        )

    choose = str(spec.get("choose", "first"))
    if choose not in {"first", "random", "longest_span", "shortest_span"}:
        raise ValueError(
            "site_selection choose must be 'first', 'random', 'longest_span', or 'shortest_span'."
        )

    anchor = dict(spec.get("anchor", {"mode": "first"}) or {"mode": "first"})
    anchor_mode = str(anchor.get("mode", "first"))
    if anchor_mode not in {"first", "random", "index"}:
        raise ValueError("site_selection anchor mode must be 'first', 'random', or 'index'.")

    strand_relation = str(spec.get("strand_relation", "any"))
    if strand_relation not in {"any", "same", "opposite"}:
        raise ValueError("site_selection strand_relation must be 'any', 'same', or 'opposite'.")

    orientation = str(spec.get("orientation", "genomic"))
    if orientation not in {"genomic", "anchor_strand"}:
        raise ValueError("site_selection orientation must be 'genomic' or 'anchor_strand'.")

    exclude = spec.get("exclude")
    exclude_dict = dict(exclude) if exclude is not None else None
    selection_seed = spec.get("selection_seed")

    return _RasterSiteSelection(
        mode=mode,
        n_windows=n_windows_effective,
        min_distance_bp=min_distance,
        max_distance_bp=max_distance,
        selection_multiplicity=multiplicity,
        choose=choose,
        selection_seed=None if selection_seed is None else int(selection_seed),
        anchor=anchor,
        strand_relation=strand_relation,
        exclude=exclude_dict,
        orientation=orientation,
        window_offsets_bp=offsets,
        primary_window_index=primary_index,
    )


def _metadata_center_bp(meta_row: dict[str, Any], distance_mode: str) -> float:
    if distance_mode == "center":
        return float((int(meta_row.get("region_start", 0)) + int(meta_row.get("region_end", 0))) // 2)
    return float(int(meta_row.get("read_start", 0)))


def _metadata_interval(meta_row: dict[str, Any]) -> tuple[str | None, int | None, int | None]:
    chrom = meta_row.get("chromosome")
    try:
        start = int(meta_row.get("region_start"))
        end = int(meta_row.get("region_end"))
    except Exception:
        return chrom, None, None
    return chrom, start, end


def _interval_overlaps_any(
    chrom: str | None,
    start: int | None,
    end: int | None,
    regions_dict: dict[str, list[tuple[int, int, str]]],
) -> bool:
    if chrom is None or start is None or end is None:
        return False
    for region_start, region_end, _ in regions_dict.get(str(chrom), []):
        if max(start, int(region_start)) < min(end, int(region_end)):
            return True
    return False


def _build_raster_site_candidates(
    meta: Sequence[dict[str, Any]],
    *,
    original_indices: np.ndarray,
    selector: _RasterSiteSelection,
    distance_mode: str,
) -> tuple[list[_RasterSiteCandidate], int]:
    excluded_original_rows = set()
    excluded_regions = None
    if selector.exclude:
        excluded_original_rows = {
            int(value) for value in selector.exclude.get("row_indices", []) or []
        }
        regions = selector.exclude.get("regions")
        if regions is not None:
            excluded_regions = utils.regions_dict_from_input(regions, window_size=None)

    candidates: list[_RasterSiteCandidate] = []
    excluded_count = 0
    for row_idx, meta_row in enumerate(meta):
        original_idx = int(original_indices[row_idx])
        chrom, start, end = _metadata_interval(meta_row)
        if original_idx in excluded_original_rows or (
            excluded_regions is not None
            and _interval_overlaps_any(chrom, start, end, excluded_regions)
        ):
            excluded_count += 1
            continue
        candidates.append(
            _RasterSiteCandidate(
                row_idx=row_idx,
                center_bp=_metadata_center_bp(meta_row, distance_mode),
                chrom=str(chrom) if chrom is not None else None,
                start_bp=start,
                end_bp=end,
                strand=_coerce_strand(meta_row.get("region_strand")),
                read_key=(meta_row.get("read_name"), meta_row.get("chromosome")),
            )
        )
    return candidates, excluded_count


def _strand_relation_ok(
    selected: Sequence[_RasterSiteCandidate],
    selector: _RasterSiteSelection,
) -> bool:
    if selector.strand_relation == "any":
        return True
    strands = [candidate.strand for candidate in selected]
    if any(strand not in {"+", "-"} for strand in strands):
        raise ValueError(
            "site_selection strand_relation requires '+'/'-' region_strand metadata for all selected sites."
        )
    if selector.strand_relation == "same":
        return len(set(strands)) == 1
    anchor_strand = strands[selector.primary_window_index]
    if len(strands) == 2:
        return strands[0] != strands[1]
    return any(strand != anchor_strand for strand in strands)


def _distances_ok(
    selected: Sequence[_RasterSiteCandidate],
    selector: _RasterSiteSelection,
) -> bool:
    centers = [candidate.center_bp for candidate in selected]
    if len(centers) < 2:
        return True
    for left, right in zip(centers, centers[1:], strict=False):
        distance = abs(float(right) - float(left))
        if distance < selector.min_distance_bp:
            return False
        if selector.max_distance_bp is not None and distance > selector.max_distance_bp:
            return False
    return True


def _apply_anchor_orientation(
    selected: Sequence[_RasterSiteCandidate],
    selector: _RasterSiteSelection,
) -> list[_RasterSiteCandidate]:
    ordered = list(selected)
    if (
        selector.orientation == "anchor_strand"
        and ordered
        and ordered[selector.primary_window_index].strand == "-"
    ):
        ordered = list(reversed(ordered))
    return ordered


def _fixed_offset_sets(
    items: Sequence[_RasterSiteCandidate],
    selector: _RasterSiteSelection,
) -> list[list[_RasterSiteCandidate]]:
    assert selector.window_offsets_bp is not None
    out: list[list[_RasterSiteCandidate]] = []
    by_center: defaultdict[float, list[_RasterSiteCandidate]] = defaultdict(list)
    for item in items:
        by_center[item.center_bp].append(item)
    for anchor in items:
        selected: list[_RasterSiteCandidate] = []
        used: set[int] = set()
        for offset in selector.window_offsets_bp:
            matches = by_center.get(anchor.center_bp + float(offset), [])
            chosen = next((candidate for candidate in matches if candidate.row_idx not in used), None)
            if chosen is None:
                selected = []
                break
            selected.append(chosen)
            used.add(chosen.row_idx)
        if selected and _strand_relation_ok(selected, selector):
            out.append(_apply_anchor_orientation(selected, selector))
    return out


def _anchor_position(
    items: Sequence[_RasterSiteCandidate],
    selector: _RasterSiteSelection,
    rng: np.random.Generator,
) -> int:
    mode = str(selector.anchor.get("mode", "first"))
    if mode == "first":
        return 0
    if mode == "random":
        return int(rng.integers(0, len(items)))
    index = int(selector.anchor.get("index", 0))
    if index < 0 or index >= len(items):
        raise ValueError("site_selection anchor index is out of range for a read.")
    return index


def _cooccurring_sets(
    items: Sequence[_RasterSiteCandidate],
    selector: _RasterSiteSelection,
    rng: np.random.Generator,
) -> list[list[_RasterSiteCandidate]]:
    out: list[list[_RasterSiteCandidate]] = []
    ordered_items = list(items)
    if selector.mode == "anchor_plus_neighbors":
        anchor_pos = _anchor_position(ordered_items, selector, rng)
        anchor = ordered_items[anchor_pos]
        combos_source = [
            [anchor, *combo]
            for combo in combinations(
                [item for idx, item in enumerate(ordered_items) if idx != anchor_pos],
                max(0, selector.n_windows - 1),
            )
        ]
    else:
        combos_source = combinations(ordered_items, selector.n_windows)

    for combo in combos_source:
        selected = sorted(combo, key=lambda candidate: candidate.center_bp)
        if not _distances_ok(selected, selector):
            continue
        if not _strand_relation_ok(selected, selector):
            continue
        out.append(_apply_anchor_orientation(selected, selector))
    return out


def _choose_site_set(
    valid_sets: Sequence[list[_RasterSiteCandidate]],
    selector: _RasterSiteSelection,
    rng: np.random.Generator,
) -> list[_RasterSiteCandidate] | None:
    if not valid_sets:
        return None
    if selector.choose == "first":
        return list(valid_sets[0])
    if selector.choose == "random":
        return list(valid_sets[int(rng.integers(0, len(valid_sets)))])

    def span(site_set: Sequence[_RasterSiteCandidate]) -> float:
        centers = [candidate.center_bp for candidate in site_set]
        return float(max(centers) - min(centers))

    if selector.choose == "longest_span":
        return list(max(valid_sets, key=span))
    return list(min(valid_sets, key=span))


def _summarize_observed_offsets(
    observed_offsets: Sequence[Sequence[float]],
) -> list[dict[str, float | int | None]]:
    summary: list[dict[str, float | int | None]] = []
    for values in observed_offsets:
        arr = np.asarray(values, dtype=float)
        if arr.size == 0:
            summary.append(
                {
                    "n": 0,
                    "min": None,
                    "median": None,
                    "max": None,
                    "unique": 0,
                }
            )
            continue
        summary.append(
            {
                "n": int(arr.size),
                "min": float(np.min(arr)),
                "median": float(np.median(arr)),
                "max": float(np.max(arr)),
                "unique": int(np.unique(arr).size),
            }
        )
    return summary


def _select_raster_site_windows(
    meta: Sequence[dict[str, Any]],
    *,
    original_indices: np.ndarray,
    selector: _RasterSiteSelection,
    distance_mode: str,
    selection_spans_all_windows,
) -> _RasterSiteSelectionResult:
    candidates, excluded_count = _build_raster_site_candidates(
        meta,
        original_indices=original_indices,
        selector=selector,
        distance_mode=distance_mode,
    )
    groups: defaultdict[tuple[Any, Any], list[_RasterSiteCandidate]] = defaultdict(list)
    for candidate in candidates:
        groups[candidate.read_key].append(candidate)

    rng = np.random.default_rng(selector.selection_seed)
    per_window_indices: list[list[int]] = [[] for _ in range(selector.n_windows)]
    per_window_centers: list[list[float]] = [[] for _ in range(selector.n_windows)]
    selected_read_keys: list[tuple[Any, Any]] = []
    dropped_for_window_span = 0
    observed_offsets: list[list[float]] = [[] for _ in range(selector.n_windows)]

    for read_key, items in groups.items():
        if len(items) < selector.n_windows:
            continue
        ordered_items = sorted(items, key=lambda candidate: candidate.center_bp)
        if selector.mode == "fixed_offsets":
            valid_sets = _fixed_offset_sets(ordered_items, selector)
        else:
            valid_sets = _cooccurring_sets(ordered_items, selector, rng)
        chosen = _choose_site_set(valid_sets, selector, rng)
        if chosen is None:
            continue
        candidate_for_span = [
            (candidate.row_idx, candidate.center_bp, window_pos)
            for window_pos, candidate in enumerate(chosen)
        ]
        if not selection_spans_all_windows(candidate_for_span):
            dropped_for_window_span += 1
            continue
        primary_center = chosen[selector.primary_window_index].center_bp
        for window_pos, candidate in enumerate(chosen):
            per_window_indices[window_pos].append(int(candidate.row_idx))
            per_window_centers[window_pos].append(float(candidate.center_bp))
            observed_offsets[window_pos].append(float(candidate.center_bp - primary_center))
        selected_read_keys.append(read_key)

    if not all(per_window_indices):
        raise ValueError("No read sets found meeting site-selection criteria.")

    stats = {
        "rows_are": "reads",
        "unique_reads": len(set(selected_read_keys)),
        "site_sets": len(selected_read_keys),
        "orientation_applied": selector.orientation,
        "observed_window_center_offsets_bp": observed_offsets,
        "observed_window_center_offsets_summary_bp": _summarize_observed_offsets(
            observed_offsets
        ),
        "site_selection": {
            "mode": selector.mode,
            "n_windows": selector.n_windows,
            "min_distance_bp": selector.min_distance_bp,
            "max_distance_bp": selector.max_distance_bp,
            "selection_multiplicity": selector.selection_multiplicity,
            "choose": selector.choose,
            "selection_seed": selector.selection_seed,
            "anchor": selector.anchor,
            "strand_relation": selector.strand_relation,
            "orientation": selector.orientation,
            "excluded_sites": excluded_count,
        },
        "dropped_for_window_span": dropped_for_window_span,
    }
    return _RasterSiteSelectionResult(
        window_indices=[np.asarray(values, dtype=int) for values in per_window_indices],
        panel_centers=[np.asarray(values, dtype=float) for values in per_window_centers],
        stats=stats,
    )


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


def _iter_region_batches(
    regions: str | Path | list[str | Path],
    batch_size: int,
    tmpdir: str,
):
    """Yield region-spec batches (<= ``batch_size`` regions each) for low-memory
    windowed extraction, preserving BED columns/strand.

    - a list/tuple of regions is sliced directly;
    - a ``.bed`` file is split into temporary sub-BEDs (all columns kept, so strand and
      orientation are unchanged);
    - a single region string (e.g. ``"chr1:100-200"``) yields one batch.
    """
    if isinstance(regions, (list, tuple)):
        seq = list(regions)
        for i in range(0, len(seq), batch_size):
            yield seq[i : i + batch_size]
        return
    if isinstance(regions, (str, Path)):
        path = Path(regions)
        if path.suffix == ".bed" and path.exists():
            lines = [
                ln
                for ln in path.read_text().splitlines()
                if ln.strip() and not ln.startswith("#")
            ]
            if not lines:
                yield [regions]
                return
            for i in range(0, len(lines), batch_size):
                sub = Path(tmpdir) / f"batch_{i // batch_size}.bed"
                sub.write_text("\n".join(lines[i : i + batch_size]) + "\n")
                yield sub
            return
    # single region string or unrecognized input -> one batch
    yield [regions]


def _extract_read_windows_batched(
    hdf5_file: str | Path,
    motifs: Sequence[str],
    regions: str | Path | list[str | Path],
    *,
    config: ReadWindowExtractionConfig | None,
    window_size: int | None,
    single_strand: bool,
    subset_parameters: dict | None,
    span_full_window: bool,
    region_batch_size: int,
    quiet: bool,
) -> ReadWindowExtractionResult:
    """Low-memory driver for :func:`extract_read_windows`.

    Loads reads in region batches instead of all at once, so only one batch's full-length
    read vectors are resident at a time (long genomic reads store their entire span, so
    materializing every selected read can exhaust memory). Each batch is windowed via
    ``extract_read_windows`` itself (identical windowing) and the results are concatenated.
    A shared ``window_size`` guarantees every batch produces the same window width. Row
    order may differ from the single-pass path; results are otherwise equivalent.
    """
    results: list[ReadWindowExtractionResult] = []
    with tempfile.TemporaryDirectory(prefix="dimelo_lowmem_") as tmp:
        for batch in _iter_region_batches(regions, region_batch_size, tmp):
            try:
                res = extract_read_windows(
                    hdf5_file,
                    motifs,
                    batch,
                    config=config,
                    window_size=window_size,
                    single_strand=single_strand,
                    subset_parameters=subset_parameters,
                    span_full_window=span_full_window,
                    quiet=quiet,
                    low_memory=False,
                )
            except ValueError as exc:
                # A batch may legitimately contain no reads that fully span the window;
                # skip only that case and re-raise anything else.
                if "No reads produced a full window" in str(exc):
                    continue
                raise
            results.append(res)

    if not results:
        raise ValueError("No reads produced a full window in any region batch.")

    data_matrix = np.vstack([r.data_matrix for r in results])
    val_present = all(r.val_matrix is not None for r in results)
    val_matrix = (
        np.vstack([r.val_matrix for r in results]) if val_present else None
    )
    metadata = [row for r in results for row in r.metadata]
    datasets = results[0].datasets
    regions_dict: dict = {}
    for r in results:
        if r.regions_dict:
            for chrom, entries in r.regions_dict.items():
                regions_dict.setdefault(chrom, []).extend(entries)
    return ReadWindowExtractionResult(
        data_matrix=data_matrix,
        val_matrix=val_matrix,
        metadata=metadata,
        datasets=datasets,
        regions_dict=regions_dict or None,
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
    low_memory: bool = False,
    region_batch_size: int = 64,
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
        low_memory: if True, load reads in region batches instead of all at once, so only
            one batch's full-length read vectors are resident at a time. Useful for large
            region sets and long genomic reads (whose full spans are stored), which can
            otherwise exhaust memory. Requires an explicit window_size (or config.window_size)
            so every batch shares a window width; row order may differ from the single-pass
            path but results are otherwise equivalent.
        region_batch_size: number of regions per batch when low_memory is True.

    Returns:
        ReadWindowExtractionResult containing thresholded/binary mod windows by default,
        plus val matrices and metadata.
    """

    if low_memory and regions is not None:
        if window_size is None and (config is None or config.window_size is None):
            raise ValueError(
                "low_memory=True requires an explicit window_size (or config.window_size) "
                "so all region batches produce the same window width."
            )
        return _extract_read_windows_batched(
            hdf5_file,
            motifs,
            regions,
            config=config,
            window_size=window_size,
            single_strand=single_strand,
            subset_parameters=subset_parameters,
            span_full_window=span_full_window,
            region_batch_size=region_batch_size,
            quiet=quiet,
        )

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


def _merge_feature_spec(feature_spec: dict[str, Any] | None) -> dict[str, Any]:
    spec = {
        key: (dict(value) if isinstance(value, dict) else value)
        for key, value in DEFAULT_RICH_FEATURE_SPEC.items()
    }
    if feature_spec:
        for key, value in feature_spec.items():
            if isinstance(value, dict) and isinstance(spec.get(key), dict):
                merged = dict(spec[key])
                merged.update(value)
                spec[key] = merged
            else:
                spec[key] = value
    return spec


def _feature_enabled(spec: dict[str, Any], key: str) -> bool:
    value = spec.get(key, {})
    if isinstance(value, dict):
        return bool(value.get("enabled", False))
    return bool(value)


def _window_mean(matrix: np.ndarray, center: int, start: int, end: int) -> np.ndarray:
    start_idx = max(0, center + int(start))
    end_idx = min(matrix.shape[1], center + int(end))
    if end_idx <= start_idx:
        return np.zeros(matrix.shape[0], dtype=float)
    return matrix[:, start_idx:end_idx].mean(axis=1)


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=float),
        where=np.abs(denominator) > 1e-12,
    )


def _fft_period_power(matrix: np.ndarray, period_bp: float) -> np.ndarray:
    if matrix.shape[1] < 3 or period_bp <= 0:
        return np.zeros(matrix.shape[0], dtype=float)
    centered = matrix - matrix.mean(axis=1, keepdims=True)
    power = np.abs(np.fft.rfft(centered, axis=1)) ** 2
    freqs = np.fft.rfftfreq(matrix.shape[1], d=1.0)
    target = 1.0 / period_bp
    idx = int(np.argmin(np.abs(freqs - target)))
    total = power[:, 1:].sum(axis=1)
    if idx == 0:
        return np.zeros(matrix.shape[0], dtype=float)
    return _safe_ratio(power[:, idx], total)


def _append_feature(
    columns: list[np.ndarray],
    names: list[str],
    rows: list[dict[str, Any]] | None,
    values: np.ndarray,
    name: str,
    *,
    family: str,
    motif: str,
    scope: str,
    window: str = "full",
) -> None:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        arr = arr[:, None]
    columns.append(arr)
    force_index = name == "pca" or name.endswith("__pca")
    if arr.shape[1] == 1 and not force_index:
        names.append(name)
        if rows is not None:
            rows.append(
                {
                    "feature_name": name,
                    "family": family,
                    "motif": motif,
                    "scope": scope,
                    "window": window,
                }
            )
    else:
        for i in range(arr.shape[1]):
            col_name = f"{name}_{i}"
            names.append(col_name)
            if rows is not None:
                rows.append(
                    {
                        "feature_name": col_name,
                        "family": family,
                        "motif": motif,
                        "scope": scope,
                        "window": window,
                    }
                )


def _motif_views(
    data_matrix: np.ndarray,
    val_matrix: np.ndarray | None,
    *,
    motif_count: int,
    motif_labels: Sequence[str] | None,
) -> list[tuple[str, np.ndarray, np.ndarray | None]]:
    motif_count = max(1, int(motif_count))
    if data_matrix.shape[1] % motif_count != 0:
        raise ValueError(
            "data_matrix width must be divisible by motif_count for motif-specific features."
        )
    width = data_matrix.shape[1] // motif_count
    labels = (
        list(motif_labels)
        if motif_labels is not None
        else [f"motif_{i}" for i in range(motif_count)]
    )
    if len(labels) != motif_count:
        raise ValueError("motif_labels length must match motif_count.")
    views = []
    for i, label in enumerate(labels):
        sl = slice(i * width, (i + 1) * width)
        views.append(
            (
                str(label),
                data_matrix[:, sl],
                val_matrix[:, sl] if val_matrix is not None else None,
            )
        )
    return views


def _append_rich_features_for_view(
    columns: list[np.ndarray],
    names: list[str],
    rows: list[dict[str, Any]] | None,
    matrix: np.ndarray,
    val_matrix: np.ndarray | None,
    spec: dict[str, Any],
    *,
    prefix: str,
    motif: str,
    scope: str,
    peak_prominence: float,
) -> None:
    n_reads, width = matrix.shape
    center = width // 2

    if _feature_enabled(spec, "summary"):
        q25 = np.percentile(matrix, 25, axis=1)
        q75 = np.percentile(matrix, 75, axis=1)
        for label, values in [
            ("mean", matrix.mean(axis=1)),
            ("var", matrix.var(axis=1)),
            ("median", np.median(matrix, axis=1)),
            ("q25", q25),
            ("q75", q75),
            ("iqr", q75 - q25),
        ]:
            _append_feature(
                columns,
                names,
                rows,
                values,
                f"{prefix}__{label}",
                family="summary",
                motif=motif,
                scope=scope,
            )
        if val_matrix is not None:
            valid_sum = val_matrix.sum(axis=1)
            frac = _safe_ratio(matrix.sum(axis=1), valid_sum)
            _append_feature(
                columns,
                names,
                rows,
                frac,
                f"{prefix}__mod_fraction",
                family="summary",
                motif=motif,
                scope=scope,
            )

    if _feature_enabled(spec, "densities"):
        windows = spec["densities"].get("windows", DEFAULT_DENSITY_WINDOWS)
        for label, start, end in windows:
            _append_feature(
                columns,
                names,
                rows,
                _window_mean(matrix, center, start, end),
                f"{prefix}__{label}",
                family="density",
                motif=motif,
                scope=scope,
                window=f"{start}..{end}",
            )

    if _feature_enabled(spec, "asymmetry"):
        for span in spec["asymmetry"].get("spans", (100, 300, 500)):
            left = _window_mean(matrix, center, -int(span), 0)
            right = _window_mean(matrix, center, 0, int(span))
            _append_feature(
                columns,
                names,
                rows,
                right - left,
                f"{prefix}__right_minus_left_pm{span}",
                family="asymmetry",
                motif=motif,
                scope=scope,
                window=f"-{span}..{span}",
            )

    if _feature_enabled(spec, "center_edge"):
        center_bp = int(spec["center_edge"].get("center_bp", 100))
        edge_bp = int(spec["center_edge"].get("edge_bp", 200))
        center_vals = _window_mean(matrix, center, -center_bp, center_bp)
        left_edge = _window_mean(matrix, center, -width // 2, -width // 2 + edge_bp)
        right_edge = _window_mean(matrix, center, width // 2 - edge_bp, width // 2)
        edge_vals = (left_edge + right_edge) / 2.0
        _append_feature(
            columns,
            names,
            rows,
            center_vals - edge_vals,
            f"{prefix}__center_minus_edge",
            family="center_edge",
            motif=motif,
            scope=scope,
            window=f"center_pm{center_bp}_edge{edge_bp}",
        )
        _append_feature(
            columns,
            names,
            rows,
            _safe_ratio(center_vals, edge_vals),
            f"{prefix}__center_edge_ratio",
            family="center_edge",
            motif=motif,
            scope=scope,
            window=f"center_pm{center_bp}_edge{edge_bp}",
        )

    if _feature_enabled(spec, "autocorr"):
        for lag in spec["autocorr"].get("lags", DEFAULT_AUTOCORR_LAGS):
            values = np.array([compute_autocorrelation_feature(row, int(lag)) for row in matrix])
            _append_feature(
                columns,
                names,
                rows,
                values,
                f"{prefix}__autocorr_{lag}",
                family="autocorr",
                motif=motif,
                scope=scope,
            )

    if _feature_enabled(spec, "fft"):
        for period in spec["fft"].get("periods_bp", (10, 50, 100, 150, 200)):
            _append_feature(
                columns,
                names,
                rows,
                _fft_period_power(matrix, float(period)),
                f"{prefix}__fft_power_period_{period}bp",
                family="fft",
                motif=motif,
                scope=scope,
            )

    if _feature_enabled(spec, "peaks") and find_peaks is not None:
        prominence = float(spec["peaks"].get("prominence", peak_prominence))
        peak_counts = np.zeros(n_reads, dtype=float)
        peak_prominences = np.zeros(n_reads, dtype=float)
        for row_idx, row in enumerate(matrix):
            peaks, props = find_peaks(row, prominence=prominence)
            peak_counts[row_idx] = len(peaks)
            if peaks.size > 0:
                peak_prominences[row_idx] = float(np.mean(props["prominences"]))
        _append_feature(
            columns,
            names,
            rows,
            peak_counts,
            f"{prefix}__peak_count",
            family="peaks",
            motif=motif,
            scope=scope,
        )
        _append_feature(
            columns,
            names,
            rows,
            peak_prominences,
            f"{prefix}__peak_prominence",
            family="peaks",
            motif=motif,
            scope=scope,
        )


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
    feature_spec: dict[str, Any] | None = None,
    motif_count: int | None = None,
    motif_labels: Sequence[str] | None = None,
    scale_features: bool = False,
    scaling_method: str = "standard",
    family_weighting: str | None = None,
    family_weights: dict[str, float] | None = None,
    unscaled_families: Sequence[str] | None = None,
    return_feature_table: bool = False,
    return_scale_table: bool = False,
) -> (
    tuple[np.ndarray, list[str]]
    | tuple[np.ndarray, list[str], pd.DataFrame]
    | tuple[np.ndarray, list[str], pd.DataFrame, pd.DataFrame]
):
    """
    Convert read windows into an augmented feature matrix including PCA components,
    autocorrelation, density summaries, and basic statistics.

    Set use_peak_features=True to append peak counts/prominences (requires SciPy).
    If require_nonzero_valid is True, rows with no valid sites (or below min_valid_fraction)
    are dropped before feature computation (requires val_matrix to be present).
    Set scale_features=True to return a scaled matrix suitable for distance-based
    methods such as k-means. Use return_scale_table=True to inspect the scaling.
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
    build_feature_table = return_feature_table or scale_features
    feature_rows: list[dict[str, Any]] | None = [] if build_feature_table else None

    if feature_spec is not None:
        spec = _merge_feature_spec(feature_spec)
        spec_motif_count = int(motif_count or spec.get("motif_count") or 1)
        spec_motif_labels = motif_labels or spec.get("motif_labels")
        views = _motif_views(
            data_matrix,
            val_matrix,
            motif_count=spec_motif_count,
            motif_labels=spec_motif_labels,
        )
        motif_mode = str(spec.get("motif_mode", "per_motif"))
        include_pooled = bool(spec.get("pooled", motif_mode in {"pooled", "both"}))

        pca_spec = spec.get("pca", {})
        if isinstance(pca_spec, dict) and pca_spec.get("enabled", False):
            pca_scope = str(pca_spec.get("scope", "pooled"))
            pca_n = int(pca_spec.get("n_components", n_pca))
            if pca_scope == "per_motif":
                pca_views = views
            else:
                pca_views = [("pooled", data_matrix, val_matrix)]
            for label, matrix, _view_val in pca_views:
                n_components = min(pca_n, n_reads, matrix.shape[1])
                if n_components <= 0:
                    continue
                from sklearn.decomposition import PCA

                pca = PCA(n_components=n_components, random_state=random_state)
                pca_vals = pca.fit_transform(matrix)
                _append_feature(
                    columns,
                    names,
                    feature_rows,
                    pca_vals,
                    f"{label}__pca",
                    family="pca",
                    motif=label,
                    scope=pca_scope,
                )

        if motif_mode in {"per_motif", "both"}:
            for label, matrix, view_val in views:
                _append_rich_features_for_view(
                    columns,
                    names,
                    feature_rows,
                    matrix,
                    view_val,
                    spec,
                    prefix=label,
                    motif=label,
                    scope="per_motif",
                    peak_prominence=peak_prominence,
                )

        if include_pooled:
            _append_rich_features_for_view(
                columns,
                names,
                feature_rows,
                data_matrix,
                val_matrix,
                spec,
                prefix="pooled",
                motif="pooled",
                scope="pooled",
                peak_prominence=peak_prominence,
            )

        if (
            spec_motif_count > 1
            and isinstance(spec.get("cross_motif"), dict)
            and spec["cross_motif"].get("enabled", False)
        ):
            for (left_label, left, _), (right_label, right, _) in combinations(views, 2):
                left_mean = left.mean(axis=1)
                right_mean = right.mean(axis=1)
                _append_feature(
                    columns,
                    names,
                    feature_rows,
                    left_mean - right_mean,
                    f"{left_label}_minus_{right_label}__mean",
                    family="cross_motif",
                    motif=f"{left_label}|{right_label}",
                    scope="cross_motif",
                )
                corrs = np.zeros(n_reads, dtype=float)
                for row_idx in range(n_reads):
                    if np.std(left[row_idx]) > 0 and np.std(right[row_idx]) > 0:
                        corrs[row_idx] = float(np.corrcoef(left[row_idx], right[row_idx])[0, 1])
                _append_feature(
                    columns,
                    names,
                    feature_rows,
                    corrs,
                    f"{left_label}_{right_label}__correlation",
                    family="cross_motif",
                    motif=f"{left_label}|{right_label}",
                    scope="cross_motif",
                )

        if not columns:
            raise ValueError("feature_spec produced no features.")
        feature_matrix = np.hstack(columns)
        table = pd.DataFrame(feature_rows)
        if scale_features:
            feature_matrix, scale_table = scale_feature_matrix(
                feature_matrix,
                names,
                feature_table=table,
                method=scaling_method,
                family_weighting=family_weighting,
                family_weights=family_weights,
                unscaled_families=unscaled_families,
            )
            if return_feature_table or return_scale_table:
                if return_scale_table:
                    return feature_matrix, names, table, scale_table
                return feature_matrix, names, table
        if return_feature_table:
            return feature_matrix, names, table
        return feature_matrix, names

    n_components = min(n_pca, n_reads, window_size)
    if n_components > 0:
        from sklearn.decomposition import PCA

        pca = PCA(n_components=n_components, random_state=random_state)
        pca_vals = pca.fit_transform(data_matrix)
        _append_feature(
            columns,
            names,
            feature_rows,
            pca_vals,
            "pca",
            family="pca",
            motif="pooled",
            scope="legacy_pooled",
        )

    # Precompute cumulative sums to accelerate density windows
    cumsum = np.cumsum(data_matrix, axis=1)
    cumsum = np.pad(cumsum, ((0, 0), (1, 0)), mode="constant", constant_values=0)

    for lag in autocorr_lags:
        values = np.array(
            [compute_autocorrelation_feature(row, lag) for row in data_matrix]
        )
        _append_feature(
            columns,
            names,
            feature_rows,
            values,
            f"autocorr_{lag}",
            family="autocorr",
            motif="pooled",
            scope="legacy_pooled",
        )

    for label, start, end in density_windows:
        start_idx = max(0, center + start)
        end_idx = min(window_size, center + end)
        if end_idx <= start_idx:
            values = np.zeros(n_reads)
        else:
            length = end_idx - start_idx
            window_sum = cumsum[:, end_idx] - cumsum[:, start_idx]
            values = window_sum / length
        _append_feature(
            columns,
            names,
            feature_rows,
            values,
            label,
            family="density",
            motif="pooled",
            scope="legacy_pooled",
            window=f"{start}..{end}",
        )

    global_mean = data_matrix.mean(axis=1)
    global_var = data_matrix.var(axis=1)
    global_median = np.median(data_matrix, axis=1)
    q25 = np.percentile(data_matrix, 25, axis=1)
    q75 = np.percentile(data_matrix, 75, axis=1)

    for label, values in [
        ("global_mean", global_mean),
        ("global_var", global_var),
        ("global_median", global_median),
        ("q25", q25),
        ("q75", q75),
        ("iqr", q75 - q25),
    ]:
        _append_feature(
            columns,
            names,
            feature_rows,
            values,
            label,
            family="summary",
            motif="pooled",
            scope="legacy_pooled",
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
        _append_feature(
            columns,
            names,
            feature_rows,
            frac,
            "global_mod_fraction",
            family="summary",
            motif="pooled",
            scope="legacy_pooled",
        )

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
        _append_feature(
            columns,
            names,
            feature_rows,
            np.array(peak_counts),
            "peak_count",
            family="peaks",
            motif="pooled",
            scope="legacy_pooled",
        )
        _append_feature(
            columns,
            names,
            feature_rows,
            np.array(peak_prominences),
            "peak_prominence",
            family="peaks",
            motif="pooled",
            scope="legacy_pooled",
        )

    feature_matrix = np.hstack(columns)
    table = pd.DataFrame(feature_rows)
    if scale_features:
        feature_matrix, scale_table = scale_feature_matrix(
            feature_matrix,
            names,
            feature_table=table,
            method=scaling_method,
            family_weighting=family_weighting,
            family_weights=family_weights,
            unscaled_families=unscaled_families,
        )
        if return_scale_table:
            return feature_matrix, names, table, scale_table
    if return_feature_table:
        return feature_matrix, names, table
    return feature_matrix, names


def summarize_feature_matrix(
    feature_matrix: np.ndarray,
    feature_names: Sequence[str],
    *,
    feature_table: pd.DataFrame | None = None,
    labels: Sequence[Any] | None = None,
    top_n_variable: int = 12,
) -> dict[str, Any]:
    """Summarize feature matrix dimensions, missingness, variance, and optional labels."""

    X = np.asarray(feature_matrix, dtype=float)
    if X.ndim != 2:
        raise ValueError("feature_matrix must be 2-dimensional.")
    if X.shape[1] != len(feature_names):
        raise ValueError("feature_names length must match feature_matrix width.")

    variance = np.nanvar(X, axis=0)
    missing_fraction = np.isnan(X).mean(axis=0)
    order = np.argsort(variance)[::-1][:top_n_variable]
    variable = pd.DataFrame(
        {
            "feature_name": [feature_names[i] for i in order],
            "variance": variance[order],
            "missing_fraction": missing_fraction[order],
        }
    )
    summary: dict[str, Any] = {
        "n_reads": int(X.shape[0]),
        "n_features": int(X.shape[1]),
        "missing_values": int(np.isnan(X).sum()),
        "zero_variance_features": int(np.sum(variance <= 1e-12)),
        "top_variable_features": variable,
    }
    if feature_table is not None and not feature_table.empty:
        summary["feature_counts_by_family"] = (
            feature_table.groupby("family", dropna=False)
            .size()
            .rename("n_features")
            .reset_index()
            .sort_values("n_features", ascending=False)
        )
        summary["feature_counts_by_motif"] = (
            feature_table.groupby("motif", dropna=False)
            .size()
            .rename("n_features")
            .reset_index()
            .sort_values("n_features", ascending=False)
        )
    if labels is not None:
        summary["label_counts"] = pd.Series(labels).value_counts(dropna=False).rename("n")
    return summary


def scale_feature_matrix(
    feature_matrix: np.ndarray,
    feature_names: Sequence[str],
    *,
    feature_table: pd.DataFrame | None = None,
    method: str = "standard",
    family_weighting: str | None = None,
    family_weights: dict[str, float] | None = None,
    unscaled_families: Sequence[str] | None = None,
    return_scaler: bool = False,
) -> tuple[np.ndarray, pd.DataFrame] | tuple[np.ndarray, pd.DataFrame, Any]:
    """
    Scale engineered features for distance-based methods.

    K-means, Ward agglomerative clustering, and Euclidean nearest-neighbor views are
    sensitive to feature scale. The default "standard" method centers each feature
    and scales it to unit variance. "robust" uses median/IQR scaling for heavier tails.
    If feature_table is provided, family_weighting="equal_family" downweights each
    column by sqrt(number of scaled columns in its family) so families with many
    columns do not dominate Euclidean distance solely by column count.
    """

    X = np.asarray(feature_matrix, dtype=float)
    if X.ndim != 2:
        raise ValueError("feature_matrix must be 2-dimensional.")
    if X.shape[1] != len(feature_names):
        raise ValueError("feature_names length must match feature_matrix width.")

    method_norm = method.lower()
    unscaled_family_set = set(unscaled_families or ())
    families = pd.Series(["feature"] * X.shape[1], dtype=object)
    if feature_table is not None and not feature_table.empty:
        if "feature_name" not in feature_table.columns or "family" not in feature_table.columns:
            raise ValueError("feature_table must include feature_name and family columns.")
        family_map = feature_table.set_index("feature_name")["family"].to_dict()
        families = pd.Series(
            [family_map.get(name, "feature") for name in feature_names],
            dtype=object,
        )
    scaled_mask = ~families.isin(unscaled_family_set).to_numpy()

    if method_norm == "standard":
        from sklearn.preprocessing import StandardScaler

        scaler = StandardScaler()
    elif method_norm == "robust":
        from sklearn.preprocessing import RobustScaler

        scaler = RobustScaler(quantile_range=(25.0, 75.0))
    elif method_norm in {"none", "identity"}:
        scaled = X.copy()
        table = pd.DataFrame(
            {
                "feature_name": list(feature_names),
                "scaling_method": "none",
                "family": families.to_numpy(),
                "family_weight": np.ones(X.shape[1], dtype=float),
                "center": np.zeros(X.shape[1], dtype=float),
                "scale": np.ones(X.shape[1], dtype=float),
                "raw_mean": np.nanmean(X, axis=0),
                "raw_std": np.nanstd(X, axis=0),
            }
        )
        if return_scaler:
            return scaled, table, None
        return scaled, table
    else:
        raise ValueError("method must be 'standard', 'robust', or 'none'.")

    scaled = X.copy()
    center = np.zeros(X.shape[1], dtype=float)
    scale = np.ones(X.shape[1], dtype=float)
    if np.any(scaled_mask):
        scaled[:, scaled_mask] = scaler.fit_transform(X[:, scaled_mask])
        center_values = getattr(
            scaler, "mean_", getattr(scaler, "center_", np.zeros(np.sum(scaled_mask)))
        )
        scale_values = getattr(scaler, "scale_", np.ones(np.sum(scaled_mask)))
        center[scaled_mask] = center_values
        scale[scaled_mask] = scale_values

    weights = np.ones(X.shape[1], dtype=float)
    if family_weighting is not None:
        weighting = family_weighting.lower()
        if weighting != "equal_family":
            raise ValueError("family_weighting must be None or 'equal_family'.")
        for _family, idx in families.groupby(families).groups.items():
            idx_arr = np.array(list(idx), dtype=int)
            scaled_idx = idx_arr[scaled_mask[idx_arr]]
            if len(scaled_idx) > 0:
                weights[scaled_idx] *= 1.0 / math.sqrt(len(scaled_idx))
    if family_weights:
        for family, weight in family_weights.items():
            weights[families.to_numpy() == family] *= float(weight)
    scaled *= weights

    table = pd.DataFrame(
        {
            "feature_name": list(feature_names),
            "scaling_method": method_norm,
            "family": families.to_numpy(),
            "family_weight": weights,
            "scaled": scaled_mask,
            "center": center,
            "scale": scale,
            "raw_mean": np.nanmean(X, axis=0),
            "raw_std": np.nanstd(X, axis=0),
            "scaled_mean": np.nanmean(scaled, axis=0),
            "scaled_std": np.nanstd(scaled, axis=0),
        }
    )
    if return_scaler:
        return scaled, table, scaler
    return scaled, table


def rank_read_features_by_group_difference(
    feature_matrix: np.ndarray,
    feature_names: Sequence[str],
    labels: Sequence[Any],
    *,
    top_n: int = 20,
) -> pd.DataFrame:
    """Rank features by standardized mean separation across two or more labels."""

    X = np.asarray(feature_matrix, dtype=float)
    y = pd.Series(labels)
    if X.ndim != 2:
        raise ValueError("feature_matrix must be 2-dimensional.")
    if X.shape[0] != len(y):
        raise ValueError("labels length must match feature_matrix rows.")
    if X.shape[1] != len(feature_names):
        raise ValueError("feature_names length must match feature_matrix width.")

    groups = list(pd.unique(y))
    if len(groups) < 2:
        raise ValueError("At least two label groups are required.")

    rows = []
    overall_std = np.nanstd(X, axis=0)
    for idx, name in enumerate(feature_names):
        means = []
        for group in groups:
            means.append(float(np.nanmean(X[y.to_numpy() == group, idx])))
        effect = (max(means) - min(means)) / (overall_std[idx] + 1e-12)
        rows.append(
            {
                "feature_name": name,
                "effect_score": float(effect),
                "min_group_mean": min(means),
                "max_group_mean": max(means),
                "overall_std": float(overall_std[idx]),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values("effect_score", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )


def raster_stats_brief(stats: dict[str, Any]) -> dict[str, Any]:
    """Return the raster stats fields most useful for notebook display."""

    keys = [
        "pairs",
        "rows_are",
        "unique_reads",
        "site_sets",
        "rows_before_downsample",
        "rows_after_downsample",
        "downsampled",
        "downsample_method",
        "coordinate_mode",
        "selection_mode",
        "window_offsets_bp",
        "window_widths_bp",
        "ml_score_thresholds",
    ]
    brief = {key: stats.get(key) for key in keys if key in stats}
    brief["observed_offsets_bp"] = stats.get("observed_window_center_offsets_summary_bp")
    brief["site_selection"] = stats.get("site_selection")
    return brief


def infer_region_source_labels(
    region_axis_table: pd.DataFrame,
    read_metadata: Sequence[dict[str, Any]],
    *,
    source_field: str = "source_label",
    unknown_label: str = "unknown",
) -> pd.DataFrame:
    """Attach the most frequent read source label to each region-axis row."""

    axis_table = region_axis_table.copy()
    if "region_id" not in axis_table.columns:
        raise ValueError("region_axis_table must include a region_id column.")
    meta_df = pd.DataFrame(read_metadata).copy()
    if meta_df.empty or not {"chromosome", "region_start", "region_end"}.issubset(
        meta_df.columns
    ):
        axis_table[source_field] = unknown_label
        return axis_table

    if "region_strand" not in meta_df.columns:
        meta_df["region_strand"] = "."
    if source_field not in meta_df.columns:
        meta_df[source_field] = unknown_label
    meta_df["region_id"] = meta_df.apply(
        lambda row: (
            f"{row['chromosome']}:{int(row['region_start'])}-"
            f"{int(row['region_end'])}:{row.get('region_strand', '.')}"
        ),
        axis=1,
    )
    source_by_region = (
        meta_df.groupby(["region_id", source_field], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["region_id", "n"], ascending=[True, False])
        .drop_duplicates(subset=["region_id"], keep="first")
        .loc[:, ["region_id", source_field]]
    )
    axis_table = axis_table.merge(source_by_region, on="region_id", how="left")
    axis_table[source_field] = axis_table[source_field].fillna(unknown_label)
    return axis_table


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
    scale: bool = False,
    **kwargs,
) -> ClusterResult:
    """
    Run clustering on a read feature matrix with support for several algorithms.
    auto_k=True grid-searches K for methods that require it using the requested score.

    scale=True applies a plain per-column StandardScaler to ``feature_matrix`` before
    clustering. This matters for the distance-based methods (kmeans, agglomerative,
    spectral, umap_kmeans): without it a high-variance feature family (e.g. the raw
    PCA/positional block) dominates the Euclidean distance and drowns out lower-variance
    families such as autocorrelation, so the partition collapses onto overall
    methylation level. For family-aware weighting use :func:`scale_feature_matrix`
    beforehand instead.
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
    if scale:
        from sklearn.preprocessing import StandardScaler

        X = StandardScaler().fit_transform(X)

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


def select_k_by_stability(
    feature_matrix: np.ndarray,
    *,
    k_grid: Sequence[int] | None = None,
    n_boot: int = 10,
    sample_frac: float = 0.8,
    stability_thresh: float = 0.7,
    margin: float = 0.1,
    random_state: int = 42,
    default: int | None = None,
) -> dict[str, Any]:
    """Choose the number of clusters by *reproducibility* rather than a single
    internal score.

    For each candidate k the full matrix is clustered (k-means), then re-clustered on
    many random ``sample_frac`` subsamples; each subsample partition is scored against
    the full-data labels with the adjusted Rand index (ARI). A k whose partition is
    recovered under resampling (mean ARI >= ``stability_thresh``) reflects real
    structure; one that is not is an artifact of the particular sample. Among the k
    that clear the threshold we keep those within ``margin`` of the peak stability and
    return the LARGEST (resolve as many robust states as the data supports, without
    over-splitting when a coarser partition is far more reproducible). Falls back to
    the most stable k, then ``default`` / the smallest grid value.

    Complements ``cluster_read_windows(auto_k=True, score="silhouette")``, which
    favours few well-separated blobs and tends to under-split biological state data.

    Returns a dict with ``best_k`` and a per-k ``stats`` list of
    ``{"k", "stability", "silhouette"}``.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score, silhouette_score

    X = np.asarray(feature_matrix)
    n = X.shape[0]
    grid = [int(k) for k in (k_grid or range(2, 11)) if 2 <= k < n]
    rng = np.random.RandomState(random_state)
    stats: list[dict[str, float]] = []
    for k in grid:
        base = KMeans(n_clusters=k, random_state=random_state, n_init=10).fit_predict(X)
        if len(np.unique(base)) < 2:
            continue
        aris = []
        for _ in range(n_boot):
            idx = rng.choice(n, int(sample_frac * n), replace=False)
            lab = KMeans(
                n_clusters=k, random_state=rng.randint(0, 1_000_000), n_init=10
            ).fit_predict(X[idx])
            aris.append(adjusted_rand_score(base[idx], lab))
        sil = float("nan")
        with contextlib.suppress(Exception):
            sil = float(silhouette_score(X, base))
        stats.append(
            {"k": int(k), "stability": float(np.mean(aris)), "silhouette": sil}
        )
    if not stats:
        fallback = default if default is not None else (grid[0] if grid else 2)
        return {"best_k": int(fallback), "stats": stats}
    stable = [s for s in stats if s["stability"] >= stability_thresh]
    if stable:
        peak = max(s["stability"] for s in stats)
        near = [s for s in stable if s["stability"] >= peak - margin]
        best = max(near, key=lambda s: s["k"])["k"]
    else:
        best = max(stats, key=lambda s: s["stability"])["k"]
    return {"best_k": int(best), "stats": stats}


def bernoulli_mixture(
    binary_matrix: np.ndarray,
    n_components: int,
    *,
    n_iter: int = 100,
    tol: float = 1e-5,
    prior_a: float = 2.0,
    prior_b: float = 2.0,
    dirichlet_conc: float = 1.05,
    random_state: int = 42,
) -> dict[str, Any]:
    """Bayesian (MAP) EM for a mixture of ``n_components`` multivariate Bernoulli
    distributions on a binary read-window matrix ``binary_matrix`` [n_reads x n_pos].

    This is the generative model matching per-position binary methylation calls:
    component c carries a profile ``mu[c]`` in (0, 1) giving P(methylated) at each
    position with mixing weight ``pi[c]``; a read is n_pos independent Bernoulli draws
    at biases ``mu[c]``. Conjugate priors make the fit Bayesian and numerically stable:

    * a per-position ``Beta(prior_a, prior_b)`` prior on each ``mu`` regularizes
      positions seen in few reads toward the prior mean instead of collapsing to 0/1
      (``a = b = 2`` is a weak "expect ~0.5, worth one pseudo-observation each way"
      prior = principled Laplace smoothing);
    * a symmetric ``Dirichlet(dirichlet_conc)`` prior on the weights shrinks empty
      components so a slightly-too-large ``n_components`` degrades gracefully.

    EM ascends the log-posterior: the E-step computes responsibilities in log-space
    (log-sum-exp); the M-step uses the posterior modes
    ``mu = (S + a - 1) / (N + a + b - 2)`` and
    ``pi = (N + conc - 1) / (n + k (conc - 1))``. The
    ``log p(x|c) = const_c + x . (log mu - log(1 - mu))`` identity avoids materializing
    ``1 - X``.

    Unlike distance-based clustering the assignment is soft (``resp``) and the
    component means ``mu`` are directly the per-cluster occupancy metaprofiles.
    Returns dict with ``labels`` (hard argmax), ``resp`` [n x k], ``mu`` [k x n_pos],
    ``pi`` [k] and the ``logpost`` trajectory.
    """
    rng = np.random.RandomState(random_state)
    X = np.asarray(binary_matrix, dtype=np.float64)
    n, d = X.shape
    k = min(int(n_components), n)
    eps = 1e-6
    mu = np.clip(X[rng.choice(n, k, replace=False)], eps, 1 - eps)
    pi = np.full(k, 1.0 / k)
    log_resp = np.empty((n, k), dtype=np.float64)
    logpost: list[float] = []
    prev = -np.inf
    resp = np.full((n, k), 1.0 / k)
    for _ in range(int(n_iter)):
        log_mu = np.log(mu)
        log_1m = np.log1p(-mu)
        for c in range(k):
            log_resp[:, c] = (
                math.log(pi[c]) + log_1m[c].sum() + X @ (log_mu[c] - log_1m[c])
            )
        row_max = log_resp.max(axis=1, keepdims=True)
        lse = row_max[:, 0] + np.log(np.exp(log_resp - row_max).sum(axis=1))
        data_ll = float(lse.sum())
        log_prior = float(
            ((prior_a - 1) * log_mu + (prior_b - 1) * log_1m).sum()
            + (dirichlet_conc - 1) * np.log(pi).sum()
        )
        current = data_ll + log_prior
        resp = np.exp(log_resp - lse[:, None])
        nk = resp.sum(axis=0)
        s = resp.T @ X
        mu = np.clip(
            (s + prior_a - 1) / (nk[:, None] + prior_a + prior_b - 2), eps, 1 - eps
        )
        pi = (nk + dirichlet_conc - 1) / (n + k * (dirichlet_conc - 1))
        logpost.append(current)
        if current - prev < tol * abs(prev if prev != 0 else 1.0):
            break
        prev = current
    return {
        "labels": resp.argmax(axis=1),
        "resp": resp,
        "mu": mu,
        "pi": pi,
        "logpost": logpost,
    }


def _n_states(labels: np.ndarray) -> int:
    """Number of non-noise clusters in a label vector (noise == -1 is excluded)."""
    ref = np.asarray(labels)
    non_noise = ref[ref >= 0]
    return int(len(np.unique(non_noise if non_noise.size else ref)))


def crossvalidate_clusters(
    feature_matrix: np.ndarray,
    labels: np.ndarray,
    *,
    covariance_type: str = "diag",
    max_extra_components: int = 3,
    n_init: int = 3,
    random_state: int = 42,
) -> dict[str, Any]:
    """Bayesian cross-check of a hard partition on the SAME feature matrix.

    Fits a variational ``BayesianGaussianMixture`` with a Dirichlet-PROCESS weight
    prior, seeded with a few more components than ``labels`` has, and lets variational
    inference drive unsupported components' weights toward ~0 (automatic relevance
    determination) — a Bayesian read on how many states the data supports. A high
    adjusted Rand index versus ``labels`` means the partition is model-agnostic rather
    than an artifact of the original algorithm's assumptions.

    Returns dict: ``ari`` vs ``labels``, ``n_effective`` (components with weight > 1%),
    ``mean_confidence`` (average top posterior; low => many ambiguous reads),
    ``weights``, per-read ``posteriors`` [n x k] and hard ``labels`` (row-aligned).
    """
    from sklearn.metrics import adjusted_rand_score
    from sklearn.mixture import BayesianGaussianMixture

    X = np.asarray(feature_matrix)
    ref = np.asarray(labels)
    k = _n_states(ref)
    if k < 2 or X.shape[0] <= k:
        return {
            "ari": float("nan"),
            "n_effective": 0,
            "mean_confidence": float("nan"),
            "weights": None,
            "posteriors": None,
            "labels": None,
        }
    upper = min(X.shape[0] - 1, k + max(0, int(max_extra_components)))
    model = BayesianGaussianMixture(
        n_components=upper,
        covariance_type=covariance_type,
        weight_concentration_prior_type="dirichlet_process",
        weight_concentration_prior=1.0 / upper,
        n_init=n_init,
        max_iter=300,
        random_state=random_state,
    )
    pred = model.fit_predict(X)
    post = model.predict_proba(X)
    return {
        "ari": float(adjusted_rand_score(ref, pred)),
        "n_effective": int(np.sum(model.weights_ > 0.01)),
        "mean_confidence": float(post.max(axis=1).mean()),
        "weights": np.asarray(model.weights_),
        "posteriors": post,
        "labels": pred,
    }


def bernoulli_crosscheck(
    binary_matrix: np.ndarray,
    labels: np.ndarray,
    *,
    cap: int = 4000,
    random_state: int = 42,
) -> dict[str, Any]:
    """Cross-check a hard partition against a Bayesian Bernoulli mixture on the RAW
    binary windows (its correct likelihood), complementing
    :func:`crossvalidate_clusters` which works in the engineered feature space.

    Fits :func:`bernoulli_mixture` with the same number of states as ``labels`` on a
    random subsample (<= ``cap`` reads) and reports agreement (adjusted Rand index).
    Returns dict: ``ari``, ``mean_confidence``, ``mu`` (state profiles) and ``n_used``.
    """
    from sklearn.metrics import adjusted_rand_score

    X = np.asarray(binary_matrix)
    ref = np.asarray(labels)
    k = _n_states(ref)
    n = X.shape[0]
    if k < 2 or n <= k:
        return {
            "ari": float("nan"),
            "mean_confidence": float("nan"),
            "mu": None,
            "n_used": 0,
        }
    rng = np.random.RandomState(random_state)
    idx = rng.choice(n, cap, replace=False) if n > cap else np.arange(n)
    fit = bernoulli_mixture(X[idx], k, random_state=random_state)
    return {
        "ari": float(adjusted_rand_score(ref[idx], fit["labels"])),
        "mean_confidence": float(fit["resp"].max(axis=1).mean()),
        "mu": fit["mu"],
        "n_used": int(len(idx)),
    }


def plot_cluster_validation(
    labels: np.ndarray,
    *,
    stability: dict[str, Any] | None = None,
    gmm: dict[str, Any] | None = None,
    bernoulli: dict[str, Any] | None = None,
    best_k: int | None = None,
    title: str = "",
    cmap_name: str = "viridis",
    smooth_window: int = 35,
    span_bp: int | None = None,
):
    """Single-page "is this clustering real?" validation dashboard (six panels):

    (a) k-selection    — stability (bootstrap ARI) and silhouette vs k, chosen k marked;
    (b) cross-model    — ARI of the partition vs the Bayesian GMM (features) and vs the
                         Bernoulli mixture (raw binary), with the 0.7 reference line;
    (c) GMM weights    — Dirichlet-process posterior component weights (bars below the
                         1% line are auto-pruned); title reports n_effective;
    (d) BMM profiles   — Bernoulli-mixture component means (per-state P(methylated));
    (e) confidence     — histogram of the top soft posterior per read (ambiguity);
    (f) overlap        — partition x GMM-component contingency (row-normalized).

    Inputs are the dicts from :func:`select_k_by_stability` (``stability``),
    :func:`crossvalidate_clusters` (``gmm``) and :func:`bernoulli_crosscheck`
    (``bernoulli``); any missing input yields an "n/a" placeholder. Returns the Figure.
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec

    ref = np.asarray(labels)
    n_states = _n_states(ref)
    cmap = plt.get_cmap(cmap_name)
    fig = plt.figure(figsize=(16, 8.5))
    grid = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.32)

    def _na(ax, msg="n/a"):
        ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=11)
        ax.set_axis_off()

    # (a) k-selection
    ax = fig.add_subplot(grid[0, 0])
    if stability and stability.get("stats"):
        stats = stability["stats"]
        ks = [s["k"] for s in stats]
        ax.plot(ks, [s["stability"] for s in stats], "o-", color="C0")
        ax.set_ylim(0, 1.03)
        ax.set_ylabel("stability (boot ARI)", color="C0")
        ax.set_xlabel("k")
        ax2 = ax.twinx()
        ax2.plot(ks, [s["silhouette"] for s in stats], "s--", color="C3")
        ax2.set_ylabel("silhouette", color="C3")
        chosen = best_k if best_k is not None else stability.get("best_k")
        if chosen is not None:
            ax.axvline(chosen, color="k", ls=":", lw=1.3)
        ax.set_title(f"k-selection (chosen k={chosen})", fontsize=10)
    else:
        _na(ax, "no k-sweep")

    # (b) cross-model ARI
    ax = fig.add_subplot(grid[0, 1])
    names, vals = [], []
    if gmm and gmm.get("ari") == gmm.get("ari"):
        names.append("Bayes GMM\n(features)")
        vals.append(gmm["ari"])
    if bernoulli and bernoulli.get("ari") == bernoulli.get("ari"):
        names.append("Bernoulli MM\n(raw binary)")
        vals.append(bernoulli["ari"])
    if vals:
        bars = ax.bar(names, vals, color=["#4c72b0", "#55a868"][: len(vals)])
        ax.set_ylim(0, 1.06)
        ax.axhline(0.7, color="grey", ls="--", lw=0.8)
        for bar, val in zip(bars, vals, strict=False):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                val + 0.02,
                f"{val:.2f}",
                ha="center",
                fontsize=9,
            )
        ax.set_ylabel("ARI vs partition")
        ax.set_title("cross-model agreement", fontsize=10)
    else:
        _na(ax)

    # (c) Bayesian-GMM Dirichlet-process weights
    ax = fig.add_subplot(grid[0, 2])
    if gmm and gmm.get("weights") is not None:
        weights = np.sort(np.asarray(gmm["weights"]))[::-1]
        ax.bar(range(len(weights)), weights, color="#4c72b0")
        ax.axhline(0.01, color="r", ls="--", lw=0.8)
        ax.set_xlabel("component (sorted)")
        ax.set_ylabel("posterior weight")
        ax.set_title(f"Bayes GMM DP prior: n_eff={gmm.get('n_effective')}", fontsize=10)
    else:
        _na(ax)

    # (d) Bernoulli-mixture state profiles
    ax = fig.add_subplot(grid[1, 0])
    if bernoulli and bernoulli.get("mu") is not None:
        mu = np.asarray(bernoulli["mu"])
        d = mu.shape[1]
        x_axis = _centered_x_axis(d, span_bp)
        for c in range(mu.shape[0]):
            ax.plot(
                x_axis,
                _smooth_profile_vector(
                    mu[c], smoothing="gaussian", smooth_win=smooth_window
                ),
                lw=1.3,
                color=cmap(c / max(1, mu.shape[0] - 1)),
                label=f"S{c}",
            )
        ax.axvline(0, color="k", lw=0.6, alpha=0.4)
        ax.set_xlabel("position (bp from center)" if span_bp else "position")
        ax.set_ylabel("P(methylated)")
        ax.set_title("Bernoulli-MM state profiles (μ)", fontsize=10)
        ax.legend(fontsize=7, ncol=2, framealpha=0.4)
    else:
        _na(ax)

    # (e) soft-assignment confidence
    ax = fig.add_subplot(grid[1, 1])
    if gmm and gmm.get("posteriors") is not None:
        conf = np.asarray(gmm["posteriors"]).max(axis=1)
        ax.hist(
            conf,
            bins=20,
            range=(1.0 / max(2, n_states), 1.0),
            color="#8172b3",
            alpha=0.85,
        )
        ax.axvline(np.median(conf), color="k", ls=":", lw=1)
        ax.set_xlabel("top posterior per read")
        ax.set_ylabel("reads")
        ax.set_title(
            f"soft-assignment confidence (median {np.median(conf):.2f})", fontsize=10
        )
    else:
        _na(ax)

    # (f) partition x GMM-component overlap
    ax = fig.add_subplot(grid[1, 2])
    if gmm and gmm.get("labels") is not None:
        gmm_labels = np.asarray(gmm["labels"])
        ku = np.sort(np.unique(ref))
        gu = np.sort(np.unique(gmm_labels))
        contingency = np.zeros((len(ku), len(gu)))
        for i, a in enumerate(ku):
            for j, b in enumerate(gu):
                contingency[i, j] = np.sum((ref == a) & (gmm_labels == b))
        row_norm = contingency / contingency.sum(axis=1, keepdims=True).clip(min=1)
        im = ax.imshow(row_norm, aspect="auto", cmap="magma", vmin=0, vmax=1)
        ax.set_xticks(range(len(gu)))
        ax.set_xticklabels(gu, fontsize=7)
        ax.set_yticks(range(len(ku)))
        ax.set_yticklabels(ku, fontsize=7)
        ax.set_xlabel("Bayes-GMM component")
        ax.set_ylabel("partition cluster")
        ax.set_title("partition → GMM overlap (row-norm)", fontsize=10)
        fig.colorbar(im, ax=ax, fraction=0.046)
    else:
        _na(ax)

    fig.suptitle(f"{title} — clustering validation", fontsize=12)
    fig.subplots_adjust(left=0.055, right=0.97, bottom=0.08, top=0.92)
    return fig


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
                color=motif_colors[motif_id],
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
    selection_mode: str = "cooccurring",
    window_offsets_bp: Sequence[int] | None = None,
    primary_window_index: int = 0,
    coordinate_mode: str | None = None,
    site_selection: dict[str, Any] | None = None,
    distance_mode: str = "center",  # or "bounds"
    window_size: int | None = None,
    motif_index: int = 0,
    motif_count: int | None = None,
    plot_all_motifs: bool = False,
    motif_labels: Sequence[str] | None = None,
    window_widths_bp: int | Sequence[int] | None = None,
    enforce_full_window_span: bool = True,
    min_read_length_bp: int | None = None,
    smoothing: str | None = "gaussian",  # None, "boxcar", "gaussian"
    smooth_win: int = 21,
    smooth_sigma_bp: float = 6.0,
    max_rows: int | None = 500,
    downsample_method: str = "auto",  # "auto" | "bin_mean" | "uniform" | "random"
    downsample_seed: int | None = None,
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
        n_windows: number of windows per read for co-occurring-region mode.
        min_separation_bp: minimum separation between selected site centers.
        selection_mode: shorthand selector mode. "fixed_offsets" requires exact
            window_offsets_bp matches per read; "cooccurring" chooses regions from the
            same read separated by min_separation_bp. The legacy alias
            "cooccurring_regions" is also accepted.
        window_offsets_bp: fixed window offsets relative to the primary window center.
            Pass one value per window; the primary window offset should usually be 0.
        primary_window_index: window used as the reference coordinate frame and default sort window.
        coordinate_mode: "relative_to_primary" preserves offsets along each read; "local_window"
            centers every panel on its own selected site.
        site_selection: optional rich selector dict. Supports modes "cooccurring",
            "fixed_offsets", and "anchor_plus_neighbors"; distance bounds, seeded
            random choice, strand filters, exclusions, and orientation by anchor strand.
            When provided, it supersedes selection_mode/window_offsets_bp/n_windows
            shorthand for site selection.
        distance_mode: "center" uses region_start/end center; "bounds" uses read_start/read_end
        window_size: half-window in bp; if provided, per-motif slice width is ``2 * window_size``.
        motif_index: reference motif slice index for single-motif plotting and ordering
        motif_count: total motif slices concatenated in data_matrix; inferred when possible
        plot_all_motifs: when True, render all motif slices in a motif x window panel grid
        motif_labels: optional motif names aligned to motif slices
        window_widths_bp: per-window display widths in bp (common width when scalar). Defaults
            to the extracted span for each panel.
        enforce_full_window_span: when True (default), drop read-sets that cannot span the
            full displayed windows based on read_start/read_end (or read_length when available).
        min_read_length_bp: optional minimum read length filter. Must be at least the required
            span implied by requested centers/separation plus plotted window width.
        smoothing: None, "boxcar", or "gaussian"
        smooth_win: smoothing window length
        smooth_sigma_bp: sigma for gaussian smoothing
        max_rows: cap rows; when exceeded, downsample according to ``downsample_method``.
            Pass None to disable downsampling.
        downsample_method: "auto" uses "uniform" for scatter plots to preserve real reads
            and "bin_mean" for heatmaps; "bin_mean" averages adjacent sorted read rows;
            "uniform" takes evenly spaced sorted rows; "random" samples sorted rows.
        downsample_seed: seed for reproducible random downsampling.
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
    if selection_mode not in {"fixed_offsets", "cooccurring_regions", "cooccurring"}:
        raise ValueError(
            "selection_mode must be 'fixed_offsets', 'cooccurring_regions', or 'cooccurring'."
        )
    if coordinate_mode is not None and coordinate_mode not in {
        "relative_to_primary",
        "local_window",
    }:
        raise ValueError(
            "coordinate_mode must be None, 'relative_to_primary', or 'local_window'."
        )
    if distance_mode not in {"center", "bounds"}:
        raise ValueError("distance_mode must be 'center' or 'bounds'.")
    if n_windows < 1:
        raise ValueError("n_windows must be >= 1.")
    if min_separation_bp < 0:
        raise ValueError("min_separation_bp must be >= 0.")
    if downsample_method not in {"auto", "bin_mean", "uniform", "random"}:
        raise ValueError(
            "downsample_method must be 'auto', 'bin_mean', 'uniform', or 'random'."
        )
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

    def _normalize_ml_scores(matrix: np.ndarray) -> np.ndarray:
        values = np.asarray(matrix, dtype=float)
        finite = values[np.isfinite(values)]
        if finite.size and float(np.nanmax(finite)) > 1.0:
            values = values / 255.0
        return np.clip(values, 0.0, 1.0)

    ML_by_motif = [_normalize_ml_scores(matrix) for matrix in X_by_motif]

    raw_selector_spec = dict(site_selection or {})
    raw_mode = raw_selector_spec.get("mode", selection_mode)
    raw_offsets = raw_selector_spec.get("window_offsets_bp", window_offsets_bp)
    if raw_mode in {"fixed_offsets"}:
        if raw_offsets is None:
            raise ValueError(
                "window_offsets_bp is required when selection_mode='fixed_offsets'."
            )
        n_windows_effective = len(list(raw_offsets))
    else:
        n_windows_effective = int(raw_selector_spec.get("n_windows", n_windows))

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

    selector = _resolve_raster_site_selection(
        site_selection,
        selection_mode=selection_mode,
        window_offsets_bp=window_offsets_bp,
        n_windows=n_windows,
        min_separation_bp=min_separation_bp,
        primary_window_index=primary_window_index,
        panel_widths_bp=panel_widths_bp,
    )
    n_windows_effective = selector.n_windows
    primary_window_index = selector.primary_window_index
    target_offsets = (
        np.asarray(selector.window_offsets_bp, dtype=float)
        if selector.window_offsets_bp is not None
        else None
    )
    effective_coordinate_mode = (
        coordinate_mode
        if coordinate_mode is not None
        else ("relative_to_primary" if selector.mode == "fixed_offsets" else "local_window")
    )

    if target_offsets is None:
        required_read_length_bp = int(
            max(
                slice_width,
                ((n_windows_effective - 1) * selector.min_distance_bp) + slice_width,
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
        if selector.mode == "fixed_offsets" and min_read_length_bp is None
        else min_read_length_bp
    )

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

        if unresolved_present and candidate_lengths:
            required_span = max(right_edges) - min(left_edges)
            if max(candidate_lengths) < required_span:
                return False
        return True

    selection_result = _select_raster_site_windows(
        meta,
        original_indices=filtered_original_indices,
        selector=selector,
        distance_mode=distance_mode,
        selection_spans_all_windows=_selection_spans_all_windows,
    )
    window_indices = selection_result.window_indices
    panel_centers = selection_result.panel_centers
    dropped_for_window_span = int(selection_result.stats["dropped_for_window_span"])
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
        sort_window_pos = int(primary_window_index)
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
    primary_centers = centers_sorted[int(primary_window_index)]
    if effective_coordinate_mode == "local_window":
        center_offsets = [
            np.zeros_like(center_values, dtype=float)
            for center_values in centers_sorted
        ]
    else:
        center_offsets = [
            center_values - primary_centers for center_values in centers_sorted
        ]
    if render_mode == "heatmap" and effective_coordinate_mode == "relative_to_primary":
        for offsets in center_offsets:
            if offsets.size > 1 and not np.allclose(
                offsets, offsets[0], equal_nan=True
            ):
                raise ValueError(
                    "Heatmap rendering with coordinate_mode='relative_to_primary' requires constant "
                    "per-read offsets within each panel because imshow uses one rectangular extent. "
                    "Use render_mode='scatter', coordinate_mode='local_window', or fixed offsets "
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
    effective_downsample_method = (
        ("uniform" if render_mode == "scatter" else "bin_mean")
        if downsample_method == "auto"
        else downsample_method
    )
    downsampled = False
    downsample_factor = 1
    if max_rows is not None and max_rows < P:
        step = math.ceil(P / max_rows)
        downsample_factor = step
        if effective_downsample_method == "bin_mean":
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
        elif effective_downsample_method == "random":
            rng = np.random.default_rng(downsample_seed)
            keep = np.sort(rng.choice(P, size=max_rows, replace=False))
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

    if effective_coordinate_mode == "relative_to_primary":
        window_plot_order = sorted(
            range(n_windows_effective),
            key=lambda idx: float(np.nanmedian(center_offsets[idx]))
            if center_offsets[idx].size
            else 0.0,
            reverse=bool(rotate and effective_axis_orientation == "position_y"),
        )
    else:
        window_plot_order = list(range(n_windows_effective))

    share_x_axes = True
    if effective_axis_orientation == "position_x" and n_windows_effective > 1:
        # When each window has a different positional frame (relative-to-primary),
        # sharing x limits forces misleading axis alignment across panels.
        share_x_axes = False
    share_y_axes = True
    if effective_axis_orientation == "position_y" and n_windows_effective > 1:
        # Same issue as above, but position is encoded on y.
        share_y_axes = False

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
                sharey=share_y_axes,
            )
            axes = np.atleast_1d(axes).reshape(1, n_windows_effective)
    else:
        if rotate:
            n_rows, n_cols = n_windows_effective, n_plot_motifs
            fig, axes = plt.subplots(
                n_rows,
                n_cols,
                figsize=(
                    max(9, 4.4 * n_cols),
                    max(
                        4.2,
                        2.7 * n_rows
                        + (0.35 if render_mode == "heatmap" else 0)
                        + (0.7 if threshold_map is not None else 0),
                    ),
                ),
                sharex=share_x_axes,
                sharey=share_y_axes,
            )
        else:
            n_rows, n_cols = n_plot_motifs, n_windows_effective
            fig, axes = plt.subplots(
                n_rows,
                n_cols,
                figsize=(
                    max(9, 4.2 * n_cols),
                    max(4.2, 2.5 * n_rows + (0.35 if render_mode == "heatmap" else 0)),
                ),
                sharey=share_y_axes,
            )
        axes = np.atleast_2d(axes)

    def _axis_for(motif_pos: int, window_pos: int):
        display_window_pos = window_plot_order.index(window_pos)
        if len(motif_ids) == 1:
            return axes[display_window_pos, 0] if rotate else axes[0, display_window_pos]
        if rotate:
            return axes[display_window_pos, motif_pos]
        return axes[motif_pos, display_window_pos]

    first_mappable = None
    first_mappable_by_motif: dict[int, Any] = {}
    threshold_legend_entries: dict[str, tuple[str, float]] = {}
    if threshold_map is not None:
        threshold_legend_entries = {
            str(motif_labels_resolved[motif_id]): (
                str(motif_color_map[motif_id]),
                float(threshold_map[motif_id]),
            )
            for motif_id in motif_ids
            if motif_id in threshold_map
        }
    for motif_pos, motif_id in enumerate(motif_ids):
        motif_panels = panels_sorted_by_motif[motif_pos]
        motif_raw_panels = panels_raw_sorted_by_motif[motif_pos]
        motif_smooth_panels = panels_smooth_sorted_by_motif[motif_pos]
        motif_ml_panels = panels_ml_sorted_by_motif[motif_pos]
        motif_label = str(motif_labels_resolved[motif_id])
        for window_pos in window_plot_order:
            display_window_pos = window_plot_order.index(window_pos)
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
                            first_mappable_by_motif.setdefault(motif_id, first_mappable)
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
                            first_mappable_by_motif.setdefault(motif_id, first_mappable)
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

            if window_pos == primary_window_index:
                for spine in ax.spines.values():
                    spine.set_linewidth(2.2)
                    spine.set_edgecolor("black")

            if render_mode == "heatmap":
                site_label = f"Site {window_pos + 1}"
                if window_pos == primary_window_index:
                    site_label += " (primary)"
                heatmap_title = (
                    f"{motif_label} | {site_label}"
                    if len(motif_ids) > 1
                    else site_label
                )
                ax.set_title(heatmap_title, fontsize=9)
            elif len(motif_ids) > 1:
                if rotate and display_window_pos == 0:
                    ax.set_title(motif_label)
                if not rotate and motif_pos == 0:
                    ax.set_title(f"Site {window_pos + 1}")

            if effective_axis_orientation == "position_x":
                if rotate and motif_pos == 0:
                    ax.set_ylabel("Reads (sorted)")
                if (not rotate) and len(motif_ids) > 1 and window_pos == 0:
                    ax.set_ylabel(f"{motif_label}\nReads")
            else:
                if rotate and motif_pos == 0:
                    ax.set_ylabel("Position (bp)")
                if (not rotate) and len(motif_ids) > 1 and window_pos == 0:
                    ax.set_ylabel(f"{motif_label}\nPosition")

    position_label = (
        "Position relative to primary center (bp)"
        if effective_coordinate_mode == "relative_to_primary"
        else "Distance from window center (bp)"
    )
    if rotate:
        bottom_window = window_plot_order[-1]
        bottom_axes = [
            _axis_for(motif_pos, bottom_window)
            for motif_pos in range(len(motif_ids))
        ]
        xlabel = (
            position_label
            if effective_axis_orientation == "position_x"
            else "Reads (sorted)"
        )
        for ax in bottom_axes:
            ax.set_xlabel(xlabel)
    else:
        for motif_pos in range(len(motif_ids)):
            left_window = window_plot_order[0]
            if effective_axis_orientation == "position_x":
                _axis_for(motif_pos, left_window).set_ylabel(
                    "Reads (sorted)"
                    if len(motif_ids) == 1
                    else _axis_for(motif_pos, left_window).get_ylabel()
                )
            else:
                _axis_for(motif_pos, left_window).set_ylabel(
                    "Position (bp)"
                    if len(motif_ids) == 1
                    else _axis_for(motif_pos, left_window).get_ylabel()
                )
            for window_pos in range(n_windows_effective):
                if effective_axis_orientation == "position_x":
                    _axis_for(motif_pos, window_pos).set_xlabel(
                        position_label
                    )
                else:
                    _axis_for(motif_pos, window_pos).set_xlabel("Reads (sorted)")

    sort_target = f"Site {sort_window_pos + 1}"
    if sort_window_pos == primary_window_index:
        sort_target += " (primary)"
    selection_label = selector.mode.replace("_", " ")
    title_parts = [
        f"{selection_label} raster",
        f"primary Site {primary_window_index + 1} shown with bold border",
        f"sorted by {sort_by} in {sort_target}",
        f"{effective_coordinate_mode}, {effective_axis_orientation}",
        f"rows {panels_sorted_by_motif[0][0].shape[0]} / {rows_before_downsample}",
    ]
    fig.suptitle(" | ".join(title_parts), fontsize=9, y=0.985)
    top_margin = 0.70 if threshold_map is not None else 0.88
    if render_mode == "heatmap":
        top_margin -= 0.03
    left_margin = 0.16 if rotate and n_windows_effective > 1 else 0.10
    fig.subplots_adjust(
        top=top_margin,
        right=0.88,
        left=left_margin,
        hspace=0.45,
        wspace=0.35,
    )

    if rotate and n_windows_effective > 1:
        for window_pos in window_plot_order:
            row_axis = _axis_for(0, window_pos)
            bbox = row_axis.get_position()
            site_label = f"Site {window_pos + 1}"
            if window_pos == primary_window_index:
                site_label += " (primary)"
            fig.text(
                max(0.01, bbox.x0 - 0.105),
                (bbox.y0 + bbox.y1) / 2,
                site_label,
                rotation=90,
                va="center",
                ha="center",
                fontsize=9,
                fontweight="bold" if window_pos == primary_window_index else "normal",
            )

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
        sigma_txt = f", σ={smooth_sigma_bp}" if smoothing == "gaussian" else ""
        if (
            render_mode == "scatter"
            and effective_scatter_color_values == "ml_score"
            and len(motif_ids) > 1
        ):
            for motif_id in motif_ids:
                motif_pos = motif_ids.index(motif_id)
                if motif_id not in first_mappable_by_motif:
                    first_mappable_by_motif[motif_id] = _axis_for(motif_pos, 0).scatter(
                        [],
                        [],
                        c=[],
                        cmap=ml_cmap_names[motif_id],
                        vmin=vmin,
                        vmax=vmax,
                    )
                cbar = fig.colorbar(
                    first_mappable_by_motif[motif_id],
                    ax=[
                        _axis_for(motif_pos, window_pos)
                        for window_pos in range(n_windows_effective)
                    ],
                    shrink=0.6,
                    pad=0.06,
                )
                cbar.set_label(
                    f"{motif_labels_resolved[motif_id]} normalized ML score (0-1)\n"
                    f"{'downsampled' if downsampled else 'full'}"
                )
        else:
            cbar = fig.colorbar(
                first_mappable, ax=np.ravel(axes).tolist(), shrink=0.6, pad=0.06
            )
            if (
                render_mode == "scatter"
                and effective_scatter_color_values == "ml_score"
            ):
                cbar.set_label(
                    f"normalized ML score (0-1)\n{'downsampled' if downsampled else 'full'}"
                )
            elif render_mode == "scatter":
                scale_label = "fraction modified signal"
                smoothing_label = (
                    f"{smoothing}, win={smooth_win}{sigma_txt}"
                    if render_value_mode == "smoothed"
                    else "none"
                )
                cbar.set_label(
                    f"{scale_label} (values={render_value_mode}, smoother={smoothing_label})\n"
                    f"{'downsampled' if downsampled else 'full'}"
                )
            else:
                scale_label = "Signal heatmap"
                smoothing_label = (
                    f"{smoothing} win={smooth_win}{sigma_txt}"
                    if render_value_mode == "smoothed"
                    else "none"
                )
                cbar.set_label(
                    f"{scale_label}\n"
                    f"{render_value_mode}; {smoothing_label}\n"
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
            loc="upper center",
            bbox_to_anchor=(0.5, 0.86),
            ncol=max(1, len(handles)),
            frameon=False,
        )

    return fig, {
        "pairs": panels_sorted_by_motif[0][0].shape[0],
        **selection_result.stats,
        "downsampled": downsampled,
        "downsample_method": (effective_downsample_method if downsampled else "none"),
        "downsample_factor": int(downsample_factor),
        "rows_before_downsample": rows_before_downsample,
        "rows_after_downsample": int(panels_sorted_by_motif[0][0].shape[0]),
        "vmin": vmin,
        "vmax": vmax,
        "render_mode": render_mode,
        "render_values": render_value_mode,
        "scatter_color_values": effective_scatter_color_values,
        "axis_orientation": effective_axis_orientation,
        "coordinate_mode": effective_coordinate_mode,
        "window_plot_order": [int(idx) for idx in window_plot_order],
        "sort_by": sort_by,
        "sort_window_index": sort_window_pos,
        "sort_descending": effective_sort_descending,
        "motifs_plotted": [motif_labels_resolved[m] for m in motif_ids],
        "n_windows": n_windows_effective,
        "selection_mode": selector.mode,
        "primary_window_index": int(primary_window_index),
        "window_offsets_bp": target_offsets.tolist() if target_offsets is not None else None,
        "window_widths_bp": [int(value) for value in panel_widths_bp],
        "required_read_length_bp": required_read_length_bp,
        "min_read_length_bp": effective_min_read_length,
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
    max_rows: int | None = 500,
    downsample_method: str = "auto",
    downsample_seed: int | None = None,
    **kwargs,
):
    """
    Convenience wrapper for two-site raster plots in primary-region coordinates.

    Args:
        read_windows: ReadWindowExtractionResult from extract_read_windows / build_multimotif_read_windows
        second_site_offset_bp: secondary site center (bp) relative to the primary site at 0
        window_width_bp: common width (bp) applied to both windows when provided
        max_rows: cap displayed reads; pass None to show all rows
        downsample_method: "auto" uses uniform read sampling for scatter and bin-mean for heatmap
        downsample_seed: seed for reproducible random downsampling
        **kwargs: forwarded to plot_multisite_read_raster
    """

    if (
        "window_offsets_bp" in kwargs
        or "selection_mode" in kwargs
        or "site_selection" in kwargs
    ):
        raise ValueError(
            "plot_two_site_read_raster manages fixed offsets internally; do not pass window_offsets_bp, selection_mode, or site_selection."
        )
    if "n_windows" in kwargs:
        raise ValueError("plot_two_site_read_raster always uses exactly two windows.")
    if "primary_window_index" in kwargs:
        raise ValueError("plot_two_site_read_raster always uses window 1 as primary.")

    forwarded: dict[str, Any] = dict(kwargs)
    forwarded["max_rows"] = max_rows
    forwarded["downsample_method"] = downsample_method
    forwarded["downsample_seed"] = downsample_seed
    forwarded["site_selection"] = {
        "mode": "fixed_offsets",
        "window_offsets_bp": [0, int(second_site_offset_bp)],
        "primary_window_index": 0,
    }
    if window_width_bp is not None:
        forwarded["window_widths_bp"] = [int(window_width_bp), int(window_width_bp)]
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
    "summarize_feature_matrix",
    "scale_feature_matrix",
    "rank_read_features_by_group_difference",
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
    "raster_stats_brief",
    "infer_region_source_labels",
]
