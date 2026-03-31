from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

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
    # Controls how windows are cut from per-read vectors
    window_size: int = 2000
    orientation_aware: bool = True
    filter_multi_region_reads: bool = False


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
    region_metadata = [
        (chromosome, start, end, strand)
        for chromosome, region_list in regions_dict.items()
        for start, end, strand in region_list
    ]
    if len(region_metadata) == 0:
        raise ValueError("No regions were found to build clustering features.")

    loader = partial(
        load_processed.pileup_vectors_from_bedmethyl,
        bedmethyl_file=bedmethyl_file,
        motif=motif,
        window_size=window_size,
        single_strand=single_strand,
        regions_5to3prime=regions_5to3prime,
        quiet=quiet,
        cores=cores,
        chunk_size=chunk_size,
    )
    per_region_vectors = load_processed.regions_to_list(
        function_handle=loader,
        regions=regions,
        window_size=window_size,
        quiet=quiet,
        cores=cores,
        split_large_regions=split_large_regions,
    )

    if len(per_region_vectors) != len(region_metadata):
        raise RuntimeError(
            "Region metadata and pileup vectors have mismatched lengths. "
            "Double-check the region inputs."
        )

    feature_vectors = []
    expected_length: int | None = None
    for modified_vec, valid_vec in per_region_vectors:
        modified = np.asarray(modified_vec, dtype=np.float64)
        valid = np.asarray(valid_vec, dtype=np.float64)
        if modified.shape != valid.shape:
            raise ValueError("Modified and valid vectors must have matching shapes.")
        if expected_length is None:
            expected_length = modified.shape[0]
        elif modified.shape[0] != expected_length:
            raise ValueError(
                "All regions must have the same length. "
                "Pass window_size to enforce length-matched regions."
            )

        denominator = valid + (2 * pseudo_count)
        with np.errstate(divide="ignore", invalid="ignore"):
            fraction = np.divide(
                modified + pseudo_count,
                denominator,
                out=np.zeros_like(modified),
                where=denominator > 0,
            )
        feature_vectors.append(fraction)

    feature_matrix = np.vstack(feature_vectors)
    return feature_matrix, region_metadata


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
    except ImportError as exc:  # pragma: no cover - exercised in environments w/o sklearn
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
            "Install xgboost to use classifier='xgboost'. "
            "Try `pip install xgboost`."
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


def _extract_window_from_record(
    record: tuple,
    idx: dict[str, int],
    window_size: int,
    orientation_aware: bool,
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
    val_vector = (
        np.asarray(record[idx["val_vector"]]) if "val_vector" in idx else None
    )

    if read_end <= read_start or region_end <= region_start:
        return None

    center = (region_start + region_end) // 2
    half = window_size // 2
    # Grab a symmetric window around the region center
    window_start = center - half
    window_end = window_start + window_size
    # Skip reads that do not fully span the requested window
    if window_start < read_start or window_end > read_end:
        return None

    slice_start = window_start - read_start
    slice_end = slice_start + window_size
    if slice_start < 0 or slice_end > len(mod_vector):
        return None

    mod_window = mod_vector[slice_start:slice_end]
    val_window = (
        val_vector[slice_start:slice_end] if val_vector is not None else None
    )

    if orientation_aware:
        region_strand = (
            _coerce_strand(record[idx["region_strand"]])
            if "region_strand" in idx
            else None
        )
        read_strand = (
            _coerce_strand(record[idx["strand"]]) if "strand" in idx else None
        )
        # If read/reference strands disagree, flip to align everything 5'->3'
        if _should_flip(region_strand, read_strand):
            mod_window = np.flip(mod_window)
            if val_window is not None:
                val_window = np.flip(val_window)

    return mod_window, val_window


def build_multimotif_read_windows(
    hdf5_file: str | Path,
    motifs: Sequence[str],
    regions: str | Path | list[str | Path] | None = None,
    *,
    window_size: int,
    orientation_aware: bool = True,
    single_strand: bool = False,
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
    require_all_motifs: bool = True,
) -> ReadWindowExtractionResult:
    """
    Group per-motif rows by read and return a combined window per read, concatenating motifs.

    Each requested motif contributes a window of length `window_size`; missing motifs are filled with zeros.
    Windows are concatenated in the order provided by `motifs`.

    Args:
        hdf5_file: path to extract .h5 file
        motifs: motifs to include (order matters for concatenation)
        regions: optional region filter
        window_size: fixed window length around region center to extract per motif
        orientation_aware: flip windows if read strand != region strand
        single_strand: passed to loader
        subset_parameters: passed to loader for subsetting
        span_full_window: passed to loader (if True, only reads spanning the region are loaded)
        require_all_motifs: if True, drop reads that are missing any requested motif

    Returns:
        ReadWindowExtractionResult with data_matrix of shape (n_reads, len(motifs) * window_size)
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
    groups: defaultdict[tuple[Any, Any, Any, Any, Any], dict[str, tuple]] = defaultdict(dict)
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

    for key, motif_map in groups.items():
        motif_windows = []
        motif_val_windows = []
        motifs_present = []
        for motif in motifs:
            if motif in motif_map:
                extracted = _extract_window_from_record(
                    motif_map[motif], idx, window_size, orientation_aware
                )
                if extracted is None:
                    # If window cannot be extracted, treat as missing
                    motif_windows.append(np.zeros(window_size, dtype=float))
                    motif_val_windows.append(
                        np.zeros(window_size, dtype=float) if has_val else None
                    )
                else:
                    mw, vw = extracted
                    motif_windows.append(np.asarray(mw, dtype=float))
                    motif_val_windows.append(
                        np.asarray(vw, dtype=float) if vw is not None else None
                    )
                    motifs_present.append(motif)
            else:
                motif_windows.append(np.zeros(window_size, dtype=float))
                motif_val_windows.append(np.zeros(window_size, dtype=float) if has_val else None)

        if require_all_motifs and len(motifs_present) < len(motifs):
            continue

        combined_windows.append(np.concatenate(motif_windows, axis=0))
        if has_val:
            combined_vals.append(
                np.concatenate(
                    [vw if vw is not None else np.zeros(window_size, dtype=float) for vw in motif_val_windows],
                    axis=0,
                )
            )
        metadata.append(
            {
                "read_name": key[0],
                "chromosome": key[1],
                "region_start": key[2],
                "region_end": key[3],
                "region_strand": key[4],
                "motifs_present": motifs_present,
            }
        )

    if not combined_windows:
        raise ValueError("No reads produced combined motif windows; check inputs and window_size.")

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
        config: controls window size and orientation handling
        window_size: overrides config.window_size when provided

    Returns:
        ReadWindowExtractionResult containing raw mod/val matrices and metadata
    """

    cfg = config or ReadWindowExtractionConfig()
    effective_window = window_size or cfg.window_size

    # Load all requested reads/vectors from extract output
    read_tuples, dataset_names, regions_dict = load_processed.read_vectors_from_hdf5(
        file=hdf5_file,
        motifs=list(motifs),
        regions=regions,
        window_size=window_size,
        single_strand=single_strand,
        subset_parameters=subset_parameters,
        span_full_window=span_full_window,
    )
    idx = _build_dataset_index(dataset_names)

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
            rec, idx, effective_window, cfg.orientation_aware
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
            "region_start",
            "region_end",
            "region_strand",
            "strand",
            "motif",
        ]
        metadata.append(
            {
                field: rec[idx[field]] if field in idx else None
                for field in meta_fields
            }
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
        values = np.array([compute_autocorrelation_feature(row, lag) for row in data_matrix])
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
    try:
        scores["silhouette"] = float(silhouette_score(X, labels))
    except Exception:  # pragma: no cover
        pass
    try:
        scores["calinski_harabasz"] = float(calinski_harabasz_score(X, labels))
    except Exception:  # pragma: no cover
        pass
    try:
        scores["davies_bouldin"] = float(davies_bouldin_score(X, labels))
    except Exception:  # pragma: no cover
        pass
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
    for raw_label, ordered_label in zip(labels_raw, labels_size_ordered):
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
        for lbl, p in zip(unique, probs):
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
        AgglomerativeClustering,
        Birch,
        DBSCAN,
        KMeans,
        MiniBatchKMeans,
        OPTICS,
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
                min_samples=kwargs.get("min_samples", None),
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
    view_bp: int | None = None,
    window_bp: int | None = None,
    motif_index: int = 0,
    motif_count: int | None = None,
    cmap_name: str = "viridis",
    invert_y: bool = True,
    point_size: float = 1.0,
    point_alpha: float = 0.01,
):
    """
    Scatter plot of per-read modification calls with per-cluster average profiles.
    Useful for quick QC after clustering.
    Set motif_count when data_matrix contains concatenated motifs to select motif_index safely.
    """

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.colors import Normalize

    X_full = np.asarray(data_matrix)
    labs = np.asarray(labels)
    if include_clusters is not None:
        mask = np.isin(labs, list(include_clusters))
        if not np.any(mask):
            raise ValueError("No reads remain after filtering include_clusters.")
        X_full = X_full[mask]
        labs = labs[mask]
        if val_matrix is not None:
            val_matrix = np.asarray(val_matrix)[mask]

    # If concatenated motifs are present, select the requested motif slice
    if motif_index is not None and X_full.ndim == 2:
        window_len = X_full.shape[1]
        # Prefer explicit motif_count; otherwise infer from window_bp/view_bp if they evenly tile the window
        inferred_count = None
        if motif_count is not None:
            inferred_count = int(max(1, motif_count))
        elif window_bp and window_bp > 0 and window_len % window_bp == 0:
            inferred_count = window_len // window_bp
        elif view_bp and view_bp > 0 and window_len % view_bp == 0:
            inferred_count = window_len // view_bp
        n_motifs = inferred_count or 1
        slice_width = window_len // n_motifs if n_motifs > 0 else window_len
        if slice_width == 0:
            slice_width = window_len
            n_motifs = 1
        start = motif_index * slice_width
        end = min(window_len, start + slice_width)
        if start >= window_len:
            raise ValueError(
                f"motif_index {motif_index} exceeds concatenated width {window_len}; "
                "pass motif_count/window_bp to disambiguate."
            )
        X = X_full[:, start:end]
        V = val_matrix[:, start:end] if val_matrix is not None else None
    else:
        X = X_full
        V = val_matrix

    n_reads, window = X.shape
    window_bp = window_bp or window
    half = window_bp // 2
    if view_bp is not None:
        half = view_bp // 2
    x_axis = np.arange(-half, -half + X.shape[1])

    mod_fraction = X.mean(axis=1)
    order = np.lexsort((-mod_fraction, labs))
    X = X[order]
    labs = labs[order]
    V = V[order] if V is not None else None

    unique_clusters = np.unique(labs)
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=int(unique_clusters.min()), vmax=int(unique_clusters.max()))

    rows, cols = np.nonzero(X)
    colors = cmap(norm(labs[rows]))

    fig = plt.figure(figsize=(12, max(4, len(unique_clusters))))
    gs = gridspec.GridSpec(nrows=len(unique_clusters), ncols=2, width_ratios=[3, 1])
    ax_left = fig.add_subplot(gs[:, 0])
    ax_left.scatter(x_axis[cols], rows, s=point_size, alpha=point_alpha, c=colors)
    ax_left.set_xlabel("Distance from region center (bp)")
    ax_left.set_ylabel("Sorted read index")
    if invert_y:
        ax_left.invert_yaxis()

    change_points = np.flatnonzero(np.diff(labs)) + 1
    for cp in change_points:
        ax_left.axhline(cp, color="0.2", linestyle="--", linewidth=0.3)

    for i, cluster_id in enumerate(unique_clusters):
        ax = fig.add_subplot(gs[i, 1])
        mask = labs == cluster_id
        mean_profile = X[mask].mean(axis=0)
        ax.plot(x_axis, mean_profile, color=cmap(norm(cluster_id)))
        if V is not None:
            mean_val = V[mask].mean(axis=0)
            ax.plot(x_axis, mean_val, color="0.35", linestyle="--", linewidth=1.0)
        ax.set_title(f"Cluster {cluster_id} (n={mask.sum()})")
        ax.set_xlim(x_axis[0], x_axis[-1])
        ax.set_ylim(0, min(1.0, max(0.05, mean_profile.max() + 0.05)))

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

    le = LabelEncoder()
    y = le.fit_transform(y_raw)
    classes = le.classes_
    if len(classes) != 2:
        raise ValueError("Binary classification requires exactly two sample labels.")

    X_train, X_test, y_train, y_test, lbl_train, lbl_test = train_test_split(
        X, y, y_raw, test_size=test_size, random_state=random_state, stratify=y
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
            max_depth=kwargs.get("max_depth", None),
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
            n_jobs=kwargs.get("n_jobs", None),
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
        }
    )

    all_preds = pd.concat([train_df, preds_df], ignore_index=True)

    return {
        "model": model,
        "metrics": metrics,
        "predictions": all_preds,
        "label_mapping": {cls: int(code) for cls, code in zip(classes, le.transform(classes))},
    }


def plot_confusion_matrices(predictions: pd.DataFrame, *, cmap: str = "Blues"):
    """
    Plot train/test confusion matrices from the predictions DataFrame returned by classify_read_features_binary.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.metrics import confusion_matrix

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    for ax, split in zip(axes, ["train", "test"]):
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


def plot_classification_profiles(
    data_matrix: np.ndarray,
    predictions: pd.DataFrame,
    *,
    val_matrix: np.ndarray | None = None,
    split: str = "test",
    motif_index: int = 0,
    view_bp: int | None = None,
    cmap_name: str = "viridis",
    invert_y: bool = True,
    point_size: float = 1.0,
    point_alpha: float = 0.01,
):
    """
    Visualize per-read modification calls colored by predicted label, with mean profiles per predicted class.

    Args:
        data_matrix: per-read windows (n_reads, window_len)
        predictions: DataFrame from classify_read_features_binary (must align row-wise with data_matrix)
        val_matrix: optional valid-site matrix to overlay
        split: which split to visualize ('train' or 'test')
    """

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.colors import Normalize

    preds = predictions[predictions["split"] == split].reset_index(drop=True)
    if len(preds) != data_matrix.shape[0]:
        raise ValueError("predictions and data_matrix must have matching rows for the chosen split.")

    labels = preds["pred_label"].to_numpy()
    X_full = np.asarray(data_matrix)
    V = np.asarray(val_matrix) if val_matrix is not None else None

    # Slice motif window if concatenated
    window_len = X_full.shape[1]
    motif_width = window_len // max(1, (window_len // (view_bp or window_len))) if window_len else window_len
    if motif_index is not None and X_full.ndim == 2:
        n_motifs = max(1, window_len // (view_bp or window_len))
        slice_width = window_len // n_motifs
        start = motif_index * slice_width
        end = start + slice_width
        X = X_full[:, start:end]
        V = V[:, start:end] if V is not None else None
    else:
        X = X_full

    n_reads, window = X.shape
    half = window // 2
    if view_bp is not None:
        half = view_bp // 2
    x_axis = np.arange(-half, -half + X.shape[1])

    mod_fraction = X.mean(axis=1)
    order = np.lexsort((-mod_fraction, labels))
    X = X[order]
    labels = labels[order]
    V = V[order] if V is not None else None

    unique_labels = np.unique(labels)
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=0, vmax=max(len(unique_labels) - 1, 1))

    rows, cols = np.nonzero(X)
    colors = cmap(norm(pd.factorize(labels)[0][rows]))

    fig = plt.figure(figsize=(12, max(4, len(unique_labels))))
    gs = gridspec.GridSpec(nrows=len(unique_labels), ncols=2, width_ratios=[3, 1])
    ax_left = fig.add_subplot(gs[:, 0])
    ax_left.scatter(x_axis[cols], rows, s=point_size, alpha=point_alpha, c=colors)
    ax_left.set_xlabel("Distance from region center (bp)")
    ax_left.set_ylabel("Sorted read index")
    if invert_y:
        ax_left.invert_yaxis()

    change_points = np.flatnonzero(np.diff(pd.factorize(labels)[0])) + 1
    for cp in change_points:
        ax_left.axhline(cp, color="0.2", linestyle="--", linewidth=0.3)

    for i, lbl in enumerate(unique_labels):
        ax = fig.add_subplot(gs[i, 1])
        mask = labels == lbl
        mean_profile = X[mask].mean(axis=0)
        ax.plot(x_axis, mean_profile, color=cmap(norm(i)), label=str(lbl))
        if V is not None:
            mean_val = V[mask].mean(axis=0)
            ax.plot(x_axis, mean_val, color="0.35", linestyle="--", linewidth=1.0)
        ax.set_title(f"{lbl} (n={mask.sum()})")
        ax.set_xlim(x_axis[0], x_axis[-1])
        ax.set_ylim(0, min(1.0, max(0.05, mean_profile.max() + 0.05)))
    fig.tight_layout()
    return fig


def plot_multisite_read_raster(
    read_windows: ReadWindowExtractionResult,
    *,
    n_windows: int = 2,
    min_separation_bp: int = 5000,
    distance_mode: str = "center",  # or "bounds"
    window_bp: int | None = None,
    motif_index: int = 0,
    smoothing: str | None = "gaussian",  # None, "boxcar", "gaussian"
    smooth_win: int = 21,
    smooth_sigma_bp: float = 6.0,
    max_rows: int = 500,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    sort_by: str = "mod_fraction",  # or "cluster"
    sort_labels: Sequence[int] | None = None,
    beds: Sequence[str | Path] | None = None,
    rotate: bool = True,
):
    """
    Plot n windows per read (e.g., paired sites) as stacked rasters, sorted by cluster or mod fraction.

    Args:
        read_windows: ReadWindowExtractionResult from extract_read_windows / build_multimotif_read_windows
        n_windows: number of windows per read to plot (e.g., 2 for pairs)
        min_separation_bp: minimum separation between site centers (distance_mode="center") or read starts (bounds)
        distance_mode: "center" uses region_start/end center; "bounds" uses read_start/read_end
        window_bp: overrides default window length inferred from data_matrix slice
        motif_index: which motif slice to use if data_matrix concatenates motifs
        smoothing: None, "boxcar", or "gaussian"
        smooth_win: smoothing window length
        smooth_sigma_bp: sigma for gaussian smoothing
        max_rows: cap rows; downsample by averaging if exceeded
        cmap: matplotlib colormap
        vmin/vmax: color scale limits; None auto-scales
        sort_by: "mod_fraction" or "cluster"
        sort_labels: cluster labels to sort by when sort_by="cluster"
        beds: optional list of BED/GTF/GFF paths; only reads intersecting all beds are retained
        rotate: if True, rotate rasters for vertical stacking
    """

    import matplotlib.pyplot as plt
    import math

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
        k = _kernel("gaussian" if smoothing == "gaussian" else "boxcar", smooth_win, smooth_sigma_bp)
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

    meta = read_windows.metadata
    X_full = read_windows.data_matrix
    # motif slicing
    window_len = X_full.shape[1]
    slice_width = window_bp or (window_len // max(1, window_len // (window_bp or window_len)))
    n_motifs = max(1, window_len // slice_width)
    start = motif_index * slice_width
    end = start + slice_width
    X = X_full[:, start:end]

    # Filter by beds if provided
    keep_idx = np.arange(len(meta))
    if bed_filters:
        keep_mask = []
        for i, m in enumerate(meta):
            try:
                chrom = m.get("chromosome")
                rs = int(m.get("region_start", -1))
                re = int(m.get("region_end", -1))
            except Exception:
                keep_mask.append(False)
                continue
            pr_row = pr.PyRanges(pd.DataFrame({"Chromosome": [chrom], "Start": [rs], "End": [re]}))
            ok = all(len(pr_row.join(bf)) > 0 for bf in bed_filters)
            keep_mask.append(ok)
        keep_idx = np.flatnonzero(keep_mask)
        X = X[keep_idx]
        meta = [meta[i] for i in keep_idx]

    # Group by read key
    groups = defaultdict(list)
    for i, m in zip(range(len(meta)), meta):
        key = (
            m.get("read_name"),
            m.get("chromosome"),
            m.get("region_strand"),
        )
        if distance_mode == "center":
            center = (int(m.get("region_start", 0)) + int(m.get("region_end", 0))) // 2
        else:
            center = int(m.get("read_start", 0))
        groups[key].append((i, m.get("chromosome"), center))

    panels = [[] for _ in range(n_windows)]
    for _, items in groups.items():
        if len(items) < n_windows:
            continue
        items = sorted(items, key=lambda x: (x[1], x[2]))
        # sliding window of n_windows with separation constraint
        for j in range(len(items) - n_windows + 1):
            centers = [items[j + k][2] for k in range(n_windows)]
            if any((centers[k + 1] - centers[k]) < min_separation_bp for k in range(n_windows - 1)):
                continue
            idxs = [items[j + k][0] for k in range(n_windows)]
            for k, idx in enumerate(idxs):
                panels[k].append(X[idx])

    if not all(len(p) > 0 for p in panels):
        raise ValueError("No read sets found meeting separation criteria.")

    panels = [np.vstack(p) for p in panels]
    # Smooth
    panels_smooth = [_smooth(p) for p in panels]

    # Sorting
    if sort_by == "mod_fraction":
        means = [p.mean(axis=1) for p in panels_smooth]
        order = np.lexsort(tuple(-m for m in reversed(means)))
    elif sort_by == "cluster" and sort_labels is not None:
        sort_labels = np.asarray(sort_labels)[keep_idx]
        # repeat labels to match panels length
        order = np.argsort(sort_labels)
    else:
        order = np.arange(panels_smooth[0].shape[0])

    panels_sorted = [p[order] for p in panels_smooth]

    # Downsample rows if too tall
    def bin_rows(M, k):
        n_bins = math.ceil(M.shape[0] / k)
        out = np.zeros((n_bins, M.shape[1]))
        for i in range(n_bins):
            s, e = i * k, min((i + 1) * k, M.shape[0])
            out[i] = M[s:e].mean(axis=0)
        return out

    P = panels_sorted[0].shape[0]
    downsampled = False
    if P > max_rows:
        step = math.ceil(P / max_rows)
        panels_sorted = [bin_rows(p, step) for p in panels_sorted]
        P = panels_sorted[0].shape[0]
        downsampled = True

    if vmin is None or vmax is None:
        vmax_auto = max(p.max() for p in panels_sorted)
        vmin = 0.0 if vmin is None else vmin
        vmax = vmax_auto if vmax is None else vmax

    # Plot
    import matplotlib.pyplot as plt

    if rotate:
        fig, axes = plt.subplots(
            n_windows, 1, figsize=(max(6, P * 0.06), 3 + n_windows * 1.5), sharex=True
        )
        axes = np.atleast_1d(axes)
        for i, ax in enumerate(axes):
            im = ax.imshow(
                panels_sorted[i].T,
                aspect="auto",
                origin="lower",
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
            )
            ax.axhline(slice_width // 2, linestyle="--", linewidth=1, color="0.7")
            ax.set_ylabel(f"Window {i+1} (bp)")
        axes[-1].set_xlabel("Reads (sorted)")
    else:
        fig, axes = plt.subplots(
            1, n_windows, figsize=(10, max(3, P * 0.06)), sharey=True
        )
        axes = np.atleast_1d(axes)
        for i, ax in enumerate(axes):
            im = ax.imshow(
                panels_sorted[i],
                aspect="auto",
                origin="upper",
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
            )
            ax.axvline(slice_width // 2 - 0.5, linestyle="--", linewidth=1, color="0.7")
            ax.set_xlabel("Position (bp)")
        axes[0].set_ylabel("Reads (sorted)")

    cbar = fig.colorbar(im, ax=np.ravel(axes).tolist(), shrink=0.6, pad=0.02)
    sigma_txt = f", σ={smooth_sigma_bp}" if smoothing == "gaussian" else ""
    cbar.set_label(
        f"Modification density (smoothed={smoothing}, win={smooth_win}{sigma_txt})\\n"
        f"{'downsampled' if downsampled else 'full'}"
    )

    return fig, {"pairs": panels_sorted[0].shape[0], "downsampled": downsampled, "vmin": vmin, "vmax": vmax}


def plot_region_cluster_profiles(
    pileup_matrix: np.ndarray,
    labels: np.ndarray,
    *,
    window_bp: int | None = None,
    motif_index: int | None = None,
    motif_count: int | None = None,
    show_all_motifs: bool = False,
    cmap_name: str = "viridis",
):
    """
    Visualize region-level clustering by overlaying per-cluster mean profiles and a heatmap.

    Args:
        pileup_matrix: 2D array (n_regions, region_length) of modification fractions
        labels: cluster assignments per region
        window_bp: optional x-axis width in bp (assumes symmetric window); defaults to range(len(profile))
        motif_index: if pileup_matrix is a concatenation of multiple motifs, select which motif slice to plot
        motif_count: total number of motifs concatenated; if None, inferred when possible
        show_all_motifs: if True and motifs are concatenated, plot mean profiles for each motif slice per cluster
        cmap_name: matplotlib colormap name for clusters
    """

    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    from matplotlib.colors import Normalize

    X_full = np.asarray(pileup_matrix)
    labs = np.asarray(labels)
    if X_full.shape[0] != labs.shape[0]:
        raise ValueError("pileup_matrix rows and labels must have the same length.")

    # Motif slicing
    region_len = X_full.shape[1]
    motif_idx = motif_index if motif_index is not None else 0
    if motif_count is None and window_bp and window_bp > 0 and region_len % window_bp == 0:
        motif_count = region_len // window_bp
    motif_count = max(1, int(motif_count or 1))
    slice_width = region_len // motif_count if motif_count > 0 else region_len
    if slice_width == 0:
        slice_width = region_len
        motif_count = 1

    start = motif_idx * slice_width
    end = min(region_len, start + slice_width)
    if start >= region_len:
        raise ValueError(
            f"motif_index {motif_idx} is out of range for concatenated length {region_len}; "
            "pass motif_count to disambiguate."
        )
    X = X_full[:, start:end]

    n_regions, region_len = X.shape
    if window_bp is None:
        x_axis = np.arange(region_len)
    else:
        half = window_bp // 2
        x_axis = np.linspace(-half, half, region_len)

    # Sort regions by cluster for heatmap readability
    order = np.lexsort((np.arange(n_regions), labs))
    X_sorted = X[order]
    labs_sorted = labs[order]

    unique_clusters = np.unique(labs_sorted)
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=int(unique_clusters.min()), vmax=int(unique_clusters.max()))

    fig = plt.figure(figsize=(10, 6 if not show_all_motifs else 8))
    rows = 2 if not show_all_motifs else (2 + max(0, motif_count - 1))
    gs = gridspec.GridSpec(rows, 1, height_ratios=[1] * (rows - 1) + [2], hspace=0.3)

    # Mean profiles (optionally per motif)
    ax0 = fig.add_subplot(gs[0])
    for cl in unique_clusters:
        mask = labs == cl
        mean_profile = X[mask].mean(axis=0)
        ax0.plot(x_axis, mean_profile, color=cmap(norm(int(cl))), label=f"C{cl} (n={mask.sum()})")
    ax0.set_ylabel("Modified fraction")
    ax0.legend(loc="upper right", ncol=min(len(unique_clusters), 4), fontsize="small")

    if show_all_motifs and motif_count and motif_count > 1:
        for mi in range(1, motif_count):
            ax_m = fig.add_subplot(gs[mi])
            start_m = mi * slice_width
            end_m = start_m + slice_width
            Xm = X_full[:, start_m:end_m]
            for cl in unique_clusters:
                mask = labs == cl
                ax_m.plot(
                    x_axis,
                    Xm[mask].mean(axis=0),
                    color=cmap(norm(int(cl))),
                    label=f"C{cl}",
                )
            ax_m.set_ylabel(f"Motif {mi}")

    # Heatmap sorted by cluster (using selected motif slice)
    ax1 = fig.add_subplot(gs[-1])
    im = ax1.imshow(
        X_sorted,
        aspect="auto",
        interpolation="none",
        extent=[x_axis[0], x_axis[-1], n_regions, 0],
        cmap="viridis",
    )
    ax1.set_xlabel("Position (bp)" if window_bp is not None else "Index")
    ax1.set_ylabel("Region (sorted by cluster)")
    fig.colorbar(im, ax=ax1, label="Modified fraction")

    # Cluster boundaries
    change_points = np.flatnonzero(np.diff(labs_sorted)) + 1
    for cp in change_points:
        ax1.axhline(cp, color="0.8", linestyle="--", linewidth=0.5)

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
    for (chrom, start, end, strand), label in zip(region_metadata, labels):
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
    cmap_name: str = "viridis",
    linewidth: float = 6.0,
    figsize_per_chrom: float = 0.4,
):
    """
    Plot cluster-labeled regions along chromosomes (ideogram-style), colored by cluster.

    Args:
        region_bed: BED file, typically from export_region_clusters_to_bed, with Name carrying cluster id
        chrom_sizes: 2-column TSV (chrom, length) e.g. ref.fasta.fai
        cmap_name: matplotlib colormap to map cluster ids to colors
        linewidth: thickness of region segments
        figsize_per_chrom: vertical size per chromosome for sizing the figure
    """
    import matplotlib.pyplot as plt

    clusters_df = pd.read_csv(
        region_bed,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["Chromosome", "Start", "End", "Name"],
    )
    chrom_df = pd.read_csv(
        chrom_sizes,
        sep="\t",
        header=None,
        usecols=[0, 1],
        names=["Chromosome", "Length"],
    )
    df = clusters_df.merge(chrom_df, on="Chromosome", how="inner")
    if df.empty:
        raise ValueError("No overlapping chromosomes between BED and chrom_sizes.")

    df["start_frac"] = df["Start"] / df["Length"]
    df["end_frac"] = df["End"] / df["Length"]
    df["cluster_id"] = (
        df["Name"].astype(str).str.extract(r"(\d+)$", expand=False).astype(int)
    )

    chrom_order = (
        df.groupby("Chromosome")["Length"].max().sort_values(ascending=False).index
    )
    cmap = plt.get_cmap(cmap_name)
    fig_height = max(3.0, figsize_per_chrom * len(chrom_order))
    fig, ax = plt.subplots(figsize=(10, fig_height))

    for i, chrom in enumerate(chrom_order):
        ax.plot([0, 1], [i, i], color="0.9", lw=linewidth)  # chromosome backbone
        sub = df[df["Chromosome"] == chrom]
        for _, row in sub.iterrows():
            ax.plot(
                [row.start_frac, row.end_frac],
                [i, i],
                color=cmap(row.cluster_id % cmap.N),
                lw=linewidth,
                solid_capstyle="butt",
            )
        ax.text(-0.02, i, chrom, ha="right", va="center", fontsize=8)

    ax.set_xlim(0, 1)
    ax.set_ylim(-1, len(chrom_order))
    ax.set_yticks([])
    ax.set_xlabel("Relative position along chromosome")
    ax.set_title("Clustered regions by chromosome")
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
    df["cluster"] = df["Name"].astype(str).str.replace(f"{name_prefix}_", "", regex=False).astype(int)
    df["feature"] = df[feature_col].astype(str)

    counts = df.groupby(["cluster", "feature"]).size().reset_index(name="count")
    totals = counts.groupby("cluster")["count"].sum().reset_index(name="total_cluster")
    merged = counts.merge(totals, on="cluster", how="left")
    merged["fraction"] = merged["count"] / merged["total_cluster"]
    return merged


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

    def region_key(meta: dict[str, Any]) -> tuple[Any, Any, Any, Any]:
        strand = meta.get("region_strand") if include_strand else None
        return (
            meta.get("chromosome"),
            meta.get("region_start"),
            meta.get("region_end"),
            strand if strand in {"+", "-"} else ".",
        )

    rows = []
    for meta, label in zip(metadata, labels):
        key = region_key(meta)
        rows.append((*key, label))

    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "cluster"])
    if df.empty:
        return df

    counts = df.groupby(["chrom", "start", "end", "strand", "cluster"]).size()
    counts = counts.reset_index(name="count")
    pivot = counts.pivot_table(
        index=["chrom", "start", "end", "strand"],
        columns="cluster",
        values="count",
        fill_value=0,
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

    dom, dom_frac, ent = zip(*pivot.apply(_dominance_and_entropy, axis=1))
    pivot["dominant_cluster"] = dom
    pivot["dominant_fraction"] = dom_frac
    pivot["entropy"] = ent

    pivot = pivot.reset_index()
    return pivot


__all__ = [
    "cluster_features",
    "read_mod_fraction_table",
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
    "summarize_read_clusters_by_region",
    "plot_cluster_karyotype",
    "plot_classification_profiles",
    "plot_confusion_matrices",
    "sample_rows",
    "cluster_label_mapping",
    "apply_cluster_label_mapping",
    "plot_multisite_read_raster",
]
