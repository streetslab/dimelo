import numpy as np
import pytest

from dimelo import cluster


def test_region_feature_matrix_from_pileup(monkeypatch):
    # Mock the region parsing and loading helpers so we can focus on the feature logic
    fake_regions_dict = {
        "chr1": [(0, 3, "+")],
        "chr2": [(5, 8, "-")],
    }
    monkeypatch.setattr(
        cluster.utils, "regions_dict_from_input", lambda *args, **kwargs: fake_regions_dict
    )
    monkeypatch.setattr(
        cluster.load_processed,
        "regions_to_list",
        lambda **kwargs: [
            (np.array([1, 1, 0], dtype=float), np.array([2, 2, 1], dtype=float)),
            (np.array([0, 2, 2], dtype=float), np.array([1, 4, 4], dtype=float)),
        ],
    )

    features, metadata = cluster.region_feature_matrix_from_pileup(
        bedmethyl_file="pileup.bed.gz",
        motif="A,0",
        regions="regions.bed",
        pseudo_count=0.0,
    )

    np.testing.assert_allclose(
        features,
        np.array(
            [
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
            ]
        ),
    )
    assert metadata == [("chr1", 0, 3, "+"), ("chr2", 5, 8, "-")]


def test_read_mod_fraction_table(monkeypatch):
    datasets = [
        "read_name",
        "chromosome",
        "motif",
        "mod_vector",
        "val_vector",
        "region_start",
        "region_end",
        "region_strand",
        "read_length",
        "A,0_mod_fraction",
        "CG,0_mod_fraction",
    ]
    read_records = [
        ("read1", "chr1", "A,0", None, None, 10, 110, "+", 100, 0.25, 0.75),
        ("read2", "chr1", "CG,0", None, None, 20, 140, "-", 120, 0.10, 0.60),
    ]
    monkeypatch.setattr(
        cluster.load_processed,
        "read_vectors_from_hdf5",
        lambda **kwargs: (read_records, datasets, {"chr1": [(0, 100, "+")]}),
    )

    features, feature_names, metadata, regions_dict = cluster.read_mod_fraction_table(
        hdf5_file="reads.h5",
        motifs=["A,0", "CG,0"],
        metadata_fields=("read_name", "chromosome"),
    )

    np.testing.assert_allclose(
        features,
        np.array(
            [
                [0.25, 0.75],
                [0.10, 0.60],
            ]
        ),
    )
    assert feature_names == ["A,0_mod_fraction", "CG,0_mod_fraction"]
    assert metadata == [
        {"read_name": "read1", "chromosome": "chr1"},
        {"read_name": "read2", "chromosome": "chr1"},
    ]
    assert regions_dict == {"chr1": [(0, 100, "+")]}


def test_cluster_features_invokes_kmeans(monkeypatch):
    class DummyKMeans:
        def __init__(self, n_clusters, random_state, **kwargs):
            self.n_clusters = n_clusters
            self.random_state = random_state
            self.kwargs = kwargs
            self.observed = None

        def fit_predict(self, matrix):
            self.observed = matrix
            return np.arange(matrix.shape[0])

    monkeypatch.setattr(cluster, "_get_kmeans", lambda: DummyKMeans)

    feature_matrix = np.array([[0.1, 0.2], [0.4, 0.6], [0.9, 0.95]])
    labels, estimator = cluster.cluster_features(
        feature_matrix, n_clusters=3, random_state=7
    )

    np.testing.assert_array_equal(labels, np.array([0, 1, 2]))
    assert isinstance(estimator, DummyKMeans)
    assert estimator.n_clusters == 3
    np.testing.assert_allclose(estimator.observed, feature_matrix)


def test_cluster_features_rejects_unknown_method():
    with pytest.raises(ValueError):
        cluster.cluster_features(np.ones((2, 2)), method="dbscan")


def test_extract_read_windows_orientation(monkeypatch):
    dataset_names = [
        "chromosome",
        "mod_vector",
        "motif",
        "read_end",
        "read_name",
        "read_start",
        "strand",
        "val_vector",
        "region_start",
        "region_end",
        "region_strand",
        "read_length",
    ]
    record = (
        "chr1",
        np.arange(10, dtype=float),
        "A,0",
        110,
        "read1",
        100,
        "-",
        np.ones(10, dtype=float),
        104,
        108,
        "+",
        10,
    )

    def fake_loader(**kwargs):
        return ([record], dataset_names, None)

    monkeypatch.setattr(cluster.load_processed, "read_vectors_from_hdf5", fake_loader)
    result = cluster.extract_read_windows(
        hdf5_file="reads.h5",
        motifs=["A,0"],
        config=cluster.ReadWindowExtractionConfig(window_size=4, orientation_aware=True),
    )

    np.testing.assert_allclose(result.data_matrix, np.array([[7.0, 6.0, 5.0, 4.0]]))
    assert result.val_matrix is not None
    assert len(result.metadata) == 1


def test_extract_read_windows_filter_multi_region(monkeypatch):
    dataset_names = [
        "chromosome",
        "mod_vector",
        "motif",
        "read_end",
        "read_name",
        "read_start",
        "strand",
        "val_vector",
        "region_start",
        "region_end",
        "region_strand",
        "read_length",
    ]
    rec1 = (
        "chr1",
        np.zeros(10, dtype=float),
        "A,0",
        110,
        "dup",
        100,
        "+",
        np.ones(10, dtype=float),
        104,
        108,
        "+",
        10,
    )
    rec2 = (
        "chr1",
        np.ones(10, dtype=float),
        "A,0",
        115,
        "dup",
        100,
        "+",
        np.ones(10, dtype=float),
        106,
        110,
        "+",
        10,
    )
    rec3 = (
        "chr1",
        np.arange(10, dtype=float),
        "A,0",
        112,
        "keep",
        100,
        "+",
        np.ones(10, dtype=float),
        104,
        108,
        "+",
        10,
    )

    def fake_loader(**kwargs):
        return ([rec1, rec2, rec3], dataset_names, None)

    monkeypatch.setattr(cluster.load_processed, "read_vectors_from_hdf5", fake_loader)
    result = cluster.extract_read_windows(
        hdf5_file="reads.h5",
        motifs=["A,0"],
        config=cluster.ReadWindowExtractionConfig(
            window_size=4, orientation_aware=False, filter_multi_region_reads=True
        ),
    )
    assert result.data_matrix.shape[0] == 1
    assert result.metadata[0]["read_name"] == "keep"


def test_read_window_feature_matrix():
    data = np.array(
        [
            [0, 0, 1, 1],
            [1, 1, 0, 0],
            [0, 1, 0, 1],
        ],
        dtype=float,
    )
    result = cluster.ReadWindowExtractionResult(
        data_matrix=data,
        val_matrix=None,
        metadata=[],
        datasets=[],
        regions_dict=None,
    )
    features, names = cluster.read_window_feature_matrix(
        result,
        n_pca=1,
        autocorr_lags=(1,),
        density_windows=(("center", -2, 2),),
        require_nonzero_valid=False,
        min_valid_fraction=0.0,
    )
    assert features.shape[0] == 3
    assert "pca_0" in names
    assert "autocorr_1" in names
    assert "center" in names
    assert "global_mean" in names
    assert "iqr" in names


def test_plot_cluster_profiles_motif_index(monkeypatch):
    X = np.hstack([np.ones((3, 2)), np.zeros((3, 2))])
    labels = np.array([0, 1, 0])
    import matplotlib

    matplotlib.use("Agg")
    fig = cluster.plot_cluster_profiles(X, labels, motif_index=1, view_bp=2)
    assert fig is not None


def test_plot_region_cluster_profiles_motif_index(monkeypatch):
    X = np.hstack([np.ones((2, 2)), np.zeros((2, 2))])
    labels = np.array([0, 1])
    import matplotlib

    matplotlib.use("Agg")
    fig = cluster.plot_region_cluster_profiles(
        X, labels, motif_index=1, motif_count=2, window_bp=2
    )
    assert fig is not None


def test_read_window_feature_matrix_filters():
    data = np.array(
        [
            [0, 0, 1, 1],
            [1, 1, 0, 0],
        ],
        dtype=float,
    )
    val = np.array(
        [
            [0, 0, 1, 1],
            [0, 0, 0, 0],
        ],
        dtype=float,
    )
    result = cluster.ReadWindowExtractionResult(
        data_matrix=data,
        val_matrix=val,
        metadata=[],
        datasets=[],
        regions_dict=None,
    )
    features, _ = cluster.read_window_feature_matrix(
        result, require_nonzero_valid=True, min_valid_fraction=0.5, n_pca=0, autocorr_lags=()
    )
    # Second row should be dropped due to zero valid fraction
    assert features.shape[0] == 1


def test_plot_multisite_read_raster_smoke(monkeypatch):
    # Build a minimal ReadWindowExtractionResult with two motifs, two reads
    data = np.hstack(
        [
            np.tile(np.array([[0, 1, 0, 1]], dtype=float), (2, 1)),
            np.tile(np.array([[1, 0, 1, 0]], dtype=float), (2, 1)),
        ]
    )
    meta = [
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 6000,
            "region_end": 6004,
            "region_strand": "+",
        },
    ]
    r = cluster.ReadWindowExtractionResult(
        data_matrix=data,
        val_matrix=None,
        metadata=meta,
        datasets=[],
        regions_dict=None,
    )
    import matplotlib

    matplotlib.use("Agg")
    fig, stats = cluster.plot_multisite_read_raster(
        r, n_windows=2, min_separation_bp=5000, motif_index=0, smoothing=None
    )
    assert fig is not None
    assert stats["pairs"] >= 1


def test_cluster_read_windows_kmeans():
    rng = np.random.default_rng(0)
    feature_matrix = rng.normal(size=(20, 5))
    result = cluster.cluster_read_windows(
        feature_matrix, method="kmeans", n_clusters=3, random_state=0
    )
    assert isinstance(result.labels_raw, np.ndarray)
    assert result.labels_raw.shape[0] == 20
    assert result.labels_size_ordered.shape == result.labels_raw.shape
    assert "silhouette" in result.metrics


def test_cluster_label_mapping_round_trips_labels():
    labels_raw = np.array([3, 3, 7, 7, 7])
    labels_size_ordered = np.array([1, 1, 0, 0, 0])

    mapping = cluster.cluster_label_mapping(labels_raw, labels_size_ordered)
    remapped = cluster.apply_cluster_label_mapping(np.array([7, 3, 9]), mapping)

    assert mapping == {3: 1, 7: 0}
    np.testing.assert_array_equal(remapped, np.array([0, 1, -1]))
