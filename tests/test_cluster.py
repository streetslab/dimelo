import numpy as np
import pandas as pd
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
            np.array([0.5, 0.5, 0.0], dtype=float),
            np.array([0.0, 0.5, 0.5], dtype=float),
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


def test_pileup_fraction_vector_from_bedmethyl(monkeypatch):
    monkeypatch.setattr(
        cluster.load_processed,
        "pileup_vectors_from_bedmethyl",
        lambda **kwargs: (
            np.array([1, 1, 0], dtype=float),
            np.array([2, 2, 1], dtype=float),
        ),
    )
    fraction = cluster._pileup_fraction_vector_from_bedmethyl(
        bedmethyl_file="pileup.bed.gz",
        motif="A,0",
        regions="chr1:0-3,+",
        pseudo_count=0.0,
    )
    np.testing.assert_allclose(fraction, np.array([0.5, 0.5, 0.0], dtype=float))


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


def test_summarize_read_cluster_region_associations_counts_and_fractions():
    metadata = [
        {"chromosome": "chr1", "region_start": 0, "region_end": 100, "region_strand": "+"},
        {"chromosome": "chr1", "region_start": 0, "region_end": 100, "region_strand": "+"},
        {"chromosome": "chr1", "region_start": 0, "region_end": 100, "region_strand": "+"},
        {"chromosome": "chr1", "region_start": 200, "region_end": 300, "region_strand": "-"},
        {"chromosome": "chr1", "region_start": 200, "region_end": 300, "region_strand": "-"},
    ]
    labels = [0, 0, 1, 1, 1]

    summary = cluster.summarize_read_cluster_region_associations(metadata, labels)

    region1 = summary[
        (summary["chrom"] == "chr1")
        & (summary["start"] == 0)
        & (summary["end"] == 100)
        & (summary["strand"] == "+")
    ].sort_values("cluster")
    assert region1["count"].tolist() == [2, 1]
    np.testing.assert_allclose(region1["fraction"].to_numpy(), np.array([2 / 3, 1 / 3]))
    assert region1["total_reads"].tolist() == [3, 3]


def test_summarize_read_cluster_region_associations_enrichment_and_multiple_testing():
    metadata = [
        {"chromosome": "chr1", "region_start": 0, "region_end": 100, "region_strand": "+"},
        {"chromosome": "chr1", "region_start": 0, "region_end": 100, "region_strand": "+"},
        {"chromosome": "chr1", "region_start": 0, "region_end": 100, "region_strand": "+"},
        {"chromosome": "chr1", "region_start": 200, "region_end": 300, "region_strand": "-"},
        {"chromosome": "chr1", "region_start": 200, "region_end": 300, "region_strand": "-"},
    ]
    labels = [0, 0, 1, 1, 1]

    summary = cluster.summarize_read_cluster_region_associations(metadata, labels)

    enriched = summary[
        (summary["chrom"] == "chr1")
        & (summary["start"] == 200)
        & (summary["end"] == 300)
        & (summary["strand"] == "-")
        & (summary["cluster"] == 1)
    ].iloc[0]
    assert enriched["log2_enrichment"] > 0
    assert 0.0 <= enriched["p_value"] <= 1.0
    assert 0.0 <= enriched["q_value"] <= 1.0

    assert set(["p_value", "q_value"]).issubset(summary.columns)
    assert summary["p_value"].between(0.0, 1.0).all()
    assert summary["q_value"].between(0.0, 1.0).all()


def test_summarize_read_cluster_region_associations_min_reads_filter():
    metadata = [
        {"chromosome": "chr1", "region_start": 0, "region_end": 100, "region_strand": "+"},
        {"chromosome": "chr1", "region_start": 200, "region_end": 300, "region_strand": "-"},
        {"chromosome": "chr1", "region_start": 200, "region_end": 300, "region_strand": "-"},
    ]
    labels = [0, 1, 1]

    summary = cluster.summarize_read_cluster_region_associations(
        metadata,
        labels,
        min_reads_per_region=2,
    )

    assert summary["start"].tolist() == [200, 200]
    assert summary["total_reads"].tolist() == [2, 2]


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
        config=cluster.ReadWindowExtractionConfig(window_size=2, orientation_aware=True),
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
            window_size=2, orientation_aware=False, filter_multi_region_reads=True
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
    fig = cluster.plot_cluster_profiles(X, labels, motif_index=1, view_window_size=1)
    assert fig is not None


def test_plot_region_cluster_profiles_motif_index(monkeypatch):
    X = np.hstack([np.ones((2, 2)), np.zeros((2, 2))])
    labels = np.array([0, 1])
    import matplotlib

    matplotlib.use("Agg")
    fig = cluster.plot_region_cluster_profiles(
        X, labels, motif_index=1, motif_count=2, window_size=2
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
    assert stats["render_mode"] == "scatter"


def test_plot_multisite_read_raster_supports_multimotif_grid():
    import matplotlib

    matplotlib.use("Agg")
    # 3 windows from one read, 2 concatenated motifs of width 4 each
    motif_a = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 1, 1, 0],
        ],
        dtype=float,
    )
    motif_c = np.array(
        [
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [1, 1, 0, 0],
        ],
        dtype=float,
    )
    data = np.hstack([motif_a, motif_c])
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
            "region_start": 2500,
            "region_end": 2504,
            "region_strand": "+",
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 5000,
            "region_end": 5004,
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
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=3,
        min_separation_bp=2000,
        motif_index=0,
        motif_count=2,
        plot_all_motifs=True,
        motif_labels=["A,0", "CG,0"],
        x_axis_mode="relative_to_primary",
        smoothing=None,
    )
    assert fig is not None
    assert stats["pairs"] >= 1
    assert stats["n_windows"] == 3
    assert stats["motifs_plotted"] == ["A,0", "CG,0"]


def test_plot_multisite_read_raster_supports_explicit_window_centers_and_length_filter():
    import matplotlib

    matplotlib.use("Agg")
    motif_a = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 1, 1, 0],
        ],
        dtype=float,
    )
    data = motif_a
    meta = [
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 2500,
            "region_end": 2504,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 5000,
            "region_end": 5004,
            "region_strand": "+",
            "read_length": 7000,
        },
    ]
    r = cluster.ReadWindowExtractionResult(
        data_matrix=data,
        val_matrix=None,
        metadata=meta,
        datasets=[],
        regions_dict=None,
    )
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        motif_index=0,
        window_center_offsets=[-2500, 0, 2500],
        x_axis_mode="relative_to_primary",
        smoothing=None,
    )
    assert fig is not None
    assert stats["n_windows"] == 3
    assert stats["window_center_offsets"] == [-2500.0, 0.0, 2500.0]
    assert stats["required_read_length_bp"] > 0
    assert stats["min_read_length_bp"] == stats["required_read_length_bp"]


def test_plot_multisite_read_raster_heatmap_centered_titles():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
        ],
        dtype=float,
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
    fig, _ = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        render_mode="heatmap",
        x_axis_mode="centered",
        rotate=False,
        smoothing=None,
    )
    titles = [axis.get_title() for axis in fig.axes if axis.get_title()]
    assert "Window 1" in titles
    assert "Window 2" in titles


def test_plot_multisite_read_raster_multimotif_nonrotated_uses_effective_window_count():
    import matplotlib

    matplotlib.use("Agg")
    motif_a = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 1, 1, 0],
        ],
        dtype=float,
    )
    motif_c = np.array(
        [
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [1, 1, 0, 0],
        ],
        dtype=float,
    )
    data = np.hstack([motif_a, motif_c])
    meta = [
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 2500,
            "region_end": 2504,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 5000,
            "region_end": 5004,
            "region_strand": "+",
            "read_length": 7000,
        },
    ]
    r = cluster.ReadWindowExtractionResult(
        data_matrix=data,
        val_matrix=None,
        metadata=meta,
        datasets=[],
        regions_dict=None,
    )
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        motif_index=0,
        motif_count=2,
        plot_all_motifs=True,
        motif_labels=["A,0", "CG,0"],
        window_center_offsets=[-2500, 0, 2500],
        center_tolerance_bp=500,
        rotate=False,
        render_mode="scatter",
        x_axis_mode="relative_to_primary",
        smoothing=None,
    )
    assert stats["n_windows"] == 3
    # 2 motifs x 3 windows + 1 colorbar axis
    assert len(fig.axes) == 7


def test_plot_multisite_read_raster_rejects_too_small_min_read_length():
    data = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 1, 1, 0],
        ],
        dtype=float,
    )
    meta = [
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 2500,
            "region_end": 2504,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 5000,
            "region_end": 5004,
            "region_strand": "+",
            "read_length": 7000,
        },
    ]
    r = cluster.ReadWindowExtractionResult(
        data_matrix=data,
        val_matrix=None,
        metadata=meta,
        datasets=[],
        regions_dict=None,
    )
    with pytest.raises(ValueError, match="min_read_length_bp"):
        cluster.plot_multisite_read_raster(
            r,
            motif_index=0,
            window_center_offsets=[-2500, 0, 2500],
            min_read_length_bp=100,
            smoothing=None,
        )


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


def test_extract_read_windows_infers_shortest_region_window(monkeypatch):
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
    records = [
        (
            "chr1",
            np.arange(20, dtype=float),
            "A,0",
            120,
            "r1",
            100,
            "+",
            np.ones(20, dtype=float),
            104,
            110,  # length 6 -> shortest
            "+",
            20,
        ),
        (
            "chr1",
            np.arange(20, dtype=float),
            "A,0",
            120,
            "r2",
            100,
            "+",
            np.ones(20, dtype=float),
            102,
            112,  # length 10
            "+",
            20,
        ),
    ]
    monkeypatch.setattr(
        cluster.load_processed, "read_vectors_from_hdf5", lambda **kwargs: (records, dataset_names, None)
    )
    result = cluster.extract_read_windows(
        hdf5_file="reads.h5",
        motifs=["A,0"],
        config=cluster.ReadWindowExtractionConfig(window_size=None, orientation_aware=False),
    )
    assert result.data_matrix.shape == (2, 6)


def test_classify_read_features_binary_includes_row_index():
    X = np.vstack([np.zeros((8, 3)), np.ones((8, 3))])
    labels = np.array(["a"] * 8 + ["b"] * 8)
    result = cluster.classify_read_features_binary(
        X,
        labels,
        classifier="logreg",
        random_state=0,
    )
    preds = result["predictions"]
    assert "row_index" in preds.columns
    assert preds["row_index"].between(0, X.shape[0] - 1).all()
    assert preds["row_index"].nunique() == X.shape[0]


def test_plot_classification_profiles_confusion_smoke():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 0, 1, 1],
            [1, 1, 0, 0],
        ],
        dtype=float,
    )
    preds = pd.DataFrame(
        {
            "true_label": ["neg", "neg", "pos", "pos"],
            "pred_label": ["neg", "pos", "neg", "pos"],
            "proba": [0.1, 0.7, 0.4, 0.8],
            "split": ["test"] * 4,
            "sample_label": ["s"] * 4,
            "row_index": [0, 1, 2, 3],
        }
    )
    fig = cluster.plot_classification_profiles(
        data,
        preds,
        split="test",
        group_by="confusion",
        show_valid_profile=False,
    )
    assert fig is not None


def test_plot_region_classification_profiles_smoke():
    import matplotlib

    matplotlib.use("Agg")
    X = np.array(
        [
            [0.1, 0.2, 0.1, 0.2],
            [0.8, 0.7, 0.8, 0.7],
            [0.2, 0.1, 0.2, 0.1],
            [0.7, 0.8, 0.7, 0.8],
        ],
        dtype=float,
    )
    preds = pd.DataFrame(
        {
            "true_label": ["neg", "pos", "neg", "pos"],
            "pred_label": ["neg", "pos", "pos", "neg"],
            "split": ["test"] * 4,
            "row_index": [0, 1, 2, 3],
        }
    )
    fig = cluster.plot_region_classification_profiles(
        X,
        preds,
        split="test",
        group_by="confusion",
    )
    assert fig is not None


def test_plot_classification_profiles_accepts_split_aligned_matrix_with_global_row_index():
    import matplotlib

    matplotlib.use("Agg")
    full_matrix = np.array(
        [
            [0.0, 0.1, 0.0, 0.1],
            [0.9, 0.8, 0.9, 0.8],
            [0.2, 0.1, 0.2, 0.1],
            [0.8, 0.9, 0.8, 0.9],
            [0.1, 0.2, 0.1, 0.2],
            [0.7, 0.6, 0.7, 0.6],
        ],
        dtype=float,
    )
    split_rows = np.array([5, 2, 4], dtype=int)
    split_matrix = full_matrix[split_rows]
    preds = pd.DataFrame(
        {
            "true_label": ["pos", "neg", "neg"],
            "pred_label": ["pos", "neg", "pos"],
            "proba": [0.9, 0.2, 0.7],
            "split": ["test"] * 3,
            "sample_label": ["s"] * 3,
            # Global row indices do not fit split_matrix bounds.
            "row_index": split_rows,
        }
    )
    fig = cluster.plot_classification_profiles(
        split_matrix,
        preds,
        split="test",
        group_by="pred_label",
    )
    assert fig is not None


def test_plot_region_classification_profiles_accepts_split_aligned_matrix_with_global_row_index():
    import matplotlib

    matplotlib.use("Agg")
    full_matrix = np.array(
        [
            [0.1, 0.2, 0.1, 0.2],
            [0.8, 0.7, 0.8, 0.7],
            [0.2, 0.1, 0.2, 0.1],
            [0.7, 0.8, 0.7, 0.8],
            [0.3, 0.2, 0.3, 0.2],
            [0.6, 0.5, 0.6, 0.5],
        ],
        dtype=float,
    )
    split_rows = np.array([4, 5, 2], dtype=int)
    split_matrix = full_matrix[split_rows]
    preds = pd.DataFrame(
        {
            "true_label": ["neg", "pos", "neg"],
            "pred_label": ["neg", "pos", "pos"],
            "split": ["test"] * 3,
            "row_index": split_rows,
        }
    )
    fig = cluster.plot_region_classification_profiles(
        split_matrix,
        preds,
        split="test",
        group_by="confusion",
    )
    assert fig is not None


def test_plot_cluster_karyotype_handles_non_numeric_cluster_names(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    bed = tmp_path / "clusters.bed"
    bed.write_text(
        "chr1\t0\t10\talpha\t0\t+\n"
        "chr1\t20\t40\tbeta\t0\t+\n"
        "chr2\t5\t15\talpha\t0\t+\n"
    )
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t100\nchr2\t80\n")
    fig = cluster.plot_cluster_karyotype(bed, sizes)
    assert fig is not None


def test_infer_shared_window_size(monkeypatch):
    def fake_regions_dict_from_input(regions, window_size=None):
        if regions == "a.bed":
            return {"chr1": [(0, 100, "+"), (50, 120, "+")]}
        if regions == "b.bed":
            return {"chr1": [(0, 30, "+")], "chr2": [(10, 70, "-")]}
        raise ValueError("unexpected input")

    monkeypatch.setattr(cluster.utils, "regions_dict_from_input", fake_regions_dict_from_input)
    shared = cluster._infer_shared_window_size(["a.bed", "b.bed"])
    assert shared == 15


def test_merge_read_window_results_center_crop():
    r1 = cluster.ReadWindowExtractionResult(
        data_matrix=np.arange(24, dtype=float).reshape(3, 8),
        val_matrix=np.ones((3, 8), dtype=float),
        metadata=[{"read_name": "r1"}, {"read_name": "r2"}, {"read_name": "r3"}],
        datasets=["read_name"],
        regions_dict=None,
    )
    r2 = cluster.ReadWindowExtractionResult(
        data_matrix=np.arange(10, dtype=float).reshape(2, 5),
        val_matrix=np.ones((2, 5), dtype=float),
        metadata=[{"read_name": "r4"}, {"read_name": "r5"}],
        datasets=["read_name"],
        regions_dict=None,
    )
    merged = cluster.merge_read_window_results(
        [r1, r2],
        source_labels=["on", "off"],
        align="center_crop",
    )
    assert merged.data_matrix.shape == (5, 5)
    assert merged.val_matrix is not None
    assert merged.val_matrix.shape == (5, 5)
    assert {m["source_label"] for m in merged.metadata} == {"on", "off"}


def test_merge_read_window_results_error_on_mismatch():
    r1 = cluster.ReadWindowExtractionResult(
        data_matrix=np.zeros((2, 4), dtype=float),
        val_matrix=None,
        metadata=[],
        datasets=[],
        regions_dict=None,
    )
    r2 = cluster.ReadWindowExtractionResult(
        data_matrix=np.zeros((2, 6), dtype=float),
        val_matrix=None,
        metadata=[],
        datasets=[],
        regions_dict=None,
    )
    with pytest.raises(ValueError, match="do not match"):
        cluster.merge_read_window_results([r1, r2], align="error")


def test_extract_read_windows_rejects_nonpositive_window_size(monkeypatch):
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
    records = [
        (
            "chr1",
            np.arange(20, dtype=float),
            "A,0",
            120,
            "r1",
            100,
            "+",
            np.ones(20, dtype=float),
            104,
            112,
            "+",
            20,
        ),
    ]
    monkeypatch.setattr(
        cluster.load_processed,
        "read_vectors_from_hdf5",
        lambda **kwargs: (records, dataset_names, None),
    )
    with pytest.raises(ValueError, match="window_size must be a positive integer"):
        cluster.extract_read_windows(
            hdf5_file="reads.h5",
            motifs=["A,0"],
            window_size=0,
        )
