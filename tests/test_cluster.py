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
        cluster.utils,
        "regions_dict_from_input",
        lambda *args, **kwargs: fake_regions_dict,
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
        {
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 100,
            "region_strand": "+",
        },
        {
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 100,
            "region_strand": "+",
        },
        {
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 100,
            "region_strand": "+",
        },
        {
            "chromosome": "chr1",
            "region_start": 200,
            "region_end": 300,
            "region_strand": "-",
        },
        {
            "chromosome": "chr1",
            "region_start": 200,
            "region_end": 300,
            "region_strand": "-",
        },
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
        {
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 100,
            "region_strand": "+",
        },
        {
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 100,
            "region_strand": "+",
        },
        {
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 100,
            "region_strand": "+",
        },
        {
            "chromosome": "chr1",
            "region_start": 200,
            "region_end": 300,
            "region_strand": "-",
        },
        {
            "chromosome": "chr1",
            "region_start": 200,
            "region_end": 300,
            "region_strand": "-",
        },
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
        {
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 100,
            "region_strand": "+",
        },
        {
            "chromosome": "chr1",
            "region_start": 200,
            "region_end": 300,
            "region_strand": "-",
        },
        {
            "chromosome": "chr1",
            "region_start": 200,
            "region_end": 300,
            "region_strand": "-",
        },
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
        config=cluster.ReadWindowExtractionConfig(
            window_size=2,
            orientation_aware=True,
            enforce_thresholded_vectors=False,
        ),
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
            window_size=2,
            orientation_aware=False,
            filter_multi_region_reads=True,
            enforce_thresholded_vectors=False,
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
        result,
        require_nonzero_valid=True,
        min_valid_fraction=0.5,
        n_pca=0,
        autocorr_lags=(),
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
    assert stats["ml_score_cmaps"] == {"A,0": "Blues", "CG,0": "Oranges"}


def test_plot_multisite_read_raster_ml_scatter_uses_sequential_cmaps_for_other_motifs():
    import matplotlib

    matplotlib.use("Agg")
    motif_1 = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
        ],
        dtype=float,
    )
    motif_2 = np.array(
        [
            [1, 0, 1, 0],
            [0, 1, 0, 1],
        ],
        dtype=float,
    )
    motif_3 = np.array(
        [
            [0, 0, 1, 1],
            [1, 1, 0, 0],
        ],
        dtype=float,
    )
    data = np.hstack([motif_1, motif_2, motif_3])
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
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        motif_count=3,
        plot_all_motifs=True,
        motif_labels=["X,0", "Y,0", "Z,0"],
        smoothing=None,
    )
    assert fig is not None
    assert stats["ml_score_cmaps"] == {"X,0": "Greens", "Y,0": "Purples", "Z,0": "Reds"}


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


def test_plot_two_site_read_raster_wrapper():
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
            "read_length": 8000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 2000,
            "region_end": 2004,
            "region_strand": "+",
            "read_length": 8000,
        },
    ]
    r = cluster.ReadWindowExtractionResult(
        data_matrix=data,
        val_matrix=None,
        metadata=meta,
        datasets=[],
        regions_dict=None,
    )
    fig, stats = cluster.plot_two_site_read_raster(
        r,
        second_site_offset_bp=2000,
        window_width_bp=4,
        motif_index=0,
        smoothing=None,
    )
    assert fig is not None
    assert stats["n_windows"] == 2
    assert stats["window_centers_bp"] == [0.0, 2000.0]
    assert stats["window_widths_bp"] == [4, 4]


def test_plot_multisite_read_raster_supports_symmetric_center_spec():
    import matplotlib

    matplotlib.use("Agg")
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
            "read_length": 9000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 2000,
            "region_end": 2004,
            "region_strand": "+",
            "read_length": 9000,
        },
        {
            "read_name": "r1",
            "chromosome": "chr1",
            "region_start": 4000,
            "region_end": 4004,
            "region_strand": "+",
            "read_length": 9000,
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
        symmetric_side_windows=1,
        symmetric_max_offset_bp=2000,
        window_match_tolerance_bp=200,
        smoothing=None,
    )
    assert fig is not None
    assert stats["n_windows"] == 3
    assert stats["window_centers_bp"] == [-2000.0, 0.0, 2000.0]


def test_plot_multisite_read_raster_heatmap_window_titles():
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
        rotate=False,
        smoothing=None,
    )
    titles = [axis.get_title() for axis in fig.axes if axis.get_title()]
    assert "Window 1" in titles
    assert "Window 2" in titles


def test_plot_multisite_read_raster_accepts_centered_axis_mode():
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
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        x_axis_mode="centered",
        smoothing=None,
    )

    assert fig is not None
    assert stats["x_axis_mode"] == "centered"
    data_axes = [ax for ax in fig.axes if ax.lines]
    assert len(data_axes) == 2
    assert [float(ax.lines[0].get_xdata()[0]) for ax in data_axes] == [0.0, 0.0]
    assert all(ax.get_xlim()[0] < 0 < ax.get_xlim()[1] for ax in data_axes)


def test_plot_multisite_read_raster_heatmap_auto_uses_raw_values():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0, 1, 0, 0],
            [0, 0, 1, 0],
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
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        render_mode="heatmap",
        render_values="auto",
        axis_orientation="position_x",
        smoothing="gaussian",
        rotate=False,
    )
    assert fig is not None
    assert stats["render_values"] == "raw"


def test_plot_multisite_read_raster_heatmap_relative_axis_uses_window_offsets():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0.0, 0.8, 0.0, 0.0],
            [0.0, 0.0, 0.7, 0.0],
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
        render_values="smoothed",
        smoothing="gaussian",
        axis_orientation="position_x",
        x_axis_mode="relative_to_primary",
        rotate=True,
    )
    data_axes = [ax for ax in fig.axes if ax.images]
    assert len(data_axes) == 2
    center_0 = float(data_axes[0].lines[0].get_xdata()[0])
    center_1 = float(data_axes[1].lines[0].get_xdata()[0])
    assert abs(center_0) < 5
    assert center_1 > 1000


def test_plot_multisite_read_raster_heatmap_relative_axis_rejects_variable_offsets():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0.0, 0.8, 0.0, 0.0],
            [0.0, 0.0, 0.7, 0.0],
            [0.0, 0.6, 0.0, 0.0],
            [0.0, 0.0, 0.5, 0.0],
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
        {
            "read_name": "r2",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
        },
        {
            "read_name": "r2",
            "chromosome": "chr1",
            "region_start": 7000,
            "region_end": 7004,
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

    with pytest.raises(
        ValueError,
        match="render_mode='scatter'.*x_axis_mode='centered'.*fixed/explicit centers",
    ):
        cluster.plot_multisite_read_raster(
            r,
            n_windows=2,
            min_separation_bp=5000,
            motif_index=0,
            render_mode="heatmap",
            x_axis_mode="relative_to_primary",
            smoothing=None,
        )


def test_plot_multisite_read_raster_heatmap_smoothing_label_mentions_per_read_axis():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0.0, 0.8, 0.0, 0.0],
            [0.0, 0.0, 0.7, 0.0],
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
        render_values="smoothed",
        smoothing="gaussian",
        axis_orientation="position_x",
        rotate=True,
    )
    colorbar_axes = [
        ax
        for ax in fig.axes
        if "per-read smoothing along position axis" in ax.get_ylabel()
    ]
    assert colorbar_axes


def test_plot_multisite_read_raster_heatmap_thresholds_still_use_colorbar():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0.0, 0.9, 0.0, 0.0],
            [0.0, 0.0, 0.4, 0.0],
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
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        render_mode="heatmap",
        ml_score_thresholds=[0.5],
        smoothing=None,
    )

    assert fig is not None
    assert stats["render_mode"] == "heatmap"
    assert stats["ml_score_thresholds"] is None
    assert any("window signal heatmap" in ax.get_ylabel() for ax in fig.axes)
    assert not any(
        legend.get_title().get_text() == "Motif thresholds" for legend in fig.legends
    )


def test_plot_multisite_read_raster_accepts_axis_orientation_and_sort_modes():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
        ],
        dtype=float,
    )
    meta = [
        {
            "read_name": "z_read",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "z_read",
            "chromosome": "chr1",
            "region_start": 3000,
            "region_end": 3004,
            "region_strand": "+",
            "read_length": 7000,
        },
        {
            "read_name": "a_read",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
            "read_length": 6500,
        },
        {
            "read_name": "a_read",
            "chromosome": "chr1",
            "region_start": 3000,
            "region_end": 3004,
            "region_strand": "+",
            "read_length": 6500,
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
        n_windows=2,
        min_separation_bp=2000,
        motif_index=0,
        render_mode="scatter",
        axis_orientation="position_y",
        sort_by="read_name",
        sort_descending=False,
        smoothing=None,
    )
    assert fig is not None
    assert stats["axis_orientation"] == "position_y"
    assert stats["sort_by"] == "read_name"
    # At least one data axis should now expose reads on x.
    assert any(ax.get_xlabel() == "Reads (sorted)" for ax in fig.axes)


def test_plot_multisite_read_raster_default_scatter_uses_ml_score_colorbar():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0.0, 0.9, 0.0, 0.0],
            [0.0, 0.0, 0.4, 0.0],
        ],
        dtype=float,
    )
    val = np.array(
        [
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
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
        val_matrix=val,
        metadata=meta,
        datasets=[],
        regions_dict=None,
    )
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        render_mode="scatter",
        smoothing=None,
    )
    assert fig is not None
    assert stats["scatter_color_values"] == "ml_score"
    assert any("normalized ML score" in ax.get_ylabel() for ax in fig.axes)
    # Ensure ML coloring tracks mod values (0.9 / 0.4), not val-mask (1.0 / 1.0).
    scatter_arrays = [
        np.asarray(collection.get_array(), dtype=float)
        for axis in fig.axes
        for collection in getattr(axis, "collections", [])
        if collection.get_array() is not None and len(collection.get_array()) > 0
    ]
    assert any(np.isclose(arr, 0.9).any() for arr in scatter_arrays)
    assert any(np.isclose(arr, 0.4).any() for arr in scatter_arrays)


def test_plot_multisite_read_raster_uniform_downsampling_stats():
    import matplotlib

    matplotlib.use("Agg")
    data = np.tile(np.array([[0.0, 0.9, 0.0, 0.0]], dtype=float), (20, 1))
    meta = []
    for idx in range(20):
        meta.append(
            {
                "read_name": f"r{idx}",
                "chromosome": "chr1",
                "region_start": 0,
                "region_end": 4,
                "region_strand": "+",
            }
        )
        meta.append(
            {
                "read_name": f"r{idx}",
                "chromosome": "chr1",
                "region_start": 6000,
                "region_end": 6004,
                "region_strand": "+",
            }
        )
    r = cluster.ReadWindowExtractionResult(
        data_matrix=np.vstack([data, data]),
        val_matrix=None,
        metadata=meta,
        datasets=[],
        regions_dict=None,
    )
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        render_mode="scatter",
        max_rows=8,
        downsample_method="uniform",
        smoothing=None,
    )
    assert fig is not None
    assert stats["downsampled"] is True
    assert stats["downsample_method"] == "uniform"
    assert stats["rows_before_downsample"] > stats["rows_after_downsample"]


def test_plot_multisite_read_raster_drops_reads_that_do_not_span_display_windows():
    import matplotlib

    matplotlib.use("Agg")
    data = np.array(
        [
            [0.0, 1.0, 0.0, 0.0],  # long read, window 1
            [0.0, 1.0, 0.0, 0.0],  # long read, window 2
            [0.0, 1.0, 0.0, 0.0],  # short read, window 1
            [0.0, 1.0, 0.0, 0.0],  # short read, window 2 (outside read_end)
        ],
        dtype=float,
    )
    meta = [
        {
            "read_name": "r_long",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
            "read_start": -10,
            "read_end": 3000,
        },
        {
            "read_name": "r_long",
            "chromosome": "chr1",
            "region_start": 2000,
            "region_end": 2004,
            "region_strand": "+",
            "read_start": -10,
            "read_end": 3000,
        },
        {
            "read_name": "r_short",
            "chromosome": "chr1",
            "region_start": 0,
            "region_end": 4,
            "region_strand": "+",
            "read_start": -10,
            "read_end": 2001,
        },
        {
            "read_name": "r_short",
            "chromosome": "chr1",
            "region_start": 2000,
            "region_end": 2004,
            "region_strand": "+",
            "read_start": -10,
            "read_end": 2001,
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
        n_windows=2,
        min_separation_bp=1000,
        motif_index=0,
        window_widths_bp=4,
        enforce_full_window_span=True,
        smoothing=None,
    )
    assert fig is not None
    assert stats["pairs"] == 1
    assert stats["dropped_for_window_span"] >= 1


def test_plot_multisite_read_raster_ml_thresholds_use_motif_colors():
    import matplotlib

    matplotlib.use("Agg")
    data = np.hstack(
        [
            np.array(
                [
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                ],
                dtype=float,
            ),
            np.array(
                [
                    [0, 1, 0, 0],
                    [0, 1, 0, 0],
                ],
                dtype=float,
            ),
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
    fig, stats = cluster.plot_multisite_read_raster(
        r,
        n_windows=2,
        min_separation_bp=5000,
        motif_index=0,
        motif_count=2,
        plot_all_motifs=True,
        motif_labels=["A,0", "CG,0"],
        render_mode="scatter",
        ml_score_thresholds=[0.5, 0.75],
        motif_colors=["#1f77b4", "#ff7f0e"],
        smoothing=None,
    )
    assert fig is not None
    assert stats["ml_score_thresholds"] == {"A,0": 0.5, "CG,0": 0.75}
    assert stats["read_order_consistent"] is True
    assert any(
        legend.get_title().get_text() == "Motif thresholds" for legend in fig.legends
    )


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
        cluster.load_processed,
        "read_vectors_from_hdf5",
        lambda **kwargs: (records, dataset_names, None),
    )
    result = cluster.extract_read_windows(
        hdf5_file="reads.h5",
        motifs=["A,0"],
        config=cluster.ReadWindowExtractionConfig(
            window_size=None,
            orientation_aware=False,
            enforce_thresholded_vectors=False,
        ),
    )
    assert result.data_matrix.shape == (2, 6)


def test_extract_read_windows_auto_thresholds_raw_vectors_by_default(monkeypatch):
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
        np.array([0.0, 0.3, 0.8, 0.2, 0.9, 0.0], dtype=float),
        "A,0",
        106,
        "read1",
        100,
        "+",
        np.ones(6, dtype=float),
        101,
        105,
        "+",
        6,
    )

    monkeypatch.setattr(
        cluster.load_processed,
        "read_vectors_from_hdf5",
        lambda **kwargs: ([record], dataset_names, None),
    )
    with pytest.warns(RuntimeWarning, match="applying automatic thresholding"):
        result = cluster.extract_read_windows(
            hdf5_file="reads.h5",
            motifs=["A,0"],
            config=cluster.ReadWindowExtractionConfig(
                window_size=2, orientation_aware=False
            ),
        )

    # default 190/255 ~= 0.745 -> [0.3, 0.8, 0.2, 0.9] => [0,1,0,1]
    np.testing.assert_allclose(result.data_matrix, np.array([[0.0, 1.0, 0.0, 1.0]]))


def test_extract_read_windows_respects_adjustable_auto_threshold(monkeypatch):
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
        np.array([0.0, 0.45, 0.8, 0.2, 0.9, 0.0], dtype=float),
        "A,0",
        106,
        "read1",
        100,
        "+",
        np.ones(6, dtype=float),
        101,
        105,
        "+",
        6,
    )

    monkeypatch.setattr(
        cluster.load_processed,
        "read_vectors_from_hdf5",
        lambda **kwargs: ([record], dataset_names, None),
    )
    with pytest.warns(RuntimeWarning):
        result = cluster.extract_read_windows(
            hdf5_file="reads.h5",
            motifs=["A,0"],
            config=cluster.ReadWindowExtractionConfig(
                window_size=2,
                orientation_aware=False,
                auto_threshold_if_raw=100,  # 100/255 ~= 0.392
            ),
        )

    # threshold 100/255 -> [0.45,0.8,0.2,0.9] => [1,1,0,1]
    np.testing.assert_allclose(result.data_matrix, np.array([[1.0, 1.0, 0.0, 1.0]]))


def test_extract_read_windows_errors_when_raw_vectors_disallowed(monkeypatch):
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
        np.array([0.0, 0.45, 0.8, 0.2, 0.9, 0.0], dtype=float),
        "A,0",
        106,
        "read1",
        100,
        "+",
        np.ones(6, dtype=float),
        101,
        105,
        "+",
        6,
    )

    monkeypatch.setattr(
        cluster.load_processed,
        "read_vectors_from_hdf5",
        lambda **kwargs: ([record], dataset_names, None),
    )
    with pytest.raises(ValueError, match="unthresholded"):
        cluster.extract_read_windows(
            hdf5_file="reads.h5",
            motifs=["A,0"],
            config=cluster.ReadWindowExtractionConfig(
                window_size=2,
                orientation_aware=False,
                auto_threshold_if_raw=None,
                enforce_thresholded_vectors=True,
            ),
        )


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
        "chr1\t0\t10\talpha\t0\t+\nchr1\t20\t40\tbeta\t0\t+\nchr2\t5\t15\talpha\t0\t+\n"
    )
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t100\nchr2\t80\n")
    fig = cluster.plot_cluster_karyotype(bed, sizes)
    assert fig is not None


def test_plot_cluster_karyotype_uses_top_bottom_coordinate_labels_and_draws_rare_cluster_on_top(
    tmp_path,
):
    import matplotlib

    matplotlib.use("Agg")
    bed = tmp_path / "clusters_overlap.bed"
    bed.write_text(
        "chr1\t10\t20\t0\t0\t+\nchr1\t10\t20\t1\t0\t+\nchr2\t5\t15\t0\t0\t+\n"
    )
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t100\nchr2\t80\n")

    fig = cluster.plot_cluster_karyotype(
        bed,
        sizes,
        chromosome_order="natural",
        min_visible_bp=1,
        min_visible_fraction=0.0,
    )
    assert fig is not None
    ax = fig.axes[0]

    text_values = [t.get_text() for t in ax.texts]
    # One "0" label per chromosome at the top.
    assert text_values.count("0") >= 2
    # Bottom end-index labels should be present.
    assert "100" in text_values
    assert "80" in text_values
    top_zero_labels = [text for text in ax.texts if text.get_text() == "0"]
    assert len(top_zero_labels) >= 2
    # Top "0" labels should sit just above chromosome starts and centered on each line.
    for text in top_zero_labels:
        x_pos, y_pos = text.get_position()
        assert x_pos == pytest.approx(round(x_pos))
        assert y_pos < 0
        assert y_pos > -5

    end_labels = [text for text in ax.texts if text.get_text() in {"80", "100"}]
    assert len(end_labels) == 2
    for text in end_labels:
        x_pos, y_pos = text.get_position()
        end_val = float(text.get_text())
        # End labels should be just below the chromosome end and centered.
        assert x_pos == pytest.approx(round(x_pos))
        assert y_pos > end_val
        assert y_pos < end_val + 6

    legend = ax.get_legend()
    assert legend is not None
    label_to_color = {
        text.get_text(): handle.get_color()
        for text, handle in zip(legend.get_texts(), legend.legend_handles, strict=False)
    }
    assert "C0" in label_to_color
    assert "C1" in label_to_color

    overlap_lines = []
    for line in ax.lines:
        x_vals = np.asarray(line.get_xdata(), dtype=float)
        y_vals = np.asarray(line.get_ydata(), dtype=float)
        if x_vals.size != 2 or y_vals.size != 2:
            continue
        if np.allclose(x_vals, [0.0, 0.0]) and np.allclose(y_vals, [10.0, 20.0]):
            overlap_lines.append(line)

    assert len(overlap_lines) == 2
    top_line = max(overlap_lines, key=lambda line: line.get_zorder())
    assert tuple(top_line.get_color()) == tuple(label_to_color["C1"])
    assert not ax.yaxis.get_visible()


def test_infer_shared_window_size(monkeypatch):
    def fake_regions_dict_from_input(regions, window_size=None):
        if regions == "a.bed":
            return {"chr1": [(0, 100, "+"), (50, 120, "+")]}
        if regions == "b.bed":
            return {"chr1": [(0, 30, "+")], "chr2": [(10, 70, "-")]}
        raise ValueError("unexpected input")

    monkeypatch.setattr(
        cluster.utils, "regions_dict_from_input", fake_regions_dict_from_input
    )
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
