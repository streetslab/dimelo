import numpy as np
import pandas as pd

from dimelo import plot_read_browser


def test_collapse_rows_preserves_expected_assignment():
    read_extent_df = pd.DataFrame(
        {
            "read_start": [0, 10, 30, 60],
            "read_end": [50, 20, 40, 70],
            "y_index": [0, 1, 2, 3],
            "read_name": ["r0", "r1", "r2", "r3"],
        }
    )

    collapsed = plot_read_browser.collapse_rows(
        read_extent_df,
        minimum_gap=0,
        meta_sort=None,
    )

    assert collapsed.to_dict() == {0: 0, 1: 1, 2: 1, 3: 0}


def test_make_browser_figure_batches_read_extents_into_single_trace():
    read_extent_df = pd.DataFrame(
        {
            "read_start": [100, 150, 200],
            "read_end": [130, 185, 240],
            "y_index": [0, 1, 2],
            "read_name": ["readA", "readB", "readC"],
        }
    )
    mod_event_df = pd.DataFrame(
        {
            "y_index": [0, 1, 2],
            "read_name": ["readA", "readB", "readC"],
            "motif": ["A,0", "A,0", "A,0"],
            "pos": [105, 160, 225],
            "prob": [0.4, 0.8, 0.6],
        }
    )

    fig = plot_read_browser.make_browser_figure(
        read_extent_df=read_extent_df,
        mod_event_df=mod_event_df,
        collapse=False,
        chrom="chr1",
        region_start=90,
        region_end=260,
        colorscales={"A,0": "Viridis"},
    )

    assert len(fig.data) == 2
    line_trace = fig.data[0]
    assert line_trace.mode == "lines"
    assert np.count_nonzero(np.array(line_trace.x, dtype=object) == None) == 3  # noqa: E711


def test_format_browser_data_preserves_event_vectors_and_duplicate_extent_handling():
    entry_labels = [
        "chromosome",
        "strand",
        "region_start",
        "region_end",
        "read_start",
        "read_end",
        "read_name",
        "motif",
        "mod_vector",
        "val_vector",
    ]
    read_tuples = [
        (
            "chr1",
            "+",
            90,
            120,
            100,
            110,
            "read1",
            "A,0",
            np.array([0.1, 0.8, 0.2]),
            np.array([0, 1, 1]),
        ),
        (
            "chr1",
            "+",
            90,
            120,
            102,
            109,
            "read1",
            "C,0",
            np.array([0.3, 0.4]),
            np.array([1, 0]),
        ),
        (
            "chr1",
            "-",
            180,
            210,
            200,
            204,
            "read2",
            "A,0",
            np.array([0.9, 0.7]),
            np.array([1, 0]),
        ),
    ]

    read_extent_df, mod_event_df = plot_read_browser.format_browser_data(
        read_tuples=read_tuples,
        entry_labels=entry_labels,
    )

    assert read_extent_df.reset_index(drop=True).to_dict("records") == [
        {"read_start": 100, "read_end": 110, "y_index": 0, "read_name": "read1"},
        {"read_start": 200, "read_end": 204, "y_index": 1, "read_name": "read2"},
    ]

    assert list(
        mod_event_df[["y_index", "read_name", "motif", "pos", "prob"]].itertuples(
            index=False, name=None
        )
    ) == [
        (0, "read1", "A,0", 101, 0.8),
        (0, "read1", "A,0", 102, 0.2),
        (0, "read1", "C,0", 102, 0.3),
        (1, "read2", "A,0", 200, 0.9),
    ]
