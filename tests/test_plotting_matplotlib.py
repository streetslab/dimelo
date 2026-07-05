from __future__ import annotations

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from dimelo import plotting_matplotlib as pm  # noqa: E402


def _new_fig():
    plt = pm._import_pyplot()
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    return fig


def test_save_figure_writes_file(tmp_path):
    fig = _new_fig()
    out = tmp_path / "nested" / "figure.png"
    returned = pm.save_figure(fig, out, dpi=72)
    assert returned == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_require_payload_keys_raises_on_missing():
    with pytest.raises(ValueError, match="missing required keys"):
        pm._require_payload_keys({"a": 1}, ("a", "b"), owner="owner")
    # all present -> no error
    pm._require_payload_keys({"a": 1, "b": 2}, ("a", "b"), owner="owner")


def test_require_payload_table_type_and_presence():
    table = pd.DataFrame({"x": [1]})
    assert pm._require_payload_table({"t": table}, "t") is table
    with pytest.raises(ValueError, match="missing required key"):
        pm._require_payload_table({}, "t")
    with pytest.raises(TypeError, match="must be a pandas DataFrame"):
        pm._require_payload_table({"t": [1, 2, 3]}, "t")


def test_make_axes_creates_and_reuses():
    plt = pm._import_pyplot()
    fig, axes = pm._make_axes(axes=None, n_axes=3)
    assert len(axes) == 3
    assert fig is not None

    existing_fig, existing_axes = plt.subplots(2, 1, squeeze=False)
    reuse_fig, reuse_axes = pm._make_axes(axes=existing_axes, n_axes=2)
    assert len(reuse_axes) == 2
    # the provided axes must actually be reused, not silently replaced
    assert reuse_fig is existing_fig
    assert reuse_axes[0] is existing_axes[0, 0]
    assert reuse_axes[1] is existing_axes[1, 0]


def test_make_axis_creates_and_reuses():
    plt = pm._import_pyplot()
    fig, ax = pm._make_axis(ax=None)
    assert fig is not None and ax is not None

    _, provided = plt.subplots()
    reuse_fig, reuse_ax = pm._make_axis(ax=provided)
    assert reuse_ax is provided


def test_plot_region_discovery_scan_renders_with_hits():
    payload = {
        "scan_table": pd.DataFrame(
            {
                "contig": ["chr1", "chr1", "chr2"],
                "window_midpoint": [100.0, 200.0, 50.0],
                "score_value": [0.2, 0.9, 0.5],
            }
        ),
        "hit_table": pd.DataFrame(
            {
                "contig": ["chr1"],
                "start": [150],
                "end": [250],
                "score_value": [0.9],
            }
        ),
        "metadata": {"contig_order": ["chr1", "chr2"], "score_column": "score_value"},
    }
    fig, axes = pm.plot_region_discovery_scan_matplotlib(payload, title="scan")
    assert fig is not None
    assert len(axes) == 2


def test_plot_region_discovery_scan_empty_contigs():
    payload = {
        "scan_table": pd.DataFrame(columns=["contig", "window_midpoint", "score_value"]),
        "hit_table": pd.DataFrame(columns=["contig", "start", "end", "score_value"]),
        "metadata": {"contig_order": [], "score_column": "score_value"},
    }
    fig, axes = pm.plot_region_discovery_scan_matplotlib(payload)
    assert fig is not None
    assert len(axes) == 1


def test_plot_region_discovery_scan_missing_key_raises():
    with pytest.raises(ValueError, match="missing required keys"):
        pm.plot_region_discovery_scan_matplotlib({"scan_table": pd.DataFrame()})
