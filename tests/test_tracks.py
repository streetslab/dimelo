from __future__ import annotations

import pandas as pd
import pytest

from dimelo import tracks


def _regions():
    return pd.DataFrame(
        {"region_id": ["r1"], "chromosome": ["chr1"], "start": [0], "end": [100]}
    )


def test_bin_regions_tiles_exactly():
    binned = tracks.bin_regions(_regions(), bins=5)
    assert len(binned) == 5
    assert list(binned["bin_start"]) == [0, 20, 40, 60, 80]
    assert list(binned["bin_end"]) == [20, 40, 60, 80, 100]
    assert binned["bin_end"].iloc[-1] == 100  # last bin absorbs the end


def test_bin_regions_validates():
    with pytest.raises(ValueError, match="requires columns"):
        tracks.bin_regions(pd.DataFrame({"chromosome": ["chr1"]}), bins=2)
    with pytest.raises(ValueError, match="bins must be"):
        tracks.bin_regions(_regions(), bins=0)


def _write_bigwig(path):
    import pyBigWig

    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([("chr1", 100)])
    # value == genomic position, so per-bin means increase across the region
    bw.addEntries("chr1", 0, values=[float(i) for i in range(100)], span=1, step=1)
    bw.close()


def test_import_bigwig_signal_increasing_bins(tmp_path):
    bw_path = tmp_path / "cov.bw"
    _write_bigwig(bw_path)
    binned = tracks.import_bigwig_signal(bw_path, _regions(), bins=5)
    assert list(binned["bin"]) == [0, 1, 2, 3, 4]
    # bin 0 covers positions 0..19 -> mean 9.5, bin 1 -> 29.5, ...
    assert binned["signal"].iloc[0] == pytest.approx(9.5)
    assert binned["signal"].iloc[1] == pytest.approx(29.5)
    assert list(binned["signal"]) == sorted(binned["signal"])  # strictly increasing


def test_import_bigwig_signal_missing_contig(tmp_path):
    bw_path = tmp_path / "cov.bw"
    _write_bigwig(bw_path)
    other = pd.DataFrame(
        {"region_id": ["x"], "chromosome": ["chr2"], "start": [0], "end": [50]}
    )
    binned = tracks.import_bigwig_signal(bw_path, other, bins=2)
    assert binned["signal"].isna().all()


def test_import_bedgraph_signal_overlap_weighted(tmp_path):
    bg_path = tmp_path / "cov.bedgraph"
    bg_path.write_text("chr1\t0\t50\t1\nchr1\t50\t100\t5\n")
    binned = tracks.import_bedgraph_signal(bg_path, _regions(), bins=5)
    signal = list(binned["signal"])
    # bins 0,1 (0-40) fully in value-1 interval; bin 2 (40-60) straddles: (10*1+10*5)/20=3
    assert signal[0] == pytest.approx(1.0)
    assert signal[1] == pytest.approx(1.0)
    assert signal[2] == pytest.approx(3.0)
    assert signal[3] == pytest.approx(5.0)
    assert signal[4] == pytest.approx(5.0)


def test_correlate_binned_signals_positive(tmp_path):
    bw_path = tmp_path / "cov.bw"
    _write_bigwig(bw_path)
    external = tracks.import_bigwig_signal(bw_path, _regions(), bins=5)
    # DiMeLo-like binned signal that increases with the external track
    dimelo = external[["region_id", "bin"]].copy()
    dimelo["signal"] = [0.1, 0.3, 0.5, 0.7, 0.9]

    paired, stats = tracks.correlate_binned_signals(dimelo, external)
    assert stats["n"] == 5
    assert stats["pearson_r"] > 0.99
    assert stats["spearman_rho"] == pytest.approx(1.0)
    assert list(paired.columns) == ["region_id", "bin", "dimelo", "external"]


def test_correlate_binned_signals_too_few_points():
    dimelo = pd.DataFrame({"region_id": ["r1"], "bin": [0], "signal": [0.5]})
    external = pd.DataFrame({"region_id": ["r1"], "bin": [0], "signal": [1.0]})
    _, stats = tracks.correlate_binned_signals(dimelo, external)
    assert stats["n"] == 1
    assert pd.isna(stats["pearson_r"])


def test_track_correlation_plotting(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    from dimelo import plotting, plotting_matplotlib

    bw_path = tmp_path / "cov.bw"
    _write_bigwig(bw_path)
    external = tracks.import_bigwig_signal(bw_path, _regions(), bins=5)
    dimelo = external[["region_id", "bin"]].copy()
    dimelo["signal"] = [0.1, 0.3, 0.5, 0.7, 0.9]
    paired, stats = tracks.correlate_binned_signals(dimelo, external)

    payload = plotting.prepare_track_correlation_data(paired=paired, stats=stats)
    fig, _ = plotting_matplotlib.plot_track_correlation_matplotlib(payload, title="corr")
    assert fig is not None
