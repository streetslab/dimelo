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


def test_bin_regions_absorbs_remainder():
    # non-divisible region: fractional linspace edges truncate to ints, last bin absorbs
    binned = tracks.bin_regions(_regions(), bins=3)
    assert len(binned) == 3
    assert list(binned["bin_start"]) == [0, 33, 66]
    assert list(binned["bin_end"]) == [33, 66, 100]
    # bins tile contiguously with no gaps/overlaps
    assert list(binned["bin_start"])[1:] == list(binned["bin_end"])[:-1]


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


def test_import_bigwig_signal_zero_width_bins_are_nan(tmp_path):
    # more bins than the region is long -> some integer edges collide (zero-width bins).
    # These must yield NaN, not crash pyBigWig.stats with 'Invalid interval bounds!'.
    bw_path = tmp_path / "cov.bw"
    _write_bigwig(bw_path)
    short = pd.DataFrame(
        {"region_id": ["s"], "chromosome": ["chr1"], "start": [0], "end": [3]}
    )
    binned = tracks.import_bigwig_signal(bw_path, short, bins=5)
    assert len(binned) == 5
    zero_width = binned["bin_end"] <= binned["bin_start"]
    assert zero_width.any()  # e.g. bins (0,0) and (1,1)
    assert binned.loc[zero_width, "signal"].isna().all()
    assert binned.loc[~zero_width, "signal"].notna().any()


def test_import_bedgraph_signal_overlap_weighted(tmp_path):
    bg_path = tmp_path / "cov.bedgraph"
    bg_path.write_text("chr1\t0\t50\t1\nchr1\t50\t100\t5\n")
    # region [0,90] in 3 bins -> edges [0,30,60,90]
    regions = pd.DataFrame(
        {"region_id": ["r"], "chromosome": ["chr1"], "start": [0], "end": [90]}
    )
    binned = tracks.import_bedgraph_signal(bg_path, regions, bins=3)
    signal = list(binned["signal"])
    assert signal[0] == pytest.approx(1.0)  # bin 0-30 fully value 1
    # bin 30-60 straddles ASYMMETRICALLY: 20bp of value 1 + 10bp of value 5 -> 70/30,
    # which a plain (unweighted) mean of {1,5}=3.0 would NOT produce.
    assert signal[1] == pytest.approx(70 / 30)
    assert signal[1] != pytest.approx(3.0)
    assert signal[2] == pytest.approx(5.0)  # bin 60-90 fully value 5


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
    fig, ax = plotting_matplotlib.plot_track_correlation_matplotlib(payload, title="corr")
    assert ax.get_xlabel() == "DiMeLo signal"
    assert ax.get_ylabel() == "external signal"
    assert ax.collections[0].get_offsets().shape[0] == 5  # 5 scatter points
    annotation = " ".join(text.get_text() for text in ax.texts)
    assert "Pearson r = 1.00" in annotation


def test_track_correlation_plotting_empty_and_nan():
    import matplotlib

    matplotlib.use("Agg")
    from dimelo import plotting, plotting_matplotlib

    # empty paired input -> early return, no scatter, no crash
    empty = plotting.prepare_track_correlation_data(
        paired=pd.DataFrame({"dimelo": [], "external": []}), stats={"n": 0}
    )
    fig, ax = plotting_matplotlib.plot_track_correlation_matplotlib(empty)
    assert len(ax.collections) == 0

    # NaN correlation must not be rendered as the literal 'nan'
    nan_payload = plotting.prepare_track_correlation_data(
        paired=pd.DataFrame({"dimelo": [0.1, 0.2], "external": [0.3, 0.4]}),
        stats={"pearson_r": float("nan"), "spearman_rho": float("nan")},
    )
    _, ax2 = plotting_matplotlib.plot_track_correlation_matplotlib(nan_payload)
    annotation = " ".join(text.get_text() for text in ax2.texts)
    assert "nan" not in annotation.lower()


# --------------------------------------------------------------------------- #
# S6b: Hi-C (cooler)                                                          #
# --------------------------------------------------------------------------- #


def _write_cool(path):
    import cooler

    bins = pd.DataFrame(
        {
            "chrom": ["chr1"] * 10,
            "start": list(range(0, 100, 10)),
            "end": list(range(10, 101, 10)),
        }
    )
    # high contact between bins {0,1}x{2,3}, medium {2,3}x{8,9}, low {0,1}x{8,9}
    pix = []
    for b1 in (0, 1):
        for b2 in (2, 3):
            pix.append((b1, b2, 100))
        for b2 in (8, 9):
            pix.append((b1, b2, 1))
    for b1 in (2, 3):
        for b2 in (8, 9):
            pix.append((b1, b2, 10))
    pixels = (
        pd.DataFrame(pix, columns=["bin1_id", "bin2_id", "count"])
        .sort_values(["bin1_id", "bin2_id"])
        .reset_index(drop=True)
    )
    cooler.create_cooler(str(path), bins, pixels)


def test_import_hic_contacts(tmp_path):
    pytest.importorskip("cooler")
    cool_path = tmp_path / "test.cool"
    _write_cool(cool_path)
    pairs = pd.DataFrame(
        {
            "pair_id": ["near", "mid", "far"],
            "region_a": ["chr1:0-20", "chr1:20-40", "chr1:0-20"],
            "region_b": ["chr1:20-40", "chr1:80-100", "chr1:80-100"],
        }
    )
    contacts = tracks.import_hic_contacts(cool_path, pairs).set_index("pair_id")
    # {0,1}x{2,3} = 4*100 ; {2,3}x{8,9} = 4*10 ; {0,1}x{8,9} = 4*1
    assert contacts.loc["near", "hic_contact"] == pytest.approx(400.0)
    assert contacts.loc["mid", "hic_contact"] == pytest.approx(40.0)
    assert contacts.loc["far", "hic_contact"] == pytest.approx(4.0)


def test_correlate_hic_vs_joint_occupancy(tmp_path):
    pytest.importorskip("cooler")
    cool_path = tmp_path / "test.cool"
    _write_cool(cool_path)
    pairs = pd.DataFrame(
        {
            "pair_id": ["near", "mid", "far"],
            "region_a": ["chr1:0-20", "chr1:20-40", "chr1:0-20"],
            "region_b": ["chr1:20-40", "chr1:80-100", "chr1:80-100"],
        }
    )
    contacts = tracks.import_hic_contacts(cool_path, pairs)
    joint = pd.DataFrame(
        {"pair_id": ["near", "mid", "far"], "log2_obs_exp": [2.0, 1.0, 0.1]}
    )
    paired, stats = tracks.correlate_hic_vs_joint_occupancy(contacts, joint)
    assert stats["n"] == 3
    # Hi-C contact (400,40,4) and joint occupancy (2.0,1.0,0.1) are rank-concordant
    assert stats["spearman_rho"] == pytest.approx(1.0)
    assert list(paired.columns) == ["pair_id", "hic_contact", "joint_occupancy"]


def test_correlate_hic_vs_joint_occupancy_validates():
    with pytest.raises(ValueError, match="hic_contacts requires"):
        tracks.correlate_hic_vs_joint_occupancy(
            pd.DataFrame({"pair_id": ["p"]}),
            pd.DataFrame({"pair_id": ["p"], "log2_obs_exp": [1.0]}),
        )
