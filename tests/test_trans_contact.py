from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dimelo import trans_contact


def test_estimate_cis_reach_half_decay():
    # density peaks at offset 0 (1.0) and decays to background 0.0 by offset 100
    profile = pd.DataFrame(
        {
            "offset_bp": [0, 25, 50, 75, 100],
            "mean_density": [1.0, 0.8, 0.4, 0.1, 0.0],
        }
    )
    # peak 1.0, background 0.0, decay_fraction 0.5 -> threshold 0.5; first offset below 0.5
    # is 50 (density 0.4)
    reach = trans_contact.estimate_cis_reach(profile)
    assert reach == 50


def test_estimate_cis_reach_never_decays_returns_max_offset():
    profile = pd.DataFrame(
        {"offset_bp": [0, 50, 100], "mean_density": [1.0, 1.0, 1.0]}
    )
    assert trans_contact.estimate_cis_reach(profile) == 100


def _positions():
    return pd.DataFrame(
        {
            "chromosome": ["chr1", "chr1", "chr1", "chr2"],
            "position": [1000, 1050, 5000, 1000],
            "signal": [0.9, 0.8, 0.7, 0.6],
            # per-position max Hi-C contact to a bound site
            "hic_contact": [10.0, 10.0, 90.0, 5.0],
        }
    )


def _bound_sites():
    return pd.DataFrame({"chromosome": ["chr1"], "position": [1000]})


def test_classify_cis_trans_assigns_three_classes():
    classified = trans_contact.classify_cis_trans(
        _positions(),
        _bound_sites(),
        cis_reach_bp=100,
        hic_contact_threshold=50.0,
    )

    def _row(chrom, pos):
        return classified[
            (classified["chromosome"] == chrom) & (classified["position"] == pos)
        ].iloc[0]

    # within 100 bp of the bound site at chr1:1000 -> cis
    assert _row("chr1", 1000)["classification"] == trans_contact.CIS
    assert _row("chr1", 1050)["classification"] == trans_contact.CIS
    # chr1:5000 is distal (4000 bp) with high Hi-C contact (90 >= 50) -> trans
    assert _row("chr1", 5000)["classification"] == trans_contact.TRANS
    # chr2:1000 has no bound site on chr2 -> distal, low Hi-C -> background
    chr2_row = _row("chr2", 1000)
    assert np.isinf(chr2_row["nearest_bound_distance"])
    assert chr2_row["classification"] == trans_contact.DISTAL_BACKGROUND


def test_classify_cis_trans_without_hic_has_no_trans():
    classified = trans_contact.classify_cis_trans(
        _positions(), _bound_sites(), cis_reach_bp=100
    )
    assert (classified["classification"] != trans_contact.TRANS).all()
    # the distal chr1:5000 becomes background rather than trans
    distal = classified[classified["position"] == 5000].iloc[0]
    assert distal["classification"] == trans_contact.DISTAL_BACKGROUND


def test_classify_cis_trans_requires_hic_column_when_thresholded():
    positions = _positions().drop(columns=["hic_contact"])
    with pytest.raises(ValueError, match="hic_contact"):
        trans_contact.classify_cis_trans(
            positions, _bound_sites(), cis_reach_bp=100, hic_contact_threshold=50.0
        )


def test_apply_trans_correction_splits_signal():
    classified = trans_contact.classify_cis_trans(
        _positions(), _bound_sites(), cis_reach_bp=100, hic_contact_threshold=50.0
    )
    corrected = trans_contact.apply_trans_correction(classified)

    def _row(chrom, pos):
        return corrected[
            (corrected["chromosome"] == chrom) & (corrected["position"] == pos)
        ].iloc[0]

    # a cis position keeps its signal in-site, contributes nothing to trans
    cis = _row("chr1", 1050)
    assert cis["in_site_signal"] == pytest.approx(0.8)
    assert cis["trans_signal"] == pytest.approx(0.0)
    # the trans position (chr1:5000) moves to the trans track only
    trans = _row("chr1", 5000)
    assert trans["in_site_signal"] == pytest.approx(0.0)
    assert trans["trans_signal"] == pytest.approx(0.7)


def test_summarize_trans_correction():
    classified = trans_contact.classify_cis_trans(
        _positions(), _bound_sites(), cis_reach_bp=100, hic_contact_threshold=50.0
    )
    summary = trans_contact.apply_trans_correction(classified)
    result = trans_contact.summarize_trans_correction(summary).iloc[0]
    assert result["n_cis"] == 2
    assert result["n_trans"] == 1
    assert result["n_distal_background"] == 1
    assert result["total_in_site_signal"] == pytest.approx(0.9 + 0.8)
    assert result["total_trans_signal"] == pytest.approx(0.7)
    # trans_fraction = 0.7 / (1.7 + 0.7)
    assert result["trans_fraction"] == pytest.approx(0.7 / (1.7 + 0.7))
