from __future__ import annotations

import numpy as np
import pytest

from dimelo import footprint
from dimelo.footprint import BernoulliHMM


def test_viterbi_recovers_obvious_path():
    # state 0 protected (p=0.02), state 1 accessible (p=0.98)
    hmm = BernoulliHMM(p_emit=(0.02, 0.98), stay=0.9)
    obs = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0])
    path = hmm.viterbi(obs)
    assert list(path) == [1, 1, 0, 0, 0, 1, 1]


def test_viterbi_handles_missing_positions():
    hmm = BernoulliHMM(p_emit=(0.02, 0.98), stay=0.9)
    obs = np.array([1.0, np.nan, 0.0, 0.0, np.nan, 1.0])
    path = hmm.viterbi(obs)
    # missing positions are uninformative; observed 1->accessible, 0->protected
    assert path[0] == 1 and path[2] == 0 and path[3] == 0 and path[5] == 1


def _footprint_reads():
    # 15 footprinted reads: methylated flanks, protected (unmethylated) center 5..9
    footprinted = np.tile(
        np.array([1.0] * 5 + [0.0] * 5 + [1.0] * 5), (15, 1)
    )
    # 5 fully-accessible reads (no footprint)
    accessible = np.ones((5, 15))
    return np.vstack([footprinted, accessible])


def test_fit_orders_protected_state_and_recovers_emissions():
    reads = _footprint_reads()
    hmm = BernoulliHMM().fit([reads[i] for i in range(reads.shape[0])])
    # state 0 is protected (lower methylation) after fit
    assert hmm.p_emit[0] < hmm.p_emit[1]
    assert hmm.p_emit[0] < 0.2
    assert hmm.p_emit[1] > 0.8


def test_call_footprints_flags_footprinted_reads_only():
    reads = _footprint_reads()
    table, hmm = footprint.call_footprints(reads, anchor_index=7, min_width=3)
    table = table.set_index("read_index")

    # the 15 footprinted reads are flagged with a protected run over the center
    footprinted = table.loc[0:14]
    assert bool(footprinted["footprint_present"].all())
    assert (footprinted["protected_start"] == 5).all()
    assert (footprinted["protected_end"] == 10).all()
    assert (footprinted["protected_width"] == 5).all()
    assert (footprinted["footprint_score"] > 0.5).all()

    # the 5 accessible reads are not
    accessible = table.loc[15:19]
    assert not bool(accessible["footprint_present"].any())


def test_call_footprints_respects_width_bounds():
    reads = _footprint_reads()
    # require a wider footprint than the 5-position protected run -> none called
    table, _ = footprint.call_footprints(
        reads, anchor_index=7, min_width=8
    )
    assert not bool(table["footprint_present"].any())


def test_call_footprints_anchor_off_the_footprint():
    reads = _footprint_reads()
    # anchor at position 0 (in the accessible flank) -> protected center does not overlap
    table, _ = footprint.call_footprints(reads, anchor_index=0, min_width=3)
    assert not bool(table.loc[table["read_index"] < 15, "footprint_present"].any())


def test_aggregate_footprint_profile_dips_at_center():
    reads = _footprint_reads()
    profile = footprint.aggregate_footprint_profile(reads).set_index("position")
    # methylation is high at the flanks, low at the protected center (positions 5-9)
    assert profile.loc[0, "mean_methylation"] == pytest.approx(1.0)
    assert profile.loc[7, "mean_methylation"] < profile.loc[0, "mean_methylation"]
    # protected posterior peaks at the center
    assert (
        profile.loc[7, "mean_protected_posterior"]
        > profile.loc[0, "mean_protected_posterior"]
    )
    assert (profile["n_observed"] == reads.shape[0]).all()


def test_call_footprints_validates_inputs():
    with pytest.raises(ValueError, match="2D"):
        footprint.call_footprints(np.array([1.0, 0.0]), anchor_index=0)
    with pytest.raises(ValueError, match="anchor_index"):
        footprint.call_footprints(np.ones((2, 5)), anchor_index=9)
    empty, _ = footprint.call_footprints(np.empty((0, 5)), anchor_index=0)
    assert empty.empty


def test_footprint_profile_plotting():
    import matplotlib

    matplotlib.use("Agg")
    from dimelo import plotting, plotting_matplotlib

    reads = _footprint_reads()
    profile = footprint.aggregate_footprint_profile(reads)
    payload = plotting.prepare_footprint_profile_data(profile=profile, anchor_index=7)
    assert list(payload["profile_table"]["position"]) == list(range(15))
    fig, ax = plotting_matplotlib.plot_footprint_profile_matplotlib(
        payload, title="footprint"
    )
    assert fig is not None
