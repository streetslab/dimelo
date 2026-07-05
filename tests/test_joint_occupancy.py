from __future__ import annotations

import pandas as pd
import pytest

from dimelo import joint_occupancy


def _read(region_id, read_id, is_sig, posterior):
    return {
        "region_id": region_id,
        "read_id": read_id,
        "is_true_signal": is_sig,
        "occupancy_posterior": posterior,
    }


def _called_reads():
    rows = []
    # Coupled pair (siteA, siteB): reads are mostly both-bound or both-unbound.
    for i in range(12):  # both
        rows.append(_read("siteA", f"c{i}", True, 0.9))
        rows.append(_read("siteB", f"c{i}", True, 0.9))
    for i in range(12, 24):  # neither
        rows.append(_read("siteA", f"c{i}", False, 0.1))
        rows.append(_read("siteB", f"c{i}", False, 0.1))
    for i in range(24, 27):  # a only
        rows.append(_read("siteA", f"c{i}", True, 0.9))
        rows.append(_read("siteB", f"c{i}", False, 0.1))
    for i in range(27, 30):  # b only
        rows.append(_read("siteA", f"c{i}", False, 0.1))
        rows.append(_read("siteB", f"c{i}", True, 0.9))
    # a read seen only at siteA (does not span both) must be ignored
    rows.append(_read("siteA", "lonely", True, 0.9))

    # Independent pair (siteC, siteD): 4 of each joint state -> OR ~ 1.
    states = [(True, True), (True, False), (False, True), (False, False)]
    idx = 0
    for sa, sb in states:
        for _ in range(4):
            rows.append(_read("siteC", f"i{idx}", sa, 0.9 if sa else 0.1))
            rows.append(_read("siteD", f"i{idx}", sb, 0.9 if sb else 0.1))
            idx += 1

    # Asymmetric coupled pair (siteE, siteF): n_a_only (5) != n_b_only (1), so an
    # a_only/b_only label or contingency transpose is detectable.
    asym = (
        [(True, True)] * 10
        + [(True, False)] * 5
        + [(False, True)] * 1
        + [(False, False)] * 4
    )
    for j, (sa, sb) in enumerate(asym):
        rows.append(_read("siteE", f"e{j}", sa, 0.9 if sa else 0.1))
        rows.append(_read("siteF", f"e{j}", sb, 0.9 if sb else 0.1))
    return pd.DataFrame(rows)


def _pairs():
    return pd.DataFrame(
        {
            "pair_id": ["coupled", "independent", "asymmetric"],
            "site_a": ["siteA", "siteC", "siteE"],
            "site_b": ["siteB", "siteD", "siteF"],
            "distance_bp": [1000, 50000, 2000],
        }
    )


def test_joint_occupancy_detects_positive_coupling():
    pair_summary, per_read = joint_occupancy.joint_occupancy(
        called_reads=_called_reads(), pairs=_pairs()
    )
    pair_summary = pair_summary.set_index("pair_id")

    coupled = pair_summary.loc["coupled"]
    assert coupled["n_reads_spanning_both"] == 30  # 'lonely' read excluded
    assert coupled["n_both"] == 12
    assert coupled["odds_ratio"] > 1.0
    assert coupled["log2_obs_exp"] > 0.0
    assert bool(coupled["significant"]) is True
    assert coupled["posterior_corr"] > 0.5
    assert coupled["distance_bp"] == 1000

    independent = pair_summary.loc["independent"]
    assert independent["n_reads_spanning_both"] == 16
    assert independent["odds_ratio"] == pytest.approx(1.0)
    assert bool(independent["significant"]) is False

    # per-read joint states are emitted for spanning reads only
    coupled_reads = per_read[per_read["pair_id"] == "coupled"]
    assert len(coupled_reads) == 30
    assert set(coupled_reads["joint_state"]) <= {"both", "a_only", "b_only", "neither"}
    assert (coupled_reads["joint_state"] == "both").sum() == 12
    assert "lonely" not in set(coupled_reads["read_id"])


def test_joint_occupancy_asymmetric_contingency_orientation():
    # asymmetric pair pins a_only (5) vs b_only (1) so a contingency/label transpose is
    # caught (a symmetric fixture would leave all four cells and the OR unchanged).
    pair_summary, per_read = joint_occupancy.joint_occupancy(
        called_reads=_called_reads(), pairs=_pairs()
    )
    asym = pair_summary.set_index("pair_id").loc["asymmetric"]
    assert asym["n_both"] == 10
    assert asym["n_a_only"] == 5
    assert asym["n_b_only"] == 1
    assert asym["n_neither"] == 4

    asym_reads = per_read[per_read["pair_id"] == "asymmetric"]
    assert (asym_reads["joint_state"] == "a_only").sum() == 5
    assert (asym_reads["joint_state"] == "b_only").sum() == 1


def test_joint_occupancy_posterior_corr_ignores_nan_reads():
    # a single read with a NaN posterior at one anchor must not null the correlation for
    # an otherwise well-sampled, correlated pair.
    rows = []
    for i in range(6):
        post = 0.9 if i % 2 == 0 else 0.1
        rows.append(_read("siteA", f"r{i}", i % 2 == 0, post))
        rows.append(_read("siteB", f"r{i}", i % 2 == 0, post))
    # one extra spanning read with a NaN posterior at siteA
    rows.append(_read("siteA", "rn", True, float("nan")))
    rows.append(_read("siteB", "rn", True, 0.9))
    pairs = pd.DataFrame({"pair_id": ["p"], "site_a": ["siteA"], "site_b": ["siteB"]})
    pair_summary, _ = joint_occupancy.joint_occupancy(
        called_reads=pd.DataFrame(rows), pairs=pairs, min_spanning_reads=3
    )
    corr = pair_summary.iloc[0]["posterior_corr"]
    assert not pd.isna(corr)
    assert corr > 0.5


def test_summarize_joint_states_fractions():
    _, per_read = joint_occupancy.joint_occupancy(
        called_reads=_called_reads(), pairs=_pairs()
    )
    summary = joint_occupancy.summarize_joint_states(per_read).set_index("pair_id")
    assert summary.loc["coupled", "both"] == 12
    assert summary.loc["coupled", "total"] == 30
    assert summary.loc["coupled", "frac_both"] == pytest.approx(12 / 30)


def test_joint_occupancy_min_spanning_reads_marks_untestable():
    called = _called_reads()
    # require more spanning reads than the independent pair has -> not tested
    pair_summary, _ = joint_occupancy.joint_occupancy(
        called_reads=called, pairs=_pairs(), min_spanning_reads=20
    )
    pair_summary = pair_summary.set_index("pair_id")
    assert pd.isna(pair_summary.loc["independent", "p_value"])
    assert pd.isna(pair_summary.loc["independent", "q_value"])
    assert bool(pair_summary.loc["independent", "significant"]) is False
    # the coupled pair (30 reads) is still tested
    assert not pd.isna(pair_summary.loc["coupled", "p_value"])


def test_joint_occupancy_validates_inputs():
    with pytest.raises(ValueError, match="called_reads requires columns"):
        joint_occupancy.joint_occupancy(
            called_reads=pd.DataFrame({"read_id": ["r"]}), pairs=_pairs()
        )
    with pytest.raises(ValueError, match="pairs requires columns"):
        joint_occupancy.joint_occupancy(
            called_reads=_called_reads(), pairs=pd.DataFrame({"site_a": ["x"]})
        )


def test_anchor_pairs_from_sequence():
    pairs = joint_occupancy.anchor_pairs_from_sequence(
        [("a", "b"), ("c", "d")], distances_bp=[100, 200]
    )
    assert list(pairs["pair_id"]) == ["a|b", "c|d"]
    assert list(pairs["site_a"]) == ["a", "c"]
    assert list(pairs["distance_bp"]) == [100, 200]


def test_joint_occupancy_plotting():
    import matplotlib

    matplotlib.use("Agg")
    from dimelo import plotting, plotting_matplotlib

    pair_summary, per_read = joint_occupancy.joint_occupancy(
        called_reads=_called_reads(), pairs=_pairs()
    )
    states = joint_occupancy.summarize_joint_states(per_read)

    dist_payload = plotting.prepare_joint_state_distribution_data(joint_states=states)
    fig1, _ = plotting_matplotlib.plot_joint_state_distribution_matplotlib(
        dist_payload, title="states"
    )
    assert fig1 is not None

    oe_payload = plotting.prepare_joint_occupancy_distance_data(pair_summary=pair_summary)
    assert len(oe_payload["distance_table"]) == 3  # all pairs have distances
    fig2, _ = plotting_matplotlib.plot_joint_occupancy_distance_matplotlib(
        oe_payload, title="distance"
    )
    assert fig2 is not None
