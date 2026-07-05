from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dimelo import background, plotting, plotting_matplotlib


def _read(region_id, condition, read_id, k, n):
    return {
        "region_id": region_id,
        "condition": condition,
        "read_id": read_id,
        "modified_count": k,
        "valid_count": n,
    }


def test_fit_beta_binomial_recovers_low_background_rate():
    # control reads clustered near a low methylation rate
    ks = np.array([0, 1, 0, 2, 1, 0, 1, 0])
    ns = np.array([20, 20, 20, 20, 20, 20, 20, 20])
    alpha, beta = background.fit_beta_binomial(ks, ns)
    assert alpha > 0 and beta > 0
    # mean = alpha/(alpha+beta) should track the pooled control rate (~5/160)
    assert alpha / (alpha + beta) == pytest.approx(ks.sum() / ns.sum(), abs=0.05)


def test_fit_beta_binomial_requires_coverage():
    with pytest.raises(ValueError, match="valid_count > 0"):
        background.fit_beta_binomial(np.array([0, 0]), np.array([0, 0]))


def _evidence():
    rows = []
    # siteA: 6 low-methylation control reads -> per-site null (min_control_reads=5)
    for i, k in enumerate([0, 1, 0, 1, 2, 0]):
        rows.append(_read("siteA", "control", f"cA{i}", k, 20))
    # siteB: only 2 control reads -> pooled fallback
    for i, k in enumerate([0, 1]):
        rows.append(_read("siteB", "control", f"cB{i}", k, 20))
    # target reads
    rows.append(_read("siteA", "target", "tA_high", 18, 20))  # clearly modified
    rows.append(_read("siteA", "target", "tA_bg", 1, 20))  # background-like
    rows.append(_read("siteB", "target", "tB_high", 19, 20))  # clearly modified
    return pd.DataFrame(rows)


def test_call_true_signal_reads_separates_signal_from_background():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        fdr=0.05,
        min_control_reads=5,
    )
    called = called.set_index("read_id")

    # only the 3 target reads are called
    assert set(called.index) == {"tA_high", "tA_bg", "tB_high"}
    # high-methylation reads are true signal; background-like is not
    assert bool(called.loc["tA_high", "is_true_signal"]) is True
    assert bool(called.loc["tB_high", "is_true_signal"]) is True
    assert bool(called.loc["tA_bg", "is_true_signal"]) is False
    assert called.loc["tA_high", "p_value"] < called.loc["tA_bg", "p_value"]
    # null source: siteA has enough control coverage (per-site), siteB falls back to pooled
    assert called.loc["tA_high", "null_source"] == "per_site"
    assert called.loc["tB_high", "null_source"] == "pooled"
    # background_rate is the fitted control mean (low)
    assert called.loc["tA_high", "background_rate"] < 0.3


def test_call_true_signal_reads_zero_coverage_is_untestable():
    evidence = _evidence()
    evidence = pd.concat(
        [evidence, pd.DataFrame([_read("siteA", "target", "tA_zero", 0, 0)])],
        ignore_index=True,
    )
    called = background.call_true_signal_reads(
        evidence=evidence,
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    ).set_index("read_id")

    assert np.isnan(called.loc["tA_zero", "p_value"])
    assert np.isnan(called.loc["tA_zero", "q_value"])
    assert bool(called.loc["tA_zero", "is_true_signal"]) is False


def test_call_true_signal_reads_requires_control_and_columns():
    with pytest.raises(ValueError, match="requires columns"):
        background.call_true_signal_reads(
            evidence=pd.DataFrame({"region_id": ["s"], "condition": ["target"]}),
            target_conditions=["target"],
            control_conditions=["control"],
        )
    with pytest.raises(ValueError, match="no control reads"):
        background.call_true_signal_reads(
            evidence=_evidence(),
            target_conditions=["target"],
            control_conditions=["not_a_condition"],
        )


def test_call_true_signal_reads_empty_target_returns_typed_frame():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["not_a_condition"],
        control_conditions=["control"],
    )
    assert called.empty
    assert "is_true_signal" in called.columns


def test_occupancy_posterior_ranks_signal_above_background():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    ).set_index("read_id")

    assert "occupancy_posterior" in called.columns
    # high-methylation reads carry higher occupancy posterior than the background-like one
    assert called.loc["tA_high", "occupancy_posterior"] > called.loc[
        "tA_bg", "occupancy_posterior"
    ]
    assert 0.0 <= called.loc["tA_bg", "occupancy_posterior"] <= 1.0
    assert 0.0 <= called.loc["tA_high", "occupancy_posterior"] <= 1.0


def test_occupancy_posterior_nan_for_zero_coverage():
    evidence = pd.concat(
        [_evidence(), pd.DataFrame([_read("siteA", "target", "tA_zero", 0, 0)])],
        ignore_index=True,
    )
    called = background.call_true_signal_reads(
        evidence=evidence,
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    ).set_index("read_id")
    assert np.isnan(called.loc["tA_zero", "occupancy_posterior"])


def test_background_removed_pileup_hard_uses_true_signal_only():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    )
    pileup = background.background_removed_pileup(called, weighting="hard").set_index(
        "region_id"
    )
    # siteA: only tA_high (18/20) is true signal (tA_bg dropped)
    assert pileup.loc["siteA", "modified_count"] == pytest.approx(18.0)
    assert pileup.loc["siteA", "valid_count"] == pytest.approx(20.0)
    assert pileup.loc["siteA", "mod_fraction"] == pytest.approx(0.9)
    # siteB: tB_high (19/20) is true signal
    assert pileup.loc["siteB", "mod_fraction"] == pytest.approx(0.95)


def test_background_removed_pileup_posterior_weights_counts():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    )
    pileup = background.background_removed_pileup(called, weighting="posterior")
    # posterior-weighted counts are bounded by the raw counts and non-negative
    assert (pileup["valid_count"] >= 0).all()
    assert (pileup["mod_fraction"] >= 0).all() and (pileup["mod_fraction"] <= 1).all()


def test_background_removed_pileup_rejects_bad_weighting():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    )
    with pytest.raises(ValueError, match="weighting must be"):
        background.background_removed_pileup(called, weighting="soft")


def test_summarize_site_occupancy_rate():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    )
    occupancy = background.summarize_site_occupancy(called).set_index("region_id")
    # siteA: 2 target reads, 1 true signal -> 0.5 ; siteB: 1 read, 1 true -> 1.0
    assert occupancy.loc["siteA", "n_reads"] == 2
    assert occupancy.loc["siteA", "n_true_signal"] == 1
    assert occupancy.loc["siteA", "occupancy_rate"] == pytest.approx(0.5)
    assert occupancy.loc["siteB", "occupancy_rate"] == pytest.approx(1.0)


# --------------------------------------------------------------------------- #
# S2 plotting                                                                 #
# --------------------------------------------------------------------------- #


def _called():
    return background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    )


def test_prepare_true_signal_read_data_orders_and_validates():
    payload = plotting.prepare_true_signal_read_data(called_reads=_called())
    table = payload["read_table"]
    assert "read_order" in table.columns
    assert set(payload["metadata"]["region_ids"]) == {"siteA", "siteB"}
    # within a region, reads are ordered by ascending mod fraction
    site_a = table.loc[table["region_id"] == "siteA"]
    assert list(site_a["read_mod_fraction"]) == sorted(site_a["read_mod_fraction"])

    with pytest.raises(ValueError, match="color_by"):
        plotting.prepare_true_signal_read_data(called_reads=_called(), color_by="nope")
    with pytest.raises(ValueError, match="No called reads"):
        plotting.prepare_true_signal_read_data(called_reads=_called(), region_id="ghost")


def test_prepare_background_removed_pileup_overlay():
    payload = plotting.prepare_background_removed_pileup_overlay_data(
        called_reads=_called()
    )
    overlay = payload["overlay_table"]
    assert set(overlay["stage"]) == {"raw", "background_removed"}
    removed = overlay.loc[overlay["stage"] == "background_removed"].set_index("region_id")
    # siteA background-removed = true-signal read only (18/20 = 0.9)
    assert removed.loc["siteA", "mod_fraction"] == pytest.approx(0.9)
    raw = overlay.loc[overlay["stage"] == "raw"].set_index("region_id")
    # siteA raw pools both target reads (18+1)/(20+20) = 0.475
    assert raw.loc["siteA", "mod_fraction"] == pytest.approx(19 / 40)


def test_s2_renderers_smoke():
    import matplotlib

    matplotlib.use("Agg")
    called = _called()
    occupancy = background.summarize_site_occupancy(called)

    fig1, _ = plotting_matplotlib.plot_true_signal_read_raster_matplotlib(
        plotting.prepare_true_signal_read_data(called_reads=called), title="raster"
    )
    assert fig1 is not None
    fig1b, _ = plotting_matplotlib.plot_true_signal_read_raster_matplotlib(
        plotting.prepare_true_signal_read_data(
            called_reads=called, color_by="occupancy_posterior"
        )
    )
    assert fig1b is not None
    fig2, _ = plotting_matplotlib.plot_occupancy_rate_track_matplotlib(
        plotting.prepare_occupancy_rate_track_data(site_occupancy=occupancy)
    )
    assert fig2 is not None
    fig3, _ = plotting_matplotlib.plot_background_removed_pileup_overlay_matplotlib(
        plotting.prepare_background_removed_pileup_overlay_data(called_reads=called)
    )
    assert fig3 is not None
