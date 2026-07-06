from __future__ import annotations

import pandas as pd
import pytest

from dimelo import background, plotting, plotting_matplotlib


def _read(region_id, condition, read_id, k, n=20):
    return {
        "region_id": region_id,
        "condition": condition,
        "read_id": read_id,
        "modified_count": k,
        "valid_count": n,
    }


def _evidence():
    rows = []
    # bound_site: control low; target mostly bound (16 high, 4 background)
    for i, k in enumerate([0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0]):
        rows.append(_read("bound_site", "control", f"cB{i}", k))
    for i in range(16):
        rows.append(_read("bound_site", "target", f"tB{i}", 18))
    for i in range(4):
        rows.append(_read("bound_site", "target", f"tBbg{i}", 1))
    # unbound_site: control low; target also background-like (no binding)
    for i, k in enumerate([0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0, 1, 0, 1, 0, 2, 1, 0]):
        rows.append(_read("unbound_site", "control", f"cU{i}", k))
    for i in range(20):
        rows.append(_read("unbound_site", "target", f"tU{i}", 1))
    return pd.DataFrame(rows)


def test_binding_strength_control_calibrated_ranking():
    strength = background.binding_strength(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        fdr=0.05,
        min_control_reads=5,
    ).set_index("region_id")

    # bound_site: 16/20 bound. occupancy is the HARD point estimate (0.8), distinct from
    # the soft Beta-posterior mean; control occupancy ~0; positive excess; significant.
    assert strength.loc["bound_site", "occupancy"] == pytest.approx(16 / 20)
    assert strength.loc["bound_site", "occupancy"] != pytest.approx(
        strength.loc["bound_site", "occupancy_posterior_mean"]
    )
    assert strength.loc["bound_site", "control_occupancy"] < 0.2
    assert strength.loc["bound_site", "occupancy_excess"] > 0.3
    assert bool(strength.loc["bound_site", "significant"]) is True
    # the credible interval separates signal from noise (not a Beta self-consistency check)
    assert strength.loc["bound_site", "occupancy_ci_low"] > 0.5
    # unbound_site: hard occupancy is exactly 0, distinct from the nonzero posterior mean
    assert strength.loc["unbound_site", "occupancy"] == pytest.approx(0.0)
    assert strength.loc["unbound_site", "occupancy_posterior_mean"] > 0.0
    assert bool(strength.loc["unbound_site", "significant"]) is False
    assert strength.loc["bound_site", "rank"] < strength.loc["unbound_site", "rank"]


def test_binding_strength_control_calibration_rejects_leaky_site():
    # A site where the CONTROL itself carries spurious high signal (few control reads ->
    # pooled low null -> control reads called true): target occupancy is high, but so is
    # control occupancy, so the control-calibrated excess is ~0 and the site is NOT
    # significant. This exercises the control-subtraction (a mutant that dropped the
    # control arm would report the raw occupancy and call it significant).
    rows = []
    for i in range(120):  # large tight-low control bulk -> low pooled null
        rows.append(_read("bulk_site", "control", f"cb{i}", 1 if i % 5 == 0 else 0))
    for i in range(20):
        rows.append(_read("bulk_site", "target", f"tb{i}", 1))
    for i in range(2):  # < min_control_reads -> uses the pooled (low) null
        rows.append(_read("leaky_site", "control", f"cL{i}", 20))
    for i in range(12):
        rows.append(_read("leaky_site", "target", f"tL{i}", 20))

    strength = background.binding_strength(
        evidence=pd.DataFrame(rows),
        target_conditions=["target"],
        control_conditions=["control"],
        fdr=0.05,
        min_control_reads=5,
    ).set_index("region_id")

    leaky = strength.loc["leaky_site"]
    assert leaky["occupancy"] > 0.8  # raw occupancy looks strong...
    assert leaky["control_occupancy"] > 0.5  # ...but the control is just as high
    assert leaky["occupancy_excess"] == pytest.approx(
        leaky["occupancy"] - leaky["control_occupancy"]
    )
    assert leaky["occupancy_excess"] < 0.2
    # control-calibration rejects it despite the high raw occupancy
    assert bool(leaky["significant"]) is False


def test_binding_strength_empty_target_returns_typed_frame():
    strength = background.binding_strength(
        evidence=_evidence(),
        target_conditions=["not_a_condition"],
        control_conditions=["control"],
    )
    assert strength.empty
    assert "occupancy_excess" in strength.columns


def test_binding_strength_consistency_qc_true_reads_carry_more_signal():
    called = background.call_true_signal_reads(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    )
    qc = background.binding_strength_consistency_qc(called).set_index("region_id")
    # at the bound site, true-signal reads carry much higher methylation than background
    assert qc.loc["bound_site", "mean_signal_true"] > qc.loc[
        "bound_site", "mean_signal_background"
    ]
    assert qc.loc["bound_site", "signal_gap"] > 0.3


def test_binding_strength_consistency_qc_requires_columns():
    with pytest.raises(ValueError, match="requires columns"):
        background.binding_strength_consistency_qc(
            pd.DataFrame({"region_id": ["s"], "is_true_signal": [True]})
        )


def test_binding_strength_plotting():
    import matplotlib

    matplotlib.use("Agg")
    strength = background.binding_strength(
        evidence=_evidence(),
        target_conditions=["target"],
        control_conditions=["control"],
        min_control_reads=5,
    )
    payload = plotting.prepare_binding_strength_data(binding_strength=strength)
    assert list(payload["binding_table"]["region_id"]) == payload["metadata"][
        "region_order"
    ]
    fig, _ = plotting_matplotlib.plot_binding_strength_matplotlib(payload, title="strength")
    assert fig is not None

    top = plotting.prepare_binding_strength_data(binding_strength=strength, top_n=1)
    assert len(top["binding_table"]) == 1
