from __future__ import annotations

import pandas as pd
import pytest

from dimelo import plotting, plotting_matplotlib, region_contrasts
from dimelo.models import ContrastSpec


def _row(
    region_id,
    sample_id,
    condition,
    modified,
    valid,
    *,
    chrom="chr1",
    start=0,
    end=10,
    strand="+",
    replicate=1,
):
    return {
        "region_id": region_id,
        "chromosome": chrom,
        "start": start,
        "end": end,
        "strand": strand,
        "sample_id": sample_id,
        "condition": condition,
        "replicate": replicate,
        "modified_count": modified,
        "valid_count": valid,
        "mod_fraction": 0.0 if valid == 0 else modified / valid,
    }


def _factors(rows):
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# S1a — global normalization helper                                           #
# --------------------------------------------------------------------------- #


def test_apply_global_normalization_shifts_and_rederives_counts():
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 800, 1000),  # mod_fraction 0.8
            _row("reg1", "d1", "control", 200, 1000),  # 0.2
        ]
    )
    factors = _factors(
        [
            {"sample_id": "n1", "motif": "A,0", "global_offset": 0.1},
            {"sample_id": "d1", "motif": "A,0", "global_offset": 0.0},
        ]
    )

    out = region_contrasts._apply_global_normalization(
        evidence, factors=factors, motif="A,0"
    )

    assert list(out["mod_fraction"]) == pytest.approx([0.7, 0.2])
    assert list(out["modified_count"]) == [700, 200]
    # valid_count is never touched
    assert list(out["valid_count"]) == [1000, 1000]


def test_apply_global_normalization_nan_or_absent_offset_is_identity():
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 800, 1000),
            _row("reg1", "d1", "control", 200, 1000),
        ]
    )
    # n1 offset is NaN (zero pooled coverage motif); d1 absent from factors entirely.
    factors = _factors(
        [{"sample_id": "n1", "motif": "A,0", "global_offset": float("nan")}]
    )

    out = region_contrasts._apply_global_normalization(
        evidence, factors=factors, motif="A,0"
    )

    assert list(out["mod_fraction"]) == pytest.approx([0.8, 0.2])
    assert list(out["modified_count"]) == [800, 200]


def test_apply_global_normalization_clamps_to_unit_interval():
    evidence = pd.DataFrame([_row("reg1", "n1", "treated", 100, 1000)])  # 0.1
    factors = _factors(
        [{"sample_id": "n1", "motif": "A,0", "global_offset": 0.5}]
    )  # 0.1 - 0.5 = -0.4 -> clamp 0
    out = region_contrasts._apply_global_normalization(
        evidence, factors=factors, motif="A,0"
    )
    assert out.loc[0, "mod_fraction"] == 0.0
    assert out.loc[0, "modified_count"] == 0


# --------------------------------------------------------------------------- #
# S1b — background pooling + correction helpers                               #
# --------------------------------------------------------------------------- #


def _bg_evidence():
    return pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 900, 1000),  # p 0.9
            _row("reg1", "d1", "control", 300, 1000),  # p 0.3
            _row("reg1", "b1", "bg", 100, 1000),
            _row("reg1", "b2", "bg", 300, 1000),  # pooled bg = 400/2000 = 0.2
        ]
    )


def test_pool_background_rate_coverage_weighted_and_present():
    # Unequal coverage so a coverage-weighted pool (400/1500 ~= 0.267) is
    # distinguishable from a naive mean-of-fractions ((0.1 + 0.6)/2 = 0.35).
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 900, 1000),
            _row("reg1", "d1", "control", 300, 1000),
            _row("reg1", "b1", "bg", 100, 1000),  # fraction 0.1, coverage 1000
            _row("reg1", "b2", "bg", 300, 500),  # fraction 0.6, coverage 500
        ]
    )
    br = region_contrasts._pool_background_rate(evidence, background_conditions=["bg"])
    assert len(br) == 1
    assert br.loc[0, "background_fraction"] == pytest.approx(400 / 1500)
    # not the unweighted mean of the two replicate fractions
    assert br.loc[0, "background_fraction"] != pytest.approx(0.35)
    assert bool(br.loc[0, "background_present"]) is True


def test_pool_background_rate_missing_condition_raises():
    with pytest.raises(ValueError, match="Missing background"):
        region_contrasts._pool_background_rate(
            _bg_evidence(), background_conditions=["not_a_condition"]
        )


def test_apply_background_correction_math_and_invariants():
    contrast = ContrastSpec(
        mode="background_adjusted",
        numerator=["treated"],
        denominator=["control"],
        background=["bg"],
    )
    evidence = _bg_evidence()
    br = region_contrasts._pool_background_rate(evidence, background_conditions=["bg"])
    out = region_contrasts._apply_background_correction(
        evidence, background_rate=br, contrast=contrast
    )

    b = 0.2
    treated = out.loc[out["condition"] == "treated"].iloc[0]
    control = out.loc[out["condition"] == "control"].iloc[0]
    assert treated["mod_fraction"] == pytest.approx((0.9 - b) / (1 - b))  # 0.875
    assert control["mod_fraction"] == pytest.approx((0.3 - b) / (1 - b))  # 0.125
    assert treated["modified_count"] == int(round(0.875 * 1000))
    # k <= n preserved, valid_count unchanged
    assert (out["modified_count"] <= out["valid_count"]).all()
    assert list(out["valid_count"]) == [1000, 1000, 1000, 1000]
    # background rows themselves are untouched
    bg_rows = out.loc[out["condition"] == "bg"]
    assert list(bg_rows["modified_count"]) == [100, 300]


def test_apply_background_correction_saturated_background_gives_zero():
    contrast = ContrastSpec(
        mode="background_adjusted",
        numerator=["treated"],
        denominator=["control"],
        background=["bg"],
    )
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 900, 1000),
            _row("reg1", "d1", "control", 300, 1000),
            _row("reg1", "b1", "bg", 1000, 1000),  # b = 1.0 saturates
        ]
    )
    br = region_contrasts._pool_background_rate(evidence, background_conditions=["bg"])
    out = region_contrasts._apply_background_correction(
        evidence, background_rate=br, contrast=contrast
    )
    corrected = out.loc[out["condition"].isin(["treated", "control"])]
    assert (corrected["mod_fraction"] == 0.0).all()
    assert (corrected["modified_count"] == 0).all()


def test_apply_background_correction_missing_region_background_is_identity():
    contrast = ContrastSpec(
        mode="background_adjusted",
        numerator=["treated"],
        denominator=["control"],
        background=["bg"],
    )
    # reg2 has no background rows; correction should be identity there.
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 900, 1000),
            _row("reg1", "b1", "bg", 200, 1000),
            _row("reg2", "n1", "treated", 700, 1000),
        ]
    )
    br = region_contrasts._pool_background_rate(evidence, background_conditions=["bg"])
    out = region_contrasts._apply_background_correction(
        evidence, background_rate=br, contrast=contrast
    )
    reg2 = out.loc[out["region_id"] == "reg2"].iloc[0]
    assert reg2["mod_fraction"] == pytest.approx(0.7)
    assert reg2["modified_count"] == 700


# --------------------------------------------------------------------------- #
# score_regions integration                                                   #
# --------------------------------------------------------------------------- #


def _patch_evidence(monkeypatch, evidence):
    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )


def test_score_regions_global_normalization_end_to_end(monkeypatch):
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 800, 1000),
            _row("reg1", "d1", "control", 200, 1000),
        ]
    )
    _patch_evidence(monkeypatch, evidence)
    factors = _factors(
        [
            {"sample_id": "n1", "motif": "A,0", "global_offset": 0.1},
            {"sample_id": "d1", "motif": "A,0", "global_offset": 0.0},
        ]
    )
    contrast = ContrastSpec(
        mode="pairwise", numerator=["treated"], denominator=["control"]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
        normalization_mode="global",
        global_normalization_factors=factors,
    )

    row = result.summary.iloc[0]
    assert row["fraction"] == pytest.approx(0.7)
    assert row["reference_fraction"] == pytest.approx(0.2)
    assert row["delta_fraction"] == pytest.approx(0.5)
    assert result.metadata["normalization_mode"] == "global"


def test_score_regions_global_requires_factors(monkeypatch):
    evidence = pd.DataFrame([_row("reg1", "n1", "treated", 800, 1000)])
    _patch_evidence(monkeypatch, evidence)
    contrast = ContrastSpec(
        mode="pairwise", numerator=["treated"], denominator=["control"]
    )
    with pytest.raises(ValueError, match="requires global_normalization_factors"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
            normalization_mode="global",
            global_normalization_factors=None,
        )


def test_score_regions_normalization_mode_metadata_fallback(monkeypatch):
    # param default 'none' but contrast.metadata sets 'global' -> global is applied
    # (proven by the requires-factors ValueError firing).
    evidence = pd.DataFrame([_row("reg1", "n1", "treated", 800, 1000)])
    _patch_evidence(monkeypatch, evidence)
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
        metadata={"normalization_mode": "global"},
    )
    with pytest.raises(ValueError, match="requires global_normalization_factors"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
        )


def test_score_regions_rejects_unsupported_normalization_mode(monkeypatch):
    evidence = pd.DataFrame([_row("reg1", "n1", "treated", 800, 1000)])
    _patch_evidence(monkeypatch, evidence)
    contrast = ContrastSpec(
        mode="pairwise", numerator=["treated"], denominator=["control"]
    )
    with pytest.raises(ValueError, match="Unsupported normalization_mode"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
            normalization_mode="bogus",
        )


def test_score_regions_rejects_normalization_for_non_ensemble_unit():
    # global normalization is only wired into the ensemble path; requesting it for a
    # single_read/cluster_occupancy unit must fail loud rather than silently mislabel.
    contrast = ContrastSpec(
        mode="pairwise", numerator=["treated"], denominator=["control"]
    )
    with pytest.raises(ValueError, match="only supported for analysis_unit"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="extract_reads",
            read_table=pd.DataFrame(),
            normalization_mode="global",
            global_normalization_factors=pd.DataFrame(),
        )


def test_score_regions_background_adjusted_end_to_end(monkeypatch):
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 900, 1000),  # p 0.9
            _row("reg1", "d1", "control", 300, 1000),  # p 0.3
            _row("reg1", "b1", "bg", 200, 1000),  # b 0.2
        ]
    )
    _patch_evidence(monkeypatch, evidence)
    contrast = ContrastSpec(
        mode="background_adjusted",
        numerator=["treated"],
        denominator=["control"],
        background=["bg"],
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="beta_binomial",
    )

    row = result.summary.iloc[0]
    # corrected fractions: (p - b)/(1 - b)
    assert row["fraction"] == pytest.approx(0.875)
    assert row["reference_fraction"] == pytest.approx(0.125)
    assert row["background_fraction"] == pytest.approx(0.2)
    assert bool(row["background_present"]) is True
    assert {"background_fraction", "background_present"}.issubset(result.summary.columns)
    assert {"p_value", "adjusted_p_value"}.issubset(result.summary.columns)

    # raw-vs-corrected overlay carries the pre-correction pooled fractions
    assert "region_effect_sizes_raw" in result.plot_data
    raw = result.plot_data["region_effect_sizes_raw"].iloc[0]
    assert raw["raw_fraction"] == pytest.approx(0.9)
    assert raw["raw_reference_fraction"] == pytest.approx(0.3)

    assert result.metadata["background_correction"] is True
    assert result.metadata["background_conditions"] == ["bg"]


def test_score_regions_global_then_background_order(monkeypatch):
    # Global normalization runs BEFORE background pooling, so the background rate is
    # computed on normalized counts. With a +0.2 offset on the bg sample, its
    # normalized fraction (and thus b) shifts down, changing the correction.
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 900, 1000),
            _row("reg1", "d1", "control", 300, 1000),
            _row("reg1", "b1", "bg", 400, 1000),  # raw b 0.4
        ]
    )
    _patch_evidence(monkeypatch, evidence)
    factors = _factors(
        [
            {"sample_id": "n1", "motif": "A,0", "global_offset": 0.0},
            {"sample_id": "d1", "motif": "A,0", "global_offset": 0.0},
            {"sample_id": "b1", "motif": "A,0", "global_offset": 0.2},
        ]
    )
    contrast = ContrastSpec(
        mode="background_adjusted",
        numerator=["treated"],
        denominator=["control"],
        background=["bg"],
    )
    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
        normalization_mode="global",
        global_normalization_factors=factors,
    )
    # normalized bg fraction = 0.4 - 0.2 = 0.2 -> b = 0.2 (not 0.4)
    row = result.summary.iloc[0]
    assert row["background_fraction"] == pytest.approx(0.2)
    assert row["fraction"] == pytest.approx((0.9 - 0.2) / (1 - 0.2))  # 0.875


def test_existing_pairwise_path_unchanged_without_normalization(monkeypatch):
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 800, 1000),
            _row("reg1", "d1", "control", 200, 1000),
        ]
    )
    _patch_evidence(monkeypatch, evidence)
    contrast = ContrastSpec(
        mode="pairwise", numerator=["treated"], denominator=["control"]
    )
    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
    )
    row = result.summary.iloc[0]
    assert row["fraction"] == pytest.approx(0.8)
    assert row["reference_fraction"] == pytest.approx(0.2)
    assert result.metadata["normalization_mode"] == "none"
    assert set(result.plot_data) == {"region_effect_sizes"}
    assert "background_fraction" not in result.summary.columns


# --------------------------------------------------------------------------- #
# S1 plotting overlay (raw vs corrected)                                      #
# --------------------------------------------------------------------------- #


def _background_adjusted_result(monkeypatch):
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 900, 1000),
            _row("reg1", "d1", "control", 300, 1000),
            _row("reg1", "b1", "bg", 200, 1000),
            _row("reg2", "n1", "treated", 500, 1000, start=100, end=110),
            _row("reg2", "d1", "control", 500, 1000, start=100, end=110),
            _row("reg2", "b1", "bg", 100, 1000, start=100, end=110),
        ]
    )
    _patch_evidence(monkeypatch, evidence)
    contrast = ContrastSpec(
        mode="background_adjusted",
        numerator=["treated"],
        denominator=["control"],
        background=["bg"],
    )
    return region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
    )


def test_prepare_correction_overlay_emits_raw_and_corrected(monkeypatch):
    result = _background_adjusted_result(monkeypatch)
    payload = plotting.prepare_region_contrast_correction_overlay_data(result=result)
    overlay = payload["overlay_table"]

    assert set(overlay["stage"]) == {"raw", "corrected"}
    assert "background_fraction" in overlay.columns

    raw = overlay[overlay["stage"] == "raw"].set_index("region_id")
    corrected = overlay[overlay["stage"] == "corrected"].set_index("region_id")
    # raw delta for reg1 = 0.9 - 0.3 = 0.6
    assert raw.loc["reg1", "delta_fraction"] == pytest.approx(0.6)
    # corrected delta for reg1 = 0.875 - 0.125 = 0.75 (b = 0.2)
    assert corrected.loc["reg1", "delta_fraction"] == pytest.approx(0.75)
    assert corrected.loc["reg1", "background_fraction"] == pytest.approx(0.2)


def test_prepare_correction_overlay_requires_background_result(monkeypatch):
    evidence = pd.DataFrame(
        [
            _row("reg1", "n1", "treated", 800, 1000),
            _row("reg1", "d1", "control", 200, 1000),
        ]
    )
    _patch_evidence(monkeypatch, evidence)
    contrast = ContrastSpec(
        mode="pairwise", numerator=["treated"], denominator=["control"]
    )
    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
    )
    with pytest.raises(ValueError, match="background_adjusted"):
        plotting.prepare_region_contrast_correction_overlay_data(result=result)


def test_prepare_correction_overlay_top_n(monkeypatch):
    result = _background_adjusted_result(monkeypatch)
    payload = plotting.prepare_region_contrast_correction_overlay_data(
        result=result, top_n=1
    )
    overlay = payload["overlay_table"]
    # only the single largest-|corrected-delta| region survives, and it must be reg1
    # (corrected delta 0.75) not reg2 (corrected delta 0.0) -- catches a sort-order swap
    assert overlay["region_id"].nunique() == 1
    assert overlay["region_id"].unique().tolist() == ["reg1"]
    assert payload["metadata"]["region_order"] == ["reg1"]


def test_correction_overlay_renderer_smoke(monkeypatch):
    import matplotlib

    matplotlib.use("Agg")
    result = _background_adjusted_result(monkeypatch)
    payload = plotting.prepare_region_contrast_correction_overlay_data(result=result)
    fig, ax = plotting_matplotlib.plot_region_contrast_correction_overlay_matplotlib(
        payload, title="overlay"
    )
    assert fig is not None
