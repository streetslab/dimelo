from __future__ import annotations

import pandas as pd
import pytest

from dimelo import dmr, plotting, plotting_matplotlib
from dimelo.models import ModkitDMRMultiResult, ModkitDMRPairResult


def _pair_result(sites: pd.DataFrame, segments: pd.DataFrame | None = None):
    # Build high_confidence_sites the same way run_dmr_pair does, so the fixture
    # stays faithful and tolerates frames that omit map_pvalue/effect_size.
    high_confidence = dmr._select_high_confidence_sites(
        sites, pvalue_max=0.01, abs_effect_size_min=0.1
    )
    return ModkitDMRPairResult(
        output_path="dmr_sites.bed",
        segment_path=None if segments is None else "dmr_segments.bed",
        command=["modkit", "dmr", "pair"],
        sites=sites,
        segments=segments,
        high_confidence_sites=high_confidence,
        metadata={"bases": ["A"]},
    )


def _sites_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "start": [10, 40, 5],
            "end": [12, 42, 7],
            "score": [0.3, 0.01, 0.9],
            "map_pvalue": [0.001, 0.6, 0.002],
            "effect_size": [0.4, 0.01, -0.5],
            "a_total": [10, 3, 8],
        }
    )


def test_prepare_dmr_site_data_computes_features_and_high_confidence_flag():
    result = _pair_result(_sites_frame())

    payload = plotting.prepare_dmr_site_data(result=result)
    table = payload["site_table"]

    assert list(table["contig"]) == ["chr1", "chr1", "chr2"]
    # midpoint / span
    row = table.loc[table["start"] == 10].iloc[0]
    assert row["midpoint"] == 11.0
    assert row["span"] == 2.0
    # -log10(0.001) == 3
    assert row["neg_log10_pvalue"] == pytest.approx(3.0)
    # two sites pass p<=0.01 & |effect|>=0.1 (chr1:10 and chr2:5)
    assert table["is_high_confidence"].sum() == 2
    assert bool(table.loc[table["start"] == 40].iloc[0]["is_high_confidence"]) is False
    assert payload["metadata"]["n_high_confidence"] == 2


def test_prepare_dmr_site_data_top_n_by_abs_effect_size():
    result = _pair_result(_sites_frame())

    payload = plotting.prepare_dmr_site_data(result=result, top_n=2)
    table = payload["site_table"]

    assert len(table) == 2
    # largest |effect| are 0.5 (chr2) and 0.4 (chr1:10); 0.01 dropped
    assert set(zip(table["contig"], table["start"], strict=True)) == {
        ("chr1", 10),
        ("chr2", 5),
    }


def test_prepare_dmr_site_data_contigs_filter_and_empty():
    result = _pair_result(_sites_frame())

    payload = plotting.prepare_dmr_site_data(result=result, contigs=["chr2"])
    assert list(payload["site_table"]["contig"]) == ["chr2"]

    with pytest.raises(ValueError, match="not present"):
        plotting.prepare_dmr_site_data(result=result, contigs=["chrX"])

    empty = _pair_result(_sites_frame().iloc[0:0])
    empty_payload = plotting.prepare_dmr_site_data(result=empty)
    assert empty_payload["site_table"].empty
    assert empty_payload["metadata"]["contig_order"] == []


def test_prepare_dmr_site_data_zero_pvalue_is_nan_not_inf():
    sites = _sites_frame()
    sites.loc[0, "map_pvalue"] = 0.0
    result = _pair_result(sites)

    table = plotting.prepare_dmr_site_data(result=result)["site_table"]
    row = table.loc[table["start"] == 10].iloc[0]
    assert pd.isna(row["neg_log10_pvalue"])
    assert not any(table["neg_log10_pvalue"].dropna().apply(lambda v: v == float("inf")))


def test_prepare_dmr_site_data_missing_pvalue_and_effect_columns():
    # regions-mode output: no map_pvalue, no effect_size.
    sites = pd.DataFrame(
        {"chrom": ["chr1", "chr2"], "start": [10, 20], "end": [12, 22]}
    )
    result = _pair_result(sites)

    payload = plotting.prepare_dmr_site_data(result=result)
    table = payload["site_table"]
    assert table["neg_log10_pvalue"].isna().all()
    # still a valid, populated table with derived coordinate features
    assert list(table["contig"]) == ["chr1", "chr2"]
    assert list(table["span"]) == [2.0, 2.0]


def test_prepare_dmr_site_data_only_high_confidence():
    result = _pair_result(_sites_frame())

    payload = plotting.prepare_dmr_site_data(result=result, only_high_confidence=True)
    table = payload["site_table"]
    # sourced from high_confidence_sites (2 rows) and all flagged True
    assert len(table) == len(result.high_confidence_sites)
    assert bool(table["is_high_confidence"].all()) is True
    assert payload["metadata"]["only_high_confidence"] is True


def test_prepare_dmr_segment_data_none_and_populated():
    no_seg = _pair_result(_sites_frame(), segments=None)
    payload = plotting.prepare_dmr_segment_data(result=no_seg)
    assert payload["segment_table"].empty
    assert payload["metadata"]["has_segments"] is False

    segments = pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "chrom_start": [1, 100],
            "chrom_end": [10, 140],
            "name": ["diff", "diff"],
            "score": [0.2, 0.7],
        }
    )
    with_seg = _pair_result(_sites_frame(), segments=segments)
    seg_payload = plotting.prepare_dmr_segment_data(result=with_seg)
    seg_table = seg_payload["segment_table"]
    assert seg_payload["metadata"]["has_segments"] is True
    assert list(seg_table["span"]) == [9.0, 40.0]
    assert list(seg_table["midpoint"]) == [5.5, 120.0]


def _write_multi_pair_files(tmp_path):
    # modkit regions / multi mode does NOT emit map_pvalue, so significance is
    # driven by the effect-size cut alone. Row 1 clears |effect|>=0.1; row 2 does not.
    good = tmp_path / "s1_vs_s2.bed"
    good.write_text(
        "#chrom\tstart\tend\tscore\teffect_size\ta_total\n"
        "chr1\t1\t2\t0.3\t0.40\t10\n"
        "chr1\t3\t4\t0.01\t0.05\t3\n"
    )
    empty = tmp_path / "s1_vs_s3.bed"
    empty.write_text("#chrom\tstart\tend\tscore\teffect_size\ta_total\n")
    pair_files = pd.DataFrame(
        {
            "pair_file": [good, empty],
            "pair_name": ["s1_vs_s2", "s1_vs_s3"],
        }
    )
    return ModkitDMRMultiResult(
        out_dir=tmp_path,
        command=["modkit", "dmr", "multi"],
        pair_files=pair_files,
        metadata={"bases": ["A"]},
    )


def test_prepare_dmr_multi_summary_reads_pair_files(tmp_path):
    result = _write_multi_pair_files(tmp_path)

    payload = plotting.prepare_dmr_multi_summary_data(result=result)
    summary = payload["summary_table"]

    s1_s2 = summary.loc[summary["pair_name"] == "s1_vs_s2"].iloc[0]
    assert s1_s2["n_sites"] == 2
    assert s1_s2["n_significant"] == 1
    assert s1_s2["significant_fraction"] == pytest.approx(0.5)
    s1_s3 = summary.loc[summary["pair_name"] == "s1_vs_s3"].iloc[0]
    assert s1_s3["n_sites"] == 0
    assert s1_s3["n_significant"] == 0


def test_prepare_dmr_multi_summary_significance_tracks_effect_threshold(tmp_path):
    # The significance count must respond to the effect-size cut specifically:
    # raising the threshold above both rows drops the count to 0; lowering it
    # below both rows raises it to 2. A swapped/dropped filter would not do this.
    result = _write_multi_pair_files(tmp_path)

    strict = plotting.prepare_dmr_multi_summary_data(
        result=result, abs_effect_size_min=0.5
    )["summary_table"]
    assert (
        strict.loc[strict["pair_name"] == "s1_vs_s2"].iloc[0]["n_significant"] == 0
    )

    loose = plotting.prepare_dmr_multi_summary_data(
        result=result, abs_effect_size_min=0.01
    )["summary_table"]
    assert (
        loose.loc[loose["pair_name"] == "s1_vs_s2"].iloc[0]["n_significant"] == 2
    )


def test_dmr_renderers_smoke(tmp_path):
    import matplotlib

    matplotlib.use("Agg")

    result = _pair_result(
        _sites_frame(),
        segments=pd.DataFrame(
            {
                "chrom": ["chr1"],
                "chrom_start": [1],
                "chrom_end": [10],
                "score": [0.2],
            }
        ),
    )

    site_payload = plotting.prepare_dmr_site_data(result=result)
    fig, axes = plotting_matplotlib.plot_dmr_site_matplotlib(
        site_payload, title="sites"
    )
    assert fig is not None

    seg_payload = plotting.prepare_dmr_segment_data(result=result)
    fig2, axes2 = plotting_matplotlib.plot_dmr_segment_matplotlib(seg_payload)
    assert fig2 is not None

    good = tmp_path / "s1_vs_s2.bed"
    good.write_text(
        "#chrom\tstart\tend\tscore\tmap_pvalue\teffect_size\ta_total\n"
        "chr1\t1\t2\t0.3\t0.001\t0.4\t10\n"
    )
    multi = ModkitDMRMultiResult(
        out_dir=tmp_path,
        command=["modkit", "dmr", "multi"],
        pair_files=pd.DataFrame(
            {"pair_file": [good], "pair_name": ["s1_vs_s2"]}
        ),
        metadata={},
    )
    multi_payload = plotting.prepare_dmr_multi_summary_data(result=multi)
    fig3, ax3 = plotting_matplotlib.plot_dmr_multi_summary_matplotlib(multi_payload)
    assert fig3 is not None
