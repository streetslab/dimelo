from __future__ import annotations

from pathlib import Path

import pandas as pd

from dimelo import dmr, run_modkit


def _fake_capabilities() -> run_modkit.ModkitCapabilities:
    return run_modkit.ModkitCapabilities(
        executable="modkit",
        version_raw="modkit 0.6.1",
        version="0.6.1",
        version_tuple=(0, 6, 1),
        supports_mod_threshold=True,
        supports_mod_thresholds=False,
        supports_modified_bases=True,
        supports_force_allow_implicit=False,
        supports_extract_subcommands=True,
        extract_supports_reference_long=True,
        extract_supports_reference_short=True,
    )


def test_run_dmr_pair_builds_command_and_parses_outputs(tmp_path, monkeypatch):
    a = tmp_path / "a.bed.gz"
    b = tmp_path / "b.bed.gz"
    ref = tmp_path / "ref.fa"
    out = tmp_path / "pair.bed"
    seg = tmp_path / "pair.segment.bed"
    for path in (a, b, ref):
        path.write_text("x")

    captured: dict[str, object] = {}

    def fake_run(cmd, check):
        captured["cmd"] = list(cmd)
        out.write_text(
            "#chrom\tstart\tend\tscore\tmap_pvalue\teffect_size\ta_total\n"
            "chr1\t1\t2\t0.3\t0.001\t0.4\t10\n"
            "chr1\t3\t4\t0.01\t0.6\t0.01\t3\n"
        )
        seg.write_text("#chrom\tchrom_start\tchrom_end\tname\tscore\nchr1\t1\t10\tdiff\t0.2\n")
        return None

    monkeypatch.setattr(dmr.subprocess, "run", fake_run)
    monkeypatch.setattr(run_modkit, "_ensure_modkit_available", lambda **_: _fake_capabilities())

    result = dmr.run_dmr_pair(
        control_bed_methyl=a,
        experiment_bed_methyl=b,
        ref_genome=ref,
        out_path=out,
        segment_path=seg,
        bases=["A"],
        threads=2,
        io_threads=1,
        min_total_coverage=5,
    )

    assert isinstance(result.sites, pd.DataFrame)
    assert isinstance(result.segments, pd.DataFrame)
    assert result.high_confidence_sites.shape[0] == 1
    command = captured["cmd"]
    assert "--segment" in command
    assert "--base" in command
    assert "--header" in command


def test_run_dmr_pair_rejects_regions_and_segment_together(tmp_path):
    a = tmp_path / "a.bed.gz"
    b = tmp_path / "b.bed.gz"
    ref = tmp_path / "ref.fa"
    regions = tmp_path / "regions.bed"
    for path in (a, b, ref, regions):
        path.write_text("x")

    try:
        dmr.run_dmr_pair(
            control_bed_methyl=a,
            experiment_bed_methyl=b,
            ref_genome=ref,
            out_path=tmp_path / "out.bed",
            regions_bed=regions,
            segment_path=tmp_path / "seg.bed",
        )
    except ValueError as exc:
        assert "does not allow --regions-bed with --segment" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("Expected ValueError for regions_bed + segment_path")


def test_run_dmr_multi_collects_pair_files(tmp_path, monkeypatch):
    sample_a = tmp_path / "a.bed.gz"
    sample_b = tmp_path / "b.bed.gz"
    regions = tmp_path / "regions.bed"
    ref = tmp_path / "ref.fa"
    out_dir = tmp_path / "dmr_multi"
    for path in (sample_a, sample_b, regions, ref):
        path.write_text("x")

    captured: dict[str, object] = {}

    def fake_run(cmd, check):
        captured["cmd"] = list(cmd)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "s1_vs_s2.bed").write_text("chr1\t1\t2\n")
        return None

    monkeypatch.setattr(dmr.subprocess, "run", fake_run)
    monkeypatch.setattr(run_modkit, "_ensure_modkit_available", lambda **_: _fake_capabilities())

    result = dmr.run_dmr_multi(
        samples={"s1": sample_a, "s2": sample_b},
        regions_bed=regions,
        ref_genome=ref,
        out_dir=out_dir,
        bases=["A"],
        threads=2,
    )

    assert result.pair_files.shape[0] == 1
    assert result.pair_files.iloc[0]["pair_name"] == "s1_vs_s2"
    command = captured["cmd"]
    assert "--sample" in command
    assert "--regions-bed" in command
