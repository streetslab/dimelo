"""End-to-end parse -> pileup -> modkit dmr -> prep/plot integration test.

Skipped automatically unless a real ``modkit`` executable AND the demo inputs
(the committed CTCF barcode BAMs plus the CHM13 reference FASTA used to align
them) are all present. In CI where modkit and/or the multi-GB reference are not
available, the whole module is skipped rather than failing.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dimelo import plotting, plotting_matplotlib

DATA_DIR = Path("dimelo/test/data")
REF_GENOME = Path("dimelo/test/output/chm13.draft_v1.0.fasta")
BARCODE17_BAM = DATA_DIR / "barcode17.merged.sorted.ctcf_demo.sorted.bam"
BARCODE18_BAM = DATA_DIR / "barcode18.merged.sorted.ctcf_demo.sorted.bam"
ANALYSIS_REGION = "chr16:63442391-63446452"


def _modkit_available() -> bool:
    try:
        from dimelo import run_modkit

        run_modkit._ensure_modkit_available(quiet=True)
        return True
    except Exception:
        return False


def _missing_inputs() -> list[str]:
    required = [
        BARCODE17_BAM,
        BARCODE17_BAM.with_suffix(BARCODE17_BAM.suffix + ".bai"),
        BARCODE18_BAM,
        BARCODE18_BAM.with_suffix(BARCODE18_BAM.suffix + ".bai"),
        REF_GENOME,
        REF_GENOME.with_suffix(REF_GENOME.suffix + ".fai"),
    ]
    return [str(path) for path in required if not path.exists()]


_SKIP_REASON: str | None = None
if not _modkit_available():
    _SKIP_REASON = "modkit executable not available"
else:
    _missing = _missing_inputs()
    if _missing:
        _SKIP_REASON = "missing demo inputs: " + ", ".join(_missing)

pytestmark = pytest.mark.skipif(
    _SKIP_REASON is not None, reason=_SKIP_REASON or "integration prerequisites met"
)


def _build_pileups(out_dir: Path) -> dict[str, Path]:
    from dimelo import parse_bam as _parse_bam

    pileups: dict[str, Path] = {}
    for sample_id, bam_path in (
        ("barcode17_targeted", BARCODE17_BAM),
        ("barcode18_control", BARCODE18_BAM),
    ):
        pileup_path, _ = _parse_bam.pileup(
            bam_path,
            output_name=f"{sample_id}_A",
            output_directory=out_dir,
            ref_genome=REF_GENOME,
            regions=[ANALYSIS_REGION],
            motifs=["A,0"],
            thresh=0.5,
            cores=2,
            log=False,
            cleanup=True,
            overwrite=True,
            quiet=True,
            override_checks=True,
        )
        pileups[sample_id] = pileup_path
    return pileups


def test_parse_pileup_dmr_pair_and_prep(tmp_path):
    from dimelo import dmr

    pileups = _build_pileups(tmp_path)
    for pileup_path in pileups.values():
        assert pileup_path.exists()
        assert pileup_path.with_suffix(pileup_path.suffix + ".tbi").exists()

    result = dmr.run_dmr_pair(
        control_bed_methyl=pileups["barcode18_control"],
        experiment_bed_methyl=pileups["barcode17_targeted"],
        ref_genome=REF_GENOME,
        out_path=tmp_path / "pair.sites.bed",
        segment_path=tmp_path / "pair.segments.bed",
        bases=["A"],
        min_valid_coverage=1,
        threads=2,
        pvalue_max=0.05,
        abs_effect_size_min=0.05,
    )

    # Real single-site modkit output carries these columns.
    assert not result.sites.empty
    assert {"map_pvalue", "effect_size"}.issubset(result.sites.columns)
    # high-confidence sites are a subset of all sites.
    assert len(result.high_confidence_sites) <= len(result.sites)
    assert result.segments is not None

    # New prep helper consumes the real modkit columns and derives features.
    site_payload = plotting.prepare_dmr_site_data(result=result)
    site_table = site_payload["site_table"]
    for column in ("contig", "midpoint", "span", "neg_log10_pvalue", "is_high_confidence"):
        assert column in site_table.columns
    assert site_table["is_high_confidence"].sum() == len(result.high_confidence_sites)
    assert site_payload["metadata"]["n_high_confidence"] == len(
        result.high_confidence_sites
    )

    segment_payload = plotting.prepare_dmr_segment_data(result=result)
    assert segment_payload["metadata"]["has_segments"] is True

    # Renderers run headlessly on the real payloads.
    import matplotlib

    matplotlib.use("Agg")
    site_fig, _ = plotting_matplotlib.plot_dmr_site_matplotlib(site_payload)
    assert site_fig is not None
    segment_fig, _ = plotting_matplotlib.plot_dmr_segment_matplotlib(segment_payload)
    assert segment_fig is not None


def test_parse_pileup_dmr_multi_and_summary(tmp_path):
    from dimelo import dmr

    pileups = _build_pileups(tmp_path)
    region_bed = tmp_path / "region.bed"
    chrom, span = ANALYSIS_REGION.split(":")
    start, end = span.split("-")
    region_bed.write_text(f"{chrom}\t{start}\t{end}\n")

    multi_result = dmr.run_dmr_multi(
        samples={
            "barcode18_control": pileups["barcode18_control"],
            "barcode17_targeted": pileups["barcode17_targeted"],
        },
        regions_bed=region_bed,
        ref_genome=REF_GENOME,
        out_dir=tmp_path / "multi",
        bases=["A"],
        min_valid_coverage=1,
        threads=2,
        prefix="barcode17_18_A",
    )

    assert not multi_result.pair_files.empty

    summary_payload = plotting.prepare_dmr_multi_summary_data(result=multi_result)
    summary = summary_payload["summary_table"]
    assert not summary.empty
    assert (summary["n_sites"] >= 0).all()
    assert {"pair_name", "n_sites", "n_significant"}.issubset(summary.columns)

    import matplotlib

    matplotlib.use("Agg")
    fig, _ = plotting_matplotlib.plot_dmr_multi_summary_matplotlib(summary_payload)
    assert fig is not None
