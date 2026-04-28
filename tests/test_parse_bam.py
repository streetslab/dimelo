import gzip
from pathlib import Path

import pytest

from dimelo import parse_bam, run_modkit


class _FakeRead:
    def __init__(self, tags=None):
        self._tags = tags or []

    def to_dict(self):
        return {"tags": self._tags}


class _FakeAlignmentFile:
    def __init__(self, read_count):
        self.read_count = read_count
        self.reads_yielded = 0

    def fetch(self):
        for _ in range(self.read_count):
            self.reads_yielded += 1
            yield _FakeRead()


class _TagSequenceAlignmentFile:
    def __init__(self, reads):
        self._reads = reads

    def fetch(self):
        yield from self._reads


def test_check_bam_format_consumes_at_most_the_first_100_reads(monkeypatch):
    fake_alignment_file = _FakeAlignmentFile(parse_bam.NUM_READS_TO_CHECK + 1)

    monkeypatch.setattr(
        parse_bam.pysam,
        "AlignmentFile",
        lambda *args, **kwargs: fake_alignment_file,
    )

    parse_bam.check_bam_format("dummy.bam", motifs=["A,0", "CG,0"])

    assert fake_alignment_file.reads_yielded == parse_bam.NUM_READS_TO_CHECK


def test_check_bam_format_raises_on_malformed_tags_after_motifs_are_found(monkeypatch):
    fake_alignment_file = _TagSequenceAlignmentFile(
        [
            _FakeRead(["MM:Z:A+a?;C+m?"]),
            _FakeRead(["Mm:Z:A+a?"]),
        ]
    )

    monkeypatch.setattr(
        parse_bam.pysam,
        "AlignmentFile",
        lambda *args, **kwargs: fake_alignment_file,
    )

    with pytest.raises(ValueError, match="Mm and Ml instead of MM and ML"):
        parse_bam.check_bam_format("dummy.bam", motifs=["A,0", "CG,0"])


def test_threads_command_list_defaults_to_all_available_cores(monkeypatch):
    monkeypatch.setattr(parse_bam.multiprocessing, "cpu_count", lambda: 12)

    command = parse_bam._threads_command_list(cores=None, quiet=True)

    assert command == ["--threads", "12"]


def test_threads_command_list_caps_requested_cores_to_available(monkeypatch):
    monkeypatch.setattr(parse_bam.multiprocessing, "cpu_count", lambda: 8)

    command = parse_bam._threads_command_list(cores=32, quiet=True)

    assert command == ["--threads", "8"]


def test_threads_command_list_uses_requested_cores_when_available(monkeypatch):
    monkeypatch.setattr(parse_bam.multiprocessing, "cpu_count", lambda: 16)

    command = parse_bam._threads_command_list(cores=6, quiet=True)

    assert command == ["--threads", "6"]


def test_parse_should_render_live_progress_respects_env_off(monkeypatch):
    monkeypatch.setenv("DIMELO_PROGRESS_MODE", "off")
    monkeypatch.setattr(parse_bam.sys.stderr, "isatty", lambda: True)

    assert parse_bam._should_render_live_progress() is False


def test_parse_should_render_live_progress_auto_disables_notebook(monkeypatch):
    monkeypatch.delenv("DIMELO_PROGRESS_MODE", raising=False)
    monkeypatch.setenv("JPY_PARENT_PID", "12345")
    monkeypatch.setattr(parse_bam.sys.stderr, "isatty", lambda: True)

    assert parse_bam._should_render_live_progress() is False


def test_create_region_command_list_returns_empty_without_regions(tmp_path):
    command, bed_path = parse_bam.create_region_command_list(
        output_path=tmp_path,
        regions=None,
        window_size=None,
    )

    assert command == []
    assert bed_path is None


def test_create_region_command_list_writes_processed_bed(tmp_path):
    command, bed_path = parse_bam.create_region_command_list(
        output_path=tmp_path,
        regions="chr1:10-20",
        window_size=None,
    )

    assert bed_path == Path(tmp_path) / "regions.processed.bed"
    assert command == ["--include-bed", str(bed_path)]
    assert bed_path.exists()


def test_build_pileup_targeting_command_list_routes_to_modified_bases_on_modkit_0_6():
    capabilities = run_modkit.ModkitCapabilities(
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

    command = parse_bam._build_pileup_targeting_command_list(
        motifs=["A,0,Y", "CG,0,Z"],
        capabilities=capabilities,
    )

    assert "--modified-bases" in command
    assert "A:Y" in command
    assert "C:Z" in command
    assert "--cpg" in command
    assert "--motif" not in command


def test_build_pileup_targeting_command_list_uses_version_fallback_for_modkit_0_6():
    capabilities = run_modkit.ModkitCapabilities(
        executable="modkit",
        version_raw="modkit 0.6.1",
        version="0.6.1",
        version_tuple=(0, 6, 1),
        supports_mod_threshold=True,
        supports_mod_thresholds=False,
        # Simulate stale/incorrect capability detection.
        supports_modified_bases=False,
        supports_force_allow_implicit=False,
        supports_extract_subcommands=True,
        extract_supports_reference_long=True,
        extract_supports_reference_short=True,
    )

    command = parse_bam._build_pileup_targeting_command_list(
        motifs=["A,0", "CG,0"],
        capabilities=capabilities,
    )

    assert "--modified-bases" in command
    assert "A:Y" in command
    assert "A:a" in command
    assert "C:Z" in command
    assert "C:m" in command
    assert "--motif" not in command


def test_build_mod_threshold_command_list_routes_flag_by_modkit_capabilities():
    legacy = run_modkit.ModkitCapabilities(
        executable="modkit",
        version_raw="modkit 0.2.4",
        version="0.2.4",
        version_tuple=(0, 2, 4),
        supports_mod_threshold=False,
        supports_mod_thresholds=True,
        supports_modified_bases=False,
        supports_force_allow_implicit=True,
        supports_extract_subcommands=False,
        extract_supports_reference_long=False,
        extract_supports_reference_short=True,
    )
    modern = run_modkit.ModkitCapabilities(
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

    legacy_flags = parse_bam._build_mod_threshold_command_list(
        motifs=["A,0,Y"],
        motif_thresholds={"A,0,Y": 0.75},
        capabilities=legacy,
    )
    modern_flags = parse_bam._build_mod_threshold_command_list(
        motifs=["A,0,Y"],
        motif_thresholds={"A,0,Y": 0.75},
        capabilities=modern,
    )

    assert legacy_flags == ["--mod-thresholds", "Y:0.75"]
    assert modern_flags == ["--mod-threshold", "Y:0.75"]


def test_build_mod_threshold_command_list_rejects_conflicting_thresholds_for_shared_mod_code():
    capabilities = run_modkit.ModkitCapabilities(
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

    with pytest.raises(
        ValueError,
        match="Cannot apply different thresholds to motifs that share mod code",
    ):
        parse_bam._build_mod_threshold_command_list(
            motifs=["A,0,a", "GATC,1,a"],
            motif_thresholds={"A,0,a": 0.7, "GATC,1,a": 0.8},
            capabilities=capabilities,
        )


def test_resolve_motif_thresholds_accepts_scalar_and_per_motif_dict():
    motifs = ["A,0", "CG,0", "GCH,1"]

    scalar = parse_bam._resolve_motif_thresholds(
        motifs=motifs,
        thresh=0.7,
        quiet=True,
    )
    assert scalar == {"A,0": 0.7, "CG,0": 0.7, "GCH,1": 0.7}

    by_motif = parse_bam._resolve_motif_thresholds(
        motifs=motifs,
        thresh={"A,0": 0.6, "CG,0": 0.8, "default": 0.9},
        quiet=True,
    )
    assert by_motif == {"A,0": 0.6, "CG,0": 0.8, "GCH,1": 0.9}


def test_resolve_motif_thresholds_supports_canonical_key_alias():
    thresholds = parse_bam._resolve_motif_thresholds(
        motifs=["A,0,a", "CG,0,m"],
        thresh={"A,0": 0.65, "CG,0": 0.75},
        quiet=True,
    )
    assert thresholds == {"A,0,a": 0.65, "CG,0,m": 0.75}


def test_build_extract_command_prefix_uses_full_subcommand_for_modkit_0_6():
    modern = run_modkit.ModkitCapabilities(
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
    legacy = run_modkit.ModkitCapabilities(
        executable="modkit",
        version_raw="modkit 0.2.4",
        version="0.2.4",
        version_tuple=(0, 2, 4),
        supports_mod_threshold=False,
        supports_mod_thresholds=True,
        supports_modified_bases=False,
        supports_force_allow_implicit=True,
        supports_extract_subcommands=False,
        extract_supports_reference_long=False,
        extract_supports_reference_short=True,
    )

    modern_prefix = parse_bam._build_extract_command_prefix(
        input_file=Path("input.bam"),
        output_txt=Path("out.txt"),
        capabilities=modern,
    )
    legacy_prefix = parse_bam._build_extract_command_prefix(
        input_file=Path("input.bam"),
        output_txt=Path("out.txt"),
        capabilities=legacy,
    )

    assert modern_prefix == [
        "modkit",
        "extract",
        "full",
        Path("input.bam"),
        Path("out.txt"),
    ]
    assert legacy_prefix == ["modkit", "extract", Path("input.bam"), Path("out.txt")]


def test_build_extract_command_prefix_uses_version_fallback_for_modkit_0_6():
    stale = run_modkit.ModkitCapabilities(
        executable="modkit",
        version_raw="modkit 0.6.1",
        version="0.6.1",
        version_tuple=(0, 6, 1),
        supports_mod_threshold=True,
        supports_mod_thresholds=False,
        supports_modified_bases=True,
        supports_force_allow_implicit=False,
        # Simulate stale capability detection.
        supports_extract_subcommands=False,
        extract_supports_reference_long=False,
        extract_supports_reference_short=False,
    )

    prefix = parse_bam._build_extract_command_prefix(
        input_file=Path("input.bam"),
        output_txt=Path("out.txt"),
        capabilities=stale,
    )
    ref_flags = parse_bam._build_extract_reference_command_list(
        ref_genome=Path("ref.fa"),
        capabilities=stale,
    )

    assert prefix == ["modkit", "extract", "full", Path("input.bam"), Path("out.txt")]
    assert ref_flags == ["--reference", Path("ref.fa")]


def test_prep_output_directory_respects_overwrite_flag(tmp_path):
    _, (output_file,) = parse_bam.prep_output_directory(
        output_directory=tmp_path,
        output_name="demo",
        input_file=tmp_path,
        output_file_names=["result.txt"],
        overwrite=True,
    )
    output_file.write_text("keep-me")

    parse_bam.prep_output_directory(
        output_directory=tmp_path,
        output_name="demo",
        input_file=tmp_path,
        output_file_names=["result.txt"],
        overwrite=False,
    )
    assert output_file.exists()

    parse_bam.prep_output_directory(
        output_directory=tmp_path,
        output_name="demo",
        input_file=tmp_path,
        output_file_names=["result.txt"],
        overwrite=True,
    )
    assert not output_file.exists()


def test_extract_overwrite_false_reuses_existing_outputs(tmp_path, monkeypatch):
    output_path = tmp_path / "existing_extract"
    output_path.mkdir(parents=True)
    output_reads_path = output_path / "reads.combined_basemods.h5"
    output_reads_path.write_bytes(b"")
    processed_regions_path = output_path / "regions.processed.bed"
    processed_regions_path.write_text("chr1\t0\t1\n")

    def _unexpected_modkit(*args, **kwargs):
        raise AssertionError("modkit should not be called when outputs are reused")

    monkeypatch.setattr(run_modkit, "_ensure_modkit_available", _unexpected_modkit)

    output_reads, processed_regions = parse_bam.extract(
        input_file=tmp_path / "missing.bam",
        output_name="existing_extract",
        ref_genome=tmp_path / "missing.fa",
        output_directory=tmp_path,
        regions="chr1:1-2",
        overwrite=False,
        quiet=True,
    )

    assert output_reads == output_reads_path
    assert processed_regions == processed_regions_path


def test_extract_overwrite_false_raises_on_partial_existing_outputs(
    tmp_path, monkeypatch
):
    output_path = tmp_path / "partial_extract"
    output_path.mkdir(parents=True)
    output_reads_path = output_path / "reads.combined_basemods.h5"
    output_reads_path.write_bytes(b"")

    def _unexpected_modkit(*args, **kwargs):
        raise AssertionError("modkit should not be called before conflict check")

    monkeypatch.setattr(run_modkit, "_ensure_modkit_available", _unexpected_modkit)

    with pytest.raises(FileExistsError, match="overwrite=False"):
        parse_bam.extract(
            input_file=tmp_path / "missing.bam",
            output_name="partial_extract",
            ref_genome=tmp_path / "missing.fa",
            output_directory=tmp_path,
            regions="chr1:1-2",
            overwrite=False,
            quiet=True,
        )


def test_pileup_overwrite_false_reuses_existing_outputs(tmp_path, monkeypatch):
    output_path = tmp_path / "existing_pileup"
    output_path.mkdir(parents=True)
    output_pileup_path = output_path / "pileup.sorted.bed.gz"
    output_pileup_path.write_bytes(b"")
    output_pileup_tbi_path = output_path / "pileup.sorted.bed.gz.tbi"
    output_pileup_tbi_path.write_bytes(b"")
    processed_regions_path = output_path / "regions.processed.bed"
    processed_regions_path.write_text("chr1\t0\t1\n")

    def _unexpected_modkit(*args, **kwargs):
        raise AssertionError("modkit should not be called when outputs are reused")

    monkeypatch.setattr(run_modkit, "_ensure_modkit_available", _unexpected_modkit)

    output_pileup, processed_regions = parse_bam.pileup(
        input_file=tmp_path / "missing.bam",
        output_name="existing_pileup",
        ref_genome=tmp_path / "missing.fa",
        output_directory=tmp_path,
        regions="chr1:1-2",
        overwrite=False,
        quiet=True,
    )

    assert output_pileup == output_pileup_path
    assert processed_regions == processed_regions_path


def test_pileup_preserves_explicit_threshold_empty_output_without_retry(
    tmp_path, monkeypatch
):
    calls = []

    capabilities = run_modkit.ModkitCapabilities(
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

    monkeypatch.setattr(
        run_modkit, "_ensure_modkit_available", lambda **_: capabilities
    )
    monkeypatch.setattr(parse_bam, "verify_inputs", lambda *args, **kwargs: None)

    def _fake_run_with_progress_bars(*, command_list, **kwargs):
        calls.append(command_list)
        output_bed = Path(command_list[3])
        output_bed.write_text("")

    monkeypatch.setattr(
        run_modkit, "run_with_progress_bars", _fake_run_with_progress_bars
    )

    output_pileup, _ = parse_bam.pileup(
        input_file=tmp_path / "missing.bam",
        output_name="test_out",
        ref_genome=tmp_path / "missing.fa",
        output_directory=tmp_path,
        motifs=["A,0"],
        thresh=190,
        regions=None,
        quiet=True,
        overwrite=True,
    )

    assert len(calls) == 1
    assert "--mod-threshold" in calls[0]

    with gzip.open(output_pileup, "rt") as handle:
        rows = [line for line in handle if line.strip()]
    assert rows == []


def test_group_motifs_for_pileup_splits_multi_context_on_modkit_0_6():
    capabilities = run_modkit.ModkitCapabilities(
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
    groups = parse_bam._group_motifs_for_pileup(["A,0", "CG,0"], capabilities)
    assert groups == [("A,0", ["A,0"]), ("CG,0", ["CG,0"])]


def test_group_motifs_for_pileup_keeps_combined_on_legacy_modkit():
    capabilities = run_modkit.ModkitCapabilities(
        executable="modkit",
        version_raw="modkit 0.2.4",
        version="0.2.4",
        version_tuple=(0, 2, 4),
        supports_mod_threshold=False,
        supports_mod_thresholds=True,
        supports_modified_bases=False,
        supports_force_allow_implicit=True,
        supports_extract_subcommands=False,
        extract_supports_reference_long=False,
        extract_supports_reference_short=True,
    )
    groups = parse_bam._group_motifs_for_pileup(["A,0", "CG,0"], capabilities)
    assert groups == [("combined", ["A,0", "CG,0"])]


def test_pileup_multi_motif_split_merges_outputs_on_modkit_0_6(tmp_path, monkeypatch):
    calls = []
    capabilities = run_modkit.ModkitCapabilities(
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

    monkeypatch.setattr(
        run_modkit, "_ensure_modkit_available", lambda **_: capabilities
    )
    monkeypatch.setattr(parse_bam, "verify_inputs", lambda *args, **kwargs: None)

    def _fake_run_with_progress_bars(*, command_list, **kwargs):
        calls.append(command_list)
        output_bed = Path(command_list[3])
        modified_bases = []
        for i, token in enumerate(command_list):
            if token == "--modified-bases":
                modified_bases.append(command_list[i + 1])
        if any(base.startswith("A:") for base in modified_bases):
            output_bed.write_text(
                "chr1\t10\t11\ta\t1\t+\t10\t11\t255,0,0\t1\t0.0\t0\t1\t0\t0\t0\t0\t0\n"
            )
        elif any(base.startswith("C:") for base in modified_bases):
            output_bed.write_text(
                "chr1\t12\t13\tm\t1\t+\t12\t13\t255,0,0\t1\t0.0\t0\t1\t0\t0\t0\t0\t0\n"
            )
        else:
            output_bed.write_text("")

    monkeypatch.setattr(
        run_modkit, "run_with_progress_bars", _fake_run_with_progress_bars
    )

    output_pileup, _ = parse_bam.pileup(
        input_file=tmp_path / "missing.bam",
        output_name="test_out_multi",
        ref_genome=tmp_path / "missing.fa",
        output_directory=tmp_path,
        motifs=["A,0", "CG,0"],
        thresh=None,
        regions=None,
        quiet=True,
        overwrite=True,
    )

    assert len(calls) == 2
    with gzip.open(output_pileup, "rt") as handle:
        rows = [line for line in handle if line.strip()]
    assert len(rows) == 2
