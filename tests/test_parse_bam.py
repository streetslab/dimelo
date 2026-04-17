import pytest
from pathlib import Path

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
        for read in self._reads:
            yield read


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

    assert modern_prefix == ["modkit", "extract", "full", Path("input.bam"), Path("out.txt")]
    assert legacy_prefix == ["modkit", "extract", Path("input.bam"), Path("out.txt")]
