import pytest
from pathlib import Path

from dimelo import parse_bam


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
