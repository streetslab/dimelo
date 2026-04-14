from __future__ import annotations

import gzip
from pathlib import Path

import h5py
import numpy as np
import pytest

from dimelo import parse_bam


class _FakeRead:
    def __init__(
        self,
        *,
        tags: dict[str, object] | None = None,
        query_sequence: str | None = None,
        aligned_pairs: list[tuple[int, int]] | None = None,
        reference_name: str = "chr1",
    ) -> None:
        self._tags = tags or {}
        self.query_sequence = query_sequence
        self._aligned_pairs = aligned_pairs or []
        self.reference_name = reference_name

    def has_tag(self, tag: str) -> bool:
        return tag in self._tags

    def get_tag(self, tag: str):
        return self._tags[tag]

    def get_aligned_pairs(self, matches_only: bool = False):
        assert matches_only is True
        return list(self._aligned_pairs)


class _FakeBam:
    def __init__(self, reads: list[_FakeRead]) -> None:
        self._reads = reads

    def fetch(self):
        return iter(self._reads)


class _FakeFasta:
    def __init__(self, mapping: dict[tuple[str, int, int], str]) -> None:
        self.mapping = mapping
        self.calls: list[tuple[str, int, int]] = []

    def fetch(self, chrom: str, start: int, end: int) -> str:
        key = (chrom, start, end)
        self.calls.append(key)
        return self.mapping[key]


def test_check_bam_format_rejects_lowercase_mm_ml(monkeypatch):
    fake_bam = _FakeBam(reads=[_FakeRead(tags={"Mm": "x"})])
    monkeypatch.setattr(parse_bam.pysam, "AlignmentFile", lambda *_args, **_kwargs: fake_bam)

    with pytest.raises(ValueError, match="Mm and Ml"):
        parse_bam.check_bam_format("ignored.bam", motifs=["A,0,a"])


def test_check_bam_format_rejects_implicit_mm_probability(monkeypatch):
    fake_bam = _FakeBam(reads=[_FakeRead(tags={"MM": "A+a,0;"})])
    monkeypatch.setattr(parse_bam.pysam, "AlignmentFile", lambda *_args, **_kwargs: fake_bam)

    with pytest.raises(ValueError, match="Need \\? or \\."):
        parse_bam.check_bam_format("ignored.bam", motifs=["A,0,a"])


def test_check_bam_format_accepts_expected_mod_code(monkeypatch):
    fake_bam = _FakeBam(reads=[_FakeRead(tags={"MM": "A+a?,0;"})])
    monkeypatch.setattr(parse_bam.pysam, "AlignmentFile", lambda *_args, **_kwargs: fake_bam)

    # Should return without raising because the expected mod code is present.
    parse_bam.check_bam_format("ignored.bam", motifs=["A,0,a"])


def test_resolve_modkit_threads_command_auto(monkeypatch):
    monkeypatch.setattr(parse_bam, "_memory_limited_modkit_threads", lambda cores: cores)
    monkeypatch.setattr(
        parse_bam.utils,
        "cores_to_run",
        lambda cores: 6 if cores is None else min(int(cores), 6),
    )
    assert parse_bam._resolve_modkit_threads_command(cores=None, quiet=True) == [
        "--threads",
        "6",
    ]


def test_resolve_modkit_threads_command_caps_requested_cores(monkeypatch):
    monkeypatch.setattr(parse_bam, "_memory_limited_modkit_threads", lambda cores: cores)
    monkeypatch.setattr(
        parse_bam.utils,
        "cores_to_run",
        lambda cores: 5 if cores is None else min(int(cores), 5),
    )
    assert parse_bam._resolve_modkit_threads_command(cores=12, quiet=True) == [
        "--threads",
        "5",
    ]


def test_get_alignment_quality_fetches_reference_window_per_read(tmp_path: Path, monkeypatch):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\nACGT\n", encoding="utf-8")
    (tmp_path / "ref.fa.fai").write_text("fake-index\n", encoding="utf-8")

    fake_reads = [
        _FakeRead(
            query_sequence="ACG",
            aligned_pairs=[(0, 100), (1, 101), (2, 102)],
            reference_name="chr1",
        ),
        _FakeRead(
            query_sequence="TT",
            aligned_pairs=[(0, 200), (1, 201)],
            reference_name="chr1",
        ),
    ]
    fake_bam = _FakeBam(fake_reads)
    fake_fasta = _FakeFasta(
        mapping={
            ("chr1", 100, 103): "ACG",
            ("chr1", 200, 202): "TA",
        }
    )

    monkeypatch.setattr(
        parse_bam.pysam, "AlignmentFile", lambda *_args, **_kwargs: fake_bam
    )
    monkeypatch.setattr(
        parse_bam.pysam, "FastaFile", lambda *_args, **_kwargs: fake_fasta
    )

    correct, total = parse_bam.get_alignment_quality("ignored.bam", ref)

    assert total == 5
    assert correct == 4
    assert fake_fasta.calls == [("chr1", 100, 103), ("chr1", 200, 202)]


def test_build_compressed_read_vectors_raw_probabilities():
    mod_compressed, valid_compressed = parse_bam._build_compressed_read_vectors(
        valid_coordinates=[0, 2],
        mod_values=[0.5, 1.0],
        read_len=4,
        thresh=None,
        compress_level=1,
    )

    mod_vector = np.frombuffer(gzip.decompress(mod_compressed.tobytes()), dtype=np.uint8)
    valid_vector = np.frombuffer(gzip.decompress(valid_compressed.tobytes()), dtype=np.uint8)

    assert mod_vector.tolist() == [128, 0, 0]
    assert valid_vector.tolist() == [1, 0, 1]


def test_build_compressed_read_vectors_thresholded_binary():
    mod_compressed, valid_compressed = parse_bam._build_compressed_read_vectors(
        valid_coordinates=[1, 3],
        mod_values=[1, 0],
        read_len=5,
        thresh=0.8,
        compress_level=1,
    )

    mod_vector = np.frombuffer(gzip.decompress(mod_compressed.tobytes()), dtype=np.uint8)
    valid_vector = np.frombuffer(gzip.decompress(valid_compressed.tobytes()), dtype=np.uint8)

    assert mod_vector.tolist() == [0, 1, 0, 0]
    assert valid_vector.tolist() == [0, 1, 0, 1]


def test_read_by_base_txt_to_hdf5_round_trip_thresholded(tmp_path: Path):
    input_txt = tmp_path / "reads.A,0,a.txt"
    output_h5 = tmp_path / "reads.h5"
    input_txt.write_text(
        "\n".join(
            [
                "read_id\tread_pos\tref_pos\tchrom\tref_kmer\tref_strand\tref_mod_strand\tfw_soft_clipped_start\tfw_soft_clipped_end\tread_length\tcall_prob\tcall_code\tbase_qual\tref_position\tref_kmer_5mer\tcanonical_base\textra_col",
                "read1\t0\t100\tchr1\tNNNNN\t+\t+\t0\t0\t4\t0.9\ta\t0\t100\tAAAAA\tA\tX",
                "read1\t1\t101\tchr1\tNNNNN\t+\t+\t0\t0\t4\t0.1\ta\t0\t101\tAAAAA\tA\tX",
                "read1\t2\t102\tchr1\tNNNNN\t+\t+\t0\t0\t4\t0.8\tm\t0\t102\tCCCCC\tC\tX",
                "read2\t0\t200\tchr2\tNNNNN\t+\t+\t0\t0\t3\t0.7\ta\t0\t200\tAAAAA\tA\tX",
                "",
            ]
        ),
        encoding="utf-8",
    )

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=input_txt,
        output_h5=output_h5,
        motif="A,0,a",
        thresh=0.5,
        quiet=True,
        chunk_size=1,
    )

    with h5py.File(output_h5, "r") as h5:
        assert h5["read_name"].shape[0] == 2
        assert h5["read_name"][0].decode("utf-8") == "read1"
        assert h5["read_name"][1].decode("utf-8") == "read2"
        assert h5["chromosome"][0].decode("utf-8") == "chr1"
        assert h5["chromosome"][1].decode("utf-8") == "chr2"
        assert h5["read_start"][0] == 100
        assert h5["read_end"][0] == 104
        assert h5["read_start"][1] == 200
        assert h5["read_end"][1] == 203

        read1_mod = np.frombuffer(
            gzip.decompress(h5["mod_vector"][0].tobytes()), dtype=np.uint8
        )
        read1_val = np.frombuffer(
            gzip.decompress(h5["val_vector"][0].tobytes()), dtype=np.uint8
        )
        read2_mod = np.frombuffer(
            gzip.decompress(h5["mod_vector"][1].tobytes()), dtype=np.uint8
        )
        read2_val = np.frombuffer(
            gzip.decompress(h5["val_vector"][1].tobytes()), dtype=np.uint8
        )

        assert read1_mod.tolist() == [1, 0]
        assert read1_val.tolist() == [1, 1]
        assert read2_mod.tolist() == [1]
        assert read2_val.tolist() == [1]


def test_read_by_base_txt_to_hdf5_appends_across_calls(tmp_path: Path):
    first_txt = tmp_path / "reads_first.A,0,a.txt"
    second_txt = tmp_path / "reads_second.A,0,a.txt"
    output_h5 = tmp_path / "reads.h5"

    header = (
        "read_id\tread_pos\tref_pos\tchrom\tref_kmer\tref_strand\tref_mod_strand\t"
        "fw_soft_clipped_start\tfw_soft_clipped_end\tread_length\tcall_prob\tcall_code\t"
        "base_qual\tref_position\tref_kmer_5mer\tcanonical_base\textra_col"
    )
    first_txt.write_text(
        "\n".join(
            [
                header,
                "read1\t0\t100\tchr1\tNNNNN\t+\t+\t0\t0\t3\t0.9\ta\t0\t100\tAAAAA\tA\tX",
                "",
            ]
        ),
        encoding="utf-8",
    )
    second_txt.write_text(
        "\n".join(
            [
                header,
                "read2\t0\t200\tchr2\tNNNNN\t+\t+\t0\t0\t3\t0.7\ta\t0\t200\tAAAAA\tA\tX",
                "read3\t0\t300\tchr3\tNNNNN\t+\t+\t0\t0\t3\t0.6\ta\t0\t300\tAAAAA\tA\tX",
                "",
            ]
        ),
        encoding="utf-8",
    )

    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=first_txt,
        output_h5=output_h5,
        motif="A,0,a",
        thresh=0.5,
        quiet=True,
    )
    parse_bam.read_by_base_txt_to_hdf5(
        input_txt=second_txt,
        output_h5=output_h5,
        motif="A,0,a",
        thresh=0.5,
        quiet=True,
    )

    with h5py.File(output_h5, "r") as h5:
        assert h5["read_name"].shape[0] == 3
        assert [value.decode("utf-8") for value in h5["read_name"][:]] == [
            "read1",
            "read2",
            "read3",
        ]


def test_read_by_base_txt_to_hdf5_require_motif_context_missing_column_raises(
    tmp_path: Path,
):
    input_txt = tmp_path / "reads_missing_context.txt"
    output_h5 = tmp_path / "reads_missing_context.h5"
    input_txt.write_text(
        "\n".join(
            [
                "read_id\tread_pos\tref_pos\tchrom\tref_kmer\tref_strand\tref_mod_strand\tfw_soft_clipped_start\tfw_soft_clipped_end\tread_length\tcall_prob\tcall_code\tbase_qual\tref_position",
                "read1\t0\t100\tchr1\tNNNNN\t+\t+\t0\t0\t4\t0.9\ta\t0\t100",
                "",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="missing the kmer context column"):
        parse_bam.read_by_base_txt_to_hdf5(
            input_txt=input_txt,
            output_h5=output_h5,
            motif="A,0,a",
            thresh=0.5,
            quiet=True,
            require_motif_context=True,
        )


def test_kmer_matches_parsed_motif_with_iupac():
    parsed = parse_bam.utils.ParsedMotif("GCH,1,m")
    assert parse_bam._kmer_matches_parsed_motif("AGCAT", parsed) is True
    assert parse_bam._kmer_matches_parsed_motif("AGCGT", parsed) is False


def test_verify_inputs_uses_cache(tmp_path: Path, monkeypatch):
    bam = tmp_path / "input.bam"
    ref = tmp_path / "ref.fa"
    bam.write_text("bam", encoding="utf-8")
    ref.write_text(">chr1\nA\n", encoding="utf-8")

    called = {"format": 0, "align": 0}

    def _fake_check(_bam, _motifs):
        called["format"] += 1

    def _fake_align(_bam, _ref):
        called["align"] += 1
        return 90, 100

    monkeypatch.setattr(parse_bam, "check_bam_format", _fake_check)
    monkeypatch.setattr(parse_bam, "get_alignment_quality", _fake_align)

    parse_bam.verify_inputs(bam, ["A,0,a"], ref, quiet=True, cache=True)
    parse_bam.verify_inputs(bam, ["A,0,a"], ref, quiet=True, cache=True)

    assert called["format"] == 1
    assert called["align"] == 1


def test_resolve_modkit_threads_command_memory_caps_auto(monkeypatch):
    monkeypatch.setattr(
        parse_bam.utils,
        "cores_to_run",
        lambda cores: 8 if cores is None else min(int(cores), 8),
    )
    monkeypatch.setattr(parse_bam, "_available_memory_bytes_for_modkit", lambda: 4 * 1024 * 1024 * 1024)
    monkeypatch.setenv("DIMELO_MODKIT_MEMORY_PER_THREAD_BYTES", str(2 * 1024 * 1024 * 1024))
    monkeypatch.setenv("DIMELO_MODKIT_MEMORY_FRACTION", "1.0")

    assert parse_bam._resolve_modkit_threads_command(cores=None, quiet=True) == [
        "--threads",
        "2",
    ]


def test_extract_runs_modkit_once_and_routes_per_motif(tmp_path: Path, monkeypatch):
    bam = tmp_path / "input.bam"
    ref = tmp_path / "ref.fa"
    bam.write_text("bam", encoding="utf-8")
    ref.write_text(">chr1\nA\n", encoding="utf-8")
    output_path = tmp_path / "out"
    output_path.mkdir(parents=True, exist_ok=True)
    output_reads = output_path / "reads.combined_basemods.h5"

    monkeypatch.setattr(
        parse_bam.utils,
        "sanitize_path_args",
        lambda *_args: (bam, ref, tmp_path),
    )
    monkeypatch.setattr(parse_bam, "verify_inputs", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        parse_bam,
        "prep_output_directory",
        lambda **_kwargs: (output_path, [output_reads]),
    )
    monkeypatch.setattr(
        parse_bam,
        "create_region_command_list",
        lambda *_args, **_kwargs: ([], None),
    )
    monkeypatch.setattr(
        parse_bam,
        "_resolve_modkit_threads_command",
        lambda **_kwargs: ["--threads", "2"],
    )

    modkit_commands: list[list[str]] = []

    def _fake_run_with_progress_bars(**kwargs):
        modkit_commands.append(list(kwargs["command_list"]))
        return ""

    routed_reads: list[tuple[str, bool, Path]] = []

    def _fake_read_to_h5(input_txt, _output_h5, motif, _thresh, **kwargs):
        routed_reads.append((motif, kwargs.get("require_motif_context", False), Path(input_txt)))
        return None

    monkeypatch.setattr(parse_bam.run_modkit, "run_with_progress_bars", _fake_run_with_progress_bars)
    monkeypatch.setattr(parse_bam, "read_by_base_txt_to_hdf5", _fake_read_to_h5)

    parse_bam.extract(
        input_file=bam,
        output_name="out",
        ref_genome=ref,
        output_directory=tmp_path,
        motifs=["A,0,a", "A,0,Y", "CG,0,m"],
        quiet=True,
        cleanup=False,
    )

    assert len(modkit_commands) == 1
    command = modkit_commands[0]
    assert command.count("--motif") == 2
    assert len(routed_reads) == 3
    assert all(require_context for _motif, require_context, _input_txt in routed_reads)
    assert len({input_txt for _motif, _require_context, input_txt in routed_reads}) == 1
