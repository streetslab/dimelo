from pathlib import Path

import pytest

from dimelo import export


def _make_bedmethyl_row(
    *,
    chrom: str,
    start: int,
    strand: str,
    modified: int,
    valid: int,
) -> str:
    end = start + 1
    other = max(valid - modified, 0)
    # modkit-compatible 18-column row
    return (
        f"{chrom}\t{start}\t{end}\tY\t{valid}\t{strand}\t{start}\t{end}\t255,0,0\t"
        f"{valid}\t0.00\t{modified}\t{other}\t0\t0\t0\t0\t0"
    )


def test_pileup_to_bigwig_sorts_and_collapses_duplicate_positions(
    monkeypatch, tmp_path
):
    rows_chr1 = [
        _make_bedmethyl_row(chrom="chr1", start=100, strand="+", modified=1, valid=2),
        _make_bedmethyl_row(chrom="chr1", start=99, strand="+", modified=1, valid=1),
        _make_bedmethyl_row(chrom="chr1", start=100, strand="-", modified=2, valid=2),
        _make_bedmethyl_row(chrom="chr1", start=101, strand="+", modified=1, valid=2),
    ]

    class _FakeTabixFile:
        def __init__(self, *_args, **_kwargs):
            self.contigs = ["chr1"]

        def fetch(self, contig):
            assert contig == "chr1"
            return iter(rows_chr1)

    class _FakeFastaFile:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def get_reference_length(self, contig):
            assert contig == "chr1"
            return 1000

    class _FakeBigWig:
        def __init__(self):
            self.header = None
            self.entries = []

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def addHeader(self, header):
            self.header = header

        def addEntries(self, contigs, starts, *, ends, values):
            self.entries.extend(zip(contigs, starts, ends, values, strict=False))

    fake_bw = _FakeBigWig()

    monkeypatch.setattr(
        export.utils,
        "sanitize_path_args",
        lambda bedmethyl_file, bigwig_file, ref_genome: (
            Path(bedmethyl_file),
            Path(bigwig_file) if bigwig_file is not None else None,
            Path(ref_genome) if ref_genome is not None else None,
        ),
    )
    monkeypatch.setattr(export.pysam, "TabixFile", _FakeTabixFile)
    monkeypatch.setattr(export.pysam, "FastaFile", _FakeFastaFile)
    monkeypatch.setattr(export.pyBigWig, "open", lambda *_args, **_kwargs: fake_bw)

    output_bigwig = tmp_path / "out.bw"
    export.pileup_to_bigwig(
        bedmethyl_file=tmp_path / "in.bed.gz",
        motif="A,0",
        bigwig_file=output_bigwig,
        ref_genome=tmp_path / "ref.fa",
        strand=".",
        chunk_size=2,
    )

    assert fake_bw.header == [("chr1", 1000)]
    starts = [entry[1] for entry in fake_bw.entries]
    assert starts == [99, 100, 101]
    values_by_start = {entry[1]: entry[3] for entry in fake_bw.entries}
    assert values_by_start[99] == pytest.approx(1.0)
    # 100 has two rows merged: (1+2)/(2+2) = 0.75
    assert values_by_start[100] == pytest.approx(0.75)
    assert values_by_start[101] == pytest.approx(0.5)


def test_pileup_to_bigwig_sorts_rows_buffered_before_late_inversion(
    monkeypatch, tmp_path
):
    rows_chr1 = [
        _make_bedmethyl_row(chrom="chr1", start=100, strand="+", modified=1, valid=2),
        _make_bedmethyl_row(chrom="chr1", start=101, strand="+", modified=1, valid=4),
        _make_bedmethyl_row(chrom="chr1", start=99, strand="+", modified=1, valid=1),
        _make_bedmethyl_row(chrom="chr1", start=100, strand="-", modified=2, valid=2),
    ]

    class _FakeTabixFile:
        def __init__(self, *_args, **_kwargs):
            self.contigs = ["chr1"]

        def fetch(self, contig):
            assert contig == "chr1"
            return iter(rows_chr1)

    class _FakeFastaFile:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def get_reference_length(self, contig):
            assert contig == "chr1"
            return 1000

    class _FakeBigWig:
        def __init__(self):
            self.header = None
            self.entries = []

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def addHeader(self, header):
            self.header = header

        def addEntries(self, contigs, starts, *, ends, values):
            self.entries.extend(zip(contigs, starts, ends, values, strict=False))

    fake_bw = _FakeBigWig()

    monkeypatch.setattr(
        export.utils,
        "sanitize_path_args",
        lambda bedmethyl_file, bigwig_file, ref_genome: (
            Path(bedmethyl_file),
            Path(bigwig_file) if bigwig_file is not None else None,
            Path(ref_genome) if ref_genome is not None else None,
        ),
    )
    monkeypatch.setattr(export.pysam, "TabixFile", _FakeTabixFile)
    monkeypatch.setattr(export.pysam, "FastaFile", _FakeFastaFile)
    monkeypatch.setattr(export.pyBigWig, "open", lambda *_args, **_kwargs: fake_bw)

    output_bigwig = tmp_path / "out.bw"
    export.pileup_to_bigwig(
        bedmethyl_file=tmp_path / "in.bed.gz",
        motif="A,0",
        bigwig_file=output_bigwig,
        ref_genome=tmp_path / "ref.fa",
        strand=".",
        chunk_size=1,
    )

    assert fake_bw.header == [("chr1", 1000)]
    starts = [entry[1] for entry in fake_bw.entries]
    assert starts == [99, 100, 101]
    values_by_start = {entry[1]: entry[3] for entry in fake_bw.entries}
    assert values_by_start[99] == pytest.approx(1.0)
    # 100 has one row written before inversion plus one later row: (1+2)/(2+2) = 0.75
    assert values_by_start[100] == pytest.approx(0.75)
    assert values_by_start[101] == pytest.approx(0.25)


def test_pileup_to_bigwig_streams_sorted_rows_before_iterator_is_exhausted(
    monkeypatch, tmp_path
):
    rows_chr1 = [
        _make_bedmethyl_row(chrom="chr1", start=100, strand="+", modified=1, valid=2),
        _make_bedmethyl_row(chrom="chr1", start=101, strand="+", modified=1, valid=4),
        _make_bedmethyl_row(chrom="chr1", start=102, strand="+", modified=3, valid=4),
    ]

    class _FakeTabixFile:
        def __init__(self, *_args, **_kwargs):
            self.contigs = ["chr1"]
            self.write_fetch_exhausted = True

        def fetch(self, contig):
            assert contig == "chr1"
            self.write_fetch_exhausted = False
            try:
                yield from rows_chr1
            finally:
                self.write_fetch_exhausted = True

    class _FakeFastaFile:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def get_reference_length(self, contig):
            assert contig == "chr1"
            return 1000

    class _FakeBigWig:
        def __init__(self, tabix):
            self.header = None
            self.entries = []
            self.add_entries_exhausted_states = []
            self.tabix = tabix

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def addHeader(self, header):
            self.header = header

        def addEntries(self, contigs, starts, *, ends, values):
            self.add_entries_exhausted_states.append(self.tabix.write_fetch_exhausted)
            self.entries.extend(zip(contigs, starts, ends, values, strict=False))

    fake_tabix = _FakeTabixFile()
    fake_bw = _FakeBigWig(fake_tabix)

    monkeypatch.setattr(
        export.utils,
        "sanitize_path_args",
        lambda bedmethyl_file, bigwig_file, ref_genome: (
            Path(bedmethyl_file),
            Path(bigwig_file) if bigwig_file is not None else None,
            Path(ref_genome) if ref_genome is not None else None,
        ),
    )
    monkeypatch.setattr(export.pysam, "TabixFile", lambda *_args, **_kwargs: fake_tabix)
    monkeypatch.setattr(export.pysam, "FastaFile", _FakeFastaFile)
    monkeypatch.setattr(export.pyBigWig, "open", lambda *_args, **_kwargs: fake_bw)

    export.pileup_to_bigwig(
        bedmethyl_file=tmp_path / "in.bed.gz",
        motif="A,0",
        bigwig_file=tmp_path / "out.bw",
        ref_genome=tmp_path / "ref.fa",
        strand=".",
        chunk_size=1,
    )

    assert fake_bw.add_entries_exhausted_states
    assert fake_bw.add_entries_exhausted_states[0] is False
    assert [entry[1] for entry in fake_bw.entries] == [100, 101, 102]


def test_pileup_to_bigwig_uses_max_end_for_header_without_reference(
    monkeypatch, tmp_path
):
    rows_chr1 = [
        _make_bedmethyl_row(chrom="chr1", start=100, strand="+", modified=1, valid=2),
        _make_bedmethyl_row(chrom="chr1", start=105, strand="+", modified=1, valid=4),
        _make_bedmethyl_row(chrom="chr1", start=99, strand="+", modified=1, valid=1),
    ]

    class _FakeTabixFile:
        def __init__(self, *_args, **_kwargs):
            self.contigs = ["chr1"]

        def fetch(self, contig):
            assert contig == "chr1"
            return iter(rows_chr1)

    class _FakeBigWig:
        def __init__(self):
            self.header = None
            self.entries = []

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def addHeader(self, header):
            self.header = header

        def addEntries(self, contigs, starts, *, ends, values):
            self.entries.extend(zip(contigs, starts, ends, values, strict=False))

    fake_bw = _FakeBigWig()

    monkeypatch.setattr(
        export.utils,
        "sanitize_path_args",
        lambda bedmethyl_file, bigwig_file, ref_genome: (
            Path(bedmethyl_file),
            Path(bigwig_file) if bigwig_file is not None else None,
            Path(ref_genome) if ref_genome is not None else None,
        ),
    )
    monkeypatch.setattr(export.pysam, "TabixFile", _FakeTabixFile)
    monkeypatch.setattr(export.pyBigWig, "open", lambda *_args, **_kwargs: fake_bw)

    export.pileup_to_bigwig(
        bedmethyl_file=tmp_path / "in.bed.gz",
        motif="A,0",
        bigwig_file=tmp_path / "out.bw",
        ref_genome=None,
        strand=".",
        chunk_size=2,
    )

    assert fake_bw.header == [("chr1", 106)]
    assert [entry[1] for entry in fake_bw.entries] == [99, 100, 105]
