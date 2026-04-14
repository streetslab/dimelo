import numpy as np

from dimelo import load_processed


def test_counts_from_pileup_alias_forwards_arguments(monkeypatch):
    captured = {}

    def _fake_counts_from_bedmethyl(**kwargs):
        captured.update(kwargs)
        return (3, 9)

    monkeypatch.setattr(
        load_processed,
        "pileup_counts_from_bedmethyl",
        _fake_counts_from_bedmethyl,
    )

    result = load_processed.counts_from_pileup(
        bedmethyl_file="reads.bed.gz",
        motif="A,0",
        regions="chr1:10-20",
        window_size=50,
        single_strand=True,
        quiet=True,
        cores=2,
        chunk_size=123,
    )

    assert result == (3, 9)
    assert captured == {
        "bedmethyl_file": "reads.bed.gz",
        "motif": "A,0",
        "regions": "chr1:10-20",
        "window_size": 50,
        "single_strand": True,
        "quiet": True,
        "cores": 2,
        "chunk_size": 123,
    }


def test_vectors_from_pileup_alias_forwards_arguments(monkeypatch):
    captured = {}
    expected = (np.array([1, 2, 3]), np.array([4, 5, 6]))

    def _fake_vectors_from_bedmethyl(**kwargs):
        captured.update(kwargs)
        return expected

    monkeypatch.setattr(
        load_processed,
        "pileup_vectors_from_bedmethyl",
        _fake_vectors_from_bedmethyl,
    )

    result = load_processed.vectors_from_pileup(
        bedmethyl_file="reads.bed.gz",
        motif="CG,0",
        regions="chr2:100-200",
        window_size=25,
        single_strand=False,
        regions_5to3prime=True,
        quiet=False,
        cores=4,
        chunk_size=999,
    )

    np.testing.assert_array_equal(result[0], expected[0])
    np.testing.assert_array_equal(result[1], expected[1])
    assert captured == {
        "bedmethyl_file": "reads.bed.gz",
        "motif": "CG,0",
        "regions": "chr2:100-200",
        "window_size": 25,
        "single_strand": False,
        "regions_5to3prime": True,
        "quiet": False,
        "cores": 4,
        "chunk_size": 999,
    }
