import pandas as pd
import pytest

from dimelo import global_analysis
from dimelo.models import SampleSpec


def test_summarize_global_samples_from_pileup(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="s2.h5",
            replicate=2,
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]

    counts = {
        ("s1.bed.gz", "A,0"): (2, 10),
        ("s1.bed.gz", "CG,0"): (1, 4),
        ("s2.bed.gz", "A,0"): (3, 6),
        ("s2.bed.gz", "CG,0"): (0, 0),
    }

    def fake_global_counts_from_bedmethyl(bedmethyl_file, motif, quiet=True):
        return counts[(bedmethyl_file, motif)]

    monkeypatch.setattr(
        global_analysis,
        "_global_counts_from_bedmethyl",
        fake_global_counts_from_bedmethyl,
    )

    result = global_analysis.summarize_global_samples(
        samples=samples,
        motifs=["A,0", "CG,0"],
    )

    expected = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "motif": "A,0",
                "modified_count": 2,
                "valid_count": 10,
                "global_fraction": 0.2,
            },
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "motif": "CG,0",
                "modified_count": 1,
                "valid_count": 4,
                "global_fraction": 0.25,
            },
            {
                "sample_id": "s2",
                "condition": "15min",
                "replicate": 2,
                "motif": "A,0",
                "modified_count": 3,
                "valid_count": 6,
                "global_fraction": 0.5,
            },
            {
                "sample_id": "s2",
                "condition": "15min",
                "replicate": 2,
                "motif": "CG,0",
                "modified_count": 0,
                "valid_count": 0,
                "global_fraction": float("nan"),
            },
        ]
    )

    pd.testing.assert_frame_equal(result, expected)


def test_global_counts_from_bedmethyl_sums_only_matching_rows(monkeypatch):
    rows_by_contig = {
        "chr1": ["row-1", "row-2"],
        "chr2": ["row-3"],
    }

    class FakeTabixFile:
        def __init__(self, path):
            assert path == "pileup.bed.gz"
            self.contigs = tuple(rows_by_contig)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def fetch(self, contig):
            return iter(rows_by_contig[contig])

    parsed_motifs = []

    def fake_parsed_motif(motif):
        parsed_motifs.append(motif)
        return {"motif": motif}

    processed_rows = []
    process_results = {
        "row-1": (True, 10, 3, 8),
        "row-2": (False, 11, 20, 30),
        "row-3": (True, 12, 5, 10),
    }

    def fake_process_pileup_row(row, parsed_motif, region_strand, single_strand=False):
        processed_rows.append((row, parsed_motif, region_strand, single_strand))
        return process_results[row]

    monkeypatch.setattr(global_analysis.pysam, "TabixFile", FakeTabixFile)
    monkeypatch.setattr(global_analysis.utils, "ParsedMotif", fake_parsed_motif)
    monkeypatch.setattr(
        global_analysis.load_processed,
        "process_pileup_row",
        fake_process_pileup_row,
    )

    modified_count, valid_count = global_analysis._global_counts_from_bedmethyl(
        bedmethyl_file="pileup.bed.gz",
        motif="A,0",
    )

    assert (modified_count, valid_count) == (8, 18)
    assert parsed_motifs == ["A,0"]
    assert processed_rows == [
        ("row-1", {"motif": "A,0"}, ".", False),
        ("row-2", {"motif": "A,0"}, ".", False),
        ("row-3", {"motif": "A,0"}, ".", False),
    ]


def test_summarize_global_samples_requires_pileup_path():
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
        )
    ]

    with pytest.raises(ValueError, match="missing metadata\\['pileup_path'\\]"):
        global_analysis.summarize_global_samples(samples=samples, motifs=["A,0"])
