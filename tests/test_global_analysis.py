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

    observed_calls = []

    def fake_global_counts_for_motifs_from_bedmethyl(
        bedmethyl_file, motifs, quiet=True
    ):
        observed_calls.append((bedmethyl_file, tuple(motifs), quiet))
        return {motif: counts[(bedmethyl_file, motif)] for motif in motifs}

    monkeypatch.setattr(
        global_analysis,
        "_global_counts_for_motifs_from_bedmethyl",
        fake_global_counts_for_motifs_from_bedmethyl,
    )
    monkeypatch.setattr(
        global_analysis,
        "_global_counts_from_bedmethyl",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("legacy single-motif counter should not be used here")
        ),
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
    assert observed_calls == [
        ("s1.bed.gz", ("A,0", "CG,0"), True),
        ("s2.bed.gz", ("A,0", "CG,0"), True),
    ]


def test_global_counts_for_motifs_from_bedmethyl_single_pass(monkeypatch):
    rows_by_contig = {
        "chr1": ["row-1", "row-2"],
        "chr2": ["row-3"],
    }

    class FakeTabixFile:
        opened_paths = []
        fetch_calls = []

        def __init__(self, path):
            self.opened_paths.append(path)
            self.contigs = tuple(rows_by_contig)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def fetch(self, contig):
            self.fetch_calls.append(contig)
            return iter(rows_by_contig[contig])

    parsed_motifs = []

    def fake_parsed_motif(motif):
        parsed_motifs.append(motif)
        return {"motif": motif}

    processed_rows = []
    process_results = {
        ("row-1", "A,0"): (True, 10, 3, 8),
        ("row-1", "CG,0"): (False, 10, 99, 99),
        ("row-2", "A,0"): (False, 11, 99, 99),
        ("row-2", "CG,0"): (True, 11, 4, 9),
        ("row-3", "A,0"): (True, 12, 5, 10),
        ("row-3", "CG,0"): (True, 12, 6, 12),
    }

    def fake_process_pileup_row(row, parsed_motif, region_strand, single_strand=False):
        motif = parsed_motif["motif"]
        processed_rows.append((row, motif, region_strand, single_strand))
        return process_results[(row, motif)]

    monkeypatch.setattr(global_analysis.pysam, "TabixFile", FakeTabixFile)
    monkeypatch.setattr(global_analysis.utils, "ParsedMotif", fake_parsed_motif)
    monkeypatch.setattr(
        global_analysis.load_processed,
        "process_pileup_row",
        fake_process_pileup_row,
    )

    counts_by_motif = global_analysis._global_counts_for_motifs_from_bedmethyl(
        bedmethyl_file="pileup.bed.gz",
        motifs=["A,0", "CG,0"],
    )

    assert counts_by_motif == {
        "A,0": (8, 18),
        "CG,0": (10, 21),
    }
    assert parsed_motifs == ["A,0", "CG,0"]
    assert FakeTabixFile.opened_paths == ["pileup.bed.gz"]
    assert FakeTabixFile.fetch_calls == ["chr1", "chr2"]
    assert processed_rows == [
        ("row-1", "A,0", ".", False),
        ("row-1", "CG,0", ".", False),
        ("row-2", "A,0", ".", False),
        ("row-2", "CG,0", ".", False),
        ("row-3", "A,0", ".", False),
        ("row-3", "CG,0", ".", False),
    ]


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


def test_tile_windows_from_genome_sizes_dict():
    windows = global_analysis.tile_genome_windows(
        genome_sizes={"chr1": 2500},
        window_size=1000,
        step_size=500,
    )

    assert windows["window_id"].tolist() == [
        "chr1:0-1000",
        "chr1:500-1500",
        "chr1:1000-2000",
        "chr1:1500-2500",
    ]


def test_build_window_summary_from_regions_to_list(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
    ]
    windows = pd.DataFrame(
        {
            "window_id": ["chr1:0-1000", "chr1:500-1500"],
            "chromosome": ["chr1", "chr1"],
            "start": [0, 500],
            "end": [1000, 1500],
            "strand": [".", "."],
        }
    )

    captured = {}

    def fake_tile_genome_windows(
        genome_sizes,
        window_size,
        step_size,
        include_contigs=None,
        exclude_contigs=None,
    ):
        captured["genome_sizes"] = genome_sizes
        captured["window_size"] = window_size
        captured["step_size"] = step_size
        captured["include_contigs"] = include_contigs
        captured["exclude_contigs"] = exclude_contigs
        return windows

    def fake_regions_to_list(
        function_handle,
        regions,
        window_size=None,
        quiet=True,
        cores=None,
        split_large_regions=False,
    ):
        captured["function_handle"] = function_handle
        captured["regions"] = list(regions)
        captured["window_size_arg"] = window_size
        captured["quiet"] = quiet
        captured["cores"] = cores
        captured["split_large_regions"] = split_large_regions
        return [(5, 10), (0, 0)]

    monkeypatch.setattr(
        global_analysis, "tile_genome_windows", fake_tile_genome_windows
    )
    monkeypatch.setattr(
        global_analysis.load_processed,
        "regions_to_list",
        fake_regions_to_list,
    )

    summary = global_analysis.build_window_summary(
        samples=samples,
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        window_size=1000,
        step_size=500,
    )

    assert captured["genome_sizes"] == {"chr1": 1500}
    assert captured["window_size"] == 1000
    assert captured["step_size"] == 500
    assert captured["regions"] == ["chr1:0-1000,.", "chr1:500-1500,."]
    assert captured["function_handle"].keywords["bedmethyl_file"] == "s1.bed.gz"
    assert captured["function_handle"].keywords["motif"] == "A,0"
    assert summary["window_fraction"].tolist() == [0.5, 0.0]


def test_compute_global_normalization_factors_from_summary():
    summary = pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "condition": ["NS", "15min"],
            "replicate": [None, None],
            "motif": ["A,0", "A,0"],
            "modified_count": [10, 30],
            "valid_count": [100, 100],
            "global_fraction": [0.1, 0.3],
        }
    )

    factors = global_analysis.compute_global_normalization_factors(summary)

    assert set(factors.columns) >= {
        "sample_id",
        "motif",
        "global_fraction",
        "reference_fraction",
        "global_offset",
    }
    assert factors.loc[factors["sample_id"] == "s1", "global_offset"].iloc[
        0
    ] == pytest.approx(-0.1)
    assert factors.loc[factors["sample_id"] == "s2", "global_offset"].iloc[
        0
    ] == pytest.approx(0.1)


def test_run_global_analysis_returns_result(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
    ]

    monkeypatch.setattr(
        global_analysis,
        "summarize_global_samples",
        lambda *, samples, motifs, quiet=True: pd.DataFrame(
            {
                "sample_id": ["s1"],
                "condition": ["NS"],
                "replicate": [None],
                "motif": ["A,0"],
                "modified_count": [10],
                "valid_count": [100],
                "global_fraction": [0.1],
            }
        ),
    )
    monkeypatch.setattr(
        global_analysis,
        "build_window_summary",
        lambda **kwargs: pd.DataFrame({"window_id": ["chr1:0-1000"]}),
    )
    monkeypatch.setattr(
        global_analysis,
        "compute_global_normalization_factors",
        lambda summary: pd.DataFrame({"sample_id": ["s1"], "global_offset": [0.0]}),
    )

    result = global_analysis.run_global_analysis(
        samples=samples,
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=1000,
        step_size=1000,
    )

    assert list(result.summary["sample_id"]) == ["s1"]
    assert list(result.windows["window_id"]) == ["chr1:0-1000"]
    assert "global_fraction_bar" in result.plot_data


def test_run_global_analysis_supports_multiple_motifs(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
    ]
    observed_build_calls = []

    monkeypatch.setattr(
        global_analysis,
        "summarize_global_samples",
        lambda *, samples, motifs, quiet=True: pd.DataFrame(
            {
                "sample_id": ["s1", "s1"],
                "condition": ["NS", "NS"],
                "replicate": [None, None],
                "motif": list(motifs),
                "modified_count": [10, 20],
                "valid_count": [100, 200],
                "global_fraction": [0.1, 0.1],
            }
        ),
    )

    def fake_build_window_summary(**kwargs):
        motifs = list(kwargs["motifs"])
        observed_build_calls.append(motifs)
        return pd.DataFrame(
            {
                "window_id": ["chr1:0-1000"],
                "motif": motifs,
            }
        )

    monkeypatch.setattr(
        global_analysis, "build_window_summary", fake_build_window_summary
    )
    monkeypatch.setattr(
        global_analysis,
        "compute_global_normalization_factors",
        lambda summary: pd.DataFrame(
            {
                "sample_id": ["s1", "s1"],
                "motif": ["A,0", "CG,0"],
                "global_offset": [0.0, 0.0],
            }
        ),
    )

    result = global_analysis.run_global_analysis(
        samples=samples,
        motifs=["A,0", "CG,0"],
        genome_sizes={"chr1": 1000},
        window_size=1000,
        step_size=1000,
    )

    assert observed_build_calls == [["A,0"], ["CG,0"]]
    assert result.windows["motif"].tolist() == ["A,0", "CG,0"]


def test_run_global_analysis_requires_at_least_one_motif():
    with pytest.raises(ValueError, match="at least one motif"):
        global_analysis.run_global_analysis(
            samples=[],
            motifs=[],
            genome_sizes={"chr1": 1000},
            window_size=1000,
            step_size=1000,
        )
