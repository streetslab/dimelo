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
                "global_fraction": 0.0,
            },
        ]
    )

    pd.testing.assert_frame_equal(result, expected)


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
