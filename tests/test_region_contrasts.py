import pandas as pd
import pytest

from dimelo import region_contrasts
from dimelo.models import SampleSpec


def test_validate_supported_v1_combination():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="ensemble_region",
        representation="modified_fraction",
        signal_source="pileup_counts",
        test="beta_binomial",
    )


def test_validate_rejects_unsupported_single_read_beta_binomial():
    with pytest.raises(ValueError, match="analysis_unit='ensemble_region'"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="pileup_counts",
            test="beta_binomial",
        )


def test_build_region_evidence_table_from_pileup_counts(monkeypatch):
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

    monkeypatch.setattr(
        region_contrasts.utils,
        "regions_dict_from_input",
        lambda regions, window_size=None: {
            "chr1": [(0, 10, "+")],
            "chr2": [(20, 30, "-")],
        },
    )

    counts_by_file = {
        "s1.bed.gz": [(2, 10), (0, 0)],
        "s2.bed.gz": [(7, 10), (8, 10)],
    }

    def fake_regions_to_list(
        function_handle,
        regions,
        window_size,
        quiet,
        cores,
        split_large_regions=False,
    ):
        return counts_by_file[function_handle.keywords["bedmethyl_file"]]

    monkeypatch.setattr(
        region_contrasts.load_processed,
        "regions_to_list",
        fake_regions_to_list,
    )

    evidence = region_contrasts.build_region_evidence_table(
        samples=samples,
        regions="matched.bed",
        motifs=["A,0"],
    )

    expected = pd.DataFrame(
        [
            {
                "region_id": "chr1:0-10,+",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "chr2:20-30,-",
                "chromosome": "chr2",
                "start": 20,
                "end": 30,
                "strand": "-",
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "modified_count": 0,
                "valid_count": 0,
                "mod_fraction": 0.0,
            },
            {
                "region_id": "chr1:0-10,+",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "s2",
                "condition": "15min",
                "replicate": 2,
                "modified_count": 7,
                "valid_count": 10,
                "mod_fraction": 0.7,
            },
            {
                "region_id": "chr2:20-30,-",
                "chromosome": "chr2",
                "start": 20,
                "end": 30,
                "strand": "-",
                "sample_id": "s2",
                "condition": "15min",
                "replicate": 2,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
        ]
    )

    pd.testing.assert_frame_equal(evidence.reset_index(drop=True), expected)
