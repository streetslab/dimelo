import pandas as pd
import pytest

from dimelo import region_contrasts
from dimelo.models import ContrastSpec, SampleSpec


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


def test_score_regions_effect_size_only_pairwise_ranks_largest_delta_first(monkeypatch):
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["control"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "n1",
                "condition": "treated",
                "replicate": 1,
                "modified_count": 6,
                "valid_count": 10,
                "mod_fraction": 0.6,
            },
            {
                "region_id": "reg2",
                "chromosome": "chr2",
                "start": 10,
                "end": 20,
                "strand": "-",
                "sample_id": "d1",
                "condition": "control",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
    )

    expected_summary = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "fraction": 0.8,
                "reference_fraction": 0.2,
                "delta_fraction": 0.6,
                "rank": 1,
            },
            {
                "region_id": "reg2",
                "fraction": 0.6,
                "reference_fraction": 0.4,
                "delta_fraction": 0.2,
                "rank": 2,
            },
        ]
    )

    actual_summary = result.summary[
        [
            "region_id",
            "fraction",
            "reference_fraction",
            "delta_fraction",
            "log2_fc",
            "rank",
        ]
    ].reset_index(drop=True)

    pd.testing.assert_frame_equal(
        actual_summary.drop(columns=["log2_fc"]),
        expected_summary,
    )
    assert actual_summary.loc[0, "log2_fc"] == pytest.approx(1.99999459)
    assert actual_summary.loc[1, "log2_fc"] == pytest.approx(0.58496124)
    assert list(result.regions["region_id"]) == ["reg1", "reg2"]
    assert set(result.plot_data) == {"region_effect_sizes"}
    assert result.metadata["contrast_mode"] == "pairwise"


def test_score_regions_effect_size_only_group_vs_group_pools_conditions(monkeypatch):
    contrast = ContrastSpec(
        mode="group_vs_group",
        numerator=["treated_a", "treated_b"],
        denominator=["control_a", "control_b"],
    )
    evidence = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n1",
                "condition": "treated_a",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "n2",
                "condition": "treated_b",
                "replicate": 1,
                "modified_count": 8,
                "valid_count": 10,
                "mod_fraction": 0.8,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d1",
                "condition": "control_a",
                "replicate": 1,
                "modified_count": 2,
                "valid_count": 10,
                "mod_fraction": 0.2,
            },
            {
                "region_id": "reg1",
                "chromosome": "chr1",
                "start": 0,
                "end": 10,
                "strand": "+",
                "sample_id": "d2",
                "condition": "control_b",
                "replicate": 1,
                "modified_count": 4,
                "valid_count": 10,
                "mod_fraction": 0.4,
            },
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: evidence.copy(),
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
    )

    summary_row = result.summary.iloc[0]
    assert summary_row["region_id"] == "reg1"
    assert summary_row["fraction"] == pytest.approx(0.6)
    assert summary_row["reference_fraction"] == pytest.approx(0.3)
    assert summary_row["delta_fraction"] == pytest.approx(0.3)
    assert summary_row["numerator_replicate_n"] == 2
    assert summary_row["denominator_replicate_n"] == 2


def test_score_regions_effect_size_only_rejects_unsupported_contrast_mode(monkeypatch):
    contrast = ContrastSpec(
        mode="matched_pairwise",
        numerator=["treated"],
        denominator=["control"],
        pairing_key="donor_id",
    )

    monkeypatch.setattr(
        region_contrasts,
        "build_region_evidence_table",
        lambda **kwargs: pd.DataFrame(),
    )

    with pytest.raises(NotImplementedError, match="matched_pairwise"):
        region_contrasts.score_regions(
            samples=[],
            regions="regions.bed",
            motifs=["A,0"],
            contrast=contrast,
            test="effect_size_only",
        )
