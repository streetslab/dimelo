from __future__ import annotations

import pandas as pd
import pytest

from dimelo import global_analysis
from dimelo.models import ContrastSpec, RegionDiscoveryResult, SampleSpec

from dimelo import region_discovery


def _mock_window_summary() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-1000",
                "chromosome": "chr1",
                "start": 0,
                "end": 1000,
                "strand": ".",
                "modified_count": 10,
                "valid_count": 20,
                "window_fraction": 0.5,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-1000",
                "chromosome": "chr1",
                "start": 0,
                "end": 1000,
                "strand": ".",
                "modified_count": 18,
                "valid_count": 20,
                "window_fraction": 0.9,
            },
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:1000-2000",
                "chromosome": "chr1",
                "start": 1000,
                "end": 2000,
                "strand": ".",
                "modified_count": 12,
                "valid_count": 20,
                "window_fraction": 0.6,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:1000-2000",
                "chromosome": "chr1",
                "start": 1000,
                "end": 2000,
                "strand": ".",
                "modified_count": 14,
                "valid_count": 20,
                "window_fraction": 0.7,
            },
        ]
    )


def test_scan_genome_basic_behavior_with_mocked_window_summary(monkeypatch):
    captured = {}

    def fake_build_window_summary(**kwargs):
        captured.update(kwargs)
        return _mock_window_summary()

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    result = region_discovery.scan_genome(
        samples=[
            SampleSpec(
                sample_id="s1",
                condition="NS",
                extract_h5="s1.h5",
                metadata={"pileup_path": "s1.bed.gz"},
            ),
            SampleSpec(
                sample_id="s2",
                condition="treated",
                extract_h5="s2.h5",
                metadata={"pileup_path": "s2.bed.gz"},
            ),
        ],
        motifs=["A,0"],
        genome_sizes={"chr1": 2000},
        window_size=1000,
        step_size=1000,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["treated"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        score="effect_size_only",
        min_coverage=5,
    )

    assert isinstance(result, RegionDiscoveryResult)
    assert captured["motifs"] == ["A,0"]
    assert captured["window_size"] == 1000
    assert captured["step_size"] == 1000
    assert list(result.windows["window_id"]) == ["chr1:0-1000", "chr1:1000-2000"]
    assert result.windows.loc[0, "score_value"] == pytest.approx(0.4)
    assert list(result.hits["rank"]) == [1, 2]
    assert list(result.hits["window_id"]) == ["chr1:0-1000", "chr1:1000-2000"]
    assert set(result.plot_data) >= {"window_score_table", "top_hits_table"}
    assert result.metadata["score"] == "effect_size_only"


def test_scan_genome_filters_low_coverage_windows(monkeypatch):
    def fake_build_window_summary(**kwargs):
        return pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 10,
                    "valid_count": 20,
                    "window_fraction": 0.5,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 18,
                    "valid_count": 20,
                    "window_fraction": 0.9,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 1,
                    "valid_count": 2,
                    "window_fraction": 0.5,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 1,
                    "valid_count": 2,
                    "window_fraction": 0.5,
                },
            ]
        )

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    result = region_discovery.scan_genome(
        samples=[
            SampleSpec(
                sample_id="s1",
                condition="NS",
                extract_h5="s1.h5",
                metadata={"pileup_path": "s1.bed.gz"},
            ),
            SampleSpec(
                sample_id="s2",
                condition="treated",
                extract_h5="s2.h5",
                metadata={"pileup_path": "s2.bed.gz"},
            ),
        ],
        motifs=["A,0"],
        genome_sizes={"chr1": 2000},
        window_size=1000,
        step_size=1000,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["treated"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        score="effect_size_only",
        min_coverage=10,
    )

    assert list(result.hits["window_id"]) == ["chr1:0-1000"]
    low_coverage = result.windows.loc[
        result.windows["window_id"] == "chr1:1000-2000", "score_value"
    ]
    assert low_coverage.isna().all()
    assert result.windows.loc[
        result.windows["window_id"] == "chr1:1000-2000", "rank"
    ].isna().all()


def test_scan_genome_hands_include_and_exclude_contigs_to_window_builder(monkeypatch):
    captured = {}

    def fake_build_window_summary(**kwargs):
        captured.update(kwargs)
        return pd.DataFrame(
            columns=[
                "sample_id",
                "condition",
                "replicate",
                "motif",
                "window_id",
                "chromosome",
                "start",
                "end",
                "strand",
                "modified_count",
                "valid_count",
                "window_fraction",
            ]
        )

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    region_discovery.scan_genome(
        samples=[],
        motifs=["A,0"],
        genome_sizes={"chr1": 2000, "chr2": 1000},
        window_size=1000,
        step_size=500,
        include_contigs=["chr1"],
        exclude_contigs=["chr2"],
        score="effect_size_only",
    )

    assert captured["include_contigs"] == ["chr1"]
    assert captured["exclude_contigs"] == ["chr2"]


def test_scan_genome_beta_binomial_adds_ranking_fields(monkeypatch):
    def fake_build_window_summary(**kwargs):
        return pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 4,
                    "valid_count": 20,
                    "window_fraction": 0.2,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 16,
                    "valid_count": 20,
                    "window_fraction": 0.8,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 10,
                    "valid_count": 20,
                    "window_fraction": 0.5,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 12,
                    "valid_count": 20,
                    "window_fraction": 0.6,
                },
            ]
        )

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    result = region_discovery.scan_genome(
        samples=[
            SampleSpec(
                sample_id="s1",
                condition="NS",
                extract_h5="s1.h5",
                metadata={"pileup_path": "s1.bed.gz"},
            ),
            SampleSpec(
                sample_id="s2",
                condition="treated",
                extract_h5="s2.h5",
                metadata={"pileup_path": "s2.bed.gz"},
            ),
        ],
        motifs=["A,0"],
        genome_sizes={"chr1": 2000},
        window_size=1000,
        step_size=1000,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["treated"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        score="beta_binomial",
        min_coverage=5,
    )

    assert {"score_value", "p_value", "adjusted_p_value", "rank"}.issubset(
        result.windows.columns
    )
    assert list(result.hits["rank"]) == [1, 2]
    assert result.windows.loc[0, "adjusted_p_value"] <= result.windows.loc[1, "adjusted_p_value"]
    assert result.hits.loc[0, "score_value"] >= result.hits.loc[1, "score_value"]


def test_scan_genome_excludes_low_coverage_before_beta_binomial_scoring(monkeypatch):
    observed_p_value_calls = []
    observed_bh_inputs = []

    def fake_build_window_summary(**kwargs):
        return pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 5,
                    "valid_count": 20,
                    "window_fraction": 0.25,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 15,
                    "valid_count": 20,
                    "window_fraction": 0.75,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 1,
                    "valid_count": 1,
                    "window_fraction": 1.0,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 0,
                    "valid_count": 1,
                    "window_fraction": 0.0,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:2000-3000",
                    "chromosome": "chr1",
                    "start": 2000,
                    "end": 3000,
                    "strand": ".",
                    "modified_count": 8,
                    "valid_count": 20,
                    "window_fraction": 0.4,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:2000-3000",
                    "chromosome": "chr1",
                    "start": 2000,
                    "end": 3000,
                    "strand": ".",
                    "modified_count": 12,
                    "valid_count": 20,
                    "window_fraction": 0.6,
                },
            ]
        )

    def fake_p_value(modified_count, valid_count, alpha, beta):
        observed_p_value_calls.append((modified_count, valid_count, alpha, beta))
        return 0.5

    def fake_bh_adjustment(p_values):
        observed_bh_inputs.append(len(p_values))
        return p_values.astype(float)

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)
    monkeypatch.setattr(region_discovery, "_beta_binomial_two_sided_p_value", fake_p_value)
    monkeypatch.setattr(region_discovery, "_adjust_p_values_bh", fake_bh_adjustment)

    result = region_discovery.scan_genome(
        samples=[],
        motifs=["A,0"],
        genome_sizes={"chr1": 3000},
        window_size=1000,
        step_size=1000,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["treated"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        score="beta_binomial",
        min_coverage=10,
    )

    assert len(observed_p_value_calls) == 2
    assert observed_bh_inputs == [2]
    assert list(result.hits["window_id"]) == ["chr1:0-1000", "chr1:2000-3000"]


def test_scan_genome_effect_size_only_ranks_by_absolute_difference(monkeypatch):
    def fake_build_window_summary(**kwargs):
        return pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 19,
                    "valid_count": 20,
                    "window_fraction": 0.95,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 1,
                    "valid_count": 20,
                    "window_fraction": 0.05,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 5,
                    "valid_count": 20,
                    "window_fraction": 0.25,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 15,
                    "valid_count": 20,
                    "window_fraction": 0.75,
                },
            ]
        )

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    result = region_discovery.scan_genome(
        samples=[],
        motifs=["A,0"],
        genome_sizes={"chr1": 2000},
        window_size=1000,
        step_size=1000,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["treated"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        score="effect_size_only",
        min_coverage=1,
    )

    assert list(result.hits["window_id"]) == ["chr1:0-1000", "chr1:1000-2000"]
    assert list(result.hits["score_value"]) == [pytest.approx(0.9), pytest.approx(0.5)]


def test_scan_genome_raises_for_missing_contrast_condition(monkeypatch):
    def fake_build_window_summary(**kwargs):
        return pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 10,
                    "valid_count": 20,
                    "window_fraction": 0.5,
                }
            ]
        )

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    with pytest.raises(ValueError, match="treated"):
        region_discovery.scan_genome(
            samples=[],
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=1000,
            step_size=1000,
            contrast=ContrastSpec(
                mode="pairwise",
                numerator=["treated"],
                denominator=["NS"],
                reference_condition="NS",
            ),
            score="effect_size_only",
            min_coverage=1,
        )


def test_scan_genome_rejects_unimplemented_contrast_modes():
    with pytest.raises(ValueError, match="pairwise/group_vs_group"):
        region_discovery.scan_genome(
            samples=[],
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=1000,
            step_size=1000,
            contrast=ContrastSpec(
                mode="matched_pairwise",
                numerator=["treated"],
                denominator=["NS"],
                pairing_key="donor",
            ),
            score="effect_size_only",
            min_coverage=1,
        )


def test_scan_genome_reranks_hits_after_min_coverage_filtering(monkeypatch):
    def fake_build_window_summary(**kwargs):
        return pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 1,
                    "valid_count": 20,
                    "window_fraction": 0.05,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": ".",
                    "modified_count": 19,
                    "valid_count": 20,
                    "window_fraction": 0.95,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 1,
                    "valid_count": 1,
                    "window_fraction": 1.0,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": ".",
                    "modified_count": 0,
                    "valid_count": 1,
                    "window_fraction": 0.0,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:2000-3000",
                    "chromosome": "chr1",
                    "start": 2000,
                    "end": 3000,
                    "strand": ".",
                    "modified_count": 6,
                    "valid_count": 20,
                    "window_fraction": 0.3,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:2000-3000",
                    "chromosome": "chr1",
                    "start": 2000,
                    "end": 3000,
                    "strand": ".",
                    "modified_count": 14,
                    "valid_count": 20,
                    "window_fraction": 0.7,
                },
            ]
        )

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    result = region_discovery.scan_genome(
        samples=[],
        motifs=["A,0"],
        genome_sizes={"chr1": 3000},
        window_size=1000,
        step_size=1000,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["treated"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        score="effect_size_only",
        min_coverage=10,
    )

    assert list(result.hits["window_id"]) == ["chr1:0-1000", "chr1:2000-3000"]
    assert list(result.hits["rank"]) == [1, 2]


def test_scan_genome_merge_hits_keeps_strand_and_recomputes_window_fields(monkeypatch):
    def fake_build_window_summary(**kwargs):
        return pd.DataFrame(
            [
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": "+",
                    "modified_count": 1,
                    "valid_count": 10,
                    "window_fraction": 0.1,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:0-1000",
                    "chromosome": "chr1",
                    "start": 0,
                    "end": 1000,
                    "strand": "+",
                    "modified_count": 1,
                    "valid_count": 10,
                    "window_fraction": 0.1,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": "+",
                    "modified_count": 4,
                    "valid_count": 10,
                    "window_fraction": 0.4,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:1000-2000",
                    "chromosome": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "strand": "+",
                    "modified_count": 4,
                    "valid_count": 10,
                    "window_fraction": 0.4,
                },
                {
                    "sample_id": "s1",
                    "condition": "NS",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:2000-3000",
                    "chromosome": "chr1",
                    "start": 2000,
                    "end": 3000,
                    "strand": "-",
                    "modified_count": 5,
                    "valid_count": 10,
                    "window_fraction": 0.5,
                },
                {
                    "sample_id": "s2",
                    "condition": "treated",
                    "replicate": 1,
                    "motif": "A,0",
                    "window_id": "chr1:2000-3000",
                    "chromosome": "chr1",
                    "start": 2000,
                    "end": 3000,
                    "strand": "-",
                    "modified_count": 5,
                    "valid_count": 10,
                    "window_fraction": 0.5,
                },
            ]
        )

    monkeypatch.setattr(global_analysis, "build_window_summary", fake_build_window_summary)

    result = region_discovery.scan_genome(
        samples=[],
        motifs=["A,0"],
        genome_sizes={"chr1": 3000},
        window_size=1000,
        step_size=1000,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["treated"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        score="effect_size_only",
        min_coverage=1,
        merge_hits=True,
        merge_distance=1000,
    )

    assert list(result.hits["window_id"]) == ["chr1:0-2000", "chr1:2000-3000"]
    assert list(result.hits["strand"]) == ["+", "-"]
    assert result.hits.loc[0, "window_fraction"] == pytest.approx(0.25)
    assert result.hits.loc[0, "merged_window_count"] == 2
    assert list(result.hits["rank"]) == [1, 2]
