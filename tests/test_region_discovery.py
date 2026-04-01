from __future__ import annotations

import pandas as pd
import pytest

from dimelo import global_analysis
from dimelo.models import ContrastSpec, RegionDiscoveryResult, SampleSpec

from dimelo import region_contrasts, region_discovery


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


def _mock_paired_window_summary() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "sample_id": "t1",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 10,
                "valid_count": 20,
                "window_fraction": 0.5,
            },
            {
                "sample_id": "t1_rep2",
                "condition": "targeting",
                "replicate": 2,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 10,
                "window_fraction": 0.4,
            },
            {
                "sample_id": "d1",
                "condition": "nontargeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 20,
                "window_fraction": 0.2,
            },
            {
                "sample_id": "t2",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 6,
                "valid_count": 20,
                "window_fraction": 0.3,
            },
            {
                "sample_id": "d2",
                "condition": "nontargeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 2,
                "valid_count": 20,
                "window_fraction": 0.1,
            },
            {
                "sample_id": "t3",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 100,
                "valid_count": 100,
                "window_fraction": 1.0,
            },
        ]
    )


def _mock_paired_pairwise_window_summary() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "sample_id": "t1",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 5,
                "window_fraction": 0.8,
            },
            {
                "sample_id": "t1_rep2",
                "condition": "targeting",
                "replicate": 2,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 5,
                "window_fraction": 0.8,
            },
            {
                "sample_id": "d1",
                "condition": "nontargeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 10,
                "window_fraction": 0.4,
            },
            {
                "sample_id": "t2",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 7,
                "valid_count": 10,
                "window_fraction": 0.7,
            },
            {
                "sample_id": "d2",
                "condition": "nontargeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 10,
                "window_fraction": 0.4,
            },
            {
                "sample_id": "t1",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:500-1000",
                "chromosome": "chr1",
                "start": 500,
                "end": 1000,
                "strand": ".",
                "modified_count": 1,
                "valid_count": 5,
                "window_fraction": 0.2,
            },
            {
                "sample_id": "t1_rep2",
                "condition": "targeting",
                "replicate": 2,
                "motif": "A,0",
                "window_id": "chr1:500-1000",
                "chromosome": "chr1",
                "start": 500,
                "end": 1000,
                "strand": ".",
                "modified_count": 1,
                "valid_count": 5,
                "window_fraction": 0.2,
            },
            {
                "sample_id": "d1",
                "condition": "nontargeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:500-1000",
                "chromosome": "chr1",
                "start": 500,
                "end": 1000,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 10,
                "window_fraction": 0.4,
            },
            {
                "sample_id": "t2",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:500-1000",
                "chromosome": "chr1",
                "start": 500,
                "end": 1000,
                "strand": ".",
                "modified_count": 3,
                "valid_count": 10,
                "window_fraction": 0.3,
            },
            {
                "sample_id": "d2",
                "condition": "nontargeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:500-1000",
                "chromosome": "chr1",
                "start": 500,
                "end": 1000,
                "strand": ".",
                "modified_count": 4,
                "valid_count": 10,
                "window_fraction": 0.4,
            },
            {
                "sample_id": "t3",
                "condition": "targeting",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-500",
                "chromosome": "chr1",
                "start": 0,
                "end": 500,
                "strand": ".",
                "modified_count": 5,
                "valid_count": 10,
                "window_fraction": 0.5,
            },
        ]
    )


def _paired_samplespecs() -> list[SampleSpec]:
    return [
        SampleSpec(
            sample_id="t1",
            condition="targeting",
            extract_h5="t1.h5",
            metadata={"pileup_path": "t1.bed.gz", "pair_id": "pair-1"},
        ),
        SampleSpec(
            sample_id="t1_rep2",
            condition="targeting",
            extract_h5="t1_rep2.h5",
            metadata={"pileup_path": "t1_rep2.bed.gz", "pair_id": "pair-1"},
        ),
        SampleSpec(
            sample_id="d1",
            condition="nontargeting",
            extract_h5="d1.h5",
            metadata={"pileup_path": "d1.bed.gz", "pair_id": "pair-1"},
        ),
        SampleSpec(
            sample_id="t2",
            condition="targeting",
            extract_h5="t2.h5",
            metadata={"pileup_path": "t2.bed.gz", "pair_id": "pair-2"},
        ),
        SampleSpec(
            sample_id="d2",
            condition="nontargeting",
            extract_h5="d2.h5",
            metadata={"pileup_path": "d2.bed.gz", "pair_id": "pair-2"},
        ),
        SampleSpec(
            sample_id="t3",
            condition="targeting",
            extract_h5="t3.h5",
            metadata={"pileup_path": "t3.bed.gz", "pair_id": "pair-3"},
        ),
    ]


def _merge_helper_hits() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "chromosome": "chr2",
                "start": 100,
                "end": 200,
                "strand": "+",
                "window_id": "chr2:100-200",
                "score_value": 0.95,
                "p_value": 0.01,
                "adjusted_p_value": 0.01,
                "rank": 1,
                "modified_count": 19,
                "valid_count": 20,
                "window_fraction": 0.95,
            },
            {
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": "+",
                "window_id": "chr1:0-100",
                "score_value": 0.2,
                "p_value": 0.04,
                "adjusted_p_value": 0.08,
                "rank": 2,
                "modified_count": 2,
                "valid_count": 10,
                "window_fraction": 0.2,
                "numerator_modified_count": 2,
                "numerator_valid_count": 4,
                "denominator_modified_count": 1,
                "denominator_valid_count": 4,
            },
            {
                "chromosome": "chr1",
                "start": 101,
                "end": 200,
                "strand": "+",
                "window_id": "chr1:101-200",
                "score_value": 0.4,
                "p_value": 0.03,
                "adjusted_p_value": 0.06,
                "rank": 3,
                "modified_count": 6,
                "valid_count": 10,
                "window_fraction": 0.6,
                "numerator_modified_count": 4,
                "numerator_valid_count": 6,
                "denominator_modified_count": 2,
                "denominator_valid_count": 6,
            },
            {
                "chromosome": "chr1",
                "start": 500,
                "end": 600,
                "strand": "-",
                "window_id": "chr1:500-600",
                "score_value": 0.8,
                "p_value": 0.05,
                "adjusted_p_value": 0.07,
                "rank": 4,
                "modified_count": 8,
                "valid_count": 10,
                "window_fraction": 0.8,
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
    assert set(result.plot_data) == {"window_score_table", "top_hits_table"}
    assert result.metadata == {
        "analysis_unit": "ensemble_region",
        "representation": "modified_fraction",
        "signal_source": "pileup_counts",
        "score": "effect_size_only",
        "window_size": 1000,
        "step_size": 1000,
        "min_coverage": 5,
        "merge_hits": False,
        "merge_distance": 0,
        "motifs": ["A,0"],
        "include_contigs": None,
        "exclude_contigs": None,
        "contrast_mode": "pairwise",
        "contrast_numerator": ["treated"],
        "contrast_denominator": ["NS"],
    }


def test_scan_genome_matched_pairwise_uses_only_complete_pairs(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_window_summary())

    result = region_discovery.scan_genome(
        samples=_paired_samplespecs(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=500,
        step_size=500,
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["targeting"],
            denominator=["nontargeting"],
            pairing_key="pair_id",
        ),
        score="effect_size_only",
    )

    assert result.metadata["pairing_policy"] == "complete_pairs_only"
    assert result.metadata["n_pairs_used"] == 2
    assert result.metadata["n_pairs_dropped"] == 1
    assert result.windows.loc[0, "modified_count"] == 26
    assert result.windows.loc[0, "valid_count"] == 90
    assert result.windows.loc[0, "window_fraction"] == pytest.approx(26 / 90)
    assert result.windows.loc[0, "score_value"] == pytest.approx(7 / 30)


def test_scan_genome_matched_pairwise_ranks_by_mean_abs_delta(monkeypatch):
    monkeypatch.setattr(
        global_analysis,
        "build_window_summary",
        lambda **_: _mock_paired_pairwise_window_summary(),
    )

    result = region_discovery.scan_genome(
        samples=_paired_samplespecs(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=500,
        step_size=500,
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["targeting"],
            denominator=["nontargeting"],
            pairing_key="pair_id",
        ),
        score="effect_size_only",
    )

    assert list(result.hits["window_id"]) == ["chr1:0-500", "chr1:500-1000"]
    assert result.hits.loc[0, "mean_delta"] == pytest.approx(0.35)
    assert result.hits.loc[0, "mean_abs_delta"] == pytest.approx(0.35)
    assert result.hits.loc[0, "sign_agreement"] == pytest.approx(1.0)
    assert result.metadata["paired_mode"] == "matched_pairwise"
    assert result.metadata["rank_by"] == "mean_abs_delta"


def test_scan_genome_matched_pairwise_errors_on_missing_pairs_in_strict_mode(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_window_summary())

    with pytest.raises(ValueError, match="incomplete matched units"):
        region_discovery.scan_genome(
            samples=_paired_samplespecs(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=500,
            step_size=500,
            contrast=ContrastSpec(
                mode="matched_pairwise",
                numerator=["targeting"],
                denominator=["nontargeting"],
                pairing_key="pair_id",
            ),
            score="effect_size_only",
            pairing_policy="error_on_missing",
        )


def test_build_paired_window_table_collapses_duplicate_rows():
    paired_table, pairing_meta = region_discovery._build_paired_window_table(
        _mock_paired_window_summary(),
        samples=_paired_samplespecs(),
        pairing_key="pair_id",
        required_conditions=["targeting", "nontargeting"],
        pairing_policy="complete_pairs_only",
    )

    collapsed = paired_table.loc[
        (paired_table["window_id"] == "chr1:0-500")
        & (paired_table["pair_id"] == "pair-1")
        & (paired_table["condition"] == "targeting")
    ]

    assert len(collapsed) == 1
    assert collapsed.iloc[0]["modified_count"] == 14
    assert collapsed.iloc[0]["valid_count"] == 30
    assert collapsed.iloc[0]["window_fraction"] == pytest.approx(14 / 30)
    assert pairing_meta == {"n_pairs_used": 2, "n_pairs_dropped": 1}
def test_scan_genome_time_course_errors_on_missing_pairing_key(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_window_summary())

    with pytest.raises(ValueError, match="pairing_key"):
        region_discovery.scan_genome(
            samples=_paired_samplespecs(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=500,
            step_size=500,
            contrast=ContrastSpec(
                mode="time_course",
                time_order=["0min", "15min"],
            ),
            score="effect_size_only",
        )


def test_merge_adjacent_hits_preserves_rank_order_and_merges_counts():
    merged = region_discovery.merge_adjacent_hits(_merge_helper_hits(), merge_distance=1)

    assert list(merged["window_id"]) == ["chr2:100-200", "chr1:0-200", "chr1:500-600"]
    assert list(merged["rank"]) == [1, 2, 4]
    assert list(merged["merged_window_count"]) == [1, 2, 1]
    assert list(merged["score_value"]) == [pytest.approx(0.95), pytest.approx(0.3), pytest.approx(0.8)]
    assert list(merged["p_value"]) == [pytest.approx(0.01), pytest.approx(0.03), pytest.approx(0.05)]
    assert list(merged["adjusted_p_value"]) == [pytest.approx(0.01), pytest.approx(0.06), pytest.approx(0.07)]
    assert list(merged["modified_count"]) == [19, 8, 8]
    assert list(merged["valid_count"]) == [20, 20, 10]
    assert list(merged["window_fraction"]) == [pytest.approx(0.95), pytest.approx(0.4), pytest.approx(0.8)]
    assert pd.isna(merged.loc[0, "numerator_modified_count"])
    assert merged.loc[1, "numerator_modified_count"] == 6
    assert pd.isna(merged.loc[2, "numerator_modified_count"])
    assert pd.isna(merged.loc[0, "numerator_valid_count"])
    assert merged.loc[1, "numerator_valid_count"] == 10
    assert pd.isna(merged.loc[2, "numerator_valid_count"])
    assert pd.isna(merged.loc[0, "denominator_modified_count"])
    assert merged.loc[1, "denominator_modified_count"] == 3
    assert pd.isna(merged.loc[2, "denominator_modified_count"])
    assert pd.isna(merged.loc[0, "denominator_valid_count"])
    assert merged.loc[1, "denominator_valid_count"] == 10
    assert pd.isna(merged.loc[2, "denominator_valid_count"])


def test_hits_to_bed_projects_required_columns_in_order():
    merged = region_discovery.merge_adjacent_hits(_merge_helper_hits(), merge_distance=1)
    bed = region_discovery.hits_to_bed(merged)

    assert list(bed.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    assert pd.api.types.is_integer_dtype(bed["score"])
    assert bed["score"].between(0, 1000).all()
    assert bed.to_dict(orient="records") == [
        {"chrom": "chr2", "start": 100, "end": 200, "name": "chr2:100-200", "score": 950, "strand": "+"},
        {"chrom": "chr1", "start": 0, "end": 200, "name": "chr1:0-200", "score": 300, "strand": "+"},
        {"chrom": "chr1", "start": 500, "end": 600, "name": "chr1:500-600", "score": 800, "strand": "-"},
    ]


def test_discovery_bed_handoff_into_region_contrasts(tmp_path, monkeypatch):
    hits = _merge_helper_hits()
    bed_df = region_discovery.hits_to_bed(hits)
    bed_path = tmp_path / "discovered_hits.bed"
    bed_df.to_csv(bed_path, sep="\t", header=False, index=False)

    counts_by_pileup = {
        "s1.bed.gz": [(1, 10), (1, 10), (8, 10)],
        "s2.bed.gz": [(9, 10), (2, 10), (1, 10)],
    }

    def fake_regions_to_list(function_handle, regions, window_size=None, quiet=True, cores=None, split_large_regions=False):
        pileup_path = function_handle.keywords["bedmethyl_file"]
        regions_dict = region_contrasts.utils.regions_dict_from_input(regions, window_size)
        n_regions = sum(len(region_list) for region_list in regions_dict.values())
        base_counts = counts_by_pileup[pileup_path]
        if len(base_counts) >= n_regions:
            return base_counts[:n_regions]
        repeats = (n_regions // len(base_counts)) + 1
        return (base_counts * repeats)[:n_regions]

    monkeypatch.setattr(region_contrasts.load_processed, "regions_to_list", fake_regions_to_list)

    samples = [
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
    ]
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["NS"],
        reference_condition="NS",
    )

    result = region_contrasts.score_regions(
        samples=samples,
        regions=bed_path,
        motifs=["A,0"],
        contrast=contrast,
        test="effect_size_only",
    )

    assert len(result.regions) == 4
    assert set(result.regions["region_id"]) == {
        "chr2:100-200,+",
        "chr1:0-100,+",
        "chr1:101-200,+",
        "chr1:500-600,-",
    }
    assert list(result.summary["rank"]) == [1, 2, 3, 4]


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
    with pytest.raises(ValueError, match="paired contrast modes"):
        region_discovery.scan_genome(
            samples=[],
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=1000,
            step_size=1000,
            contrast=ContrastSpec(
                mode="single_dataset",
            ),
            score="effect_size_only",
            min_coverage=1,
        )


def test_scan_genome_beta_binomial_requires_explicit_contrast():
    with pytest.raises(ValueError, match="requires an explicit contrast"):
        region_discovery.scan_genome(
            samples=[],
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=1000,
            step_size=1000,
            score="beta_binomial",
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
    assert list(result.hits["rank"]) == [1, 3]
