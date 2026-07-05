from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dimelo import utils
from dimelo.utils import ParsedMotif


def test_cores_to_run_clamps_and_defaults():
    avail = utils._effective_cpu_count()
    assert utils.cores_to_run(None) == avail
    assert utils.cores_to_run(1) == 1
    # request more than available -> clamped to available
    assert utils.cores_to_run(avail + 1000) == avail
    # request <= 0 -> floored to 1
    assert utils.cores_to_run(0) == 1
    assert utils.cores_to_run(-5) == 1


def test_parsed_motif_two_and_three_parts():
    two = ParsedMotif("A,0")
    assert two.motif_seq == "A"
    assert two.modified_pos == 0
    assert two.modified_base == "A"
    assert two.mod_codes == {"a", "Y"}

    cpg = ParsedMotif("CG,0")
    assert cpg.modified_base == "C"
    assert cpg.mod_codes == {"m", "Z"}

    three = ParsedMotif("A,0,a")
    assert three.mod_codes == {"a"}


def test_parsed_motif_invalid_inputs():
    with pytest.raises(ValueError, match="out-of-range"):
        ParsedMotif("A,2")
    with pytest.raises(ValueError, match="2 or 3 comma-separated"):
        ParsedMotif("A")


def test_parse_region_string_basic_strand_and_window():
    chrom, (start, end, strand) = utils.parse_region_string("chr1:100-200", None)
    assert chrom == "chr1"
    assert (start, end, strand) == (100, 200, ".")

    _, (_, _, strand2) = utils.parse_region_string("chr1:100-200,+", None)
    assert strand2 == "+"

    # window_size recenters around the midpoint (150 +/- 25)
    _, (ws_start, ws_end, _) = utils.parse_region_string("chr1:100-200", 25)
    assert (ws_start, ws_end) == (125, 175)


def test_parse_region_string_invalid():
    with pytest.raises(ValueError, match="Invalid region string"):
        utils.parse_region_string("not-a-region", None)


def test_sanitize_path_args_and_check_len_equal():
    a, b, c = utils.sanitize_path_args("x/y", None, Path("z"))
    assert a == Path("x/y")
    assert b is None
    assert c == Path("z")

    assert utils.check_len_equal([1, 2], [3, 4], [5, 6]) is True
    assert utils.check_len_equal([1, 2], [3]) is False


def test_adjust_threshold_scales_255_and_passes_fraction():
    # values in 0-1 pass through
    assert utils.adjust_threshold(0.7) == 0.7
    # values > 1 assumed 0-255 and scaled
    assert utils.adjust_threshold(255) == pytest.approx(1.0)
    assert utils.adjust_threshold(128) == pytest.approx(128 / 255)
    # zero/negative pass through unchanged
    assert utils.adjust_threshold(0) == 0


def test_regions_dict_from_input_region_strings_sorted():
    regions_dict = utils.regions_dict_from_input(
        ["chr1:300-400", "chr1:100-200,+", "chr2:5-9"]
    )
    assert set(regions_dict) == {"chr1", "chr2"}
    # sorted by start within a contig
    assert regions_dict["chr1"] == [(100, 200, "+"), (300, 400, ".")]
    assert regions_dict["chr2"] == [(5, 9, ".")]


def test_regions_dict_from_input_rejects_nonpositive_window():
    with pytest.raises(ValueError, match="Invalid window_size"):
        utils.regions_dict_from_input(["chr1:100-200"], window_size=0)


def test_add_region_to_dict_reads_bed_with_and_without_strand(tmp_path):
    bed = tmp_path / "regions.bed"
    bed.write_text("chr1\t10\t20\nchr2\t30\t40\tname\t0\t-\n")
    regions_dict = utils.regions_dict_from_input(bed)
    assert regions_dict["chr1"] == [(10, 20, ".")]
    assert regions_dict["chr2"] == [(30, 40, "-")]


def test_add_region_to_dict_rejects_non_bed_path(tmp_path):
    bad = tmp_path / "regions.txt"
    bad.write_text("chr1\t10\t20\n")
    with pytest.raises(ValueError, match="not pointing to a .bed file"):
        utils.regions_dict_from_input(bad)


def test_smooth_rolling_mean_centered():
    smoothed = utils.smooth_rolling_mean(np.array([0.0, 3.0, 6.0]), window=3)
    # centered window, min_periods=1: edges average 2 values, middle averages 3
    assert smoothed[0] == pytest.approx(1.5)
    assert smoothed[1] == pytest.approx(3.0)
    assert smoothed[2] == pytest.approx(4.5)


def test_random_sample_is_reproducible_with_seed():
    array = np.arange(100)
    first = utils.random_sample(array, n=10, seed=1234)
    second = utils.random_sample(array, n=10, seed=1234)
    assert np.array_equal(first, second)
    assert len(first) == 10
    assert set(first.tolist()).issubset(set(array.tolist()))
