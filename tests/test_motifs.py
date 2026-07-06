from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from dimelo import jaspar, motifs


def _ctcf_like_counts():
    # a short, strongly-informative synthetic motif: positions favor A, C, G, T in turn
    return pd.DataFrame(
        {
            "position": [1, 2, 3, 4],
            "A": [40, 0, 0, 0],
            "C": [0, 40, 0, 0],
            "G": [0, 0, 40, 0],
            "T": [0, 0, 0, 40],
        }
    )


# --------------------------------------------------------------------------- #
# JASPAR client (offline)                                                     #
# --------------------------------------------------------------------------- #


def test_jaspar_matrix_from_json_parses_pfm():
    payload = {
        "matrix_id": "MA0139.1",
        "name": "CTCF",
        "pfm": {"A": [1, 2], "C": [3, 4], "G": [5, 6], "T": [7, 8]},
    }
    frame = jaspar.jaspar_matrix_from_json(payload)
    assert list(frame["position"]) == [1, 2]
    assert list(frame["A"]) == [1.0, 2.0]
    assert frame.attrs["name"] == "CTCF"


def test_jaspar_matrix_from_json_validates():
    with pytest.raises(ValueError, match="pfm"):
        jaspar.jaspar_matrix_from_json({"matrix_id": "x"})
    with pytest.raises(ValueError, match="missing base"):
        jaspar.jaspar_matrix_from_json({"pfm": {"A": [1], "C": [1], "G": [1]}})


def test_fetch_jaspar_matrix_uses_cache_without_network(tmp_path):
    # pre-populate the cache so no network request is made
    payload = {
        "matrix_id": "MA0139.1",
        "name": "CTCF",
        "pfm": {"A": [1], "C": [2], "G": [3], "T": [4]},
    }
    (tmp_path / "MA0139.1.json").write_text(json.dumps(payload))
    frame = jaspar.fetch_jaspar_matrix("MA0139.1", cache_dir=tmp_path)
    assert frame.attrs["matrix_id"] == "MA0139.1"
    assert list(frame["T"]) == [4.0]


# --------------------------------------------------------------------------- #
# motif analysis (offline)                                                    #
# --------------------------------------------------------------------------- #


def test_counts_to_probability_rows_sum_to_one():
    prob = motifs.counts_to_probability(_ctcf_like_counts(), pseudocount=0.8)
    row_sums = prob[["A", "C", "G", "T"]].sum(axis=1)
    assert np.allclose(row_sums, 1.0)
    # the favored base still dominates its position
    assert prob.loc[0, "A"] > 0.9


def test_information_content_peaks_at_conserved_positions():
    prob = motifs.counts_to_probability(_ctcf_like_counts(), pseudocount=0.1)
    ic = motifs.information_content(prob)
    per_position = ic.groupby("position")["bits"].sum()
    # near-deterministic positions carry close to 2 bits
    assert (per_position > 1.5).all()
    assert per_position.max() <= 2.0 + 1e-9


def test_scan_sequences_scores_matching_window_highest():
    prob = motifs.counts_to_probability(_ctcf_like_counts(), pseudocount=0.1)
    # the motif consensus is ACGT; embed it with flanks
    result = motifs.scan_sequences(["TTACGTTT", "TTTTTTTT"], prob).set_index(
        "sequence_index"
    )
    assert result.loc[0, "best_offset"] == 2  # 'ACGT' starts at index 2
    assert result.loc[0, "best_score"] > result.loc[1, "best_score"]


def test_observed_pfm_and_compare_matches_reference():
    # sequences that all match the ACGT consensus -> observed motif == reference motif
    sequences = ["ACGT"] * 20
    observed = motifs.observed_pfm(sequences)
    assert list(observed["A"]) == [20.0, 0.0, 0.0, 0.0]

    observed_prob = motifs.counts_to_probability(observed, pseudocount=0.1)
    reference_prob = motifs.counts_to_probability(_ctcf_like_counts(), pseudocount=0.1)
    comparison = motifs.compare_motifs(observed_prob, reference_prob)
    assert comparison["pearson_r"] > 0.99
    assert comparison["mean_per_position_r"] > 0.99


def test_observed_pfm_requires_equal_length():
    with pytest.raises(ValueError, match="equal-length"):
        motifs.observed_pfm(["ACGT", "AC"])


def test_compare_motifs_requires_equal_shape():
    a = motifs.counts_to_probability(_ctcf_like_counts())
    b = a.iloc[:2].copy()
    with pytest.raises(ValueError, match="equal-length"):
        motifs.compare_motifs(a, b)


def test_extract_site_sequences_from_fasta(tmp_path):
    import pysam

    fasta_path = tmp_path / "mini.fa"
    fasta_path.write_text(">chr1\n" + "ACGT" * 25 + "\n")  # 100 bp
    pysam.faidx(str(fasta_path))
    sites = pd.DataFrame({"chromosome": ["chr1"], "center": [50]})
    seqs = motifs.extract_site_sequences(fasta_path, sites, width=8)
    assert len(seqs) == 1
    assert len(seqs[0]) == 8
    assert set(seqs[0]) <= set("ACGT")


def test_sequence_logo_plotting():
    import matplotlib

    matplotlib.use("Agg")
    from dimelo import plotting, plotting_matplotlib

    prob = motifs.counts_to_probability(_ctcf_like_counts(), pseudocount=0.1)
    ic = motifs.information_content(prob)
    payload = plotting.prepare_sequence_logo_data(information_content=ic)
    assert payload["metadata"]["positions"] == [1, 2, 3, 4]
    fig, _ = plotting_matplotlib.plot_sequence_logo_matplotlib(payload, title="logo")
    assert fig is not None
