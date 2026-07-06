from __future__ import annotations

import pandas as pd
import pytest

from dimelo import cluster_reference, motifs


def _sites():
    # cluster C0 sites sit on chr1:0-3000; cluster C1 sites on chr1:10000-13000
    rows = []
    for i in range(10):
        rows.append(("C0", "chr1", i * 300, i * 300 + 100))
    for i in range(10):
        rows.append(("C1", "chr1", 10000 + i * 300, 10000 + i * 300 + 100))
    return pd.DataFrame(rows, columns=["cluster", "chromosome", "start", "end"])


def _ctcf_peaks():
    # peaks overlap all of C0's sites, none of C1's
    return pd.DataFrame(
        {
            "chromosome": ["chr1"] * 10,
            "start": [i * 300 for i in range(10)],
            "end": [i * 300 + 100 for i in range(10)],
        }
    )


def test_site_overlaps_group_flags_overlap():
    with_flag = cluster_reference.site_overlaps_group(_sites(), _ctcf_peaks())
    # all 10 C0 sites overlap, no C1 sites overlap
    c0 = with_flag[with_flag["cluster"] == "C0"]
    c1 = with_flag[with_flag["cluster"] == "C1"]
    assert bool(c0["overlaps"].all())
    assert not bool(c1["overlaps"].any())


def test_site_overlaps_group_validates():
    with pytest.raises(ValueError, match="requires columns"):
        cluster_reference.site_overlaps_group(
            pd.DataFrame({"chromosome": ["chr1"]}), _ctcf_peaks()
        )


def test_site_overlaps_group_half_open_boundaries():
    site = pd.DataFrame({"chromosome": ["chr1"], "start": [100], "end": [200]})

    def _overlaps(ref_start, ref_end):
        ref = pd.DataFrame(
            {"chromosome": ["chr1"], "start": [ref_start], "end": [ref_end]}
        )
        return bool(cluster_reference.site_overlaps_group(site, ref)["overlaps"].iloc[0])

    assert _overlaps(0, 100) is False  # touching at ref.end == site.start -> no overlap
    assert _overlaps(200, 300) is False  # touching at ref.start == site.end -> no overlap
    assert _overlaps(199, 300) is True  # 1 bp overlap
    assert _overlaps(100, 101) is True  # 1 bp overlap at the left edge


def test_cluster_feature_enrichment_finds_associated_cluster():
    with_flag = cluster_reference.site_overlaps_group(_sites(), _ctcf_peaks())
    enrichment = cluster_reference.cluster_feature_enrichment(
        with_flag, feature_column="overlaps"
    ).set_index("cluster")

    # C0 is fully overlapping and enriched vs the rest; C1 is not
    assert enrichment.loc["C0", "feature_fraction"] == pytest.approx(1.0)
    assert enrichment.loc["C1", "feature_fraction"] == pytest.approx(0.0)
    assert enrichment.loc["C0", "odds_ratio"] > 1.0
    assert bool(enrichment.loc["C0", "significant"]) is True
    assert bool(enrichment.loc["C1", "significant"]) is False


def test_cluster_feature_enrichment_ignores_unlabeled_sites():
    # unlabeled (NaN-cluster) sites must not be tested NOR leak into the 2x2 background
    sites = pd.DataFrame(
        {
            "cluster": ["C0", "C0", "C1", "C1", float("nan"), float("nan")],
            "overlaps": [True, True, False, False, True, True],
        }
    )
    enrichment = cluster_reference.cluster_feature_enrichment(
        sites, feature_column="overlaps"
    ).set_index("cluster")
    assert set(enrichment.index) == {"C0", "C1"}
    assert enrichment["n_sites"].sum() == 4  # the 2 unlabeled sites are excluded
    assert enrichment.loc["C0", "feature_fraction"] == pytest.approx(1.0)
    assert enrichment.loc["C1", "feature_fraction"] == pytest.approx(0.0)


def test_annotate_clusters_by_site_groups_labels_association():
    groups = {"CTCF_peaks": _ctcf_peaks()}
    enrichment = cluster_reference.annotate_clusters_by_site_groups(_sites(), groups)
    assert set(enrichment["group"]) == {"CTCF_peaks"}

    summary = cluster_reference.summarize_cluster_associations(enrichment).set_index(
        "cluster"
    )
    # C0's sites are associated with the CTCF peak group; C1's are not
    assert summary.loc["C0", "associated_groups"] == "CTCF_peaks"
    assert summary.loc["C0", "n_significant_groups"] == 1
    assert summary.loc["C1", "associated_groups"] == ""


def test_annotate_clusters_by_motif():
    # C0 sites carry the ACGT consensus; C1 sites do not
    rows = []
    for _ in range(10):
        rows.append(("C0", "AACGTA"))
    for _ in range(10):
        rows.append(("C1", "TTTTTT"))
    sites = pd.DataFrame(rows, columns=["cluster", "sequence"])
    counts = pd.DataFrame(
        {
            "position": [1, 2, 3, 4],
            "A": [40, 0, 0, 0],
            "C": [0, 40, 0, 0],
            "G": [0, 0, 40, 0],
            "T": [0, 0, 0, 40],
        }
    )
    prob = motifs.counts_to_probability(counts, pseudocount=0.1)
    enrichment = cluster_reference.annotate_clusters_by_motif(
        sites, prob, sequence_column="sequence", score_threshold=0.0
    ).set_index("cluster")

    assert enrichment.loc["C0", "feature_fraction"] == pytest.approx(1.0)
    assert enrichment.loc["C1", "feature_fraction"] == pytest.approx(0.0)
    assert bool(enrichment.loc["C0", "significant"]) is True
