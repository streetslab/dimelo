"""Characterize single-molecule clusters against reference site-groups (Q8/Q10 x clustering).

Answers "what is this cluster's sites associated with?" by testing, for each cluster, whether
its member sites are enriched for a reference feature — overlap with a group of reference
sites (e.g. CTCF ChIP-Atlas/ENCODE peaks, enhancers, a histone-mark domain), or presence of
a sequence motif (JASPAR PWM via ``dimelo.motifs``) — relative to the other clusters.

The common primitive is a per-site binary feature; ``cluster_feature_enrichment`` then runs a
per-cluster 2x2 Fisher enrichment (this cluster vs the rest) with BH across clusters.
``annotate_clusters_by_site_groups`` and ``annotate_clusters_by_motif`` compute the feature
(peak overlap / motif hit) and enrich in one call, so a cluster can be labeled with the
reference groups it is significantly associated with.
"""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
from scipy import stats

from .background import _benjamini_hochberg

_REGION_COLUMNS = ("chromosome", "start", "end")


def site_overlaps_group(
    sites: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    out_column: str = "overlaps",
) -> pd.DataFrame:
    """Add a boolean ``out_column``: does each site overlap any interval in ``reference``?

    Both frames need ``chromosome, start, end`` (half-open). Overlap is
    ``ref.start < site.end and ref.end > site.start`` on the same chromosome.
    """
    for name, frame in (("sites", sites), ("reference", reference)):
        missing = set(_REGION_COLUMNS) - set(frame.columns)
        if missing:
            raise ValueError(f"{name} requires columns: {', '.join(sorted(missing))}.")
    out = sites.copy().reset_index(drop=True)
    reference_by_chrom = {
        str(chrom): group[["start", "end"]].to_numpy(dtype=float)
        for chrom, group in reference.groupby("chromosome")
    }
    overlaps = np.zeros(len(out), dtype=bool)
    for i, (chrom, start, end) in enumerate(
        zip(
            out["chromosome"].astype(str),
            out["start"].astype(float),
            out["end"].astype(float),
            strict=True,
        )
    ):
        intervals = reference_by_chrom.get(chrom)
        if intervals is None or len(intervals) == 0:
            continue
        overlaps[i] = bool(
            np.any((intervals[:, 0] < end) & (intervals[:, 1] > start))
        )
    out[out_column] = overlaps
    return out


def cluster_feature_enrichment(
    sites: pd.DataFrame,
    *,
    feature_column: str,
    cluster_column: str = "cluster",
    fdr: float = 0.05,
) -> pd.DataFrame:
    """Per-cluster enrichment of a binary site feature vs all other clusters.

    For each cluster, a one-sided Fisher's exact test on the 2x2 table
    ``[[in-cluster & feature, in-cluster & not], [other & feature, other & not]]``
    (``alternative='greater'`` — over-representation), with BH across clusters. Returns
    ``cluster, n_sites, n_feature, feature_fraction, odds_ratio, p_value, q_value,
    significant``.
    """
    for column in (feature_column, cluster_column):
        if column not in sites.columns:
            raise ValueError(f"sites requires column: {column}.")
    columns = [
        "cluster",
        "n_sites",
        "n_feature",
        "feature_fraction",
        "odds_ratio",
        "p_value",
        "q_value",
        "significant",
    ]
    # Drop sites with no cluster label so the 2x2 "other clusters" background and the
    # tested clusters share one population (groupby drops NaN groups; the totals must too,
    # else unassigned/noise sites silently inflate every cluster's background row).
    sites = sites[sites[cluster_column].notna()]
    if sites.empty:
        return pd.DataFrame(columns=columns)

    feature = sites[feature_column].astype(bool)
    total = len(sites)
    total_feature = int(feature.sum())

    rows: list[dict[str, object]] = []
    for cluster, group in sites.groupby(cluster_column, sort=False):
        n_sites = len(group)
        n_feature = int(group[feature_column].astype(bool).sum())
        other = total - n_sites
        other_feature = total_feature - n_feature
        table = [
            [n_feature, n_sites - n_feature],
            [other_feature, other - other_feature],
        ]
        odds_ratio, p_value = stats.fisher_exact(table, alternative="greater")
        rows.append(
            {
                "cluster": cluster,
                "n_sites": n_sites,
                "n_feature": n_feature,
                "feature_fraction": n_feature / n_sites if n_sites else float("nan"),
                "odds_ratio": float(odds_ratio),
                "p_value": float(p_value),
            }
        )

    result = pd.DataFrame(rows)
    result["q_value"] = _benjamini_hochberg(result["p_value"])
    result["significant"] = (result["q_value"] <= fdr).fillna(False)
    return result.loc[:, columns]


def annotate_clusters_by_site_groups(
    sites: pd.DataFrame,
    reference_groups: Mapping[str, pd.DataFrame],
    *,
    cluster_column: str = "cluster",
    fdr: float = 0.05,
) -> pd.DataFrame:
    """For each named reference site-group, test per-cluster overlap enrichment.

    ``sites`` needs ``chromosome, start, end`` and ``cluster_column``; ``reference_groups``
    maps a group name (e.g. ``"CTCF_peaks"``) to a ``chromosome/start/end`` frame. Returns
    the per-cluster enrichment table (one block per group) with a ``group`` column, so each
    cluster is characterized by the groups it is significantly associated with.
    """
    if not reference_groups:
        raise ValueError("annotate_clusters_by_site_groups requires at least one group.")
    blocks: list[pd.DataFrame] = []
    for group_name, reference in reference_groups.items():
        with_feature = site_overlaps_group(
            sites, reference, out_column="_overlaps_group"
        )
        enrichment = cluster_feature_enrichment(
            with_feature,
            feature_column="_overlaps_group",
            cluster_column=cluster_column,
            fdr=fdr,
        )
        enrichment.insert(0, "group", group_name)
        blocks.append(enrichment)
    return pd.concat(blocks, ignore_index=True)


def annotate_clusters_by_motif(
    sites: pd.DataFrame,
    probability: pd.DataFrame,
    *,
    sequence_column: str = "sequence",
    cluster_column: str = "cluster",
    score_threshold: float,
    background: dict[str, float] | None = None,
    fdr: float = 0.05,
) -> pd.DataFrame:
    """Test per-cluster enrichment of a sequence motif at the cluster's sites.

    ``sites`` needs ``cluster_column`` and a ``sequence_column`` (the sequence under each
    site). A site "has the motif" when its best PSSM score (``motifs.scan_sequences``) meets
    ``score_threshold``; enrichment is then per-cluster vs the rest. Returns the
    ``cluster_feature_enrichment`` table.
    """
    from . import motifs

    for column in (sequence_column, cluster_column):
        if column not in sites.columns:
            raise ValueError(f"sites requires column: {column}.")
    scanned = motifs.scan_sequences(
        list(sites[sequence_column].astype(str)),
        probability,
        background=background,
    )
    with_feature = sites.copy().reset_index(drop=True)
    with_feature["has_motif"] = (
        scanned["best_score"].to_numpy() >= score_threshold
    )
    return cluster_feature_enrichment(
        with_feature,
        feature_column="has_motif",
        cluster_column=cluster_column,
        fdr=fdr,
    )


def summarize_cluster_associations(
    enrichment: pd.DataFrame,
    *,
    group_column: str = "group",
    cluster_column: str = "cluster",
) -> pd.DataFrame:
    """Collapse a multi-group cluster enrichment table to one row per cluster listing the
    significantly-associated groups (by ascending q-value).

    Returns ``cluster, n_significant_groups, associated_groups`` (comma-joined group names,
    empty when none are significant).
    """
    required = {group_column, cluster_column, "significant", "q_value"}
    missing = required - set(enrichment.columns)
    if missing:
        raise ValueError(
            f"summarize_cluster_associations requires columns: {', '.join(sorted(missing))}."
        )
    rows: list[dict[str, object]] = []
    for cluster, group in enrichment.groupby(cluster_column, sort=False):
        significant = group.loc[group["significant"].astype(bool)].sort_values("q_value")
        rows.append(
            {
                "cluster": cluster,
                "n_significant_groups": int(len(significant)),
                "associated_groups": ", ".join(
                    significant[group_column].astype(str).tolist()
                ),
            }
        )
    return pd.DataFrame(
        rows, columns=["cluster", "n_significant_groups", "associated_groups"]
    )
