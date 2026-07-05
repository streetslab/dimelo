"""Joint occupancy on single long reads (Q7).

For a user-supplied set of anchor pairs (bedpe / ChIA-PET-style, or promoter-enhancer
pairs), this quantifies whether the two anchors are co-occupied on the *same* single
molecule more (or less) often than expected if binding were independent. Only reads that
span BOTH anchors contribute.

Per spanning read, the joint state (both / a_only / b_only / neither) comes from the Q3
single-molecule ``is_true_signal`` call at each anchor (``background.call_true_signal_reads``).
Across the spanning reads of a pair: a 2x2 contingency table -> Fisher's exact independence
test, an odds ratio, and ``log2`` observed/expected co-occupancy; plus a Spearman
correlation of the soft ``occupancy_posterior`` between the anchors (continuous coupling).
Benjamini-Hochberg controls the FDR across pairs.

Positive coupling — especially at short 1D distance — may reflect a trans / long-range
contact artifact (Q4) rather than independent binding; ``pair_distance_bp`` is reported so
that confound can be assessed (and corrected once Q4 lands).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .background import _benjamini_hochberg

_READ_COLUMNS = ("read_id", "region_id", "is_true_signal")


def _log2_observed_over_expected(p_both: float, p_a: float, p_b: float) -> float:
    """log2 of observed co-occupancy over the independence expectation, with a small
    pseudocount so an all-zero cell yields a finite (very negative) value."""
    eps = 1e-9
    expected = p_a * p_b
    return float(np.log2((p_both + eps) / (expected + eps)))


def joint_occupancy(
    *,
    called_reads: pd.DataFrame,
    pairs: pd.DataFrame,
    min_spanning_reads: int = 5,
    multiple_testing: str = "fdr_bh",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Quantify single-molecule joint occupancy for anchor pairs.

    Parameters
    ----------
    called_reads:
        Per-read Q3 calls with columns ``read_id, region_id, is_true_signal`` (and
        optionally ``occupancy_posterior`` for the soft correlation), e.g. from
        ``background.call_true_signal_reads``. ``region_id`` is the anchor id; a read that
        spans two anchors appears once per anchor with a shared ``read_id``.
    pairs:
        Anchor pairs with columns ``site_a, site_b`` (matching ``region_id`` values) and
        optionally ``pair_id`` and ``distance_bp``.
    min_spanning_reads:
        Pairs with fewer spanning reads are reported but not tested (``p_value`` NaN,
        excluded from the BH total).

    Returns
    -------
    ``(pair_summary, per_read_states)``:
    - ``pair_summary``: one row per pair with ``pair_id, site_a, site_b,
      n_reads_spanning_both, n_both, n_a_only, n_b_only, n_neither, frac_both, odds_ratio,
      log2_obs_exp, p_value, q_value, significant, posterior_corr, posterior_corr_p,
      distance_bp``.
    - ``per_read_states``: ``pair_id, site_a, site_b, read_id, joint_state, is_signal_a,
      is_signal_b, posterior_a, posterior_b``.
    """
    if multiple_testing != "fdr_bh":
        raise ValueError("joint_occupancy currently requires multiple_testing='fdr_bh'.")
    missing = set(_READ_COLUMNS) - set(called_reads.columns)
    if missing:
        raise ValueError(
            f"joint_occupancy called_reads requires columns: {', '.join(sorted(missing))}."
        )
    if not {"site_a", "site_b"}.issubset(pairs.columns):
        raise ValueError("joint_occupancy pairs requires columns: site_a, site_b.")

    has_posterior = "occupancy_posterior" in called_reads.columns
    read_cols = ["read_id", "is_true_signal"] + (
        ["occupancy_posterior"] if has_posterior else []
    )

    pair_columns = [
        "pair_id",
        "site_a",
        "site_b",
        "n_reads_spanning_both",
        "n_both",
        "n_a_only",
        "n_b_only",
        "n_neither",
        "frac_both",
        "odds_ratio",
        "log2_obs_exp",
        "p_value",
        "q_value",
        "significant",
        "posterior_corr",
        "posterior_corr_p",
        "distance_bp",
    ]
    read_state_columns = [
        "pair_id",
        "site_a",
        "site_b",
        "read_id",
        "joint_state",
        "is_signal_a",
        "is_signal_b",
        "posterior_a",
        "posterior_b",
    ]
    if pairs.empty:
        return (
            pd.DataFrame(columns=pair_columns),
            pd.DataFrame(columns=read_state_columns),
        )

    pair_rows: list[dict[str, object]] = []
    read_rows: list[dict[str, object]] = []
    for _, pair in pairs.reset_index(drop=True).iterrows():
        site_a = pair["site_a"]
        site_b = pair["site_b"]
        pair_id = pair["pair_id"] if "pair_id" in pairs.columns else f"{site_a}|{site_b}"
        distance_bp = pair["distance_bp"] if "distance_bp" in pairs.columns else pd.NA

        side_a = called_reads.loc[called_reads["region_id"] == site_a, read_cols]
        side_b = called_reads.loc[called_reads["region_id"] == site_b, read_cols]
        spanning = side_a.merge(side_b, on="read_id", suffixes=("_a", "_b"))
        n = len(spanning)

        is_a = spanning["is_true_signal_a"].astype(bool) if n else pd.Series(dtype=bool)
        is_b = spanning["is_true_signal_b"].astype(bool) if n else pd.Series(dtype=bool)
        n_both = int((is_a & is_b).sum())
        n_a_only = int((is_a & ~is_b).sum())
        n_b_only = int((~is_a & is_b).sum())
        n_neither = int((~is_a & ~is_b).sum())

        # Emit per-read joint states by iterating the merged frame once (no per-read
        # re-lookup, so duplicate read_ids are handled row-by-row).
        if n:
            posterior_a_values = (
                spanning["occupancy_posterior_a"].to_numpy(dtype=float)
                if has_posterior
                else np.full(n, np.nan)
            )
            posterior_b_values = (
                spanning["occupancy_posterior_b"].to_numpy(dtype=float)
                if has_posterior
                else np.full(n, np.nan)
            )
            for read_id, sa, sb, post_a, post_b in zip(
                spanning["read_id"],
                is_a,
                is_b,
                posterior_a_values,
                posterior_b_values,
                strict=True,
            ):
                joint_state = (
                    "both"
                    if sa and sb
                    else "a_only"
                    if sa
                    else "b_only"
                    if sb
                    else "neither"
                )
                read_rows.append(
                    {
                        "pair_id": pair_id,
                        "site_a": site_a,
                        "site_b": site_b,
                        "read_id": read_id,
                        "joint_state": joint_state,
                        "is_signal_a": bool(sa),
                        "is_signal_b": bool(sb),
                        "posterior_a": float(post_a),
                        "posterior_b": float(post_b),
                    }
                )

        testable = n >= min_spanning_reads
        if n > 0:
            frac_both = n_both / n
            p_a = (n_both + n_a_only) / n
            p_b = (n_both + n_b_only) / n
            log2_obs_exp = _log2_observed_over_expected(frac_both, p_a, p_b)
        else:
            frac_both = float("nan")
            log2_obs_exp = float("nan")

        if testable:
            odds_ratio, p_value = stats.fisher_exact(
                [[n_both, n_a_only], [n_b_only, n_neither]]
            )
            odds_ratio = float(odds_ratio)
            p_value = float(p_value)
        else:
            odds_ratio = float("nan")
            p_value = float("nan")

        posterior_corr = float("nan")
        posterior_corr_p = float("nan")
        if has_posterior and n >= 2:
            post_a = spanning["occupancy_posterior_a"].to_numpy(dtype=float)
            post_b = spanning["occupancy_posterior_b"].to_numpy(dtype=float)
            if np.nanstd(post_a) > 0 and np.nanstd(post_b) > 0:
                rho, corr_p = stats.spearmanr(post_a, post_b)
                posterior_corr = float(rho)
                posterior_corr_p = float(corr_p)

        pair_rows.append(
            {
                "pair_id": pair_id,
                "site_a": site_a,
                "site_b": site_b,
                "n_reads_spanning_both": n,
                "n_both": n_both,
                "n_a_only": n_a_only,
                "n_b_only": n_b_only,
                "n_neither": n_neither,
                "frac_both": frac_both,
                "odds_ratio": odds_ratio,
                "log2_obs_exp": log2_obs_exp,
                "p_value": p_value,
                "distance_bp": distance_bp,
                "posterior_corr": posterior_corr,
                "posterior_corr_p": posterior_corr_p,
                "_testable": testable,
            }
        )

    pair_summary = pd.DataFrame(pair_rows)
    pair_summary["q_value"] = _benjamini_hochberg(
        pair_summary["p_value"], testable=pair_summary["_testable"]
    )
    pair_summary["significant"] = (pair_summary["q_value"] <= 0.05).fillna(False)
    pair_summary = pair_summary.drop(columns="_testable").loc[:, pair_columns]

    per_read_states = pd.DataFrame(read_rows, columns=read_state_columns)
    return pair_summary, per_read_states


def summarize_joint_states(per_read_states: pd.DataFrame) -> pd.DataFrame:
    """Per-pair joint-state counts and fractions (for a stacked-bar view)."""
    required = {"pair_id", "joint_state"}
    missing = required - set(per_read_states.columns)
    if missing:
        raise ValueError(
            f"summarize_joint_states requires columns: {', '.join(sorted(missing))}."
        )
    states = ["both", "a_only", "b_only", "neither"]
    columns = ["pair_id", *states, "total", *[f"frac_{s}" for s in states]]
    if per_read_states.empty:
        return pd.DataFrame(columns=columns)
    counts = (
        per_read_states.groupby(["pair_id", "joint_state"]).size().unstack(fill_value=0)
    )
    for state in states:
        if state not in counts.columns:
            counts[state] = 0
    counts = counts.loc[:, states]
    counts["total"] = counts.sum(axis=1)
    for state in states:
        counts[f"frac_{state}"] = counts[state] / counts["total"].where(
            counts["total"] != 0
        )
    return counts.reset_index().loc[:, columns]


def anchor_pairs_from_sequence(
    pairs: Sequence[tuple[str, str]],
    *,
    distances_bp: Sequence[int] | None = None,
) -> pd.DataFrame:
    """Build a ``pairs`` DataFrame from a sequence of ``(site_a, site_b)`` tuples."""
    frame = pd.DataFrame(list(pairs), columns=["site_a", "site_b"])
    frame.insert(0, "pair_id", [f"{a}|{b}" for a, b in zip(frame["site_a"], frame["site_b"], strict=True)])
    if distances_bp is not None:
        frame["distance_bp"] = list(distances_bp)
    return frame
