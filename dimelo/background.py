"""Negative-control-calibrated single-molecule true-signal calling (Q3 / gap G1).

Given per-read modification counts ``(modified_count k, valid_count n)`` in a window at
each site, for a *target* sample and one or more *negative-control* samples, this module
builds a beta-binomial background null from the control reads and calls each target read
as true signal or background at a controlled FDR.

Null model (per the finalized design):
- **Per-site** beta-binomial fit to the control reads at that site, when the site has at
  least ``min_control_reads`` control reads with coverage.
- **Pooled fallback** beta-binomial fit to *all* control reads, used where a site lacks
  enough control coverage for its own fit. (A context-conditioned pooled fallback keyed on
  coverage / local motif density is a planned refinement; the current fallback is a single
  pooled null.)

Per read, the one-sided upper-tail probability ``P(K >= k | n, null)`` is the p-value that
the read carries *more* modification than background; Benjamini-Hochberg across the tested
target reads yields ``q_value`` and the binary ``is_true_signal`` call at FDR ``alpha``.

The per-read evidence frame is the output of
``region_contrasts.build_single_read_mod_fraction_evidence_table`` (columns
``region_id, condition, read_id, modified_count, valid_count``); ``region_id`` is the site
key and ``condition`` selects target vs control reads.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize
from scipy.special import betaln

_REQUIRED_COLUMNS = ("region_id", "condition", "read_id", "modified_count", "valid_count")


def _beta_binomial_nll(mu: float, s: float, ks: np.ndarray, ns: np.ndarray) -> float:
    """Negative log-likelihood of BetaBinom((mu*s), ((1-mu)*s)) over paired (k, n)."""
    mu = min(max(mu, 1e-6), 1 - 1e-6)
    s = max(s, 1e-6)
    alpha = mu * s
    beta = (1 - mu) * s
    return -float(np.sum(betaln(ks + alpha, ns - ks + beta) - betaln(alpha, beta)))


def fit_beta_binomial(ks, ns) -> tuple[float, float]:
    """MLE ``(alpha, beta)`` of a beta-binomial fit to paired ``(k, n)`` counts.

    Parameterized through ``(mu, s)`` (``alpha = mu*s``, ``beta = (1-mu)*s``) with a
    multi-start over the concentration ``s`` to avoid poor local optima. Reads with
    ``valid_count == 0`` contribute nothing to the likelihood.
    """
    ks = np.asarray(ks, dtype=float)
    ns = np.asarray(ns, dtype=float)
    if ks.shape != ns.shape:
        raise ValueError("fit_beta_binomial requires ks and ns of equal length.")
    if len(ks) == 0 or ns.sum() <= 0:
        raise ValueError(
            "fit_beta_binomial requires at least one read with valid_count > 0."
        )
    mu0 = float(ks.sum() / max(ns.sum(), 1.0))
    mu0 = min(max(mu0, 1e-4), 1 - 1e-4)
    best = None
    for s0 in (5.0, 25.0, 100.0):
        result = minimize(
            lambda p: _beta_binomial_nll(p[0], p[1], ks, ns),
            x0=[mu0, s0],
            method="L-BFGS-B",
            bounds=[(1e-6, 1 - 1e-6), (1e-3, 1e6)],
        )
        if best is None or result.fun < best.fun:
            best = result
    mu, s = float(best.x[0]), float(best.x[1])
    return mu * s, (1 - mu) * s


def _fit_beta_binomial_weighted(
    ks: np.ndarray,
    ns: np.ndarray,
    weights: np.ndarray,
    *,
    mu_lower: float = 1e-6,
) -> tuple[float, float]:
    """Weighted beta-binomial MLE (the M-step of the occupancy mixture). ``mu_lower``
    lower-bounds the component mean so the signal component stays above the background
    (prevents label-switching)."""
    total_weight = float(np.sum(weights))
    if total_weight <= 0 or ns.sum() <= 0:
        # Degenerate weights: fall back to a high-methylation anchor.
        return 0.9, 0.1
    weighted_mean = float(np.sum(weights * ks) / max(np.sum(weights * ns), 1e-9))
    mu0 = min(max(weighted_mean, max(mu_lower, 1e-4)), 1 - 1e-4)

    def nll(p: np.ndarray) -> float:
        mu = min(max(p[0], 1e-6), 1 - 1e-6)
        s = max(p[1], 1e-6)
        alpha = mu * s
        beta = (1 - mu) * s
        per_read = betaln(ks + alpha, ns - ks + beta) - betaln(alpha, beta)
        return -float(np.sum(weights * per_read))

    best = None
    for s0 in (5.0, 25.0, 100.0):
        result = minimize(
            nll,
            x0=[mu0, s0],
            method="L-BFGS-B",
            bounds=[(max(mu_lower, 1e-6), 1 - 1e-6), (1e-3, 1e6)],
        )
        if best is None or result.fun < best.fun:
            best = result
    mu, s = float(best.x[0]), float(best.x[1])
    return mu * s, (1 - mu) * s


def _occupancy_posterior_em(
    ks: np.ndarray,
    ns: np.ndarray,
    bg_alphas: np.ndarray,
    bg_betas: np.ndarray,
    *,
    max_iter: int = 100,
    tol: float = 1e-4,
) -> np.ndarray:
    """Per-read posterior P(signal | k, n) from a 2-component beta-binomial mixture.

    Component 0 (background) is fixed to each read's assigned control null; component 1
    (signal, constrained to a mean above the background) and the mixing weight are fit by
    EM. Reads with ``n == 0`` receive NaN. Returns an array aligned to the inputs.
    """
    posterior = np.full(len(ks), np.nan)
    testable = ns > 0
    if testable.sum() == 0:
        return posterior

    k_t = ks[testable].astype(float)
    n_t = ns[testable].astype(float)
    a0 = bg_alphas[testable].astype(float)
    b0 = bg_betas[testable].astype(float)

    f0 = stats.betabinom.pmf(k_t, n_t, a0, b0)
    background_mean = float(np.mean(a0 / (a0 + b0)))
    # signal component starts above the background mean
    mu1 = min(0.95, max(background_mean + 0.2, 0.5))
    s1 = 20.0
    alpha1, beta1 = mu1 * s1, (1 - mu1) * s1
    pi = 0.5

    try:
        for _ in range(max_iter):
            f1 = stats.betabinom.pmf(k_t, n_t, alpha1, beta1)
            numerator = pi * f1
            responsibilities = numerator / (numerator + (1 - pi) * f0 + 1e-300)
            pi_new = float(np.mean(responsibilities))
            alpha1, beta1 = _fit_beta_binomial_weighted(
                k_t, n_t, responsibilities, mu_lower=background_mean
            )
            converged = abs(pi_new - pi) < tol
            pi = pi_new
            if converged:
                break
        f1 = stats.betabinom.pmf(k_t, n_t, alpha1, beta1)
        numerator = pi * f1
        responsibilities = numerator / (numerator + (1 - pi) * f0 + 1e-300)
    except (ValueError, FloatingPointError):
        return posterior

    posterior[testable] = responsibilities
    return posterior


def _upper_tail_pvalue(k: int, n: int, alpha: float, beta: float) -> float:
    """``P(K >= k | n, BetaBinom(alpha, beta))`` — one-sided test for more modification
    than the control-derived background. Returns 1.0 when ``n <= 0``."""
    if n <= 0:
        return 1.0
    return float(stats.betabinom.sf(k - 1, n, alpha, beta))


def _benjamini_hochberg(p_values, *, testable=None) -> np.ndarray:
    """BH q-values; the multiple-testing total ``m`` counts only testable, finite rows.
    Non-testable rows return NaN (mirrors ``region_contrasts._adjust_p_values_bh``)."""
    p = np.asarray(p_values, dtype=float)
    q = np.full(len(p), np.nan)
    if testable is None:
        mask = np.isfinite(p)
    else:
        mask = np.asarray(testable, dtype=bool) & np.isfinite(p)
    positions = np.where(mask)[0]
    m = len(positions)
    if m == 0:
        return q
    ordered = positions[np.argsort(p[positions])]
    running_min = 1.0
    for rank_from_end, position in enumerate(ordered[::-1], start=1):
        rank = m - rank_from_end + 1
        running_min = min(running_min, float(p[position]) * m / rank)
        q[position] = running_min
    return q


def build_control_null(
    control_reads: pd.DataFrame,
    *,
    min_control_reads: int = 20,
) -> dict[str, object]:
    """Fit the background null from negative-control reads.

    Returns ``{"pooled": (alpha, beta), "per_site": {region_id: (alpha, beta)}}``. A
    per-site null is fit only where the site has ``>= min_control_reads`` control reads
    with positive total coverage; all other sites fall back to the pooled null.
    """
    if control_reads.empty or control_reads["valid_count"].sum() <= 0:
        raise ValueError(
            "build_control_null requires control reads with valid_count > 0."
        )
    pooled = fit_beta_binomial(
        control_reads["modified_count"], control_reads["valid_count"]
    )
    per_site: dict[str, tuple[float, float]] = {}
    for site, rows in control_reads.groupby("region_id", sort=False):
        if len(rows) >= min_control_reads and rows["valid_count"].sum() > 0:
            try:
                per_site[site] = fit_beta_binomial(
                    rows["modified_count"], rows["valid_count"]
                )
            except ValueError:
                continue
    return {"pooled": pooled, "per_site": per_site}


def call_true_signal_reads(
    *,
    evidence: pd.DataFrame,
    target_conditions: Sequence[str],
    control_conditions: Sequence[str],
    fdr: float = 0.05,
    min_control_reads: int = 20,
) -> pd.DataFrame:
    """Call each target read as true signal or background against the control null.

    Parameters
    ----------
    evidence:
        Per-read frame with columns ``region_id, condition, read_id, modified_count,
        valid_count`` (e.g. from
        ``region_contrasts.build_single_read_mod_fraction_evidence_table``).
    target_conditions / control_conditions:
        Condition labels selecting the target reads to call and the negative-control
        reads that define the background null.
    fdr:
        Benjamini-Hochberg FDR threshold for ``is_true_signal``.
    min_control_reads:
        Minimum control reads at a site to fit a per-site null (else pooled fallback).

    Returns
    -------
    Per target read with ``region_id, read_id, condition, modified_count, valid_count,
    read_mod_fraction, background_rate, null_source, p_value, q_value, is_true_signal``.
    Reads with ``valid_count == 0`` are untestable (``p_value``/``q_value`` NaN,
    ``is_true_signal`` False).
    """
    missing = set(_REQUIRED_COLUMNS) - set(evidence.columns)
    if missing:
        raise ValueError(
            "call_true_signal_reads requires columns: "
            f"{', '.join(sorted(missing))}."
        )
    if not 0.0 < fdr < 1.0:
        raise ValueError("fdr must be in the open interval (0, 1).")

    control_conditions = list(control_conditions)
    target_conditions = list(target_conditions)
    control = evidence.loc[evidence["condition"].isin(control_conditions)]
    target = evidence.loc[evidence["condition"].isin(target_conditions)].copy()

    output_columns = [
        "region_id",
        "read_id",
        "condition",
        "modified_count",
        "valid_count",
        "read_mod_fraction",
        "background_rate",
        "null_source",
        "p_value",
        "q_value",
        "is_true_signal",
        "occupancy_posterior",
    ]
    if control.empty:
        raise ValueError(
            "call_true_signal_reads found no control reads for the requested "
            f"control_conditions: {control_conditions}."
        )
    if target.empty:
        return pd.DataFrame(columns=output_columns)

    null = build_control_null(control, min_control_reads=min_control_reads)
    pooled_alpha, pooled_beta = null["pooled"]
    per_site = null["per_site"]

    p_values: list[float] = []
    background_rates: list[float] = []
    null_sources: list[str] = []
    alphas: list[float] = []
    betas: list[float] = []
    for site, modified_count, valid_count in zip(
        target["region_id"],
        target["modified_count"].astype(int),
        target["valid_count"].astype(int),
        strict=True,
    ):
        if site in per_site:
            alpha, beta = per_site[site]
            source = "per_site"
        else:
            alpha, beta = pooled_alpha, pooled_beta
            source = "pooled"
        alphas.append(alpha)
        betas.append(beta)
        background_rates.append(alpha / (alpha + beta))
        null_sources.append(source)
        if valid_count <= 0:
            p_values.append(float("nan"))
        else:
            p_values.append(_upper_tail_pvalue(modified_count, valid_count, alpha, beta))

    target["background_rate"] = background_rates
    target["null_source"] = null_sources
    target["p_value"] = p_values
    target["read_mod_fraction"] = (
        target["modified_count"]
        .astype(float)
        .div(target["valid_count"].where(target["valid_count"] != 0))
        .fillna(0.0)
    )
    testable = target["valid_count"] > 0
    target["q_value"] = _benjamini_hochberg(target["p_value"], testable=testable)
    target["is_true_signal"] = (target["q_value"] <= fdr).fillna(False)
    target["occupancy_posterior"] = _occupancy_posterior_em(
        target["modified_count"].to_numpy(dtype=float),
        target["valid_count"].to_numpy(dtype=float),
        np.asarray(alphas, dtype=float),
        np.asarray(betas, dtype=float),
    )

    return target.loc[:, output_columns].reset_index(drop=True)


def summarize_site_occupancy(called_reads: pd.DataFrame) -> pd.DataFrame:
    """Per-site control-calibrated occupancy rate = fraction of true-signal reads.

    Returns ``region_id, n_reads, n_true_signal, occupancy_rate`` (feeds Q6 binding
    strength).
    """
    required = {"region_id", "read_id", "is_true_signal"}
    missing = required - set(called_reads.columns)
    if missing:
        raise ValueError(
            f"summarize_site_occupancy requires columns: {', '.join(sorted(missing))}."
        )
    if called_reads.empty:
        return pd.DataFrame(
            columns=["region_id", "n_reads", "n_true_signal", "occupancy_rate"]
        )
    grouped = called_reads.groupby("region_id", sort=False, as_index=False).agg(
        n_reads=("read_id", "size"),
        n_true_signal=("is_true_signal", "sum"),
    )
    grouped["n_true_signal"] = grouped["n_true_signal"].astype(int)
    grouped["occupancy_rate"] = grouped["n_true_signal"] / grouped["n_reads"].where(
        grouped["n_reads"] != 0
    )
    return grouped


def background_removed_pileup(
    called_reads: pd.DataFrame,
    *,
    weighting: str = "hard",
) -> pd.DataFrame:
    """Per-site background-removed pileup from called reads.

    ``weighting="hard"`` sums counts over ``is_true_signal`` reads only.
    ``weighting="posterior"`` weights each read's counts by ``occupancy_posterior``
    (reads with NaN posterior contribute 0). Returns ``region_id, modified_count,
    valid_count, mod_fraction`` — the background-removed signal track (links to the Q2
    ensemble ``background_adjusted`` correction).
    """
    if weighting not in {"hard", "posterior"}:
        raise ValueError("weighting must be 'hard' or 'posterior'.")
    required = {"region_id", "modified_count", "valid_count"}
    if weighting == "hard":
        required.add("is_true_signal")
    else:
        required.add("occupancy_posterior")
    missing = required - set(called_reads.columns)
    if missing:
        raise ValueError(
            f"background_removed_pileup requires columns: {', '.join(sorted(missing))}."
        )

    columns = ["region_id", "modified_count", "valid_count", "mod_fraction"]
    if called_reads.empty:
        return pd.DataFrame(columns=columns)

    frame = called_reads.copy()
    if weighting == "hard":
        weights = frame["is_true_signal"].astype(bool).astype(float)
    else:
        weights = pd.to_numeric(frame["occupancy_posterior"], errors="coerce").fillna(
            0.0
        )
    frame["_w_modified"] = weights * frame["modified_count"].astype(float)
    frame["_w_valid"] = weights * frame["valid_count"].astype(float)

    pooled = frame.groupby("region_id", sort=False, as_index=False).agg(
        modified_count=("_w_modified", "sum"),
        valid_count=("_w_valid", "sum"),
    )
    pooled["mod_fraction"] = pooled["modified_count"] / pooled["valid_count"].where(
        pooled["valid_count"] != 0
    )
    pooled["mod_fraction"] = pooled["mod_fraction"].fillna(0.0)
    return pooled.loc[:, columns]
