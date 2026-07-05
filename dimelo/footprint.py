"""Single-molecule footprint scoring via a 2-state Bernoulli HMM (Q5).

A protein bound to DNA occludes the tethered methyltransferase, so a bound read shows a
central run of *protected* (unmethylated) modifiable positions flanked by *accessible*
(densely methylated) positions. This module segments each read's binary methylation vector
with a lightweight, hand-rolled 2-state Bernoulli HMM (no heavy dependency) and calls a
footprint = a protected run of roughly the expected factor width near the motif anchor.

The HMM has two states — accessible (high methylation probability) and protected (low) —
with Bernoulli emissions over modifiable positions. Parameters are fit by Baum-Welch
(pooled across reads) in log space; per read, Viterbi decodes the state path and
forward-backward gives the per-position protected posterior used to score the footprint.

Footprint presence is a chromatin-state feature (Q5) and a *confidence input* to the Q3
bound/unbound call — it is deliberately NOT part of the Q6 binding-strength metric.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import logsumexp

_EPS = 1e-6

ACCESSIBLE = "accessible"
PROTECTED = "protected"


def _safe_log(x: np.ndarray) -> np.ndarray:
    return np.log(np.clip(np.asarray(x, dtype=float), 1e-300, None))


class BernoulliHMM:
    """2-state Bernoulli HMM over binary methylation vectors (1 = methylated).

    States are ordered so state 0 is the lower-methylation (protected) state after
    ``fit``. ``nan`` observations are treated as missing (uninformative emissions).
    """

    def __init__(
        self,
        *,
        p_emit: tuple[float, float] = (0.05, 0.6),
        stay: float = 0.9,
        start: tuple[float, float] | None = None,
    ) -> None:
        self.p_emit = np.asarray(p_emit, dtype=float)
        self.n_states = len(self.p_emit)
        off = (1.0 - stay) / (self.n_states - 1)
        trans = np.full((self.n_states, self.n_states), off)
        np.fill_diagonal(trans, stay)
        self.trans = trans
        self.start = (
            np.full(self.n_states, 1.0 / self.n_states)
            if start is None
            else np.asarray(start, dtype=float)
        )

    # ------------------------------------------------------------------ #
    def _log_emission(self, obs: np.ndarray) -> np.ndarray:
        """(T,) obs with nan for missing -> (T, n_states) log emission probabilities."""
        p = np.clip(self.p_emit, _EPS, 1 - _EPS)
        log_p = np.log(p)
        log_1mp = np.log(1 - p)
        log_e = np.zeros((len(obs), self.n_states))
        valid = ~np.isnan(obs)
        y = obs[valid][:, None]
        log_e[valid] = y * log_p[None, :] + (1 - y) * log_1mp[None, :]
        return log_e

    def _forward(self, log_e: np.ndarray) -> np.ndarray:
        n = len(log_e)
        log_trans = _safe_log(self.trans)
        log_a = np.full((n, self.n_states), -np.inf)
        log_a[0] = _safe_log(self.start) + log_e[0]
        for t in range(1, n):
            log_a[t] = (
                logsumexp(log_a[t - 1][:, None] + log_trans, axis=0) + log_e[t]
            )
        return log_a

    def _backward(self, log_e: np.ndarray) -> np.ndarray:
        n = len(log_e)
        log_trans = _safe_log(self.trans)
        log_b = np.full((n, self.n_states), -np.inf)
        log_b[-1] = 0.0
        for t in range(n - 2, -1, -1):
            log_b[t] = logsumexp(
                log_trans + (log_e[t + 1] + log_b[t + 1])[None, :], axis=1
            )
        return log_b

    def posterior(self, obs: np.ndarray) -> np.ndarray:
        """Per-position state posteriors ``gamma`` (T, n_states) via forward-backward."""
        log_e = self._log_emission(obs)
        log_a = self._forward(log_e)
        log_b = self._backward(log_e)
        log_gamma = log_a + log_b
        log_gamma -= logsumexp(log_gamma, axis=1, keepdims=True)
        return np.exp(log_gamma)

    def viterbi(self, obs: np.ndarray) -> np.ndarray:
        """Most likely state path (T,) of integer state indices."""
        log_e = self._log_emission(obs)
        n = len(log_e)
        log_trans = _safe_log(self.trans)
        delta = np.full((n, self.n_states), -np.inf)
        back = np.zeros((n, self.n_states), dtype=int)
        delta[0] = _safe_log(self.start) + log_e[0]
        for t in range(1, n):
            scores = delta[t - 1][:, None] + log_trans
            back[t] = np.argmax(scores, axis=0)
            delta[t] = scores[back[t], np.arange(self.n_states)] + log_e[t]
        path = np.zeros(n, dtype=int)
        path[-1] = int(np.argmax(delta[-1]))
        for t in range(n - 2, -1, -1):
            path[t] = back[t + 1, path[t + 1]]
        return path

    def fit(
        self,
        sequences: list[np.ndarray],
        *,
        n_iter: int = 50,
        tol: float = 1e-4,
    ) -> BernoulliHMM:
        """Baum-Welch over pooled read sequences. States are re-ordered so state 0 is the
        protected (lower-methylation) state on return."""
        usable = [np.asarray(s, dtype=float) for s in sequences if len(s) > 0]
        if not usable:
            raise ValueError("BernoulliHMM.fit requires at least one non-empty sequence.")

        prev_ll = -np.inf
        for _ in range(n_iter):
            trans_num = np.zeros((self.n_states, self.n_states))
            start_acc = np.zeros(self.n_states)
            emit_num = np.zeros(self.n_states)
            emit_den = np.zeros(self.n_states)
            gamma_den = np.zeros(self.n_states)
            total_ll = 0.0
            log_trans = _safe_log(self.trans)

            for obs in usable:
                log_e = self._log_emission(obs)
                log_a = self._forward(log_e)
                log_b = self._backward(log_e)
                ll = logsumexp(log_a[-1])
                total_ll += ll
                log_gamma = log_a + log_b - ll
                gamma = np.exp(log_gamma)
                start_acc += gamma[0]
                # transitions
                for t in range(len(obs) - 1):
                    log_xi = (
                        log_a[t][:, None]
                        + log_trans
                        + (log_e[t + 1] + log_b[t + 1])[None, :]
                        - ll
                    )
                    trans_num += np.exp(log_xi)
                    gamma_den += gamma[t]
                # emissions (only observed positions)
                valid = ~np.isnan(obs)
                emit_num += (gamma[valid] * obs[valid][:, None]).sum(axis=0)
                emit_den += gamma[valid].sum(axis=0)

            self.start = start_acc / max(start_acc.sum(), _EPS)
            row_sums = np.clip(gamma_den[:, None], _EPS, None)
            self.trans = trans_num / row_sums
            self.trans /= self.trans.sum(axis=1, keepdims=True)
            self.p_emit = np.clip(emit_num / np.clip(emit_den, _EPS, None), _EPS, 1 - _EPS)

            if abs(total_ll - prev_ll) < tol:
                break
            prev_ll = total_ll

        # order states so state 0 is protected (lower methylation)
        if self.p_emit[0] > self.p_emit[1]:
            order = np.argsort(self.p_emit)
            self.p_emit = self.p_emit[order]
            self.start = self.start[order]
            self.trans = self.trans[np.ix_(order, order)]
        return self


def _protected_runs(path: np.ndarray, protected_state: int) -> list[tuple[int, int]]:
    """Maximal [start, end) index runs where the path is in the protected state."""
    runs: list[tuple[int, int]] = []
    in_run = False
    start = 0
    for i, state in enumerate(path):
        if state == protected_state and not in_run:
            in_run = True
            start = i
        elif state != protected_state and in_run:
            in_run = False
            runs.append((start, i))
    if in_run:
        runs.append((start, len(path)))
    return runs


def call_footprints(
    reads: np.ndarray,
    *,
    anchor_index: int,
    hmm: BernoulliHMM | None = None,
    min_width: int = 3,
    max_width: int | None = None,
    anchor_tolerance: int = 2,
) -> tuple[pd.DataFrame, BernoulliHMM]:
    """Call a per-read footprint from binary methylation vectors.

    ``reads`` is a ``(n_reads, n_positions)`` array of aligned modifiable positions with
    values in ``{0, 1}`` and ``nan`` for missing calls; ``anchor_index`` is the column of
    the motif center. A footprint is a protected run (Viterbi state 0) that overlaps the
    anchor within ``anchor_tolerance`` and whose width is in ``[min_width, max_width]``.

    Returns ``(per_read_table, fitted_hmm)`` where the table has ``read_index,
    footprint_present, protected_start, protected_end, protected_width, footprint_score``
    (``footprint_score`` = mean protected posterior over the called run). If no ``hmm`` is
    given, one is fit by Baum-Welch on ``reads``.
    """
    reads = np.asarray(reads, dtype=float)
    if reads.ndim != 2:
        raise ValueError("reads must be a 2D (n_reads, n_positions) array.")
    n_reads, n_positions = reads.shape
    if not 0 <= anchor_index < n_positions:
        raise ValueError("anchor_index must be a valid column of reads.")
    columns = [
        "read_index",
        "footprint_present",
        "protected_start",
        "protected_end",
        "protected_width",
        "footprint_score",
    ]
    if n_reads == 0:
        return pd.DataFrame(columns=columns), (hmm or BernoulliHMM())

    if hmm is None:
        hmm = BernoulliHMM().fit([reads[i] for i in range(n_reads)])
    protected_state = int(np.argmin(hmm.p_emit))
    width_ceiling = max_width if max_width is not None else n_positions

    rows: list[dict[str, object]] = []
    for i in range(n_reads):
        obs = reads[i]
        path = hmm.viterbi(obs)
        gamma = hmm.posterior(obs)
        runs = _protected_runs(path, protected_state)
        best = None
        for start, end in runs:
            width = end - start
            overlaps_anchor = (
                start - anchor_tolerance <= anchor_index < end + anchor_tolerance
            )
            if overlaps_anchor and min_width <= width <= width_ceiling:
                score = float(gamma[start:end, protected_state].mean())
                if best is None or width > best["protected_width"]:
                    best = {
                        "protected_start": int(start),
                        "protected_end": int(end),
                        "protected_width": int(width),
                        "footprint_score": score,
                    }
        if best is None:
            rows.append(
                {
                    "read_index": i,
                    "footprint_present": False,
                    "protected_start": pd.NA,
                    "protected_end": pd.NA,
                    "protected_width": 0,
                    "footprint_score": float("nan"),
                }
            )
        else:
            rows.append({"read_index": i, "footprint_present": True, **best})
    return pd.DataFrame(rows, columns=columns), hmm


def aggregate_footprint_profile(
    reads: np.ndarray,
    *,
    hmm: BernoulliHMM | None = None,
) -> pd.DataFrame:
    """Per-position aggregate footprint profile across reads.

    Returns ``position, mean_methylation, mean_protected_posterior, n_observed`` — the
    methylation density and mean protected-state posterior at each aligned position (the
    aggregate footprint shows up as a dip in methylation / rise in protected posterior).
    """
    reads = np.asarray(reads, dtype=float)
    if reads.ndim != 2:
        raise ValueError("reads must be a 2D (n_reads, n_positions) array.")
    n_reads, n_positions = reads.shape
    columns = ["position", "mean_methylation", "mean_protected_posterior", "n_observed"]
    if n_reads == 0:
        return pd.DataFrame(columns=columns)

    if hmm is None:
        hmm = BernoulliHMM().fit([reads[i] for i in range(n_reads)])
    protected_state = int(np.argmin(hmm.p_emit))

    protected_posteriors = np.vstack(
        [hmm.posterior(reads[i])[:, protected_state] for i in range(n_reads)]
    )
    with np.errstate(invalid="ignore"):
        mean_methylation = np.nanmean(reads, axis=0)
    n_observed = (~np.isnan(reads)).sum(axis=0)
    return pd.DataFrame(
        {
            "position": np.arange(n_positions),
            "mean_methylation": mean_methylation,
            "mean_protected_posterior": protected_posteriors.mean(axis=0),
            "n_observed": n_observed,
        },
        columns=columns,
    )
