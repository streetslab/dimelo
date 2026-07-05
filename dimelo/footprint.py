"""Single-molecule footprint scoring via a 2-state Bernoulli HMM (Q5).

A protein bound to DNA occludes the tethered methyltransferase, so a bound read shows a
central run of *protected* (unmethylated) modifiable positions flanked by *accessible*
(densely methylated) positions. This module segments each read's binary methylation vector
with a 2-state HMM (accessible / protected, Bernoulli emissions over modifiable positions)
and calls a footprint = a protected run of roughly the expected factor width near the motif
anchor.

The HMM is backed by ``hmmlearn``'s ``CategoricalHMM`` (two symbols: 0 = unmethylated,
1 = methylated). Parameters are fit by Baum-Welch (pooled across reads); per read, Viterbi
decodes the accessible/protected path and the forward-backward posterior gives the
per-position protected probability used to score the footprint. ``nan`` observations are
treated as missing (dropped for fitting/decoding, then carried forward so the decoded path
stays aligned to the input columns).

Footprint presence is a chromatin-state feature (Q5) and a *confidence input* to the Q3
bound/unbound call — it is deliberately NOT part of the Q6 binding-strength metric.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from hmmlearn.hmm import CategoricalHMM

_EPS = 1e-6

ACCESSIBLE = "accessible"
PROTECTED = "protected"


class BernoulliHMM:
    """2-state Bernoulli HMM over binary methylation vectors, backed by
    ``hmmlearn.CategoricalHMM``. State 0 is the protected (lower-methylation) state after
    ``fit``. ``nan`` observations are missing (dropped for fit/decode, carried forward)."""

    def __init__(
        self,
        *,
        p_emit: tuple[float, float] = (0.05, 0.6),
        stay: float = 0.9,
        start: tuple[float, float] | None = None,
        n_iter: int = 50,
        tol: float = 1e-4,
    ) -> None:
        self.n_states = len(p_emit)
        self.n_iter = int(n_iter)
        self.tol = float(tol)
        self._model = self._build_model(
            np.asarray(p_emit, dtype=float),
            self._default_trans(stay),
            self._default_start(start),
        )

    # ------------------------------------------------------------------ #
    def _default_trans(self, stay: float) -> np.ndarray:
        off = (1.0 - stay) / (self.n_states - 1)
        trans = np.full((self.n_states, self.n_states), off)
        np.fill_diagonal(trans, stay)
        return trans

    def _default_start(self, start: tuple[float, float] | None) -> np.ndarray:
        if start is None:
            return np.full(self.n_states, 1.0 / self.n_states)
        return np.asarray(start, dtype=float)

    def _build_model(
        self, p_emit: np.ndarray, trans: np.ndarray, start: np.ndarray
    ) -> CategoricalHMM:
        model = CategoricalHMM(
            n_components=self.n_states,
            n_features=2,
            n_iter=self.n_iter,
            tol=self.tol,
            random_state=0,
        )
        model.startprob_ = start / start.sum()
        model.transmat_ = trans / trans.sum(axis=1, keepdims=True)
        p = np.clip(p_emit, _EPS, 1 - _EPS)
        # symbol 0 = unmethylated, symbol 1 = methylated
        model.emissionprob_ = np.column_stack([1 - p, p])
        return model

    @property
    def p_emit(self) -> np.ndarray:
        return self._model.emissionprob_[:, 1]

    @property
    def trans(self) -> np.ndarray:
        return self._model.transmat_

    @property
    def start(self) -> np.ndarray:
        return self._model.startprob_

    # ------------------------------------------------------------------ #
    def fit(
        self,
        sequences: list[np.ndarray],
        *,
        n_iter: int | None = None,
        tol: float | None = None,
    ) -> BernoulliHMM:
        """Baum-Welch over pooled read sequences (nan dropped). States are re-ordered so
        state 0 is protected (lower methylation) on return."""
        symbol_blocks: list[np.ndarray] = []
        lengths: list[int] = []
        for sequence in sequences:
            values = np.asarray(sequence, dtype=float)
            observed = values[~np.isnan(values)]
            if len(observed) == 0:
                continue
            symbol_blocks.append(observed.astype(int).reshape(-1, 1))
            lengths.append(len(observed))
        if not symbol_blocks:
            raise ValueError("BernoulliHMM.fit requires at least one non-empty sequence.")

        # Informed initialization (init_params="" -> keep these as the EM starting point):
        # spread the states' methylation emission low->high so Baum-Welch separates a
        # protected (low) from an accessible (high) state instead of collapsing from a
        # poor random init.
        model = CategoricalHMM(
            n_components=self.n_states,
            n_features=2,
            n_iter=self.n_iter if n_iter is None else int(n_iter),
            tol=self.tol if tol is None else float(tol),
            random_state=0,
            init_params="",
            params="ste",
        )
        model.startprob_ = np.full(self.n_states, 1.0 / self.n_states)
        model.transmat_ = self._default_trans(0.9)
        init_p = np.clip(np.linspace(0.1, 0.9, self.n_states), _EPS, 1 - _EPS)
        model.emissionprob_ = np.column_stack([1 - init_p, init_p])
        model.fit(np.vstack(symbol_blocks), lengths)
        # order states so state 0 is protected (lowest methylation emission)
        order = np.argsort(model.emissionprob_[:, 1])
        model.startprob_ = model.startprob_[order]
        model.transmat_ = model.transmat_[np.ix_(order, order)]
        model.emissionprob_ = model.emissionprob_[order]
        self._model = model
        return self

    # ------------------------------------------------------------------ #
    def _decode(self, obs: np.ndarray, *, proba: bool):
        obs = np.asarray(obs, dtype=float)
        n = len(obs)
        valid = ~np.isnan(obs)
        if valid.sum() == 0:
            if proba:
                return np.full((n, self.n_states), 1.0 / self.n_states)
            return np.zeros(n, dtype=int)
        symbols = obs[valid].astype(int).reshape(-1, 1)
        if proba:
            values = self._model.predict_proba(symbols)
            full = np.full((n, self.n_states), np.nan)
            full[valid] = values
            return pd.DataFrame(full).ffill().bfill().to_numpy()
        states = self._model.predict(symbols)
        full = np.full(n, np.nan)
        full[valid] = states
        return pd.Series(full).ffill().bfill().fillna(0).to_numpy().astype(int)

    def viterbi(self, obs: np.ndarray) -> np.ndarray:
        """Most likely state path (T,) aligned to the input columns (missing carried)."""
        return self._decode(obs, proba=False)

    def posterior(self, obs: np.ndarray) -> np.ndarray:
        """Per-position state posteriors (T, n_states) aligned to the input columns."""
        return self._decode(obs, proba=True)


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
        "n_observed",
    ]
    if n_reads == 0:
        return pd.DataFrame(columns=columns), (hmm or BernoulliHMM())

    if hmm is None:
        hmm = BernoulliHMM().fit([reads[i] for i in range(n_reads)])
    protected_state = int(np.argmin(hmm.p_emit))
    width_ceiling = max_width if max_width is not None else n_positions
    anchor_lo = max(0, anchor_index - anchor_tolerance)
    anchor_hi = min(n_positions, anchor_index + anchor_tolerance + 1)

    def _absent_row(read_index: int, n_observed: int) -> dict[str, object]:
        return {
            "read_index": read_index,
            "footprint_present": False,
            "protected_start": pd.NA,
            "protected_end": pd.NA,
            "protected_width": 0,
            "footprint_score": float("nan"),
            "n_observed": n_observed,
        }

    rows: list[dict[str, object]] = []
    for i in range(n_reads):
        obs = reads[i]
        valid = ~np.isnan(obs)
        n_observed = int(valid.sum())
        # A footprint call requires actual methylation data near the anchor; otherwise
        # Viterbi decodes from the transition/start prior over missing positions and would
        # emit a spurious full-width footprint (S4 review finding).
        if not valid[anchor_lo:anchor_hi].any():
            rows.append(_absent_row(i, n_observed))
            continue
        path = hmm.viterbi(obs)
        gamma = hmm.posterior(obs)
        runs = _protected_runs(path, protected_state)
        best = None
        for start, end in runs:
            width = end - start
            overlaps_anchor = (
                start - anchor_tolerance <= anchor_index < end + anchor_tolerance
            )
            # the run must be backed by at least one observed position (not a pure
            # prior-fill over missing calls)
            if (
                overlaps_anchor
                and min_width <= width <= width_ceiling
                and valid[start:end].any()
            ):
                score = float(gamma[start:end, protected_state].mean())
                if best is None or width > best["protected_width"]:
                    best = {
                        "protected_start": int(start),
                        "protected_end": int(end),
                        "protected_width": int(width),
                        "footprint_score": score,
                    }
        if best is None:
            rows.append(_absent_row(i, n_observed))
        else:
            rows.append(
                {"read_index": i, "footprint_present": True, "n_observed": n_observed, **best}
            )
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
