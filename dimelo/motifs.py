"""Sequence-motif analysis at DiMeLo sites and comparison to annotated motifs (Q8).

Given a position frequency matrix (PFM) — e.g. from ``dimelo.jaspar`` — and the genomic
sequences under a set of DiMeLo sites, this module:

- scans each site sequence for the motif (log-odds PSSM, best hit),
- builds the *observed* motif from the site sequences (a position frequency matrix), and
- compares the observed motif to the annotated PFM (per-position and overall correlation),

so you can ask "do my DiMeLo sites actually carry the annotated binding motif, and does the
observed logo match?". Sequence extraction from a FASTA uses ``pysam`` (an existing
dependency); all scoring/logo/comparison functions work offline on plain sequences.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

_BASES = ("A", "C", "G", "T")
_BASE_INDEX = {base: index for index, base in enumerate(_BASES)}


def counts_to_probability(
    counts: pd.DataFrame, *, pseudocount: float = 0.8
) -> pd.DataFrame:
    """Convert a per-position ``A/C/G/T`` count table to a probability matrix (rows sum to
    1) with additive ``pseudocount`` smoothing."""
    matrix = counts.loc[:, list(_BASES)].to_numpy(dtype=float) + pseudocount
    probability = matrix / matrix.sum(axis=1, keepdims=True)
    frame = pd.DataFrame(probability, columns=list(_BASES))
    if "position" in counts.columns:
        frame.insert(0, "position", counts["position"].to_numpy())
    else:
        frame.insert(0, "position", range(1, len(frame) + 1))
    return frame


def information_content(probability: pd.DataFrame) -> pd.DataFrame:
    """Per-(position, base) information-content heights (bits) for a sequence logo.

    Height = base probability x position information content (``2 - entropy``, capped at 2
    bits for DNA). Returns ``position, base, bits``.
    """
    prob = probability.loc[:, list(_BASES)].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        per_base = np.where(prob > 0, prob * np.log2(prob), 0.0)
    entropy = -np.sum(per_base, axis=1)
    info = np.log2(len(_BASES)) - entropy  # bits at this position (max 2 for DNA)
    heights = prob * info[:, None]
    positions = (
        probability["position"].to_numpy()
        if "position" in probability.columns
        else np.arange(1, prob.shape[0] + 1)
    )
    rows: list[dict[str, object]] = []
    for row_index in range(prob.shape[0]):
        for base_index, base in enumerate(_BASES):
            rows.append(
                {
                    "position": int(positions[row_index]),
                    "base": base,
                    "bits": float(heights[row_index, base_index]),
                }
            )
    return pd.DataFrame(rows, columns=["position", "base", "bits"])


def _log_odds(probability: pd.DataFrame, background: dict[str, float]) -> np.ndarray:
    prob = np.clip(probability.loc[:, list(_BASES)].to_numpy(dtype=float), 1e-9, None)
    bg = np.array([background[base] for base in _BASES], dtype=float)
    return np.log2(prob / bg[None, :])


def scan_sequences(
    sequences: Sequence[str],
    probability: pd.DataFrame,
    *,
    background: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Best forward-strand log-odds (PSSM) match per sequence.

    Slides the motif over each sequence (all offsets), scoring only windows of full ACGT
    (any other character voids that offset). Returns ``sequence_index, best_score,
    best_offset`` (``best_score`` is ``-inf`` and ``best_offset`` ``-1`` when no scorable
    window exists).
    """
    background = background or {base: 0.25 for base in _BASES}
    pssm = _log_odds(probability, background)
    motif_length = pssm.shape[0]

    rows: list[dict[str, object]] = []
    for sequence_index, sequence in enumerate(sequences):
        seq = str(sequence).upper()
        best_score = float("-inf")
        best_offset = -1
        for offset in range(0, len(seq) - motif_length + 1):
            window = seq[offset : offset + motif_length]
            score = 0.0
            scorable = True
            for column, char in enumerate(window):
                base_index = _BASE_INDEX.get(char)
                if base_index is None:
                    scorable = False
                    break
                score += pssm[column, base_index]
            if scorable and score > best_score:
                best_score = score
                best_offset = offset
        rows.append(
            {
                "sequence_index": sequence_index,
                "best_score": best_score,
                "best_offset": best_offset,
            }
        )
    return pd.DataFrame(rows, columns=["sequence_index", "best_score", "best_offset"])


def observed_pfm(sequences: Sequence[str]) -> pd.DataFrame:
    """Position frequency matrix (counts) from equal-length sequences (the observed motif).

    Non-ACGT characters are ignored. Returns ``position, A, C, G, T``.
    """
    seqs = [str(s).upper() for s in sequences]
    if not seqs:
        raise ValueError("observed_pfm requires at least one sequence.")
    length = len(seqs[0])
    if any(len(s) != length for s in seqs):
        raise ValueError("observed_pfm requires equal-length sequences.")
    counts = np.zeros((length, len(_BASES)), dtype=float)
    for seq in seqs:
        for column, char in enumerate(seq):
            base_index = _BASE_INDEX.get(char)
            if base_index is not None:
                counts[column, base_index] += 1
    frame = pd.DataFrame(counts, columns=list(_BASES))
    frame.insert(0, "position", range(1, length + 1))
    return frame


def compare_motifs(
    observed_probability: pd.DataFrame, reference_probability: pd.DataFrame
) -> dict[str, float]:
    """Similarity of two equal-length probability matrices.

    Returns ``pearson_r`` / ``pearson_p`` over the flattened matrices and
    ``mean_per_position_r`` (mean of per-position Pearson correlations, ignoring positions
    where either matrix is constant).
    """
    from scipy.stats import pearsonr

    observed = observed_probability.loc[:, list(_BASES)].to_numpy(dtype=float)
    reference = reference_probability.loc[:, list(_BASES)].to_numpy(dtype=float)
    if observed.shape != reference.shape:
        raise ValueError("compare_motifs requires equal-length probability matrices.")
    pearson_r, pearson_p = pearsonr(observed.ravel(), reference.ravel())
    per_position = []
    for row_index in range(observed.shape[0]):
        a = observed[row_index]
        b = reference[row_index]
        if np.std(a) > 0 and np.std(b) > 0:
            per_position.append(float(pearsonr(a, b)[0]))
    mean_per_position = float(np.nanmean(per_position)) if per_position else float("nan")
    return {
        "pearson_r": float(pearson_r),
        "pearson_p": float(pearson_p),
        "mean_per_position_r": mean_per_position,
    }


def extract_site_sequences(
    fasta_path: str | Path,
    sites: pd.DataFrame,
    *,
    width: int,
    chromosome_column: str = "chromosome",
    center_column: str = "center",
) -> list[str]:
    """Extract fixed-``width`` sequences centered on each site from a FASTA (via ``pysam``).

    ``sites`` needs ``chromosome_column`` and ``center_column``. Returns uppercase strings
    of length ``width`` (windows falling off a contig end are skipped).
    """
    import pysam

    half = width // 2
    fasta = pysam.FastaFile(str(fasta_path))
    try:
        contig_lengths = dict(zip(fasta.references, fasta.lengths, strict=True))
        sequences: list[str] = []
        for _, site in sites.iterrows():
            chromosome = str(site[chromosome_column])
            center = int(site[center_column])
            start = center - half
            end = start + width
            if chromosome not in contig_lengths or start < 0 or end > contig_lengths[chromosome]:
                continue
            sequences.append(fasta.fetch(chromosome, start, end).upper())
    finally:
        fasta.close()
    return sequences
