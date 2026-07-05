"""Import external genomics tracks and compare them to DiMeLo signal (Q8, coverage).

DiMeLo data becomes far more interpretable next to orthogonal assays — MNase / ATAC /
CUT&Tag / ChIP coverage (bigWig or bedGraph). This module bins an external coverage track
over the same regions used for DiMeLo analysis and correlates the two binned signals
(Pearson + Spearman) so a locus- or genome-wide comparison is a couple of calls.

Coverage import (bigWig via ``pyBigWig``, bedGraph via overlap-weighted binning) is the
first modality; Hi-C (``cooler``) and the Hi-C-vs-joint-occupancy view (Q7 bridge) are a
follow-up that adds a heavier dependency.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as _scipy_stats

_REGION_COLUMNS = ("chromosome", "start", "end")


def _region_id(chromosome: str, start: int, end: int) -> str:
    return f"{chromosome}:{start}-{end}"


def bin_regions(regions: pd.DataFrame, *, bins: int = 1) -> pd.DataFrame:
    """Split each region into ``bins`` equal-width bins.

    ``regions`` needs ``chromosome, start, end`` (and optionally ``region_id``). Returns
    ``region_id, chromosome, bin, bin_start, bin_end`` (bin index 0..bins-1). The last bin
    absorbs any remainder so the bins exactly tile ``[start, end)``.
    """
    missing = set(_REGION_COLUMNS) - set(regions.columns)
    if missing:
        raise ValueError(f"regions requires columns: {', '.join(sorted(missing))}.")
    if bins < 1:
        raise ValueError("bins must be >= 1.")

    rows: list[dict[str, object]] = []
    for _, region in regions.iterrows():
        chromosome = str(region["chromosome"])
        start = int(region["start"])
        end = int(region["end"])
        region_id = (
            str(region["region_id"])
            if "region_id" in regions.columns and not pd.isna(region.get("region_id"))
            else _region_id(chromosome, start, end)
        )
        if end <= start:
            raise ValueError(f"region {region_id} has end <= start.")
        edges = np.linspace(start, end, bins + 1).astype(int)
        edges[-1] = end
        for b in range(bins):
            rows.append(
                {
                    "region_id": region_id,
                    "chromosome": chromosome,
                    "bin": b,
                    "bin_start": int(edges[b]),
                    "bin_end": int(edges[b + 1]),
                }
            )
    return pd.DataFrame(
        rows, columns=["region_id", "chromosome", "bin", "bin_start", "bin_end"]
    )


def import_bigwig_signal(
    bigwig_path: str | Path,
    regions: pd.DataFrame,
    *,
    bins: int = 1,
    stat: str = "mean",
) -> pd.DataFrame:
    """Binned external coverage from a bigWig aligned to ``regions`` (via ``pyBigWig``).

    Returns ``bin_regions`` output plus a ``signal`` column (the per-bin ``stat`` over the
    bigWig; ``NaN`` where the track has no data). Bins with no coverage are ``NaN``.
    """
    import pyBigWig

    binned = bin_regions(regions, bins=bins)
    if binned.empty:
        binned["signal"] = pd.Series(dtype=float)
        return binned

    bigwig = pyBigWig.open(str(bigwig_path))
    try:
        available = set(bigwig.chroms().keys())
        signals: list[float] = []
        for _, region in binned.iterrows():
            chromosome = region["chromosome"]
            if chromosome not in available:
                signals.append(float("nan"))
                continue
            values = bigwig.stats(
                chromosome,
                int(region["bin_start"]),
                int(region["bin_end"]),
                type=stat,
            )
            value = values[0] if values else None
            signals.append(float("nan") if value is None else float(value))
    finally:
        bigwig.close()
    binned["signal"] = signals
    return binned


def import_bedgraph_signal(
    bedgraph_path: str | Path,
    regions: pd.DataFrame,
    *,
    bins: int = 1,
) -> pd.DataFrame:
    """Binned external coverage from a bedGraph aligned to ``regions``.

    Each bin's signal is the overlap-length-weighted mean of the bedGraph intervals it
    overlaps (``NaN`` where nothing overlaps).
    """
    binned = bin_regions(regions, bins=bins)
    if binned.empty:
        binned["signal"] = pd.Series(dtype=float)
        return binned

    bedgraph = pd.read_csv(
        bedgraph_path,
        sep="\t",
        header=None,
        comment="#",
        names=["chromosome", "start", "end", "value"],
    )
    bedgraph["chromosome"] = bedgraph["chromosome"].astype(str)

    signals: list[float] = []
    for _, region in binned.iterrows():
        chromosome = region["chromosome"]
        bin_start = int(region["bin_start"])
        bin_end = int(region["bin_end"])
        overlapping = bedgraph.loc[
            (bedgraph["chromosome"] == chromosome)
            & (bedgraph["end"] > bin_start)
            & (bedgraph["start"] < bin_end)
        ]
        if overlapping.empty:
            signals.append(float("nan"))
            continue
        overlap_len = (
            np.minimum(overlapping["end"], bin_end)
            - np.maximum(overlapping["start"], bin_start)
        ).clip(lower=0)
        total = float(overlap_len.sum())
        if total <= 0:
            signals.append(float("nan"))
            continue
        signals.append(float((overlapping["value"] * overlap_len).sum() / total))
    binned["signal"] = signals
    return binned


def correlate_binned_signals(
    dimelo_signal: pd.DataFrame,
    external_signal: pd.DataFrame,
    *,
    keys: tuple[str, ...] = ("region_id", "bin"),
    dimelo_value: str = "signal",
    external_value: str = "signal",
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Correlate two binned signals matched on ``keys``.

    Returns ``(paired, stats)`` where ``paired`` has ``*keys, dimelo, external`` over bins
    present in both (NaN dropped), and ``stats`` has ``n, pearson_r, pearson_p,
    spearman_rho, spearman_p`` (correlations NaN when fewer than 2 paired points or a
    signal is constant).
    """
    keys = tuple(keys)
    left = dimelo_signal.loc[:, [*keys, dimelo_value]].rename(
        columns={dimelo_value: "dimelo"}
    )
    right = external_signal.loc[:, [*keys, external_value]].rename(
        columns={external_value: "external"}
    )
    paired = (
        left.merge(right, on=list(keys), how="inner")
        .dropna(subset=["dimelo", "external"])
        .reset_index(drop=True)
    )
    stats: dict[str, float] = {
        "n": int(len(paired)),
        "pearson_r": float("nan"),
        "pearson_p": float("nan"),
        "spearman_rho": float("nan"),
        "spearman_p": float("nan"),
    }
    if len(paired) >= 2 and paired["dimelo"].std() > 0 and paired["external"].std() > 0:
        pearson_r, pearson_p = _scipy_stats.pearsonr(paired["dimelo"], paired["external"])
        spearman_rho, spearman_p = _scipy_stats.spearmanr(
            paired["dimelo"], paired["external"]
        )
        stats.update(
            pearson_r=float(pearson_r),
            pearson_p=float(pearson_p),
            spearman_rho=float(spearman_rho),
            spearman_p=float(spearman_p),
        )
    return paired, stats
