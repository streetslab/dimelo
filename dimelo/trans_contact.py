"""Trans / long-range-contact correction (Q4) — exploratory research spike.

DiMeLo's tethered methyltransferase can methylate DNA that is far in 1D sequence but close
in 3D (a trans / long-range contact), which inflates apparent in-site signal. Because trans
methylation is target-driven it is *not* removed by the negative-control model (Q3). This
module separates in-site (cis) methylation from contact-driven (trans) methylation using two
pieces of independent evidence:

1. **Cis reach** — methylation within the enzyme's 1D reach of a bound-site footprint is
   cis. The reach is estimated empirically from the methylation-density decay around
   isolated, strongly-bound sites (``estimate_cis_reach``), or supplied by the user.
2. **Hi-C-informed trans flag** — a methylated position *outside* the cis reach is flagged
   trans when it has high Hi-C contact to a Q3-called bound site (3D-proximal, 1D-distal).

Outputs: a per-position cis / trans / distal-background classification, an in-site-only
corrected signal, and a separate trans track. This is a spike — validate that trans-flagging
reduces spurious Q7 joint occupancy at high-Hi-C distal pairs before relying on it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

CIS = "cis"
TRANS = "trans"
DISTAL_BACKGROUND = "distal_background"


def estimate_cis_reach(
    profile: pd.DataFrame,
    *,
    offset_column: str = "offset_bp",
    density_column: str = "mean_density",
    background_level: float | None = None,
    decay_fraction: float = 0.5,
) -> int:
    """Estimate the enzyme's 1D cis reach from a methylation-density decay profile.

    ``profile`` is methylation density vs absolute offset from isolated, strongly-bound
    sites (``offset_bp >= 0``, e.g. an aggregate around footprints). The reach is the first
    offset (moving out from 0) at which density has decayed from its peak toward background
    by ``decay_fraction`` — i.e. crosses ``background + (1 - decay_fraction) * (peak -
    background)``. ``background_level`` defaults to the density at the largest offset.
    """
    _require = {offset_column, density_column}
    missing = _require - set(profile.columns)
    if missing:
        raise ValueError(f"profile requires columns: {', '.join(sorted(missing))}.")
    if not 0.0 < decay_fraction < 1.0:
        raise ValueError("decay_fraction must be in (0, 1).")
    ordered = profile.sort_values(offset_column).reset_index(drop=True)
    offsets = ordered[offset_column].to_numpy(dtype=float)
    density = ordered[density_column].to_numpy(dtype=float)
    if len(ordered) == 0:
        raise ValueError("profile is empty.")

    peak = float(np.nanmax(density))
    background = (
        float(background_level) if background_level is not None else float(density[-1])
    )
    if peak <= background:
        return int(offsets[-1])
    threshold = background + (1.0 - decay_fraction) * (peak - background)
    # first offset where density falls below the threshold
    below = np.where(density < threshold)[0]
    if len(below) == 0:
        return int(offsets[-1])
    return int(offsets[below[0]])


def classify_cis_trans(
    positions: pd.DataFrame,
    bound_sites: pd.DataFrame,
    *,
    cis_reach_bp: int,
    hic_contact_threshold: float | None = None,
) -> pd.DataFrame:
    """Classify each methylated position as cis / trans / distal-background.

    ``positions`` needs ``chromosome, position`` (and typically a ``signal`` column and, for
    trans calling, a per-position ``hic_contact`` = max Hi-C contact to a bound site).
    ``bound_sites`` needs ``chromosome, position`` (Q3-called bound anchors). A position
    within ``cis_reach_bp`` of a bound site on the same chromosome is ``cis``; a more distal
    position with ``hic_contact >= hic_contact_threshold`` is ``trans``; otherwise
    ``distal_background``. Returns ``positions`` with ``nearest_bound_distance`` and
    ``classification`` columns.
    """
    for column in ("chromosome", "position"):
        if column not in positions.columns:
            raise ValueError(f"positions requires column: {column}.")
        if column not in bound_sites.columns:
            raise ValueError(f"bound_sites requires column: {column}.")
    if hic_contact_threshold is not None and "hic_contact" not in positions.columns:
        raise ValueError(
            "positions requires a 'hic_contact' column when hic_contact_threshold is set."
        )

    out = positions.copy().reset_index(drop=True)
    bound_by_chrom = {
        str(chrom): np.sort(group["position"].to_numpy(dtype=float))
        for chrom, group in bound_sites.groupby("chromosome")
    }

    distances = np.full(len(out), np.inf)
    for i, (chrom, pos) in enumerate(
        zip(out["chromosome"].astype(str), out["position"].astype(float), strict=True)
    ):
        anchors = bound_by_chrom.get(chrom)
        if anchors is None or len(anchors) == 0:
            continue
        distances[i] = float(np.min(np.abs(anchors - pos)))
    out["nearest_bound_distance"] = distances

    is_cis = distances <= cis_reach_bp
    if hic_contact_threshold is not None:
        hic = pd.to_numeric(out["hic_contact"], errors="coerce").to_numpy()
        is_trans = (~is_cis) & (hic >= hic_contact_threshold)
    else:
        is_trans = np.zeros(len(out), dtype=bool)

    classification = np.full(len(out), DISTAL_BACKGROUND, dtype=object)
    classification[is_cis] = CIS
    classification[is_trans] = TRANS
    out["classification"] = classification
    return out


def apply_trans_correction(
    classified: pd.DataFrame,
    *,
    signal_column: str = "signal",
) -> pd.DataFrame:
    """Split the signal into in-site (cis) and trans tracks.

    Adds ``in_site_signal`` (the signal at ``cis`` positions, 0 elsewhere) and
    ``trans_signal`` (the signal at ``trans`` positions, 0 elsewhere). Distal-background
    positions contribute to neither. The in-site track is the trans-corrected signal; the
    trans track is the removed contact-driven signal (useful for architecture / Q10).
    """
    if "classification" not in classified.columns:
        raise ValueError("classified requires a 'classification' column.")
    if signal_column not in classified.columns:
        raise ValueError(f"classified requires the signal column: {signal_column}.")
    out = classified.copy()
    signal = pd.to_numeric(out[signal_column], errors="coerce").fillna(0.0)
    out["in_site_signal"] = np.where(out["classification"] == CIS, signal, 0.0)
    out["trans_signal"] = np.where(out["classification"] == TRANS, signal, 0.0)
    return out


def summarize_trans_correction(corrected: pd.DataFrame) -> pd.DataFrame:
    """Totals of the raw / in-site / trans signal and per-class position counts.

    Returns a one-row frame: ``n_positions, n_cis, n_trans, n_distal_background,
    total_in_site_signal, total_trans_signal, trans_fraction`` (trans_fraction = trans
    signal / (in-site + trans signal))."""
    required = {"classification", "in_site_signal", "trans_signal"}
    missing = required - set(corrected.columns)
    if missing:
        raise ValueError(
            f"summarize_trans_correction requires columns: {', '.join(sorted(missing))}."
        )
    in_site_total = float(corrected["in_site_signal"].sum())
    trans_total = float(corrected["trans_signal"].sum())
    denom = in_site_total + trans_total
    counts = corrected["classification"].value_counts()
    return pd.DataFrame(
        [
            {
                "n_positions": int(len(corrected)),
                "n_cis": int(counts.get(CIS, 0)),
                "n_trans": int(counts.get(TRANS, 0)),
                "n_distal_background": int(counts.get(DISTAL_BACKGROUND, 0)),
                "total_in_site_signal": in_site_total,
                "total_trans_signal": trans_total,
                "trans_fraction": (trans_total / denom) if denom > 0 else 0.0,
            }
        ]
    )
