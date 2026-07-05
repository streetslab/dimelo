"""JASPAR motif (PWM) fetching — a thin cached client (Q8 reference connectors).

Pulls a position frequency matrix (PFM) for a JASPAR matrix id (e.g. ``MA0139.1`` for CTCF)
from the JASPAR REST API and caches the raw JSON to disk, so repeated use and offline work
hit the cache instead of the network. Parsing (``jaspar_matrix_from_json``) is separated
from fetching so it is fully testable without network access.

The returned PFM (a per-position ``A/C/G/T`` count table) feeds ``dimelo.motifs`` to scan
DiMeLo-site sequences, build an observed sequence logo, and compare it to the annotated
motif.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

JASPAR_MATRIX_URL = "https://jaspar.genereg.net/api/v1/matrix/{matrix_id}/"
_BASES = ("A", "C", "G", "T")


def jaspar_matrix_from_json(payload: dict) -> pd.DataFrame:
    """Parse a JASPAR matrix-API JSON payload into a per-position count table.

    Returns a ``DataFrame`` with ``position, A, C, G, T`` (1-based positions) and
    ``matrix_id`` / ``name`` recorded in ``.attrs``.
    """
    if "pfm" not in payload:
        raise ValueError("JASPAR payload is missing the 'pfm' field.")
    pfm = payload["pfm"]
    missing = [base for base in _BASES if base not in pfm]
    if missing:
        raise ValueError(f"JASPAR pfm is missing base(s): {', '.join(missing)}.")
    lengths = {len(pfm[base]) for base in _BASES}
    if len(lengths) != 1:
        raise ValueError("JASPAR pfm columns have inconsistent lengths.")
    frame = pd.DataFrame({base: [float(x) for x in pfm[base]] for base in _BASES})
    frame.insert(0, "position", range(1, len(frame) + 1))
    frame.attrs["matrix_id"] = payload.get("matrix_id")
    frame.attrs["name"] = payload.get("name")
    return frame


def fetch_jaspar_matrix(
    matrix_id: str,
    *,
    cache_dir: str | Path | None = None,
    timeout: int = 30,
    refresh: bool = False,
) -> pd.DataFrame:
    """Fetch a JASPAR PFM by matrix id, caching the raw JSON under ``cache_dir``.

    On a cache hit (and ``refresh=False``) no network request is made. Returns the parsed
    count table from :func:`jaspar_matrix_from_json`.
    """
    cache_root = (
        Path(cache_dir)
        if cache_dir is not None
        else Path.home() / ".cache" / "dimelo" / "jaspar"
    )
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_file = cache_root / f"{matrix_id}.json"

    if cache_file.exists() and not refresh:
        try:
            return jaspar_matrix_from_json(json.loads(cache_file.read_text()))
        except (ValueError, json.JSONDecodeError):
            # corrupt / previously-poisoned cache -> fall through and refetch
            pass

    import requests

    response = requests.get(
        JASPAR_MATRIX_URL.format(matrix_id=matrix_id), timeout=timeout
    )
    response.raise_for_status()
    payload = response.json()
    # validate BEFORE writing the cache so a malformed 200 never poisons it
    frame = jaspar_matrix_from_json(payload)
    cache_file.write_text(json.dumps(payload))
    return frame
