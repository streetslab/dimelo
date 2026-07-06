"""ENCODE data-portal connector (Q8 reference connectors).

Search the ENCODE portal for experiments (by target / assay / biosample / assembly), list
their files, and download a file (cached). The fetched artifacts route into the comparison
machinery already built:

- ChIP/TF peak beds -> ``cluster_reference`` / ``regulatory_enrichment`` (validate binding
  sites, cluster association),
- signal bigWigs (histone marks, ATAC, MNase) -> ``tracks.import_bigwig_signal`` + correlate,
- Hi-C ``.cool`` -> ``tracks.import_hic_contacts``.

Network fetching (``search_*`` / ``fetch_*`` / ``download_*``, which lazily import
``requests``) is kept separate from the offline parsers (``parse_encode_search`` /
``parse_encode_files`` / ``load_encode_peaks``) so parsing and comparison are fully testable
without network access. Downloads are cached to disk.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ENCODE_BASE_URL = "https://www.encodeproject.org"
ENCODE_SEARCH_URL = f"{ENCODE_BASE_URL}/search/"
_JSON_HEADERS = {"Accept": "application/json"}


def parse_encode_search(payload: dict) -> pd.DataFrame:
    """Parse an ENCODE ``/search/`` JSON payload into an experiment table.

    Returns ``accession, target, assay_title, biosample, assembly, description`` (one row
    per experiment in ``@graph``).
    """
    graph = payload.get("@graph", [])
    rows: list[dict[str, object]] = []
    for experiment in graph:
        target = experiment.get("target") or {}
        biosample = experiment.get("biosample_ontology") or {}
        assembly = experiment.get("assembly")
        rows.append(
            {
                "accession": experiment.get("accession"),
                "target": target.get("label"),
                "assay_title": experiment.get("assay_title"),
                "biosample": biosample.get("term_name"),
                "assembly": assembly[0] if isinstance(assembly, list) and assembly else assembly,
                "description": experiment.get("description"),
            }
        )
    return pd.DataFrame(
        rows,
        columns=["accession", "target", "assay_title", "biosample", "assembly", "description"],
    )


def parse_encode_files(
    payload: dict,
    *,
    file_format: str | None = None,
    output_type: str | None = None,
    assembly: str | None = None,
) -> pd.DataFrame:
    """Parse an ENCODE experiment JSON payload into a downloadable-files table.

    Returns ``accession, file_format, output_type, assembly, file_size, href, url``
    (``url`` is the absolute download URL). Optional filters narrow by ``file_format``
    (e.g. ``"bigWig"``, ``"bed"``), ``output_type`` (e.g. ``"fold change over control"``,
    ``"IDR thresholded peaks"``), and ``assembly`` (e.g. ``"GRCh38"``).
    """
    files = payload.get("files", [])
    rows: list[dict[str, object]] = []
    for file_record in files:
        href = file_record.get("href")
        rows.append(
            {
                "accession": file_record.get("accession"),
                "file_format": file_record.get("file_format"),
                "output_type": file_record.get("output_type"),
                "assembly": file_record.get("assembly"),
                "file_size": file_record.get("file_size"),
                "href": href,
                "url": f"{ENCODE_BASE_URL}{href}" if href else None,
            }
        )
    table = pd.DataFrame(
        rows,
        columns=[
            "accession",
            "file_format",
            "output_type",
            "assembly",
            "file_size",
            "href",
            "url",
        ],
    )
    if file_format is not None:
        table = table.loc[table["file_format"] == file_format]
    if output_type is not None:
        table = table.loc[table["output_type"] == output_type]
    if assembly is not None:
        table = table.loc[table["assembly"] == assembly]
    return table.reset_index(drop=True)


def load_encode_peaks(path: str | Path, *, max_peaks: int | None = None) -> pd.DataFrame:
    """Load a BED / narrowPeak file (``.gz`` ok) into a ``chromosome, start, end`` table
    ready for overlap (``cluster_reference`` / ``tracks``). ``max_peaks`` keeps the first N."""
    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        comment="#",
        compression="infer",
        usecols=[0, 1, 2],
        names=["chromosome", "start", "end"],
    )
    # Coerce coords and drop non-data rows: UCSC 'track'/'browser' header lines are not
    # '#'-commented, so read_csv keeps them; they coerce to NaN in start/end and are
    # dropped here (rather than crashing astype(int) with IntCastingNaNError).
    table["start"] = pd.to_numeric(table["start"], errors="coerce")
    table["end"] = pd.to_numeric(table["end"], errors="coerce")
    table = table.dropna(subset=["start", "end"]).copy()
    table["chromosome"] = table["chromosome"].astype(str)
    table["start"] = table["start"].astype(int)
    table["end"] = table["end"].astype(int)
    if max_peaks is not None:
        table = table.head(max_peaks)
    return table.reset_index(drop=True)


def search_encode_experiments(
    *,
    target: str | None = None,
    assay_title: str | None = None,
    biosample: str | None = None,
    assembly: str | None = None,
    limit: int = 25,
    timeout: int = 30,
) -> pd.DataFrame:
    """Search ENCODE for experiments (network). Returns :func:`parse_encode_search`."""
    import requests

    params: dict[str, str] = {"type": "Experiment", "format": "json", "limit": str(limit)}
    if target is not None:
        params["target.label"] = target
    if assay_title is not None:
        params["assay_title"] = assay_title
    if biosample is not None:
        params["biosample_ontology.term_name"] = biosample
    if assembly is not None:
        params["assembly"] = assembly
    response = requests.get(
        ENCODE_SEARCH_URL, params=params, headers=_JSON_HEADERS, timeout=timeout
    )
    response.raise_for_status()
    return parse_encode_search(response.json())


def fetch_encode_experiment_files(
    accession: str,
    *,
    file_format: str | None = None,
    output_type: str | None = None,
    assembly: str | None = None,
    timeout: int = 30,
) -> pd.DataFrame:
    """Fetch an ENCODE experiment's file list (network). Returns
    :func:`parse_encode_files` with the given filters applied."""
    import requests

    response = requests.get(
        f"{ENCODE_BASE_URL}/experiments/{accession}/",
        headers=_JSON_HEADERS,
        timeout=timeout,
    )
    response.raise_for_status()
    return parse_encode_files(
        response.json(),
        file_format=file_format,
        output_type=output_type,
        assembly=assembly,
    )


def download_encode_file(
    url: str,
    *,
    cache_dir: str | Path | None = None,
    timeout: int = 120,
    refresh: bool = False,
) -> Path:
    """Download an ENCODE file to ``cache_dir`` (streamed), returning the local path.

    A cache hit (and ``refresh=False``) makes no network request. ``url`` may be an
    absolute ENCODE URL or a portal ``href`` (``/files/.../@@download/NAME``).
    """
    if url.startswith("/"):
        url = f"{ENCODE_BASE_URL}{url}"
    cache_root = (
        Path(cache_dir)
        if cache_dir is not None
        else Path.home() / ".cache" / "dimelo" / "encode"
    )
    cache_root.mkdir(parents=True, exist_ok=True)
    filename = url.rstrip("/").split("/")[-1]
    if not filename:
        raise ValueError(f"Could not derive a filename from ENCODE url: {url!r}.")
    destination = cache_root / filename
    if destination.exists() and not refresh:
        return destination

    import requests

    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        # write to a temp file first so an interrupted download doesn't leave a bad cache
        partial = destination.with_suffix(destination.suffix + ".partial")
        with open(partial, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1 << 20):
                if chunk:
                    handle.write(chunk)
        partial.replace(destination)
    return destination
