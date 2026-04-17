from __future__ import annotations

from dataclasses import dataclass, field
import html
import http.cookiejar
from pathlib import Path
import re
import shutil
import subprocess
import time
from typing import Any, Iterable
from urllib import error, parse, request

from .models import UniBindJobResult

_UNIBIND_SUPPORTED_SPECIES = (
    "arabidopsis_thaliana",
    "caenorhabditis_elegans",
    "danio_rerio",
    "drosophila_melanogaster",
    "homo_sapiens",
    "mus_musculus",
    "rattus_norvegicus",
    "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe",
)

_SCREEN_SUPPORTED_SPECIES = {
    "homo_sapiens",
    "mus_musculus",
}

_DEFAULT_SPECIES = "homo_sapiens"
_DEFAULT_TARGET_GENOME_BY_SPECIES = {
    "homo_sapiens": "hg38",
    "mus_musculus": "mm10",
}

_PROVIDER_ALIASES = {
    "screen": "screen",
    "screen2": "screen",
    "screen_2": "screen",
    "screen2_0": "screen",
    "screen_2_0": "screen",
    "unibind": "unibind",
    "unibind_2021": "unibind",
    "unibind2021": "unibind",
}

_SPECIES_ALIASES = {
    "arabidopsis_thaliana": "arabidopsis_thaliana",
    "arabidopsis": "arabidopsis_thaliana",
    "3702": "arabidopsis_thaliana",
    "caenorhabditis_elegans": "caenorhabditis_elegans",
    "c_elegans": "caenorhabditis_elegans",
    "celegans": "caenorhabditis_elegans",
    "worm": "caenorhabditis_elegans",
    "6239": "caenorhabditis_elegans",
    "danio_rerio": "danio_rerio",
    "zebrafish": "danio_rerio",
    "7955": "danio_rerio",
    "drosophila_melanogaster": "drosophila_melanogaster",
    "d_melanogaster": "drosophila_melanogaster",
    "fly": "drosophila_melanogaster",
    "fruit_fly": "drosophila_melanogaster",
    "7227": "drosophila_melanogaster",
    "homo_sapiens": "homo_sapiens",
    "human": "homo_sapiens",
    "9606": "homo_sapiens",
    "mus_musculus": "mus_musculus",
    "mouse": "mus_musculus",
    "10090": "mus_musculus",
    "rattus_norvegicus": "rattus_norvegicus",
    "rat": "rattus_norvegicus",
    "10116": "rattus_norvegicus",
    "saccharomyces_cerevisiae": "saccharomyces_cerevisiae",
    "s_cerevisiae": "saccharomyces_cerevisiae",
    "yeast": "saccharomyces_cerevisiae",
    "4932": "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe": "schizosaccharomyces_pombe",
    "s_pombe": "schizosaccharomyces_pombe",
    "4896": "schizosaccharomyces_pombe",
}

_GENOME_SPECIES_HINTS = {
    "hg38": "homo_sapiens",
    "grch38": "homo_sapiens",
    "hg19": "homo_sapiens",
    "grch37": "homo_sapiens",
    "hs1": "homo_sapiens",
    "chm13": "homo_sapiens",
    "t2t_chm13v2_0": "homo_sapiens",
    "t2t_chm13v2": "homo_sapiens",
    "mm10": "mus_musculus",
    "grcm38": "mus_musculus",
    "mm39": "mus_musculus",
    "grcm39": "mus_musculus",
}

_UNIBIND_TRACKHUB_URLS = {
    "robust": "https://unibind.uio.no/static/data/latest/UniBind_hubs_Robust/UCSC/hub.txt",
    "permissive": "https://unibind.uio.no/static/data/latest/UniBind_hubs_Permissive/UCSC/hub.txt",
}

DEFAULT_UNIBIND_TFBS_EXTRACTION_URL = "https://unibind.uio.no/TFBS_extraction/"
DEFAULT_UNIBIND_ENRICHMENT_URL = "https://unibind.uio.no/enrichment/"

_UNIBIND_FORM_SPECIES_CODES = {
    "arabidopsis_thaliana": "3702",
    "caenorhabditis_elegans": "6239",
    "danio_rerio": "7955",
    "drosophila_melanogaster": "7227",
    "homo_sapiens": "9606",
    "mus_musculus": "10090",
    "rattus_norvegicus": "10116",
    "saccharomyces_cerevisiae": "4932",
    "schizosaccharomyces_pombe": "4896",
}
_UNIBIND_FORM_SPECIES_BY_CODE = {v: k for k, v in _UNIBIND_FORM_SPECIES_CODES.items()}
_UNIBIND_COLLECTION_ALIASES = {
    "robust": "Robust",
    "permissive": "Permissive",
}
_UNIBIND_ENRICHMENT_ANALYSIS_ALIASES = {
    "onesetbg": "oneSetBg",
    "one_set_bg": "oneSetBg",
    "with_background": "oneSetBg",
    "twosets": "twoSets",
    "two_sets": "twoSets",
    "differential": "twoSets",
}
_UNIBIND_STATUS_ALIASES = {
    "queued": "queued",
    "pending": "queued",
    "running": "running",
    "submitted": "queued",
    "processing": "running",
    "completed": "completed",
    "finished": "completed",
    "done": "completed",
    "success": "completed",
    "failed": "error",
    "failure": "error",
    "error": "error",
    "cancelled": "cancelled",
    "canceled": "cancelled",
}
_UNIBIND_TERMINAL_STATUSES = {"completed", "error", "cancelled"}
_UNIBIND_SUCCESS_STATUSES = {"completed"}


class RegulatoryEnrichmentSpecError(ValueError):
    """Invalid regulatory enrichment provider configuration."""


def _normalize_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def _normalize_assembly_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.strip().lower())


def normalize_species_name(species: str) -> str:
    token = _normalize_token(species)
    if not token:
        raise RegulatoryEnrichmentSpecError("Species cannot be empty.")
    canonical = _SPECIES_ALIASES.get(token, token)
    if canonical not in _UNIBIND_SUPPORTED_SPECIES:
        supported = ", ".join(_UNIBIND_SUPPORTED_SPECIES)
        raise RegulatoryEnrichmentSpecError(
            "Unsupported species for UniBind-backed workflows: "
            f"{species!r}. Supported species: {supported}."
        )
    return canonical


def _infer_species_from_genomes(*, reference_genome: str | None, target_genome: str | None) -> str | None:
    for genome in (target_genome, reference_genome):
        if genome is None:
            continue
        token = _normalize_token(str(genome))
        if token in _GENOME_SPECIES_HINTS:
            return _GENOME_SPECIES_HINTS[token]
    return None


def _normalize_provider(provider: str) -> str:
    token = _normalize_token(provider)
    if token not in _PROVIDER_ALIASES:
        supported = ", ".join(sorted(set(_PROVIDER_ALIASES.values())))
        raise RegulatoryEnrichmentSpecError(
            f"Unsupported provider {provider!r}. Supported providers: {supported}."
        )
    return _PROVIDER_ALIASES[token]


def _unique_preserve_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _normalize_unibind_species_code(species: str) -> str:
    raw = str(species).strip()
    if not raw:
        raise ValueError("UniBind species cannot be empty.")
    if raw in _UNIBIND_FORM_SPECIES_BY_CODE:
        return raw
    token = _normalize_token(raw)
    if token in _UNIBIND_FORM_SPECIES_BY_CODE:
        return token
    canonical = normalize_species_name(raw)
    return _UNIBIND_FORM_SPECIES_CODES[canonical]


def _normalize_unibind_collection(collection: str) -> str:
    token = _normalize_token(collection)
    if token not in _UNIBIND_COLLECTION_ALIASES:
        allowed = ", ".join(sorted(_UNIBIND_COLLECTION_ALIASES.values()))
        raise ValueError(
            f"Unsupported UniBind collection {collection!r}. Supported: {allowed}."
        )
    return _UNIBIND_COLLECTION_ALIASES[token]


def _normalize_unibind_enrichment_analysis_type(analysis_type: str) -> str:
    token = _normalize_token(analysis_type)
    if token not in _UNIBIND_ENRICHMENT_ANALYSIS_ALIASES:
        allowed = ", ".join(sorted(set(_UNIBIND_ENRICHMENT_ANALYSIS_ALIASES.values())))
        raise ValueError(
            "Unsupported UniBind enrichment analysis_type "
            f"{analysis_type!r}. Supported: {allowed}."
        )
    return _UNIBIND_ENRICHMENT_ANALYSIS_ALIASES[token]


def _normalize_unibind_status(status: str | None) -> str:
    if status is None:
        return "unknown"
    token = _normalize_token(str(status))
    if not token:
        return "unknown"
    return _UNIBIND_STATUS_ALIASES.get(token, token)


def _regions_to_bed_bytes(regions: Any) -> bytes:
    from . import chip_atlas

    text = chip_atlas._regions_to_bed_text(regions)
    return text.encode("utf-8")


def _coerce_optional_text_file_payload(
    value: str | Path | Iterable[str] | None,
    *,
    default_filename: str,
) -> tuple[str, bytes, str] | None:
    if value is None:
        return None

    filename = default_filename
    if isinstance(value, Path):
        payload = value.read_bytes()
        filename = value.name
    elif isinstance(value, str):
        candidate = Path(value).expanduser()
        if "\n" not in value and "\t" not in value and candidate.exists():
            payload = candidate.read_bytes()
            filename = candidate.name
        else:
            payload = value.encode("utf-8")
    else:
        lines = [str(item).strip() for item in value if str(item).strip()]
        payload = ("\n".join(lines) + ("\n" if lines else "")).encode("utf-8")

    if not payload:
        return None
    return filename, payload, "text/plain"


def _build_multipart_form(
    *,
    fields: dict[str, str],
    files: dict[str, tuple[str, bytes, str]],
) -> tuple[bytes, str]:
    boundary = f"----dimelo-{int(time.time() * 1000)}"
    chunks: list[bytes] = []
    for key, value in fields.items():
        chunks.append(f"--{boundary}\r\n".encode("utf-8"))
        chunks.append(
            f'Content-Disposition: form-data; name="{key}"\r\n\r\n'.encode("utf-8")
        )
        chunks.append(str(value).encode("utf-8"))
        chunks.append(b"\r\n")
    for key, (filename, payload, content_type) in files.items():
        chunks.append(f"--{boundary}\r\n".encode("utf-8"))
        chunks.append(
            (
                f'Content-Disposition: form-data; name="{key}"; '
                f'filename="{filename}"\r\n'
            ).encode("utf-8")
        )
        chunks.append(f"Content-Type: {content_type}\r\n\r\n".encode("utf-8"))
        chunks.append(payload)
        chunks.append(b"\r\n")
    chunks.append(f"--{boundary}--\r\n".encode("utf-8"))
    return b"".join(chunks), boundary


def _extract_csrf_token(page_html: str) -> str:
    match = re.search(
        r'name="csrfmiddlewaretoken"\s+value="([^"]+)"',
        page_html,
        flags=re.IGNORECASE,
    )
    if match is None:
        raise RuntimeError("UniBind form did not include a CSRF token.")
    return html.unescape(match.group(1))


def _extract_unibind_results_url(job_html: str, *, base_url: str) -> str | None:
    match = re.search(
        r"<th>\s*Results URL:\s*</th>\s*<td>\s*<a href=\"([^\"]+)\"",
        job_html,
        flags=re.IGNORECASE,
    )
    if match is None:
        return None
    return parse.urljoin(base_url, html.unescape(match.group(1)).strip())


def _extract_unibind_job_status(job_html: str) -> str:
    match = re.search(
        r"<th>\s*Job status:\s*</th>\s*<td>\s*([^<]+)\s*</td>",
        job_html,
        flags=re.IGNORECASE,
    )
    if match is None:
        return "unknown"
    raw = html.unescape(match.group(1)).strip()
    return _normalize_unibind_status(raw)


def _looks_like_unibind_output_url(url: str) -> bool:
    parsed = parse.urlparse(url)
    path = parsed.path.lower()
    if "/temp/" in path:
        return True
    return path.endswith(
        (
            ".bed",
            ".bed.gz",
            ".tsv",
            ".csv",
            ".txt",
            ".zip",
            ".pdf",
        )
    )


def _extract_unibind_output_urls(job_html: str, *, base_url: str) -> list[str]:
    links: list[str] = []
    for href in re.findall(r'href="([^"]+)"', job_html, flags=re.IGNORECASE):
        resolved = parse.urljoin(base_url, html.unescape(href).strip())
        if not _looks_like_unibind_output_url(resolved):
            continue
        links.append(resolved)
    return _unique_preserve_order(links)


def _extract_unibind_job_id(job_url: str) -> str:
    parts = [segment for segment in parse.urlparse(job_url).path.split("/") if segment]
    if not parts:
        raise ValueError(f"Could not determine UniBind job id from URL: {job_url}")
    return parts[-1]


def _submit_unibind_form(
    *,
    endpoint_url: str,
    fields: dict[str, str],
    files: dict[str, tuple[str, bytes, str]],
    timeout_seconds: float,
    query: dict[str, Any],
) -> UniBindJobResult:
    cookie_jar = http.cookiejar.CookieJar()
    opener = request.build_opener(request.HTTPCookieProcessor(cookie_jar))

    get_req = request.Request(
        endpoint_url,
        method="GET",
        headers={"User-Agent": "dimelo-toolkit/1.0"},
    )
    with opener.open(get_req, timeout=timeout_seconds) as response:
        landing_html = response.read().decode("utf-8", errors="replace")

    csrf_token = _extract_csrf_token(landing_html)
    fields_with_csrf = dict(fields)
    fields_with_csrf["csrfmiddlewaretoken"] = csrf_token
    body, boundary = _build_multipart_form(fields=fields_with_csrf, files=files)

    post_req = request.Request(
        endpoint_url,
        data=body,
        method="POST",
        headers={
            "User-Agent": "dimelo-toolkit/1.0",
            "Content-Type": f"multipart/form-data; boundary={boundary}",
            "Referer": endpoint_url,
        },
    )
    with opener.open(post_req, timeout=timeout_seconds) as response:
        final_url = str(response.geturl())
        job_html = response.read().decode("utf-8", errors="replace")

    results_url = _extract_unibind_results_url(job_html, base_url=final_url)
    job_url = results_url or final_url
    status = _extract_unibind_job_status(job_html)
    output_urls = _extract_unibind_output_urls(job_html, base_url=job_url)
    job_id = _extract_unibind_job_id(job_url)

    return UniBindJobResult(
        job_id=job_id,
        status=status,
        job_url=job_url,
        endpoint_url=endpoint_url,
        results_url=results_url,
        download_urls=output_urls,
        query=query,
        status_history=[
            {
                "status": status,
                "timestamp": time.time(),
                "stage": "submitted",
                "url": job_url,
            }
        ],
        metadata={},
    )


def _fetch_bytes(url: str, *, timeout_seconds: float = 60.0) -> bytes:
    req = request.Request(url, headers={"User-Agent": "dimelo-toolkit/1.0"})
    try:
        with request.urlopen(req, timeout=timeout_seconds) as response:
            return response.read()
    except error.HTTPError as exc:
        raise RuntimeError(
            f"Track hub request failed with HTTP {exc.code}: {url}"
        ) from exc
    except error.URLError as exc:
        raise RuntimeError(f"Track hub request failed for {url}: {exc.reason}") from exc


def _fetch_text(url: str, *, timeout_seconds: float = 60.0) -> str:
    return _fetch_bytes(url, timeout_seconds=timeout_seconds).decode("utf-8", errors="replace")


def _parse_kv_blocks(text: str) -> list[dict[str, str]]:
    blocks: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            if current:
                blocks.append(current)
                current = {}
            continue
        parts = line.split(None, 1)
        key = parts[0]
        value = parts[1].strip() if len(parts) > 1 else ""
        current[key] = value
    if current:
        blocks.append(current)
    return blocks


def _load_trackdb_stanzas(
    url: str,
    *,
    timeout_seconds: float = 60.0,
    _visited: set[str] | None = None,
) -> list[dict[str, str]]:
    visited = _visited if _visited is not None else set()
    if url in visited:
        return []
    visited.add(url)

    stanzas: list[dict[str, str]] = []
    current: dict[str, str] = {}
    text = _fetch_text(url, timeout_seconds=timeout_seconds)

    def flush_current() -> None:
        nonlocal current
        if current:
            stanzas.append(current)
            current = {}

    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            flush_current()
            continue

        if line.lower().startswith("include "):
            flush_current()
            include_ref = line.split(None, 1)[1].strip().split()[0]
            include_url = parse.urljoin(url, include_ref)
            stanzas.extend(
                _load_trackdb_stanzas(
                    include_url,
                    timeout_seconds=timeout_seconds,
                    _visited=visited,
                )
            )
            continue

        parts = line.split(None, 1)
        key = parts[0]
        value = parts[1].strip() if len(parts) > 1 else ""
        current[key] = value

    flush_current()
    return stanzas


def _select_trackdb_url(
    *,
    hub_url: str,
    assembly: str | None,
    timeout_seconds: float,
) -> tuple[str, str]:
    hub_text = _fetch_text(hub_url, timeout_seconds=timeout_seconds)
    hub_meta = _parse_kv_blocks(hub_text)
    if not hub_meta:
        raise RuntimeError(f"Track hub descriptor is empty: {hub_url}")

    genomes_file = hub_meta[0].get("genomesFile")
    if not genomes_file:
        raise RuntimeError(
            "Track hub descriptor did not include 'genomesFile'. "
            "Only multi-file hubs are currently supported."
        )

    genomes_url = parse.urljoin(hub_url, genomes_file)
    genome_blocks = _parse_kv_blocks(_fetch_text(genomes_url, timeout_seconds=timeout_seconds))
    candidates: list[tuple[str, str]] = []
    for block in genome_blocks:
        genome = block.get("genome")
        trackdb = block.get("trackDb")
        if genome and trackdb:
            candidates.append((genome, parse.urljoin(genomes_url, trackdb)))
    if not candidates:
        raise RuntimeError(f"No genome/trackDb entries found in {genomes_url}")

    if assembly is None:
        genome, trackdb_url = candidates[0]
        return genome, trackdb_url

    requested = _normalize_assembly_token(assembly)
    for genome, trackdb_url in candidates:
        if _normalize_assembly_token(genome) == requested:
            return genome, trackdb_url

    available = ", ".join(genome for genome, _ in candidates)
    raise ValueError(
        f"Assembly {assembly!r} was not found in track hub. Available assemblies: {available}"
    )


def _is_bigbed_path(path: Path) -> bool:
    lower_name = path.name.lower()
    return lower_name.endswith(".bb") or lower_name.endswith(".bigbed")


def _download_track_file(
    *,
    url: str,
    cache_dir: Path,
    allow_cached: bool,
    timeout_seconds: float,
) -> Path:
    filename = Path(parse.urlparse(url).path).name
    if not filename:
        filename = "track_file"
    output = cache_dir / filename

    if allow_cached and output.exists() and output.stat().st_size > 0:
        return output

    payload = _fetch_bytes(url, timeout_seconds=timeout_seconds)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = output.with_suffix(output.suffix + ".tmp")
    tmp_output.write_bytes(payload)
    tmp_output.replace(output)
    return output


def _convert_bigbed_to_bed(path: Path) -> Path:
    converter = shutil.which("bigBedToBed")
    if converter is None:
        raise RuntimeError(
            "Requested bigBed conversion but 'bigBedToBed' was not found on PATH."
        )

    converted = path.with_suffix(".bed")
    proc = subprocess.run(
        [converter, str(path), str(converted)],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"bigBedToBed failed for {path}: {proc.stderr.strip()[:400]}"
        )
    return converted


def unibind_trackhub_url(collection: str = "robust") -> str:
    key = _normalize_token(collection)
    if key not in _UNIBIND_TRACKHUB_URLS:
        available = ", ".join(sorted(_UNIBIND_TRACKHUB_URLS))
        raise ValueError(
            f"Unsupported UniBind track hub collection {collection!r}. Available: {available}."
        )
    return _UNIBIND_TRACKHUB_URLS[key]


def search_unibind_trackhub_tracks(
    *,
    trackhub_url: str | None = None,
    collection: str = "robust",
    assembly: str | None = None,
    search_terms: Iterable[str] | None = None,
    timeout_seconds: float = 60.0,
) -> list[dict[str, str]]:
    resolved_hub_url = trackhub_url if trackhub_url is not None else unibind_trackhub_url(collection)
    resolved_assembly, trackdb_url = _select_trackdb_url(
        hub_url=resolved_hub_url,
        assembly=assembly,
        timeout_seconds=timeout_seconds,
    )
    stanzas = _load_trackdb_stanzas(trackdb_url, timeout_seconds=timeout_seconds)

    terms = [str(term).strip().lower() for term in (search_terms or []) if str(term).strip()]

    rows: list[dict[str, str]] = []
    for stanza in stanzas:
        big_data = stanza.get("bigDataUrl")
        if not big_data:
            continue

        track_name = stanza.get("track", "")
        short_label = stanza.get("shortLabel", "")
        long_label = stanza.get("longLabel", "")
        track_type = stanza.get("type", "")
        big_data_url = parse.urljoin(trackdb_url, big_data)

        if terms:
            haystack = " ".join(
                [track_name, short_label, long_label, big_data_url]
            ).lower()
            if not any(term in haystack for term in terms):
                continue

        rows.append(
            {
                "assembly": resolved_assembly,
                "track": track_name,
                "short_label": short_label,
                "long_label": long_label,
                "type": track_type,
                "url": big_data_url,
            }
        )

    return rows


def resolve_unibind_track_paths(
    *,
    track_paths: Iterable[str | Path] | None = None,
    trackhub_url: str | None = None,
    collection: str = "robust",
    assembly: str | None = None,
    search_terms: Iterable[str] | None = None,
    cache_dir: str | Path = "cache/unibind_tracks",
    max_tracks: int | None = None,
    allow_cached: bool = True,
    timeout_seconds: float = 60.0,
    convert_bigbed_to_bed: bool = False,
) -> list[Path]:
    if track_paths is not None:
        resolved_paths: list[Path] = []
        for raw in track_paths:
            candidate = Path(raw).expanduser().resolve()
            if not candidate.exists():
                raise FileNotFoundError(f"UniBind track path does not exist: {candidate}")
            resolved_paths.append(candidate)
        if not resolved_paths:
            raise ValueError("track_paths was provided but empty.")
        return resolved_paths

    discovered = search_unibind_trackhub_tracks(
        trackhub_url=trackhub_url,
        collection=collection,
        assembly=assembly,
        search_terms=search_terms,
        timeout_seconds=timeout_seconds,
    )
    if not discovered:
        raise ValueError(
            "No UniBind tracks matched the requested track hub search criteria."
        )

    if max_tracks is not None:
        discovered = discovered[: int(max_tracks)]

    cache_root = Path(cache_dir)
    cache_root.mkdir(parents=True, exist_ok=True)

    outputs: list[Path] = []
    for row in discovered:
        downloaded = _download_track_file(
            url=row["url"],
            cache_dir=cache_root,
            allow_cached=allow_cached,
            timeout_seconds=timeout_seconds,
        )
        if convert_bigbed_to_bed and _is_bigbed_path(downloaded):
            downloaded = _convert_bigbed_to_bed(downloaded)
        outputs.append(downloaded)

    return outputs


def submit_unibind_tfbs_extraction(
    *,
    regions: Any,
    species: str = _DEFAULT_SPECIES,
    collection: str = "robust",
    tf_list: str | Path | Iterable[str] | None = None,
    experiment_ids: str | Path | Iterable[str] | None = None,
    name: str | None = None,
    email: str | None = None,
    endpoint_url: str = DEFAULT_UNIBIND_TFBS_EXTRACTION_URL,
    timeout_seconds: float = 120.0,
) -> UniBindJobResult:
    species_code = _normalize_unibind_species_code(species)
    resolved_collection = _normalize_unibind_collection(collection)
    bed_payload = _regions_to_bed_bytes(regions)

    files: dict[str, tuple[str, bytes, str]] = {
        "bed_file": ("regions.bed", bed_payload, "text/plain"),
    }
    tf_payload = _coerce_optional_text_file_payload(
        tf_list,
        default_filename="tf_list.txt",
    )
    if tf_payload is not None:
        files["TFs_file"] = tf_payload
    experiment_payload = _coerce_optional_text_file_payload(
        experiment_ids,
        default_filename="experiment_ids.txt",
    )
    if experiment_payload is not None:
        files["experiments_file"] = experiment_payload

    fields = {
        "name": (name or "").strip(),
        "email": (email or "").strip(),
        "species": species_code,
        "collection": resolved_collection,
        "performTFBSextraction": "1",
    }
    return _submit_unibind_form(
        endpoint_url=endpoint_url,
        fields=fields,
        files=files,
        timeout_seconds=timeout_seconds,
        query={
            "species": _UNIBIND_FORM_SPECIES_BY_CODE.get(species_code, species_code),
            "species_code": species_code,
            "collection": resolved_collection,
            "name": fields["name"],
            "email": fields["email"],
            "has_tf_filter": tf_payload is not None,
            "has_experiment_filter": experiment_payload is not None,
        },
    )


def submit_unibind_enrichment(
    *,
    regions: Any,
    analysis_type: str = "oneSetBg",
    background_regions: Any | None = None,
    comparison_regions: Any | None = None,
    species: str = _DEFAULT_SPECIES,
    collection: str = "robust",
    name: str | None = None,
    email: str | None = None,
    endpoint_url: str = DEFAULT_UNIBIND_ENRICHMENT_URL,
    timeout_seconds: float = 120.0,
) -> UniBindJobResult:
    resolved_analysis_type = _normalize_unibind_enrichment_analysis_type(analysis_type)
    species_code = _normalize_unibind_species_code(species)
    resolved_collection = _normalize_unibind_collection(collection)
    bed_payload = _regions_to_bed_bytes(regions)

    files: dict[str, tuple[str, bytes, str]] = {
        "bed_file_1": ("regions_1.bed", bed_payload, "text/plain"),
    }
    if resolved_analysis_type == "oneSetBg":
        if background_regions is None:
            raise ValueError(
                "UniBind enrichment with analysis_type='oneSetBg' requires "
                "background_regions."
            )
        files["bed_file_background"] = (
            "background_regions.bed",
            _regions_to_bed_bytes(background_regions),
            "text/plain",
        )
    elif resolved_analysis_type == "twoSets":
        if comparison_regions is None:
            raise ValueError(
                "UniBind enrichment with analysis_type='twoSets' requires "
                "comparison_regions."
            )
        files["bed_file_2"] = (
            "regions_2.bed",
            _regions_to_bed_bytes(comparison_regions),
            "text/plain",
        )

    fields = {
        "analysis_type": resolved_analysis_type,
        "name": (name or "").strip(),
        "email": (email or "").strip(),
        "species": species_code,
        "collection": resolved_collection,
        "performEnrichment": "1",
    }
    return _submit_unibind_form(
        endpoint_url=endpoint_url,
        fields=fields,
        files=files,
        timeout_seconds=timeout_seconds,
        query={
            "analysis_type": resolved_analysis_type,
            "species": _UNIBIND_FORM_SPECIES_BY_CODE.get(species_code, species_code),
            "species_code": species_code,
            "collection": resolved_collection,
            "name": fields["name"],
            "email": fields["email"],
            "has_background": background_regions is not None,
            "has_comparison": comparison_regions is not None,
        },
    )


def poll_unibind_job(
    job: UniBindJobResult | str,
    *,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 1200.0,
    request_timeout_seconds: float = 60.0,
) -> UniBindJobResult:
    if isinstance(job, UniBindJobResult):
        current = job
    else:
        parsed = str(job).strip()
        if not parsed:
            raise ValueError("job URL cannot be empty.")
        current = UniBindJobResult(
            job_id=_extract_unibind_job_id(parsed),
            status="unknown",
            job_url=parsed,
            endpoint_url=parse.urljoin(parsed, "."),
            query={},
            status_history=[],
            metadata={},
        )

    history = list(current.status_history)
    started = time.monotonic()
    status = current.status
    results_url = current.results_url
    output_urls = list(current.download_urls)

    while True:
        page_html = _fetch_text(current.job_url, timeout_seconds=request_timeout_seconds)
        status = _extract_unibind_job_status(page_html)
        extracted_results_url = _extract_unibind_results_url(page_html, base_url=current.job_url)
        results_url = extracted_results_url or results_url or current.job_url
        output_urls = _extract_unibind_output_urls(page_html, base_url=current.job_url)

        history.append(
            {
                "status": status,
                "timestamp": time.time(),
                "stage": "poll",
                "url": current.job_url,
            }
        )
        if status in _UNIBIND_TERMINAL_STATUSES:
            break

        elapsed = time.monotonic() - started
        if elapsed > timeout_seconds:
            raise TimeoutError(
                "Timed out waiting for UniBind job completion: "
                f"{current.job_url} (last status={status!r})."
            )
        time.sleep(max(0.1, float(poll_interval_seconds)))

    return UniBindJobResult(
        job_id=current.job_id,
        status=status,
        job_url=current.job_url,
        endpoint_url=current.endpoint_url,
        results_url=results_url,
        download_urls=output_urls,
        query=dict(current.query),
        status_history=history,
        metadata=dict(current.metadata),
    )


def download_unibind_job_outputs(
    *,
    job: UniBindJobResult,
    output_dir: str | Path = "cache/unibind_jobs",
    allow_cached: bool = True,
    timeout_seconds: float = 120.0,
) -> list[Path]:
    if not job.download_urls:
        return []

    root = Path(output_dir).expanduser().resolve() / job.job_id
    root.mkdir(parents=True, exist_ok=True)
    downloaded: list[Path] = []

    for index, url in enumerate(job.download_urls):
        name = Path(parse.urlparse(url).path).name or f"output_{index}"
        destination = root / name
        if destination.exists() and destination.stat().st_size > 0 and allow_cached:
            downloaded.append(destination)
            continue

        payload = _fetch_bytes(url, timeout_seconds=timeout_seconds)
        tmp_path = destination.with_suffix(destination.suffix + ".tmp")
        tmp_path.write_bytes(payload)
        tmp_path.replace(destination)
        downloaded.append(destination)
    return downloaded


def run_unibind_tfbs_extraction(
    *,
    regions: Any,
    species: str = _DEFAULT_SPECIES,
    collection: str = "robust",
    tf_list: str | Path | Iterable[str] | None = None,
    experiment_ids: str | Path | Iterable[str] | None = None,
    name: str | None = None,
    email: str | None = None,
    endpoint_url: str = DEFAULT_UNIBIND_TFBS_EXTRACTION_URL,
    submit_timeout_seconds: float = 120.0,
    wait: bool = True,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 1200.0,
    download_outputs: bool = False,
    output_dir: str | Path = "cache/unibind_jobs",
    allow_cached_downloads: bool = True,
    download_timeout_seconds: float = 120.0,
) -> UniBindJobResult:
    result = submit_unibind_tfbs_extraction(
        regions=regions,
        species=species,
        collection=collection,
        tf_list=tf_list,
        experiment_ids=experiment_ids,
        name=name,
        email=email,
        endpoint_url=endpoint_url,
        timeout_seconds=submit_timeout_seconds,
    )
    if wait:
        result = poll_unibind_job(
            result,
            poll_interval_seconds=poll_interval_seconds,
            timeout_seconds=timeout_seconds,
        )
    if download_outputs:
        if not wait:
            raise ValueError("download_outputs=True requires wait=True.")
        downloaded = download_unibind_job_outputs(
            job=result,
            output_dir=output_dir,
            allow_cached=allow_cached_downloads,
            timeout_seconds=download_timeout_seconds,
        )
        result.metadata["downloaded_outputs"] = [str(path) for path in downloaded]
    return result


def run_unibind_enrichment(
    *,
    regions: Any,
    analysis_type: str = "oneSetBg",
    background_regions: Any | None = None,
    comparison_regions: Any | None = None,
    species: str = _DEFAULT_SPECIES,
    collection: str = "robust",
    name: str | None = None,
    email: str | None = None,
    endpoint_url: str = DEFAULT_UNIBIND_ENRICHMENT_URL,
    submit_timeout_seconds: float = 120.0,
    wait: bool = True,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 1200.0,
    download_outputs: bool = False,
    output_dir: str | Path = "cache/unibind_jobs",
    allow_cached_downloads: bool = True,
    download_timeout_seconds: float = 120.0,
) -> UniBindJobResult:
    result = submit_unibind_enrichment(
        regions=regions,
        analysis_type=analysis_type,
        background_regions=background_regions,
        comparison_regions=comparison_regions,
        species=species,
        collection=collection,
        name=name,
        email=email,
        endpoint_url=endpoint_url,
        timeout_seconds=submit_timeout_seconds,
    )
    if wait:
        result = poll_unibind_job(
            result,
            poll_interval_seconds=poll_interval_seconds,
            timeout_seconds=timeout_seconds,
        )
    if download_outputs:
        if not wait:
            raise ValueError("download_outputs=True requires wait=True.")
        downloaded = download_unibind_job_outputs(
            job=result,
            output_dir=output_dir,
            allow_cached=allow_cached_downloads,
            timeout_seconds=download_timeout_seconds,
        )
        result.metadata["downloaded_outputs"] = [str(path) for path in downloaded]
    return result


def unibind_supported_species() -> tuple[str, ...]:
    return _UNIBIND_SUPPORTED_SPECIES


def screen_supported_species() -> tuple[str, ...]:
    return tuple(sorted(_SCREEN_SUPPORTED_SPECIES))


@dataclass
class RegulatoryEnrichmentSpec:
    providers: tuple[str, ...] = ("screen", "unibind")
    species: str | None = None
    reference_genome: str | None = None
    target_genome: str | None = None
    crossmap_chain_file: str | Path | None = None
    crossmap_chain_url: str | None = None
    crossmap_chain_cache_dir: str | Path | None = None
    crossmap_executable: str | None = "CrossMap.py"
    strict_provider_support: bool = False
    enabled_providers: tuple[str, ...] = field(init=False)
    provider_notes: dict[str, str] = field(init=False, default_factory=dict)

    def __post_init__(self) -> None:
        requested = self.providers if self.providers else ("screen", "unibind")
        normalized_providers = _unique_preserve_order(_normalize_provider(p) for p in requested)

        inferred = _infer_species_from_genomes(
            reference_genome=self.reference_genome,
            target_genome=self.target_genome,
        )
        resolved_species = self.species if self.species is not None else (inferred or _DEFAULT_SPECIES)
        canonical_species = normalize_species_name(str(resolved_species))

        enabled: list[str] = []
        notes: dict[str, str] = {}
        for provider in normalized_providers:
            if provider == "screen" and canonical_species not in _SCREEN_SUPPORTED_SPECIES:
                message = (
                    "SCREEN disabled for species "
                    f"{canonical_species!r}; SCREEN supports only "
                    "homo_sapiens and mus_musculus."
                )
                if self.strict_provider_support:
                    raise RegulatoryEnrichmentSpecError(message)
                notes[provider] = message
                continue
            enabled.append(provider)

        if not enabled:
            raise RegulatoryEnrichmentSpecError(
                "No providers remain after species-based filtering. "
                f"Requested providers={normalized_providers}, species={canonical_species!r}."
            )

        if self.target_genome is None and canonical_species in _DEFAULT_TARGET_GENOME_BY_SPECIES:
            self.target_genome = _DEFAULT_TARGET_GENOME_BY_SPECIES[canonical_species]

        self.providers = tuple(normalized_providers)
        self.species = canonical_species
        self.enabled_providers = tuple(enabled)
        self.provider_notes = notes

    def as_dict(self) -> dict[str, Any]:
        return {
            "providers": list(self.providers),
            "enabled_providers": list(self.enabled_providers),
            "species": self.species,
            "reference_genome": self.reference_genome,
            "target_genome": self.target_genome,
            "crossmap_chain_file": (
                str(self.crossmap_chain_file)
                if isinstance(self.crossmap_chain_file, Path)
                else self.crossmap_chain_file
            ),
            "crossmap_chain_url": self.crossmap_chain_url,
            "crossmap_chain_cache_dir": (
                str(self.crossmap_chain_cache_dir)
                if isinstance(self.crossmap_chain_cache_dir, Path)
                else self.crossmap_chain_cache_dir
            ),
            "crossmap_executable": self.crossmap_executable,
            "provider_notes": dict(self.provider_notes),
        }


__all__ = [
    "DEFAULT_UNIBIND_ENRICHMENT_URL",
    "DEFAULT_UNIBIND_TFBS_EXTRACTION_URL",
    "RegulatoryEnrichmentSpec",
    "RegulatoryEnrichmentSpecError",
    "download_unibind_job_outputs",
    "normalize_species_name",
    "poll_unibind_job",
    "resolve_unibind_track_paths",
    "run_unibind_enrichment",
    "run_unibind_tfbs_extraction",
    "screen_supported_species",
    "search_unibind_trackhub_tracks",
    "submit_unibind_enrichment",
    "submit_unibind_tfbs_extraction",
    "unibind_supported_species",
    "unibind_trackhub_url",
]
