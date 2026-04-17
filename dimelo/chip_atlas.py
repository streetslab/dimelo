from __future__ import annotations

import csv
import io
import json
import re
import shutil
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path
from typing import Any, Iterable, Mapping
from urllib import error, parse, request

import numpy as np
import pandas as pd

from .models import ChipAtlasEnrichmentResult

DEFAULT_SUBMIT_URL = "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/"
DEFAULT_STATUS_URL = DEFAULT_SUBMIT_URL
DEFAULT_RESULT_URL = DEFAULT_SUBMIT_URL
DEFAULT_CHAIN_CACHE_DIR = Path("cache/chains")
DEFAULT_METADATA_CACHE_DIR = Path("cache/chip_atlas")
DEFAULT_EXPERIMENT_LIST_URL = (
    "https://dbarchive.biosciencedbc.jp/data/chip-atlas/LATEST/chip_atlas_experiment_list.zip"
)

_SUCCESS_STATUSES = {"finished"}
_TERMINAL_STATUSES = {"finished", "error", "cancelled", "unknown"}
_STATUS_ALIASES = {
    "finished": "finished",
    "complete": "finished",
    "completed": "finished",
    "done": "finished",
    "success": "finished",
    "ok": "finished",
    "error": "error",
    "failed": "error",
    "failure": "error",
    "cancelled": "cancelled",
    "canceled": "cancelled",
    "running": "running",
    "preparing": "running",
    "queued": "queued",
    "pending": "queued",
    "submitted": "queued",
}

_GENOME_ALIASES = {
    "chm13": "hs1",
    "t2t-chm13v2.0": "hs1",
    "t2t_chm13v2.0": "hs1",
    "t2t_chm13v2": "hs1",
    "hs1": "hs1",
    "hg38": "hg38",
    "grch38": "hg38",
    "hg19": "hg19",
    "grch37": "hg19",
    "mm10": "mm10",
    "mm39": "mm39",
}

_CHAIN_URL_OVERRIDES = {
    ("hs1", "hg38"): "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg38.over.chain.gz",
    ("hs1", "hg19"): "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/liftOver/hs1ToHg19.over.chain.gz",
}

_BED_THRESHOLD_COLUMN_MAP = {
    "05": "Peak-call (BED) (q < 1E-05)",
    "5": "Peak-call (BED) (q < 1E-05)",
    "q05": "Peak-call (BED) (q < 1E-05)",
    "1e-05": "Peak-call (BED) (q < 1E-05)",
    "10": "Peak-call (BED) (q < 1E-10)",
    "q10": "Peak-call (BED) (q < 1E-10)",
    "1e-10": "Peak-call (BED) (q < 1E-10)",
    "20": "Peak-call (BED) (q < 1E-20)",
    "q20": "Peak-call (BED) (q < 1E-20)",
    "1e-20": "Peak-call (BED) (q < 1E-20)",
}

_STRATIFY_TO_BIN_COUNT = {
    "quartile": 4,
    "quartiles": 4,
    "quintile": 5,
    "quintiles": 5,
    "decile": 10,
    "deciles": 10,
}


def _http_request(
    *,
    url: str,
    method: str = "GET",
    data: bytes | None = None,
    timeout_seconds: float = 60.0,
    headers: Mapping[str, str] | None = None,
) -> tuple[int, str, str]:
    request_headers = {"Accept": "application/json,text/plain,text/html,*/*"}
    if headers:
        request_headers.update(dict(headers))
    if method.upper() == "POST":
        request_headers.setdefault("Content-Type", "application/x-www-form-urlencoded")

    req = request.Request(
        url=url,
        data=data,
        method=method.upper(),
        headers=request_headers,
    )
    try:
        with request.urlopen(req, timeout=timeout_seconds) as response:
            status_code = int(response.getcode())
            content_type = response.headers.get("Content-Type", "")
            body = response.read().decode("utf-8", errors="replace")
            return status_code, content_type, body
    except error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        hint = ""
        lowered = detail.lower()
        if exc.code == 404 and ("chip-atlas: 404" in lowered or "<!doctype html>" in lowered):
            hint = (
                " The endpoint may be outdated. The current API base URL is "
                "'https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/'."
            )
        raise RuntimeError(
            f"ChIP-Atlas API request failed with HTTP {exc.code}: {detail[:400]}{hint}"
        ) from exc
    except error.URLError as exc:
        raise RuntimeError(f"ChIP-Atlas API request failed: {exc.reason}") from exc


def _normalize_genome_name(genome: str) -> str:
    token = str(genome).strip()
    if not token:
        raise ValueError("Genome name cannot be empty.")
    lowered = token.lower()
    return _GENOME_ALIASES.get(lowered, lowered)


def _chain_url_for_pair(source_genome: str, target_genome: str) -> str:
    source = _normalize_genome_name(source_genome)
    target = _normalize_genome_name(target_genome)
    override = _CHAIN_URL_OVERRIDES.get((source, target))
    if override:
        return override
    target_token = target[0].upper() + target[1:] if target else target
    return (
        f"https://hgdownload.soe.ucsc.edu/goldenPath/{source}/liftOver/"
        f"{source}To{target_token}.over.chain.gz"
    )


def _resolve_chain_file(
    *,
    source_genome: str,
    target_genome: str,
    chain_file: str | Path | None = None,
    chain_url: str | None = None,
    chain_cache_dir: str | Path | None = None,
    timeout_seconds: float = 120.0,
) -> Path:
    if chain_file is not None:
        resolved = Path(chain_file).expanduser().resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"CrossMap chain file not found: {resolved}")
        return resolved

    source = _normalize_genome_name(source_genome)
    target = _normalize_genome_name(target_genome)
    cache_root = Path(chain_cache_dir) if chain_cache_dir is not None else DEFAULT_CHAIN_CACHE_DIR
    cache_root.mkdir(parents=True, exist_ok=True)
    cached_path = cache_root / f"{source}_to_{target}.over.chain.gz"
    if cached_path.exists() and cached_path.stat().st_size > 0:
        return cached_path

    url = chain_url if chain_url is not None else _chain_url_for_pair(source, target)
    request_obj = request.Request(url, headers={"User-Agent": "dimelo-toolkit/1.0"})
    try:
        with request.urlopen(request_obj, timeout=timeout_seconds) as response:
            with cached_path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
    except Exception as exc:
        raise RuntimeError(
            "Failed to download CrossMap chain file for "
            f"{source_genome!r} -> {target_genome!r} from {url!r}. "
            "Provide an explicit chain_file path to continue."
        ) from exc
    return cached_path


def cache_chain_files(
    *,
    source_genome: str = "chm13",
    target_genomes: Iterable[str] = ("hg38", "hg19"),
    chain_cache_dir: str | Path | None = None,
    timeout_seconds: float = 120.0,
) -> dict[str, Path]:
    cached: dict[str, Path] = {}
    for target in target_genomes:
        cached[str(target)] = _resolve_chain_file(
            source_genome=source_genome,
            target_genome=str(target),
            chain_cache_dir=chain_cache_dir,
            timeout_seconds=timeout_seconds,
        )
    return cached


def _resolve_crossmap_executable(preferred: str | None = None) -> str:
    candidates = []
    if preferred:
        candidates.append(str(preferred))
    for fallback in ("CrossMap.py", "CrossMap"):
        if fallback not in candidates:
            candidates.append(fallback)
    for candidate in candidates:
        if shutil.which(candidate):
            return candidate
    raise RuntimeError(
        "CrossMap executable was not found on PATH. Install CrossMap and ensure "
        "either 'CrossMap.py' or 'CrossMap' is available."
    )


def _normalize_text_token(value: str | None) -> str:
    if value is None:
        return ""
    return str(value).strip().lower()


def _normalize_bed_threshold(threshold: str = "05") -> str:
    token = _normalize_text_token(threshold)
    if token not in _BED_THRESHOLD_COLUMN_MAP:
        supported = ", ".join(sorted(_BED_THRESHOLD_COLUMN_MAP))
        raise ValueError(
            f"Unsupported bed threshold {threshold!r}. Supported aliases: {supported}."
        )
    return token


def _resolve_stratification_bins(stratify: str | int | None) -> int | None:
    if stratify is None:
        return None
    if isinstance(stratify, int):
        if stratify < 2:
            raise ValueError("stratify must be >= 2 when provided as an integer.")
        return int(stratify)
    token = _normalize_text_token(stratify)
    if token in _STRATIFY_TO_BIN_COUNT:
        return _STRATIFY_TO_BIN_COUNT[token]
    if token.isdigit():
        parsed = int(token)
        if parsed < 2:
            raise ValueError("stratify must be >= 2.")
        return parsed
    allowed = ", ".join(sorted(_STRATIFY_TO_BIN_COUNT))
    raise ValueError(
        "Unsupported stratify value "
        f"{stratify!r}. Use an integer >=2 or one of: {allowed}."
    )


def _download_bytes(
    *,
    url: str,
    timeout_seconds: float = 120.0,
    headers: Mapping[str, str] | None = None,
) -> bytes:
    req_headers = {"User-Agent": "dimelo-toolkit/1.0"}
    if headers:
        req_headers.update(dict(headers))
    req = request.Request(url, headers=req_headers)
    try:
        with request.urlopen(req, timeout=timeout_seconds) as response:
            return response.read()
    except error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(
            f"ChIP-Atlas download failed with HTTP {exc.code}: {url}\n{detail[:300]}"
        ) from exc
    except error.URLError as exc:
        raise RuntimeError(f"ChIP-Atlas download failed for {url}: {exc.reason}") from exc


def _ensure_experiment_list_zip(
    *,
    experiment_list_url: str,
    cache_dir: str | Path,
    allow_cached: bool = True,
    timeout_seconds: float = 120.0,
) -> Path:
    cache_root = Path(cache_dir).expanduser().resolve()
    cache_root.mkdir(parents=True, exist_ok=True)
    output_path = cache_root / "chip_atlas_experiment_list.zip"
    if allow_cached and output_path.exists() and output_path.stat().st_size > 0:
        return output_path
    payload = _download_bytes(url=experiment_list_url, timeout_seconds=timeout_seconds)
    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    tmp_path.write_bytes(payload)
    tmp_path.replace(output_path)
    return output_path


def _iter_experiment_rows(experiment_list_zip: Path) -> Iterable[dict[str, str]]:
    with zipfile.ZipFile(experiment_list_zip) as archive:
        member = next((name for name in archive.namelist() if name.lower().endswith(".csv")), None)
        if member is None:
            raise RuntimeError(
                "ChIP-Atlas experiment-list archive did not contain a CSV file."
            )
        with archive.open(member, "r") as handle:
            reader = csv.DictReader(io.TextIOWrapper(handle, encoding="utf-8", errors="replace"))
            for row in reader:
                yield {str(k): ("" if v is None else str(v)) for k, v in row.items()}


def _score_series_for_bed(table: pd.DataFrame) -> pd.Series:
    if table.shape[1] > 4:
        return pd.to_numeric(table.iloc[:, 4], errors="coerce")
    return pd.Series(np.zeros(table.shape[0]), index=table.index, dtype=float)


def _sorted_bed_by_score(table: pd.DataFrame) -> pd.DataFrame:
    scores = _score_series_for_bed(table).fillna(float("-inf"))
    order = scores.sort_values(ascending=False, kind="stable").index
    return table.loc[order].reset_index(drop=True)


def _coerce_bed6(table: pd.DataFrame) -> pd.DataFrame:
    if table.shape[1] < 3:
        raise ValueError("BED table requires at least 3 columns.")
    out = pd.DataFrame(
        {
            "chrom": table.iloc[:, 0].astype(str),
            "start": pd.to_numeric(table.iloc[:, 1], errors="coerce").fillna(-1).astype(int),
            "end": pd.to_numeric(table.iloc[:, 2], errors="coerce").fillna(-1).astype(int),
            "name": table.iloc[:, 3].astype(str) if table.shape[1] > 3 else ".",
            "score": table.iloc[:, 4] if table.shape[1] > 4 else 0,
            "strand": table.iloc[:, 5].astype(str) if table.shape[1] > 5 else ".",
        }
    )
    out = out.loc[out["end"] > out["start"]].reset_index(drop=True)
    return out


def _extract_request_id(body: str) -> str:
    stripped = body.strip()
    if not stripped:
        raise ValueError("ChIP-Atlas submit response was empty; no request id found.")

    try:
        parsed = json.loads(stripped)
    except json.JSONDecodeError:
        parsed = None

    if isinstance(parsed, dict):
        for key in ("id", "requestId", "request_id"):
            value = parsed.get(key)
            if value:
                return str(value)

    if "?" in stripped and "id=" in stripped:
        parsed_url = parse.urlparse(stripped)
        query = parse.parse_qs(parsed_url.query)
        request_ids = query.get("id")
        if request_ids:
            return str(request_ids[0])

    id_match = re.search(r"\b(?:id|requestId|request_id)\s*[:=]\s*['\"]?([A-Za-z0-9._-]+)", stripped)
    if id_match:
        return id_match.group(1)

    if re.fullmatch(r"[A-Za-z0-9._-]+", stripped):
        return stripped

    raise ValueError(
        "Could not parse ChIP-Atlas request id from submit response. "
        f"Response snippet: {stripped[:160]}"
    )


def _normalize_status(status: str | None) -> str:
    if status is None:
        return "unknown"
    normalized = str(status).strip().lower()
    if not normalized:
        return "unknown"
    return _STATUS_ALIASES.get(normalized, normalized)


def _extract_status_value(body: str, parsed_body: Any) -> str | None:
    if isinstance(parsed_body, dict):
        for key in ("status", "state", "job_status", "phase"):
            value = parsed_body.get(key)
            if value is not None:
                return str(value)

    if isinstance(parsed_body, str) and parsed_body.strip():
        return parsed_body.strip().splitlines()[0].strip()

    lowered = body.lower()
    for token in sorted(_STATUS_ALIASES, key=len, reverse=True):
        if re.search(rf"\b{re.escape(token)}\b", lowered):
            return token
    return None


def _parse_region_id(region_id: str) -> tuple[str, int, int, str]:
    raw = region_id.strip()
    if not raw:
        raise ValueError("Empty region id is not allowed.")
    if ":" not in raw or "-" not in raw:
        raise ValueError(
            "Region ids must look like 'chr:start-end' or 'chr:start-end,strand'. "
            f"Received: {region_id!r}"
        )

    chrom, remainder = raw.split(":", 1)
    strand = "."
    if "," in remainder:
        core, strand_token = remainder.rsplit(",", 1)
        if strand_token in {"+", "-", "."}:
            strand = strand_token
        remainder = core
    elif remainder.count(":") >= 1:
        core, suffix = remainder.rsplit(":", 1)
        if suffix in {"+", "-", "."}:
            strand = suffix
            remainder = core

    start_text, end_text = remainder.split("-", 1)
    start = int(start_text)
    end = int(end_text)
    if start < 0 or end < 0 or end <= start:
        raise ValueError(f"Invalid region coordinates in {region_id!r}")
    return chrom, start, end, strand


def region_ids_to_bed_dataframe(region_ids: Iterable[str]) -> pd.DataFrame:
    rows = []
    for index, region_id in enumerate(region_ids):
        chrom, start, end, strand = _parse_region_id(str(region_id))
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": f"region_{index}",
                "score": 0,
                "strand": strand,
            }
        )
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "score", "strand"])


def _regions_dataframe_from_input(
    regions: pd.DataFrame | str | Path | Iterable[str],
) -> pd.DataFrame:
    if isinstance(regions, pd.DataFrame):
        source = regions.copy()
    elif isinstance(regions, Path):
        source = pd.read_csv(regions, sep="\t", header=None, comment="#")
    elif isinstance(regions, str):
        candidate = Path(regions)
        if "\n" not in regions and "\t" not in regions and candidate.exists():
            source = pd.read_csv(candidate, sep="\t", header=None, comment="#")
        else:
            lines = [line.strip() for line in regions.splitlines() if line.strip()]
            return _regions_dataframe_from_lines(lines)
    else:
        lines = [str(line).strip() for line in regions if str(line).strip()]
        return _regions_dataframe_from_lines(lines)

    if not isinstance(source, pd.DataFrame):
        raise TypeError("Could not coerce ChIP-Atlas regions to a DataFrame.")
    if source.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "name", "score", "strand"])

    if source.columns.dtype == object and {"start", "end"}.issubset(source.columns):
        chrom_column = "chrom" if "chrom" in source.columns else "chromosome"
        if chrom_column not in source.columns:
            raise ValueError(
                "Region DataFrame must include 'start', 'end', and 'chrom' or 'chromosome'."
            )
        normalized = pd.DataFrame(
            {
                "chrom": source[chrom_column].astype(str),
                "start": source["start"].astype(int),
                "end": source["end"].astype(int),
                "name": source.get("name", pd.Series(".", index=source.index)).astype(str),
                "score": source.get("score", pd.Series(0, index=source.index)),
                "strand": source.get("strand", pd.Series(".", index=source.index)).astype(str),
            }
        )
        return normalized

    if len(source.columns) < 3:
        raise ValueError("BED-like input requires at least 3 columns.")
    normalized = pd.DataFrame(
        {
            "chrom": source.iloc[:, 0].astype(str),
            "start": source.iloc[:, 1].astype(int),
            "end": source.iloc[:, 2].astype(int),
            "name": source.iloc[:, 3].astype(str) if len(source.columns) > 3 else ".",
            "score": source.iloc[:, 4] if len(source.columns) > 4 else 0,
            "strand": source.iloc[:, 5].astype(str) if len(source.columns) > 5 else ".",
        }
    )
    return normalized


def _regions_dataframe_from_lines(lines: Iterable[str]) -> pd.DataFrame:
    rows = []
    for index, line in enumerate(lines):
        if "\t" in line:
            tokens = line.split("\t")
            chrom = tokens[0]
            start = int(tokens[1])
            end = int(tokens[2])
            name = tokens[3] if len(tokens) > 3 else f"region_{index}"
            score = tokens[4] if len(tokens) > 4 else 0
            strand = tokens[5] if len(tokens) > 5 else "."
        else:
            chrom, start, end, strand = _parse_region_id(line)
            name = f"region_{index}"
            score = 0
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": name,
                "score": score,
                "strand": strand if strand in {"+", "-", "."} else ".",
            }
        )
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "name", "score", "strand"])


def _regions_to_bed_text(regions: pd.DataFrame | str | Path | Iterable[str]) -> str:
    region_frame = _regions_dataframe_from_input(regions)
    if region_frame.empty:
        raise ValueError("ChIP-Atlas enrichment requires at least one region.")
    lines = []
    for row in region_frame.itertuples(index=False):
        lines.append(
            "\t".join(
                [
                    str(row.chrom),
                    str(int(row.start)),
                    str(int(row.end)),
                    str(row.name),
                    str(row.score),
                    str(row.strand if row.strand in {"+", "-", "."} else "."),
                ]
            )
        )
    return "\n".join(lines)


def convert_regions_with_crossmap(
    *,
    regions: pd.DataFrame | str | Path | Iterable[str],
    source_genome: str,
    target_genome: str,
    chain_file: str | Path | None = None,
    chain_url: str | None = None,
    chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
    timeout_seconds: float = 300.0,
) -> pd.DataFrame:
    source = _normalize_genome_name(source_genome)
    target = _normalize_genome_name(target_genome)
    region_frame = _regions_dataframe_from_input(regions)
    if source == target:
        return region_frame
    if region_frame.empty:
        return region_frame

    resolved_chain = _resolve_chain_file(
        source_genome=source,
        target_genome=target,
        chain_file=chain_file,
        chain_url=chain_url,
        chain_cache_dir=chain_cache_dir,
        timeout_seconds=min(timeout_seconds, 120.0),
    )
    with tempfile.TemporaryDirectory(prefix="dimelo-crossmap-") as tmpdir:
        tmp_path = Path(tmpdir)
        input_bed = tmp_path / "input.bed"
        output_bed = tmp_path / "mapped.bed"
        region_frame.to_csv(
            input_bed,
            sep="\t",
            index=False,
            header=False,
            columns=["chrom", "start", "end", "name", "score", "strand"],
        )
        resolved_executable = _resolve_crossmap_executable(crossmap_executable)
        cmd = [
            resolved_executable,
            "bed",
            str(resolved_chain),
            str(input_bed),
            str(output_bed),
        ]
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
                check=False,
            )
        except FileNotFoundError as exc:
            raise RuntimeError(
                f"CrossMap executable {resolved_executable!r} was not found. "
                "Install CrossMap and ensure it is available on PATH."
            ) from exc
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError(
                "CrossMap conversion timed out while converting regions from "
                f"{source_genome!r} to {target_genome!r}."
            ) from exc
        if proc.returncode != 0:
            raise RuntimeError(
                "CrossMap conversion failed with exit code "
                f"{proc.returncode}. stderr: {proc.stderr.strip()[:500]}"
            )
        if not output_bed.exists() or output_bed.stat().st_size == 0:
            raise RuntimeError(
                "CrossMap conversion produced no mapped regions. "
                f"source_genome={source_genome!r}, target_genome={target_genome!r}"
            )

        mapped = pd.read_csv(output_bed, sep="\t", header=None, comment="#")
        if mapped.empty:
            raise RuntimeError(
                "CrossMap conversion returned an empty mapped region table."
            )
        normalized = pd.DataFrame(
            {
                "chrom": mapped.iloc[:, 0].astype(str),
                "start": mapped.iloc[:, 1].astype(int),
                "end": mapped.iloc[:, 2].astype(int),
                "name": mapped.iloc[:, 3].astype(str) if len(mapped.columns) > 3 else ".",
                "score": mapped.iloc[:, 4] if len(mapped.columns) > 4 else 0,
                "strand": mapped.iloc[:, 5].astype(str) if len(mapped.columns) > 5 else ".",
            }
        )
    normalized = normalized[normalized["end"] > normalized["start"]].reset_index(drop=True)
    if normalized.empty:
        raise RuntimeError(
            "CrossMap conversion filtered all regions after coordinate validation."
        )
    return normalized


def _build_url(base_url: str, params_dict: Mapping[str, Any]) -> str:
    query = parse.urlencode({key: value for key, value in params_dict.items() if value is not None})
    return f"{base_url}?{query}" if query else base_url


def _request_url_with_request_id(base_url: str, request_id: str) -> str:
    base = str(base_url).strip()
    if not base:
        raise ValueError("ChIP-Atlas URL cannot be empty.")
    rid = parse.quote(str(request_id).strip(), safe="")
    if not rid:
        raise ValueError("request_id cannot be empty.")
    if "{request_id}" in base:
        return base.format(request_id=rid)
    if base.rstrip("/").endswith(str(request_id).strip()):
        return base
    return f"{base.rstrip('/')}/{rid}"


def submit_enrichment(
    *,
    regions: pd.DataFrame | str | Path | Iterable[str],
    genome: str = "hg38",
    regions_genome: str | None = None,
    antigen_class: str | None = "TFs and others",
    antigen: str | None = None,
    cell_type_class: str | None = "No description",
    cell_type: str | None = None,
    distance: int | None = None,
    threshold: str | None = "100",
    crossmap_chain_file: str | Path | None = None,
    crossmap_chain_url: str | None = None,
    crossmap_chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
    params: Mapping[str, Any] | None = None,
    submit_url: str = DEFAULT_SUBMIT_URL,
    timeout_seconds: float = 60.0,
) -> dict[str, Any]:
    region_payload: pd.DataFrame | str | Path | Iterable[str] = regions
    conversion_metadata: dict[str, Any] = {"applied": False}
    if regions_genome is not None:
        input_build = _normalize_genome_name(regions_genome)
        atlas_build = _normalize_genome_name(genome)
        if input_build != atlas_build:
            resolved_chain = _resolve_chain_file(
                source_genome=input_build,
                target_genome=atlas_build,
                chain_file=crossmap_chain_file,
                chain_url=crossmap_chain_url,
                chain_cache_dir=crossmap_chain_cache_dir,
                timeout_seconds=max(timeout_seconds, 120.0),
            )
            region_payload = convert_regions_with_crossmap(
                regions=regions,
                source_genome=input_build,
                target_genome=atlas_build,
                chain_file=resolved_chain,
                chain_url=crossmap_chain_url,
                chain_cache_dir=crossmap_chain_cache_dir,
                crossmap_executable=crossmap_executable,
                timeout_seconds=max(timeout_seconds, 120.0),
            )
            conversion_metadata = {
                "applied": True,
                "source_genome": input_build,
                "target_genome": atlas_build,
                "chain_file": str(resolved_chain),
            }
        else:
            conversion_metadata = {
                "applied": False,
                "source_genome": input_build,
                "target_genome": atlas_build,
            }

    bed_text = _regions_to_bed_text(region_payload)
    payload: dict[str, Any] = {
        "format": "text",
        "result": "www",
        "genome": genome,
        "typeA": "bed",
        "bedAFile": bed_text,
        "typeB": "rnd",
        "bedBFile": "empty",
        "permTime": 10,
        "descriptionA": "query_regions",
        "descriptionB": "random_permutation_background",
        "title": "dimelo_toolkit_enrichment",
    }
    if antigen_class is not None:
        payload["antigenClass"] = antigen_class
    if antigen is not None:
        payload["antigen"] = antigen
    if cell_type_class is not None:
        payload["cellClass"] = cell_type_class
        payload["celltypeClass"] = cell_type_class
    if cell_type is not None:
        payload["cellType"] = cell_type
        payload["celltype"] = cell_type
    if distance is not None:
        payload["distanceUp"] = int(distance)
        payload["distanceDown"] = int(distance)
    if threshold is not None:
        payload["threshold"] = threshold
    if params:
        payload.update(dict(params))
    encoded = parse.urlencode(
        {
            key: ",".join(map(str, value)) if isinstance(value, (list, tuple, set)) else str(value)
            for key, value in payload.items()
            if value is not None
        }
    ).encode("utf-8")

    status_code, content_type, body = _http_request(
        url=submit_url,
        method="POST",
        data=encoded,
        timeout_seconds=timeout_seconds,
    )
    request_id = _extract_request_id(body)
    return {
        "request_id": request_id,
        "status_code": status_code,
        "content_type": content_type,
        "response": body,
        "submit_url": submit_url,
        "query": {
            "genome": genome,
            "regions_genome": regions_genome,
            "antigen_class": antigen_class,
            "antigen": antigen,
            "cell_type_class": cell_type_class,
            "cell_type": cell_type,
            "distance": distance,
            "threshold": threshold,
            "crossmap": conversion_metadata,
            "params": dict(params or {}),
        },
    }


def get_status(
    request_id: str,
    *,
    status_url: str = DEFAULT_STATUS_URL,
    params: Mapping[str, Any] | None = None,
    timeout_seconds: float = 30.0,
) -> dict[str, Any]:
    url = _request_url_with_request_id(status_url, request_id)
    query_params: dict[str, Any] = {}
    if params:
        query_params.update(dict(params))
    if query_params:
        url = _build_url(url, query_params)
    status_code, content_type, body = _http_request(
        url=url,
        method="GET",
        timeout_seconds=timeout_seconds,
    )
    try:
        parsed_body = json.loads(body)
    except json.JSONDecodeError:
        parsed_body = None

    raw_status = _extract_status_value(body, parsed_body)
    status = _normalize_status(raw_status)
    return {
        "request_id": request_id,
        "status": status,
        "raw_status": raw_status,
        "status_code": status_code,
        "content_type": content_type,
        "response": parsed_body if parsed_body is not None else body,
        "status_url": url,
    }


def poll_request(
    request_id: str,
    *,
    status_url: str = DEFAULT_STATUS_URL,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 600.0,
    params: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    history: list[dict[str, Any]] = []
    deadline = time.monotonic() + timeout_seconds
    while True:
        snapshot = get_status(
            request_id,
            status_url=status_url,
            params=params,
            timeout_seconds=min(30.0, max(poll_interval_seconds, 1.0)),
        )
        history.append(
            {
                "time_epoch_s": time.time(),
                "status": snapshot["status"],
                "raw_status": snapshot.get("raw_status"),
                "status_code": snapshot.get("status_code"),
            }
        )
        if snapshot["status"] in _TERMINAL_STATUSES:
            return {
                "request_id": request_id,
                "status": snapshot["status"],
                "history": history,
                "last_snapshot": snapshot,
            }
        if time.monotonic() >= deadline:
            raise TimeoutError(
                f"Timed out waiting for ChIP-Atlas request {request_id!r} after {timeout_seconds} seconds."
            )
        time.sleep(max(poll_interval_seconds, 0.0))


def fetch_result(
    request_id: str,
    *,
    result_url: str = DEFAULT_RESULT_URL,
    params: Mapping[str, Any] | None = None,
    timeout_seconds: float = 120.0,
    as_dataframe: bool = True,
    sep: str = "\t",
) -> pd.DataFrame | str:
    url = _request_url_with_request_id(result_url, request_id)
    query_params: dict[str, Any] = {"info": "result", "format": "tsv"}
    if params:
        query_params.update(dict(params))
    url = _build_url(url, query_params)
    _, _, body = _http_request(
        url=url,
        method="GET",
        timeout_seconds=timeout_seconds,
    )
    if not as_dataframe:
        return body
    return pd.read_csv(io.StringIO(body), sep=sep)


def run_enrichment(
    *,
    regions: pd.DataFrame | str | Path | Iterable[str],
    genome: str = "hg38",
    regions_genome: str | None = None,
    antigen_class: str | None = "TFs and others",
    antigen: str | None = None,
    cell_type_class: str | None = "No description",
    cell_type: str | None = None,
    distance: int | None = None,
    threshold: str | None = "100",
    crossmap_chain_file: str | Path | None = None,
    crossmap_chain_url: str | None = None,
    crossmap_chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
    params: Mapping[str, Any] | None = None,
    wait: bool = True,
    fetch_results: bool = True,
    raise_on_failure: bool = True,
    submit_url: str = DEFAULT_SUBMIT_URL,
    status_url: str = DEFAULT_STATUS_URL,
    result_url: str = DEFAULT_RESULT_URL,
    poll_interval_seconds: float = 5.0,
    timeout_seconds: float = 600.0,
) -> ChipAtlasEnrichmentResult:
    submission = submit_enrichment(
        regions=regions,
        genome=genome,
        regions_genome=regions_genome,
        antigen_class=antigen_class,
        antigen=antigen,
        cell_type_class=cell_type_class,
        cell_type=cell_type,
        distance=distance,
        threshold=threshold,
        crossmap_chain_file=crossmap_chain_file,
        crossmap_chain_url=crossmap_chain_url,
        crossmap_chain_cache_dir=crossmap_chain_cache_dir,
        crossmap_executable=crossmap_executable,
        params=params,
        submit_url=submit_url,
    )
    request_id = str(submission["request_id"])
    history: list[dict[str, Any]] = []
    final_status = "queued"

    if wait:
        poll = poll_request(
            request_id,
            status_url=status_url,
            poll_interval_seconds=poll_interval_seconds,
            timeout_seconds=timeout_seconds,
        )
        final_status = poll["status"]
        history = poll["history"]

    if raise_on_failure and wait and final_status not in _SUCCESS_STATUSES:
        raise RuntimeError(
            f"ChIP-Atlas request {request_id!r} finished with status {final_status!r}."
        )

    result_table: pd.DataFrame | None = None
    if fetch_results and wait and final_status in _SUCCESS_STATUSES:
        fetched = fetch_result(
            request_id,
            result_url=result_url,
            timeout_seconds=timeout_seconds,
            as_dataframe=True,
        )
        result_table = fetched if isinstance(fetched, pd.DataFrame) else None

    return ChipAtlasEnrichmentResult(
        request_id=request_id,
        status=final_status,
        results=result_table,
        query=submission["query"],
        status_history=history,
        submit_url=submit_url,
        status_url=status_url,
        result_url=result_url,
        metadata={"submission": {k: v for k, v in submission.items() if k != "response"}},
    )


def search_peak_datasets(
    *,
    antigen: str,
    genome: str = "hg38",
    cell_type: str | None = None,
    antigen_class: str | None = None,
    cell_type_class: str | None = None,
    threshold: str = "05",
    match_mode: str = "contains",
    max_results: int | None = None,
    experiment_list_url: str = DEFAULT_EXPERIMENT_LIST_URL,
    cache_dir: str | Path = DEFAULT_METADATA_CACHE_DIR,
    allow_cached_metadata: bool = True,
    timeout_seconds: float = 120.0,
) -> pd.DataFrame:
    """
    Search ChIP-Atlas experiment-level peak datasets by antigen and optional metadata filters.

    This reads the published ChIP-Atlas experiment metadata table and returns candidate datasets
    with direct BED URLs for the requested peak threshold.
    """

    antigen_query = str(antigen).strip()
    if not antigen_query:
        raise ValueError("antigen is required.")
    if match_mode not in {"contains", "exact"}:
        raise ValueError("match_mode must be 'contains' or 'exact'.")
    if max_results is not None and int(max_results) <= 0:
        raise ValueError("max_results must be positive when provided.")

    threshold_key = _normalize_bed_threshold(threshold)
    bed_column = _BED_THRESHOLD_COLUMN_MAP[threshold_key]
    experiment_zip = _ensure_experiment_list_zip(
        experiment_list_url=experiment_list_url,
        cache_dir=cache_dir,
        allow_cached=allow_cached_metadata,
        timeout_seconds=timeout_seconds,
    )

    genome_query = _normalize_text_token(genome)
    antigen_class_query = _normalize_text_token(antigen_class)
    cell_type_class_query = _normalize_text_token(cell_type_class)
    cell_type_query = _normalize_text_token(cell_type)
    antigen_norm = _normalize_text_token(antigen_query)

    def _matches(value: str, query: str) -> bool:
        value_norm = _normalize_text_token(value)
        if not query:
            return True
        if match_mode == "exact":
            return value_norm == query
        return query in value_norm

    rows: list[dict[str, Any]] = []
    for raw in _iter_experiment_rows(experiment_zip):
        row_genome = _normalize_text_token(raw.get("Genome assembly", ""))
        if genome_query and row_genome != genome_query:
            continue
        if not _matches(raw.get("Antigen", ""), antigen_norm):
            continue
        if antigen_class_query and not _matches(raw.get("Antigen class", ""), antigen_class_query):
            continue
        if cell_type_class_query and not _matches(raw.get("Cell type class", ""), cell_type_class_query):
            continue
        if cell_type_query and not _matches(raw.get("Cell type", ""), cell_type_query):
            continue

        bed_url = str(raw.get(bed_column, "")).strip()
        if not bed_url:
            continue
        rows.append(
            {
                "dataset_id": str(raw.get("Experimental ID", "")).strip(),
                "genome_assembly": str(raw.get("Genome assembly", "")).strip(),
                "antigen_class": str(raw.get("Antigen class", "")).strip(),
                "antigen": str(raw.get("Antigen", "")).strip(),
                "cell_type_class": str(raw.get("Cell type class", "")).strip(),
                "cell_type": str(raw.get("Cell type", "")).strip(),
                "cell_type_description": str(raw.get("Cell type description", "")).strip(),
                "processing_logs": str(raw.get("Processing logs", "")).strip(),
                "title": str(raw.get("Title", "")).strip(),
                "metadata": str(raw.get("Meta data", "")).strip(),
                "bigwig_url": str(raw.get("BigWig", "")).strip(),
                "bed_url": bed_url,
                "bigbed_url": str(
                    raw.get(
                        bed_column.replace("Peak-call (BED)", "Peak-call (BigBed)"),
                        "",
                    )
                ).strip(),
                "threshold_key": threshold_key,
                "threshold_column": bed_column,
            }
        )
        if max_results is not None and len(rows) >= int(max_results):
            break

    return pd.DataFrame(rows)


def download_peak_datasets(
    *,
    datasets: pd.DataFrame,
    dataset_ids: Iterable[str] | None = None,
    output_dir: str | Path = "cache/chip_atlas/peak_sets",
    include_complete_sorted: bool = True,
    include_top_n: int | None = 3000,
    include_bottom_n: int | None = 3000,
    stratify: str | int | None = None,
    allow_cached: bool = True,
    timeout_seconds: float = 180.0,
    crossmap_target_genome: str | None = None,
    crossmap_chain_file: str | Path | None = None,
    crossmap_chain_url: str | None = None,
    crossmap_chain_cache_dir: str | Path | None = None,
    crossmap_executable: str | None = "CrossMap.py",
) -> pd.DataFrame:
    """
    Download selected ChIP-Atlas BED datasets and optionally generate sorted/stratified subsets.

    Outputs are written as BED files under ``output_dir/<dataset_id>/`` and returned in a manifest.
    """

    if not isinstance(datasets, pd.DataFrame):
        raise TypeError("datasets must be a pandas DataFrame.")
    required = {"dataset_id", "bed_url", "genome_assembly"}
    missing = required - set(datasets.columns)
    if missing:
        raise ValueError(
            "datasets is missing required columns: "
            f"{', '.join(sorted(missing))}"
        )
    if include_top_n is not None and int(include_top_n) < 0:
        raise ValueError("include_top_n must be non-negative.")
    if include_bottom_n is not None and int(include_bottom_n) < 0:
        raise ValueError("include_bottom_n must be non-negative.")

    stratify_bins = _resolve_stratification_bins(stratify)
    selected = datasets.copy()
    if dataset_ids is not None:
        wanted = {str(dataset_id).strip() for dataset_id in dataset_ids if str(dataset_id).strip()}
        selected = selected.loc[selected["dataset_id"].astype(str).isin(wanted)].copy()
    if selected.empty:
        return pd.DataFrame(
            columns=[
                "dataset_id",
                "variant",
                "path",
                "n_rows",
                "source_genome",
                "target_genome",
                "crossmapped",
                "bed_url",
            ]
        )

    output_root = Path(output_dir).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[dict[str, Any]] = []

    def _write_variant(path: Path, frame: pd.DataFrame) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, sep="\t", index=False, header=False)
        return path

    for row in selected.itertuples(index=False):
        dataset_id = str(getattr(row, "dataset_id"))
        bed_url = str(getattr(row, "bed_url"))
        source_genome = str(getattr(row, "genome_assembly"))
        dataset_dir = output_root / dataset_id
        raw_path = dataset_dir / "raw.bed"

        if not (allow_cached and raw_path.exists() and raw_path.stat().st_size > 0):
            payload = _download_bytes(url=bed_url, timeout_seconds=timeout_seconds)
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            tmp_path = raw_path.with_suffix(".bed.tmp")
            tmp_path.write_bytes(payload)
            tmp_path.replace(raw_path)

        bed_table = pd.read_csv(raw_path, sep="\t", header=None, comment="#")
        if bed_table.empty:
            continue
        sorted_table = _sorted_bed_by_score(bed_table)
        variants: list[tuple[str, pd.DataFrame]] = []
        if include_complete_sorted:
            variants.append(("full_sorted", sorted_table))
        if include_top_n is not None and int(include_top_n) > 0:
            variants.append((f"top_{int(include_top_n)}", sorted_table.head(int(include_top_n)).copy()))
        if include_bottom_n is not None and int(include_bottom_n) > 0:
            variants.append((f"bottom_{int(include_bottom_n)}", sorted_table.tail(int(include_bottom_n)).copy()))
        if stratify_bins is not None and len(sorted_table) > 0:
            chunks = np.array_split(np.arange(len(sorted_table)), stratify_bins)
            for i, idx in enumerate(chunks, start=1):
                if len(idx) == 0:
                    continue
                variants.append(
                    (
                        f"quantile_{i}_of_{stratify_bins}",
                        sorted_table.iloc[idx].copy(),
                    )
                )

        target_genome = source_genome
        if crossmap_target_genome is not None:
            target_genome = str(crossmap_target_genome)

        for variant_name, variant_table in variants:
            variant_path = dataset_dir / f"{variant_name}.bed"
            _write_variant(variant_path, variant_table)
            manifest_rows.append(
                {
                    "dataset_id": dataset_id,
                    "variant": variant_name,
                    "path": str(variant_path),
                    "n_rows": int(len(variant_table)),
                    "source_genome": source_genome,
                    "target_genome": source_genome,
                    "crossmapped": False,
                    "bed_url": bed_url,
                }
            )
            if crossmap_target_genome is None:
                continue
            if _normalize_genome_name(source_genome) == _normalize_genome_name(target_genome):
                continue
            bed6 = _coerce_bed6(variant_table)
            if bed6.empty:
                continue
            mapped = convert_regions_with_crossmap(
                regions=bed6,
                source_genome=source_genome,
                target_genome=target_genome,
                chain_file=crossmap_chain_file,
                chain_url=crossmap_chain_url,
                chain_cache_dir=crossmap_chain_cache_dir,
                crossmap_executable=crossmap_executable,
            )
            mapped_path = dataset_dir / f"{variant_name}.{target_genome}.bed"
            mapped.to_csv(mapped_path, sep="\t", index=False, header=False)
            manifest_rows.append(
                {
                    "dataset_id": dataset_id,
                    "variant": f"{variant_name}_crossmapped",
                    "path": str(mapped_path),
                    "n_rows": int(len(mapped)),
                    "source_genome": source_genome,
                    "target_genome": target_genome,
                    "crossmapped": True,
                    "bed_url": bed_url,
                }
            )

    return pd.DataFrame(manifest_rows)


__all__ = [
    "DEFAULT_EXPERIMENT_LIST_URL",
    "DEFAULT_CHAIN_CACHE_DIR",
    "DEFAULT_METADATA_CACHE_DIR",
    "DEFAULT_RESULT_URL",
    "DEFAULT_STATUS_URL",
    "DEFAULT_SUBMIT_URL",
    "cache_chain_files",
    "convert_regions_with_crossmap",
    "download_peak_datasets",
    "fetch_result",
    "get_status",
    "poll_request",
    "region_ids_to_bed_dataframe",
    "run_enrichment",
    "search_peak_datasets",
    "submit_enrichment",
]
