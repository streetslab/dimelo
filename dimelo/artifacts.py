from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable

from dimelo.models import DatasetArtifact

COMPATIBILITY_PROVENANCE_KEYS = frozenset(
    {
        "normalization_mode",
        "feature_scaling",
        "cluster_basis",
        "motifs",
        "window_size",
        "region_source",
        "region_digest",
        "pipeline",
    }
)


def _params_hash(params: dict[str, object]) -> str:
    payload = json.dumps(params, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _normalize_sequence(values: object) -> tuple[object, ...]:
    if values is None:
        return None
    if isinstance(values, str):
        return (values,)
    return tuple(values)


def _normalize_source_files(values: object) -> tuple[str, ...] | None:
    normalized = _normalize_sequence(values)
    if normalized is None:
        return None
    return tuple(sorted(str(value) for value in normalized))


def _normalize_source_fingerprints(
    values: object,
) -> tuple[dict[str, object], ...] | None:
    if values is None:
        return None
    normalized = [dict(value) for value in values]
    normalized.sort(
        key=lambda value: json.dumps(value, sort_keys=True, separators=(",", ":"))
    )
    return tuple(normalized)


def _mapping_subset_matches(
    requested: dict[str, object],
    candidate: dict[str, object],
) -> bool:
    return all(item in candidate.items() for item in requested.items())


def _requested_params_hash_matches(
    requested: dict[str, object],
    candidate: dict[str, object],
) -> bool:
    if not _mapping_subset_matches(requested, candidate):
        return False
    candidate_subset = {key: candidate[key] for key in requested}
    return _params_hash(requested) == _params_hash(candidate_subset)


def _compatibility_provenance(requested: dict[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in requested.items()
        if key in COMPATIBILITY_PROVENANCE_KEYS
    }


def artifact_fingerprint(artifact: DatasetArtifact) -> dict[str, object]:
    return {
        "schema_version": artifact.metadata.get("schema_version"),
        "package_version": artifact.metadata.get("package_version"),
        "source_files": _normalize_source_files(
            artifact.provenance.get(
                "source_files", artifact.metadata.get("source_files")
            )
        ),
        "source_fingerprints": _normalize_source_fingerprints(
            artifact.provenance.get(
                "source_fingerprints", artifact.metadata.get("source_fingerprints")
            )
        ),
        "upstream_lineage": _normalize_sequence(
            artifact.provenance.get(
                "upstream_lineage", artifact.metadata.get("upstream_lineage")
            )
        ),
        "params_hash": _params_hash(artifact.params),
    }


def _has_required_fingerprint_fields(fingerprint: dict[str, object]) -> bool:
    required_fields = (
        "schema_version",
        "package_version",
        "source_files",
        "source_fingerprints",
        "upstream_lineage",
    )
    if any(fingerprint[field] is None for field in required_fields):
        return False
    return bool(fingerprint["source_files"]) and bool(
        fingerprint["source_fingerprints"]
    )


def artifact_is_compatible(
    requested: DatasetArtifact,
    candidate: DatasetArtifact,
) -> bool:
    requested_fingerprint = artifact_fingerprint(requested)
    candidate_fingerprint = artifact_fingerprint(candidate)
    if not _has_required_fingerprint_fields(requested_fingerprint):
        return False
    if not _has_required_fingerprint_fields(candidate_fingerprint):
        return False
    if requested.sample_id != candidate.sample_id:
        return False
    if requested.artifact_type != candidate.artifact_type:
        return False
    if any(
        requested_fingerprint[field] != candidate_fingerprint[field]
        for field in (
            "schema_version",
            "package_version",
            "source_files",
            "source_fingerprints",
            "upstream_lineage",
        )
    ):
        return False
    if not _requested_params_hash_matches(requested.params, candidate.params):
        return False
    return _mapping_subset_matches(
        _compatibility_provenance(requested.provenance),
        candidate.provenance,
    )


def resolve_artifact(
    requested: DatasetArtifact,
    candidates: Iterable[DatasetArtifact],
    artifact_policy: str = "prefer_cached",
) -> DatasetArtifact | None:
    if artifact_policy == "rebuild":
        return None

    if artifact_policy not in {"prefer_cached", "require_cached"}:
        raise ValueError(f"Unknown artifact_policy: {artifact_policy}")

    for candidate in candidates:
        if artifact_is_compatible(requested, candidate):
            return candidate

    if artifact_policy == "prefer_cached":
        return None
    raise FileNotFoundError("No compatible cached artifact found for require_cached")
