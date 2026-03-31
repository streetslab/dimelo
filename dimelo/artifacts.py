from __future__ import annotations

from collections.abc import Iterable

from dimelo.models import DatasetArtifact

_FINGERPRINT_FIELDS = (
    "sample_id",
    "artifact_type",
    "format",
    "params",
    "provenance",
)


def artifact_fingerprint(artifact: DatasetArtifact) -> dict[str, object]:
    return {field: getattr(artifact, field) for field in _FINGERPRINT_FIELDS}


def artifact_is_compatible(
    requested: DatasetArtifact,
    candidate: DatasetArtifact,
) -> bool:
    return all(
        getattr(requested, field) == getattr(candidate, field)
        for field in _FINGERPRINT_FIELDS
    )


def resolve_artifact(
    requested: DatasetArtifact,
    candidates: Iterable[DatasetArtifact],
    artifact_policy: str = "prefer_cached",
) -> DatasetArtifact:
    if artifact_policy == "prefer_cached":
        for candidate in candidates:
            if artifact_is_compatible(requested, candidate):
                return candidate
    return requested
