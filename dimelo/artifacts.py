from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable

from dimelo.models import DatasetArtifact


def _params_hash(params: dict[str, object]) -> str:
    payload = json.dumps(params, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def artifact_fingerprint(artifact: DatasetArtifact) -> dict[str, object]:
    return {
        "schema": artifact.metadata.get("schema"),
        "package": artifact.metadata.get("package"),
        "source": artifact.provenance.get("source"),
        "params_hash": _params_hash(artifact.params),
    }


def artifact_is_compatible(
    requested: DatasetArtifact,
    candidate: DatasetArtifact,
) -> bool:
    requested_fingerprint = artifact_fingerprint(requested)
    candidate_fingerprint = artifact_fingerprint(candidate)
    if any(
        requested_fingerprint[field] != candidate_fingerprint[field]
        for field in ("schema", "package")
    ):
        return False
    if (
        requested_fingerprint["source"] is not None
        and requested_fingerprint["source"] != candidate_fingerprint["source"]
    ):
        return False

    return all(item in candidate.params.items() for item in requested.params.items()) and all(
        item in candidate.provenance.items() for item in requested.provenance.items()
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
    raise LookupError("No compatible cached artifact found for require_cached")
