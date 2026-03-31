import hashlib
import json
from pathlib import Path

from dimelo.artifacts import artifact_fingerprint, artifact_is_compatible, resolve_artifact
from dimelo.models import DatasetArtifact


def make_artifact(
    *,
    path: str,
    params: dict[str, object],
    provenance: dict[str, object],
    metadata: dict[str, object] | None = None,
) -> DatasetArtifact:
    return DatasetArtifact(
        sample_id="sample-1",
        artifact_type="extract",
        path=Path(path),
        format="hdf5",
        params=params,
        provenance=provenance,
        metadata=metadata or {},
    )


def test_artifact_fingerprint_includes_required_keys():
    artifact = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam", "source": "modkit"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )

    fingerprint = artifact_fingerprint(artifact)

    assert set(fingerprint) == {
        "schema",
        "package",
        "source",
        "params_hash",
    }
    expected_params_hash = hashlib.sha256(
        json.dumps({"window_size": 200}, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert fingerprint == {
        "schema": "artifact-v1",
        "package": "dimelo",
        "source": "modkit",
        "params_hash": expected_params_hash,
    }


def test_artifact_is_compatible_accepts_subset_matches():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )
    candidate = make_artifact(
        path="cached.h5",
        params={"window_size": 200, "threads": 8},
        provenance={"pipeline": "parse_bam", "source": "modkit"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )

    assert artifact_is_compatible(requested, candidate)


def test_artifact_is_compatible_rejects_parameter_mismatch():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )
    candidate = make_artifact(
        path="cached.h5",
        params={"window_size": 100, "threads": 8},
        provenance={"pipeline": "parse_bam", "source": "modkit"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )

    assert not artifact_is_compatible(requested, candidate)


def test_resolve_artifact_prefers_matching_cached_artifact():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )
    cached_match = make_artifact(
        path="cached.h5",
        params={"window_size": 200, "threads": 8},
        provenance={"pipeline": "parse_bam", "source": "modkit"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )
    cached_mismatch = make_artifact(
        path="other-cached.h5",
        params={"window_size": 100},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )

    resolved = resolve_artifact(
        requested,
        [cached_mismatch, cached_match],
        artifact_policy="prefer_cached",
    )

    assert resolved is cached_match


def test_resolve_artifact_returns_none_on_prefer_cached_miss():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )

    resolved = resolve_artifact(
        requested,
        [],
        artifact_policy="prefer_cached",
    )

    assert resolved is None


def test_resolve_artifact_raises_on_require_cached_miss():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )

    try:
        resolve_artifact(
            requested,
            [],
            artifact_policy="require_cached",
        )
    except LookupError as exc:
        assert "require_cached" in str(exc)
    else:
        raise AssertionError("Expected resolve_artifact to raise LookupError")


def test_resolve_artifact_returns_none_for_rebuild():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )
    cached_match = make_artifact(
        path="cached.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
        metadata={"schema": "artifact-v1", "package": "dimelo"},
    )

    resolved = resolve_artifact(
        requested,
        [cached_match],
        artifact_policy="rebuild",
    )

    assert resolved is None
