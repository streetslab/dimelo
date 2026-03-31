import hashlib
import json
from pathlib import Path

import pytest

from dimelo.artifacts import artifact_fingerprint, artifact_is_compatible, resolve_artifact
from dimelo.models import DatasetArtifact


def make_artifact(
    *,
    sample_id: str = "sample-1",
    artifact_type: str = "extract",
    path: str,
    params: dict[str, object],
    provenance: dict[str, object],
    metadata: dict[str, object] | None = None,
) -> DatasetArtifact:
    return DatasetArtifact(
        sample_id=sample_id,
        artifact_type=artifact_type,
        path=Path(path),
        format="hdf5",
        params=params,
        provenance=provenance,
        metadata=metadata or {},
    )


def test_artifact_fingerprint_includes_required_keys():
    artifact = make_artifact(
        path="requested.h5",
        params={"window_size": 200, "min_coverage": 5},
        provenance={
            "source_files": ["reads-a.bam", "reads-b.bam"],
            "source_fingerprints": [
                {"path": "reads-a.bam", "size": 100, "mtime": 10, "digest": "aaa"},
                {"path": "reads-b.bam", "size": 200, "mtime": 20, "digest": "bbb"},
            ],
            "upstream_lineage": ["parse_bam", "normalize"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )

    fingerprint = artifact_fingerprint(artifact)

    assert set(fingerprint) == {
        "schema_version",
        "package_version",
        "source_files",
        "source_fingerprints",
        "upstream_lineage",
        "params_hash",
    }
    assert fingerprint["schema_version"] == "artifact-v1"
    assert fingerprint["package_version"] == "1.2.3"
    assert fingerprint["source_files"] == ("reads-a.bam", "reads-b.bam")
    assert fingerprint["source_fingerprints"] == (
        {"path": "reads-a.bam", "size": 100, "mtime": 10, "digest": "aaa"},
        {"path": "reads-b.bam", "size": 200, "mtime": 20, "digest": "bbb"},
    )
    assert fingerprint["upstream_lineage"] == ("parse_bam", "normalize")
    expected_params_hash = hashlib.sha256(
        json.dumps(
            {"min_coverage": 5, "window_size": 200},
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()
    assert fingerprint["params_hash"] == expected_params_hash


def test_artifact_is_compatible_rejects_sample_id_mismatch():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )
    candidate = make_artifact(
        sample_id="sample-2",
        path="cached.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )

    assert not artifact_is_compatible(requested, candidate)


def test_artifact_is_compatible_rejects_artifact_type_mismatch():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )
    candidate = make_artifact(
        artifact_type="other",
        path="cached.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )

    assert not artifact_is_compatible(requested, candidate)


def test_resolve_artifact_prefers_matching_cached_artifact():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )
    cached_match = make_artifact(
        path="cached.h5",
        params={"window_size": 200, "threads": 8},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )

    resolved = resolve_artifact(
        requested,
        [cached_match],
        artifact_policy="prefer_cached",
    )

    assert resolved is cached_match


def test_resolve_artifact_returns_none_on_prefer_cached_miss():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )
    candidate = make_artifact(
        sample_id="sample-2",
        path="cached.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )

    resolved = resolve_artifact(
        requested,
        [candidate],
        artifact_policy="prefer_cached",
    )

    assert resolved is None


def test_resolve_artifact_raises_on_require_cached_miss():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )

    with pytest.raises(LookupError, match="require_cached"):
        resolve_artifact(
            requested,
            [],
            artifact_policy="require_cached",
        )


def test_resolve_artifact_returns_none_for_rebuild():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )
    candidate = make_artifact(
        path="cached.h5",
        params={"window_size": 200},
        provenance={
            "source_files": ["reads-a.bam"],
            "source_fingerprints": [{"path": "reads-a.bam", "size": 100, "mtime": 10}],
            "upstream_lineage": ["parse_bam"],
        },
        metadata={"schema_version": "artifact-v1", "package_version": "1.2.3"},
    )

    resolved = resolve_artifact(
        requested,
        [candidate],
        artifact_policy="rebuild",
    )

    assert resolved is None
