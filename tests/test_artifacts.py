from pathlib import Path

from dimelo.artifacts import artifact_fingerprint, artifact_is_compatible, resolve_artifact
from dimelo.models import DatasetArtifact


def make_artifact(
    *,
    path: str,
    params: dict[str, object],
    provenance: dict[str, object],
) -> DatasetArtifact:
    return DatasetArtifact(
        sample_id="sample-1",
        artifact_type="extract",
        path=Path(path),
        format="hdf5",
        params=params,
        provenance=provenance,
    )


def test_artifact_fingerprint_includes_required_keys():
    artifact = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
    )

    fingerprint = artifact_fingerprint(artifact)

    assert set(fingerprint) == {
        "sample_id",
        "artifact_type",
        "format",
        "params",
        "provenance",
    }
    assert fingerprint["params"] == {"window_size": 200}
    assert fingerprint["provenance"] == {"pipeline": "parse_bam"}


def test_artifact_is_compatible_rejects_parameter_mismatch():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
    )
    candidate = make_artifact(
        path="cached.h5",
        params={"window_size": 100},
        provenance={"pipeline": "parse_bam"},
    )

    assert not artifact_is_compatible(requested, candidate)


def test_resolve_artifact_prefers_matching_cached_artifact():
    requested = make_artifact(
        path="requested.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
    )
    cached_match = make_artifact(
        path="cached.h5",
        params={"window_size": 200},
        provenance={"pipeline": "parse_bam"},
    )
    cached_mismatch = make_artifact(
        path="other-cached.h5",
        params={"window_size": 100},
        provenance={"pipeline": "parse_bam"},
    )

    resolved = resolve_artifact(
        requested,
        [cached_mismatch, cached_match],
        artifact_policy="prefer_cached",
    )

    assert resolved is cached_match
