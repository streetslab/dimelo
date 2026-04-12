from __future__ import annotations

from dimelo.models import ContrastSpec, SharedClusterContrastResult, SharedClusterResult


def _require_supported_shared_cluster_mode(contrast: ContrastSpec) -> None:
    supported = {"pairwise", "matched_pairwise", "group_vs_group", "time_course"}
    if contrast.mode not in supported:
        raise NotImplementedError(
            f"Shared cluster tests are not implemented for contrast mode '{contrast.mode}'."
        )


def shared_cluster_tests(
    *,
    result: SharedClusterResult,
    contrast: ContrastSpec,
    test: str = "permutation",
    multiple_testing: str = "fdr_bh",
    n_permutations: int = 1000,
    random_state: int | None = 42,
    include_pairwise: bool = False,
) -> SharedClusterContrastResult:
    _require_supported_shared_cluster_mode(contrast)
    raise NotImplementedError("shared_cluster_tests is not implemented yet.")
