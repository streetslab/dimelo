# Shared Cluster Tests Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a new result-centric inference layer for global shared-cluster composition that consumes `SharedClusterResult`, reuses `ContrastSpec`, supports replicate-aware permutation testing plus optional pooled chi-squared / G-test screening, and returns contrast-style result tables plus plot-ready payloads.

**Architecture:** Keep `workflows.shared_cluster_distribution(...)` descriptive and add a separate `dimelo.shared_cluster_tests` module. The new module should consume existing shared-clustering tables, derive sample-level evidence, and return a `SharedClusterContrastResult` with `summary`, `details`, optional `pairwise`, `plot_data`, and `metadata`. Region-level occupancy inference remains in `dimelo.region_contrasts`, while a small workflow tweak propagates `SampleSpec.metadata` into shared-cluster assignments so `ContrastSpec.pairing_key` has a reliable metadata source.

**Tech Stack:** Python 3.11, pandas, numpy, scipy.stats, pytest, existing `dimelo.models`, `dimelo.workflows`, `dimelo.region_contrasts`, `dimelo.plotting`-style result conventions

---

## File Map

- Create: `dimelo/shared_cluster_tests.py`
  - public `shared_cluster_tests(...)`
  - contrast-mode validation
  - sample-level evidence builders
  - permutation, chi-squared, and G-test scoring helpers
  - result assembly and `plot_data`
- Modify: `dimelo/models.py`
  - add `SharedClusterContrastResult`
- Modify: `dimelo/__init__.py`
  - export `shared_cluster_tests`
- Modify: `dimelo/workflows.py`
  - propagate `SampleSpec.metadata` keys into shared-cluster assignment rows
- Create: `tests/test_shared_cluster_tests.py`
  - focused inference-layer coverage
- Modify: `tests/test_models.py`
  - result-model validation coverage
- Modify: `tests/test_workflows.py`
  - workflow metadata propagation coverage for paired designs
- Modify: `README.md`
  - mention `shared_cluster_tests(...)` and the boundary vs `region_contrasts`
- Modify: `docs/shared-clustering.md`
  - add user-facing examples and inference-boundary guidance

## Task 1: Add Failing Tests For The Result Model And Module Skeleton

**Files:**
- Create: `tests/test_shared_cluster_tests.py`
- Modify: `tests/test_models.py`
- Reference: `dimelo/shared_cluster_tests.py`
- Reference: `dimelo/models.py`

- [ ] **Step 1: Write the failing module-import and result-model tests**

Add this near the existing result-model coverage in `tests/test_models.py`:

```python
from dimelo.models import SharedClusterContrastResult


def test_shared_cluster_contrast_result_requires_summary_details_and_plot_data():
    with pytest.raises(ValueError, match="SharedClusterContrastResult requires non-None values"):
        SharedClusterContrastResult(
            summary=None,
            details=pd.DataFrame(),
            plot_data={},
        )
```

Create `tests/test_shared_cluster_tests.py` with:

```python
import pandas as pd


def test_shared_cluster_tests_module_exports_entry_point():
    from dimelo import shared_cluster_tests

    assert hasattr(shared_cluster_tests, "shared_cluster_tests")


def test_shared_cluster_tests_rejects_unsupported_contrast_mode():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = _make_shared_cluster_test_result()

    with pytest.raises(NotImplementedError, match="background_adjusted"):
        shared_cluster_tests.shared_cluster_tests(
            result=result,
            contrast=ContrastSpec(
                mode="background_adjusted",
                numerator=["treated"],
                denominator=["NS"],
                background=["bg"],
            ),
        )
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py tests/test_shared_cluster_tests.py -k "shared_cluster_contrast_result or shared_cluster_tests_module_exports_entry_point or shared_cluster_tests_rejects_unsupported_contrast_mode" -q
```

Expected: FAIL because `SharedClusterContrastResult` and `dimelo.shared_cluster_tests` do not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_models.py tests/test_shared_cluster_tests.py
git commit -m "test: add shared cluster tests skeleton coverage"
```

## Task 2: Add `SharedClusterContrastResult` And The Module Skeleton

**Files:**
- Create: `dimelo/shared_cluster_tests.py`
- Modify: `dimelo/models.py`
- Modify: `dimelo/__init__.py`
- Modify: `tests/test_models.py`
- Modify: `tests/test_shared_cluster_tests.py`

- [ ] **Step 1: Add the new result model**

Add this dataclass in `dimelo/models.py` near the other result models:

```python
@dataclass
class SharedClusterContrastResult:
    summary: pd.DataFrame
    details: pd.DataFrame
    pairwise: pd.DataFrame | None = None
    plot_data: dict[str, pd.DataFrame | dict[str, Any]] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "summary": self.summary,
            "details": self.details,
            "plot_data": self.plot_data,
            "metadata": self.metadata,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "SharedClusterContrastResult requires non-None values for: "
                f"{', '.join(missing)}"
            )
```

- [ ] **Step 2: Create the module skeleton and export it**

Create `dimelo/shared_cluster_tests.py` with:

```python
from __future__ import annotations

from dimelo.models import ContrastSpec, SharedClusterContrastResult


def _require_supported_shared_cluster_mode(contrast: ContrastSpec) -> None:
    supported = {"pairwise", "matched_pairwise", "group_vs_group", "time_course"}
    if contrast.mode not in supported:
        raise NotImplementedError(
            f"Shared cluster tests are not implemented for contrast mode '{contrast.mode}'."
        )


def shared_cluster_tests(
    *,
    result,
    contrast: ContrastSpec,
    test: str = "permutation",
    multiple_testing: str = "fdr_bh",
    n_permutations: int = 1000,
    random_state: int | None = 42,
    include_pairwise: bool = False,
) -> SharedClusterContrastResult:
    _require_supported_shared_cluster_mode(contrast)
    raise NotImplementedError("shared_cluster_tests is not implemented yet.")
```

Update `dimelo/__init__.py`:

```python
from . import (
    ...
    shared_cluster_tests,
    ...
)

__all__ = [
    ...
    "shared_cluster_tests",
    ...
]
```

- [ ] **Step 3: Run the focused tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py tests/test_shared_cluster_tests.py -k "shared_cluster_contrast_result or shared_cluster_tests_module_exports_entry_point or shared_cluster_tests_rejects_unsupported_contrast_mode" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the skeleton**

```bash
git add dimelo/models.py dimelo/shared_cluster_tests.py dimelo/__init__.py tests/test_models.py tests/test_shared_cluster_tests.py
git commit -m "feat: add shared cluster tests module skeleton"
```

## Task 3: Add Failing Tests For Pairing-Metadata Propagation From Workflows

**Files:**
- Modify: `tests/test_workflows.py`
- Modify: `dimelo/workflows.py`
- Reference: `dimelo/shared_cluster_tests.py`

- [ ] **Step 1: Write the failing workflow metadata propagation test**

Add this near the existing shared-clustering workflow tests in `tests/test_workflows.py`:

```python
def test_shared_cluster_distribution_propagates_sample_metadata_into_assignments(monkeypatch):
    monkeypatch.setattr(workflows.cluster, "extract_read_windows", _fake_extract_result)
    monkeypatch.setattr(workflows.cluster, "read_window_feature_matrix", _fake_feature_matrix)
    monkeypatch.setattr(workflows.cluster, "cluster_read_windows", _fake_cluster_fit_result)

    result = workflows.shared_cluster_distribution(
        samples=[
            SampleSpec(
                sample_id="s1",
                condition="NS",
                extract_h5="s1.h5",
                metadata={"subject_id": "mouse-1"},
            ),
            SampleSpec(
                sample_id="s2",
                condition="treated",
                extract_h5="s2.h5",
                metadata={"subject_id": "mouse-2"},
            ),
        ],
        mode="read_global",
        motifs=["A,0"],
        n_clusters=2,
    )

    assert "subject_id" in result.assignments.columns
    assert set(result.assignments["subject_id"]) == {"mouse-1", "mouse-2"}
```

- [ ] **Step 2: Run the focused test to verify it fails**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_workflows.py -k "propagates_sample_metadata_into_assignments" -q
```

Expected: FAIL because `SampleSpec.metadata` is not yet propagated into shared-cluster assignments.

- [ ] **Step 3: Implement the minimal workflow tweak**

Update the metadata-row construction in `dimelo/workflows.py`:

```python
        for metadata in extracted.metadata:
            row = {
                "sample_id": sample.sample_id,
                "condition": sample.condition,
                "replicate": sample.replicate,
            }
            if sample.metadata:
                row.update(sample.metadata)
            row.update(metadata)
            metadata_rows.append(row)
```

For the region-anchored path, add the same metadata merge before
`_build_shared_cluster_result(...)` by enriching `metadata_rows` by `sample_id`.

- [ ] **Step 4: Run the focused test to verify it passes**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_workflows.py -k "propagates_sample_metadata_into_assignments" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the workflow metadata support**

```bash
git add dimelo/workflows.py tests/test_workflows.py
git commit -m "feat: propagate shared cluster sample metadata"
```

## Task 4: Add Failing Tests For Pairwise And Group Permutation Inference

**Files:**
- Modify: `tests/test_shared_cluster_tests.py`
- Reference: `dimelo/shared_cluster_tests.py`

- [ ] **Step 1: Add a compact shared-cluster inference fixture**

Add this helper in `tests/test_shared_cluster_tests.py`:

```python
def _make_shared_cluster_test_result() -> SharedClusterResult:
    return SharedClusterResult(
        model=SharedClusterModel(
            mode="read_global",
            motifs=["A,0"],
            feature_names=["f0", "f1"],
            preprocessing={"signal_normalization": "none"},
            estimator=object(),
            cluster_labels=["C0", "C1"],
            fit_metadata={"clusterer": "minibatch_kmeans", "n_clusters": 2},
        ),
        assignments=pd.DataFrame(
            {
                "sample_id": ["s1", "s2", "s3", "s4"],
                "condition": ["NS", "NS", "treated", "treated"],
                "subject_id": ["p1", "p2", "p1", "p2"],
                "cluster": ["C0", "C0", "C1", "C1"],
            }
        ),
        cluster_distribution=pd.DataFrame(
            [
                {"sample_id": "s1", "condition": "NS", "cluster": "C0", "count": 80, "fraction": 0.80},
                {"sample_id": "s1", "condition": "NS", "cluster": "C1", "count": 20, "fraction": 0.20},
                {"sample_id": "s2", "condition": "NS", "cluster": "C0", "count": 75, "fraction": 0.75},
                {"sample_id": "s2", "condition": "NS", "cluster": "C1", "count": 25, "fraction": 0.25},
                {"sample_id": "s3", "condition": "treated", "cluster": "C0", "count": 30, "fraction": 0.30},
                {"sample_id": "s3", "condition": "treated", "cluster": "C1", "count": 70, "fraction": 0.70},
                {"sample_id": "s4", "condition": "treated", "cluster": "C0", "count": 25, "fraction": 0.25},
                {"sample_id": "s4", "condition": "treated", "cluster": "C1", "count": 75, "fraction": 0.75},
            ]
        ),
        condition_distribution=pd.DataFrame(
            [
                {"condition": "NS", "cluster": "C0", "count": 155, "fraction": 0.775, "replicate_n": 2},
                {"condition": "NS", "cluster": "C1", "count": 45, "fraction": 0.225, "replicate_n": 2},
                {"condition": "treated", "cluster": "C0", "count": 55, "fraction": 0.275, "replicate_n": 2},
                {"condition": "treated", "cluster": "C1", "count": 145, "fraction": 0.725, "replicate_n": 2},
            ]
        ),
        distribution_change=None,
        cluster_profiles=pd.DataFrame(columns=["cluster", "count", "f0", "f1"]),
        region_summaries=None,
        plot_data={},
        metadata={},
    )
```

- [ ] **Step 2: Write failing pairwise and matched-pairwise tests**

Add these tests:

```python
def test_shared_cluster_tests_pairwise_returns_summary_details_and_plot_data():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
        test="permutation",
        n_permutations=50,
        random_state=7,
    )

    assert set(result.summary.columns) >= {
        "contrast_id", "composition_effect_size", "omnibus_p_value", "top_cluster"
    }
    assert set(result.details.columns) >= {
        "cluster", "fraction", "reference_fraction", "delta_fraction", "p_value", "adjusted_p_value"
    }
    assert set(result.plot_data) >= {"summary_table", "cluster_effect_table"}


def test_shared_cluster_tests_matched_pairwise_uses_contrast_pairing_key():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["treated"],
            denominator=["NS"],
            pairing_key="subject_id",
        ),
        test="permutation",
        n_permutations=50,
        random_state=7,
    )

    assert result.metadata["paired"] is True
    assert result.metadata["pairing_key"] == "subject_id"
```

- [ ] **Step 3: Run the focused tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "pairwise_returns_summary_details_and_plot_data or matched_pairwise_uses_contrast_pairing_key" -q
```

Expected: FAIL because `shared_cluster_tests(...)` is still a stub.

- [ ] **Step 4: Implement pairwise/group permutation helpers**

In `dimelo/shared_cluster_tests.py`, add the following minimal helper set:

```python
def _sample_fraction_table(result) -> pd.DataFrame:
    required = {"sample_id", "condition", "cluster", "fraction"}
    missing = required - set(result.cluster_distribution.columns)
    if missing:
        raise ValueError(f"Shared cluster tests require cluster_distribution columns: {', '.join(sorted(missing))}.")
    return result.cluster_distribution.loc[:, ["sample_id", "condition", "cluster", "fraction"]].copy()


def _cluster_order(result) -> list[str]:
    return list(result.model.cluster_labels)


def _composition_effect_size(summary_table: pd.DataFrame) -> float:
    return float(summary_table["delta_fraction"].abs().sum() / 2.0)
```

Then implement the pairwise/group path to:

- build sample-level cluster fraction evidence
- compute condition means
- compute per-cluster `fraction`, `reference_fraction`, `delta_fraction`, `log2_fc`
- compute omnibus statistic from the L1 / total-variation shift
- compute permutation p-values
- add BH-adjusted per-cluster p-values
- assemble `summary`, `details`, `plot_data`, and `metadata`

- [ ] **Step 5: Run the focused tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "pairwise_returns_summary_details_and_plot_data or matched_pairwise_uses_contrast_pairing_key" -q
```

Expected: PASS.

- [ ] **Step 6: Commit the pairwise/group inference path**

```bash
git add dimelo/shared_cluster_tests.py tests/test_shared_cluster_tests.py
git commit -m "feat: add shared cluster pairwise inference"
```

## Task 5: Add Failing Tests For Pooled Screening Tests

**Files:**
- Modify: `tests/test_shared_cluster_tests.py`
- Reference: `dimelo/shared_cluster_tests.py`

- [ ] **Step 1: Write failing chi-squared and G-test coverage**

Add:

```python
def test_shared_cluster_tests_supports_chi_squared_screen():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
        test="chi_squared",
    )

    assert result.metadata["inference_level"] == "pooled_screen"
    assert result.summary.loc[0, "test"] == "chi_squared"


def test_shared_cluster_tests_supports_g_test_screen():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_test_result(),
        contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
        test="g_test",
    )

    assert result.summary.loc[0, "test"] == "g_test"
    assert 0.0 <= result.summary.loc[0, "omnibus_p_value"] <= 1.0
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "chi_squared_screen or g_test_screen" -q
```

Expected: FAIL because pooled screening tests are not yet implemented.

- [ ] **Step 3: Implement pooled count screening**

Add helpers in `dimelo/shared_cluster_tests.py`:

```python
def _pooled_count_table(result, contrast: ContrastSpec) -> pd.DataFrame:
    table = result.cluster_distribution.loc[:, ["sample_id", "condition", "cluster", "count"]].copy()
    return table
```

Then implement:

- pooled numerator and denominator cluster-count vectors
- `scipy.stats.chi2_contingency(...)` for `chi_squared`
- `scipy.stats.power_divergence(..., lambda_="log-likelihood")` or equivalent for `g_test`
- summary row with pooled-screen metadata
- per-cluster descriptive `details` rows without replicate-aware p-values

- [ ] **Step 4: Run the focused tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "chi_squared_screen or g_test_screen" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the pooled screening path**

```bash
git add dimelo/shared_cluster_tests.py tests/test_shared_cluster_tests.py
git commit -m "feat: add shared cluster pooled screening tests"
```

## Task 6: Add Failing Tests For Time-Course Inference

**Files:**
- Modify: `tests/test_shared_cluster_tests.py`
- Reference: `dimelo/shared_cluster_tests.py`

- [ ] **Step 1: Add a compact time-course fixture and failing tests**

Add:

```python
def _make_shared_cluster_time_course_result() -> SharedClusterResult:
    result = _make_shared_cluster_test_result()
    result.cluster_distribution = pd.DataFrame(
        [
            {"sample_id": "t0_a", "condition": "t0", "cluster": "C0", "count": 80, "fraction": 0.80},
            {"sample_id": "t0_a", "condition": "t0", "cluster": "C1", "count": 20, "fraction": 0.20},
            {"sample_id": "t1_a", "condition": "t1", "cluster": "C0", "count": 55, "fraction": 0.55},
            {"sample_id": "t1_a", "condition": "t1", "cluster": "C1", "count": 45, "fraction": 0.45},
            {"sample_id": "t2_a", "condition": "t2", "cluster": "C0", "count": 25, "fraction": 0.25},
            {"sample_id": "t2_a", "condition": "t2", "cluster": "C1", "count": 75, "fraction": 0.75},
        ]
    )
    return result


def test_shared_cluster_tests_time_course_returns_omnibus_and_trend_outputs():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_time_course_result(),
        contrast=ContrastSpec(mode="time_course", time_order=["t0", "t1", "t2"]),
        test="permutation",
        n_permutations=50,
        random_state=7,
    )

    assert {"omnibus_p_value", "trend_p_value", "composition_effect_size"} <= set(result.summary.columns)
    assert "time_course_table" in result.plot_data


def test_shared_cluster_tests_time_course_optional_pairwise_follow_up():
    from dimelo import shared_cluster_tests
    from dimelo.models import ContrastSpec

    result = shared_cluster_tests.shared_cluster_tests(
        result=_make_shared_cluster_time_course_result(),
        contrast=ContrastSpec(mode="time_course", time_order=["t0", "t1", "t2"]),
        test="permutation",
        include_pairwise=True,
        n_permutations=20,
        random_state=7,
    )

    assert result.pairwise is not None
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "time_course_returns_omnibus_and_trend_outputs or time_course_optional_pairwise_follow_up" -q
```

Expected: FAIL because time-course inference is not yet implemented.

- [ ] **Step 3: Implement time-course support**

In `dimelo/shared_cluster_tests.py`, add:

- a time-order validator
- an omnibus permutation statistic across timepoints
- a trend statistic using ordered time indices
- optional pairwise follow-up expansion when `include_pairwise=True`

Use helper signatures like:

```python
def _require_time_order_present(table: pd.DataFrame, time_order: list[str]) -> None:
    ...


def _score_time_course(
    *,
    sample_table: pd.DataFrame,
    contrast: ContrastSpec,
    n_permutations: int,
    random_state: int | None,
    include_pairwise: bool,
) -> SharedClusterContrastResult:
    ...
```

- [ ] **Step 4: Run the focused tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "time_course_returns_omnibus_and_trend_outputs or time_course_optional_pairwise_follow_up" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the time-course path**

```bash
git add dimelo/shared_cluster_tests.py tests/test_shared_cluster_tests.py
git commit -m "feat: add shared cluster time course inference"
```

## Task 7: Add User-Facing Docs

**Files:**
- Modify: `README.md`
- Modify: `docs/shared-clustering.md`

- [ ] **Step 1: Add a short README note**

Add near the existing shared-clustering analysis guide text:

```md
Global shared-cluster composition inference is available through
`dimelo.shared_cluster_tests.shared_cluster_tests(...)`, which consumes a
`SharedClusterResult` plus a `ContrastSpec`. Per-region occupancy inference
continues to live in `region_contrasts.score_regions(..., analysis_unit="cluster_occupancy")`.
```

- [ ] **Step 2: Add shared-clustering guide examples**

Add examples in `docs/shared-clustering.md`:

```python
from dimelo import shared_cluster_tests
from dimelo.models import ContrastSpec

contrast_result = shared_cluster_tests.shared_cluster_tests(
    result=result,
    contrast=ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["NS"],
    ),
)
```

And a paired example:

```python
contrast_result = shared_cluster_tests.shared_cluster_tests(
    result=result,
    contrast=ContrastSpec(
        mode="matched_pairwise",
        numerator=["treated"],
        denominator=["NS"],
        pairing_key="subject_id",
    ),
)
```

Also add boundary text that:

- `shared_cluster_tests(...)` is for global composition
- `region_contrasts.score_regions(...)` remains the route for region occupancy
- pooled chi-squared / G-test is screening-oriented

- [ ] **Step 3: Commit the docs**

```bash
git add README.md docs/shared-clustering.md
git commit -m "docs: add shared cluster tests guide"
```

## Task 8: Final Verification And Self-Review

**Files:**
- Reference: `dimelo/shared_cluster_tests.py`
- Reference: `dimelo/models.py`
- Reference: `dimelo/workflows.py`
- Reference: `tests/test_models.py`
- Reference: `tests/test_shared_cluster_tests.py`
- Reference: `tests/test_workflows.py`
- Reference: `README.md`
- Reference: `docs/shared-clustering.md`

- [ ] **Step 1: Run the focused shared-cluster inference suite**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -q
```

Expected: PASS.

- [ ] **Step 2: Run the model and workflow regression subset**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py tests/test_workflows.py -k "shared_cluster_contrast_result or shared_cluster_distribution" -q
```

Expected: PASS.

- [ ] **Step 3: Run a clean status check**

Run:

```bash
git status --short
```

Expected: clean working tree.

- [ ] **Step 4: Self-review against the spec**

Confirm explicitly that:

- `shared_cluster_distribution(...)` remains descriptive
- `shared_cluster_tests(...)` consumes `SharedClusterResult`
- `ContrastSpec` is reused rather than replaced
- replicate-aware permutation is the default path
- pooled chi-squared / G-test is explicit and labeled as screening
- pairwise, matched pairwise, group-vs-group, and time-course are covered
- `region_contrasts` remains the occupancy-inference route
- `plot_data` exists for both summary and per-cluster follow-up

- [ ] **Step 5: Commit a small follow-up only if verification finds a real issue**

```bash
git add dimelo/shared_cluster_tests.py dimelo/models.py dimelo/workflows.py tests/test_models.py tests/test_shared_cluster_tests.py tests/test_workflows.py README.md docs/shared-clustering.md
git commit -m "fix: tighten shared cluster tests support"
```

Only make this commit if final verification or self-review uncovers a concrete follow-up issue.
