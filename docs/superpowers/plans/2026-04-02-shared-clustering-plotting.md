# Shared Clustering Plotting Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add renderer-neutral plotting-prep helpers for `dimelo.shared_clustering` distribution/change views, cluster-profile views, and region-occupancy views.

**Architecture:** Extend the shared plotting module with three result-centric helpers that consume `SharedClusterResult` and emit stable payload dictionaries. Keep the implementation additive and data-first: preserve the existing lightweight `result.plot_data`, avoid clustering logic changes, and only aggregate from canonical `SharedClusterResult` tables already produced by `dimelo.workflows`.

**Tech Stack:** Python 3.11, pandas, pytest, existing `dimelo.plotting`, `dimelo.workflows`, `dimelo.distribution`, and `dimelo.models` modules

---

## File Map

- Modify: `dimelo/plotting.py`
  - add `SharedClusterResult` validation/helpers
  - add `prepare_shared_cluster_distribution_data(...)`
  - add `prepare_shared_cluster_profile_data(...)`
  - add `prepare_shared_cluster_region_data(...)`
- Modify: `tests/test_plotting.py`
  - add focused plotting helper coverage using synthetic `SharedClusterResult` fixtures
- Modify: `tests/test_workflows.py`
  - verify shared-clustering workflow outputs remain compatible with the new helper layer
- Modify: `docs/shared-clustering.md`
  - document the new plotting-prep helpers and how they relate to existing `result.plot_data`
- Modify: `README.md`
  - add one short note about shared-clustering plotting-prep helpers

## Task 1: Add Failing Tests For Distribution And Change Plot Prep

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/models.py`
- Reference: `dimelo/plotting.py`

- [ ] **Step 1: Write the failing distribution/change tests**

Add these tests near the other plotting helper tests in `tests/test_plotting.py`:

```python
from dimelo.models import SharedClusterModel, SharedClusterResult


def _make_shared_cluster_result() -> SharedClusterResult:
    model = SharedClusterModel(
        mode="region_anchored",
        motifs=["A,0"],
        feature_names=["f0", "f1"],
        preprocessing={"signal_normalization": "none"},
        estimator=object(),
        cluster_labels=["C0", "C1"],
        fit_metadata={"clusterer": "minibatch_kmeans", "n_clusters": 2},
    )
    assignments = pd.DataFrame(
        [
            {"sample_id": "s1", "condition": "NS", "cluster": "C0"},
            {"sample_id": "s1", "condition": "NS", "cluster": "C1"},
            {"sample_id": "s2", "condition": "treated", "cluster": "C1"},
        ]
    )
    cluster_distribution = pd.DataFrame(
        [
            {"sample_id": "s1", "condition": "NS", "cluster": "C0", "count": 2, "fraction": 2 / 3},
            {"sample_id": "s1", "condition": "NS", "cluster": "C1", "count": 1, "fraction": 1 / 3},
            {"sample_id": "s2", "condition": "treated", "cluster": "C0", "count": 1, "fraction": 1 / 4},
            {"sample_id": "s2", "condition": "treated", "cluster": "C1", "count": 3, "fraction": 3 / 4},
        ]
    )
    condition_distribution = pd.DataFrame(
        [
            {"condition": "NS", "cluster": "C0", "count": 2, "fraction": 2 / 3, "replicate_n": 1},
            {"condition": "NS", "cluster": "C1", "count": 1, "fraction": 1 / 3, "replicate_n": 1},
            {"condition": "treated", "cluster": "C0", "count": 1, "fraction": 1 / 4, "replicate_n": 1},
            {"condition": "treated", "cluster": "C1", "count": 3, "fraction": 3 / 4, "replicate_n": 1},
        ]
    )
    distribution_change = pd.DataFrame(
        [
            {
                "condition": "treated",
                "cluster": "C0",
                "count": 1,
                "fraction": 1 / 4,
                "replicate_n": 1,
                "reference_fraction": 2 / 3,
                "delta_fraction": -5 / 12,
                "log2_fc": -1.415037499278844,
            },
            {
                "condition": "treated",
                "cluster": "C1",
                "count": 3,
                "fraction": 3 / 4,
                "replicate_n": 1,
                "reference_fraction": 1 / 3,
                "delta_fraction": 5 / 12,
                "log2_fc": 1.1699250014423124,
            },
        ]
    )
    cluster_profiles = pd.DataFrame(
        [
            {"cluster": "C0", "count": 3, "f0": 0.1, "f1": 0.2},
            {"cluster": "C1", "count": 4, "f0": 0.8, "f1": 0.9},
        ]
    )
    region_summaries = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "cluster": "C0", "count": 2, "fraction": 2 / 3},
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "cluster": "C1", "count": 1, "fraction": 1 / 3},
            {"region_id": "reg1", "sample_id": "s2", "condition": "treated", "cluster": "C0", "count": 1, "fraction": 1 / 4},
            {"region_id": "reg1", "sample_id": "s2", "condition": "treated", "cluster": "C1", "count": 3, "fraction": 3 / 4},
        ]
    )
    return SharedClusterResult(
        model=model,
        assignments=assignments,
        cluster_distribution=cluster_distribution,
        condition_distribution=condition_distribution,
        distribution_change=distribution_change,
        cluster_profiles=cluster_profiles,
        region_summaries=region_summaries,
        plot_data={
            "cluster_distribution_bar": cluster_distribution.copy(),
            "cluster_distribution_heatmap": pd.DataFrame(
                [{"condition": "NS", "C0": 2 / 3, "C1": 1 / 3}, {"condition": "treated", "C0": 1 / 4, "C1": 3 / 4}]
            ),
        },
        metadata={"mode": "region_anchored"},
    )


def test_prepare_shared_cluster_distribution_data_returns_distribution_payload():
    result = _make_shared_cluster_result()

    payload = plotting.prepare_shared_cluster_distribution_data(result=result)

    assert set(payload) == {
        "sample_distribution",
        "condition_distribution",
        "distribution_change",
        "metadata",
    }
    assert payload["sample_distribution"]["sample_id"].tolist() == ["s1", "s1", "s2", "s2"]
    assert payload["condition_distribution"]["condition"].tolist() == ["NS", "NS", "treated", "treated"]
    assert payload["distribution_change"]["condition"].tolist() == ["treated", "treated"]
    assert payload["metadata"]["mode"] == "region_anchored"


def test_prepare_shared_cluster_distribution_data_handles_missing_change_table():
    result = _make_shared_cluster_result()
    result.distribution_change = None

    payload = plotting.prepare_shared_cluster_distribution_data(result=result)

    assert list(payload["distribution_change"].columns) == [
        "condition",
        "cluster",
        "count",
        "fraction",
        "replicate_n",
        "reference_fraction",
        "delta_fraction",
        "log2_fc",
    ]
    assert payload["distribution_change"].empty
```

- [ ] **Step 2: Run the new distribution/change tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "shared_cluster_distribution_data" -q
```

Expected: FAIL because `prepare_shared_cluster_distribution_data(...)` does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add shared clustering distribution plotting coverage"
```

## Task 2: Implement Distribution And Change Plot Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add the internal shared-clustering validation helpers**

Add these helpers in `dimelo/plotting.py` near the other result-validation utilities:

```python
def _validate_shared_cluster_result(result) -> None:
    if result is None:
        raise ValueError("plotting helpers require a SharedClusterResult.")
    required_attrs = (
        "cluster_distribution",
        "condition_distribution",
        "cluster_profiles",
        "plot_data",
        "metadata",
    )
    if not all(hasattr(result, attr) for attr in required_attrs):
        raise TypeError("plotting helpers require a SharedClusterResult-like object.")


def _empty_distribution_change_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "condition",
            "cluster",
            "count",
            "fraction",
            "replicate_n",
            "reference_fraction",
            "delta_fraction",
            "log2_fc",
        ]
    )
```

- [ ] **Step 2: Add `prepare_shared_cluster_distribution_data(...)`**

Add this public function in `dimelo/plotting.py` near the other public plotting helpers:

```python
def prepare_shared_cluster_distribution_data(
    *,
    result,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_shared_cluster_result(result)

    sample_distribution = prepare_cluster_distribution_bar_data(result.cluster_distribution)
    condition_distribution = (
        result.condition_distribution.loc[
            :, ["condition", "cluster", "count", "fraction", "replicate_n"]
        ]
        .sort_values(["condition", "cluster"], kind="stable")
        .reset_index(drop=True)
    )
    distribution_change = (
        _empty_distribution_change_table()
        if result.distribution_change is None
        else result.distribution_change.sort_values(["condition", "cluster"], kind="stable").reset_index(drop=True)
    )
    metadata = {
        "mode": result.model.mode,
        "cluster_labels": list(result.model.cluster_labels),
        "has_distribution_change": not distribution_change.empty,
    }
    return {
        "sample_distribution": sample_distribution,
        "condition_distribution": condition_distribution,
        "distribution_change": distribution_change,
        "metadata": metadata,
    }
```

- [ ] **Step 3: Run the distribution/change tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "shared_cluster_distribution_data" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the implementation**

```bash
git add dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add shared clustering distribution plotting prep"
```

## Task 3: Add Failing Tests For Cluster Profile Plot Prep

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting.py`

- [ ] **Step 1: Write the failing cluster-profile tests**

Add these tests in `tests/test_plotting.py` below the distribution helper coverage:

```python
def test_prepare_shared_cluster_profile_data_returns_long_form_profiles():
    result = _make_shared_cluster_result()

    payload = plotting.prepare_shared_cluster_profile_data(result=result)

    assert set(payload) == {"profile_table", "metadata"}
    assert set(payload["profile_table"].columns) == {"cluster", "feature", "value", "count"}
    assert payload["profile_table"]["cluster"].tolist() == ["C0", "C0", "C1", "C1"]
    assert payload["profile_table"]["feature"].tolist() == ["f0", "f1", "f0", "f1"]
    assert payload["metadata"]["feature_names"] == ["f0", "f1"]


def test_prepare_shared_cluster_profile_data_respects_feature_subset():
    result = _make_shared_cluster_result()

    payload = plotting.prepare_shared_cluster_profile_data(result=result, features=["f1"])

    assert payload["profile_table"]["feature"].unique().tolist() == ["f1"]
    assert payload["profile_table"]["cluster"].tolist() == ["C0", "C1"]
```

- [ ] **Step 2: Run the new cluster-profile tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "shared_cluster_profile_data" -q
```

Expected: FAIL because `prepare_shared_cluster_profile_data(...)` does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add shared clustering profile plotting coverage"
```

## Task 4: Implement Cluster Profile Plot Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add the cluster-profile helper**

Add this public function in `dimelo/plotting.py`:

```python
def prepare_shared_cluster_profile_data(
    *,
    result,
    features: list[str] | None = None,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_shared_cluster_result(result)
    cluster_profiles = result.cluster_profiles.copy()
    if cluster_profiles.empty:
        return {
            "profile_table": pd.DataFrame(columns=["cluster", "feature", "value", "count"]),
            "metadata": {"feature_names": [], "cluster_labels": list(result.model.cluster_labels)},
        }

    feature_names = [column for column in cluster_profiles.columns if column not in {"cluster", "count"}]
    if features is not None:
        missing = [feature for feature in features if feature not in feature_names]
        if missing:
            raise ValueError(f"Requested features are not present in cluster_profiles: {', '.join(missing)}")
        feature_names = [feature for feature in feature_names if feature in features]

    profile_table = (
        cluster_profiles.loc[:, ["cluster", "count", *feature_names]]
        .melt(id_vars=["cluster", "count"], var_name="feature", value_name="value")
        .sort_values(["cluster", "feature"], kind="stable")
        .reset_index(drop=True)
    )
    return {
        "profile_table": profile_table,
        "metadata": {
            "feature_names": feature_names,
            "cluster_labels": list(result.model.cluster_labels),
        },
    }
```

- [ ] **Step 2: Run the cluster-profile tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "shared_cluster_profile_data" -q
```

Expected: PASS.

- [ ] **Step 3: Commit the implementation**

```bash
git add dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add shared clustering profile plotting prep"
```

## Task 5: Add Failing Tests For Region Occupancy Plot Prep

**Files:**
- Modify: `tests/test_plotting.py`
- Modify: `tests/test_workflows.py`
- Reference: `dimelo/workflows.py`

- [ ] **Step 1: Write the failing region-occupancy tests**

Add these tests in `tests/test_plotting.py`:

```python
def test_prepare_shared_cluster_region_data_returns_sample_and_condition_tables():
    result = _make_shared_cluster_result()

    payload = plotting.prepare_shared_cluster_region_data(result=result)

    assert set(payload) == {"region_table", "condition_region_table", "metadata"}
    assert payload["region_table"]["sample_id"].tolist() == ["s1", "s1", "s2", "s2"]
    assert payload["condition_region_table"]["condition"].tolist() == ["NS", "NS", "treated", "treated"]
    assert payload["metadata"]["mode"] == "region_anchored"


def test_prepare_shared_cluster_region_data_can_disable_condition_aggregation():
    result = _make_shared_cluster_result()

    payload = plotting.prepare_shared_cluster_region_data(
        result=result,
        aggregate_conditions=False,
    )

    assert payload["condition_region_table"].empty


def test_prepare_shared_cluster_region_data_rejects_missing_region_summaries():
    result = _make_shared_cluster_result()
    result.region_summaries = None

    with pytest.raises(ValueError, match="region_summaries"):
        plotting.prepare_shared_cluster_region_data(result=result)
```

Add this workflow compatibility test in `tests/test_workflows.py` near the existing shared-clustering workflow tests:

```python
def test_shared_cluster_distribution_region_anchored_region_summaries_feed_region_plotting(monkeypatch):
    fake_samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            regions_bed="r1.bed",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="treated",
            extract_h5="s2.h5",
            regions_bed="r2.bed",
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]

    def fake_region_table(*args, **kwargs):
        return np.array([[0.1, 0.9], [0.9, 0.1]]), [
            {
                "region_id": "reg1",
                "sample_id": "s1",
                "condition": "NS",
                "replicate": None,
                "chromosome": "chr1",
                "start": 0,
                "end": 2,
                "strand": "+",
            },
            {
                "region_id": "reg1",
                "sample_id": "s2",
                "condition": "treated",
                "replicate": None,
                "chromosome": "chr1",
                "start": 0,
                "end": 2,
                "strand": "+",
            },
        ]

    monkeypatch.setattr(workflows.region_analysis, "build_region_feature_table", fake_region_table)

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="region_anchored",
        motifs=["A,0"],
        matched_regions="matched.bed",
        n_clusters=2,
        make_plots=False,
    )

    payload = plotting.prepare_shared_cluster_region_data(result=result)

    assert not payload["region_table"].empty
    assert list(payload["condition_region_table"]["condition"].unique()) == ["NS", "treated"]
```

- [ ] **Step 2: Run the new region-occupancy tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py tests/test_workflows.py -k "shared_cluster_region_data" -q
```

Expected: FAIL because `prepare_shared_cluster_region_data(...)` does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py tests/test_workflows.py
git commit -m "test: add shared clustering region plotting coverage"
```

## Task 6: Implement Region Occupancy Plot Prep And Workflow Alignment

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`
- Modify: `tests/test_workflows.py`

- [ ] **Step 1: Add the region-occupancy helper**

Add this public function in `dimelo/plotting.py`:

```python
def prepare_shared_cluster_region_data(
    *,
    result,
    aggregate_conditions: bool = True,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_shared_cluster_result(result)
    if result.region_summaries is None:
        raise ValueError("SharedClusterResult.region_summaries is required for region plotting.")

    region_table = result.region_summaries.sort_values(
        ["region_id", "sample_id", "cluster"],
        kind="stable",
    ).reset_index(drop=True)

    if aggregate_conditions and not region_table.empty:
        condition_region_table = (
            region_table.groupby(["region_id", "condition", "cluster"], as_index=False, sort=False)
            .agg(
                count=("count", "sum"),
                fraction_mean=("fraction", "mean"),
                fraction_median=("fraction", "median"),
                sample_n=("sample_id", "nunique"),
            )
            .sort_values(["region_id", "condition", "cluster"], kind="stable")
            .reset_index(drop=True)
        )
    else:
        condition_region_table = pd.DataFrame(
            columns=[
                "region_id",
                "condition",
                "cluster",
                "count",
                "fraction_mean",
                "fraction_median",
                "sample_n",
            ]
        )

    return {
        "region_table": region_table,
        "condition_region_table": condition_region_table,
        "metadata": {
            "mode": result.model.mode,
            "cluster_labels": list(result.model.cluster_labels),
            "has_condition_aggregation": aggregate_conditions,
        },
    }
```

- [ ] **Step 2: Keep workflow compatibility explicit in tests**

Do not change `dimelo/workflows.py` unless a test reveals a real mismatch. If the new helper works directly on existing `SharedClusterResult.region_summaries`, keep the workflow untouched and preserve the additive boundary.

- [ ] **Step 3: Run the region-occupancy tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py tests/test_workflows.py -k "shared_cluster_region_data or region_summaries_feed_region_plotting" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the implementation**

```bash
git add dimelo/plotting.py tests/test_plotting.py tests/test_workflows.py
git commit -m "feat: add shared clustering region plotting prep"
```

## Task 7: Document The New Shared Clustering Plotting Helpers

**Files:**
- Modify: `docs/shared-clustering.md`
- Modify: `README.md`

- [ ] **Step 1: Update the shared clustering guide**

Add a short “Plotting Prep Helpers” section to `docs/shared-clustering.md` that documents the three new helpers and keeps the continuity note about existing lightweight `result.plot_data`:

```md
## Plotting Prep Helpers

`SharedClusterResult` now has three renderer-neutral plotting helpers in `dimelo.plotting`:

- `prepare_shared_cluster_distribution_data(result=...)`
- `prepare_shared_cluster_profile_data(result=...)`
- `prepare_shared_cluster_region_data(result=...)`

These helpers sit on top of the canonical result tables:

- `result.cluster_distribution`
- `result.condition_distribution`
- `result.distribution_change`
- `result.cluster_profiles`
- `result.region_summaries`

The older lightweight payloads in `result.plot_data["cluster_distribution_bar"]` and
`result.plot_data["cluster_distribution_heatmap"]` remain supported for backward familiarity.
```

- [ ] **Step 2: Add a short README note**

Add one short paragraph in `README.md` near the other workflow/plotting summaries:

```md
Shared clustering also has renderer-neutral plotting prep in `dimelo.plotting`. Use
`prepare_shared_cluster_distribution_data(...)` for sample/condition cluster fractions,
`prepare_shared_cluster_profile_data(...)` for cluster feature summaries, and
`prepare_shared_cluster_region_data(...)` for region-level occupancy tables.
```

- [ ] **Step 3: Commit the docs**

```bash
git add docs/shared-clustering.md README.md
git commit -m "docs: add shared clustering plotting guide"
```

## Task 8: Run Final Verification

**Files:**
- Reference: `dimelo/plotting.py`
- Reference: `tests/test_plotting.py`
- Reference: `tests/test_workflows.py`
- Reference: `docs/shared-clustering.md`
- Reference: `README.md`

- [ ] **Step 1: Run the focused plotting and workflow tests**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/xdg-cache PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py tests/test_workflows.py -q
```

Expected: PASS with no new failures. Existing non-blocking warnings are acceptable if they are unchanged.

- [ ] **Step 2: Run a quick git status check**

Run:

```bash
git status --short
```

Expected: clean working tree.

- [ ] **Step 3: Self-review against the spec**

Check the completed work against `docs/superpowers/specs/2026-04-02-shared-clustering-plotting-design.md`:

- distribution/change helper present and additive
- cluster-profile helper present and data-first
- region-occupancy helper present with optional condition aggregation
- existing `result.plot_data` preserved
- docs updated without claiming renderer support

- [ ] **Step 4: Commit any final doc/test touch-ups if needed**

```bash
git add dimelo/plotting.py tests/test_plotting.py tests/test_workflows.py docs/shared-clustering.md README.md
git commit -m "fix: tighten shared clustering plotting prep"
```

Only make this commit if the self-review or verification uncovered a small follow-up patch.
