# Cluster Occupancy Region Contrasts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `dimelo.region_contrasts` to support region-level cluster occupancy contrasts over sample-level occupancy summaries, including per-cluster fractions plus descriptive dominant-cluster and entropy outputs, with inferential support for per-cluster fractions in `pairwise` and `group_vs_group` modes.

**Architecture:** Keep the current `score_regions(...)` entry point and extend it with a second validated path for `analysis_unit="cluster_occupancy"`. Build one occupancy evidence table from clustering-derived region summaries, then score per-cluster fraction contrasts across `region x sample` observations while also emitting descriptive dominant-cluster and entropy summaries in the result tables.

**Tech Stack:** Python, pandas, existing dimelo contrast models, pytest

---

### Task 1: Extend Validation And Occupancy Evidence Builders

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/region_contrasts.py`
- Test: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_contrasts.py`

- [ ] **Step 1: Write the failing validation and evidence-table tests**

```python
def test_validate_region_contrast_request_accepts_cluster_occupancy_fraction_mode():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="effect_size_only",
    )


def test_validate_region_contrast_request_rejects_beta_binomial_for_cluster_occupancy():
    with pytest.raises(ValueError, match="cluster_occupancy"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="cluster_occupancy",
            representation="cluster_fraction",
            signal_source="cluster_occupancy",
            test="beta_binomial",
        )


def test_build_cluster_occupancy_evidence_table_summarizes_region_sample_clusters():
    evidence = region_contrasts.build_cluster_occupancy_evidence_table(
        region_summaries=_mock_region_summaries(),
    )
    assert {"region_id", "sample_id", "condition", "cluster", "fraction"} <= set(evidence.columns)
    assert "dominant_cluster" in evidence.columns
    assert "cluster_entropy" in evidence.columns
```

- [ ] **Step 2: Run the targeted tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -q`
Expected: FAIL with missing builder / validation support

- [ ] **Step 3: Extend request validation**

```python
def validate_region_contrast_request(...):
    if analysis_unit == "ensemble_region":
        ...
        return
    if analysis_unit == "cluster_occupancy":
        if signal_source != "cluster_occupancy":
            raise ValueError(...)
        if representation not in {"cluster_fraction", "dominant_cluster", "cluster_entropy"}:
            raise ValueError(...)
        if test not in {"effect_size_only", "fraction_test"}:
            raise ValueError(...)
        return
    raise ValueError(...)
```

- [ ] **Step 4: Add the occupancy evidence builder**

```python
def build_cluster_occupancy_evidence_table(*, region_summaries: pd.DataFrame) -> pd.DataFrame:
    required = {"region_id", "sample_id", "condition", "cluster", "count", "fraction"}
    missing = required - set(region_summaries.columns)
    if missing:
        raise ValueError(...)

    evidence = region_summaries.copy()
    totals = evidence.groupby(["region_id", "sample_id", "condition"])["count"].transform("sum")
    evidence["fraction"] = evidence["count"].div(totals.where(totals != 0), fill_value=0).fillna(0.0)

    dominant = (
        evidence.sort_values(["region_id", "sample_id", "condition", "fraction", "cluster"],
                             ascending=[True, True, True, False, True])
        .drop_duplicates(["region_id", "sample_id", "condition"])
        .loc[:, ["region_id", "sample_id", "condition", "cluster"]]
        .rename(columns={"cluster": "dominant_cluster"})
    )
    entropy = (
        evidence.groupby(["region_id", "sample_id", "condition"], as_index=False)
        .agg(cluster_entropy=("fraction", lambda values: float(-(values * np.log2(values.where(values > 0, 1))).sum())))
    )
    return evidence.merge(dominant, ...).merge(entropy, ...)
```

- [ ] **Step 5: Run the targeted tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add cluster occupancy contrast evidence builders"
```

### Task 2: Add Sample-Level Occupancy Scoring

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/region_contrasts.py`
- Test: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_contrasts.py`

- [ ] **Step 1: Write the failing scoring tests**

```python
def test_score_regions_cluster_fraction_effect_size_only_ranks_largest_fraction_shift_first():
    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        motifs=[],
        contrast=ContrastSpec(mode="pairwise", numerator=["15min"], denominator=["NS"]),
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="effect_size_only",
        occupancy_table=_mock_cluster_occupancy_evidence(),
    )
    assert result.regions.iloc[0]["cluster"] == "C1"
    assert "delta_fraction" in result.regions.columns


def test_score_regions_cluster_fraction_fraction_test_adds_p_values():
    result = region_contrasts.score_regions(
        ...,
        analysis_unit="cluster_occupancy",
        representation="cluster_fraction",
        signal_source="cluster_occupancy",
        test="fraction_test",
        occupancy_table=_mock_cluster_occupancy_evidence(),
    )
    assert {"p_value", "adjusted_p_value"} <= set(result.regions.columns)


def test_score_regions_cluster_entropy_returns_descriptive_summary_only():
    result = region_contrasts.score_regions(
        ...,
        analysis_unit="cluster_occupancy",
        representation="cluster_entropy",
        signal_source="cluster_occupancy",
        test="effect_size_only",
        occupancy_table=_mock_cluster_occupancy_evidence(),
    )
    assert "cluster_entropy" in result.regions.columns
    assert "p_value" not in result.regions.columns
```

- [ ] **Step 2: Run the targeted tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -q`
Expected: FAIL with unsupported cluster occupancy path

- [ ] **Step 3: Implement occupancy scoring helpers**

```python
def _pool_cluster_occupancy_groups(...):
    ...


def _add_fraction_test_scores(regions_table: pd.DataFrame, *, multiple_testing: str) -> pd.DataFrame:
    ...


def _score_cluster_occupancy(...):
    if representation == "cluster_fraction":
        ...
    elif representation == "dominant_cluster":
        ...
    else:
        ...
```

Rules for v1:
- inferential support only for `representation="cluster_fraction"`
- only `pairwise` and `group_vs_group`
- sample-level region occupancy is the unit; do not pool reads as pseudo-replicates
- `dominant_cluster` and `cluster_entropy` remain descriptive

- [ ] **Step 4: Extend score_regions(...) with the new path**

```python
def score_regions(..., occupancy_table: pd.DataFrame | None = None, ...):
    validate_region_contrast_request(...)
    if analysis_unit == "cluster_occupancy":
        evidence = occupancy_table if occupancy_table is not None else build_cluster_occupancy_evidence_table(...)
        ...
        return RegionContrastResult(...)
    ...
```

- [ ] **Step 5: Add rejection-path tests**

```python
def test_score_regions_cluster_occupancy_rejects_matched_pairwise():
    ...


def test_score_regions_cluster_occupancy_rejects_beta_binomial():
    ...


def test_score_regions_cluster_occupancy_rejects_missing_occupancy_columns():
    ...
```

- [ ] **Step 6: Run the targeted tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -q`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add cluster occupancy region scoring"
```

### Task 3: Document Occupancy Contrasts And Verify The Slice

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-contrasts.md`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/shared-clustering.md`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md`
- Test: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_contrasts.py`

- [ ] **Step 1: Add one occupancy-driven example**

```python
result = region_contrasts.score_regions(
    samples=samples,
    regions=selected_regions,
    motifs=["A,0"],
    contrast=ContrastSpec(
        mode="pairwise",
        numerator=["15min"],
        denominator=["NS"],
    ),
    analysis_unit="cluster_occupancy",
    representation="cluster_fraction",
    signal_source="cluster_occupancy",
    test="fraction_test",
    occupancy_table=occupancy_table,
)
```

- [ ] **Step 2: Make the v1 occupancy contract explicit**

```markdown
- sample-level `region x sample` occupancy summaries are the unit of comparison
- `cluster_fraction` supports `effect_size_only` and `fraction_test`
- `dominant_cluster` and `cluster_entropy` are descriptive-only in v1
- supported contrast modes are `pairwise` and `group_vs_group`
```

- [ ] **Step 3: Add/adjust one metadata test if needed**

```python
def test_score_regions_cluster_occupancy_metadata_marks_signal_source():
    ...
    assert result.metadata["analysis_unit"] == "cluster_occupancy"
    assert result.metadata["signal_source"] == "cluster_occupancy"
```

- [ ] **Step 4: Run the verification slice**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py tests/test_workflows.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add docs/region-contrasts.md docs/shared-clustering.md README.md tests/test_region_contrasts.py
git commit -m "docs: add cluster occupancy contrast guide"
```

---

## Self-Review

- Spec coverage: this plan adds the occupancy evidence builder, scoring path, inference/descriptive split, and docs for the new cluster-occupancy contrasts.
- Placeholder scan: every task names exact files, tests, and commands; there are no TODO placeholders.
- Type consistency: the same API terms are used throughout: `analysis_unit="cluster_occupancy"`, `representation in {"cluster_fraction", "dominant_cluster", "cluster_entropy"}`, `signal_source="cluster_occupancy"`, and `test in {"effect_size_only", "fraction_test"}`.
