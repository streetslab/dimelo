# Single-Read Region Contrasts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `dimelo.region_contrasts` with a first `analysis_unit="single_read"` family covering `read_mod_fraction` and `read_window_features` for defined-region contrasts.

**Architecture:** Keep all new behavior inside `region_contrasts.score_regions(...)` by adding explicit request validation, read-level evidence builders, and sample-aware scoring paths. Build `read_mod_fraction` first, then add `read_window_features` using the current built-in feature set by default and a validated user-supplied `feature_table` override, while keeping the inferential boundary at the sample level.

**Tech Stack:** Python 3.11, pandas, scipy, pytest, existing `dimelo.region_contrasts`, `dimelo.region_analysis`, `dimelo.cluster`, `dimelo.models`, and docs modules

---

## File Map

- Modify: `dimelo/region_contrasts.py`
  - extend request validation for `analysis_unit="single_read"`
  - add read-level evidence builders
  - add sample-aware scoring paths for `read_mod_fraction`
  - add feature-table validation and scoring for `read_window_features`
- Modify: `tests/test_region_contrasts.py`
  - add validation coverage
  - add synthetic evidence/scoring coverage for both single-read representations
  - add matched-pairwise coverage
- Modify: `docs/region-contrasts.md`
  - document `single_read` semantics, supported modes, and v1 constraints
- Modify: `README.md`
  - add one short note describing the new `single_read` contrast family

## Task 1: Add Failing Validation Tests For `analysis_unit=\"single_read\"`

**Files:**
- Modify: `tests/test_region_contrasts.py`
- Reference: `dimelo/region_contrasts.py`

- [ ] **Step 1: Write the failing validation tests**

Add these tests near the existing validation coverage in `tests/test_region_contrasts.py`:

```python
def test_validate_region_contrast_request_accepts_single_read_mod_fraction():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="single_read",
        representation="read_mod_fraction",
        signal_source="extract_reads",
        test="sample_distribution_shift",
    )


def test_validate_region_contrast_request_accepts_single_read_window_features():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="feature_summary_shift",
    )


def test_validate_region_contrast_request_rejects_single_read_wrong_signal_source():
    with pytest.raises(ValueError, match="extract_reads"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="pileup_counts",
            test="sample_distribution_shift",
        )


def test_validate_region_contrast_request_rejects_single_read_unknown_representation():
    with pytest.raises(ValueError, match="read_mod_fraction"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_shape",
            signal_source="extract_reads",
            test="sample_distribution_shift",
        )


def test_validate_region_contrast_request_rejects_single_read_unknown_test():
    with pytest.raises(ValueError, match="sample_distribution_shift"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="extract_reads",
            test="beta_binomial",
        )
```

- [ ] **Step 2: Run the validation tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -k "single_read and validate_region_contrast_request" -q
```

Expected: FAIL because `single_read` is not yet accepted by `validate_region_contrast_request(...)`.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_region_contrasts.py
git commit -m "test: add single-read contrast validation coverage"
```

## Task 2: Implement `single_read` Request Validation

**Files:**
- Modify: `dimelo/region_contrasts.py`
- Modify: `tests/test_region_contrasts.py`

- [ ] **Step 1: Extend `validate_region_contrast_request(...)`**

Add a `single_read` branch in `dimelo/region_contrasts.py`:

```python
    if analysis_unit == "single_read":
        if representation == "read_mod_fraction":
            if signal_source != "extract_reads":
                raise ValueError(
                    "Current single_read read_mod_fraction support requires "
                    "signal_source='extract_reads'."
                )
            if test not in {"effect_size_only", "sample_distribution_shift"}:
                raise ValueError(
                    "Current single_read read_mod_fraction scoring support requires "
                    "test in {'effect_size_only', 'sample_distribution_shift'}."
                )
            return

        if representation == "read_window_features":
            if signal_source != "extract_features":
                raise ValueError(
                    "Current single_read read_window_features support requires "
                    "signal_source='extract_features'."
                )
            if test not in {"effect_size_only", "feature_summary_shift"}:
                raise ValueError(
                    "Current single_read read_window_features support requires "
                    "test in {'effect_size_only', 'feature_summary_shift'}."
                )
            return

        raise ValueError(
            "Current single_read support requires representation='read_mod_fraction' "
            "or 'read_window_features'."
        )
```

Update the final error message so unsupported units mention `single_read`:

```python
    if analysis_unit not in {"ensemble_region", "cluster_occupancy", "single_read"}:
        raise ValueError(
            "V1 region_contrasts inference requires analysis_unit='ensemble_region', "
            "'cluster_occupancy', or 'single_read'."
        )
```

- [ ] **Step 2: Run the validation tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -k "single_read and validate_region_contrast_request" -q
```

Expected: PASS.

- [ ] **Step 3: Commit the validation implementation**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add single-read contrast request validation"
```

## Task 3: Add Failing Tests For `read_mod_fraction` Evidence Building

**Files:**
- Modify: `tests/test_region_contrasts.py`
- Reference: `dimelo/region_contrasts.py`

- [ ] **Step 1: Write the failing `read_mod_fraction` evidence tests**

Add a small synthetic helper and tests in `tests/test_region_contrasts.py`:

```python
def _single_read_samples():
    return [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            replicate=1,
            metadata={"extract_path": "s1.tsv"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="treated",
            extract_h5="s2.h5",
            replicate=1,
            metadata={"extract_path": "s2.tsv"},
        ),
    ]


def test_build_single_read_mod_fraction_evidence_table():
    samples = _single_read_samples()
    extract_table = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r1", "modified_count": 2, "valid_count": 4},
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r2", "modified_count": 1, "valid_count": 4},
            {"region_id": "reg1", "sample_id": "s2", "condition": "treated", "read_id": "r3", "modified_count": 4, "valid_count": 4},
        ]
    )

    result = region_contrasts.build_single_read_mod_fraction_evidence_table(
        extract_table=extract_table
    )

    assert result.to_dict("records") == [
        {
            "region_id": "reg1",
            "sample_id": "s1",
            "condition": "NS",
            "read_id": "r1",
            "modified_count": 2,
            "valid_count": 4,
            "read_mod_fraction": pytest.approx(0.5),
        },
        {
            "region_id": "reg1",
            "sample_id": "s1",
            "condition": "NS",
            "read_id": "r2",
            "modified_count": 1,
            "valid_count": 4,
            "read_mod_fraction": pytest.approx(0.25),
        },
        {
            "region_id": "reg1",
            "sample_id": "s2",
            "condition": "treated",
            "read_id": "r3",
            "modified_count": 4,
            "valid_count": 4,
            "read_mod_fraction": pytest.approx(1.0),
        },
    ]


def test_build_single_read_mod_fraction_evidence_table_rejects_missing_columns():
    with pytest.raises(ValueError, match="modified_count"):
        region_contrasts.build_single_read_mod_fraction_evidence_table(
            extract_table=pd.DataFrame([{"region_id": "reg1", "sample_id": "s1"}])
        )
```

- [ ] **Step 2: Run the evidence tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -k "single_read_mod_fraction_evidence" -q
```

Expected: FAIL because the evidence builder does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_region_contrasts.py
git commit -m "test: add single-read mod-fraction evidence coverage"
```

## Task 4: Implement `read_mod_fraction` Evidence And Sample-Aware Scoring

**Files:**
- Modify: `dimelo/region_contrasts.py`
- Modify: `tests/test_region_contrasts.py`

- [ ] **Step 1: Add the evidence builder**

Add this helper in `dimelo/region_contrasts.py` near the other evidence builders:

```python
def build_single_read_mod_fraction_evidence_table(
    *,
    extract_table: pd.DataFrame,
) -> pd.DataFrame:
    required_columns = {
        "region_id",
        "sample_id",
        "condition",
        "read_id",
        "modified_count",
        "valid_count",
    }
    missing_columns = required_columns - set(extract_table.columns)
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "build_single_read_mod_fraction_evidence_table requires columns: "
            f"{missing_display}."
        )

    evidence = extract_table.loc[
        :,
        [
            "region_id",
            "sample_id",
            "condition",
            "read_id",
            "modified_count",
            "valid_count",
        ],
    ].copy()
    evidence["read_mod_fraction"] = evidence["modified_count"].div(
        evidence["valid_count"].where(evidence["valid_count"] != 0),
        fill_value=0,
    ).fillna(0.0)
    return evidence
```

- [ ] **Step 2: Add failing scoring tests for `read_mod_fraction`**

Add these tests in `tests/test_region_contrasts.py`:

```python
def test_score_regions_single_read_mod_fraction_effect_size_only():
    evidence = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r1", "modified_count": 1, "valid_count": 4, "read_mod_fraction": 0.25},
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r2", "modified_count": 2, "valid_count": 4, "read_mod_fraction": 0.50},
            {"region_id": "reg1", "sample_id": "s2", "condition": "treated", "read_id": "r3", "modified_count": 4, "valid_count": 4, "read_mod_fraction": 1.00},
            {"region_id": "reg1", "sample_id": "s3", "condition": "treated", "read_id": "r4", "modified_count": 3, "valid_count": 4, "read_mod_fraction": 0.75},
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        contrast=ContrastSpec(mode="group_vs_group", numerator=["treated"], denominator=["NS"]),
        analysis_unit="single_read",
        representation="read_mod_fraction",
        signal_source="extract_reads",
        test="effect_size_only",
        read_table=evidence,
    )

    row = result.summary.iloc[0]
    assert row["region_id"] == "reg1"
    assert row["sample_summary_numerator_mean"] == pytest.approx(0.875)
    assert row["sample_summary_denominator_mean"] == pytest.approx(0.375)
    assert row["delta_summary_mean"] == pytest.approx(0.5)
```

- [ ] **Step 3: Implement the sample-aware scoring path**

Add internal helpers in `dimelo/region_contrasts.py`:

```python
def _summarize_single_read_mod_fraction_by_sample(evidence: pd.DataFrame) -> pd.DataFrame:
    return (
        evidence.groupby(["region_id", "sample_id", "condition"], as_index=False, sort=False)
        .agg(
            read_n=("read_id", "nunique"),
            sample_summary_mean=("read_mod_fraction", "mean"),
            sample_summary_median=("read_mod_fraction", "median"),
            sample_summary_var=("read_mod_fraction", "var"),
        )
        .fillna({"sample_summary_var": 0.0})
    )


def _score_single_read_mod_fraction(
    *,
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
    test: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_summary = _summarize_single_read_mod_fraction_by_sample(evidence)
    pooled = _pool_region_groups(
        sample_summary.rename(columns={"sample_summary_mean": "mod_fraction"})[
            ["region_id", "sample_id", "condition", "mod_fraction"]
        ],
        contrast=contrast,
    )
    pooled = pooled.rename(
        columns={
            "fraction": "sample_summary_mean",
            "reference_fraction": "sample_summary_reference_mean",
            "delta_fraction": "delta_summary_mean",
        }
    )
    summary = pooled.sort_values("region_id", kind="stable").reset_index(drop=True)
    return evidence, summary
```

Wire it into `score_regions(...)` behind:

```python
if analysis_unit == "single_read" and representation == "read_mod_fraction":
    if read_table is None:
        raise ValueError("single_read read_mod_fraction scoring currently requires read_table.")
    evidence = build_single_read_mod_fraction_evidence_table(extract_table=read_table)
    regions_table, summary = _score_single_read_mod_fraction(
        evidence=evidence,
        contrast=contrast,
        test=test,
    )
```

Keep `test="sample_distribution_shift"` mapped to the same first sample-aware summary path in v1.

- [ ] **Step 4: Run the `read_mod_fraction` tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -k "single_read_mod_fraction" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the implementation**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add single-read mod-fraction contrasts"
```

## Task 5: Add Failing Tests For `read_window_features`

**Files:**
- Modify: `tests/test_region_contrasts.py`
- Reference: `dimelo/region_contrasts.py`

- [ ] **Step 1: Write the failing feature-evidence and scoring tests**

Add these tests:

```python
def test_build_single_read_feature_evidence_table_accepts_user_features():
    feature_table = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r1", "f0": 0.1, "f1": 0.2},
            {"region_id": "reg1", "sample_id": "s2", "condition": "treated", "read_id": "r2", "f0": 0.8, "f1": 0.9},
        ]
    )

    result = region_contrasts.build_single_read_feature_evidence_table(feature_table=feature_table)

    assert result.to_dict("records") == feature_table.to_dict("records")


def test_build_single_read_feature_evidence_table_rejects_missing_read_id():
    with pytest.raises(ValueError, match="read_id"):
        region_contrasts.build_single_read_feature_evidence_table(
            feature_table=pd.DataFrame([{"region_id": "reg1", "sample_id": "s1", "condition": "NS", "f0": 0.1}])
        )


def test_score_regions_single_read_window_features_effect_size_only():
    feature_table = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r1", "f0": 0.1, "f1": 0.2},
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r2", "f0": 0.2, "f1": 0.3},
            {"region_id": "reg1", "sample_id": "s2", "condition": "treated", "read_id": "r3", "f0": 0.8, "f1": 0.9},
            {"region_id": "reg1", "sample_id": "s3", "condition": "treated", "read_id": "r4", "f0": 0.7, "f1": 0.8},
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        contrast=ContrastSpec(mode="group_vs_group", numerator=["treated"], denominator=["NS"]),
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="effect_size_only",
        feature_table=feature_table,
    )

    row = result.summary.iloc[0]
    assert row["region_id"] == "reg1"
    assert row["f0_delta_mean"] == pytest.approx(0.6)
    assert row["f1_delta_mean"] == pytest.approx(0.6)
```

- [ ] **Step 2: Run the feature tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -k "single_read_feature or single_read_window_features" -q
```

Expected: FAIL because the feature evidence builder and scoring path do not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_region_contrasts.py
git commit -m "test: add single-read feature contrast coverage"
```

## Task 6: Implement `read_window_features` With Built-In Default And User Override

**Files:**
- Modify: `dimelo/region_contrasts.py`
- Modify: `tests/test_region_contrasts.py`

- [ ] **Step 1: Add the user-feature evidence builder**

Add this helper:

```python
def build_single_read_feature_evidence_table(
    *,
    feature_table: pd.DataFrame,
) -> pd.DataFrame:
    required_columns = {"region_id", "sample_id", "condition", "read_id"}
    missing_columns = required_columns - set(feature_table.columns)
    if missing_columns:
        missing_display = ", ".join(sorted(missing_columns))
        raise ValueError(
            "build_single_read_feature_evidence_table requires columns: "
            f"{missing_display}."
        )

    feature_columns = [
        column for column in feature_table.columns if column not in required_columns
    ]
    if not feature_columns:
        raise ValueError(
            "build_single_read_feature_evidence_table requires at least one feature column."
        )
    return feature_table.loc[:, ["region_id", "sample_id", "condition", "read_id", *feature_columns]].copy()
```

- [ ] **Step 2: Add a built-in feature-source shim**

Add a narrow internal function that is easy to replace later with the real built-in feature extraction path:

```python
def _load_builtin_single_read_feature_table(*, samples, regions, motifs):
    raise NotImplementedError(
        "Built-in single_read feature loading is not implemented yet."
    )
```

Then wire `score_regions(...)` so:

- if `representation="read_window_features"` and `feature_table is not None`, use the user table
- if `feature_table is None`, call `_load_builtin_single_read_feature_table(...)`

This keeps the built-in default contract explicit in the code and easy to patch in tests.

- [ ] **Step 3: Implement feature-summary scoring**

Add helpers:

```python
def _summarize_single_read_features_by_sample(evidence: pd.DataFrame) -> pd.DataFrame:
    group_columns = ["region_id", "sample_id", "condition"]
    feature_columns = [
        column for column in evidence.columns if column not in {"region_id", "sample_id", "condition", "read_id"}
    ]
    return (
        evidence.groupby(group_columns, as_index=False, sort=False)[feature_columns]
        .mean()
    )


def _score_single_read_features(
    *,
    evidence: pd.DataFrame,
    contrast: ContrastSpec,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_summary = _summarize_single_read_features_by_sample(evidence)
    feature_columns = [
        column for column in sample_summary.columns if column not in {"region_id", "sample_id", "condition"}
    ]
    summary_rows = []
    for region_id, region_table in sample_summary.groupby("region_id", sort=False):
        numerator_table = region_table[region_table["condition"].isin(contrast.numerator)]
        denominator_table = region_table[region_table["condition"].isin(contrast.denominator)]
        row = {"region_id": region_id}
        for feature_name in feature_columns:
            numerator_mean = numerator_table[feature_name].mean()
            denominator_mean = denominator_table[feature_name].mean()
            row[f"{feature_name}_numerator_mean"] = numerator_mean
            row[f"{feature_name}_denominator_mean"] = denominator_mean
            row[f"{feature_name}_delta_mean"] = numerator_mean - denominator_mean
        summary_rows.append(row)
    return evidence, pd.DataFrame(summary_rows)
```

- [ ] **Step 4: Add a failing built-in-default test, then implement it with a patched loader**

Add this test:

```python
def test_score_regions_single_read_window_features_uses_builtin_loader(monkeypatch):
    feature_table = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS", "read_id": "r1", "f0": 0.1},
            {"region_id": "reg1", "sample_id": "s2", "condition": "treated", "read_id": "r2", "f0": 0.9},
        ]
    )

    monkeypatch.setattr(
        region_contrasts,
        "_load_builtin_single_read_feature_table",
        lambda **kwargs: feature_table,
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions="regions.bed",
        contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="feature_summary_shift",
        motifs=["A,0"],
    )

    assert result.summary.iloc[0]["f0_delta_mean"] == pytest.approx(0.8)
```

- [ ] **Step 5: Run the feature tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -k "single_read_feature or single_read_window_features" -q
```

Expected: PASS.

- [ ] **Step 6: Commit the feature implementation**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add single-read feature contrasts"
```

## Task 7: Add `matched_pairwise` Coverage For Both Tracks

**Files:**
- Modify: `tests/test_region_contrasts.py`
- Modify: `dimelo/region_contrasts.py`

- [ ] **Step 1: Write the failing matched-pairwise tests**

Add these tests:

```python
def test_score_regions_single_read_mod_fraction_supports_matched_pairwise():
    evidence = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1_before", "condition": "before", "read_id": "r1", "modified_count": 1, "valid_count": 4, "read_mod_fraction": 0.25, "pair_id": "p1"},
            {"region_id": "reg1", "sample_id": "s1_after", "condition": "after", "read_id": "r2", "modified_count": 4, "valid_count": 4, "read_mod_fraction": 1.0, "pair_id": "p1"},
            {"region_id": "reg1", "sample_id": "s2_before", "condition": "before", "read_id": "r3", "modified_count": 2, "valid_count": 4, "read_mod_fraction": 0.5, "pair_id": "p2"},
            {"region_id": "reg1", "sample_id": "s2_after", "condition": "after", "read_id": "r4", "modified_count": 3, "valid_count": 4, "read_mod_fraction": 0.75, "pair_id": "p2"},
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["after"],
            denominator=["before"],
            pairing_key="pair_id",
        ),
        analysis_unit="single_read",
        representation="read_mod_fraction",
        signal_source="extract_reads",
        test="sample_distribution_shift",
        read_table=evidence,
    )

    assert result.summary.iloc[0]["delta_summary_mean"] > 0


def test_score_regions_single_read_window_features_supports_matched_pairwise():
    feature_table = pd.DataFrame(
        [
            {"region_id": "reg1", "sample_id": "s1_before", "condition": "before", "read_id": "r1", "pair_id": "p1", "f0": 0.1},
            {"region_id": "reg1", "sample_id": "s1_after", "condition": "after", "read_id": "r2", "pair_id": "p1", "f0": 0.9},
            {"region_id": "reg1", "sample_id": "s2_before", "condition": "before", "read_id": "r3", "pair_id": "p2", "f0": 0.2},
            {"region_id": "reg1", "sample_id": "s2_after", "condition": "after", "read_id": "r4", "pair_id": "p2", "f0": 0.8},
        ]
    )

    result = region_contrasts.score_regions(
        samples=[],
        regions=None,
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["after"],
            denominator=["before"],
            pairing_key="pair_id",
        ),
        analysis_unit="single_read",
        representation="read_window_features",
        signal_source="extract_features",
        test="feature_summary_shift",
        feature_table=feature_table,
    )

    assert result.summary.iloc[0]["f0_delta_mean"] > 0
```

- [ ] **Step 2: Implement matched-pairwise support conservatively**

Use sample-summary level tables only. Do not add pooled-read logic. The minimal first pass should:

- preserve `pairing_key` columns from the incoming read or feature table
- carry them through sample summaries
- compute deltas using only complete pairs

Add a narrow helper if needed:

```python
def _require_pairing_column(frame: pd.DataFrame, pairing_key: str) -> None:
    if pairing_key not in frame.columns:
        raise ValueError(
            f"matched_pairwise single_read scoring requires column {pairing_key!r}."
        )
```

- [ ] **Step 3: Run the matched-pairwise tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -k "matched_pairwise and single_read" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the matched-pairwise support**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add matched single-read contrasts"
```

## Task 8: Document The New `single_read` Contrast Family

**Files:**
- Modify: `docs/region-contrasts.md`
- Modify: `README.md`

- [ ] **Step 1: Update the region contrasts guide**

Add a short section to `docs/region-contrasts.md`:

```md
## Single-Read Contrasts

`region_contrasts` now supports `analysis_unit="single_read"` for defined-region
comparison on extract-backed read data.

First-version representations:

- `representation="read_mod_fraction"` with `signal_source="extract_reads"`
- `representation="read_window_features"` with `signal_source="extract_features"`

Supported contrast modes in v1:

- `pairwise`
- `group_vs_group`
- `matched_pairwise`

Important: reads remain observational units, but samples remain the inferential units.
The v1 single-read paths summarize reads within each `region x sample` before comparing
conditions.
```

- [ ] **Step 2: Add a short README note**

Add one short paragraph near the `region_contrasts` overview in `README.md`:

```md
`region_contrasts` also now supports `analysis_unit="single_read"` for extract-backed
defined-region comparison. Use `representation="read_mod_fraction"` for per-read motif
fractions or `representation="read_window_features"` for feature-vector summaries, with
sample-aware summaries rather than pooled-read inference.
```

- [ ] **Step 3: Commit the docs**

```bash
git add docs/region-contrasts.md README.md
git commit -m "docs: add single-read contrast guide"
```

## Task 9: Run Final Verification

**Files:**
- Reference: `dimelo/region_contrasts.py`
- Reference: `tests/test_region_contrasts.py`
- Reference: `docs/region-contrasts.md`
- Reference: `README.md`

- [ ] **Step 1: Run the full region contrast test file**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_region_contrasts.py -q
```

Expected: PASS.

- [ ] **Step 2: Run a targeted workflow/plotting regression subset**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/xdg-cache PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py tests/test_workflows.py -q
```

Expected: PASS with no new failures.

- [ ] **Step 3: Run a quick git status check**

Run:

```bash
git status --short
```

Expected: clean working tree.

- [ ] **Step 4: Self-review against the spec**

Check the completed work against `docs/superpowers/specs/2026-04-03-single-read-region-contrasts-design.md`:

- `analysis_unit="single_read"` added inside `region_contrasts`
- `read_mod_fraction` implemented
- `read_window_features` implemented
- built-in feature default and user override supported
- sample-aware summaries preserved
- v1 support remains bounded to `pairwise`, `group_vs_group`, and `matched_pairwise`
- no time-course or background-adjusted work accidentally added

- [ ] **Step 5: Commit any final touch-ups if needed**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py docs/region-contrasts.md README.md
git commit -m "fix: tighten single-read contrast support"
```

Only make this commit if verification or self-review uncovered a small follow-up patch.
