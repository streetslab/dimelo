# Paired Region Discovery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `dimelo.region_discovery.scan_genome(...)` with explicit paired discovery support for matched two-condition scans and ordered paired time-course scans.

**Architecture:** Keep the current tiled pileup-backed discovery path intact, then add a shared pairing layer inside `dimelo.region_discovery` that resolves complete matched units and builds an internal per-window paired table. Implement `matched_pairwise` first and `time_course` second, both returning the existing `RegionDiscoveryResult` shape with paired metadata and paired score columns.

**Tech Stack:** Python 3.11, pandas, pytest

---

## File Map

- `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/region_discovery.py`
  Add pairing-resolution helpers, paired-window table helpers, `matched_pairwise` scoring, paired `time_course` scoring, and paired metadata emission.
- `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_discovery.py`
  Add regression tests for complete-pair filtering, strict missing-pair errors, paired ranking columns, paired metadata, time-order validation, and downstream handoff.
- `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-discovery.md`
  Document when to use pooled discovery vs paired discovery, required sample metadata, and paired examples.
- `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md`
  Update the discovery guide link text to mention paired discovery support.

Scope note:

- This plan only extends `region_discovery`.
- It does not add paired inferential models.
- It does not change `parse_bam`, `global_analysis`, `region_contrasts`, or existing pooled discovery semantics.
- It keeps `RegionDiscoveryResult` data-first and plotting-optional.

### Task 1: Add Pair Resolution And Paired-Window Helpers

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/region_discovery.py`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_discovery.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_scan_genome_matched_pairwise_uses_only_complete_pairs(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_window_summary())

    result = region_discovery.scan_genome(
        samples=_paired_samplespecs(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=500,
        step_size=500,
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["targeting"],
            denominator=["nontargeting"],
            pairing_key="pair_id",
        ),
        score="effect_size_only",
    )

    assert result.metadata["pairing_policy"] == "complete_pairs_only"
    assert result.metadata["n_pairs_used"] == 2
    assert result.metadata["n_pairs_dropped"] == 1


def test_scan_genome_time_course_errors_on_missing_pairing_key(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_window_summary())

    with pytest.raises(ValueError, match="pairing_key"):
        region_discovery.scan_genome(
            samples=_paired_samplespecs(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=500,
            step_size=500,
            contrast=ContrastSpec(
                mode="time_course",
                time_order=["0min", "15min"],
            ),
            score="effect_size_only",
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
pytest tests/test_region_discovery.py::test_scan_genome_matched_pairwise_uses_only_complete_pairs tests/test_region_discovery.py::test_scan_genome_time_course_errors_on_missing_pairing_key -q
```

Expected: FAIL because `scan_genome(...)` still rejects paired contrast modes and has no paired metadata or pairing validation path.

- [ ] **Step 3: Write minimal implementation**

```python
def _pairing_policy_value(pairing_policy: str | None) -> str:
    return pairing_policy or "complete_pairs_only"


def _is_paired_contrast(contrast: ContrastSpec | None) -> bool:
    return contrast is not None and contrast.mode in {"matched_pairwise", "time_course"}


def _resolve_pair_ids(samples, pairing_key: str) -> dict[str, object]:
    pair_ids: dict[str, object] = {}
    for sample in samples:
        metadata = sample.metadata or {}
        if pairing_key not in metadata:
            raise ValueError(
                f"scan_genome paired discovery requires sample.metadata['{pairing_key}'] for every sample."
            )
        pair_ids[sample.sample_id] = metadata[pairing_key]
    return pair_ids


def _build_paired_window_table(
    summary: pd.DataFrame,
    *,
    samples,
    pairing_key: str,
    required_conditions: list[str],
    pairing_policy: str,
) -> tuple[pd.DataFrame, dict[str, int]]:
    sample_to_pair = _resolve_pair_ids(samples, pairing_key)
    paired = summary.copy()
    paired["pair_id"] = paired["sample_id"].map(sample_to_pair)
    paired = paired.dropna(subset=["pair_id"])

    present = (
        paired.loc[:, ["pair_id", "condition"]]
        .drop_duplicates()
        .groupby("pair_id")["condition"]
        .agg(lambda values: set(values))
    )
    complete_pair_ids = sorted(
        pair_id for pair_id, values in present.items() if set(required_conditions).issubset(values)
    )
    dropped_pair_count = int(len(present) - len(complete_pair_ids))

    if pairing_policy == "error_on_missing" and dropped_pair_count:
        raise ValueError("scan_genome paired discovery found incomplete matched units.")
    if not complete_pair_ids:
        raise ValueError("scan_genome paired discovery found no complete matched units.")

    paired = paired.loc[paired["pair_id"].isin(complete_pair_ids)].copy()
    aggregated = (
        paired.groupby(_WINDOW_KEY_COLUMNS + ["pair_id", "condition"], as_index=False, sort=False)
        .agg(modified_count=("modified_count", "sum"), valid_count=("valid_count", "sum"))
        .copy()
    )
    aggregated["window_fraction"] = _safe_fraction(
        aggregated["modified_count"], aggregated["valid_count"]
    )
    return aggregated, {
        "n_pairs_used": len(complete_pair_ids),
        "n_pairs_dropped": dropped_pair_count,
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run:

```bash
pytest tests/test_region_discovery.py::test_scan_genome_matched_pairwise_uses_only_complete_pairs tests/test_region_discovery.py::test_scan_genome_time_course_errors_on_missing_pairing_key -q
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/test_region_discovery.py dimelo/region_discovery.py
git commit -m "feat: add paired region discovery pair resolution"
```

### Task 2: Add `matched_pairwise` Scoring And Paired Metadata

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/region_discovery.py`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_discovery.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_scan_genome_matched_pairwise_ranks_by_mean_abs_delta(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_pairwise_window_summary())

    result = region_discovery.scan_genome(
        samples=_paired_samplespecs(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=500,
        step_size=500,
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["targeting"],
            denominator=["nontargeting"],
            pairing_key="pair_id",
        ),
        score="effect_size_only",
    )

    assert list(result.hits["window_id"]) == ["chr1:0-500", "chr1:500-1000"]
    assert result.hits.loc[0, "mean_delta"] == pytest.approx(0.35)
    assert result.hits.loc[0, "mean_abs_delta"] == pytest.approx(0.35)
    assert result.hits.loc[0, "sign_agreement"] == pytest.approx(1.0)
    assert result.metadata["paired_mode"] == "matched_pairwise"
    assert result.metadata["rank_by"] == "mean_abs_delta"


def test_scan_genome_matched_pairwise_errors_on_missing_pairs_in_strict_mode(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_window_summary())

    with pytest.raises(ValueError, match="incomplete matched units"):
        region_discovery.scan_genome(
            samples=_paired_samplespecs(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=500,
            step_size=500,
            contrast=ContrastSpec(
                mode="matched_pairwise",
                numerator=["targeting"],
                denominator=["nontargeting"],
                pairing_key="pair_id",
            ),
            score="effect_size_only",
            pairing_policy="error_on_missing",
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
pytest tests/test_region_discovery.py::test_scan_genome_matched_pairwise_ranks_by_mean_abs_delta tests/test_region_discovery.py::test_scan_genome_matched_pairwise_errors_on_missing_pairs_in_strict_mode -q
```

Expected: FAIL because paired windows are not yet summarized into `mean_delta`, `mean_abs_delta`, `delta_sd`, or `sign_agreement`.

- [ ] **Step 3: Write minimal implementation**

```python
def _score_matched_pairwise(
    paired_window_table: pd.DataFrame,
    *,
    contrast: ContrastSpec,
) -> pd.DataFrame:
    numerator = paired_window_table.loc[
        paired_window_table["condition"].isin(contrast.numerator)
    ].copy()
    denominator = paired_window_table.loc[
        paired_window_table["condition"].isin(contrast.denominator)
    ].copy()

    numerator = numerator.rename(
        columns={
            "modified_count": "numerator_modified_count",
            "valid_count": "numerator_valid_count",
            "window_fraction": "numerator_fraction",
        }
    )
    denominator = denominator.rename(
        columns={
            "modified_count": "denominator_modified_count",
            "valid_count": "denominator_valid_count",
            "window_fraction": "denominator_fraction",
        }
    )

    merged = numerator.merge(
        denominator[
            _WINDOW_KEY_COLUMNS
            + ["pair_id", "denominator_modified_count", "denominator_valid_count", "denominator_fraction"]
        ],
        on=_WINDOW_KEY_COLUMNS + ["pair_id"],
        how="inner",
    )
    merged["delta"] = merged["numerator_fraction"] - merged["denominator_fraction"]

    scored = (
        merged.groupby(_WINDOW_KEY_COLUMNS, as_index=False, sort=False)
        .agg(
            mean_delta=("delta", "mean"),
            mean_abs_delta=("delta", lambda values: float(values.abs().mean())),
            delta_sd=("delta", lambda values: float(values.std(ddof=0))),
            sign_agreement=(
                "delta",
                lambda values: float(
                    max((values.gt(0)).mean(), (values.lt(0)).mean()) if len(values) else 0.0
                ),
            ),
            n_pairs_used=("pair_id", "nunique"),
        )
        .copy()
    )
    scored["score_value"] = scored["mean_abs_delta"]
    scored["p_value"] = pd.NA
    scored["adjusted_p_value"] = pd.NA
    return scored


def scan_genome(..., pairing_policy: str | None = None, rank_by: str | None = None, ...):
    ...
    if contrast is not None and contrast.mode == "matched_pairwise":
        paired_table, pairing_meta = _build_paired_window_table(
            window_summary,
            samples=samples,
            pairing_key=contrast.pairing_key,
            required_conditions=list(contrast.numerator or []) + list(contrast.denominator or []),
            pairing_policy=_pairing_policy_value(pairing_policy),
        )
        ranked = _rank_windows(
            _score_matched_pairwise(paired_table, contrast=contrast),
            score="effect_size_only",
            primary_score_column="mean_abs_delta",
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run:

```bash
pytest tests/test_region_discovery.py::test_scan_genome_matched_pairwise_ranks_by_mean_abs_delta tests/test_region_discovery.py::test_scan_genome_matched_pairwise_errors_on_missing_pairs_in_strict_mode -q
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/test_region_discovery.py dimelo/region_discovery.py
git commit -m "feat: add matched pairwise region discovery scoring"
```

### Task 3: Add Paired Ordered `time_course` Discovery

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/region_discovery.py`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_discovery.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_scan_genome_time_course_ranks_by_trajectory_amplitude_mean(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_time_course_window_summary())

    result = region_discovery.scan_genome(
        samples=_paired_time_course_samplespecs(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=500,
        step_size=500,
        contrast=ContrastSpec(
            mode="time_course",
            time_order=["0min", "15min", "30min"],
            pairing_key="pair_id",
        ),
        score="effect_size_only",
    )

    assert result.hits.loc[0, "trajectory_amplitude_mean"] == pytest.approx(0.55)
    assert result.hits.loc[0, "trajectory_amplitude_sd"] == pytest.approx(0.05)
    assert result.metadata["paired_mode"] == "time_course"
    assert result.metadata["time_order"] == ["0min", "15min", "30min"]
    assert result.metadata["rank_by"] == "trajectory_amplitude_mean"


def test_scan_genome_time_course_errors_when_time_order_conditions_are_missing(monkeypatch):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_time_course_window_summary())

    with pytest.raises(ValueError, match="time_order"):
        region_discovery.scan_genome(
            samples=_paired_time_course_samplespecs(),
            motifs=["A,0"],
            genome_sizes={"chr1": 1000},
            window_size=500,
            step_size=500,
            contrast=ContrastSpec(
                mode="time_course",
                time_order=["0min", "15min", "60min"],
                pairing_key="pair_id",
            ),
            score="effect_size_only",
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
pytest tests/test_region_discovery.py::test_scan_genome_time_course_ranks_by_trajectory_amplitude_mean tests/test_region_discovery.py::test_scan_genome_time_course_errors_when_time_order_conditions_are_missing -q
```

Expected: FAIL because `time_course` is still rejected by `scan_genome(...)`.

- [ ] **Step 3: Write minimal implementation**

```python
def _validate_time_order(paired_window_table: pd.DataFrame, time_order: list[str]) -> None:
    available = set(paired_window_table["condition"].dropna().tolist())
    missing = [condition for condition in time_order if condition not in available]
    if missing:
        raise ValueError(
            "scan_genome paired time_course requested missing time_order condition(s): "
            + ", ".join(missing)
        )


def _score_paired_time_course(
    paired_window_table: pd.DataFrame,
    *,
    time_order: list[str],
) -> pd.DataFrame:
    _validate_time_order(paired_window_table, time_order)
    ordered = paired_window_table.copy()
    ordered["condition"] = pd.Categorical(
        ordered["condition"],
        categories=time_order,
        ordered=True,
    )
    ordered = ordered.sort_values(_WINDOW_KEY_COLUMNS + ["pair_id", "condition"])

    per_pair = (
        ordered.groupby(_WINDOW_KEY_COLUMNS + ["pair_id"], as_index=False, sort=False)
        .agg(
            trajectory_amplitude=("window_fraction", lambda values: float(values.max() - values.min())),
        )
        .copy()
    )
    scored = (
        per_pair.groupby(_WINDOW_KEY_COLUMNS, as_index=False, sort=False)
        .agg(
            trajectory_amplitude_mean=("trajectory_amplitude", "mean"),
            trajectory_amplitude_median=("trajectory_amplitude", "median"),
            trajectory_amplitude_sd=("trajectory_amplitude", lambda values: float(values.std(ddof=0))),
            n_pairs_used=("pair_id", "nunique"),
        )
        .copy()
    )
    scored["score_value"] = scored["trajectory_amplitude_mean"]
    scored["p_value"] = pd.NA
    scored["adjusted_p_value"] = pd.NA
    return scored
```

- [ ] **Step 4: Run test to verify it passes**

Run:

```bash
pytest tests/test_region_discovery.py::test_scan_genome_time_course_ranks_by_trajectory_amplitude_mean tests/test_region_discovery.py::test_scan_genome_time_course_errors_when_time_order_conditions_are_missing -q
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/test_region_discovery.py dimelo/region_discovery.py
git commit -m "feat: add paired time-course region discovery"
```

### Task 4: Add Downstream Handoff Coverage And User Docs

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_discovery.py`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-discovery.md`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md`

- [ ] **Step 1: Write the failing tests**

```python
def test_paired_discovery_hits_to_bed_and_region_contrasts_handoff(monkeypatch, tmp_path):
    monkeypatch.setattr(global_analysis, "build_window_summary", lambda **_: _mock_paired_pairwise_window_summary())

    discovery = region_discovery.scan_genome(
        samples=_paired_samplespecs(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=500,
        step_size=500,
        contrast=ContrastSpec(
            mode="matched_pairwise",
            numerator=["targeting"],
            denominator=["nontargeting"],
            pairing_key="pair_id",
        ),
        score="effect_size_only",
    )

    bed = region_discovery.hits_to_bed(discovery.hits.head(1))
    assert list(bed.columns) == ["chrom", "start", "end", "name", "score", "strand"]
    assert bed.loc[0, "name"] == "chr1:0-500"
```

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
pytest tests/test_region_discovery.py::test_paired_discovery_hits_to_bed_and_region_contrasts_handoff -q
```

Expected: FAIL because the paired hit columns and metadata are not yet exercised through the downstream handoff path.

- [ ] **Step 3: Implement wrapper polish and docs**

```python
metadata.update(
    {
        "pairing_key": contrast.pairing_key,
        "pairing_policy": active_pairing_policy,
        "n_pairs_used": pairing_meta["n_pairs_used"],
        "n_pairs_dropped": pairing_meta["n_pairs_dropped"],
        "paired_mode": contrast.mode,
        "rank_by": active_rank_by,
    }
)
if contrast.mode == "time_course":
    metadata["time_order"] = list(contrast.time_order or [])
```

```markdown
## Paired discovery

Use paired discovery when your samples are matched across conditions and pooling would hide the design:

- nontargeting versus targeting with matched replicates
- before/after on the same sample
- ordered time courses across the same matched units

`scan_genome(...)` now supports:

- `ContrastSpec(mode="matched_pairwise", ...)`
- `ContrastSpec(mode="time_course", ...)`

For paired modes, every `SampleSpec.metadata` entry must include the requested `pairing_key`.
```

- [ ] **Step 4: Run test to verify it passes**

Run:

```bash
pytest tests/test_region_discovery.py tests/test_global_analysis.py tests/test_region_contrasts.py -q
```

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add tests/test_region_discovery.py docs/region-discovery.md README.md dimelo/region_discovery.py
git commit -m "docs: add paired region discovery guide"
```

---

## Verification And Review Gates

- After each task:
  - Run the task-focused pytest subset.
  - Confirm existing pooled discovery tests still pass.
  - Review for continuity with current `scan_genome(...)` behavior.
- After Task 4:
  - Run:

```bash
pytest tests/test_region_discovery.py tests/test_global_analysis.py tests/test_region_contrasts.py tests/test_models.py tests/test_workflows.py -q
```

- If the full subset fails:
  - fix regressions before proceeding
  - do not merge a paired mode that changes current pooled semantics

## V1 Defaults To Preserve

- pooled `pairwise` and `group_vs_group` behavior stays unchanged
- paired discovery requires explicit `pairing_key`
- default `pairing_policy="complete_pairs_only"`
- strict `pairing_policy="error_on_missing"` is opt-in
- `matched_pairwise` ranks by `mean_abs_delta`
- paired `time_course` ranks by `trajectory_amplitude_mean`
- no silent fallback from paired discovery to pooled discovery
- `RegionDiscoveryResult` remains data-first with optional figures
