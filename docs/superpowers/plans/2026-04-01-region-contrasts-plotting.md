# Region Contrasts Plotting Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add fixed-window, data-first plotting helpers for `dimelo.region_contrasts` that prepare profile and heatmap payloads from a `RegionContrastResult` plus an explicit positional source table.

**Architecture:** Extend [plotting.py](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py) with a small internal join/validation layer for contrast results and positional source tables, then expose two public helpers: one for aggregate profiles and one for heatmaps. Keep `region_contrasts` unchanged; plotting prep consumes existing `RegionContrastResult` tables and parse/load-derived positional summaries, uses the shared `AxisSpec` / `AggregationSpec` fixed-window path, and returns renderer-neutral payloads.

**Tech Stack:** Python 3.11, pandas, pytest, existing `dimelo.plotting` shared axis helpers, existing `RegionContrastResult` models

---

## File Structure

- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py`
  - Add internal helpers for validating/joining `RegionContrastResult` against a positional source table.
  - Add `prepare_region_contrast_profile_data(...)`.
  - Add `prepare_region_contrast_heatmap_data(...)`.
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py`
  - Add regression tests for the new plot-prep helpers.
  - Cover both numerator/denominator/delta payload generation and failure modes.
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-contrasts.md`
  - Document the new helper APIs and the required `position_table` contract.
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md`
  - Add a short compatibility note pointing users from `region_contrasts` into the new plotting helpers.

### Task 1: Add minimal failing tests for contrast-to-position joins

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py`
- Reference: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/models.py`
- Reference: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py`

- [ ] **Step 1: Write the failing tests for a minimal profile payload**

Add these tests near the other plotting helper tests in `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py`:

```python
from dimelo.models import ContrastSpec, RegionContrastResult


def _minimal_region_contrast_result() -> RegionContrastResult:
    regions = pd.DataFrame(
        [
            {"region_id": "chr1:90-110,+", "condition": "NS", "fraction": 0.20, "rank": 2},
            {"region_id": "chr1:90-110,+", "condition": "15min", "fraction": 0.55, "rank": 2},
            {"region_id": "chr1:190-210,-", "condition": "NS", "fraction": 0.30, "rank": 1},
            {"region_id": "chr1:190-210,-", "condition": "15min", "fraction": 0.70, "rank": 1},
        ]
    )
    summary = pd.DataFrame(
        [
            {
                "region_id": "chr1:90-110,+",
                "fraction": 0.55,
                "reference_fraction": 0.20,
                "delta_fraction": 0.35,
                "rank": 2,
            },
            {
                "region_id": "chr1:190-210,-",
                "fraction": 0.70,
                "reference_fraction": 0.30,
                "delta_fraction": 0.40,
                "rank": 1,
            },
        ]
    )
    return RegionContrastResult(
        regions=regions,
        summary=summary,
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["15min"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        metadata={
            "analysis_unit": "ensemble_region",
            "representation": "modified_fraction",
            "signal_source": "pileup_counts",
            "test": "effect_size_only",
        },
        plot_data={},
    )


def test_prepare_region_contrast_profile_data_returns_all_value_modes():
    result = _minimal_region_contrast_result()
    position_table = pd.DataFrame(
        [
            {"region_id": "chr1:90-110,+", "condition": "NS", "position": 95, "anchor": 100, "value": 0.1, "region_strand": "+"},
            {"region_id": "chr1:90-110,+", "condition": "15min", "position": 95, "anchor": 100, "value": 0.6, "region_strand": "+"},
            {"region_id": "chr1:190-210,-", "condition": "NS", "position": 205, "anchor": 200, "value": 0.2, "region_strand": "-"},
            {"region_id": "chr1:190-210,-", "condition": "15min", "position": 205, "anchor": 200, "value": 0.8, "region_strand": "-"},
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=20,
        downstream_bp=20,
    )
    aggregation = plotting.AggregationSpec()

    payload = plotting.prepare_region_contrast_profile_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )

    assert set(payload["plot_table"]["value_mode"]) == {"numerator", "denominator", "delta"}
    assert payload["metadata"]["plot_family"] == "region_contrast_profile"
```

- [ ] **Step 2: Add the failing heatmap and validation tests**

Append these tests in the same file:

```python
def test_prepare_region_contrast_heatmap_data_orders_rows_by_rank():
    result = _minimal_region_contrast_result()
    position_table = pd.DataFrame(
        [
            {"region_id": "chr1:90-110,+", "condition": "NS", "position": 95, "anchor": 100, "value": 0.1, "region_strand": "+"},
            {"region_id": "chr1:90-110,+", "condition": "15min", "position": 95, "anchor": 100, "value": 0.6, "region_strand": "+"},
            {"region_id": "chr1:190-210,-", "condition": "NS", "position": 205, "anchor": 200, "value": 0.2, "region_strand": "-"},
            {"region_id": "chr1:190-210,-", "condition": "15min", "position": 205, "anchor": 200, "value": 0.8, "region_strand": "-"},
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=20,
        downstream_bp=20,
    )
    aggregation = plotting.AggregationSpec()

    payload = plotting.prepare_region_contrast_heatmap_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode="all",
    )

    rank_rows = payload["plot_table"].loc[:, ["region_id", "row_order"]].drop_duplicates()
    assert list(rank_rows.sort_values("row_order")["region_id"]) == [
        "chr1:190-210,-",
        "chr1:90-110,+",
    ]
    assert payload["metadata"]["plot_family"] == "region_contrast_heatmap"


def test_prepare_region_contrast_profile_data_requires_joinable_grouping_key():
    result = _minimal_region_contrast_result()
    position_table = pd.DataFrame(
        [
            {"region_id": "chr1:90-110,+", "position": 95, "anchor": 100, "value": 0.1, "region_strand": "+"}
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=20,
        downstream_bp=20,
    )
    aggregation = plotting.AggregationSpec()

    with pytest.raises(ValueError, match="sample_id or condition"):
        plotting.prepare_region_contrast_profile_data(
            result=result,
            position_table=position_table,
            axis=axis,
            aggregation=aggregation,
            value_mode="all",
        )
```

- [ ] **Step 3: Run the new tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py -k "region_contrast_profile_data or region_contrast_heatmap_data" -q
```

Expected: FAIL with `AttributeError` or `NameError` because the new plotting helpers do not exist yet.

- [ ] **Step 4: Commit the failing-test scaffold**

```bash
git add /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py
git commit -m "test: add region contrast plotting coverage"
```

### Task 2: Implement shared join and normalization helpers

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py`
- Test: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py`

- [ ] **Step 1: Add a contrast-result metadata reader and join-key validator**

Insert small internal helpers in `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py` below the existing validation helpers:

```python
def _region_contrast_grouping_key(result, position_table: pd.DataFrame) -> str:
    available = [column for column in ("sample_id", "condition") if column in position_table.columns]
    if not available:
        raise ValueError(
            "region contrast plotting requires position_table to include sample_id or condition."
        )
    if "condition" in available:
        return "condition"
    return "sample_id"


def _region_contrast_metadata(result) -> dict[str, object]:
    return {
        "analysis_unit": result.metadata.get("analysis_unit"),
        "representation": result.metadata.get("representation"),
        "signal_source": result.metadata.get("signal_source"),
        "test": result.metadata.get("test"),
    }
```

- [ ] **Step 2: Add a helper that validates and filters a positional source table**

Add this helper in the same file:

```python
def _prepare_region_contrast_position_table(
    *,
    result,
    position_table: pd.DataFrame,
    grouping_key: str,
) -> pd.DataFrame:
    _require_columns(
        position_table,
        ("region_id", grouping_key, "position", "value", "region_strand"),
        "position_table",
    )
    filtered = position_table[position_table["region_id"].isin(result.summary["region_id"])].copy()
    if filtered.empty:
        raise ValueError("position_table does not contain any region_id values present in the contrast result.")
    return filtered.reset_index(drop=True)
```

- [ ] **Step 3: Add a side-merging helper for numerator, denominator, and delta**

Add this helper in the same file:

```python
def _prepare_region_contrast_value_modes(
    *,
    result,
    position_table: pd.DataFrame,
    grouping_key: str,
) -> pd.DataFrame:
    contrast = result.contrast
    numerator_mask = position_table[grouping_key].isin(contrast.numerator or [])
    denominator_mask = position_table[grouping_key].isin(contrast.denominator or [])

    numerator = position_table.loc[numerator_mask].copy()
    denominator = position_table.loc[denominator_mask].copy()
    if numerator.empty or denominator.empty:
        raise ValueError("position_table does not contain rows for both contrast sides.")

    numerator["value_mode"] = "numerator"
    denominator["value_mode"] = "denominator"

    delta = numerator.merge(
        denominator,
        on=["region_id", "position", "region_strand"],
        suffixes=("_numerator", "_denominator"),
        how="inner",
    )
    if delta.empty:
        raise ValueError("Unable to compute delta because numerator and denominator positions do not align.")

    delta = delta.rename(
        columns={
            "position": "position",
            "region_strand": "region_strand",
        }
    )
    delta["value"] = delta["value_numerator"] - delta["value_denominator"]
    delta[grouping_key] = "delta"
    delta["value_mode"] = "delta"
    delta = delta.loc[:, ["region_id", grouping_key, "position", "value", "region_strand"]]

    return pd.concat(
        [
            numerator.loc[:, ["region_id", grouping_key, "position", "value", "region_strand", "value_mode"]],
            denominator.loc[:, ["region_id", grouping_key, "position", "value", "region_strand", "value_mode"]],
            delta,
        ],
        ignore_index=True,
    )
```

- [ ] **Step 4: Run the targeted tests to verify the join layer works**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py -k "region_contrast_profile_data_requires_joinable_grouping_key" -q
```

Expected: PASS for the validation test, while the helper-existence tests still fail.

- [ ] **Step 5: Commit the join/validation layer**

```bash
git add /Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py
git commit -m "feat: add region contrast plotting join helpers"
```

### Task 3: Implement `prepare_region_contrast_profile_data(...)`

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py`
- Test: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py`

- [ ] **Step 1: Add the public helper using the shared fixed-window prep**

Add this function in `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py` near the other public plotting helpers:

```python
def prepare_region_contrast_profile_data(
    *,
    result,
    position_table: pd.DataFrame,
    axis: AxisSpec,
    aggregation: AggregationSpec,
    value_mode: str = "all",
) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family="aggregate_profile")
    validate_aggregation_spec(aggregation)
    grouping_key = _region_contrast_grouping_key(result, position_table)
    filtered = _prepare_region_contrast_position_table(
        result=result,
        position_table=position_table,
        grouping_key=grouping_key,
    )
    valued = _prepare_region_contrast_value_modes(
        result=result,
        position_table=filtered,
        grouping_key=grouping_key,
    )
    if value_mode != "all":
        valued = valued.loc[valued["value_mode"] == value_mode].copy()

    valued["anchor"] = valued["position"]
    prepared = prepare_aggregate_plot_data(
        valued,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="value",
        position_column="position",
        anchor_column="anchor",
        region_strand_column="region_strand",
    )
    prepared["metadata"] = {
        **prepared["metadata"],
        **_region_contrast_metadata(result),
        "plot_family": "region_contrast_profile",
        "value_modes": sorted(valued["value_mode"].unique()),
    }
    prepared["summary_table"] = (
        prepared["plot_table"]
        .groupby(["value_mode", "plot_x"], as_index=False)["value"]
        .mean()
    )
    return prepared
```

- [ ] **Step 2: Fix the anchor handling so the profile prep uses the real fixed-window contract**

Replace the temporary `valued["anchor"] = valued["position"]` line with:

```python
_require_columns(filtered, ("anchor",), "position_table")
valued = valued.merge(
    filtered.loc[:, ["region_id", grouping_key, "position", "anchor"]].drop_duplicates(),
    on=["region_id", grouping_key, "position"],
    how="left",
)
```

This keeps the helper aligned with the actual fixed-window substrate instead of collapsing all `plot_x` values to zero.

- [ ] **Step 3: Run the targeted profile tests**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py -k "prepare_region_contrast_profile_data" -q
```

Expected: PASS for the profile tests.

- [ ] **Step 4: Commit the profile helper**

```bash
git add /Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py
git commit -m "feat: add region contrast profile plotting prep"
```

### Task 4: Implement `prepare_region_contrast_heatmap_data(...)`

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py`
- Test: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py`

- [ ] **Step 1: Add the heatmap helper on top of the same substrate**

Add this function in `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py`:

```python
def prepare_region_contrast_heatmap_data(
    *,
    result,
    position_table: pd.DataFrame,
    axis: AxisSpec,
    aggregation: AggregationSpec,
    value_mode: str = "all",
) -> dict[str, pd.DataFrame | dict[str, object]]:
    profile_payload = prepare_region_contrast_profile_data(
        result=result,
        position_table=position_table,
        axis=axis,
        aggregation=aggregation,
        value_mode=value_mode,
    )
    rank_lookup = (
        result.summary.loc[:, ["region_id", "rank"]]
        .drop_duplicates()
        .sort_values(["rank", "region_id"], kind="stable")
        .reset_index(drop=True)
    )
    rank_lookup["row_order"] = range(len(rank_lookup))
    plot_table = profile_payload["plot_table"].merge(rank_lookup, on="region_id", how="left")
    profile_payload["plot_table"] = plot_table
    profile_payload["metadata"] = {
        **profile_payload["metadata"],
        "plot_family": "region_contrast_heatmap",
    }
    profile_payload["summary_table"] = rank_lookup
    return profile_payload
```

- [ ] **Step 2: Add a failure-mode test for missing contrast-side rows**

Extend `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py` with:

```python
def test_prepare_region_contrast_heatmap_data_requires_both_contrast_sides():
    result = _minimal_region_contrast_result()
    position_table = pd.DataFrame(
        [
            {"region_id": "chr1:90-110,+", "condition": "NS", "position": 95, "anchor": 100, "value": 0.1, "region_strand": "+"}
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=20,
        downstream_bp=20,
    )
    aggregation = plotting.AggregationSpec()

    with pytest.raises(ValueError, match="both contrast sides"):
        plotting.prepare_region_contrast_heatmap_data(
            result=result,
            position_table=position_table,
            axis=axis,
            aggregation=aggregation,
            value_mode="all",
        )
```

- [ ] **Step 3: Run the full plotting suite**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py -q
```

Expected: PASS, including the new heatmap coverage.

- [ ] **Step 4: Commit the heatmap helper**

```bash
git add /Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py
git commit -m "feat: add region contrast heatmap plotting prep"
```

### Task 5: Document the new plotting helpers and verify the slice

**Files:**
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-contrasts.md`
- Modify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md`
- Verify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py`
- Verify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_contrasts.py`
- Verify: `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_workflows.py`

- [ ] **Step 1: Add a region-contrasts plotting example**

Append a short section to `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-contrasts.md`:

```python
profile_payload = plotting.prepare_region_contrast_profile_data(
    result=result,
    position_table=position_table,
    axis=plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=2000,
        downstream_bp=2000,
    ),
    aggregation=plotting.AggregationSpec(),
    value_mode="all",
)
```

Add one sentence clarifying that `position_table` must already be an aggregated positional substrate from the parsing/loading layer.

- [ ] **Step 2: Add a short README compatibility note**

Add one short paragraph to `/Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md` under the plotting section explaining:

```text
Newer region-contrast plotting helpers accept a RegionContrastResult plus an explicit position_table from the parsing/loading layer. This keeps region scoring and positional extraction separate while still using the shared plotting-axis system.
```

- [ ] **Step 3: Run the final verification set**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_region_contrasts.py /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_workflows.py -q
```

Expected: PASS. If warnings appear, note them in the summary only if they are new or blocking.

- [ ] **Step 4: Commit the docs and final verification slice**

```bash
git add /Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/region-contrasts.md /Users/ngamarra/Documents/GitHub/dimelo-toolkit/README.md /Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py /Users/ngamarra/Documents/GitHub/dimelo-toolkit/tests/test_plotting.py
git commit -m "docs: add region contrast plotting guide"
```

## Self-Review

- Spec coverage:
  - public helpers: Task 3 and Task 4
  - explicit `position_table`: Task 1 and Task 2
  - numerator/denominator/delta payloads: Task 2 and Task 3
  - fixed-window aggregate-only scope: Task 3 and Task 4
  - compatibility/docs: Task 5
- Placeholder scan:
  - no `TODO`, `TBD`, or deferred implementation steps are left in the task bodies
- Type consistency:
  - `prepare_region_contrast_profile_data(...)` and `prepare_region_contrast_heatmap_data(...)` are used consistently across tests, code steps, and docs

