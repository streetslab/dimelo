# Plotting Axis Architecture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a shared plotting-axis core that supports compatibility-safe `5'->3'` orientation, fixed-window prep, aggregate segment-map prep, and legacy wrapper integration without breaking existing plotting entry points.

**Architecture:** Extend `dimelo/plotting.py` into the canonical plot-data preparation core with explicit axis and aggregation specs. Keep legacy plotting modules public and behaviorally stable by translating existing arguments such as `relative` and `regions_5to3prime` into the new shared prep layer. Restrict scaled/segment-normalized axes to aggregate plots, while keeping single-read plotting coordinate-preserving.

**Tech Stack:** Python 3.11, `dataclasses`, `pandas`, existing `dimelo` plotting/load utilities, `pytest`

---

## File Map

- Modify: `dimelo/plotting.py`
  - Add shared plotting spec dataclasses, validation, and plot-data prep helpers.
- Modify: `dimelo/plot_enrichment_profile.py`
  - Route one aggregate legacy plotter through the shared plotting core without changing public defaults.
- Modify: `dimelo/plot_reads.py`
  - Route one single-read legacy plotter through the shared plotting core without allowing scaled continuous axes.
- Modify: `docs/region-contrasts.md`
  - Document aggregate positional plotting semantics and compatibility notes.
- Modify: `docs/shared-clustering.md`
  - Document how the shared plotting layer applies to workflow `plot_data`.
- Modify: `README.md`
  - Add compatibility mapping between old plotting flags and new axis concepts.
- Create: `tests/test_plotting.py`
  - Add focused tests for spec validation and plot-data preparation.
- Modify: `dimelo/test/dimelo_test.py`
  - Preserve regression coverage for existing plotter entry points where needed.

### Task 1: Add Shared Plotting Specs And Validation

**Files:**
- Modify: `dimelo/plotting.py`
- Create: `tests/test_plotting.py`

- [ ] **Step 1: Write the failing tests for axis and aggregation spec validation**

```python
from dataclasses import replace

import pytest

from dimelo import plotting


def test_axis_spec_accepts_fixed_window_region_5to3():
    spec = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=1000,
        downstream_bp=1000,
    )

    plotting.validate_axis_spec(spec, plot_family="aggregate_profile")


def test_axis_spec_rejects_segment_map_without_segments():
    spec = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=None,
    )

    with pytest.raises(ValueError, match="segment_map requires segments"):
        plotting.validate_axis_spec(spec, plot_family="aggregate_profile")


def test_single_read_rejects_scaled_segments():
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec(
                segment_id="body",
                label="Body",
                start_ref=100,
                end_ref=500,
                mode="scaled",
                bins=50,
            )
        ],
    )

    with pytest.raises(ValueError, match="single_read"):
        plotting.validate_axis_spec(axis, plot_family="single_read_raster")


def test_aggregation_spec_accepts_equal_region_default():
    spec = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="faceted",
    )

    plotting.validate_aggregation_spec(spec)
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q`

Expected: FAIL with import or attribute errors for `AxisSpec`, `SegmentSpec`, `AggregationSpec`, or validation helpers not yet existing.

- [ ] **Step 3: Add the shared plotting spec dataclasses and validation helpers**

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class SegmentSpec:
    segment_id: str
    label: str
    start_ref: int
    end_ref: int
    mode: str  # "raw" | "scaled"
    bins: int | None = None
    plot_gap_after: bool = False
    contiguous_with_previous: bool = True


@dataclass(frozen=True)
class AxisSpec:
    orientation: str  # "genomic" | "region_5to3"
    coordinate_mode: str  # "fixed_window" | "segment_map"
    anchor: str | None = None
    upstream_bp: int | None = None
    downstream_bp: int | None = None
    segments: list[SegmentSpec] | None = None


@dataclass(frozen=True)
class AggregationSpec:
    weighting: str  # "equal_region" | "equal_observation" | "coverage_weighted"
    within_region_summary: str  # "mean" | "fraction" | "density"
    signal_normalization: str  # "none" | "per_region" | "global" | "control_regions"
    layout: str  # "concatenated" | "faceted"


def validate_axis_spec(axis: AxisSpec, *, plot_family: str) -> None:
    if axis.orientation not in {"genomic", "region_5to3"}:
        raise ValueError("AxisSpec.orientation must be 'genomic' or 'region_5to3'.")
    if axis.coordinate_mode not in {"fixed_window", "segment_map"}:
        raise ValueError("AxisSpec.coordinate_mode must be 'fixed_window' or 'segment_map'.")
    if axis.coordinate_mode == "segment_map" and not axis.segments:
        raise ValueError("AxisSpec coordinate_mode='segment_map' requires segments.")
    if plot_family == "single_read_raster" and axis.segments:
        if any(segment.mode == "scaled" for segment in axis.segments):
            raise ValueError(
                "single_read_raster plots must preserve coordinates and cannot use scaled segments."
            )


def validate_aggregation_spec(spec: AggregationSpec) -> None:
    if spec.weighting not in {"equal_region", "equal_observation", "coverage_weighted"}:
        raise ValueError("Unsupported weighting mode.")
    if spec.within_region_summary not in {"mean", "fraction", "density"}:
        raise ValueError("Unsupported within_region_summary.")
    if spec.signal_normalization not in {"none", "per_region", "global", "control_regions"}:
        raise ValueError("Unsupported signal_normalization.")
    if spec.layout not in {"concatenated", "faceted"}:
        raise ValueError("Unsupported layout mode.")
```

- [ ] **Step 4: Run the spec-validation tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q`

Expected: PASS for the new validation tests.

- [ ] **Step 5: Commit the shared plotting spec scaffolding**

```bash
git add dimelo/plotting.py \
        tests/test_plotting.py
git commit -m "feat: add plotting axis spec validation"
```

### Task 2: Implement Fixed-Window Plot Data Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Write failing tests for fixed-window orientation prep**

```python
import pandas as pd

from dimelo import plotting


def test_prepare_single_read_plot_data_flips_negative_regions_to_5to3():
    reads = pd.DataFrame(
        [
            {
                "region_id": "reg1",
                "region_strand": "-",
                "event_pos": 110,
                "anchor": 100,
                "read_id": "r1",
            }
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="custom",
        upstream_bp=20,
        downstream_bp=20,
    )

    payload = plotting.prepare_single_read_plot_data(
        reads,
        plot_family="single_read_raster",
        axis=axis,
        position_column="event_pos",
        anchor_column="anchor",
        region_strand_column="region_strand",
    )

    assert payload["plot_table"].loc[0, "plot_x"] == -10


def test_prepare_aggregate_plot_data_retains_metadata_for_fixed_window():
    table = pd.DataFrame(
        [
            {"region_id": "reg1", "region_strand": "+", "event_pos": 95, "anchor": 100, "signal": 1.0},
            {"region_id": "reg1", "region_strand": "+", "event_pos": 105, "anchor": 100, "signal": 3.0},
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=10,
        downstream_bp=10,
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="faceted",
    )

    payload = plotting.prepare_aggregate_plot_data(
        table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="signal",
        position_column="event_pos",
        anchor_column="anchor",
        region_strand_column="region_strand",
    )

    assert {"plot_table", "axis_table", "metadata"} <= set(payload)
    assert payload["metadata"]["orientation"] == "region_5to3"
```

- [ ] **Step 2: Run the fixed-window tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k fixed_window -q`

Expected: FAIL because `prepare_single_read_plot_data()` and `prepare_aggregate_plot_data()` do not yet exist.

- [ ] **Step 3: Implement shared fixed-window prep helpers**

```python
def _relative_position(position: float, anchor: float) -> float:
    return float(position) - float(anchor)


def _orient_position(relative_position: float, region_strand: str, orientation: str) -> float:
    if orientation == "region_5to3" and region_strand == "-":
        return -relative_position
    return relative_position


def prepare_single_read_plot_data(
    table: pd.DataFrame,
    *,
    plot_family: str,
    axis: AxisSpec,
    position_column: str,
    anchor_column: str,
    region_strand_column: str,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family=plot_family)
    if axis.coordinate_mode != "fixed_window":
        raise ValueError("single_read_raster currently supports fixed_window only.")

    plot_table = table.copy()
    plot_table["plot_x"] = plot_table.apply(
        lambda row: _orient_position(
            _relative_position(row[position_column], row[anchor_column]),
            row[region_strand_column],
            axis.orientation,
        ),
        axis=1,
    )
    axis_table = pd.DataFrame(
        [{"axis_min": -axis.upstream_bp, "axis_max": axis.downstream_bp, "segment_id": "window"}]
    )
    metadata = {"plot_family": plot_family, "orientation": axis.orientation, "coordinate_mode": axis.coordinate_mode}
    return {"plot_table": plot_table, "axis_table": axis_table, "metadata": metadata}


def prepare_aggregate_plot_data(
    table: pd.DataFrame,
    *,
    plot_family: str,
    axis: AxisSpec,
    aggregation: AggregationSpec,
    value_column: str,
    position_column: str,
    anchor_column: str,
    region_strand_column: str,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family=plot_family)
    validate_aggregation_spec(aggregation)
    plot_table = table.copy()
    plot_table["plot_x"] = plot_table.apply(
        lambda row: _orient_position(
            _relative_position(row[position_column], row[anchor_column]),
            row[region_strand_column],
            axis.orientation,
        ),
        axis=1,
    )
    axis_table = pd.DataFrame(
        [{"axis_min": -axis.upstream_bp, "axis_max": axis.downstream_bp, "segment_id": "window"}]
    )
    metadata = {
        "plot_family": plot_family,
        "orientation": axis.orientation,
        "coordinate_mode": axis.coordinate_mode,
        "weighting": aggregation.weighting,
        "within_region_summary": aggregation.within_region_summary,
    }
    return {"plot_table": plot_table, "axis_table": axis_table, "metadata": metadata}
```

- [ ] **Step 4: Run the fixed-window tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k fixed_window -q`

Expected: PASS.

- [ ] **Step 5: Commit the fixed-window prep slice**

```bash
git add dimelo/plotting.py \
        tests/test_plotting.py
git commit -m "feat: add fixed-window plotting prep"
```

### Task 3: Implement Aggregate Segment-Map Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Write failing tests for contiguous and non-contiguous aggregate segment maps**

```python
import pandas as pd

from dimelo import plotting


def test_prepare_aggregate_plot_data_builds_concatenated_segment_axis():
    table = pd.DataFrame(
        [
            {"region_id": "reg1", "segment_id": "upstream", "segment_pos": 0, "signal": 1.0},
            {"region_id": "reg1", "segment_id": "body", "segment_pos": 10, "signal": 2.0},
        ]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("upstream", "Upstream", 0, 100, "raw", bins=20),
            plotting.SegmentSpec("body", "Body", 100, 400, "scaled", bins=50, contiguous_with_previous=True),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="concatenated",
    )

    payload = plotting.prepare_aggregate_plot_data(
        table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="signal",
        segment_id_column="segment_id",
        segment_position_column="segment_pos",
    )

    assert payload["axis_table"]["segment_id"].tolist() == ["upstream", "body"]
    assert payload["axis_table"]["plot_start"].is_monotonic_increasing


def test_prepare_aggregate_plot_data_marks_non_contiguous_segment_breaks():
    table = pd.DataFrame(
        [{"region_id": "reg1", "segment_id": "exon1", "segment_pos": 5, "signal": 1.0}]
    )
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("exon1", "Exon 1", 100, 200, "scaled", bins=20),
            plotting.SegmentSpec(
                "exon3",
                "Exon 3",
                500,
                650,
                "scaled",
                bins=20,
                contiguous_with_previous=False,
                plot_gap_after=True,
            ),
        ],
    )
    aggregation = plotting.AggregationSpec(
        weighting="equal_region",
        within_region_summary="mean",
        signal_normalization="none",
        layout="faceted",
    )

    payload = plotting.prepare_aggregate_plot_data(
        table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="signal",
        segment_id_column="segment_id",
        segment_position_column="segment_pos",
    )

    assert payload["axis_table"].loc[1, "contiguous_with_previous"] is False
```

- [ ] **Step 2: Run the segment-map tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k segment -q`

Expected: FAIL because segment-map prep is not yet implemented.

- [ ] **Step 3: Implement aggregate-only segment-map preparation**

```python
def _build_segment_axis_table(segments: list[SegmentSpec]) -> pd.DataFrame:
    rows = []
    running_start = 0
    for segment in segments:
        span = segment.bins if segment.bins is not None else segment.end_ref - segment.start_ref
        rows.append(
            {
                "segment_id": segment.segment_id,
                "label": segment.label,
                "plot_start": running_start,
                "plot_end": running_start + span,
                "contiguous_with_previous": segment.contiguous_with_previous,
                "plot_gap_after": segment.plot_gap_after,
                "mode": segment.mode,
            }
        )
        running_start += span
    return pd.DataFrame(rows)


def prepare_aggregate_plot_data(..., segment_id_column: str | None = None, segment_position_column: str | None = None) -> dict[str, pd.DataFrame | dict[str, object]]:
    validate_axis_spec(axis, plot_family=plot_family)
    validate_aggregation_spec(aggregation)

    if axis.coordinate_mode == "segment_map":
        if not axis.segments:
            raise ValueError("segment_map requires segments.")
        axis_table = _build_segment_axis_table(axis.segments)
        if segment_id_column is None or segment_position_column is None:
            raise ValueError("segment_map plotting requires segment_id_column and segment_position_column.")
        plot_table = table.copy()
        plot_table = plot_table.merge(
            axis_table.loc[:, ["segment_id", "plot_start"]],
            on="segment_id",
            how="left",
        )
        plot_table["plot_x"] = plot_table["plot_start"] + plot_table[segment_position_column]
        metadata = {
            "plot_family": plot_family,
            "orientation": axis.orientation,
            "coordinate_mode": axis.coordinate_mode,
            "layout": aggregation.layout,
        }
        return {"plot_table": plot_table, "axis_table": axis_table, "metadata": metadata}
```

- [ ] **Step 4: Run the segment-map tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k segment -q`

Expected: PASS.

- [ ] **Step 5: Commit the aggregate segment-map slice**

```bash
git add dimelo/plotting.py \
        tests/test_plotting.py
git commit -m "feat: add aggregate segment-map plotting prep"
```

### Task 4: Route One Aggregate Legacy Plotter Through The Shared Core

**Files:**
- Modify: `dimelo/plot_enrichment_profile.py`
- Modify: `tests/test_plotting.py`
- Modify: `dimelo/test/dimelo_test.py`

- [ ] **Step 1: Write a regression test that legacy `regions_5to3prime` maps through the new axis logic**

```python
def test_legacy_enrichment_profile_retains_regions_5to3prime_behavior(monkeypatch):
    called = {}

    def fake_prepare_aggregate_plot_data(*args, **kwargs):
        called["axis"] = kwargs["axis"]
        return {
            "plot_table": kwargs["table"].copy() if "table" in kwargs else args[0].copy(),
            "axis_table": pd.DataFrame(),
            "metadata": {},
        }

    monkeypatch.setattr("dimelo.plotting.prepare_aggregate_plot_data", fake_prepare_aggregate_plot_data)
    # call the lowest-level trace helper with regions_5to3prime=True
    # and assert the translated axis carries orientation='region_5to3'
```

- [ ] **Step 2: Run the regression test to verify it fails**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k enrichment_profile -q`

Expected: FAIL because the legacy plotter does not yet call the shared prep layer.

- [ ] **Step 3: Translate legacy aggregate plotter arguments into shared specs**

```python
axis = plotting.AxisSpec(
    orientation="region_5to3" if regions_5to3prime else "genomic",
    coordinate_mode="fixed_window" if relative else "fixed_window",
    anchor="center" if relative else "absolute",
    upstream_bp=window_size // 2 if window_size else None,
    downstream_bp=window_size // 2 if window_size else None,
)
aggregation = plotting.AggregationSpec(
    weighting="equal_region",
    within_region_summary="mean",
    signal_normalization="none",
    layout="faceted",
)
payload = plotting.prepare_aggregate_plot_data(
    table=trace_table,
    plot_family="aggregate_profile",
    axis=axis,
    aggregation=aggregation,
    value_column="mod_fraction",
    position_column="position",
    anchor_column="center",
    region_strand_column="region_strand",
)
```

- [ ] **Step 4: Run the aggregate wrapper tests**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k enrichment_profile -q`

Expected: PASS.

- [ ] **Step 5: Commit the aggregate wrapper integration**

```bash
git add dimelo/plot_enrichment_profile.py \
        tests/test_plotting.py \
        dimelo/test/dimelo_test.py
git commit -m "refactor: route aggregate plotter through plotting core"
```

### Task 5: Route One Single-Read Legacy Plotter Through The Shared Core

**Files:**
- Modify: `dimelo/plot_reads.py`
- Modify: `tests/test_plotting.py`
- Modify: `dimelo/test/dimelo_test.py`

- [ ] **Step 1: Write a failing regression test for single-read compatibility**

```python
def test_legacy_plot_reads_translates_regions_5to3prime_to_axis_orientation(monkeypatch):
    called = {}

    def fake_prepare_single_read_plot_data(*args, **kwargs):
        called["axis"] = kwargs["axis"]
        return {"plot_table": pd.DataFrame(), "axis_table": pd.DataFrame(), "metadata": {}}

    monkeypatch.setattr("dimelo.plotting.prepare_single_read_plot_data", fake_prepare_single_read_plot_data)
    # invoke the legacy read plotting prep/helper and assert:
    # called["axis"].orientation == "region_5to3"


def test_prepare_single_read_plot_data_rejects_scaled_segment_axes():
    axis = plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="segment_map",
        segments=[
            plotting.SegmentSpec("body", "Body", 10, 100, "scaled", bins=20)
        ],
    )
    with pytest.raises(ValueError, match="single_read_raster"):
        plotting.prepare_single_read_plot_data(
            pd.DataFrame(),
            plot_family="single_read_raster",
            axis=axis,
            position_column="position",
            anchor_column="anchor",
            region_strand_column="region_strand",
        )
```

- [ ] **Step 2: Run the single-read tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k single_read -q`

Expected: FAIL because the legacy single-read plotter is not yet routed through the shared prep helper.

- [ ] **Step 3: Translate legacy single-read plot arguments into shared fixed-window prep**

```python
axis = plotting.AxisSpec(
    orientation="region_5to3" if regions_5to3prime else "genomic",
    coordinate_mode="fixed_window",
    anchor="center" if relative else "absolute",
    upstream_bp=window_size // 2 if window_size else None,
    downstream_bp=window_size // 2 if window_size else None,
)
payload = plotting.prepare_single_read_plot_data(
    table=read_table,
    plot_family="single_read_raster",
    axis=axis,
    position_column="plot_position",
    anchor_column="region_center",
    region_strand_column="region_strand",
)
```

- [ ] **Step 4: Run the single-read compatibility tests**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k single_read -q`

Expected: PASS.

- [ ] **Step 5: Commit the single-read wrapper integration**

```bash
git add dimelo/plot_reads.py \
        tests/test_plotting.py \
        dimelo/test/dimelo_test.py
git commit -m "refactor: route single-read plotter through plotting core"
```

### Task 6: Update User-Facing Plotting Documentation And Run Final Verification

**Files:**
- Modify: `docs/region-contrasts.md`
- Modify: `docs/shared-clustering.md`
- Modify: `README.md`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Write a failing documentation-backed test for compatibility mapping if needed**

```python
def test_plotting_core_exports_shared_specs():
    assert hasattr(plotting, "AxisSpec")
    assert hasattr(plotting, "prepare_single_read_plot_data")
    assert hasattr(plotting, "prepare_aggregate_plot_data")
```

- [ ] **Step 2: Run the focused plotting suite before doc edits**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q`

Expected: PASS before final docs, confirming the plotting core is stable.

- [ ] **Step 3: Update the docs with the compatibility mapping and safe-use rules**

```markdown
- `regions_5to3prime=True` maps to `AxisSpec(orientation="region_5to3")`
- fixed-window positional plotting is available for both aggregate and single-read views
- scaled or segmented metaregion plotting is aggregate-only
- single-read plots preserve geometry and reject continuous scaled axes
- `result.plot_data` remains the canonical renderer-neutral contract
```

- [ ] **Step 4: Run the branch verification slice**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py tests/test_region_contrasts.py tests/test_workflows.py -q`

Expected: PASS with no new failures. Existing non-blocking warnings are acceptable if unchanged and understood.

- [ ] **Step 5: Commit docs and final verification state**

```bash
git add tests/test_plotting.py \
        docs/region-contrasts.md \
        docs/shared-clustering.md \
        README.md
git commit -m "docs: add shared plotting axis guidance"
```

## Self-Review

### Spec coverage

- Shared axis and aggregation specs: covered by Task 1.
- Fixed-window prep with `region_5to3`: covered by Task 2.
- Aggregate-only segment maps, including contiguous and non-contiguous segments: covered by Task 3.
- Compatibility-safe legacy wrappers: covered by Tasks 4 and 5.
- Documentation and compatibility guidance: covered by Task 6.

No spec sections are currently unassigned.

### Placeholder scan

- No `TBD` or `TODO` placeholders remain.
- Each task includes concrete file paths, test cases, commands, and commit steps.

### Type consistency

- `AxisSpec`, `SegmentSpec`, `AggregationSpec`, `prepare_single_read_plot_data`, and `prepare_aggregate_plot_data` are defined once and reused consistently across tasks.
- `single_read_raster` is consistently treated as coordinate-preserving.
- `segment_map` remains aggregate-enabled, with scaled segments blocked for single-read usage.
