# Region Contrasts Plotting Design

## Purpose

This design defines the first analysis-specific plotting helpers built on top of the shared plotting-axis architecture. The goal is to support positional plotting for `dimelo.region_contrasts` while preserving the package's current separation of responsibilities:

- `region_contrasts` scores known loci
- parse/load workflows provide positional substrates
- plotting helpers join scored loci to positional summaries and return data-first plot payloads

This is an additive design. It must not replace existing plotting entry points or change the semantics of the current parsing layer.

## Design Boundary

The first plotting slice is for `region_contrasts` only.

It should support:

- aggregate positional profiles across known loci
- per-region positional heatmaps
- numerator, denominator, and delta payloads in one prepared result
- both pileup-backed abundance contrasts and cluster-occupancy follow-up, as long as a compatible positional source table is supplied

It should not require `RegionContrastResult` to carry heavy positional tables internally.

## Core Rule

`RegionContrastResult` is the scored-region contract.

Positional information comes from an explicit `position_table` derived from parsing/loading workflows. Plotting helpers merge the two.

This keeps responsibilities clean:

- scoring stays in `dimelo.region_contrasts`
- positional extraction stays in the existing parse/load layer
- plotting prep stays in `dimelo.plotting`

## Public API

Add two public helpers in [plotting.py](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/dimelo/plotting.py):

- `prepare_region_contrast_profile_data(...)`
- `prepare_region_contrast_heatmap_data(...)`

Recommended signature shape:

```python
payload = plotting.prepare_region_contrast_profile_data(
    result=contrast_result,
    position_table=position_table,
    axis=axis_spec,
    aggregation=aggregation_spec,
    value_mode="all",
)
```

```python
payload = plotting.prepare_region_contrast_heatmap_data(
    result=contrast_result,
    position_table=position_table,
    axis=axis_spec,
    aggregation=aggregation_spec,
    value_mode="all",
)
```

Where:

- `result` is a `RegionContrastResult`
- `position_table` is an already-aggregated positional summary table
- `axis` is the shared `AxisSpec`
- `aggregation` is the shared `AggregationSpec`
- `value_mode` controls whether the payload exposes numerator, denominator, delta, or all three

## Positional Source Contract

The first implementation should consume already-aggregated per-region positional summaries, not per-read raw events.

The minimal required columns are:

- `region_id`
- `position`
- `value`
- `region_strand`

And at least one grouping key that can be related to the contrast result:

- `sample_id`, or
- `condition`

Optional columns that should be preserved when present:

- `chrom`
- `region_start`
- `region_end`
- `segment_id`
- `weight`
- representation-specific columns such as `cluster`

The plotting layer should not recompute low-level pileup or extract summaries. It only validates, aligns, aggregates, and annotates.

## Supported V1 Data Paths

### 1. Pileup-backed positional abundance plots

This is the primary first path.

The positional table represents aggregate region-aligned signal over known loci, for example:

- per-region modification fraction by position
- per-region enrichment/depth-like positional summaries

These are joined to `RegionContrastResult` for:

- numerator profile
- denominator profile
- delta profile
- per-region heatmaps of the same views

### 2. Cluster-occupancy follow-up

This is allowed when the caller provides a compatible positional source table.

The plotting helpers should not assume that occupancy contrasts are inherently positional. They only support them when the positional table carries the needed signal.

This keeps occupancy plotting opt-in and data-driven rather than pretending all occupancy summaries have x-axis structure.

## Axis Rules

These helpers must use the shared plotting-axis model already established in [2026-04-01-plotting-axis-architecture-design.md](/Users/ngamarra/Documents/GitHub/dimelo-toolkit/docs/superpowers/specs/2026-04-01-plotting-axis-architecture-design.md).

### V1 axis support

Required:

- `fixed_window`
- `orientation="genomic"`
- `orientation="region_5to3"`

Deferred but designed-for:

- aggregate `segment_map`

Single-read constraints do not apply here because these helpers are aggregate-only.

## Aggregation Rules

These helpers must use the shared `AggregationSpec`.

### V1 requirements

Required:

- `weighting="equal_region"`
- `layout="faceted"` as the default

Supported if straightforward from the current substrate:

- `equal_observation`
- `coverage_weighted`

The implementation should keep weighting explicit in metadata even if only a subset is active in the first code slice.

## Value Modes

The helpers should expose all three comparison views from one prepared payload:

- `numerator`
- `denominator`
- `delta`

Recommended behavior:

- `value_mode="all"` returns a long-form table with a `value_mode` column
- narrower modes may filter that table without changing the rest of the payload contract

This allows one prep call to feed:

- overlay plots
- side-by-side facets
- delta-only tracks
- numerator / denominator / delta heatmaps

## Output Contract

Each helper should return a renderer-neutral payload with:

- `plot_table`
- `axis_table`
- `metadata`
- `summary_table`

### `plot_table`

Long-form rows suitable for direct plotting. It should include:

- mapped `plot_x`
- `region_id`
- `value`
- `value_mode`
- grouping keys such as `condition` or `sample_id`

For heatmaps, the table may also include row-order metadata.

### `axis_table`

Shared axis metadata from the plotting core.

### `metadata`

Must include at least:

- `plot_family`
- `representation`
- `analysis_unit`
- `orientation`
- `coordinate_mode`
- `weighting`
- `layout`
- `value_modes`

### `summary_table`

Optional compact aggregate tables useful for default Matplotlib wrappers or quick inspection.

## Joining Rules

The plotting prep must join positional rows only for regions represented in the contrast result unless the caller explicitly asks otherwise later.

V1 behavior:

- filter `position_table` to contrast-supported `region_id`s
- join on the grouping keys available from the result and positional source
- fail clearly if the positional source cannot be matched to the contrast result

The helpers should prefer explicit failure over silently mixing incompatible region or condition scopes.

## Heatmap Semantics

Heatmaps should be aggregate region heatmaps, not single-read matrices.

V1 row model:

- one row per region for a chosen `value_mode`
- optional faceting by numerator / denominator / delta

Possible row ordering options later:

- by contrast rank
- by delta magnitude
- by supplied region order

The first implementation can default to contrast rank order.

## Compatibility Constraint

This design must preserve continuity with the existing package.

That means:

- legacy plotters remain supported
- new helpers are additive
- old flags such as `regions_5to3prime` continue to map onto the shared orientation model
- users can keep plotting directly from old parse/load outputs if they prefer

The new helpers are for users who want the newer analysis-layer workflows without losing flexibility over rendering.

## Failure Modes

The helpers should fail clearly on:

- missing required positional columns
- incompatible `region_id` or grouping keys
- unsupported axis modes for the current slice
- ambiguous joins between contrast rows and positional rows
- positional tables that are not already aggregated to the expected resolution

## Tests

The implementation plan should include tests for:

- successful profile prep from a minimal pileup-backed positional table
- successful heatmap prep from the same substrate
- `value_mode="all"` includes numerator, denominator, and delta
- `region_5to3` orientation uses the shared fixed-window axis rules
- rejection of positional tables missing required join keys
- rejection of incompatible contrast / positional scopes
- compatibility with `cluster_occupancy` only when a positional source table is explicitly supplied

## Recommended Build Order

1. Add shared internal helpers in `dimelo.plotting` for joining contrast results to positional tables.
2. Implement `prepare_region_contrast_profile_data(...)` for fixed-window aggregate inputs.
3. Implement `prepare_region_contrast_heatmap_data(...)` on the same substrate.
4. Add tests for pileup-backed known-region contrasts.
5. Document how these helpers relate to old plotting entry points and the newer renderer-neutral payload model.

