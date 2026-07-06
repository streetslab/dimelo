# Shared Clustering Matplotlib Renderers Design

## Summary

Add a `shared_clustering` Matplotlib renderer slice on top of the existing
renderer-neutral payloads in `dimelo.plotting`.

The goal is to close the built-in renderer gap that remains after the first
Matplotlib slice landed for:

- `region_contrasts`
- `region_discovery`
- `global_analysis`

This slice should keep the same architecture:

- `dimelo.plotting` remains the source of truth for payload preparation
- `dimelo.plotting_matplotlib` remains a thin figure-construction layer
- only small prep-layer adjustments are allowed when current payloads are
  awkward to render cleanly

## Goals

- Add built-in Matplotlib renderers for the three main `shared_clustering`
  plotting families:
  - distribution/change
  - cluster profiles
  - region occupancy
- Keep the renderer layer thin and payload-driven.
- Preserve compatibility with the existing lightweight `result.plot_data`
  tables.
- Make condition-level occupancy the default visual story while preserving
  sample-level detail in the payload and renderer options.

## Non-Goals

- Add new clustering outputs or new clustering algorithms.
- Add new occupancy metrics or new feature extraction logic.
- Rework `SharedClusterResult` semantics.
- Replace the renderer-neutral prep helpers.
- Add non-Matplotlib backends in this slice.

## Why This Slice Now

`shared_clustering` now has the renderer-neutral prep layer that the other
analysis families already use:

- `prepare_shared_cluster_distribution_data(...)`
- `prepare_shared_cluster_profile_data(...)`
- `prepare_shared_cluster_region_data(...)`

What it does not yet have is the same built-in plotting path that now exists
for the other major result objects. That makes it the most obvious plotting
consistency gap in the repo.

## Public API

Add these functions to `dimelo.plotting_matplotlib`:

```python
plot_shared_cluster_distribution_matplotlib(payload, *, level="condition", ax=None, title=None)
plot_shared_cluster_change_matplotlib(payload, *, ax=None, title=None)
plot_shared_cluster_profile_heatmap_matplotlib(payload, *, ax=None, title=None)
plot_shared_cluster_profile_series_matplotlib(payload, *, ax=None, title=None)
plot_shared_cluster_region_matplotlib(payload, *, level="condition", ax=None, title=None)
```

Public contract:

- renderer takes one prepared payload
- optional Matplotlib axis can be injected
- renderer returns `(fig, ax)`
- renderers do not mutate payloads
- renderers do not recompute clustering results

## Renderer Families

### Distribution Renderer

Input:

- payload from `prepare_shared_cluster_distribution_data(...)`

Responsibilities:

- render sample-level or condition-level cluster fractions
- default to `level="condition"`
- use condition/sample labels from the prepared payload

This renderer should expose the current cluster-distribution story in a built-in
bar-like form without changing the payload contract.

### Change Renderer

Input:

- payload from `prepare_shared_cluster_distribution_data(...)`

Responsibilities:

- render `distribution_change` as a cluster-by-condition-change heatmap or
  equivalent matrix-style view
- preserve stable cluster ordering from the payload

This should give users an immediate built-in view of condition-level occupancy
shifts.

### Profile Heatmap Renderer

Input:

- payload from `prepare_shared_cluster_profile_data(...)`

Responsibilities:

- render cluster-by-feature summaries as a heatmap
- act as the default first-class profile view
- preserve cluster and feature ordering from the payload

This is the best default because profile tables are naturally matrix-shaped and
scale better than many overlaid series.

### Profile Series Renderer

Input:

- payload from `prepare_shared_cluster_profile_data(...)`

Responsibilities:

- render a simple per-cluster profile series across features
- serve as the lighter alternate view for small feature sets

This renderer should stay simple and not attempt to solve every profile layout.

### Region Occupancy Renderer

Input:

- payload from `prepare_shared_cluster_region_data(...)`

Responsibilities:

- default to `level="condition"`
- render condition-aggregated region occupancy as a heatmap
- support `level="sample"` to render the sample-level region table instead

This preserves the recommended interpretation split:

- default figure tells the condition-level biological story
- sample-level occupancy remains available for replicate inspection and QC

## Allowed Prep-Layer Adjustments

Small prep changes are allowed only if rendering exposes a real payload friction.
Examples of acceptable changes:

- explicit metadata for default sort or label order
- stable ordering columns that are already implicit in the payload
- a small long-to-wide convenience field when it prevents renderer-side data
  reshaping from becoming too opaque

Not acceptable:

- new clustering outputs
- new summary statistics
- new biological interpretation layers

## Testing Strategy

Tests should verify:

- each renderer returns Matplotlib figure/axis objects
- each renderer consumes the expected payload family
- default level selection behaves as intended
- simple sample/condition switching works where supported
- empty or minimal valid payloads are handled where the prep layer allows them

Tests should not try to enforce pixel-perfect output.

## Documentation Expectations

User-facing docs should describe these renderers as:

- optional built-in Matplotlib views on top of prepared clustering payloads
- consistent with the earlier Matplotlib renderer slice
- complementary to the older lightweight `result.plot_data` tables, not a
  replacement for them

## Recommendation

Implement this as one bounded slice:

1. add failing renderer tests
2. allow only narrowly justified prep tweaks if the payloads are awkward
3. add the five renderers above
4. document the new built-in views in `docs/shared-clustering.md` and
   `README.md`

This keeps the plotting architecture coherent while filling the last major
built-in renderer gap among the current high-level analysis families.
