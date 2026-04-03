# Shared Clustering Plotting Design

## Goal

Add a richer, renderer-neutral plotting-prep layer for `dimelo.workflows.shared_cluster_distribution()` that standardizes the main `SharedClusterResult` plotting families:

- distribution and change views
- cluster profile views
- region-level occupancy views

This should deepen the current lightweight `result.plot_data` contract without replacing it.

## Scope

This design covers:

- plotting-prep helpers for `SharedClusterResult`
- distribution and change payloads from:
  - `cluster_distribution`
  - `condition_distribution`
  - `distribution_change`
- cluster profile payloads from `cluster_profiles`
- region-level occupancy payloads from `region_summaries`
- optional condition-level aggregation for region occupancy

This design does not cover:

- renderer implementations
- new clustering logic
- new region statistics
- new occupancy models
- normalized region-axis plotting

## User-Facing API

Add three public helpers in [plotting.py](../../../dimelo/plotting.py):

- `prepare_shared_cluster_distribution_data(...)`
- `prepare_shared_cluster_profile_data(...)`
- `prepare_shared_cluster_region_data(...)`

All three should accept a [SharedClusterResult](../../../dimelo/models.py) and return renderer-neutral payload dictionaries with plot-ready tables plus metadata.

## Design Principles

### Keep `SharedClusterResult` authoritative

The helpers should use the existing `SharedClusterResult` tables as the source of truth. They should not derive new biological interpretations or replace the workflow’s current `plot_data`.

### Keep the helpers family-specific

The plotting families are different enough that one generic mode-switched helper would become opaque. The API should follow the same pattern already used for `region_contrasts`, `region_discovery`, and `global_analysis`: multiple small public helpers with shared internals.

### Preserve the existing lightweight `plot_data`

The current `SharedClusterResult.plot_data` should remain available. The new helper layer is additive and richer, not a replacement.

## Plot Family 1: Distribution And Change Views

### Purpose

Prepare plot-ready tables for sample-level distributions, condition-level distributions, and condition-change summaries.

### Public Entry Point

```python
payload = plotting.prepare_shared_cluster_distribution_data(result=result)
```

### Outputs

Return a dictionary with:

- `sample_distribution`
  - rows from `result.cluster_distribution`
- `condition_distribution`
  - rows from `result.condition_distribution`
- `distribution_change`
  - rows from `result.distribution_change`
  - may be empty if no reference condition/change table exists
- `metadata`
  - includes:
    - `mode`
    - `motifs`
    - `reference_condition`

### Intended use

This should support:

- cluster-fraction bars by sample
- condition-level heatmaps
- change/delta heatmaps or scatter summaries

## Plot Family 2: Cluster Profile Views

### Purpose

Prepare plot-ready tables for visualizing the feature profile of each cluster.

### Public Entry Point

```python
payload = plotting.prepare_shared_cluster_profile_data(result=result)
```

### Outputs

Return a dictionary with:

- `profile_table`
  - rows from `result.cluster_profiles`
  - standardized so cluster identifier and feature columns are easy to consume
- `metadata`
  - includes:
    - `mode`
    - `motifs`
    - `feature_names`

### Intended use

This should support:

- per-cluster profile panels
- line/scatter plots over feature dimensions
- summary comparisons of cluster centroids or average patterns

The helper should not assume any one renderer or feature geometry.

## Plot Family 3: Region Occupancy Views

### Purpose

Prepare plot-ready tables for region-level occupancy inspection from `result.region_summaries`.

### Public Entry Point

```python
payload = plotting.prepare_shared_cluster_region_data(
    result=result,
    aggregate_conditions=True,
)
```

### Inputs

- `result: SharedClusterResult`
- `aggregate_conditions: bool`
  - when `True`, provide condition-level region summaries in addition to sample-level rows

### Outputs

Return a dictionary with:

- `region_table`
  - sample-level rows from `result.region_summaries`
- `condition_region_table`
  - optional condition-level aggregation over `region_table`
  - may be empty when `aggregate_conditions=False`
- `metadata`
  - includes:
    - `mode`
    - `motifs`
    - `aggregate_conditions`

### Aggregation rules

This helper should only aggregate what `result.region_summaries` already contains. It should not invent a new occupancy model.

When condition aggregation is enabled, aggregate from sample-level region rows first so the summary semantics match the reported counts.

## Defaults And Validation

### Result validation

Helpers should validate they received a `SharedClusterResult`-like object with the required attributes for the requested family.

Examples:

- distribution helper requires:
  - `cluster_distribution`
  - `condition_distribution`
  - `metadata`
- profile helper requires:
  - `cluster_profiles`
  - `model`
- region helper requires:
  - `region_summaries`
  - `metadata`

### Empty tables

If a relevant table is empty or absent:

- return stable empty tables with expected columns
- preserve metadata where possible
- do not raise unless required structure is missing for the helper contract

### Region helper mode handling

If `result.region_summaries` is `None`, the region helper should return stable empty outputs rather than failing, because `read_global` clustering does not always produce region summaries.

## Relationship To Existing `plot_data`

`SharedClusterResult.plot_data` currently includes lightweight payloads such as:

- `cluster_distribution_bar`
- `cluster_distribution_heatmap`

Those should remain in place.

The new helpers should:

- preserve compatibility with the existing workflow outputs
- provide a more structured, result-centric plotting-prep interface
- avoid depending on the existing lightweight `plot_data` keys as their public contract

## Backward Compatibility

This design is additive.

- it does not change `shared_cluster_distribution()`
- it does not change `SharedClusterResult`
- it does not remove existing `result.plot_data`

## Testing

Add focused tests in [tests/test_plotting.py](../../../tests/test_plotting.py) for:

- stable distribution payload keys and table propagation
- empty change-table behavior
- profile-table passthrough and metadata
- region helper empty behavior when `region_summaries` is `None`
- condition-level region aggregation behavior

Use synthetic `SharedClusterResult` fixtures rather than end-to-end clustering runs for most helper tests.

## Documentation

Update:

- [docs/shared-clustering.md](../../../docs/shared-clustering.md)
- [README.md](../../../README.md)

Document:

- the three new helpers
- the distinction between the old lightweight `result.plot_data` and the richer helper layer
- the fact that `shared_clustering` now has dedicated renderer-neutral helper prep, not just raw workflow payloads

## Recommended Implementation Order

1. Add plotting helper tests for distribution/change prep
2. Implement `prepare_shared_cluster_distribution_data(...)`
3. Add plotting helper tests for profile prep
4. Implement `prepare_shared_cluster_profile_data(...)`
5. Add plotting helper tests for region prep
6. Implement `prepare_shared_cluster_region_data(...)`
7. Update docs

## Recommendation

Build `shared_clustering` plotting as a result-centric, data-first layer with three helper families:

- `prepare_shared_cluster_distribution_data(...)`
- `prepare_shared_cluster_profile_data(...)`
- `prepare_shared_cluster_region_data(...)`

That brings `shared_clustering` into line with the richer plotting-prep coverage now available for `region_contrasts`, `region_discovery`, and `global_analysis`, while preserving the existing lightweight `result.plot_data` contract.
