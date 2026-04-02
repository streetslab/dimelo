# Region Discovery Plotting Design

## Goal

Add a small, data-first plotting-prep layer for `dimelo.region_discovery` that makes discovery results easier to inspect without coupling the analysis API to a specific renderer. The first slice should help users answer two questions:

- where are the interesting windows across the genome scan?
- what does the local scored neighborhood around the top hits look like?

The new helpers should sit on top of the existing `RegionDiscoveryResult` contract rather than expanding `scan_genome()` itself.

## Scope

This design covers:

- renderer-neutral plotting prep for discovery scan overviews
- renderer-neutral plotting prep for local top-hit context views
- per-contig organization as the default scan layout
- continued use of `RegionDiscoveryResult.windows` and `RegionDiscoveryResult.hits` as the source of truth

This design does not cover:

- Matplotlib or seaborn renderer functions
- cumulative genome axes across contigs
- normalized region-axis plotting
- direct single-read plotting
- segment-map or metaregion plotting semantics

## User-Facing API

Add two public helpers in [plotting.py](../../../dimelo/plotting.py):

- `prepare_region_discovery_scan_data(...)`
- `prepare_region_discovery_hit_context_data(...)`

These should accept a [RegionDiscoveryResult](../../../dimelo/models.py) and return renderer-neutral payload dictionaries with plot-ready tables plus metadata.

## Design Principles

### Keep discovery plotting genome-native

Discovery is still a tiled genome scan. The first plotting helpers should preserve genomic structure rather than forcing the result into a normalized region axis. The primary organization is:

- per-contig scan windows
- ranked hit overlays
- local window neighborhoods around selected hits

### Keep `scan_genome()` focused on analysis

`scan_genome()` already returns `result.windows`, `result.hits`, and lightweight `plot_data`. The richer plotting helpers should live in [plotting.py](../../../dimelo/plotting.py), not inside [region_discovery.py](../../../dimelo/region_discovery.py). That keeps analysis and plotting concerns separate and matches the newer `region_contrasts` plotting architecture.

### Prefer result-centric helpers

The primary public input should be `RegionDiscoveryResult`, not raw tables. Internally, the helpers can still work against `result.windows` and `result.hits`, but the public API should stay aligned with the analysis workflow output.

## Plot Family 1: Scan Overview

### Purpose

Prepare data for genome-scan views that show the scored window landscape and optionally overlay ranked discovery hits.

### Public Entry Point

```python
payload = plotting.prepare_region_discovery_scan_data(
    result=result,
    contigs=None,
    top_n_hits=100,
    score_column=None,
    include_all_windows=True,
)
```

### Inputs

- `result: RegionDiscoveryResult`
- `contigs: list[str] | None`
  - optional filter and ordering override for contigs to include
- `top_n_hits: int | None`
  - optional cap on hit overlays used for display
- `score_column: str | None`
  - defaults to the canonical discovery score field for the result shape
- `include_all_windows: bool`
  - when `True`, return the full filtered `scan_table`
  - when `False`, allow callers to request a lighter payload later

### Outputs

Return a dictionary with:

- `scan_table`
  - one row per scanned window
  - minimally includes:
    - `contig`
    - `start`
    - `end`
    - `score_value` or the chosen score field
    - any rank/effect columns already present in `result.windows`
    - `is_hit`
- `hit_table`
  - one row per displayed hit
  - top-N filtered if requested
- `metadata`
  - includes:
    - `contig_order`
    - `score_column`
    - `score_mode`
    - `contrast_mode`
    - `merge_hits`
    - `top_n_hits`

### Default layout semantics

The helper should default to per-contig organization, not a cumulative genome coordinate. Renderers are expected to facet or panel by contig by default.

### Derived columns

The helper should add:

- `contig`
  - normalized alias of the discovery `chromosome` field for plotting clarity
- `window_midpoint`
  - useful for line/scatter/track renderers
- `is_hit`
  - `True` for windows that correspond to the overlay hit set

The helper should not mutate or reinterpret the discovery score itself.

## Plot Family 2: Hit Context

### Purpose

Prepare local drill-down tables around selected hits so users can inspect the scored neighborhood before exporting BEDs or moving into clustering and contrasts.

### Public Entry Point

```python
payload = plotting.prepare_region_discovery_hit_context_data(
    result=result,
    top_n=12,
    hit_ids=None,
    padding_windows=5,
    padding_bp=None,
)
```

### Inputs

- `result: RegionDiscoveryResult`
- `top_n: int | None`
  - number of ranked hits to include by default
- `hit_ids: list[str] | None`
  - explicit hit/window selection override
- `padding_windows: int | None`
  - include this many neighboring windows on each side of each selected hit
- `padding_bp: int | None`
  - optional physical-distance expansion alternative

The first implementation should support `top_n` and `padding_windows` first. `hit_ids` and `padding_bp` can be part of the public contract now if the implementation still keeps the first slice modest.

### Outputs

Return a dictionary with:

- `context_table`
  - local windows surrounding each selected hit
  - includes:
    - `contig`
    - `start`
    - `end`
    - `window_midpoint`
    - score columns
    - `selected_hit_id`
    - `selected_hit_rank`
    - `relative_window_offset`
    - `is_selected_hit`
- `selected_hits`
  - the anchor hits used for the context extraction
- `metadata`
  - includes:
    - `selection_mode`
    - `top_n`
    - `padding_windows`
    - `padding_bp`
    - `score_column`

### Layout semantics

This payload is designed for small multiples or per-hit panels. Each selected hit becomes a local plotting context, usually faceted by `selected_hit_id` or `selected_hit_rank`.

## Defaults And Validation

### Score column selection

If `score_column` is omitted, use:

- `score_value` when present
- otherwise the canonical primary score column for the result shape

If no usable score column exists, raise a clear error instead of guessing silently.

### Empty results

If `result.windows` or `result.hits` is empty:

- return empty plot tables with stable columns
- preserve metadata
- do not raise unless required input structure is missing

### Contig filtering

If `contigs` is provided:

- filter both windows and hits to the requested contigs
- preserve the given order in `metadata["contig_order"]`
- error only if none of the requested contigs are present

### Hit selection

For hit-context prep:

- `top_n` should use the ranked order already present in `result.hits`
- if both `top_n` and explicit IDs are provided later, explicit IDs should win
- if no hits are available, return an empty stable payload

## Relationship To Existing `plot_data`

The helpers should treat the discovery result tables as the canonical source of truth:

- `result.windows`
- `result.hits`

Existing lightweight `result.plot_data["window_score_table"]` and `result.plot_data["top_hits_table"]` can remain available, but they should not block the helper design. The new helpers may reuse them internally if that is clean, but the public contract should not depend on those specific precomputed payload names.

## Backward Compatibility

This design is additive.

- it does not change `scan_genome()` signatures
- it does not change `RegionDiscoveryResult`
- it does not remove existing `plot_data`
- it does not change older plotting modules

This keeps the discovery workflow stable while making downstream plotting more coherent.

## Testing

Add focused tests in [tests/test_plotting.py](../../../tests/test_plotting.py) for:

- stable scan payload columns
- per-contig ordering behavior
- top-N hit overlay filtering
- local hit-context extraction using neighboring windows
- empty-result behavior
- explicit score-column validation

Use mocked or synthetic `RegionDiscoveryResult` objects rather than full end-to-end discovery runs for most helper tests.

## Documentation

Update:

- [docs/region-discovery.md](../../../docs/region-discovery.md)
- [README.md](../../../README.md)

Document:

- the two new helpers
- the per-contig default
- the difference between scan overview prep and hit-context prep
- that discovery plotting remains data-first and renderer-neutral

## Recommended Implementation Order

1. Add plotting helper tests for scan overview prep
2. Implement `prepare_region_discovery_scan_data(...)`
3. Add plotting helper tests for hit-context prep
4. Implement `prepare_region_discovery_hit_context_data(...)`
5. Update docs

## Recommendation

Build `region_discovery` plotting as a result-centric, data-prep-only layer with two helpers:

- `prepare_region_discovery_scan_data(...)`
- `prepare_region_discovery_hit_context_data(...)`

Keep the first slice genome-native, per-contig by default, and centered on discovery windows and hits rather than normalized positional region views.
