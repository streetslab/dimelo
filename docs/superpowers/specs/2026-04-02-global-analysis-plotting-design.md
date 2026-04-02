# Global Analysis Plotting Design

## Goal

Add a renderer-neutral plotting-prep layer for `dimelo.global_analysis` that makes both whole-sample summaries and tiled-window outputs easier to visualize without coupling the analysis API to a specific plotting backend.

The first slice should help users answer two questions:

- how do global motif fractions and normalization offsets vary across samples and conditions?
- how do tiled genome windows vary across contigs at broad scale?

The plotting layer should remain additive and data-first, using the existing `GlobalAnalysisResult` outputs as the source of truth.

## Scope

This design covers:

- plotting-prep helpers for whole-sample global summaries
- plotting-prep helpers for tiled-window broad genome views
- per-contig organization as the default window layout
- condition-level aggregation alongside per-sample summaries

This design does not cover:

- renderer functions
- cumulative genome axes across contigs
- normalized region-axis or segment-map plotting
- direct integration with discovery hit overlays
- any new normalization method beyond what `GlobalAnalysisResult` already provides

## User-Facing API

Add two public helpers in [plotting.py](../../../dimelo/plotting.py):

- `prepare_global_analysis_summary_data(...)`
- `prepare_global_analysis_window_data(...)`

These helpers should consume a [GlobalAnalysisResult](../../../dimelo/models.py) and return renderer-neutral payload dictionaries with plot-ready tables plus metadata.

## Design Principles

### Keep the analysis layer authoritative

The plotting helpers should expose and reshape the values already computed by [global_analysis.py](../../../dimelo/global_analysis.py). They should not invent new biological normalization logic or re-derive alternative summary quantities beyond simple aggregations.

### Separate summary and window concerns

`GlobalAnalysisResult` naturally splits into:

- global per-sample summary values and normalization factors
- tiled broad-window values across contigs

The plotting API should mirror that split rather than combining both families behind a single mode switch.

### Keep window views genome-native

Broad window plots should stay per-contig by default. The helper should not flatten contigs into one synthetic cumulative axis in the first slice.

## Plot Family 1: Global Summary

### Purpose

Prepare whole-sample and condition-level plotting payloads for global motif fractions and normalization offsets.

### Public Entry Point

```python
payload = plotting.prepare_global_analysis_summary_data(
    result=result,
    motifs=None,
    aggregate_conditions=True,
)
```

### Inputs

- `result: GlobalAnalysisResult`
- `motifs: list[str] | None`
  - optional motif filter
- `aggregate_conditions: bool`
  - when `True`, prepare condition-level summaries in addition to sample-level rows

### Outputs

Return a dictionary with:

- `sample_summary`
  - one row per sample and motif
  - includes:
    - `sample_id`
    - `condition`
    - `replicate`
    - `motif`
    - `global_fraction`
- `condition_summary`
  - one row per condition and motif
  - includes aggregated summaries such as:
    - `condition`
    - `motif`
    - `global_fraction_mean`
    - `global_fraction_median`
    - `sample_n`
  - this table may be empty when `aggregate_conditions=False`
- `normalization_table`
  - sample-level normalization factors
  - includes:
    - `sample_id`
    - `condition`
    - `replicate`
    - `motif`
    - `global_fraction`
    - `reference_fraction`
    - `global_offset`
- `metadata`
  - includes:
    - `motifs`
    - `aggregate_conditions`

### Aggregation rules

The helper should:

- preserve per-sample rows directly from `result.summary`
- preserve sample-level normalization rows directly from `result.normalization_factors`
- compute condition-level summaries by grouping per sample, not by inventing new weighted schemes

The first implementation should use simple unweighted condition-level summaries, with mean, median, and sample count as the primary outputs.

## Plot Family 2: Broad Window Views

### Purpose

Prepare tiled-window plotting payloads for broad genome inspection by sample or condition.

### Public Entry Point

```python
payload = plotting.prepare_global_analysis_window_data(
    result=result,
    contigs=None,
    motifs=None,
    aggregate_conditions=False,
)
```

### Inputs

- `result: GlobalAnalysisResult`
- `contigs: list[str] | None`
  - optional filter and ordering override
- `motifs: list[str] | None`
  - optional motif filter
- `aggregate_conditions: bool`
  - when `True`, provide condition-level window summaries in addition to sample-level windows

### Outputs

Return a dictionary with:

- `window_table`
  - sample-level tiled windows
  - includes:
    - `sample_id`
    - `condition`
    - `replicate`
    - `motif`
    - `contig`
    - `start`
    - `end`
    - `window_midpoint`
    - `window_fraction`
- `condition_window_table`
  - optional condition-level aggregation of windows
  - includes:
    - `condition`
    - `motif`
    - `contig`
    - `start`
    - `end`
    - `window_midpoint`
    - `window_fraction_mean`
    - `window_fraction_median`
    - `sample_n`
  - may be empty when `aggregate_conditions=False`
- `metadata`
  - includes:
    - `contig_order`
    - `motifs`
    - `aggregate_conditions`

### Default layout semantics

The helper should default to per-contig organization, not cumulative genome coordinates. Renderers are expected to facet or panel by contig by default.

### Derived columns

The helper should add:

- `contig`
  - normalized alias of the `chromosome` field
- `window_midpoint`
  - `(start + end) / 2`

The helper should not alter the underlying window values beyond optional condition-level aggregation.

## Defaults And Validation

### Result validation

Helpers should validate they received a `GlobalAnalysisResult`-like object with the required attributes:

- `summary`
- `windows`
- `normalization_factors`
- `metadata`

### Motif filtering

If `motifs` is provided:

- filter the relevant tables to those motifs
- raise a clear error if none of the requested motifs are present

### Contig filtering

If `contigs` is provided:

- filter the window table and condition-window table to those contigs
- preserve the given order in `metadata["contig_order"]`
- raise only if none of the requested contigs are present

### Empty results

If the underlying tables are empty:

- return empty stable tables with expected columns
- preserve metadata where possible
- do not raise unless required structure is missing

## Relationship To Existing `plot_data`

`GlobalAnalysisResult.plot_data` currently includes lightweight payloads:

- `global_fraction_bar`
- `window_fraction_table`

The new helpers should not depend on those names as their public contract. They may reuse the underlying tables internally if clean, but the canonical inputs remain:

- `result.summary`
- `result.windows`
- `result.normalization_factors`

## Backward Compatibility

This design is additive.

- it does not change `run_global_analysis()`
- it does not change `GlobalAnalysisResult`
- it does not remove existing `plot_data`
- it does not introduce renderer-specific coupling

## Testing

Add focused tests in [tests/test_plotting.py](../../../tests/test_plotting.py) for:

- stable summary payload columns
- condition-level aggregation behavior for summary data
- stable window payload columns
- per-contig ordering behavior
- motif and contig filtering
- empty-result behavior

Use synthetic `GlobalAnalysisResult` fixtures for most helper tests rather than full end-to-end global-analysis runs.

## Documentation

Update:

- [docs/global-analysis.md](../../../docs/global-analysis.md)
- [README.md](../../../README.md)

Document:

- the two new helpers
- the distinction between sample-level and condition-level summaries
- the per-contig default for window payloads
- that the helpers remain data-prep only and renderer-neutral

## Recommended Implementation Order

1. Add plotting helper tests for global summary prep
2. Implement `prepare_global_analysis_summary_data(...)`
3. Add plotting helper tests for broad window prep
4. Implement `prepare_global_analysis_window_data(...)`
5. Update docs

## Recommendation

Build `global_analysis` plotting as a simple, result-centric extension of the current data-first architecture:

- `prepare_global_analysis_summary_data(...)`
- `prepare_global_analysis_window_data(...)`

Keep the first slice additive, per-contig by default for window views, and explicitly expose both per-sample and condition-level summaries so users can move quickly between QC-style inspection and broader comparison views.
