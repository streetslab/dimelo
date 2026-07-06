# Plotting Axis Architecture Design

## Purpose

This design defines a compatibility-safe plotting architecture for `dimelo-toolkit` that:

- preserves the existing plotting entry points and familiar arguments
- introduces a shared internal model for positional plotting
- supports `5'->3'` orientation consistently across region-aligned plots
- supports both fixed-window and variable-length aggregate region normalization
- keeps single-read plots coordinate-preserving
- keeps plotting data-first so users can use Matplotlib, seaborn, Plotly, Altair, or their own stack

This is an additive architecture around the current package. It must not require a rewrite of the parsing layer or force existing users to adopt a new plotting API immediately.

## Compatibility Constraint

Backward continuity is a hard requirement.

The update must:

- keep existing plotters working
- preserve current public plotting signatures unless a change is strictly additive
- preserve familiar flags such as `relative`, `single_strand`, and `regions_5to3prime`
- preserve current behavior by default for legacy plotting entry points
- keep Matplotlib as the built-in default renderer

The new plotting system becomes the canonical internal model. Existing plotting functions remain supported as wrappers or adapters over that model.

## Problem Statement

The package currently has several plotting paths with overlapping ideas:

- region-relative plotting
- strand-aware plotting
- `regions_5to3prime` flipping
- single-read displays
- aggregate enrichment and depth profiles

These concepts already exist, but they are distributed across multiple modules and are not expressed as one shared plotting contract.

The new analysis layers also need more explicit plotting semantics:

- global summaries
- region discovery scan plots
- known-region contrasts
- clustering summaries
- occupancy follow-up

Users also need more flexible positional plotting for region-aligned analyses, especially:

- fixed centered windows
- variable-length region normalization for aggregate plots
- arbitrary segment maps for reference-defined feature structure
- support for both contiguous and non-contiguous feature segments

The design must make those capabilities available without making single-read plots scientifically misleading.

## Core Design Rule

There are three separate concerns that must stay independent:

1. `plot_family`
   What kind of visualization is being produced.

2. `axis_spec`
   How x coordinates are defined and aligned.

3. `aggregation_spec`
   How values are summarized across reads, regions, or conditions.

The term “normalization” must not be used ambiguously across these concerns.

## Safe Semantic Boundary

Single-read plots and aggregate plots must not share the same x-axis freedoms.

### Single-read plots

Single-read plots must remain coordinate-preserving.

Allowed:

- genomic coordinates
- fixed-window relative coordinates
- `5'->3'` orientation flips for strand-consistent viewing
- segment-based faceting or panel splits

Not allowed:

- stretched or continuously scaled x-axes across variable-length regions
- metaregion-style body rescaling on a continuous single-read canvas

This is required to avoid misleading geometric interpretations of read structure.

### Aggregate plots

Aggregate positional plots may use normalized axes because they explicitly represent a summarized metaregion view rather than raw per-read geometry.

Allowed:

- genomic or region-oriented axes
- fixed-window aggregation
- scaled contiguous region bodies
- arbitrary segment maps
- contiguous and non-contiguous segment composition
- concatenated or faceted segment layouts

## Canonical Axis Model

The canonical internal model should be a shared positional axis system used by plot-data preparation.

### Axis families

There are three user-facing axis families:

1. `fixed_window`
   One continuous anchored window around a position, center, or landmark.

2. `scaled_region`
   One contiguous variable-length region body, optionally with fixed flanks.

3. `segment_map`
   One plotted axis assembled from explicit per-region reference segments.

Internally:

- `fixed_window` compiles to a one-segment representation
- `scaled_region` compiles to a structured segment representation
- `segment_map` is the canonical general form

### Orientation

Orientation is independent of coordinate family.

Supported values:

- `genomic`
- `region_5to3`

`region_5to3` means region-aligned plots are oriented along the reference feature direction. For negative-strand regions, this implies flipping positional coordinates in the derived plot view.

This is conceptually the generalization of today’s `regions_5to3prime` behavior.

## Segment Model

The plotting system must support arbitrary segments keyed directly in reference coordinates.

The base contract is explicit per-region segment definitions. Feature-derived helpers are convenience layers, not the core contract.

### Segment fields

Each segment should support:

- `segment_id`
- `label`
- `region_id`
- `start_ref`
- `end_ref`
- `mode`
  - `raw`
  - `scaled`
- `bins` or `bin_width`
- `contiguous_with_previous`
- `plot_gap_after`

### Contiguous vs non-contiguous segments

Non-contiguous segments are not a separate analysis mode. They are represented by segment metadata.

This allows one common substrate to support:

- `upstream | body | downstream`
- exon-only plots
- motif neighborhoods
- arbitrary user-specified feature composites

Renderers must be able to express discontinuity either through:

- separators or gaps on one concatenated axis
- faceted small multiples by segment

## Aggregation Model

Aggregate positional plotting must make weighting explicit.

This is separate from x-axis normalization.

### Weighting options

The design should support:

- `equal_region`
  Each region contributes equally after within-region summarization.

- `equal_observation`
  Each observation contributes equally.

- `coverage_weighted`
  Regions are weighted by coverage or another confidence-like quantity.

The recommended default is `equal_region`.

### Within-region summary

The design should support at least:

- `mean`
- `fraction`
- `density`

This defines how a region is summarized internally before any cross-region aggregation.

### Signal normalization

Signal normalization must be separate from x-axis normalization and weighting.

The plotting contract should support metadata for at least:

- `none`
- `per_region`
- `global`
- `control_regions`

Whether all of these are implemented in the first slice is a separate scope decision. The design requires the conceptual separation now.

### Order-of-operations clarity

The plotting system must keep the following choices explicit:

- normalize each region then average
- average raw signal then normalize once

These are not interchangeable and should not be hidden behind one plotting flag.

## Plot Families

The plotting system should define explicit plot families rather than one generic plotting function.

### 1. `single_read_raster`

For read-level displays.

Typical uses:

- read rasters
- read event/lollipop plots
- browser-like region views

Constraints:

- coordinate-preserving only
- may use `genomic` or `region_5to3`
- may use `fixed_window`
- may use segmented faceting
- must not use continuous scaled metaregion axes

### 2. `aggregate_profile`

For aggregated enrichment, depth, abundance, or occupancy-like profiles over aligned regions.

Typical uses:

- promoter-centered metaprofiles
- gene-body plots
- upstream/body/downstream summaries
- segment-based aggregate traces

Capabilities:

- `fixed_window`
- `segment_map`
- `region_5to3`
- aggregate-only scaled segments
- concatenated or faceted segment layout

### 3. `aggregate_heatmap`

For region-by-position summaries.

Capabilities:

- same axis semantics as `aggregate_profile`
- rows can represent regions, conditions, clusters, or ranked loci
- supports faceted or concatenated segment layouts

### 4. `distribution_summary`

For non-positional summary plots.

Typical uses:

- cluster distribution bars
- condition heatmaps
- volcano plots
- rank plots
- occupancy summary charts

No positional axis normalization applies.

### 5. `track_scan`

For genome scan or local discovery tracks.

Typical uses:

- region discovery scan tracks
- hit-context plots
- chromosome-level summaries

These are primarily genomic-coordinate plots.

## Analysis-to-Plot Mapping

### `global_analysis`

Primary plot families:

- `distribution_summary`
- `track_scan`
- `aggregate_profile` where region-aligned summaries exist

### `region_discovery`

Primary plot families:

- `track_scan`
- `distribution_summary`
- local `aggregate_profile` views for top hits

### `region_contrasts`

Primary plot families:

- `distribution_summary`
- `aggregate_profile`
- `aggregate_heatmap`

### `shared_clustering`

Primary plot families:

- `distribution_summary`
- positional `aggregate_profile` later where region-aligned occupancy summaries are meaningful

### Legacy read/profile plotting

Legacy plotting modules should gradually route through the same shared positional data-prep semantics where feasible.

## API Model

The advanced plotting API should live primarily in a shared plotting module, with thin wrappers elsewhere.

Recommended shared objects:

```python
@dataclass
class AxisSpec:
    orientation: str  # "genomic" | "region_5to3"
    coordinate_mode: str  # "fixed_window" | "segment_map"
    anchor: str | None = None
    upstream_bp: int | None = None
    downstream_bp: int | None = None
    segments: list["SegmentSpec"] | None = None


@dataclass
class SegmentSpec:
    segment_id: str
    label: str
    start_ref: int
    end_ref: int
    mode: str  # "raw" | "scaled"
    bins: int | None = None
    plot_gap_after: bool = False
    contiguous_with_previous: bool = True


@dataclass
class AggregationSpec:
    weighting: str  # "equal_region" | "equal_observation" | "coverage_weighted"
    within_region_summary: str  # "mean" | "fraction" | "density"
    signal_normalization: str  # "none" | "per_region" | "global" | "control_regions"
    layout: str  # "concatenated" | "faceted"


@dataclass
class PlotRequest:
    plot_family: str
    axis: AxisSpec | None = None
    aggregation: AggregationSpec | None = None
    renderer: str = "matplotlib"
```

Recommended shared prep entry points:

- `prepare_single_read_plot_data(...)`
- `prepare_aggregate_plot_data(...)`

These prepare data and metadata only. Rendering is a separate thin layer.

## Data Contract

Every plot-prep function should return structured plotting payloads rather than figures as the canonical contract.

Recommended outputs:

- `plot_table`
  Long-form plotting data.

- `axis_table`
  Segment boundaries, axis labels, gap markers, and orientation metadata.

- `metadata`
  Plot semantics, units, normalization mode, weighting mode, and family identifiers.

- `summary_table`
  Optional pre-aggregated values when useful.

This keeps the plotting system backend-neutral and aligned with the existing data-first direction of the newer workflows.

## Legacy Wrapper Strategy

The new plotting architecture should not replace current plotters at once.

Instead:

- legacy functions remain public
- legacy flags are translated into the new shared specs internally
- defaults remain behaviorally close to the current package
- advanced users can opt into the richer plotting API directly

Examples of translation:

- `regions_5to3prime=True`
  -> `AxisSpec(orientation="region_5to3", ...)`

- `relative=True` with centered windows
  -> `AxisSpec(coordinate_mode="fixed_window", ...)`

- older profile plotting
  -> `plot_family="aggregate_profile"`

### User-facing guidance

The docs should explicitly map old concepts to new ones:

- old `regions_5to3prime`
  -> new `orientation="region_5to3"`

- old relative region profile plotting
  -> new `fixed_window`

- new scaled or segmented metaregion plotting
  -> aggregate-only advanced capability

- single-read plots remain geometry-preserving

## Renderer Strategy

Built-in rendering should remain Matplotlib-first for compatibility.

The package should still encourage users to plot returned tables with their own preferred stack.

So the plotting stack becomes:

1. canonical plot-data preparation
2. thin built-in Matplotlib rendering
3. optional user-managed rendering elsewhere

This is consistent with the current direction of `plot_data` in newer workflow results.

## Module Layout

Recommended direction:

- expand [plotting.py](../../../dimelo/plotting.py)
  to hold the shared plot-data preparation core

- keep existing modules such as:
  - [plot_enrichment_profile.py](../../../dimelo/plot_enrichment_profile.py)
  - [plot_depth_profile.py](../../../dimelo/plot_depth_profile.py)
  - [plot_reads.py](../../../dimelo/plot_reads.py)
  as wrappers or adapters during the transition

- route newer workflow `plot_data` preparation through the shared plotting core where it makes sense

## Failure and Validation Rules

The shared plotting layer should validate semantic misuse early.

Required validations:

- reject continuous scaled metaregion axes for single-read plots
- reject positional axis specs on non-positional plot families where not applicable
- reject malformed segment maps
- require explicit segment discontinuity metadata for non-contiguous segment layouts
- preserve legacy defaults when wrappers are used without advanced options

## Testing Requirements

The implementation plan for this design must include:

- regression coverage proving existing plotters still work with legacy arguments
- unit tests for shared `AxisSpec` and segment validation
- unit tests for `region_5to3` flipping in plot-data prep
- unit tests proving single-read plotting rejects scaled continuous metaregion axes
- unit tests for aggregate segment-map prep with:
  - contiguous segments
  - non-contiguous segments
  - concatenated layout metadata
  - faceted layout metadata
- tests for weighting modes in aggregate profile prep

## Implementation Order

Recommended implementation sequence:

1. Add shared plotting spec dataclasses and validation in `plotting.py`.
2. Add shared positional plot-data prep for `fixed_window`.
3. Add aggregate-only segment-map prep with contiguous and non-contiguous support.
4. Add `region_5to3` orientation handling in the shared prep layer.
5. Route one aggregate legacy plotter through the shared core.
6. Route one single-read legacy plotter through the shared core.
7. Extend newer workflow `plot_data` helpers where useful.
8. Add compatibility docs mapping old plotting flags to the new model.

## Recommendation

The package should standardize on a shared plotting axis model while keeping plotting behavior safe and scientifically honest:

- orientation is general
- scaled axes are aggregate-only
- single-read plots preserve geometry
- legacy plotting calls remain supported
- plotting data stays the canonical contract

This gives users more expressive aggregate plotting for variable-length features and arbitrary reference-defined segments without breaking the mental model of existing plotting routines.
