# Matplotlib Renderers Design

## Summary

Add a first built-in Matplotlib renderer layer on top of the existing renderer-neutral
plot-preparation payloads.

The goal is to give users a familiar, backward-compatible plotting path without
breaking the data-first architecture that now underpins `region_contrasts`,
`region_discovery`, and `global_analysis`.

This renderer slice should:

- keep `dimelo.plotting` as the renderer-neutral prep layer
- add a separate Matplotlib-specific module for figure construction
- support optional file export helpers
- avoid moving biological aggregation logic back into renderer functions

## Goals

- Add thin Matplotlib renderers for the most mature plotting payloads.
- Preserve the data-first contract: prep helpers remain the source of truth.
- Provide a simple built-in plotting path for users who prefer Matplotlib.
- Support convenient figure export to image or PDF files.
- Keep the renderer implementation small, explicit, and easy to test.

## Non-Goals

- Build a backend plugin/registry system in this slice.
- Add Plotly, Altair, or other non-Matplotlib renderers in this slice.
- Rework the renderer-neutral payload schemas.
- Add new biological summaries or plotting prep behavior just to support rendering.
- Cover every analysis family at once.

## Why A Separate Matplotlib Module

The repo now has a clear prep layer in `dimelo.plotting`:

- prep helpers validate inputs
- prep helpers align coordinates and aggregation semantics
- prep helpers return payload dictionaries/tables

If Matplotlib rendering is added into the same module, that clean separation will erode
quickly. A separate module keeps responsibilities explicit:

- `dimelo.plotting`: data prep and renderer-neutral payloads
- `dimelo.plotting_matplotlib`: figure construction and file export

This also keeps future renderer additions straightforward.

## First-Version Scope

The first renderer slice should cover:

- `region_contrasts`
  - profile renderer
  - heatmap renderer
- `region_discovery`
  - scan renderer
  - hit-context renderer
- `global_analysis`
  - summary renderer
  - window renderer

It should not yet cover:

- `shared_clustering`
- full legacy wrapper migration
- new single-read Matplotlib renderers beyond what already exists in legacy plotting modules

## Public API

Add a new module:

- `dimelo.plotting_matplotlib`

It should expose:

```python
plot_region_contrast_profile_matplotlib(payload, *, value_mode="delta", ax=None, title=None)
plot_region_contrast_heatmap_matplotlib(payload, *, value_mode="delta", ax=None, title=None)
plot_region_discovery_scan_matplotlib(payload, *, axes=None, title=None)
plot_region_discovery_hit_context_matplotlib(payload, *, axes=None, title=None)
plot_global_analysis_summary_matplotlib(payload, *, level="condition", ax=None, title=None)
plot_global_analysis_window_matplotlib(payload, *, axes=None, title=None)
save_figure(fig, path, *, dpi=300, bbox_inches="tight", transparent=False)
```

Exact argument names may vary slightly, but the public contract should stay:

- renderer takes one prepared payload
- optional Matplotlib axes can be injected
- renderer returns figure/axes
- export helper takes a figure and destination path

## Renderer Contracts

### General Rules

Each renderer should:

- accept a prepared payload dict from `dimelo.plotting`
- validate only the keys/shapes it needs
- avoid mutating the payload
- avoid recomputing biological summaries already expressed in payload tables
- return `(fig, ax)` or `(fig, axes)` depending on the layout

Each renderer should not:

- call `prepare_*` internally from result objects
- infer new biological quantities
- write files automatically

### Region Contrast Profile Renderer

Input:

- payload from `prepare_region_contrast_profile_data(...)`

Responsibilities:

- render numerator / denominator / delta views from prepared profile rows
- respect axis labels and metadata from the payload
- support a `value_mode` switch:
  - `numerator`
  - `denominator`
  - `delta`

The renderer should assume the payload already contains the aligned plot table.

### Region Contrast Heatmap Renderer

Input:

- payload from `prepare_region_contrast_heatmap_data(...)`

Responsibilities:

- render a per-region heatmap for one `value_mode`
- preserve region rank order from the payload
- use payload axis metadata for x labels or range annotations where appropriate

### Region Discovery Scan Renderer

Input:

- payload from `prepare_region_discovery_scan_data(...)`

Responsibilities:

- render per-contig score tracks by default
- overlay or highlight hits when `hit_table` is non-empty
- keep contigs separated rather than forcing a cumulative pseudo-genome axis

### Region Discovery Hit-Context Renderer

Input:

- payload from `prepare_region_discovery_hit_context_data(...)`

Responsibilities:

- render local scored neighborhoods around selected hits
- support small-multiple style layouts by hit

### Global Analysis Summary Renderer

Input:

- payload from `prepare_global_analysis_summary_data(...)`

Responsibilities:

- render sample-level or condition-level global summaries
- expose `level="sample"` or `level="condition"`
- optionally show normalization-side quantities from the payload tables

### Global Analysis Window Renderer

Input:

- payload from `prepare_global_analysis_window_data(...)`

Responsibilities:

- render per-contig broad-window summaries
- support both sample-level and condition-level window tables where present

## File Export Helper

Add:

```python
save_figure(fig, path, *, dpi=300, bbox_inches="tight", transparent=False)
```

Responsibilities:

- create parent directories if needed
- call `fig.savefig(...)`
- leave closing policy to the caller

This helper is intentionally thin. It exists for convenience and consistency, not to
become an output manager.

## Testing Strategy

Tests should verify:

- renderer returns Matplotlib figure/axes objects
- renderer handles empty or minimal valid payloads where appropriate
- renderer uses the expected payload table for the selected mode/level
- `save_figure(...)` writes a file successfully

Tests should not try to enforce pixel-perfect rendering.

## Documentation Expectations

User-facing docs should describe the new layer as:

- optional
- Matplotlib-specific
- built on top of the existing prep payloads

The docs should preserve the current message:

- prep helpers are still the canonical interface
- users can still ignore built-in renderers and use their own plotting stack

## Compatibility

This slice should be additive:

- existing payload-prep helpers stay unchanged
- existing legacy Matplotlib plotters remain available
- the new module adds a more coherent built-in rendering surface without removing older paths

## Implementation Order

Recommended order:

1. add `dimelo.plotting_matplotlib` module
2. implement `save_figure(...)`
3. implement region-contrast renderers
4. implement region-discovery renderers
5. implement global-analysis renderers
6. add focused renderer tests
7. update docs and README

