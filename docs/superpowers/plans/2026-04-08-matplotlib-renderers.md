# Matplotlib Renderers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a first built-in Matplotlib renderer layer and figure export helper on top of the existing renderer-neutral plotting payloads for `region_contrasts`, `region_discovery`, and `global_analysis`.

**Architecture:** Keep `dimelo.plotting` as the renderer-neutral prep layer and add a new `dimelo.plotting_matplotlib` module that only consumes prepared payloads and returns Matplotlib figure objects. The renderer layer stays thin: no new biological aggregation, no result-to-payload shortcuts, and only small optional export helpers for writing PNG/PDF outputs.

**Tech Stack:** Python 3.11, matplotlib, pandas, pathlib, pytest, existing `dimelo.plotting`, user docs in `docs/`

---

## File Map

- Create: `dimelo/plotting_matplotlib.py`
  - thin Matplotlib renderers over existing payload dictionaries
  - `save_figure(...)` convenience helper
- Modify: `dimelo/__init__.py`
  - export `plotting_matplotlib`
- Modify: `tests/test_plotting.py`
  - add focused renderer and export coverage
- Modify: `README.md`
  - document the optional Matplotlib renderer layer
- Modify: `docs/region-contrasts.md`
  - add short renderer example for contrast payloads
- Modify: `docs/region-discovery.md`
  - add short renderer example for discovery payloads
- Modify: `docs/global-analysis.md`
  - add short renderer example for global-analysis payloads

## Task 1: Add Failing Tests For The New Matplotlib Module Skeleton

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting_matplotlib.py`

- [ ] **Step 1: Write failing tests for module import and `save_figure(...)`**

Add these tests near the existing plotting coverage in `tests/test_plotting.py`:

```python
def test_plotting_matplotlib_module_exports_save_figure():
    from dimelo import plotting_matplotlib

    assert hasattr(plotting_matplotlib, "save_figure")


def test_save_figure_writes_png(tmp_path):
    from matplotlib import pyplot as plt
    from dimelo import plotting_matplotlib

    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])

    output_path = tmp_path / "figure.png"
    plotting_matplotlib.save_figure(fig, output_path)

    assert output_path.exists()
    assert output_path.stat().st_size > 0
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plotting_matplotlib_module_exports_save_figure or save_figure_writes_png" -q
```

Expected: FAIL because `dimelo.plotting_matplotlib` does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add matplotlib renderer export coverage"
```

## Task 2: Create `dimelo.plotting_matplotlib` And `save_figure(...)`

**Files:**
- Create: `dimelo/plotting_matplotlib.py`
- Modify: `dimelo/__init__.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Create the new module with a thin export helper**

Create `dimelo/plotting_matplotlib.py` with:

```python
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt


def save_figure(
    fig,
    path,
    *,
    dpi: int = 300,
    bbox_inches: str = "tight",
    transparent: bool = False,
) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches=bbox_inches,
        transparent=transparent,
    )
    return output_path
```

- [ ] **Step 2: Export the new module from the package**

Update `dimelo/__init__.py` to import and export `plotting_matplotlib`:

```python
from . import (
    ...
    plotting,
    plotting_matplotlib,
    ...
)

__all__ = [
    ...
    "plotting",
    "plotting_matplotlib",
    ...
]
```

- [ ] **Step 3: Run the focused tests to verify they pass**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plotting_matplotlib_module_exports_save_figure or save_figure_writes_png" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the module skeleton**

```bash
git add dimelo/plotting_matplotlib.py dimelo/__init__.py tests/test_plotting.py
git commit -m "feat: add matplotlib renderer module"
```

## Task 3: Add Failing Tests For Region Contrast Matplotlib Renderers

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting.py`
- Reference: `dimelo/plotting_matplotlib.py`

- [ ] **Step 1: Add failing tests for profile and heatmap renderers**

Add tests using the existing `prepare_region_contrast_profile_data(...)` and
`prepare_region_contrast_heatmap_data(...)` fixtures/patterns:

```python
def test_plot_region_contrast_profile_matplotlib_returns_figure_and_axis():
    from dimelo import plotting, plotting_matplotlib

    result = _mock_region_contrast_result_for_plotting()
    position_table = _mock_region_contrast_position_table()
    payload = plotting.prepare_region_contrast_profile_data(
        result=result,
        position_table=position_table,
        axis=plotting.AxisSpec(
            orientation="region_5to3",
            coordinate_mode="fixed_window",
            anchor="center",
            upstream_bp=100,
            downstream_bp=100,
        ),
        aggregation=plotting.AggregationSpec(),
    )

    fig, ax = plotting_matplotlib.plot_region_contrast_profile_matplotlib(payload)

    assert fig is not None
    assert ax is not None


def test_plot_region_contrast_heatmap_matplotlib_returns_figure_and_axis():
    from dimelo import plotting, plotting_matplotlib

    result = _mock_region_contrast_result_for_plotting()
    position_table = _mock_region_contrast_position_table()
    payload = plotting.prepare_region_contrast_heatmap_data(
        result=result,
        position_table=position_table,
        axis=plotting.AxisSpec(
            orientation="region_5to3",
            coordinate_mode="fixed_window",
            anchor="center",
            upstream_bp=100,
            downstream_bp=100,
        ),
        aggregation=plotting.AggregationSpec(),
    )

    fig, ax = plotting_matplotlib.plot_region_contrast_heatmap_matplotlib(payload)

    assert fig is not None
    assert ax is not None
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_region_contrast_profile_matplotlib or plot_region_contrast_heatmap_matplotlib" -q
```

Expected: FAIL because the renderers do not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add region contrast matplotlib coverage"
```

## Task 4: Implement Region Contrast Matplotlib Renderers

**Files:**
- Modify: `dimelo/plotting_matplotlib.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add minimal helper utilities for Matplotlib payload validation**

Add small internal helpers in `dimelo/plotting_matplotlib.py`:

```python
def _require_payload_keys(payload: dict, keys: tuple[str, ...], *, owner: str) -> None:
    missing = [key for key in keys if key not in payload]
    if missing:
        missing_display = ", ".join(missing)
        raise ValueError(f"{owner} requires payload keys: {missing_display}.")


def _make_axis(ax=None, *, figsize=(6, 4)):
    if ax is not None:
        return ax.figure, ax
    fig, created_ax = plt.subplots(figsize=figsize)
    return fig, created_ax
```

- [ ] **Step 2: Implement `plot_region_contrast_profile_matplotlib(...)`**

Implement a thin line renderer over `payload["plot_table"]`:

```python
def plot_region_contrast_profile_matplotlib(
    payload,
    *,
    value_mode: str = "delta",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(payload, ("plot_table", "metadata"), owner="plot_region_contrast_profile_matplotlib")
    plot_table = payload["plot_table"]
    if value_mode not in {"numerator", "denominator", "delta"}:
        raise ValueError("value_mode must be 'numerator', 'denominator', or 'delta'.")
    value_column = {
        "numerator": "numerator_value",
        "denominator": "denominator_value",
        "delta": "delta_value",
    }[value_mode]
    fig, ax = _make_axis(ax=ax)
    if not plot_table.empty:
        grouped = plot_table.groupby("plot_position", sort=True)[value_column].mean().reset_index()
        ax.plot(grouped["plot_position"], grouped[value_column])
    ax.set_xlabel("Position")
    ax.set_ylabel(value_mode.replace("_", " ").title())
    if title:
        ax.set_title(title)
    return fig, ax
```

- [ ] **Step 3: Implement `plot_region_contrast_heatmap_matplotlib(...)`**

Implement a thin heatmap renderer using `imshow`:

```python
def plot_region_contrast_heatmap_matplotlib(
    payload,
    *,
    value_mode: str = "delta",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(payload, ("plot_table", "metadata"), owner="plot_region_contrast_heatmap_matplotlib")
    plot_table = payload["plot_table"]
    if value_mode not in {"numerator", "denominator", "delta"}:
        raise ValueError("value_mode must be 'numerator', 'denominator', or 'delta'.")
    value_column = {
        "numerator": "numerator_value",
        "denominator": "denominator_value",
        "delta": "delta_value",
    }[value_mode]
    fig, ax = _make_axis(ax=ax, figsize=(7, 5))
    if not plot_table.empty:
        matrix = plot_table.pivot(index="region_id", columns="plot_position", values=value_column)
        image = ax.imshow(matrix.to_numpy(), aspect="auto", interpolation="nearest")
        fig.colorbar(image, ax=ax)
    ax.set_xlabel("Position")
    ax.set_ylabel("Region")
    if title:
        ax.set_title(title)
    return fig, ax
```

- [ ] **Step 4: Run the focused tests to verify they pass**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_region_contrast_profile_matplotlib or plot_region_contrast_heatmap_matplotlib" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the region-contrast renderers**

```bash
git add dimelo/plotting_matplotlib.py tests/test_plotting.py
git commit -m "feat: add region contrast matplotlib renderers"
```

## Task 5: Add Discovery Renderer Tests And Implementations

**Files:**
- Modify: `tests/test_plotting.py`
- Modify: `dimelo/plotting_matplotlib.py`

- [ ] **Step 1: Add failing tests for discovery scan and hit-context renderers**

Add:

```python
def test_plot_region_discovery_scan_matplotlib_returns_figure_and_axes():
    from dimelo import plotting, plotting_matplotlib

    payload = plotting.prepare_region_discovery_scan_data(_mock_region_discovery_result())
    fig, axes = plotting_matplotlib.plot_region_discovery_scan_matplotlib(payload)

    assert fig is not None
    assert axes is not None


def test_plot_region_discovery_hit_context_matplotlib_returns_figure_and_axes():
    from dimelo import plotting, plotting_matplotlib

    payload = plotting.prepare_region_discovery_hit_context_data(
        _mock_region_discovery_result(),
        top_n=2,
    )
    fig, axes = plotting_matplotlib.plot_region_discovery_hit_context_matplotlib(payload)

    assert fig is not None
    assert axes is not None
```

- [ ] **Step 2: Run the focused discovery renderer tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_region_discovery_scan_matplotlib or plot_region_discovery_hit_context_matplotlib" -q
```

Expected: FAIL.

- [ ] **Step 3: Implement the two discovery renderers**

Add thin per-contig and small-multiple renderers:

```python
def plot_region_discovery_scan_matplotlib(payload, *, axes=None, title: str | None = None):
    _require_payload_keys(payload, ("scan_table", "hit_table", "metadata"), owner="plot_region_discovery_scan_matplotlib")
    scan_table = payload["scan_table"]
    contigs = payload["metadata"].get("contig_order") or list(scan_table["contig"].dropna().unique())
    if axes is None:
        fig, axes = plt.subplots(len(contigs) or 1, 1, figsize=(8, max(3, 2 * max(len(contigs), 1))), squeeze=False)
        axes = axes.ravel()
    else:
        fig = axes[0].figure if isinstance(axes, (list, tuple)) else axes.figure
        axes = axes if isinstance(axes, (list, tuple)) else [axes]
    for ax, contig in zip(axes, contigs):
        contig_table = scan_table.loc[scan_table["contig"] == contig]
        ax.plot(contig_table["window_midpoint"], contig_table[payload["metadata"]["score_column"]])
        ax.set_title(str(contig))
    if title:
        fig.suptitle(title)
    return fig, axes


def plot_region_discovery_hit_context_matplotlib(payload, *, axes=None, title: str | None = None):
    _require_payload_keys(payload, ("context_table", "selected_hits", "metadata"), owner="plot_region_discovery_hit_context_matplotlib")
    context_table = payload["context_table"]
    hit_ids = list(context_table["selected_region_id"].dropna().unique()) if not context_table.empty else []
    if axes is None:
        fig, axes = plt.subplots(len(hit_ids) or 1, 1, figsize=(8, max(3, 2 * max(len(hit_ids), 1))), squeeze=False)
        axes = axes.ravel()
    else:
        fig = axes[0].figure if isinstance(axes, (list, tuple)) else axes.figure
        axes = axes if isinstance(axes, (list, tuple)) else [axes]
    score_column = payload["metadata"]["score_column"]
    for ax, hit_id in zip(axes, hit_ids):
        hit_table = context_table.loc[context_table["selected_region_id"] == hit_id]
        ax.plot(hit_table["window_midpoint"], hit_table[score_column])
        ax.set_title(str(hit_id))
    if title:
        fig.suptitle(title)
    return fig, axes
```

- [ ] **Step 4: Run the focused discovery tests to verify they pass**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_region_discovery_scan_matplotlib or plot_region_discovery_hit_context_matplotlib" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the discovery renderers**

```bash
git add dimelo/plotting_matplotlib.py tests/test_plotting.py
git commit -m "feat: add discovery matplotlib renderers"
```

## Task 6: Add Global Analysis Renderer Tests And Implementations

**Files:**
- Modify: `tests/test_plotting.py`
- Modify: `dimelo/plotting_matplotlib.py`

- [ ] **Step 1: Add failing tests for summary and window renderers**

Add:

```python
def test_plot_global_analysis_summary_matplotlib_returns_figure_and_axis():
    from dimelo import plotting, plotting_matplotlib

    payload = plotting.prepare_global_analysis_summary_data(_mock_global_analysis_result())
    fig, ax = plotting_matplotlib.plot_global_analysis_summary_matplotlib(payload)

    assert fig is not None
    assert ax is not None


def test_plot_global_analysis_window_matplotlib_returns_figure_and_axes():
    from dimelo import plotting, plotting_matplotlib

    payload = plotting.prepare_global_analysis_window_data(_mock_global_analysis_result())
    fig, axes = plotting_matplotlib.plot_global_analysis_window_matplotlib(payload)

    assert fig is not None
    assert axes is not None
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_global_analysis_summary_matplotlib or plot_global_analysis_window_matplotlib" -q
```

Expected: FAIL.

- [ ] **Step 3: Implement the two global-analysis renderers**

Add:

```python
def plot_global_analysis_summary_matplotlib(
    payload,
    *,
    level: str = "condition",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(payload, ("sample_summary", "condition_summary", "normalization_table", "metadata"), owner="plot_global_analysis_summary_matplotlib")
    if level not in {"sample", "condition"}:
        raise ValueError("level must be 'sample' or 'condition'.")
    table = payload["sample_summary"] if level == "sample" else payload["condition_summary"]
    fig, ax = _make_axis(ax=ax)
    if not table.empty:
        x_column = "sample_id" if level == "sample" else "condition"
        ax.bar(table[x_column].astype(str), table["global_fraction"])
        ax.tick_params(axis="x", rotation=45)
    ax.set_ylabel("Global Fraction")
    if title:
        ax.set_title(title)
    return fig, ax


def plot_global_analysis_window_matplotlib(payload, *, axes=None, title: str | None = None):
    _require_payload_keys(payload, ("window_table", "metadata"), owner="plot_global_analysis_window_matplotlib")
    window_table = payload["window_table"]
    contigs = payload["metadata"].get("contig_order") or list(window_table["contig"].dropna().unique())
    if axes is None:
        fig, axes = plt.subplots(len(contigs) or 1, 1, figsize=(8, max(3, 2 * max(len(contigs), 1))), squeeze=False)
        axes = axes.ravel()
    else:
        fig = axes[0].figure if isinstance(axes, (list, tuple)) else axes.figure
        axes = axes if isinstance(axes, (list, tuple)) else [axes]
    for ax, contig in zip(axes, contigs):
        contig_table = window_table.loc[window_table["contig"] == contig]
        ax.plot(contig_table["window_midpoint"], contig_table["window_fraction"])
        ax.set_title(str(contig))
    if title:
        fig.suptitle(title)
    return fig, axes
```

- [ ] **Step 4: Run the focused tests to verify they pass**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_global_analysis_summary_matplotlib or plot_global_analysis_window_matplotlib" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the global-analysis renderers**

```bash
git add dimelo/plotting_matplotlib.py tests/test_plotting.py
git commit -m "feat: add global analysis matplotlib renderers"
```

## Task 7: Add User-Facing Docs For The Renderer Layer

**Files:**
- Modify: `README.md`
- Modify: `docs/region-contrasts.md`
- Modify: `docs/region-discovery.md`
- Modify: `docs/global-analysis.md`

- [ ] **Step 1: Add a short README paragraph**

Add a short note near the plotting/helper overview:

```md
For users who want built-in figures instead of only payload tables, `dimelo.plotting_matplotlib`
provides thin Matplotlib renderers on top of the prepared payloads plus `save_figure(...)`
for writing PNG/PDF outputs. The payload-prep helpers remain the canonical interface.
```

- [ ] **Step 2: Add one short example to each analysis guide**

Add renderer examples like:

```python
from dimelo import plotting, plotting_matplotlib

payload = plotting.prepare_region_discovery_scan_data(result)
fig, axes = plotting_matplotlib.plot_region_discovery_scan_matplotlib(payload)
plotting_matplotlib.save_figure(fig, "discovery-scan.png")
```

Repeat the same pattern for:
- `region-contrasts.md`
- `region-discovery.md`
- `global-analysis.md`

- [ ] **Step 3: Commit the docs**

```bash
git add README.md docs/region-contrasts.md docs/region-discovery.md docs/global-analysis.md
git commit -m "docs: add matplotlib renderer examples"
```

## Task 8: Final Verification

**Files:**
- Reference: `dimelo/plotting_matplotlib.py`
- Reference: `tests/test_plotting.py`
- Reference: `README.md`
- Reference: `docs/region-contrasts.md`
- Reference: `docs/region-discovery.md`
- Reference: `docs/global-analysis.md`

- [ ] **Step 1: Run plotting coverage**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/xdg-cache PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

Expected: PASS.

- [ ] **Step 2: Run workflow regression subset**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/xdg-cache PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_workflows.py -q
```

Expected: PASS with no new failures.

- [ ] **Step 3: Run a clean status check**

Run:

```bash
git status --short
```

Expected: clean working tree.

- [ ] **Step 4: Self-review against the spec**

Check the completed work against `docs/superpowers/specs/2026-04-08-matplotlib-renderers-design.md`:

- separate `plotting_matplotlib` module added
- prep layer still lives in `dimelo.plotting`
- renderers cover `region_contrasts`, `region_discovery`, and `global_analysis`
- export helper exists
- no new biological prep logic was moved into renderer functions
- docs describe the renderer layer as optional and Matplotlib-specific

- [ ] **Step 5: Commit any final touch-up if needed**

```bash
git add dimelo/plotting_matplotlib.py dimelo/__init__.py tests/test_plotting.py README.md docs/region-contrasts.md docs/region-discovery.md docs/global-analysis.md
git commit -m "fix: tighten matplotlib renderer support"
```

Only make this commit if verification or self-review uncovers a small follow-up patch.

