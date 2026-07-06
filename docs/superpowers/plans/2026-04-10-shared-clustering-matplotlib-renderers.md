# Shared Clustering Matplotlib Renderers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add built-in Matplotlib renderers for `shared_clustering` distribution/change, cluster profiles, and region occupancy on top of the existing renderer-neutral payloads.

**Architecture:** Keep `dimelo.plotting` as the source of truth for clustering plot payloads and extend `dimelo.plotting_matplotlib` with thin renderers that only consume prepared tables. Allow only narrowly scoped prep-layer adjustments if current payloads are awkward to render cleanly; do not add new clustering outputs or new analysis behavior.

**Tech Stack:** Python 3.11, matplotlib, pandas, pytest, existing `dimelo.plotting`, existing `dimelo.plotting_matplotlib`, user docs in `docs/`

---

## File Map

- Modify: `dimelo/plotting_matplotlib.py`
  - add the five `shared_clustering` Matplotlib renderers
  - reuse existing axis/payload helpers where possible
- Modify: `dimelo/plotting.py`
  - only if a small payload-shape tweak is needed for renderability
- Modify: `tests/test_plotting.py`
  - add focused renderer coverage for the new functions
- Modify: `README.md`
  - document the optional shared-clustering Matplotlib renderer layer
- Modify: `docs/shared-clustering.md`
  - add short renderer examples and clarify default condition-level occupancy view

## Task 1: Add Failing Tests For Shared-Cluster Distribution And Change Renderers

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting.py`
- Reference: `dimelo/plotting_matplotlib.py`

- [ ] **Step 1: Write failing tests for distribution and change renderers**

Add tests near the existing `prepare_shared_cluster_distribution_data(...)` coverage:

```python
def test_plot_shared_cluster_distribution_matplotlib_returns_figure_and_axis():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_distribution_data(
        result=_make_shared_cluster_result()
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_distribution_matplotlib(payload)

    assert fig is not None
    assert ax is not None


def test_plot_shared_cluster_change_matplotlib_returns_figure_and_axis():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_distribution_data(
        result=_make_shared_cluster_result()
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_change_matplotlib(payload)

    assert fig is not None
    assert ax is not None
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_shared_cluster_distribution_matplotlib or plot_shared_cluster_change_matplotlib" -q
```

Expected: FAIL because the two renderer functions do not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add shared clustering distribution renderer coverage"
```

## Task 2: Implement Shared-Cluster Distribution And Change Renderers

**Files:**
- Modify: `dimelo/plotting_matplotlib.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add thin distribution and change renderers**

Add the new public functions in `dimelo/plotting_matplotlib.py`:

```python
def plot_shared_cluster_distribution_matplotlib(
    payload: Mapping[str, object],
    *,
    level: str = "condition",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("sample_distribution", "condition_distribution", "metadata"),
        owner="plot_shared_cluster_distribution_matplotlib",
    )
    if level not in {"sample", "condition"}:
        raise ValueError("level must be 'sample' or 'condition'.")

    table = _require_payload_table(
        payload,
        "sample_distribution" if level == "sample" else "condition_distribution",
    )
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not table.empty:
        x_column = "sample_id" if level == "sample" else "condition"
        pivot = (
            table.pivot_table(
                index=x_column,
                columns="cluster",
                values="fraction",
                aggfunc="mean",
                fill_value=0.0,
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        pivot.plot(kind="bar", stacked=True, ax=ax, legend=True)
        ax.tick_params(axis="x", rotation=45)

    ax.set_xlabel("Sample" if level == "sample" else "Condition")
    ax.set_ylabel("Fraction")
    ax.set_title(title or "Shared cluster distribution")
    return fig, ax


def plot_shared_cluster_change_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("distribution_change", "metadata"),
        owner="plot_shared_cluster_change_matplotlib",
    )
    change_table = _require_payload_table(payload, "distribution_change")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not change_table.empty:
        matrix = (
            change_table.pivot_table(
                index="condition",
                columns="cluster",
                values="delta_fraction",
                aggfunc="mean",
                fill_value=0.0,
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        image = ax.imshow(
            matrix.to_numpy(),
            aspect="auto",
            origin="lower",
            interpolation="nearest",
        )
        ax.figure.colorbar(image, ax=ax, label="Delta Fraction")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels([str(value) for value in matrix.columns.tolist()], rotation=45, ha="right")
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels([str(value) for value in matrix.index.tolist()])

    ax.set_xlabel("Cluster")
    ax.set_ylabel("Condition")
    ax.set_title(title or "Shared cluster change")
    return fig, ax
```

- [ ] **Step 2: Run the focused tests to verify they pass**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_shared_cluster_distribution_matplotlib or plot_shared_cluster_change_matplotlib" -q
```

Expected: PASS.

- [ ] **Step 3: Commit the distribution and change renderers**

```bash
git add dimelo/plotting_matplotlib.py tests/test_plotting.py
git commit -m "feat: add shared clustering distribution renderers"
```

## Task 3: Add Failing Tests For Shared-Cluster Profile Renderers

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting.py`
- Reference: `dimelo/plotting_matplotlib.py`

- [ ] **Step 1: Write failing tests for profile heatmap and profile-series renderers**

Add tests near the existing `prepare_shared_cluster_profile_data(...)` coverage:

```python
def test_plot_shared_cluster_profile_heatmap_matplotlib_returns_figure_and_axis():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_profile_data(
        result=_make_shared_cluster_profile_result()
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_profile_heatmap_matplotlib(payload)

    assert fig is not None
    assert ax is not None


def test_plot_shared_cluster_profile_series_matplotlib_returns_figure_and_axis():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_profile_data(
        result=_make_shared_cluster_profile_result()
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_profile_series_matplotlib(payload)

    assert fig is not None
    assert ax is not None
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_shared_cluster_profile_heatmap_matplotlib or plot_shared_cluster_profile_series_matplotlib" -q
```

Expected: FAIL because the profile renderer functions do not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add shared clustering profile renderer coverage"
```

## Task 4: Implement Shared-Cluster Profile Renderers

**Files:**
- Modify: `dimelo/plotting_matplotlib.py`
- Modify: `dimelo/plotting.py` only if a small metadata tweak is clearly needed
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add the profile heatmap and series renderers**

Implement the two profile renderers in `dimelo/plotting_matplotlib.py`:

```python
def plot_shared_cluster_profile_heatmap_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("profile_table", "metadata"),
        owner="plot_shared_cluster_profile_heatmap_matplotlib",
    )
    profile_table = _require_payload_table(payload, "profile_table")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not profile_table.empty:
        matrix = (
            profile_table.pivot_table(
                index="cluster",
                columns="feature",
                values="value",
                aggfunc="mean",
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        image = ax.imshow(
            matrix.to_numpy(),
            aspect="auto",
            origin="lower",
            interpolation="nearest",
        )
        ax.figure.colorbar(image, ax=ax, label="Value")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels([str(value) for value in matrix.columns.tolist()], rotation=45, ha="right")
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels([str(value) for value in matrix.index.tolist()])

    ax.set_xlabel("Feature")
    ax.set_ylabel("Cluster")
    ax.set_title(title or "Shared cluster profiles")
    return fig, ax


def plot_shared_cluster_profile_series_matplotlib(
    payload: Mapping[str, object],
    *,
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("profile_table", "metadata"),
        owner="plot_shared_cluster_profile_series_matplotlib",
    )
    profile_table = _require_payload_table(payload, "profile_table")
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not profile_table.empty:
        feature_order = list(payload.get("metadata", {}).get("feature_names", []))
        if not feature_order:
            feature_order = sorted(profile_table["feature"].dropna().astype(str).unique())
        feature_lookup = {feature: index for index, feature in enumerate(feature_order)}

        for cluster, cluster_table in profile_table.groupby("cluster", sort=False):
            cluster_table = cluster_table.copy()
            cluster_table["feature_order"] = cluster_table["feature"].map(feature_lookup)
            cluster_table = cluster_table.sort_values("feature_order", kind="stable")
            ax.plot(
                cluster_table["feature_order"],
                cluster_table["value"],
                marker="o",
                linewidth=1.5,
                label=str(cluster),
            )

        ax.set_xticks(range(len(feature_order)))
        ax.set_xticklabels(feature_order, rotation=45, ha="right")
        ax.legend(title="Cluster")

    ax.set_xlabel("Feature")
    ax.set_ylabel("Value")
    ax.set_title(title or "Shared cluster profile series")
    return fig, ax
```

- [ ] **Step 2: Run the focused tests to verify they pass**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_shared_cluster_profile_heatmap_matplotlib or plot_shared_cluster_profile_series_matplotlib" -q
```

Expected: PASS.

- [ ] **Step 3: Commit the profile renderers**

```bash
git add dimelo/plotting_matplotlib.py dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add shared clustering profile renderers"
```

## Task 5: Add Failing Tests For Shared-Cluster Region Renderers

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting.py`
- Reference: `dimelo/plotting_matplotlib.py`

- [ ] **Step 1: Write failing tests for the region occupancy renderer**

Add tests near the existing `prepare_shared_cluster_region_data(...)` coverage:

```python
def test_plot_shared_cluster_region_matplotlib_defaults_to_condition_level():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_region_data(
        result=_make_shared_cluster_region_result()
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_region_matplotlib(payload)

    assert fig is not None
    assert ax is not None


def test_plot_shared_cluster_region_matplotlib_supports_sample_level():
    from dimelo import plotting_matplotlib

    payload = plotting.prepare_shared_cluster_region_data(
        result=_make_shared_cluster_region_result()
    )

    fig, ax = plotting_matplotlib.plot_shared_cluster_region_matplotlib(
        payload,
        level="sample",
    )

    assert fig is not None
    assert ax is not None
```

- [ ] **Step 2: Run the focused tests to verify they fail**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_shared_cluster_region_matplotlib" -q
```

Expected: FAIL because the region renderer does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add shared clustering region renderer coverage"
```

## Task 6: Implement Shared-Cluster Region Renderer And Any Minimal Prep Fix

**Files:**
- Modify: `dimelo/plotting_matplotlib.py`
- Modify: `dimelo/plotting.py` only if a minimal payload tweak is needed
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Implement the condition-default region occupancy renderer**

Implement the renderer in `dimelo/plotting_matplotlib.py`:

```python
def plot_shared_cluster_region_matplotlib(
    payload: Mapping[str, object],
    *,
    level: str = "condition",
    ax=None,
    title: str | None = None,
):
    _require_payload_keys(
        payload,
        ("region_table", "condition_region_table", "metadata"),
        owner="plot_shared_cluster_region_matplotlib",
    )
    if level not in {"sample", "condition"}:
        raise ValueError("level must be 'sample' or 'condition'.")

    table = _require_payload_table(
        payload,
        "condition_region_table" if level == "condition" else "region_table",
    )
    fig, ax = _make_axis(ax=ax, figsize=(8, 4))

    if not table.empty:
        value_column = "fraction_mean" if level == "condition" else "fraction"
        row_columns = ["region_id", "condition"] if level == "condition" else ["region_id", "sample_id"]
        heatmap = table.copy()
        heatmap["row_id"] = heatmap[row_columns].astype(str).agg(" | ".join, axis=1)
        matrix = (
            heatmap.pivot_table(
                index="row_id",
                columns="cluster",
                values=value_column,
                aggfunc="mean",
                fill_value=0.0,
            )
            .sort_index(axis=0)
            .sort_index(axis=1)
        )
        image = ax.imshow(
            matrix.to_numpy(),
            aspect="auto",
            origin="lower",
            interpolation="nearest",
        )
        ax.figure.colorbar(image, ax=ax, label=value_column.replace("_", " ").title())
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels([str(value) for value in matrix.columns.tolist()], rotation=45, ha="right")
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels([str(value) for value in matrix.index.tolist()])

    ax.set_xlabel("Cluster")
    ax.set_ylabel("Region / Condition" if level == "condition" else "Region / Sample")
    ax.set_title(title or "Shared cluster region occupancy")
    return fig, ax
```

- [ ] **Step 2: Add only a minimal prep-layer tweak if the renderer work proves one is necessary**

If the renderer implementation exposes a real payload friction, keep any prep-layer change small and explicit. An acceptable example is adding a stable metadata ordering field in `prepare_shared_cluster_region_data(...)`:

```python
return {
    "region_table": region_table,
    "condition_region_table": condition_region_table,
    "metadata": {
        "mode": result.model.mode,
        "cluster_labels": list(result.model.cluster_labels),
        "has_condition_aggregation": aggregate_conditions,
        "region_order": region_table["region_id"].drop_duplicates().tolist(),
    },
}
```

Do not add new region statistics or new occupancy summaries here.

- [ ] **Step 3: Run the focused tests to verify they pass**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_shared_cluster_region_matplotlib" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the region renderer**

```bash
git add dimelo/plotting_matplotlib.py dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add shared clustering region renderer"
```

## Task 7: Add User-Facing Docs For Shared-Clustering Matplotlib Renderers

**Files:**
- Modify: `README.md`
- Modify: `docs/shared-clustering.md`

- [ ] **Step 1: Add a short README note for shared-clustering renderers**

Add a short note near the existing Matplotlib renderer overview:

```md
`shared_clustering` also supports built-in Matplotlib rendering through
`dimelo.plotting_matplotlib`, including cluster-distribution, profile, and
region-occupancy views on top of the prepared clustering payloads.
```

- [ ] **Step 2: Add short examples to `docs/shared-clustering.md`**

Add examples like:

```python
from dimelo import plotting, plotting_matplotlib

distribution_payload = plotting.prepare_shared_cluster_distribution_data(result=result)
fig, ax = plotting_matplotlib.plot_shared_cluster_distribution_matplotlib(distribution_payload)
plotting_matplotlib.save_figure(fig, "shared-cluster-distribution.png")
```

```python
profile_payload = plotting.prepare_shared_cluster_profile_data(result=result)
fig, ax = plotting_matplotlib.plot_shared_cluster_profile_heatmap_matplotlib(profile_payload)
```

```python
region_payload = plotting.prepare_shared_cluster_region_data(result=result)
fig, ax = plotting_matplotlib.plot_shared_cluster_region_matplotlib(
    region_payload,
    level="condition",
)
```

Keep the docs explicit that:

- condition-level region occupancy is the default rendering view
- sample-level occupancy is still available with `level="sample"`
- the older `result.plot_data` tables remain supported

- [ ] **Step 3: Commit the docs**

```bash
git add README.md docs/shared-clustering.md
git commit -m "docs: add shared clustering renderer examples"
```

## Task 8: Final Verification And Self-Review

**Files:**
- Reference: `dimelo/plotting_matplotlib.py`
- Reference: `dimelo/plotting.py`
- Reference: `tests/test_plotting.py`
- Reference: `README.md`
- Reference: `docs/shared-clustering.md`

- [ ] **Step 1: Run the full plotting coverage**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/xdg-cache PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

Expected: PASS.

- [ ] **Step 2: Run the workflow regression subset**

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

Check the completed work against `docs/superpowers/specs/2026-04-10-shared-clustering-matplotlib-renderers-design.md`:

- separate `shared_clustering` Matplotlib renderers exist in `dimelo.plotting_matplotlib`
- renderers cover distribution/change, profiles, and region occupancy
- default region occupancy view is condition-level
- only minimal prep-layer changes, if any, were introduced
- docs describe the renderer layer as optional and payload-driven

- [ ] **Step 5: Commit a small follow-up only if verification finds a real issue**

```bash
git add dimelo/plotting_matplotlib.py dimelo/plotting.py tests/test_plotting.py README.md docs/shared-clustering.md
git commit -m "fix: tighten shared clustering renderer support"
```

Only make this commit if verification or self-review uncovers a small real follow-up patch.
