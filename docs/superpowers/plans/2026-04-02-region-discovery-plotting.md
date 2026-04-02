# Region Discovery Plotting Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add renderer-neutral plotting-prep helpers for `dimelo.region_discovery` scan overviews and local top-hit context views.

**Architecture:** Extend the shared plotting module with two result-centric helpers that consume `RegionDiscoveryResult` and emit stable payload dictionaries. Keep the implementation additive: no `scan_genome()` signature changes, no renderer coupling, and no normalized region-axis logic for this slice.

**Tech Stack:** Python 3.11, pandas, pytest, existing `dimelo.plotting` and `dimelo.models` modules

---

## File Map

- Modify: `dimelo/plotting.py`
  - add internal validation and extraction helpers for `RegionDiscoveryResult`
  - add `prepare_region_discovery_scan_data(...)`
  - add `prepare_region_discovery_hit_context_data(...)`
- Modify: `tests/test_plotting.py`
  - add focused discovery plotting helper coverage using synthetic `RegionDiscoveryResult` objects
- Modify: `docs/region-discovery.md`
  - document the new plotting-prep helpers and their data-first outputs
- Modify: `README.md`
  - add one short note under plotting/analysis guides for discovery plotting helpers

## Task 1: Add Failing Tests For Discovery Scan Plot Prep

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/models.py`
- Reference: `dimelo/plotting.py`

- [ ] **Step 1: Write the failing scan-prep tests**

Add these tests near the other plotting helper tests in `tests/test_plotting.py`:

```python
from dimelo.models import RegionDiscoveryResult


def _make_region_discovery_result() -> RegionDiscoveryResult:
    windows = pd.DataFrame(
        [
            {
                "window_id": "chr1:0-100:+",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": "+",
                "score_value": 0.2,
                "rank": 3,
            },
            {
                "window_id": "chr1:100-200:+",
                "chromosome": "chr1",
                "start": 100,
                "end": 200,
                "strand": "+",
                "score_value": 0.9,
                "rank": 1,
            },
            {
                "window_id": "chr2:0-100:+",
                "chromosome": "chr2",
                "start": 0,
                "end": 100,
                "strand": "+",
                "score_value": 0.5,
                "rank": 2,
            },
        ]
    )
    hits = windows.loc[windows["rank"] <= 2].copy()
    return RegionDiscoveryResult(
        windows=windows,
        hits=hits,
        plot_data={
            "window_score_table": windows.copy(),
            "top_hits_table": hits.copy(),
        },
        metadata={
            "score": "effect_size_only",
            "contrast_mode": "pairwise",
            "merge_hits": False,
        },
    )


def test_prepare_region_discovery_scan_data_returns_expected_tables():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_scan_data(result=result)

    assert set(payload) == {"scan_table", "hit_table", "metadata"}
    assert list(payload["scan_table"]["contig"]) == ["chr1", "chr1", "chr2"]
    assert list(payload["hit_table"]["rank"]) == [1, 2]
    assert payload["scan_table"]["is_hit"].tolist() == [False, True, True]
    assert payload["metadata"]["contig_order"] == ["chr1", "chr2"]
    assert payload["metadata"]["score_column"] == "score_value"


def test_prepare_region_discovery_scan_data_filters_contigs_in_requested_order():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_scan_data(
        result=result,
        contigs=["chr2", "chr1"],
    )

    assert payload["metadata"]["contig_order"] == ["chr2", "chr1"]
    assert payload["scan_table"]["contig"].tolist() == ["chr2", "chr1", "chr1"]


def test_prepare_region_discovery_scan_data_limits_hit_overlay():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_scan_data(
        result=result,
        top_n_hits=1,
    )

    assert payload["hit_table"]["rank"].tolist() == [1]
    assert payload["scan_table"]["is_hit"].tolist() == [False, True, False]
```

- [ ] **Step 2: Run the new scan-prep tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "region_discovery_scan_data" -q
```

Expected: FAIL because `prepare_region_discovery_scan_data(...)` does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add region discovery scan plotting coverage"
```

## Task 2: Implement Discovery Scan Plot Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add the internal discovery scan helpers**

Add these internal helpers in `dimelo/plotting.py` near the other plot-data validation utilities:

```python
def _validate_region_discovery_result(result):
    if result is None:
        raise ValueError("plotting helpers require a RegionDiscoveryResult.")
    if not hasattr(result, "windows") or not hasattr(result, "hits") or not hasattr(result, "metadata"):
        raise TypeError("plotting helpers require a RegionDiscoveryResult-like object.")


def _select_discovery_score_column(windows: pd.DataFrame, score_column: str | None) -> str:
    if score_column is not None:
        if score_column not in windows.columns:
            raise ValueError(f"Unknown discovery score column: {score_column}")
        return score_column
    if "score_value" in windows.columns:
        return "score_value"
    raise ValueError("Could not infer a discovery score column from RegionDiscoveryResult.windows.")


def _empty_discovery_scan_table(score_column: str) -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "window_id",
            "contig",
            "start",
            "end",
            "strand",
            score_column,
            "window_midpoint",
            "is_hit",
        ]
    )
```

- [ ] **Step 2: Add `prepare_region_discovery_scan_data(...)`**

Add this public function in `dimelo/plotting.py` near the other public plotting helpers:

```python
def prepare_region_discovery_scan_data(
    *,
    result,
    contigs: list[str] | None = None,
    top_n_hits: int | None = 100,
    score_column: str | None = None,
    include_all_windows: bool = True,
):
    _validate_region_discovery_result(result)

    windows = result.windows.copy()
    hits = result.hits.copy()
    active_score_column = _select_discovery_score_column(windows, score_column)

    if windows.empty:
        empty_scan = _empty_discovery_scan_table(active_score_column)
        empty_hits = pd.DataFrame(columns=["window_id", "contig", "start", "end", "strand", active_score_column])
        return {
            "scan_table": empty_scan,
            "hit_table": empty_hits,
            "metadata": {
                "contig_order": [],
                "score_column": active_score_column,
                "score_mode": result.metadata.get("score"),
                "contrast_mode": result.metadata.get("contrast_mode"),
                "merge_hits": result.metadata.get("merge_hits"),
                "top_n_hits": top_n_hits,
            },
        }

    windows["contig"] = windows["chromosome"]
    hits["contig"] = hits["chromosome"]

    if contigs is not None:
        windows = windows.loc[windows["contig"].isin(contigs)].copy()
        hits = hits.loc[hits["contig"].isin(contigs)].copy()
        if windows.empty:
            raise ValueError("Requested contigs are not present in RegionDiscoveryResult.windows.")
        contig_order = list(contigs)
    else:
        contig_order = windows["contig"].drop_duplicates().tolist()

    if top_n_hits is not None and not hits.empty and "rank" in hits.columns:
        hits = hits.sort_values("rank", kind="stable").head(top_n_hits).copy()

    hit_window_ids = set(hits.get("window_id", pd.Series(dtype="object")).tolist())
    windows["window_midpoint"] = (windows["start"] + windows["end"]) / 2.0
    windows["is_hit"] = windows["window_id"].isin(hit_window_ids)

    scan_table = windows if include_all_windows else windows.loc[windows["is_hit"]].copy()
    scan_table = scan_table.sort_values(["contig", "start", "end"], kind="stable").reset_index(drop=True)
    hit_table = hits.sort_values(["rank", "contig", "start", "end"], kind="stable").reset_index(drop=True)

    return {
        "scan_table": scan_table,
        "hit_table": hit_table,
        "metadata": {
            "contig_order": contig_order,
            "score_column": active_score_column,
            "score_mode": result.metadata.get("score"),
            "contrast_mode": result.metadata.get("contrast_mode"),
            "merge_hits": result.metadata.get("merge_hits"),
            "top_n_hits": top_n_hits,
        },
    }
```

- [ ] **Step 3: Run the scan-prep tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "region_discovery_scan_data" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the scan-prep implementation**

```bash
git add dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add region discovery scan plotting prep"
```

## Task 3: Add Failing Tests For Hit-Context Plot Prep

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting.py`

- [ ] **Step 1: Add the failing hit-context tests**

Append these tests in `tests/test_plotting.py` near the new discovery scan tests:

```python
def test_prepare_region_discovery_hit_context_data_uses_top_ranked_hits():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_hit_context_data(
        result=result,
        top_n=1,
        padding_windows=1,
    )

    assert set(payload) == {"context_table", "selected_hits", "metadata"}
    assert payload["selected_hits"]["rank"].tolist() == [1]
    assert payload["context_table"]["selected_hit_rank"].nunique() == 1
    assert payload["context_table"]["selected_hit_rank"].iloc[0] == 1
    assert payload["context_table"]["is_selected_hit"].sum() == 1


def test_prepare_region_discovery_hit_context_data_adds_relative_window_offsets():
    result = _make_region_discovery_result()

    payload = plotting.prepare_region_discovery_hit_context_data(
        result=result,
        top_n=1,
        padding_windows=1,
    )

    assert sorted(payload["context_table"]["relative_window_offset"].tolist()) == [-1, 0]


def test_prepare_region_discovery_hit_context_data_returns_empty_payload_for_no_hits():
    result = RegionDiscoveryResult(
        windows=_make_region_discovery_result().windows,
        hits=pd.DataFrame(columns=_make_region_discovery_result().hits.columns),
        plot_data={},
        metadata={"score": "effect_size_only", "contrast_mode": "pairwise", "merge_hits": False},
    )

    payload = plotting.prepare_region_discovery_hit_context_data(result=result)

    assert payload["context_table"].empty
    assert payload["selected_hits"].empty
    assert payload["metadata"]["selection_mode"] == "top_n"
```

- [ ] **Step 2: Run the hit-context tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "region_discovery_hit_context_data" -q
```

Expected: FAIL because `prepare_region_discovery_hit_context_data(...)` does not yet exist.

- [ ] **Step 3: Commit the failing hit-context tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add region discovery hit context coverage"
```

## Task 4: Implement Discovery Hit-Context Plot Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add the hit selection helper**

Insert this helper in `dimelo/plotting.py` near the discovery scan helpers:

```python
def _select_region_discovery_hits(
    hits: pd.DataFrame,
    *,
    top_n: int | None,
):
    if hits.empty:
        return hits.copy(), "top_n"
    if "rank" in hits.columns:
        selected = hits.sort_values("rank", kind="stable")
    else:
        selected = hits.copy()
    if top_n is not None:
        selected = selected.head(top_n)
    return selected.copy(), "top_n"
```

- [ ] **Step 2: Add `prepare_region_discovery_hit_context_data(...)`**

Add this public function in `dimelo/plotting.py` near the scan helper:

```python
def prepare_region_discovery_hit_context_data(
    *,
    result,
    top_n: int | None = 12,
    hit_ids: list[str] | None = None,
    padding_windows: int | None = 5,
    padding_bp: int | None = None,
    score_column: str | None = None,
):
    _validate_region_discovery_result(result)

    windows = result.windows.copy()
    hits = result.hits.copy()
    active_score_column = _select_discovery_score_column(windows, score_column)

    if hits.empty or windows.empty:
        return {
            "context_table": pd.DataFrame(
                columns=[
                    "window_id",
                    "contig",
                    "start",
                    "end",
                    "strand",
                    active_score_column,
                    "window_midpoint",
                    "selected_hit_id",
                    "selected_hit_rank",
                    "relative_window_offset",
                    "is_selected_hit",
                ]
            ),
            "selected_hits": hits.copy(),
            "metadata": {
                "selection_mode": "top_n" if hit_ids is None else "explicit_ids",
                "top_n": top_n,
                "padding_windows": padding_windows,
                "padding_bp": padding_bp,
                "score_column": active_score_column,
            },
        }

    windows = windows.sort_values(["chromosome", "start", "end"], kind="stable").reset_index(drop=True)
    windows["contig"] = windows["chromosome"]
    windows["window_midpoint"] = (windows["start"] + windows["end"]) / 2.0

    if hit_ids is not None:
        selected_hits = hits.loc[hits["window_id"].isin(hit_ids)].copy()
        selection_mode = "explicit_ids"
    else:
        selected_hits, selection_mode = _select_region_discovery_hits(hits, top_n=top_n)

    if selected_hits.empty:
        return {
            "context_table": pd.DataFrame(columns=list(windows.columns) + ["selected_hit_id", "selected_hit_rank", "relative_window_offset", "is_selected_hit"]),
            "selected_hits": selected_hits,
            "metadata": {
                "selection_mode": selection_mode,
                "top_n": top_n,
                "padding_windows": padding_windows,
                "padding_bp": padding_bp,
                "score_column": active_score_column,
            },
        }

    context_frames = []
    for _, hit in selected_hits.sort_values("rank", kind="stable").iterrows():
        same_contig = windows.loc[windows["contig"] == hit["chromosome"]].reset_index(drop=True)
        hit_positions = same_contig.index[same_contig["window_id"] == hit["window_id"]].tolist()
        if not hit_positions:
            continue
        hit_position = hit_positions[0]
        start_index = max(0, hit_position - int(padding_windows or 0))
        end_index = min(len(same_contig), hit_position + int(padding_windows or 0) + 1)
        context = same_contig.iloc[start_index:end_index].copy()
        context["selected_hit_id"] = hit["window_id"]
        context["selected_hit_rank"] = hit.get("rank")
        context["relative_window_offset"] = range(start_index - hit_position, end_index - hit_position)
        context["is_selected_hit"] = context["window_id"] == hit["window_id"]
        context_frames.append(context)

    context_table = (
        pd.concat(context_frames, ignore_index=True)
        if context_frames
        else pd.DataFrame(columns=list(windows.columns) + ["selected_hit_id", "selected_hit_rank", "relative_window_offset", "is_selected_hit"])
    )

    return {
        "context_table": context_table,
        "selected_hits": selected_hits.sort_values("rank", kind="stable").reset_index(drop=True),
        "metadata": {
            "selection_mode": selection_mode,
            "top_n": top_n,
            "padding_windows": padding_windows,
            "padding_bp": padding_bp,
            "score_column": active_score_column,
        },
    }
```

- [ ] **Step 3: Run the hit-context tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "region_discovery_hit_context_data" -q
```

Expected: PASS.

- [ ] **Step 4: Run the combined discovery plotting subset**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "region_discovery_scan_data or region_discovery_hit_context_data" -q
```

Expected: PASS.

- [ ] **Step 5: Commit the hit-context implementation**

```bash
git add dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add region discovery hit context plotting prep"
```

## Task 5: Document The New Discovery Plotting Helpers

**Files:**
- Modify: `docs/region-discovery.md`
- Modify: `README.md`
- Verify: `tests/test_plotting.py`

- [ ] **Step 1: Add a short discovery plotting section to the guide**

Append a short section to `docs/region-discovery.md` with an example like:

```python
from dimelo import plotting

scan_payload = plotting.prepare_region_discovery_scan_data(
    result=result,
    top_n_hits=50,
)

context_payload = plotting.prepare_region_discovery_hit_context_data(
    result=result,
    top_n=12,
    padding_windows=5,
)
```

Explain briefly:

- scan payloads stay per-contig by default
- hit-context payloads are intended for small multiple panels
- both helpers are data-prep only and renderer-neutral

- [ ] **Step 2: Add one short README note**

Add one short paragraph to `README.md` under the plotting or analysis-guides area explaining:

- `region_discovery` now has renderer-neutral scan and hit-context plotting helpers
- they consume `RegionDiscoveryResult`
- they prepare per-contig scan tables rather than cumulative genome axes by default

- [ ] **Step 3: Run the touched plotting test file**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

Expected: PASS.

- [ ] **Step 4: Commit the docs pass**

```bash
git add docs/region-discovery.md README.md tests/test_plotting.py dimelo/plotting.py
git commit -m "docs: add region discovery plotting guide"
```

## Final Verification

**Files:**
- Verify: `dimelo/plotting.py`
- Verify: `tests/test_plotting.py`
- Verify: `docs/region-discovery.md`
- Verify: `README.md`

- [ ] **Step 1: Run the focused plotting verification**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

Expected: PASS.

- [ ] **Step 2: Confirm the working tree is clean**

Run:

```bash
git status --short
```

Expected: no output.

- [ ] **Step 3: Summarize outcomes**

Capture in the handoff:

- which public helpers were added
- which tests were added
- that discovery plotting remains data-first and per-contig by default

## Self-Review

Spec coverage check:

- scan-overview helper: covered in Tasks 1-2
- hit-context helper: covered in Tasks 3-4
- per-contig default and metadata: covered in Tasks 1-2 and Task 5 docs
- additive/no-renderer design: preserved throughout, with no `scan_genome()` changes

Placeholder scan:

- no `TODO`, `TBD`, or undefined implementation steps remain
- each task includes exact files, code snippets, commands, and expected outcomes

Type consistency:

- both helpers are named consistently with the approved spec
- `RegionDiscoveryResult` remains the public input
- payload keys stay aligned with the spec: `scan_table`, `hit_table`, `context_table`, `selected_hits`, `metadata`
