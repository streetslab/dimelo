# Global Analysis Plotting Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add renderer-neutral plotting-prep helpers for `dimelo.global_analysis` summary and broad-window views.

**Architecture:** Extend the shared plotting module with two result-centric helpers that consume `GlobalAnalysisResult` and emit stable payload dictionaries. Keep the implementation additive and data-first: no `run_global_analysis()` changes, no renderer coupling, and per-contig broad-window organization by default.

**Tech Stack:** Python 3.11, pandas, pytest, existing `dimelo.plotting`, `dimelo.global_analysis`, and `dimelo.models` modules

---

## File Map

- Modify: `dimelo/plotting.py`
  - add internal validation/helpers for `GlobalAnalysisResult`
  - add `prepare_global_analysis_summary_data(...)`
  - add `prepare_global_analysis_window_data(...)`
- Modify: `tests/test_plotting.py`
  - add focused plotting helper coverage using synthetic `GlobalAnalysisResult` fixtures
- Modify: `docs/global-analysis.md`
  - document the new plotting-prep helpers and their data-first outputs
- Modify: `README.md`
  - add one short note about the new global-analysis plotting helpers

## Task 1: Add Failing Tests For Global Summary Plot Prep

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/models.py`
- Reference: `dimelo/plotting.py`

- [ ] **Step 1: Write the failing summary-prep tests**

Add these tests near the other plotting helper tests in `tests/test_plotting.py`:

```python
from dimelo.models import GlobalAnalysisResult


def _make_global_analysis_result() -> GlobalAnalysisResult:
    summary = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "modified_count": 10,
                "valid_count": 20,
                "global_fraction": 0.5,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "modified_count": 16,
                "valid_count": 20,
                "global_fraction": 0.8,
            },
            {
                "sample_id": "s3",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "modified_count": 12,
                "valid_count": 20,
                "global_fraction": 0.6,
            },
        ]
    )
    normalization_factors = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "global_fraction": 0.5,
                "reference_fraction": 0.6333333333,
                "global_offset": -0.1333333333,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "global_fraction": 0.8,
                "reference_fraction": 0.6333333333,
                "global_offset": 0.1666666667,
            },
            {
                "sample_id": "s3",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "global_fraction": 0.6,
                "reference_fraction": 0.6333333333,
                "global_offset": -0.0333333333,
            },
        ]
    )
    windows = pd.DataFrame(
        columns=[
            "sample_id",
            "condition",
            "replicate",
            "motif",
            "window_id",
            "chromosome",
            "start",
            "end",
            "strand",
            "modified_count",
            "valid_count",
            "window_fraction",
        ]
    )
    return GlobalAnalysisResult(
        summary=summary,
        windows=windows,
        normalization_factors=normalization_factors,
        plot_data={
            "global_fraction_bar": summary.copy(),
            "window_fraction_table": windows.copy(),
        },
        metadata={"window_size": 100000, "step_size": 50000},
    )


def test_prepare_global_analysis_summary_data_returns_expected_tables():
    result = _make_global_analysis_result()

    payload = plotting.prepare_global_analysis_summary_data(result=result)

    assert set(payload) == {"sample_summary", "condition_summary", "normalization_table", "metadata"}
    assert payload["sample_summary"]["sample_id"].tolist() == ["s1", "s2", "s3"]
    assert payload["normalization_table"]["sample_id"].tolist() == ["s1", "s2", "s3"]
    assert payload["condition_summary"]["condition"].tolist() == ["NS", "treated"]
    assert payload["condition_summary"]["sample_n"].tolist() == [1, 2]
    assert payload["metadata"]["motifs"] == ["A,0"]


def test_prepare_global_analysis_summary_data_computes_condition_means():
    result = _make_global_analysis_result()

    payload = plotting.prepare_global_analysis_summary_data(result=result)
    condition_summary = payload["condition_summary"].set_index("condition")

    assert condition_summary.loc["NS", "global_fraction_mean"] == pytest.approx(0.5)
    assert condition_summary.loc["treated", "global_fraction_mean"] == pytest.approx(0.7)
    assert condition_summary.loc["treated", "global_fraction_median"] == pytest.approx(0.7)


def test_prepare_global_analysis_summary_data_filters_motifs():
    result = _make_global_analysis_result()

    payload = plotting.prepare_global_analysis_summary_data(
        result=result,
        motifs=["A,0"],
    )

    assert payload["sample_summary"]["motif"].unique().tolist() == ["A,0"]
    assert payload["normalization_table"]["motif"].unique().tolist() == ["A,0"]
```

- [ ] **Step 2: Run the new summary-prep tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "global_analysis_summary_data" -q
```

Expected: FAIL because `prepare_global_analysis_summary_data(...)` does not yet exist.

- [ ] **Step 3: Commit the failing tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add global analysis summary plotting coverage"
```

## Task 2: Implement Global Summary Plot Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add the internal global-analysis helpers**

Add these internal helpers in `dimelo/plotting.py` near the other result-validation utilities:

```python
def _validate_global_analysis_result(result) -> None:
    if result is None:
        raise ValueError("plotting helpers require a GlobalAnalysisResult.")
    required_attrs = ("summary", "windows", "normalization_factors", "metadata")
    if not all(hasattr(result, attr) for attr in required_attrs):
        raise TypeError("plotting helpers require a GlobalAnalysisResult-like object.")


def _filter_motif_table(table: pd.DataFrame, motifs: list[str] | None, *, owner: str) -> pd.DataFrame:
    if motifs is None:
        return table.copy()
    filtered = table.loc[table["motif"].isin(motifs)].copy()
    if filtered.empty:
        raise ValueError(f"Requested motifs are not present in {owner}.")
    return filtered
```

- [ ] **Step 2: Add `prepare_global_analysis_summary_data(...)`**

Add this public function in `dimelo/plotting.py` near the other public plotting helpers:

```python
def prepare_global_analysis_summary_data(
    *,
    result,
    motifs: list[str] | None = None,
    aggregate_conditions: bool = True,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_global_analysis_result(result)

    sample_summary = _filter_motif_table(result.summary, motifs, owner="result.summary")
    normalization_table = _filter_motif_table(
        result.normalization_factors,
        motifs,
        owner="result.normalization_factors",
    )

    if sample_summary.empty:
        condition_summary = pd.DataFrame(
            columns=[
                "condition",
                "motif",
                "global_fraction_mean",
                "global_fraction_median",
                "sample_n",
            ]
        )
        motif_values: list[str] = [] if motifs is None else list(motifs)
    else:
        motif_values = sample_summary["motif"].drop_duplicates().tolist()
        if aggregate_conditions:
            condition_summary = (
                sample_summary.groupby(["condition", "motif"], as_index=False, sort=False)
                .agg(
                    global_fraction_mean=("global_fraction", "mean"),
                    global_fraction_median=("global_fraction", "median"),
                    sample_n=("sample_id", "nunique"),
                )
                .copy()
            )
        else:
            condition_summary = pd.DataFrame(
                columns=[
                    "condition",
                    "motif",
                    "global_fraction_mean",
                    "global_fraction_median",
                    "sample_n",
                ]
            )

    return {
        "sample_summary": sample_summary.reset_index(drop=True),
        "condition_summary": condition_summary.reset_index(drop=True),
        "normalization_table": normalization_table.reset_index(drop=True),
        "metadata": {
            "motifs": motif_values,
            "aggregate_conditions": aggregate_conditions,
        },
    }
```

- [ ] **Step 3: Run the summary-prep tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "global_analysis_summary_data" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the summary-prep implementation**

```bash
git add dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add global analysis summary plotting prep"
```

## Task 3: Add Failing Tests For Global Window Plot Prep

**Files:**
- Modify: `tests/test_plotting.py`
- Reference: `dimelo/plotting.py`

- [ ] **Step 1: Add the failing window-prep tests**

Append these tests in `tests/test_plotting.py` near the new global summary tests:

```python
def _make_global_window_result() -> GlobalAnalysisResult:
    result = _make_global_analysis_result()
    result.windows = pd.DataFrame(
        [
            {
                "sample_id": "s1",
                "condition": "NS",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-100",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 5,
                "valid_count": 10,
                "window_fraction": 0.5,
            },
            {
                "sample_id": "s2",
                "condition": "treated",
                "replicate": 1,
                "motif": "A,0",
                "window_id": "chr1:0-100",
                "chromosome": "chr1",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 8,
                "valid_count": 10,
                "window_fraction": 0.8,
            },
            {
                "sample_id": "s3",
                "condition": "treated",
                "replicate": 2,
                "motif": "A,0",
                "window_id": "chr2:0-100",
                "chromosome": "chr2",
                "start": 0,
                "end": 100,
                "strand": ".",
                "modified_count": 6,
                "valid_count": 10,
                "window_fraction": 0.6,
            },
        ]
    )
    result.plot_data["window_fraction_table"] = result.windows.copy()
    return result


def test_prepare_global_analysis_window_data_returns_expected_tables():
    result = _make_global_window_result()

    payload = plotting.prepare_global_analysis_window_data(result=result)

    assert set(payload) == {"window_table", "condition_window_table", "metadata"}
    assert payload["window_table"]["contig"].tolist() == ["chr1", "chr1", "chr2"]
    assert payload["window_table"]["window_midpoint"].tolist() == [50.0, 50.0, 50.0]
    assert payload["metadata"]["contig_order"] == ["chr1", "chr2"]


def test_prepare_global_analysis_window_data_filters_contigs_in_requested_order():
    result = _make_global_window_result()

    payload = plotting.prepare_global_analysis_window_data(
        result=result,
        contigs=["chr2"],
    )

    assert payload["metadata"]["contig_order"] == ["chr2"]
    assert payload["window_table"]["contig"].tolist() == ["chr2"]


def test_prepare_global_analysis_window_data_aggregates_conditions():
    result = _make_global_window_result()

    payload = plotting.prepare_global_analysis_window_data(
        result=result,
        aggregate_conditions=True,
    )

    condition_rows = payload["condition_window_table"]
    ns_row = condition_rows.loc[condition_rows["condition"] == "NS"].iloc[0]
    treated_chr1 = condition_rows.loc[
        (condition_rows["condition"] == "treated")
        & (condition_rows["contig"] == "chr1")
    ].iloc[0]

    assert ns_row["window_fraction_mean"] == pytest.approx(0.5)
    assert treated_chr1["window_fraction_mean"] == pytest.approx(0.8)
```

- [ ] **Step 2: Run the new window-prep tests to verify they fail**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "global_analysis_window_data" -q
```

Expected: FAIL because `prepare_global_analysis_window_data(...)` does not yet exist.

- [ ] **Step 3: Commit the failing window-prep tests**

```bash
git add tests/test_plotting.py
git commit -m "test: add global analysis window plotting coverage"
```

## Task 4: Implement Global Window Plot Prep

**Files:**
- Modify: `dimelo/plotting.py`
- Modify: `tests/test_plotting.py`

- [ ] **Step 1: Add `prepare_global_analysis_window_data(...)`**

Add this public function in `dimelo/plotting.py` near the new summary helper:

```python
def prepare_global_analysis_window_data(
    *,
    result,
    contigs: list[str] | None = None,
    motifs: list[str] | None = None,
    aggregate_conditions: bool = False,
) -> dict[str, pd.DataFrame | dict[str, object]]:
    _validate_global_analysis_result(result)

    window_table = _filter_motif_table(result.windows, motifs, owner="result.windows")

    if window_table.empty:
        condition_window_table = pd.DataFrame(
            columns=[
                "condition",
                "motif",
                "contig",
                "start",
                "end",
                "window_midpoint",
                "window_fraction_mean",
                "window_fraction_median",
                "sample_n",
            ]
        )
        return {
            "window_table": pd.DataFrame(
                columns=[
                    "sample_id",
                    "condition",
                    "replicate",
                    "motif",
                    "window_id",
                    "contig",
                    "start",
                    "end",
                    "strand",
                    "window_fraction",
                    "window_midpoint",
                ]
            ),
            "condition_window_table": condition_window_table,
            "metadata": {
                "contig_order": [],
                "motifs": [] if motifs is None else list(motifs),
                "aggregate_conditions": aggregate_conditions,
            },
        }

    window_table = window_table.copy()
    window_table["contig"] = window_table["chromosome"]

    if contigs is not None:
        window_table = window_table.loc[window_table["contig"].isin(contigs)].copy()
        if window_table.empty:
            raise ValueError("Requested contigs are not present in result.windows.")
        contig_order = list(contigs)
    else:
        contig_order = window_table["contig"].drop_duplicates().tolist()

    window_table["window_midpoint"] = (window_table["start"] + window_table["end"]) / 2.0
    window_table = _sort_discovery_table(
        window_table,
        contig_order=contig_order,
        sort_columns=["start", "end", "sample_id"],
    )

    if aggregate_conditions:
        condition_window_table = (
            window_table.groupby(
                ["condition", "motif", "contig", "start", "end", "window_midpoint"],
                as_index=False,
                sort=False,
            )
            .agg(
                window_fraction_mean=("window_fraction", "mean"),
                window_fraction_median=("window_fraction", "median"),
                sample_n=("sample_id", "nunique"),
            )
            .copy()
        )
        condition_window_table = _sort_discovery_table(
            condition_window_table,
            contig_order=contig_order,
            sort_columns=["start", "end", "condition"],
        )
    else:
        condition_window_table = pd.DataFrame(
            columns=[
                "condition",
                "motif",
                "contig",
                "start",
                "end",
                "window_midpoint",
                "window_fraction_mean",
                "window_fraction_median",
                "sample_n",
            ]
        )

    return {
        "window_table": window_table.reset_index(drop=True),
        "condition_window_table": condition_window_table.reset_index(drop=True),
        "metadata": {
            "contig_order": contig_order,
            "motifs": window_table["motif"].drop_duplicates().tolist(),
            "aggregate_conditions": aggregate_conditions,
        },
    }
```

- [ ] **Step 2: Run the window-prep tests to verify they pass**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "global_analysis_window_data" -q
```

Expected: PASS.

- [ ] **Step 3: Run the combined global-analysis plotting subset**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "global_analysis_summary_data or global_analysis_window_data" -q
```

Expected: PASS.

- [ ] **Step 4: Commit the window-prep implementation**

```bash
git add dimelo/plotting.py tests/test_plotting.py
git commit -m "feat: add global analysis window plotting prep"
```

## Task 5: Document The New Global Analysis Plotting Helpers

**Files:**
- Modify: `docs/global-analysis.md`
- Modify: `README.md`
- Verify: `tests/test_plotting.py`

- [ ] **Step 1: Add a short global-analysis plotting section to the guide**

Append a short section to `docs/global-analysis.md` with an example like:

```python
from dimelo import plotting

summary_payload = plotting.prepare_global_analysis_summary_data(
    result=result,
    aggregate_conditions=True,
)

window_payload = plotting.prepare_global_analysis_window_data(
    result=result,
    aggregate_conditions=True,
)
```

Explain briefly:

- summary payloads include sample-level and condition-level views
- normalization values are exposed through `normalization_table`
- window payloads stay per-contig by default
- both helpers are data-prep only and renderer-neutral

- [ ] **Step 2: Add one short README note**

Add one short paragraph to `README.md` near the analysis guides or plotting section explaining:

- `global_analysis` now has renderer-neutral summary and broad-window plotting helpers
- they consume `GlobalAnalysisResult`
- broad-window payloads stay per-contig by default

- [ ] **Step 3: Run the touched plotting test file**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

Expected: PASS.

- [ ] **Step 4: Commit the docs pass**

```bash
git add docs/global-analysis.md README.md tests/test_plotting.py dimelo/plotting.py
git commit -m "docs: add global analysis plotting guide"
```

## Final Verification

**Files:**
- Verify: `dimelo/plotting.py`
- Verify: `tests/test_plotting.py`
- Verify: `docs/global-analysis.md`
- Verify: `README.md`

- [ ] **Step 1: Run the focused plotting verification**

Run:

```bash
MPLCONFIGDIR=/tmp/mpl XDG_CACHE_HOME=/tmp/xdg-cache PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
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
- which payload tables they return
- that summary helpers expose both sample-level and condition-level views
- that broad-window helpers stay per-contig by default

## Self-Review

Spec coverage check:

- summary helper with sample and condition views: covered in Tasks 1-2
- normalization-table exposure: covered in Task 2 and Task 5 docs
- broad-window helper with per-contig default: covered in Tasks 3-4 and Task 5 docs
- additive/data-first behavior: preserved throughout with no `run_global_analysis()` changes

Placeholder scan:

- no `TODO`, `TBD`, or undefined steps remain
- each task includes exact files, commands, and implementation guidance

Type consistency:

- helper names match the approved spec
- payload keys remain aligned with the design:
  - `sample_summary`
  - `condition_summary`
  - `normalization_table`
  - `window_table`
  - `condition_window_table`
  - `metadata`
