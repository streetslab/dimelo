# Global Analysis Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a first `global_analysis` layer that produces pileup-backed global summaries, window-level summaries, and reusable normalization factors without changing the current parsing architecture.

**Architecture:** Build `dimelo.global_analysis` directly on top of `parse_bam.pileup()` outputs and existing `load_processed` pileup parsing helpers. Keep the v1 scope narrow: one result type, one global-summary path, one tiled-window path, and one normalization workflow that later modules can consume.

**Tech Stack:** Python 3.11, pandas, numpy, pysam, pytest

---

## File Map

- `dimelo/models.py`
  Add `GlobalAnalysisResult` as the canonical data-first output contract.
- `dimelo/global_analysis.py`
  New module for whole-pileup summaries, genome tiling, window summaries, and normalization-factor computation.
- `dimelo/__init__.py`
  Export `global_analysis` once the module exists.
- `tests/test_models.py`
  Add model-level validation for `GlobalAnalysisResult`.
- `tests/test_global_analysis.py`
  Add focused tests for summary, tiling, window aggregation, and normalization behavior.
- `docs/global-analysis.md`
  Add the user-facing guide for broad summaries and preprocessing handoff.
- `README.md`
  Link the new guide from the analysis section.

Scope note:

- This plan intentionally covers `global_analysis` foundations only.
- It does not implement `region_discovery`.
- It does not change `parse_bam.pileup()` or `load_processed.pileup_counts_from_bedmethyl()` semantics.
- It keeps plotting data-first and renderer-thin.

### Task 1: Add Global Analysis Result Model

**Files:**
- Modify: `dimelo/models.py`
- Modify: `tests/test_models.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_models.py
from dimelo.models import GlobalAnalysisResult


def test_global_analysis_result_supports_summary_windows_and_normalization():
    result = GlobalAnalysisResult(
        summary=pd.DataFrame({"sample_id": ["s1"]}),
        windows=pd.DataFrame({"window_id": ["chr1:0-1000"]}),
        normalization_factors=pd.DataFrame({"sample_id": ["s1"], "global_offset": [0.1]}),
        plot_data={"global_fraction_bar": pd.DataFrame({"sample_id": ["s1"]})},
        metadata={"normalization_mode": "per_sample_global"},
        figures={},
    )

    assert list(result.summary["sample_id"]) == ["s1"]
    assert list(result.windows["window_id"]) == ["chr1:0-1000"]
    assert list(result.normalization_factors["sample_id"]) == ["s1"]


def test_global_analysis_result_rejects_none_core_outputs():
    with pytest.raises(ValueError, match="summary, normalization_factors, plot_data"):
        GlobalAnalysisResult(
            summary=None,
            windows=None,
            normalization_factors=None,
            plot_data=None,
            metadata={},
            figures={},
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_models.py::test_global_analysis_result_supports_summary_windows_and_normalization tests/test_models.py::test_global_analysis_result_rejects_none_core_outputs -q`
Expected: FAIL with `ImportError` or `AttributeError` because `GlobalAnalysisResult` does not exist yet

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/models.py
@dataclass
class GlobalAnalysisResult:
    summary: pd.DataFrame
    windows: pd.DataFrame | None
    normalization_factors: pd.DataFrame
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "summary": self.summary,
            "normalization_factors": self.normalization_factors,
            "plot_data": self.plot_data,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "GlobalAnalysisResult requires non-None values for: "
                f"{', '.join(missing)}"
            )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_models.py::test_global_analysis_result_supports_summary_windows_and_normalization tests/test_models.py::test_global_analysis_result_rejects_none_core_outputs -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/models.py tests/test_models.py
git commit -m "feat: add global analysis result model"
```

### Task 2: Add Whole-Pileup Global Summary Functions

**Files:**
- Create: `dimelo/global_analysis.py`
- Modify: `dimelo/__init__.py`
- Create: `tests/test_global_analysis.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_global_analysis.py
import pandas as pd

from dimelo import global_analysis
from dimelo.models import SampleSpec


def test_summarize_global_samples_from_pileup(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="s2.h5",
            metadata={"pileup_path": "s2.bed.gz"},
        ),
    ]

    counts = {
        ("s1.bed.gz", "A,0"): (10, 100),
        ("s1.bed.gz", "CG,0"): (20, 200),
        ("s2.bed.gz", "A,0"): (30, 100),
        ("s2.bed.gz", "CG,0"): (40, 200),
    }

    monkeypatch.setattr(
        global_analysis,
        "_global_counts_from_bedmethyl",
        lambda bedmethyl_file, motif, quiet=True: counts[(bedmethyl_file, motif)],
    )

    summary = global_analysis.summarize_global_samples(
        samples=samples,
        motifs=["A,0", "CG,0"],
    )

    assert set(summary.columns) >= {
        "sample_id",
        "condition",
        "motif",
        "modified_count",
        "valid_count",
        "global_fraction",
    }
    assert summary.loc[(summary["sample_id"] == "s2") & (summary["motif"] == "A,0"), "global_fraction"].iloc[0] == 0.3


def test_summarize_global_samples_requires_pileup_path():
    with pytest.raises(ValueError, match="pileup_path"):
        global_analysis.summarize_global_samples(
            samples=[SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5")],
            motifs=["A,0"],
        )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_global_analysis.py::test_summarize_global_samples_from_pileup tests/test_global_analysis.py::test_summarize_global_samples_requires_pileup_path -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'dimelo.global_analysis'`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/global_analysis.py
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd
import pysam

from . import load_processed, utils


def _global_counts_from_bedmethyl(
    bedmethyl_file: str | Path,
    motif: str,
    quiet: bool = True,
) -> tuple[int, int]:
    parsed_motif = utils.ParsedMotif(motif)
    tabix = pysam.TabixFile(str(bedmethyl_file))

    modified_count = 0
    valid_count = 0
    for contig in tabix.contigs:
        for row in tabix.fetch(contig):
            keep_basemod, _, modified_in_row, valid_in_row = load_processed.process_pileup_row(
                row,
                parsed_motif,
                ".",
                False,
            )
            if keep_basemod:
                modified_count += modified_in_row
                valid_count += valid_in_row
    return modified_count, valid_count


def summarize_global_samples(
    *,
    samples,
    motifs: Iterable[str],
    quiet: bool = True,
) -> pd.DataFrame:
    rows = []
    for sample in samples:
        metadata = sample.metadata or {}
        if "pileup_path" not in metadata:
            raise ValueError(
                f"Sample {sample.sample_id} is missing metadata['pileup_path']."
            )
        for motif in motifs:
            modified_count, valid_count = _global_counts_from_bedmethyl(
                metadata["pileup_path"],
                motif,
                quiet=quiet,
            )
            global_fraction = 0.0 if valid_count == 0 else modified_count / valid_count
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "motif": motif,
                    "modified_count": modified_count,
                    "valid_count": valid_count,
                    "global_fraction": global_fraction,
                }
            )
    return pd.DataFrame(rows)
```

```python
# dimelo/__init__.py
from . import global_analysis

__all__.append("global_analysis")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_global_analysis.py::test_summarize_global_samples_from_pileup tests/test_global_analysis.py::test_summarize_global_samples_requires_pileup_path -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/global_analysis.py dimelo/__init__.py tests/test_global_analysis.py
git commit -m "feat: add global pileup summary workflow"
```

### Task 3: Add Window Summaries And Normalization Factors

**Files:**
- Modify: `dimelo/global_analysis.py`
- Modify: `tests/test_global_analysis.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_global_analysis.py
def test_tile_windows_from_genome_sizes_dict():
    windows = global_analysis.tile_genome_windows(
        genome_sizes={"chr1": 2500},
        window_size=1000,
        step_size=500,
    )

    assert windows["window_id"].tolist() == [
        "chr1:0-1000",
        "chr1:500-1500",
        "chr1:1000-2000",
        "chr1:1500-2500",
    ]


def test_build_window_summary_from_regions_to_list(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
    ]
    windows = pd.DataFrame(
        {
            "window_id": ["chr1:0-1000", "chr1:500-1500"],
            "chromosome": ["chr1", "chr1"],
            "start": [0, 500],
            "end": [1000, 1500],
            "strand": [".", "."],
        }
    )

    monkeypatch.setattr(
        global_analysis,
        "tile_genome_windows",
        lambda genome_sizes, window_size, step_size, include_contigs=None, exclude_contigs=None: windows,
    )
    monkeypatch.setattr(
        global_analysis.load_processed,
        "regions_to_list",
        lambda function_handle, regions, window_size, quiet, cores, split_large_regions=False: [(5, 10), (1, 10)],
    )

    summary = global_analysis.build_window_summary(
        samples=samples,
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        window_size=1000,
        step_size=500,
    )

    assert summary["window_fraction"].tolist() == [0.5, 0.1]


def test_compute_global_normalization_factors_from_summary():
    summary = pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "condition": ["NS", "15min"],
            "replicate": [None, None],
            "motif": ["A,0", "A,0"],
            "modified_count": [10, 30],
            "valid_count": [100, 100],
            "global_fraction": [0.1, 0.3],
        }
    )

    factors = global_analysis.compute_global_normalization_factors(summary)

    assert set(factors.columns) >= {
        "sample_id",
        "motif",
        "global_fraction",
        "reference_fraction",
        "global_offset",
    }
    assert factors.loc[factors["sample_id"] == "s1", "global_offset"].iloc[0] == pytest.approx(-0.1)
    assert factors.loc[factors["sample_id"] == "s2", "global_offset"].iloc[0] == pytest.approx(0.1)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_global_analysis.py::test_tile_windows_from_genome_sizes_dict tests/test_global_analysis.py::test_build_window_summary_from_regions_to_list tests/test_global_analysis.py::test_compute_global_normalization_factors_from_summary -q`
Expected: FAIL because the functions do not exist yet

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/global_analysis.py
from functools import partial

import numpy as np


def tile_genome_windows(
    *,
    genome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    include_contigs: Iterable[str] | None = None,
    exclude_contigs: Iterable[str] | None = None,
) -> pd.DataFrame:
    if window_size <= 0 or step_size <= 0:
        raise ValueError("window_size and step_size must be positive.")

    include = set(include_contigs or genome_sizes.keys())
    exclude = set(exclude_contigs or [])
    rows = []
    for chromosome, contig_length in genome_sizes.items():
        if chromosome not in include or chromosome in exclude:
            continue
        start = 0
        while start < contig_length:
            end = min(start + window_size, contig_length)
            rows.append(
                {
                    "window_id": f"{chromosome}:{start}-{end}",
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "strand": ".",
                }
            )
            if end == contig_length:
                break
            start += step_size
    return pd.DataFrame(rows)


def build_window_summary(
    *,
    samples,
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    include_contigs: Iterable[str] | None = None,
    exclude_contigs: Iterable[str] | None = None,
    quiet: bool = True,
    cores: int | None = None,
) -> pd.DataFrame:
    motifs = list(motifs)
    if len(motifs) != 1:
        raise ValueError("build_window_summary currently supports exactly one motif.")

    motif = motifs[0]
    windows = tile_genome_windows(
        genome_sizes=genome_sizes,
        window_size=window_size,
        step_size=step_size,
        include_contigs=include_contigs,
        exclude_contigs=exclude_contigs,
    )
    region_strings = [
        f"{row.chromosome}:{row.start}-{row.end},{row.strand}"
        for row in windows.itertuples(index=False)
    ]

    rows = []
    for sample in samples:
        metadata = sample.metadata or {}
        if "pileup_path" not in metadata:
            raise ValueError(
                f"Sample {sample.sample_id} is missing metadata['pileup_path']."
            )
        pileup_loader = partial(
            load_processed.pileup_counts_from_bedmethyl,
            bedmethyl_file=metadata["pileup_path"],
            motif=motif,
        )
        counts_by_window = load_processed.regions_to_list(
            function_handle=pileup_loader,
            regions=region_strings,
            window_size=None,
            quiet=quiet,
            cores=cores,
        )
        for window_row, (modified_count, valid_count) in zip(
            windows.itertuples(index=False),
            counts_by_window,
        ):
            window_fraction = 0.0 if valid_count == 0 else modified_count / valid_count
            rows.append(
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "motif": motif,
                    "window_id": window_row.window_id,
                    "chromosome": window_row.chromosome,
                    "start": window_row.start,
                    "end": window_row.end,
                    "modified_count": modified_count,
                    "valid_count": valid_count,
                    "window_fraction": window_fraction,
                }
            )
    return pd.DataFrame(rows)


def compute_global_normalization_factors(summary: pd.DataFrame) -> pd.DataFrame:
    grouped_reference = (
        summary.groupby("motif", sort=True)["global_fraction"]
        .mean()
        .rename("reference_fraction")
        .reset_index()
    )
    merged = summary.merge(grouped_reference, on="motif", how="left")
    merged["global_offset"] = merged["global_fraction"] - merged["reference_fraction"]
    return merged.loc[
        :,
        ["sample_id", "condition", "replicate", "motif", "global_fraction", "reference_fraction", "global_offset"],
    ].copy()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_global_analysis.py::test_tile_windows_from_genome_sizes_dict tests/test_global_analysis.py::test_build_window_summary_from_regions_to_list tests/test_global_analysis.py::test_compute_global_normalization_factors_from_summary -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/global_analysis.py tests/test_global_analysis.py
git commit -m "feat: add global window summaries and normalization"
```

### Task 4: Add Global Analysis Result Builder And Docs

**Files:**
- Modify: `dimelo/global_analysis.py`
- Modify: `tests/test_global_analysis.py`
- Create: `docs/global-analysis.md`
- Modify: `README.md`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_global_analysis.py
def test_run_global_analysis_returns_result(monkeypatch):
    samples = [
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="s1.h5",
            metadata={"pileup_path": "s1.bed.gz"},
        ),
    ]

    monkeypatch.setattr(
        global_analysis,
        "summarize_global_samples",
        lambda samples, motifs, quiet=True: pd.DataFrame(
            {
                "sample_id": ["s1"],
                "condition": ["NS"],
                "replicate": [None],
                "motif": ["A,0"],
                "modified_count": [10],
                "valid_count": [100],
                "global_fraction": [0.1],
            }
        ),
    )
    monkeypatch.setattr(
        global_analysis,
        "build_window_summary",
        lambda **kwargs: pd.DataFrame({"window_id": ["chr1:0-1000"]}),
    )

    result = global_analysis.run_global_analysis(
        samples=samples,
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        window_size=1000,
        step_size=1000,
    )

    assert list(result.summary["sample_id"]) == ["s1"]
    assert list(result.windows["window_id"]) == ["chr1:0-1000"]
    assert "global_fraction_bar" in result.plot_data
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_global_analysis.py::test_run_global_analysis_returns_result -q`
Expected: FAIL with `AttributeError` because `run_global_analysis` does not exist yet

- [ ] **Step 3: Write minimal implementation and docs**

```python
# dimelo/global_analysis.py
from .models import GlobalAnalysisResult


def run_global_analysis(
    *,
    samples,
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    window_size: int,
    step_size: int,
    include_contigs: Iterable[str] | None = None,
    exclude_contigs: Iterable[str] | None = None,
    quiet: bool = True,
    cores: int | None = None,
) -> GlobalAnalysisResult:
    summary = summarize_global_samples(samples=samples, motifs=motifs, quiet=quiet)
    windows = build_window_summary(
        samples=samples,
        motifs=motifs,
        genome_sizes=genome_sizes,
        window_size=window_size,
        step_size=step_size,
        include_contigs=include_contigs,
        exclude_contigs=exclude_contigs,
        quiet=quiet,
        cores=cores,
    )
    normalization_factors = compute_global_normalization_factors(summary)
    plot_data = {
        "global_fraction_bar": summary.loc[
            :,
            ["sample_id", "condition", "motif", "global_fraction"],
        ].copy(),
        "window_fraction_table": windows.copy(),
    }
    return GlobalAnalysisResult(
        summary=summary,
        windows=windows,
        normalization_factors=normalization_factors,
        plot_data=plot_data,
        metadata={
            "normalization_mode": "per_sample_global",
            "window_size": window_size,
            "step_size": step_size,
            "motifs": list(motifs),
        },
        figures={},
    )
```

```markdown
# Global Analysis

`dimelo.global_analysis` provides broad pileup-backed summaries on top of `parse_bam.pileup()` outputs.

## What It Covers

- whole-pileup motif abundance summaries
- tiled window summaries across broad genomic space
- reusable per-sample global normalization offsets

## Example

```python
from dimelo import global_analysis
from dimelo.models import SampleSpec

result = global_analysis.run_global_analysis(
    samples=[
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="output/s1.extract.h5",
            metadata={"pileup_path": "output/s1.pileup.sorted.bed.gz"},
        )
    ],
    motifs=["A,0"],
    genome_sizes={"chr1": 248956422},
    window_size=100000,
    step_size=50000,
)
```

The canonical outputs are `result.summary`, `result.windows`, and `result.normalization_factors`.
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_global_analysis.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/global_analysis.py tests/test_global_analysis.py docs/global-analysis.md README.md
git commit -m "docs: add global analysis workflow guide"
```
