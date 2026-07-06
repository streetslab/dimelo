# Region Contrasts Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a stable, backward-compatible `region_contrasts` layer for known-region comparisons over pileup-backed counts, including explicit contrast models, effect-size summaries, and a first beta-binomial testing path.

**Architecture:** Build `dimelo.region_contrasts` on top of the existing `parse_bam.pileup()` and `load_processed.pileup_counts_from_bedmethyl()` flow without changing parsing semantics. Keep the v1 scope intentionally narrow: fully support only `analysis_unit="ensemble_region"` + `signal_source="pileup_counts"` with `effect_size_only` and `beta_binomial`, while validating and rejecting unsupported combinations explicitly.

**Tech Stack:** Python 3.11, pandas, numpy, scipy, pytest

---

Scope note:

- This plan intentionally covers the first executable `region_contrasts` slice only.
- It does not implement `region_discovery`, single-read contrasts, cluster-occupancy contrasts, or time-course inference.
- It preserves continuity by consuming existing pileup outputs rather than altering `parse_bam` or `load_processed` semantics.

### Task 1: Add Contrast Models And Result Types

**Files:**
- Modify: `dimelo/models.py`
- Modify: `tests/test_models.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_models.py
from dimelo.models import ContrastSpec, RegionContrastResult


def test_contrast_spec_accepts_pairwise_mode():
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["15min"],
        denominator=["NS"],
        reference_condition="NS",
    )
    assert contrast.mode == "pairwise"
    assert contrast.numerator == ["15min"]


def test_contrast_spec_rejects_missing_groups_for_pairwise():
    try:
        ContrastSpec(mode="pairwise")
    except ValueError as exc:
        assert "numerator" in str(exc)
    else:
        raise AssertionError("Expected ContrastSpec to validate pairwise groups")


def test_region_contrast_result_requires_regions_and_summary():
    contrast = ContrastSpec(
        mode="pairwise",
        numerator=["15min"],
        denominator=["NS"],
        reference_condition="NS",
    )
    try:
        RegionContrastResult(
            regions=None,
            summary=None,
            contrast=contrast,
            plot_data={},
            metadata={},
            figures={},
        )
    except ValueError as exc:
        assert "regions" in str(exc)
        assert "summary" in str(exc)
    else:
        raise AssertionError("Expected RegionContrastResult to require core tables")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_models.py -q`
Expected: FAIL with `ImportError` or `AttributeError` for missing `ContrastSpec` / `RegionContrastResult`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/models.py
from dataclasses import dataclass, field
from typing import Any


@dataclass
class ContrastSpec:
    mode: str
    numerator: list[str] | None = None
    denominator: list[str] | None = None
    background: list[str] | None = None
    time_order: list[str] | None = None
    pairing_key: str | None = None
    reference_condition: str | None = None
    metadata: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if self.mode not in {
            "single_dataset",
            "pairwise",
            "matched_pairwise",
            "group_vs_group",
            "background_adjusted",
            "time_course",
        }:
            raise ValueError(f"Unsupported contrast mode: {self.mode}")
        if self.mode in {"pairwise", "group_vs_group"}:
            if not self.numerator or not self.denominator:
                raise ValueError(
                    "ContrastSpec pairwise/group_vs_group modes require numerator and denominator."
                )
        if self.mode == "time_course" and not self.time_order:
            raise ValueError("ContrastSpec time_course mode requires time_order.")


@dataclass
class RegionContrastResult:
    regions: pd.DataFrame
    summary: pd.DataFrame
    contrast: ContrastSpec
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "regions": self.regions,
            "summary": self.summary,
            "plot_data": self.plot_data,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "RegionContrastResult requires non-None values for: "
                f"{', '.join(missing)}"
            )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_models.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/models.py tests/test_models.py
git commit -m "feat: add region contrast models"
```

### Task 2: Add Region Contrast Validation And Evidence Table Construction

**Files:**
- Create: `dimelo/region_contrasts.py`
- Create: `tests/test_region_contrasts.py`
- Modify: `dimelo/__init__.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_region_contrasts.py
import pandas as pd
import pytest

from dimelo.models import ContrastSpec, SampleSpec
from dimelo import region_contrasts


def test_validate_supported_v1_combination():
    region_contrasts.validate_region_contrast_request(
        analysis_unit="ensemble_region",
        representation="modified_fraction",
        signal_source="pileup_counts",
        test="beta_binomial",
    )


def test_validate_rejects_unsupported_single_read_beta_binomial():
    with pytest.raises(ValueError, match="ensemble_region"):
        region_contrasts.validate_region_contrast_request(
            analysis_unit="single_read",
            representation="read_mod_fraction",
            signal_source="pileup_counts",
            test="beta_binomial",
        )


def test_build_region_evidence_table_from_pileup_counts(monkeypatch):
    fake_samples = [
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5", metadata={"pileup_path": "s1.bed.gz"}),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5", metadata={"pileup_path": "s2.bed.gz"}),
    ]

    monkeypatch.setattr(
        region_contrasts.utils,
        "regions_dict_from_input",
        lambda regions, window_size=None: {"chr1": [(0, 10, "+")], "chr2": [(20, 30, "-")]},
    )
    counts_by_file = {
        "s1.bed.gz": [(2, 10), (6, 10)],
        "s2.bed.gz": [(7, 10), (8, 10)],
    }
    monkeypatch.setattr(
        region_contrasts.load_processed,
        "regions_to_list",
        lambda function_handle, regions, window_size, quiet, cores, split_large_regions=False: counts_by_file[function_handle.keywords["bedmethyl_file"]],
    )

    evidence = region_contrasts.build_region_evidence_table(
        samples=fake_samples,
        regions="matched.bed",
        motifs=["A,0"],
    )

    assert list(evidence.columns[:7]) == [
        "region_id",
        "chromosome",
        "start",
        "end",
        "strand",
        "sample_id",
        "condition",
    ]
    assert set(evidence["modified_count"]) == {2, 6, 7, 8}
    assert set(evidence["valid_count"]) == {10}
    assert "mod_fraction" in evidence.columns
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_region_contrasts.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'dimelo.region_contrasts'`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/region_contrasts.py
from __future__ import annotations

from functools import partial

import pandas as pd

from . import load_processed, utils


def validate_region_contrast_request(
    *,
    analysis_unit: str,
    representation: str,
    signal_source: str,
    test: str,
) -> None:
    if analysis_unit != "ensemble_region":
        raise ValueError("V1 inferential region_contrasts requires analysis_unit='ensemble_region'.")
    if signal_source != "pileup_counts":
        raise ValueError("V1 inferential region_contrasts requires signal_source='pileup_counts'.")
    if representation not in {"modified_fraction", "modified_count"}:
        raise ValueError("V1 inferential region_contrasts requires modified_fraction or modified_count.")
    if test not in {"beta_binomial", "effect_size_only"}:
        raise ValueError("Unsupported region contrast test.")


def build_region_evidence_table(
    *,
    samples,
    regions,
    motifs,
    window_size=None,
    quiet: bool = True,
    cores=None,
):
    if len(motifs) != 1:
        raise ValueError("build_region_evidence_table currently supports exactly one motif.")
    motif = motifs[0]
    regions_dict = utils.regions_dict_from_input(regions, window_size)
    region_metadata = [
        (chromosome, start, end, strand)
        for chromosome, region_list in regions_dict.items()
        for start, end, strand in region_list
    ]
    rows = []
    for sample in samples:
        if not sample.metadata or "pileup_path" not in sample.metadata:
            raise ValueError(f"Sample {sample.sample_id!r} is missing metadata['pileup_path'].")
        loader = partial(
            load_processed.pileup_counts_from_bedmethyl,
            bedmethyl_file=sample.metadata["pileup_path"],
            motif=motif,
            window_size=window_size,
            quiet=quiet,
            cores=cores,
        )
        per_region_counts = load_processed.regions_to_list(
            function_handle=loader,
            regions=regions,
            window_size=window_size,
            quiet=quiet,
            cores=cores,
            split_large_regions=False,
        )
        for (chrom, start, end, strand), (modified_count, valid_count) in zip(
            region_metadata, per_region_counts
        ):
            mod_fraction = 0.0 if valid_count == 0 else modified_count / valid_count
            rows.append(
                {
                    "region_id": f"{chrom}:{start}-{end}:{strand}",
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "modified_count": int(modified_count),
                    "valid_count": int(valid_count),
                    "mod_fraction": float(mod_fraction),
                }
            )
    return pd.DataFrame(rows)
```

```python
# dimelo/__init__.py
from . import (
    cluster,
    distribution,
    export,
    load_processed,
    models,
    parse_bam,
    plotting,
    plot_depth_histogram,
    plot_depth_profile,
    plot_enrichment,
    plot_enrichment_profile,
    plot_read_browser,
    plot_reads,
    region_analysis,
    region_contrasts,
    workflows,
)

__all__ = [
    "cluster",
    "distribution",
    "export",
    "load_processed",
    "models",
    "parse_bam",
    "plotting",
    "plot_depth_histogram",
    "plot_depth_profile",
    "plot_enrichment",
    "plot_enrichment_profile",
    "plot_read_browser",
    "plot_reads",
    "region_analysis",
    "region_contrasts",
    "workflows",
]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_region_contrasts.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/region_contrasts.py dimelo/__init__.py tests/test_region_contrasts.py
git commit -m "feat: add region contrast evidence builders"
```

### Task 3: Add Effect-Size Region Scoring Workflow

**Files:**
- Modify: `dimelo/region_contrasts.py`
- Modify: `tests/test_region_contrasts.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_region_contrasts.py
from dimelo.models import ContrastSpec


def test_score_regions_effect_size_pairwise(monkeypatch):
    evidence = pd.DataFrame(
        {
            "region_id": ["reg1", "reg1", "reg2", "reg2"],
            "chromosome": ["chr1", "chr1", "chr2", "chr2"],
            "start": [0, 0, 20, 20],
            "end": [10, 10, 30, 30],
            "strand": ["+", "+", "-", "-"],
            "sample_id": ["s1", "s2", "s1", "s2"],
            "condition": ["NS", "15min", "NS", "15min"],
            "replicate": [1, 1, 1, 1],
            "modified_count": [2, 8, 6, 6],
            "valid_count": [10, 10, 10, 10],
            "mod_fraction": [0.2, 0.8, 0.6, 0.6],
        }
    )
    monkeypatch.setattr(region_contrasts, "build_region_evidence_table", lambda **kwargs: evidence)

    result = region_contrasts.score_regions(
        samples=[],
        regions="matched.bed",
        motifs=["A,0"],
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["15min"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        analysis_unit="ensemble_region",
        representation="modified_fraction",
        signal_source="pileup_counts",
        test="effect_size_only",
    )

    assert set(result.regions.columns) >= {
        "region_id",
        "fraction",
        "reference_fraction",
        "delta_fraction",
        "log2_fc",
        "rank",
    }
    top = result.regions.sort_values("rank").iloc[0]
    assert top["region_id"] == "reg1"
    assert top["delta_fraction"] > 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_region_contrasts.py::test_score_regions_effect_size_pairwise -q`
Expected: FAIL with `AttributeError` for missing `score_regions`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/region_contrasts.py
from .models import ContrastSpec, RegionContrastResult


def _pool_region_groups(evidence: pd.DataFrame, contrast: ContrastSpec) -> pd.DataFrame:
    if contrast.mode not in {"pairwise", "group_vs_group"}:
        raise NotImplementedError(
            "V1 effect-size region_contrasts supports pairwise and group_vs_group only."
        )
    condition_groups = {
        "numerator": set(contrast.numerator or []),
        "denominator": set(contrast.denominator or []),
    }
    pooled_rows = []
    for label, conditions in condition_groups.items():
        subset = evidence[evidence["condition"].isin(conditions)].copy()
        grouped = (
            subset.groupby(["region_id", "chromosome", "start", "end", "strand"], sort=True)
            .agg(
                modified_count=("modified_count", "sum"),
                valid_count=("valid_count", "sum"),
                replicate_n=("sample_id", "nunique"),
            )
            .reset_index()
        )
        grouped["group"] = label
        grouped["fraction"] = grouped["modified_count"] / grouped["valid_count"].where(
            grouped["valid_count"] > 0,
            1,
        )
        pooled_rows.append(grouped)
    return pd.concat(pooled_rows, ignore_index=True)


def score_regions(
    *,
    samples,
    regions,
    motifs,
    contrast: ContrastSpec,
    analysis_unit: str = "ensemble_region",
    representation: str = "modified_fraction",
    signal_source: str = "pileup_counts",
    test: str = "beta_binomial",
    multiple_testing: str = "fdr_bh",
):
    validate_region_contrast_request(
        analysis_unit=analysis_unit,
        representation=representation,
        signal_source=signal_source,
        test=test,
    )
    evidence = build_region_evidence_table(samples=samples, regions=regions, motifs=motifs)
    pooled = _pool_region_groups(evidence, contrast)
    numerator = pooled[pooled["group"] == "numerator"].rename(
        columns={
            "modified_count": "modified_count",
            "valid_count": "valid_count",
            "fraction": "fraction",
            "replicate_n": "replicate_n",
        }
    )
    denominator = pooled[pooled["group"] == "denominator"].rename(
        columns={
            "modified_count": "reference_modified_count",
            "valid_count": "reference_valid_count",
            "fraction": "reference_fraction",
            "replicate_n": "reference_replicate_n",
        }
    )
    merged = numerator.merge(
        denominator,
        on=["region_id", "chromosome", "start", "end", "strand"],
        how="left",
    )
    merged["delta_fraction"] = merged["fraction"] - merged["reference_fraction"]
    merged["log2_fc"] = np.log2(
        (merged["fraction"] + 1e-6) / (merged["reference_fraction"] + 1e-6)
    )
    merged["rank"] = (
        merged["delta_fraction"].abs().rank(method="first", ascending=False).astype(int)
    )
    summary = merged.loc[:, ["region_id", "fraction", "reference_fraction", "delta_fraction", "log2_fc", "rank"]]
    return RegionContrastResult(
        regions=merged.sort_values("rank", kind="stable").reset_index(drop=True),
        summary=summary.sort_values("rank", kind="stable").reset_index(drop=True),
        contrast=contrast,
        plot_data={"region_effect_sizes": summary.copy()},
        metadata={
            "contrast_mode": contrast.mode,
            "analysis_unit": analysis_unit,
            "representation": representation,
            "signal_source": signal_source,
            "test": test,
            "normalization_mode": "none",
            "biological_interpretation": "Region-level motif abundance contrast.",
            "renderer": None,
        },
        figures={},
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_region_contrasts.py::test_score_regions_effect_size_pairwise -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add effect-size region contrast workflow"
```

### Task 4: Add Beta-Binomial Testing For Pooled Region Counts

**Files:**
- Modify: `dimelo/region_contrasts.py`
- Modify: `tests/test_region_contrasts.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_region_contrasts.py

def test_score_regions_beta_binomial_adds_pvalues(monkeypatch):
    evidence = pd.DataFrame(
        {
            "region_id": ["reg1", "reg1", "reg2", "reg2"],
            "chromosome": ["chr1", "chr1", "chr2", "chr2"],
            "start": [0, 0, 20, 20],
            "end": [10, 10, 30, 30],
            "strand": ["+", "+", "-", "-"],
            "sample_id": ["s1", "s2", "s1", "s2"],
            "condition": ["NS", "15min", "NS", "15min"],
            "replicate": [1, 1, 1, 1],
            "modified_count": [1, 9, 5, 6],
            "valid_count": [10, 10, 10, 10],
            "mod_fraction": [0.1, 0.9, 0.5, 0.6],
        }
    )
    monkeypatch.setattr(region_contrasts, "build_region_evidence_table", lambda **kwargs: evidence)

    result = region_contrasts.score_regions(
        samples=[],
        regions="matched.bed",
        motifs=["A,0"],
        contrast=ContrastSpec(
            mode="pairwise",
            numerator=["15min"],
            denominator=["NS"],
            reference_condition="NS",
        ),
        analysis_unit="ensemble_region",
        representation="modified_fraction",
        signal_source="pileup_counts",
        test="beta_binomial",
    )

    assert "p_value" in result.regions.columns
    assert "adjusted_p_value" in result.regions.columns
    assert result.regions["p_value"].between(0, 1).all()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_region_contrasts.py::test_score_regions_beta_binomial_adds_pvalues -q`
Expected: FAIL because p-value columns do not exist

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/region_contrasts.py
from scipy.stats import betabinom


def _estimate_beta_binomial_prior(reference_modified: pd.Series, reference_valid: pd.Series) -> tuple[float, float]:
    fractions = np.divide(
        reference_modified,
        reference_valid,
        out=np.zeros_like(reference_modified, dtype=float),
        where=reference_valid > 0,
    )
    mean = float(fractions.mean())
    var = float(fractions.var(ddof=0))
    if var <= 0:
        concentration = 100.0
    else:
        max_var = max(mean * (1.0 - mean), 1e-9)
        shrink = min(max(var / max_var, 1e-6), 0.99)
        concentration = max((1.0 / shrink) - 1.0, 1.0)
    alpha = max(mean * concentration, 1e-6)
    beta = max((1.0 - mean) * concentration, 1e-6)
    return alpha, beta


def _adjust_pvalues_bh(p_values: pd.Series) -> pd.Series:
    order = np.argsort(p_values.to_numpy())
    ranked = p_values.to_numpy()[order]
    n = len(ranked)
    adjusted = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        running = min(running, ranked[i] * n / rank)
        adjusted[i] = running
    out = np.empty(n, dtype=float)
    out[order] = adjusted
    return pd.Series(out, index=p_values.index)


def _beta_binomial_pvalues(merged: pd.DataFrame) -> pd.Series:
    alpha, beta = _estimate_beta_binomial_prior(
        merged["reference_modified_count"],
        merged["reference_valid_count"],
    )
    pmf = betabinom.pmf(merged["modified_count"], merged["valid_count"], alpha, beta)
    p_values = 1.0 - betabinom.cdf(
        merged["modified_count"] - 1,
        merged["valid_count"],
        alpha,
        beta,
    )
    return pd.Series(np.clip(np.maximum(pmf, p_values), 0.0, 1.0), index=merged.index)
```

Add to `score_regions(...)`:

```python
    if test == "beta_binomial":
        merged["p_value"] = _beta_binomial_pvalues(merged)
        if multiple_testing != "fdr_bh":
            raise ValueError("V1 beta-binomial region_contrasts supports multiple_testing='fdr_bh' only.")
        merged["adjusted_p_value"] = _adjust_pvalues_bh(merged["p_value"])
    else:
        merged["p_value"] = pd.NA
        merged["adjusted_p_value"] = pd.NA
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_region_contrasts.py::test_score_regions_beta_binomial_adds_pvalues -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/region_contrasts.py tests/test_region_contrasts.py
git commit -m "feat: add beta-binomial region contrast scoring"
```

### Task 5: Add User-Facing Docs For Known-Region Contrasts

**Files:**
- Create: `docs/region-contrasts.md`
- Modify: `README.md`

- [ ] **Step 1: Write the docs content**

```markdown
# Region Contrasts

`dimelo.region_contrasts` scores known regions from pileup-backed inputs.

## When To Use It

- use this when you already know the regions you want to compare
- use `region_discovery` later for unknown loci
- use `cluster` first if you want read-state or cluster-occupancy summaries instead of average motif abundance

## V1 Supported Path

- `analysis_unit="ensemble_region"`
- `representation="modified_fraction"` or `"modified_count"`
- `signal_source="pileup_counts"`
- `test="effect_size_only"` or `"beta_binomial"`

## Example

```python
from dimelo import region_contrasts
from dimelo.models import ContrastSpec, SampleSpec

result = region_contrasts.score_regions(
    samples=[
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5", metadata={"pileup_path": "s1.bed.gz"}),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5", metadata={"pileup_path": "s2.bed.gz"}),
    ],
    regions="matched_regions.bed",
    motifs=["A,0"],
    contrast=ContrastSpec(
        mode="pairwise",
        numerator=["15min"],
        denominator=["NS"],
        reference_condition="NS",
    ),
    analysis_unit="ensemble_region",
    representation="modified_fraction",
    signal_source="pileup_counts",
    test="beta_binomial",
)
```

The canonical outputs are `result.regions`, `result.summary`, and `result.plot_data`.
```

- [ ] **Step 2: Save the docs and link them from the README clustering/analysis section**

Run: `sed -n '1,220p' docs/region-contrasts.md`
Expected: shows the new guide with the v1 support matrix and example

- [ ] **Step 3: Commit**

```bash
git add docs/region-contrasts.md README.md
git commit -m "docs: add region contrast workflow guide"
```
