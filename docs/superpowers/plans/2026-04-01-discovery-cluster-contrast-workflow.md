# Discovery Cluster Contrast Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one end-to-end workflow that discovers loci, clusters the selected loci, runs defined-region contrasts on that same selected set, and returns all four outputs together with the full scan table as context.

**Architecture:** Extend the current `discovery_cluster_workflow(...)` path instead of inventing a parallel orchestration stack. Add one wrapper result model for the combined outputs, then add one workflow entry point in `dimelo.workflows` that composes the existing `region_discovery.scan_genome(...)`, `shared_cluster_distribution(...)`, and `region_contrasts.score_regions(...)` functions in that order, using the selected loci as the default contrast region set.

**Tech Stack:** Python dataclasses, pandas, existing dimelo discovery/clustering/contrast modules, pytest

---

### Task 1: Add The Combined Result Model

**Files:**
- Modify: `dimelo/models.py`
- Modify: `dimelo/__init__.py`
- Test: `tests/test_models.py`

- [ ] **Step 1: Write the failing model tests**

```python
def test_region_discovery_cluster_contrast_result_requires_non_none_values():
    with pytest.raises(ValueError, match="discovery"):
        RegionDiscoveryClusterContrastResult(
            discovery=None,
            clustering=_dummy_shared_cluster_result(),
            contrasts=_dummy_region_contrast_result(),
            selected_regions=pd.DataFrame([{"chrom": "chr1", "start": 0, "end": 100}]),
            metadata={},
        )


def test_region_discovery_cluster_contrast_result_validates_wrapped_types():
    with pytest.raises(TypeError, match="RegionDiscoveryResult"):
        RegionDiscoveryClusterContrastResult(
            discovery={},
            clustering=_dummy_shared_cluster_result(),
            contrasts=_dummy_region_contrast_result(),
            selected_regions=pd.DataFrame([{"chrom": "chr1", "start": 0, "end": 100}]),
            metadata={},
        )


def test_region_discovery_cluster_contrast_result_accepts_valid_payloads():
    result = RegionDiscoveryClusterContrastResult(
        discovery=_dummy_region_discovery_result(),
        clustering=_dummy_shared_cluster_result(),
        contrasts=_dummy_region_contrast_result(),
        selected_regions=pd.DataFrame([{"chrom": "chr1", "start": 0, "end": 100}]),
        metadata={"contrast_scope": "selected"},
    )
    assert result.metadata["contrast_scope"] == "selected"
```

- [ ] **Step 2: Run the model tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py -q`
Expected: FAIL with `NameError` / import errors for `RegionDiscoveryClusterContrastResult`

- [ ] **Step 3: Write the minimal model implementation**

```python
@dataclass
class RegionDiscoveryClusterContrastResult:
    discovery: RegionDiscoveryResult
    clustering: SharedClusterResult
    contrasts: RegionContrastResult
    selected_regions: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "discovery": self.discovery,
            "clustering": self.clustering,
            "contrasts": self.contrasts,
            "selected_regions": self.selected_regions,
            "metadata": self.metadata,
            "figures": self.figures,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "RegionDiscoveryClusterContrastResult requires non-None values for: "
                f"{', '.join(missing)}"
            )
        if not isinstance(self.discovery, RegionDiscoveryResult):
            raise TypeError(
                "RegionDiscoveryClusterContrastResult.discovery must be a RegionDiscoveryResult"
            )
        if not isinstance(self.clustering, SharedClusterResult):
            raise TypeError(
                "RegionDiscoveryClusterContrastResult.clustering must be a SharedClusterResult"
            )
        if not isinstance(self.contrasts, RegionContrastResult):
            raise TypeError(
                "RegionDiscoveryClusterContrastResult.contrasts must be a RegionContrastResult"
            )
```

- [ ] **Step 4: Export the new model from the package root**

```python
from .models import RegionDiscoveryClusterContrastResult
...
__all__ = [
    ...
    "RegionDiscoveryClusterContrastResult",
]
```

- [ ] **Step 5: Run the model tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add dimelo/models.py dimelo/__init__.py tests/test_models.py
git commit -m "feat: add discovery cluster contrast result"
```

### Task 2: Add The End-To-End Workflow

**Files:**
- Modify: `dimelo/workflows.py`
- Test: `tests/test_workflows.py`

- [ ] **Step 1: Write the failing workflow tests**

```python
def test_discovery_cluster_contrast_workflow_returns_all_results(monkeypatch):
    monkeypatch.setattr(workflows.region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)
    monkeypatch.setattr(workflows.region_contrasts, "score_regions", _mock_region_contrast_result)

    result = workflows.discovery_cluster_contrast_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1500},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
        contrasts={
            "contrast": ContrastSpec(
                mode="pairwise",
                numerator=["15min"],
                denominator=["NS"],
            ),
            "test": "effect_size_only",
        },
    )

    assert result.discovery.hits.shape[0] == 3
    assert result.clustering.model.mode == "region_anchored"
    assert result.contrasts.metadata["test"] == "effect_size_only"
    assert result.metadata["contrast_scope"] == "selected"
```

```python
def test_discovery_cluster_contrast_workflow_scores_selected_regions_by_default(monkeypatch):
    captured = {}
    ...
    def fake_score_regions(**kwargs):
        captured["regions"] = kwargs["regions"]
        return _mock_region_contrast_result(**kwargs)

    ...
    result = workflows.discovery_cluster_contrast_workflow(...)
    assert captured["regions"] == ["chr1:0-500,+", "chr1:500-1000,-"]
    assert result.selected_regions["name"].tolist() == ["chr1:0-500", "chr1:500-1000"]
```

```python
def test_discovery_cluster_contrast_workflow_rejects_missing_contrast_config(monkeypatch):
    ...
    with pytest.raises(ValueError, match="requires contrasts\\['contrast'\\]"):
        workflows.discovery_cluster_contrast_workflow(...)
```

- [ ] **Step 2: Run the workflow tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_workflows.py -q`
Expected: FAIL with missing `discovery_cluster_contrast_workflow`

- [ ] **Step 3: Implement the orchestration workflow**

```python
def discovery_cluster_contrast_workflow(
    *,
    samples: Iterable[SampleSpec],
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    discovery: dict[str, Any],
    clustering: dict[str, Any],
    contrasts: dict[str, Any],
    selection: dict[str, Any] | None = None,
) -> RegionDiscoveryClusterContrastResult:
    sample_list = list(samples)
    motif_list = list(motifs)
    contrast_config = dict(contrasts)
    if "contrast" not in contrast_config:
        raise ValueError("discovery_cluster_contrast_workflow requires contrasts['contrast'].")

    discovery_cluster_result = discovery_cluster_workflow(
        samples=sample_list,
        motifs=motif_list,
        genome_sizes=genome_sizes,
        discovery=discovery,
        clustering=clustering,
        selection=selection,
    )
    selected_region_spec = _selected_regions_to_region_spec(
        discovery_cluster_result.selected_regions
    )
    contrast_result = region_contrasts.score_regions(
        samples=sample_list,
        regions=selected_region_spec,
        motifs=motif_list,
        **contrast_config,
    )
    return RegionDiscoveryClusterContrastResult(
        discovery=discovery_cluster_result.discovery,
        clustering=discovery_cluster_result.clustering,
        contrasts=contrast_result,
        selected_regions=discovery_cluster_result.selected_regions,
        metadata={
            **discovery_cluster_result.metadata,
            "contrast_scope": "selected",
            "full_scan_windows": discovery_cluster_result.discovery.windows.copy(),
        },
    )
```

- [ ] **Step 4: Add the remaining workflow tests**

```python
def test_discovery_cluster_contrast_workflow_preserves_full_scan_windows_context(monkeypatch):
    ...
    result = workflows.discovery_cluster_contrast_workflow(...)
    assert result.metadata["full_scan_windows"].equals(result.discovery.windows)
```

```python
def test_discovery_cluster_contrast_workflow_fast_fails_invalid_contrast_before_reuse(monkeypatch):
    called = {"scan": False}
    def fake_scan_genome(**kwargs):
        called["scan"] = True
        return _mock_discovery_result()

    monkeypatch.setattr(workflows.region_discovery, "scan_genome", fake_scan_genome)
    with pytest.raises(ValueError, match="analysis_unit='ensemble_region'"):
        workflows.discovery_cluster_contrast_workflow(
            ...,
            contrasts={
                "contrast": ContrastSpec(mode="pairwise", numerator=["15min"], denominator=["NS"]),
                "analysis_unit": "single_read",
            },
        )
    assert called["scan"] is False
```

- [ ] **Step 5: Run the workflow tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_workflows.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add dimelo/workflows.py tests/test_workflows.py
git commit -m "feat: add discovery cluster contrast workflow"
```

### Task 3: Document The New Handoff And Run The Broader Slice

**Files:**
- Modify: `docs/region-contrasts.md`
- Modify: `docs/region-discovery.md`
- Modify: `docs/shared-clustering.md`
- Modify: `README.md`
- Test: `tests/test_workflows.py`

- [ ] **Step 1: Add one end-to-end doc example**

```python
result = workflows.discovery_cluster_contrast_workflow(
    samples=samples,
    motifs=["A,0"],
    genome_sizes={"chr1": 248_956_422},
    discovery={
        "window_size": 2000,
        "step_size": 500,
        "score": "beta_binomial",
        "contrast": discovery_contrast,
    },
    clustering={
        "mode": "region_anchored",
        "n_clusters": 6,
    },
    contrasts={
        "contrast": region_contrast,
        "test": "beta_binomial",
    },
    selection={"mode": "top_n", "top_n": 250},
)

selected = result.selected_regions
clustered = result.clustering.assignments
scored = result.contrasts.regions
windows = result.metadata["full_scan_windows"]
```

- [ ] **Step 2: Clarify the contract in docs**

```markdown
- `result.selected_regions` is the BED-style selected follow-up set
- clustering receives a serializable region-spec derived from those rows
- `result.contrasts` scores the same selected loci by default
- `result.metadata["full_scan_windows"]` carries the full discovery scan for context
```

- [ ] **Step 3: Add/adjust one test if docs exposed a missing behavior**

```python
def test_discovery_cluster_contrast_workflow_metadata_marks_selected_scope(monkeypatch):
    ...
    result = workflows.discovery_cluster_contrast_workflow(...)
    assert result.metadata["contrast_scope"] == "selected"
```

- [ ] **Step 4: Run the broader verification slice**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py tests/test_region_contrasts.py tests/test_region_discovery.py tests/test_workflows.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add docs/region-contrasts.md docs/region-discovery.md docs/shared-clustering.md README.md tests/test_workflows.py
git commit -m "docs: add discovery cluster contrast workflow guide"
```

---

## Self-Review

- Spec coverage: this plan adds the combined result object, the selected-set-first orchestration workflow, validation for contrast config, and docs for the end-to-end discovery -> cluster -> contrast path.
- Placeholder scan: every task includes concrete files, tests, commands, and implementation snippets; there are no TODO placeholders.
- Type consistency: the plan consistently uses `RegionDiscoveryClusterContrastResult`, `discovery_cluster_contrast_workflow`, `selected_regions`, and `full_scan_windows` across tasks.
