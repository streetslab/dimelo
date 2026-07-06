# Region Discovery To Cluster Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add one end-to-end workflow that runs region discovery, selects discovered loci, clusters those loci with the shared clustering engine, and returns both the discovery and clustering results.

**Architecture:** Build a thin orchestration layer in `dimelo.workflows` that composes the existing `dimelo.region_discovery.scan_genome(...)` and `dimelo.workflows.shared_cluster_distribution(...)` entry points. Keep the result data-first by introducing a small wrapper model that stores both outputs plus the exact discovered-region subset passed into clustering.

**Tech Stack:** Python dataclasses, pandas, existing dimelo pileup/extract workflows, pytest

---

### Task 1: Add The Combined Result Model

**Files:**
- Modify: `dimelo/models.py`
- Modify: `dimelo/__init__.py`
- Test: `tests/test_models.py`

- [ ] **Step 1: Write the failing tests**

```python
def test_region_discovery_cluster_result_requires_non_none_values():
    with pytest.raises(ValueError, match="discovery"):
        RegionDiscoveryClusterResult(
            discovery=None,
            clustering=_dummy_shared_cluster_result(),
            selected_regions=pd.DataFrame(),
            metadata={},
        )


def test_region_discovery_cluster_result_accepts_valid_payloads():
    result = RegionDiscoveryClusterResult(
        discovery=_dummy_region_discovery_result(),
        clustering=_dummy_shared_cluster_result(),
        selected_regions=pd.DataFrame([{"chromosome": "chr1", "start": 0, "end": 100}]),
        metadata={"selection_mode": "top_n"},
    )
    assert list(result.selected_regions.columns) == ["chromosome", "start", "end"]
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py -q`
Expected: FAIL with `NameError` or import errors for `RegionDiscoveryClusterResult`

- [ ] **Step 3: Write the minimal model implementation**

```python
@dataclass
class RegionDiscoveryClusterResult:
    discovery: RegionDiscoveryResult
    clustering: SharedClusterResult
    selected_regions: pd.DataFrame
    metadata: dict[str, Any] = field(default_factory=dict)
    figures: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        required_fields = {
            "discovery": self.discovery,
            "clustering": self.clustering,
            "selected_regions": self.selected_regions,
            "metadata": self.metadata,
            "figures": self.figures,
        }
        missing = [name for name, value in required_fields.items() if value is None]
        if missing:
            raise ValueError(
                "RegionDiscoveryClusterResult requires non-None values for: "
                f"{', '.join(missing)}"
            )
```

- [ ] **Step 4: Export the model**

```python
from .models import (
    ...
    RegionDiscoveryClusterResult,
)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add dimelo/models.py dimelo/__init__.py tests/test_models.py
git commit -m "feat: add discovery to cluster workflow result model"
```

### Task 2: Implement The Discovery To Cluster Workflow

**Files:**
- Modify: `dimelo/workflows.py`
- Test: `tests/test_workflows.py`

- [ ] **Step 1: Write the failing workflow tests**

```python
def test_region_discovery_cluster_workflow_returns_both_results(monkeypatch):
    monkeypatch.setattr(region_discovery, "scan_genome", _mock_discovery_result)
    monkeypatch.setattr(workflows, "shared_cluster_distribution", _mock_cluster_result)

    result = workflows.discovery_cluster_workflow(
        samples=_workflow_samples(),
        motifs=["A,0"],
        genome_sizes={"chr1": 1000},
        discovery={"window_size": 500, "step_size": 500},
        clustering={"mode": "region_anchored", "n_clusters": 2},
    )

    assert result.discovery.hits.shape[0] == 2
    assert result.clustering.model.mode == "region_anchored"
    assert list(result.selected_regions.columns) == ["chromosome", "start", "end", "name", "score", "strand"]


def test_region_discovery_cluster_workflow_selects_top_n_hits(monkeypatch):
    ...
    result = workflows.discovery_cluster_workflow(
        ...,
        selection={"mode": "top_n", "top_n": 1},
    )
    assert result.selected_regions["name"].tolist() == ["chr1:0-500"]
```

- [ ] **Step 2: Run the targeted tests to verify they fail**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_workflows.py -q`
Expected: FAIL with missing `discovery_cluster_workflow`

- [ ] **Step 3: Implement selection helpers and the orchestration entry point**

```python
def _select_discovery_hits(
    discovery: RegionDiscoveryResult,
    *,
    selection_mode: str,
    top_n: int | None,
    score_threshold: float | None,
    adjusted_p_value_max: float | None,
) -> pd.DataFrame:
    hits = discovery.hits.copy()
    ...
    return region_discovery.hits_to_bed(selected_hits)


def discovery_cluster_workflow(
    *,
    samples: Iterable[SampleSpec],
    motifs: Iterable[str],
    genome_sizes: dict[str, int],
    discovery: dict[str, Any],
    clustering: dict[str, Any],
    selection: dict[str, Any] | None = None,
) -> RegionDiscoveryClusterResult:
    discovery_result = region_discovery.scan_genome(...)
    selected_regions = _select_discovery_hits(...)
    region_spec = _selected_regions_to_region_spec(selected_regions)
    clustering_result = shared_cluster_distribution(
        samples=samples,
        motifs=motifs,
        matched_regions=region_spec,
        **clustering,
    )
    return RegionDiscoveryClusterResult(
        discovery=discovery_result,
        clustering=clustering_result,
        selected_regions=selected_regions,
        metadata={...},
    )
```

- [ ] **Step 4: Cover error handling in tests**

```python
def test_region_discovery_cluster_workflow_errors_when_no_hits_survive_selection(monkeypatch):
    ...
    with pytest.raises(ValueError, match="No discovery hits remained after selection"):
        workflows.discovery_cluster_workflow(...)


def test_region_discovery_cluster_workflow_rejects_unknown_selection_mode(monkeypatch):
    ...
    with pytest.raises(ValueError, match="Unsupported selection mode"):
        workflows.discovery_cluster_workflow(...)
```

- [ ] **Step 5: Run targeted workflow tests**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_workflows.py -q`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add dimelo/workflows.py tests/test_workflows.py
git commit -m "feat: add discovery to cluster workflow"
```

### Task 3: Document And Verify The Integrated User Flow

**Files:**
- Modify: `docs/region-discovery.md`
- Modify: `docs/shared-clustering.md`
- Modify: `README.md`
- Test: `tests/test_workflows.py`

- [ ] **Step 1: Add one integration test that uses real discovery hit rows as the clustering region input**

```python
def test_region_discovery_cluster_workflow_passes_selected_regions_into_clustering(monkeypatch):
    captured = {}
    ...
    def _fake_cluster(**kwargs):
        captured["matched_regions"] = kwargs["matched_regions"]
        return _mock_cluster_result(**kwargs)

    monkeypatch.setattr(workflows, "shared_cluster_distribution", _fake_cluster)
    result = workflows.discovery_cluster_workflow(...)

    assert captured["matched_regions"] == ["chr1:0-500,+", "chr1:500-1000,-"]
    assert result.selected_regions["name"].tolist() == ["chr1:0-500", "chr1:500-1000"]
```

- [ ] **Step 2: Document the default behavior**

```markdown
from dimelo import workflows

result = workflows.discovery_cluster_workflow(
    samples=samples,
    motifs=["A,0"],
    genome_sizes={"chr1": 248956422},
    discovery={"window_size": 2000, "step_size": 500, "score": "beta_binomial", "contrast": contrast},
    clustering={"mode": "region_anchored", "n_clusters": 6},
    selection={"mode": "top_n", "top_n": 250},
)
```

Explain:
- `discovery` config is passed to `scan_genome(...)`
- `selection` defaults to top-ranked hits
- `clustering` config is passed to `shared_cluster_distribution(...)`
- the workflow returns the BED-style `selected_regions` table, but passes a serializable region-spec derived from that table into `matched_regions`
- the result keeps `discovery`, `selected_regions`, and `clustering` together

- [ ] **Step 3: Run the full verification slice**

Run: `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_models.py tests/test_region_discovery.py tests/test_workflows.py -q`
Expected: PASS

- [ ] **Step 4: Commit**

```bash
git add docs/region-discovery.md docs/shared-clustering.md README.md tests/test_workflows.py
git commit -m "docs: add discovery to cluster workflow guide"
```

---

## Self-Review

- Spec coverage: this plan adds the combined return type, the orchestration workflow, selection policy, error handling, and the user-facing docs for the one-call discovery-to-cluster path.
- Placeholder scan: every task names exact files, targeted tests, and concrete implementation points; there are no TODO placeholders.
- Type consistency: `RegionDiscoveryClusterResult`, `discovery_cluster_workflow`, and `selected_regions` are used consistently across tasks.
