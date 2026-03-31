# Shared Clustering Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a stable, backward-compatible foundation for shared clustering across datasets, including result models, artifact resolution, data-first plotting payloads, and a usable `shared_cluster_distribution()` workflow covering `read_global` and `region_anchored`.

**Architecture:** Build new orchestration and data-model modules around the existing `parse_bam`, `load_processed`, and `cluster` code without changing their semantics. Keep canonical outputs as DataFrames plus `plot_data`, with Matplotlib remaining a thin compatibility renderer instead of the core contract.

**Tech Stack:** Python 3.11, pandas, numpy, scikit-learn, matplotlib, pytest

---

Scope note:

- This plan intentionally covers the first executable slice only:
  shared clustering foundations and orchestration.
- Follow-on plans should handle:
  `region_contrasts`,
  `region_discovery`,
  and `global_analysis`.

### Task 1: Add Shared Clustering Models And Public Exports

**Files:**
- Create: `dimelo/models.py`
- Create: `tests/test_models.py`
- Modify: `dimelo/__init__.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_models.py
from pathlib import Path

from dimelo.models import (
    CohortSpec,
    DatasetArtifact,
    SampleSpec,
    SharedClusterModel,
    SharedClusterResult,
)


def test_sample_spec_basic_fields():
    sample = SampleSpec(
        sample_id="s1",
        condition="NS",
        extract_h5=Path("reads.h5"),
        regions_bed=Path("regions.bed"),
        replicate=1,
    )
    assert sample.sample_id == "s1"
    assert sample.condition == "NS"
    assert sample.extract_h5 == Path("reads.h5")


def test_dataset_artifact_carries_core_metadata():
    artifact = DatasetArtifact(
        sample_id="s1",
        artifact_type="read_features",
        path=Path("features.parquet"),
        format="parquet",
        params={"window_size": 1000},
        provenance={"schema_version": 1, "package_version": "1.0.0"},
    )
    assert artifact.artifact_type == "read_features"
    assert artifact.provenance["schema_version"] == 1


def test_shared_cluster_result_holds_plot_data():
    model = SharedClusterModel(
        mode="read_global",
        motifs=["A,0"],
        feature_names=["feat1"],
        preprocessing={},
        estimator=None,
        cluster_labels=["C0"],
        fit_metadata={},
    )
    result = SharedClusterResult(
        model=model,
        assignments=None,
        cluster_distribution=None,
        condition_distribution=None,
        distribution_change=None,
        cluster_profiles=None,
        region_summaries=None,
        plot_data={"cluster_distribution_bar": {"kind": "bar"}},
        figures={},
        metadata={},
    )
    assert "cluster_distribution_bar" in result.plot_data


def test_cohort_spec_tracks_workflow_and_params():
    cohort = CohortSpec(
        cohort_id="cohort-1",
        sample_ids=["s1", "s2"],
        workflow="shared_cluster_distribution",
        params={"mode": "read_global"},
    )
    assert cohort.workflow == "shared_cluster_distribution"
    assert cohort.params["mode"] == "read_global"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_models.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'dimelo.models'`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/models.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class SampleSpec:
    sample_id: str
    condition: str
    extract_h5: str | Path
    regions_bed: str | Path | None = None
    replicate: int | None = None
    metadata: dict[str, Any] | None = None


@dataclass
class DatasetArtifact:
    sample_id: str
    artifact_type: str
    path: str | Path
    format: str
    params: dict[str, Any]
    provenance: dict[str, Any]


@dataclass
class SharedClusterModel:
    mode: str
    motifs: list[str]
    feature_names: list[str]
    preprocessing: dict[str, Any]
    estimator: Any
    cluster_labels: list[str]
    fit_metadata: dict[str, Any]


@dataclass
class SharedClusterResult:
    model: SharedClusterModel
    assignments: pd.DataFrame | None
    cluster_distribution: pd.DataFrame | None
    condition_distribution: pd.DataFrame | None
    distribution_change: pd.DataFrame | None
    cluster_profiles: pd.DataFrame | None
    region_summaries: pd.DataFrame | None
    plot_data: dict[str, pd.DataFrame | dict[str, Any]]
    figures: dict[str, Any] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class CohortSpec:
    cohort_id: str
    sample_ids: list[str]
    workflow: str
    params: dict[str, Any]
    metadata: dict[str, Any] | None = None


@dataclass
class BatchJob:
    job_id: str
    workflow: str
    cohorts: list[CohortSpec]
    artifact_policy: str = "prefer_cached"
    metadata: dict[str, Any] | None = None
```

```python
# dimelo/__init__.py
from . import (
    cluster,
    export,
    load_processed,
    models,
    parse_bam,
    plot_depth_histogram,
    plot_depth_profile,
    plot_enrichment,
    plot_enrichment_profile,
    plot_read_browser,
    plot_reads,
)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_models.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/models.py dimelo/__init__.py tests/test_models.py
git commit -m "feat: add shared clustering model types"
```

### Task 2: Add Artifact Compatibility And Resolver Foundations

**Files:**
- Create: `dimelo/artifacts.py`
- Create: `tests/test_artifacts.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_artifacts.py
from pathlib import Path

from dimelo.artifacts import artifact_fingerprint, artifact_is_compatible, resolve_artifact
from dimelo.models import DatasetArtifact


def test_artifact_fingerprint_includes_required_keys():
    fp = artifact_fingerprint(
        schema_version=1,
        package_version="1.0.0",
        source_files=[Path("reads.h5")],
        params={"window_size": 1000, "motifs": ["A,0"]},
    )
    assert fp["schema_version"] == 1
    assert "params_hash" in fp
    assert "source_files" in fp


def test_artifact_is_compatible_rejects_param_mismatch():
    artifact = DatasetArtifact(
        sample_id="s1",
        artifact_type="read_features",
        path=Path("features.parquet"),
        format="parquet",
        params={"window_size": 1000, "motifs": ["A,0"]},
        provenance={"schema_version": 1, "package_version": "1.0.0"},
    )
    assert not artifact_is_compatible(
        artifact,
        required_params={"window_size": 2000, "motifs": ["A,0"]},
        required_provenance={"schema_version": 1},
    )


def test_resolve_artifact_prefers_matching_artifact():
    artifact = DatasetArtifact(
        sample_id="s1",
        artifact_type="read_features",
        path=Path("features.parquet"),
        format="parquet",
        params={"window_size": 1000},
        provenance={"schema_version": 1, "package_version": "1.0.0"},
    )
    resolved = resolve_artifact(
        artifacts=[artifact],
        artifact_type="read_features",
        sample_id="s1",
        required_params={"window_size": 1000},
        required_provenance={"schema_version": 1},
        policy="prefer_cached",
    )
    assert resolved is artifact
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_artifacts.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'dimelo.artifacts'`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/artifacts.py
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from .models import DatasetArtifact


def _stable_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def artifact_fingerprint(
    *,
    schema_version: int,
    package_version: str,
    source_files: Iterable[str | Path],
    params: dict[str, Any],
) -> dict[str, Any]:
    source_list = [str(Path(p)) for p in source_files]
    params_hash = hashlib.sha256(_stable_json(params).encode()).hexdigest()
    return {
        "schema_version": schema_version,
        "package_version": package_version,
        "source_files": source_list,
        "params_hash": params_hash,
    }


def artifact_is_compatible(
    artifact: DatasetArtifact,
    *,
    required_params: dict[str, Any],
    required_provenance: dict[str, Any],
) -> bool:
    for key, value in required_params.items():
        if artifact.params.get(key) != value:
            return False
    for key, value in required_provenance.items():
        if artifact.provenance.get(key) != value:
            return False
    return True


def resolve_artifact(
    *,
    artifacts: list[DatasetArtifact],
    artifact_type: str,
    sample_id: str,
    required_params: dict[str, Any],
    required_provenance: dict[str, Any],
    policy: str = "prefer_cached",
) -> DatasetArtifact | None:
    matches = [
        artifact
        for artifact in artifacts
        if artifact.sample_id == sample_id and artifact.artifact_type == artifact_type
    ]
    if policy == "rebuild":
        return None
    compatible = [
        artifact
        for artifact in matches
        if artifact_is_compatible(
            artifact,
            required_params=required_params,
            required_provenance=required_provenance,
        )
    ]
    if compatible:
        return compatible[0]
    if policy == "require_cached":
        raise FileNotFoundError(
            f"No compatible cached artifact found for sample={sample_id}, type={artifact_type}"
        )
    return None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_artifacts.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/artifacts.py tests/test_artifacts.py
git commit -m "feat: add artifact compatibility helpers"
```

### Task 3: Add Distribution Summaries And Data-First Plot Payloads

**Files:**
- Create: `dimelo/distribution.py`
- Create: `dimelo/plotting.py`
- Create: `tests/test_distribution.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_distribution.py
import pandas as pd

from dimelo.distribution import (
    build_cluster_distribution,
    build_condition_distribution,
    build_distribution_change,
)
from dimelo.plotting import (
    prepare_cluster_distribution_bar_data,
    prepare_cluster_distribution_heatmap_data,
)


def test_build_cluster_distribution_counts_and_fractions():
    assignments = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s1", "s2"],
            "condition": ["NS", "NS", "NS", "15min"],
            "cluster": ["C0", "C0", "C1", "C1"],
        }
    )
    result = build_cluster_distribution(assignments)
    s1 = result[result["sample_id"] == "s1"].sort_values("cluster").reset_index(drop=True)
    assert list(s1["count"]) == [2, 1]
    assert list(s1["fraction"]) == [2 / 3, 1 / 3]


def test_build_distribution_change_against_reference():
    condition_distribution = pd.DataFrame(
        {
            "condition": ["NS", "NS", "15min", "15min"],
            "cluster": ["C0", "C1", "C0", "C1"],
            "fraction": [0.75, 0.25, 0.25, 0.75],
        }
    )
    result = build_distribution_change(condition_distribution, reference_condition="NS")
    row = result[(result["condition"] == "15min") & (result["cluster"] == "C1")].iloc[0]
    assert row["delta_fraction"] == 0.5


def test_plot_data_helpers_return_dataframes():
    cluster_distribution = pd.DataFrame(
        {
            "sample_id": ["s1", "s1"],
            "condition": ["NS", "NS"],
            "cluster": ["C0", "C1"],
            "count": [2, 1],
            "fraction": [2 / 3, 1 / 3],
        }
    )
    bar_data = prepare_cluster_distribution_bar_data(cluster_distribution)
    heatmap_data = prepare_cluster_distribution_heatmap_data(cluster_distribution)
    assert not bar_data.empty
    assert not heatmap_data.empty
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_distribution.py -q`
Expected: FAIL with `ModuleNotFoundError` for `dimelo.distribution`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/distribution.py
from __future__ import annotations

import numpy as np
import pandas as pd


def build_cluster_distribution(assignments: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        assignments.groupby(["sample_id", "condition", "cluster"])
        .size()
        .reset_index(name="count")
    )
    totals = grouped.groupby("sample_id")["count"].transform("sum")
    grouped["fraction"] = grouped["count"] / totals
    return grouped


def build_condition_distribution(cluster_distribution: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        cluster_distribution.groupby(["condition", "cluster"])["count"]
        .sum()
        .reset_index()
    )
    totals = grouped.groupby("condition")["count"].transform("sum")
    grouped["fraction"] = grouped["count"] / totals
    grouped["replicate_n"] = (
        cluster_distribution.groupby("condition")["sample_id"].nunique().reindex(grouped["condition"]).to_numpy()
    )
    return grouped


def build_distribution_change(
    condition_distribution: pd.DataFrame,
    *,
    reference_condition: str,
    pseudocount: float = 1e-6,
) -> pd.DataFrame:
    ref = (
        condition_distribution[condition_distribution["condition"] == reference_condition]
        .loc[:, ["cluster", "fraction"]]
        .rename(columns={"fraction": "reference_fraction"})
    )
    merged = condition_distribution.merge(ref, on="cluster", how="left")
    merged = merged[merged["condition"] != reference_condition].copy()
    merged["delta_fraction"] = merged["fraction"] - merged["reference_fraction"]
    merged["log2_fc"] = np.log2(
        (merged["fraction"] + pseudocount) / (merged["reference_fraction"] + pseudocount)
    )
    return merged
```

```python
# dimelo/plotting.py
from __future__ import annotations

import pandas as pd


def prepare_cluster_distribution_bar_data(cluster_distribution: pd.DataFrame) -> pd.DataFrame:
    return cluster_distribution.copy()


def prepare_cluster_distribution_heatmap_data(cluster_distribution: pd.DataFrame) -> pd.DataFrame:
    return cluster_distribution.pivot_table(
        index="condition",
        columns="cluster",
        values="fraction",
        fill_value=0.0,
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_distribution.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/distribution.py dimelo/plotting.py tests/test_distribution.py
git commit -m "feat: add shared cluster summary and plot payload helpers"
```

### Task 4: Implement `shared_cluster_distribution()` For `read_global`

**Files:**
- Create: `dimelo/workflows.py`
- Modify: `dimelo/cluster.py`
- Modify: `dimelo/__init__.py`
- Create: `tests/test_workflows.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_workflows.py
import numpy as np
import pandas as pd

from dimelo.models import SampleSpec
from dimelo import workflows


def test_shared_cluster_distribution_read_global(monkeypatch):
    fake_samples = [
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5"),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5"),
    ]

    def fake_extract(*args, **kwargs):
        class R:
            data_matrix = np.array([[0, 1, 0, 1], [1, 0, 1, 0]], dtype=float)
            val_matrix = np.ones((2, 4), dtype=float)
            metadata = [
                {"read_name": "r1", "chromosome": "chr1", "region_start": 0, "region_end": 4},
                {"read_name": "r2", "chromosome": "chr1", "region_start": 10, "region_end": 14},
            ]
            datasets = []
            regions_dict = None

        return R()

    def fake_features(result, **kwargs):
        return result.data_matrix, ["f0", "f1", "f2", "f3"]

    monkeypatch.setattr(workflows.cluster, "extract_read_windows", fake_extract)
    monkeypatch.setattr(workflows.cluster, "read_window_feature_matrix", fake_features)

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="read_global",
        motifs=["A,0"],
        n_clusters=2,
        training_sample_per_dataset=2,
        make_plots=False,
    )

    assert not result.assignments.empty
    assert set(result.assignments["sample_id"]) == {"s1", "s2"}
    assert "cluster_distribution_bar" in result.plot_data
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_workflows.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'dimelo.workflows'`

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/workflows.py
from __future__ import annotations

import numpy as np
import pandas as pd

from . import cluster, distribution, plotting
from .models import SharedClusterModel, SharedClusterResult


def shared_cluster_distribution(
    *,
    samples,
    mode: str,
    motifs,
    matched_regions=None,
    signal_normalization: str = "none",
    feature_scaling: str = "robust_zscore",
    cluster_basis: str = "shape_plus_level",
    clusterer: str = "minibatch_kmeans",
    n_clusters: int = 8,
    training_sample_per_dataset: int = 100_000,
    artifact_policy: str = "prefer_cached",
    random_state: int = 42,
    make_plots: bool = True,
) -> SharedClusterResult:
    if mode != "read_global":
        raise NotImplementedError("First slice implements read_global only.")

    assignments = []
    feature_names = None
    matrices = []
    sample_meta = []

    for sample in samples:
        extracted = cluster.extract_read_windows(
            hdf5_file=sample.extract_h5,
            motifs=motifs,
        )
        feature_matrix, feature_names = cluster.read_window_feature_matrix(extracted)
        matrices.append(feature_matrix)
        sample_meta.extend(
            [
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    **meta,
                }
                for meta in extracted.metadata
            ]
        )

    all_features = pd.DataFrame(
        data=np.vstack(matrices),
        columns=feature_names,
    )
    result = cluster.cluster_read_windows(
        all_features.to_numpy(),
        method=clusterer,
        n_clusters=n_clusters,
        random_state=random_state,
    )

    for meta, label in zip(sample_meta, result.labels_size_ordered):
        assignments.append({**meta, "cluster": f"C{int(label)}"})

    assignments_df = pd.DataFrame(assignments)
    cluster_distribution = distribution.build_cluster_distribution(assignments_df)
    condition_distribution = distribution.build_condition_distribution(cluster_distribution)
    plot_data = {
        "cluster_distribution_bar": plotting.prepare_cluster_distribution_bar_data(
            cluster_distribution
        ),
        "cluster_distribution_heatmap": plotting.prepare_cluster_distribution_heatmap_data(
            cluster_distribution
        ),
    }
    model = SharedClusterModel(
        mode=mode,
        motifs=list(motifs),
        feature_names=list(feature_names or []),
        preprocessing={
            "signal_normalization": signal_normalization,
            "feature_scaling": feature_scaling,
            "cluster_basis": cluster_basis,
        },
        estimator=result.model,
        cluster_labels=sorted(assignments_df["cluster"].unique()),
        fit_metadata={
            "clusterer": clusterer,
            "n_clusters": n_clusters,
            "artifact_policy": artifact_policy,
            "training_sample_per_dataset": training_sample_per_dataset,
        },
    )
    return SharedClusterResult(
        model=model,
        assignments=assignments_df,
        cluster_distribution=cluster_distribution,
        condition_distribution=condition_distribution,
        distribution_change=None,
        cluster_profiles=None,
        region_summaries=None,
        plot_data=plot_data,
        figures={},
        metadata={"mode": mode},
    )
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
    workflows,
    plot_depth_histogram,
    plot_depth_profile,
    plot_enrichment,
    plot_enrichment_profile,
    plot_read_browser,
    plot_reads,
)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_workflows.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/workflows.py tests/test_workflows.py
git commit -m "feat: add read-global shared clustering workflow"
```

### Task 5: Add `region_anchored` Support And Backward-Compatibility Regression Coverage

**Files:**
- Create: `dimelo/region_analysis.py`
- Modify: `dimelo/workflows.py`
- Modify: `tests/test_workflows.py`
- Modify: `tests/test_cluster.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_workflows.py
import numpy as np

from dimelo.models import SampleSpec
from dimelo import workflows


def test_shared_cluster_distribution_region_anchored(monkeypatch):
    fake_samples = [
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5", regions_bed="r1.bed"),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5", regions_bed="r2.bed"),
    ]

    def fake_region_table(*args, **kwargs):
        return np.array([[0.2, 0.8], [0.7, 0.3]]), [
            {"region_id": "reg1", "sample_id": "s1", "condition": "NS"},
            {"region_id": "reg1", "sample_id": "s2", "condition": "15min"},
        ]

    monkeypatch.setattr(workflows.region_analysis, "build_region_feature_table", fake_region_table)

    result = workflows.shared_cluster_distribution(
        samples=fake_samples,
        mode="region_anchored",
        motifs=["A,0"],
        matched_regions="matched.bed",
        n_clusters=2,
        make_plots=False,
    )
    assert not result.assignments.empty
    assert "region_id" in result.assignments.columns
```

```python
# tests/test_cluster.py
def test_existing_cluster_result_api_still_returns_metrics():
    rng = np.random.default_rng(0)
    feature_matrix = rng.normal(size=(10, 4))
    result = cluster.cluster_read_windows(feature_matrix, method="kmeans", n_clusters=2, random_state=0)
    assert "silhouette" in result.metrics
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_workflows.py tests/test_cluster.py -q`
Expected: FAIL with `AttributeError` for missing `build_region_feature_table` or missing `region_anchored` branch

- [ ] **Step 3: Write minimal implementation**

```python
# dimelo/region_analysis.py
from __future__ import annotations

from typing import Any, Sequence

import numpy as np
import pandas as pd

from . import cluster


def build_region_feature_table(
    *,
    samples,
    motifs: Sequence[str],
    matched_regions,
):
    matrices = []
    metadata = []
    for sample in samples:
        matrix, region_meta = cluster.region_feature_matrix_from_pileup(
            bedmethyl_file=sample.metadata["pileup_path"],
            motif=motifs[0],
            regions=matched_regions or sample.regions_bed,
        )
        matrices.append(matrix)
        metadata.extend(
            [
                {
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "replicate": sample.replicate,
                    "region_id": f"{chrom}:{start}-{end}",
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                }
                for chrom, start, end, strand in region_meta
            ]
        )
    return np.vstack(matrices), metadata
```

```python
# dimelo/workflows.py
    if mode == "region_anchored":
        feature_matrix, metadata_rows = region_analysis.build_region_feature_table(
            samples=samples,
            motifs=motifs,
            matched_regions=matched_regions,
        )
        result = cluster.cluster_read_windows(
            feature_matrix,
            method=clusterer,
            n_clusters=n_clusters,
            random_state=random_state,
        )
        assignments_df = pd.DataFrame(metadata_rows)
        assignments_df["cluster"] = [f"C{int(x)}" for x in result.labels_size_ordered]
        cluster_distribution = distribution.build_cluster_distribution(assignments_df)
        condition_distribution = distribution.build_condition_distribution(cluster_distribution)
        return SharedClusterResult(
            model=SharedClusterModel(
                mode=mode,
                motifs=list(motifs),
                feature_names=[f"pos_{i}" for i in range(feature_matrix.shape[1])],
                preprocessing={
                    "signal_normalization": signal_normalization,
                    "feature_scaling": feature_scaling,
                    "cluster_basis": cluster_basis,
                },
                estimator=result.model,
                cluster_labels=sorted(assignments_df["cluster"].unique()),
                fit_metadata={"clusterer": clusterer, "n_clusters": n_clusters},
            ),
            assignments=assignments_df,
            cluster_distribution=cluster_distribution,
            condition_distribution=condition_distribution,
            distribution_change=None,
            cluster_profiles=None,
            region_summaries=assignments_df,
            plot_data={
                "cluster_distribution_bar": plotting.prepare_cluster_distribution_bar_data(cluster_distribution),
                "cluster_distribution_heatmap": plotting.prepare_cluster_distribution_heatmap_data(cluster_distribution),
            },
            figures={},
            metadata={"mode": mode},
        )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_workflows.py tests/test_cluster.py -q`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add dimelo/region_analysis.py dimelo/workflows.py tests/test_workflows.py tests/test_cluster.py
git commit -m "feat: add region-anchored shared clustering workflow"
```

### Task 6: Add Usage Docs For Backward-Compatible Workflow Entry

**Files:**
- Create: `docs/shared-clustering.md`
- Modify: `README.md`

- [ ] **Step 1: Write the failing docs check**

```markdown
<!-- docs/shared-clustering.md -->
# Shared Clustering

This file should explain:
- when to run `parse_bam.extract()`
- when to run `parse_bam.pileup()`
- when to run both
- how to call `shared_cluster_distribution()`
- how to plot from returned tables without package renderers
```

- [ ] **Step 2: Verify the docs do not yet exist**

Run: `test -f docs/shared-clustering.md`
Expected: non-zero exit status

- [ ] **Step 3: Write the docs**

```markdown
# Shared Clustering

## Preprocessing Choices

- Run `parse_bam.extract()` when you want read-level clustering or single-read pattern analysis.
- Run `parse_bam.pileup()` when you want region-anchored clustering from pileup-derived summaries.
- Run both when you want abundance-style region summaries plus deeper read-level follow-up.

## Workflow Entry Point

```python
from dimelo import workflows
from dimelo.models import SampleSpec

result = workflows.shared_cluster_distribution(
    samples=[
        SampleSpec(sample_id="s1", condition="NS", extract_h5="s1.h5"),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="s2.h5"),
    ],
    mode="read_global",
    motifs=["A,0"],
    n_clusters=8,
)
```

## Custom Plotting

The canonical outputs are:

- `result.assignments`
- `result.cluster_distribution`
- `result.condition_distribution`
- `result.plot_data`

You can use these directly in Matplotlib, seaborn, Plotly, Altair, or your own plotting stack.
```
```

- [ ] **Step 4: Review the docs render and match the implemented API**

Run: `sed -n '1,220p' docs/shared-clustering.md`
Expected: shows preprocessing guidance, workflow entry point, and custom plotting guidance

- [ ] **Step 5: Commit**

```bash
git add docs/shared-clustering.md README.md
git commit -m "docs: add shared clustering usage guide"
```
