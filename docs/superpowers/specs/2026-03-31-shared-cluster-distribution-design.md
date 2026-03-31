# Shared Cluster Distribution Design

## Summary

Add a package-level workflow that defines a single set of clustering boundaries across multiple matched datasets, applies those boundaries consistently to every dataset, and reports how cluster occupancy changes across conditions. The workflow must support two related use cases:

1. `read_global`: cluster reads pooled across datasets independent of locus, then project cluster assignments back onto matched regions.
2. `region_anchored`: build matched region-level feature rows across datasets, cluster those rows once, and compare condition-level shifts in region-cluster occupancy.

The workflow should scale from thousands to millions of reads or regions per dataset, with roughly a dozen datasets per run.

## Goals

- Define cluster boundaries once and reuse them across all datasets in a run.
- Make the shared-boundary workflow available as a stable package API rather than notebook-only code.
- Preserve compatibility with existing `dimelo` extract outputs and existing feature-engineering code in `dimelo.cluster`.
- Support matched-region analyses while also allowing read clustering that is independent of region.
- Return structured tabular results suitable for downstream notebook analysis, plus a small standard plot set for quick interpretation.
- Use a performant default that remains usable at million-row scale.
- Introduce a hybrid batch model that supports both independent dataset processing and flexible cohort-level shared workflows.
- Make downstream workflows consume reusable artifacts when available, while still being able to rebuild from raw inputs when necessary.
- Preserve continuity with prior package versions by keeping the existing low-level parsing layer stable and building new workflows on top of it.

## Non-Goals

- Reproduce every ad hoc plot or helper from the collaborator notebook as package API.
- Ship experiment-specific hard-coded condition schemes as core logic.
- Recluster separately for each condition or each region.
- Store large raw matrices on result objects by default.
- Rewrite or substantially alter the behavior of the existing parsing internals unless behavior-preserving refactors are required.

## User-Facing Workflow

### Primary Entry Point

```python
from dimelo.models import SampleSpec
from dimelo import workflows

samples = [
    SampleSpec(
        sample_id="sample_1",
        condition="NS",
        extract_h5="path/to/sample_1.h5",
        regions_bed="path/to/sample_1_regions.bed",
        replicate=1,
    ),
    SampleSpec(
        sample_id="sample_2",
        condition="15min",
        extract_h5="path/to/sample_2.h5",
        regions_bed="path/to/sample_2_regions.bed",
        replicate=1,
    ),
]

result = workflows.shared_cluster_distribution(
    samples=samples,
    mode="read_global",
    motifs=["A,0"],
    matched_regions=None,
    signal_normalization="none",
    feature_scaling="robust_zscore",
    cluster_basis="shape_plus_level",
    clusterer="minibatch_kmeans",
    n_clusters=8,
    training_sample_per_dataset=100_000,
    random_state=42,
    make_plots=True,
)
```

### Supported Modes

#### `read_global`

- Extract read-level windows and read-level feature matrices from every dataset.
- Build one balanced pooled training set across datasets.
- Fit one shared clustering model.
- Assign cluster labels to every read in every dataset.
- Optionally map those read labels back to matched regions for per-region occupancy summaries.

#### `region_anchored`

- Build one feature row per matched region per dataset from aggregated read-level metrics.
- Fit one shared clustering model across the pooled region rows.
- Assign cluster labels to every matched region row.
- Summarize cluster distribution shifts across conditions.

## Recommended Default Algorithm

### Default Strategy

Use a balanced pooled reference with `MiniBatchKMeans`.

1. Extract features separately per dataset.
2. Sample up to `training_sample_per_dataset` rows from each dataset.
3. Concatenate those samples into a pooled training matrix.
4. Fit one preprocessing pipeline on the pooled matrix.
5. Fit one shared clustering model on the pooled matrix.
6. Apply the frozen preprocessing and clustering model to the full datasets in chunks.
7. Summarize cluster distributions globally and by matched region.

The preprocessing pipeline in this algorithm must distinguish between:

- biological signal normalization before feature construction
- feature scaling after feature construction

### Why This Default

- More fair than anchoring cluster boundaries on one reference dataset.
- More performant than full pooled `KMeans` when rows scale into the millions.
- Simpler and more interpretable than per-region or per-condition reclustering.
- Produces fixed centroids, which makes prediction for full datasets cheap.

### Alternatives Kept Open

- `KMeans` for smaller runs.
- `GaussianMixture`, `Birch`, or `HDBSCAN` if later needed, but not as the initial default.
- Reference-dataset anchoring as an explicit option rather than the default path.

## Architecture

### Module Layout

#### `dimelo/models.py`

Owns the workflow dataclasses:

- `SampleSpec`
- `SharedClusterModel`
- `SharedClusterResult`

#### `dimelo/distribution.py`

Owns generic summary logic:

- counts per dataset and condition
- cluster fractions
- baseline deltas
- log2 fold changes
- replicate-aware aggregations

#### `dimelo/region_analysis.py`

Owns region-aware helpers:

- matched-region validation
- read-to-region projection
- region-by-condition occupancy summaries
- matched region feature table construction for `region_anchored`

#### `dimelo/workflows.py`

Owns orchestration:

- input validation
- dataset iteration
- balanced pooled training subset construction
- fit-once / assign-many execution
- standard plot creation
- result packaging
- artifact resolution for workflow inputs

#### `dimelo/artifacts.py`

Owns reusable dataset-level derived outputs and lookup helpers:

- artifact key definitions
- artifact metadata and provenance
- cache hit / miss resolution
- compatibility checks between requested parameters and stored artifacts

#### `dimelo/batch.py`

Owns hybrid batch orchestration:

- independent dataset processing jobs
- named cohort jobs
- workflow recipes
- batch manifests for downstream workflow execution

#### `dimelo/cluster.py`

Retains low-level primitives already present on the branch:

- `extract_read_windows`
- `read_window_feature_matrix`
- `cluster_read_windows`
- generic cluster profile plotting
- generic region-cluster plotting and summaries

The new workflow should call these lower-level functions rather than duplicate them.

## Compatibility And Continuity Requirements

The shared clustering work must be additive around the existing parsing architecture.

Required constraints:

- keep `parse_bam.pileup()` and `parse_bam.extract()` as the core preprocessing entry points
- keep existing modkit-backed parsing logic intact unless a change is behavior-preserving and covered by regression tests
- avoid breaking existing output concepts and file semantics
- prefer additive metadata or manifest layers over changing current parsed output formats
- ensure new workflow wrappers call into existing parsing outputs rather than replacing them

User-facing consequence:

- previous preprocessing habits should still work
- existing pileup or extract outputs should remain usable in the updated package
- users should learn new downstream workflows without needing to relearn the parsing layer

## Data Model

### `SampleSpec`

```python
@dataclass
class SampleSpec:
    sample_id: str
    condition: str
    extract_h5: str | Path
    regions_bed: str | Path | None = None
    replicate: int | None = None
    metadata: dict[str, Any] | None = None
```

### `SharedClusterModel`

```python
@dataclass
class SharedClusterModel:
    mode: str
    motifs: list[str]
    feature_names: list[str]
    preprocessing: dict[str, Any]
    estimator: Any
    cluster_labels: list[str]
    fit_metadata: dict[str, Any]
```

### `SharedClusterResult`

```python
@dataclass
class SharedClusterResult:
    model: SharedClusterModel
    assignments: pd.DataFrame
    cluster_distribution: pd.DataFrame
    condition_distribution: pd.DataFrame
    distribution_change: pd.DataFrame | None
    cluster_profiles: pd.DataFrame
    region_summaries: pd.DataFrame | None
    figures: dict[str, Any]
    metadata: dict[str, Any]
```

### `DatasetArtifact`

```python
@dataclass
class DatasetArtifact:
    sample_id: str
    artifact_type: str
    path: str | Path
    format: str
    params: dict[str, Any]
    provenance: dict[str, Any]
```

### `CohortSpec`

```python
@dataclass
class CohortSpec:
    cohort_id: str
    sample_ids: list[str]
    workflow: str
    params: dict[str, Any]
    metadata: dict[str, Any] | None = None
```

### `BatchJob`

```python
@dataclass
class BatchJob:
    job_id: str
    workflow: str
    cohorts: list[CohortSpec]
    artifact_policy: str = "prefer_cached"
    metadata: dict[str, Any] | None = None
```

## Hybrid Batch Model

The package should support two complementary batch concepts rather than forcing one batch unit everywhere.

### Dataset-Centric Processing

Each dataset can be processed independently into a standard artifact bundle. This is the high-throughput, cache-friendly layer.

Examples of dataset-level artifacts:

- QC summaries
- read-level feature tables
- region-mapping tables
- cluster-ready metadata tables
- cached transformed feature matrices

This layer is best for:

- parallel processing of many unrelated datasets
- retry and resume at per-dataset granularity
- feeding many downstream workflows from the same reusable intermediates

### Cohort-Centric Workflows

Named cohorts should define which datasets participate in a shared analysis recipe.

Examples:

- shared-boundary clustering over six matched samples
- region occupancy analysis over one subset
- a different recipe over an unrelated subset processed in parallel

This layer is best for:

- pooled fitting
- consistent cluster boundary definition
- any workflow where the group of samples matters scientifically

### Recommended Hybrid Structure

The system should use both layers:

1. dataset processing produces reusable standardized artifacts
2. cohort workflows consume those artifacts to run shared analyses

This preserves flexibility while avoiding repeated raw-input processing.

## Artifact Resolution Policy

Downstream workflows should prefer cached artifacts but fall back to raw inputs automatically when required artifacts are missing or incompatible.

### Default Policy

Use `artifact_policy="prefer_cached"` as the default.

Resolution order:

1. look for a compatible artifact
2. if found, load it
3. if missing or incompatible, rebuild from raw inputs
4. record whether the workflow used cache hits or rebuilt inputs

### Additional Policies

- `prefer_cached`
  default pragmatic mode
- `require_cached`
  fail if required artifacts are unavailable
- `rebuild`
  ignore cached artifacts and regenerate needed intermediates

### Provenance Requirements

Every workflow result should record:

- which artifacts were loaded from cache
- which artifacts were rebuilt
- the parameters used to determine compatibility
- the raw input paths used when fallback occurred

This keeps the fallback behavior auditable rather than implicit.

## Shared Workflow Input Contract

Shared comparison workflows should consume artifact manifests or dataset specs, not notebook-era ad hoc dictionaries.

### Accepted Input Forms

1. `SampleSpec` objects with raw input paths
2. artifact manifests derived from those samples
3. `CohortSpec` objects that name a workflow and a set of samples

The shared workflow should resolve the required inputs internally through the artifact layer.

### Why This Contract Matters

- the same sample can participate in multiple shared workflows
- expensive preprocessing is reused
- workflows remain runnable even when some artifacts are absent
- batch orchestration stays decoupled from notebook-specific path plumbing

## Detailed Workflow Behavior

## Normalization And Scaling Model

The workflow must treat three concepts independently:

1. biological signal normalization
2. feature scaling for model fitting
3. cluster basis selection

These must not be collapsed into one generic "scale" option, because they answer different questions.

### Biological Signal Normalization

This operates on the raw or aggregated assay signal before clustering features are built.

Supported modes:

- `none`
  Preserve absolute global shifts across datasets.
  Use when broad changes in accessibility or compaction may be biologically real.

- `per_sample_global`
  Normalize each dataset by a global per-sample factor.
  Useful when technical global shifts are suspected, but this can remove true biology.

- `control_regions`
  Normalize using designated control loci expected to be stable across conditions.
  This is the preferred correction mode when suitable controls are available.

- `quantile`
  Force distributions into alignment across datasets.
  This is a strong technical correction mode and should not be the default.

### Feature Scaling

This operates only on the feature matrix passed to PCA or clustering. It does not redefine the underlying biological signal in reported outputs.

Supported modes:

- `none`
- `zscore`
- `robust_zscore`

Default: `robust_zscore`

### Cluster Basis Selection

This determines whether shared cluster boundaries should use absolute level, local pattern shape, or both.

Supported modes:

- `level_only`
  Shared boundaries are driven primarily by absolute accessibility level.
  This is sensitive to broad global shifts.

- `shape_only`
  Shared boundaries are driven by local pattern shape after removing row-level global offset.
  This is useful when asking whether local patterns change independently of global level.

- `shape_plus_level`
  Use both shape-derived features and explicit global-level features.
  This is the recommended default because it helps disentangle global compaction from local pattern reorganization.

### Recommended Default Interpretation

The first implementation should default to:

- `signal_normalization="none"`
- `feature_scaling="robust_zscore"`
- `cluster_basis="shape_plus_level"`

This preserves potential biological global shifts while still making the model numerically stable and keeping both local-pattern and global-level information available for clustering.

### Comparison-Oriented Helper Workflow

The package should eventually support a helper that runs the same shared-boundary workflow under multiple normalization or basis settings and compares the outputs.

Example shape:

```python
result = workflows.compare_shared_cluster_modes(
    samples=samples,
    modes=[
        {
            "signal_normalization": "none",
            "feature_scaling": "robust_zscore",
            "cluster_basis": "shape_plus_level",
        },
        {
            "signal_normalization": "control_regions",
            "feature_scaling": "robust_zscore",
            "cluster_basis": "shape_plus_level",
        },
        {
            "signal_normalization": "none",
            "feature_scaling": "robust_zscore",
            "cluster_basis": "shape_only",
        },
    ],
)
```

This comparison-oriented wrapper is not required for the first implementation, but the shared clustering workflow should be designed so that it can support this later without API breakage.

### Shared Validation

Before any feature extraction:

- Require at least two datasets.
- Require `sample_id`, `condition`, and `extract_h5` for every sample.
- Require consistent `motifs` and compatible extracted vector lengths across datasets.
- In `region_anchored`, require resolvable matched-region inputs for every dataset.
- Record dataset sizes and rows remaining after filtering.
- Resolve required artifacts using the selected artifact policy before rebuilding from raw inputs.
- Validate `signal_normalization`, `feature_scaling`, and `cluster_basis` as independent parameters.
- If `signal_normalization="control_regions"`, require control-region inputs or a resolvable control-region artifact.

### `read_global` Flow

1. For each sample, load read windows from extract output using existing cluster helpers.
2. Apply biological signal normalization if requested.
3. Generate read-level feature rows.
4. Apply cluster-basis transformation and feature scaling.
5. Build a balanced pooled training subset across samples.
6. Fit preprocessing on the pooled subset:
   - optional valid-site filtering
   - feature scaling
   - PCA or other dimensionality reduction
7. Fit `MiniBatchKMeans` on the pooled subset.
8. Transform and assign the full dataset for each sample in chunks.
9. Build a long-form assignments table with:
   - `sample_id`
   - `condition`
   - `replicate`
   - `cluster`
   - read metadata columns
   - optional region keys if present
10. Summarize:
   - global cluster occupancy per sample
   - cluster occupancy per condition
   - read-cluster occupancy by region when matched-region keys are available

### `region_anchored` Flow

1. For each sample, project reads onto matched regions.
2. Apply biological signal normalization if requested.
3. Aggregate region-level metrics per matched region.
4. Build one aligned feature row per region per dataset.
5. Apply cluster-basis transformation and feature scaling.
6. Build a balanced pooled training subset across datasets.
7. Fit preprocessing and one shared clustering model.
8. Assign cluster labels to all region rows.
9. Summarize region-cluster occupancy shifts across conditions.

## Performance Strategy

### Core Rules

- Fit on bounded pooled subsets, not on all rows.
- Use `MiniBatchKMeans` as the default estimator.
- Use chunked prediction for full assignment.
- Convert feature matrices to `float32` after extraction when precision permits.
- Avoid storing raw read-window matrices in the returned object unless explicitly requested.

### Sampling Strategy

The pooled training set should be balanced across datasets by default:

- choose `min(dataset_size, training_sample_per_dataset)` rows per dataset
- concatenate sampled rows
- record actual sampled counts in result metadata

Optional future extension:

- stratified sampling by region or prior label if metadata is available

### Scaling And Normalization Performance Notes

- `signal_normalization="none"` should be the cheapest mode because it avoids extra correction passes.
- `per_sample_global` should remain lightweight because it only requires one global factor per dataset.
- `control_regions` will require extra region-level resolution work but is still acceptable if control-region inputs are cached as artifacts.
- `quantile` is likely the most expensive normalization mode and should not be used implicitly.
- `shape_only` can often be implemented as a per-row centering or offset-removal transform and should remain computationally cheap.

## Outputs

### `assignments`

Long-form row-level output:

- one row per read in `read_global`
- one row per region in `region_anchored`

Minimum columns:

- `sample_id`
- `condition`
- `replicate`
- `cluster`
- mode-specific identifiers

### `cluster_distribution`

Per-sample cluster summaries:

- `sample_id`
- `condition`
- `cluster`
- `count`
- `fraction`

### `condition_distribution`

Condition-level aggregate summary:

- `condition`
- `cluster`
- `count`
- `fraction`
- `replicate_n`

### `distribution_change`

Optional comparison against a reference condition:

- `condition`
- `cluster`
- `fraction`
- `delta_fraction`
- `log2_fc`

### `cluster_profiles`

Interpretability table:

- mean feature values per cluster
- optional mean profile vectors summarized into compact columns

### `region_summaries`

Only when region mapping is possible:

- `region_id`
- `sample_id`
- `condition`
- `cluster`
- `count`
- `fraction`

### `metadata`

The top-level metadata payload should also include:

- artifact policy used
- artifact cache hits and misses
- per-sample rebuild decisions
- cohort identifier when the workflow is cohort-driven
- signal normalization mode
- feature scaling mode
- cluster basis mode
- any per-sample normalization factors used

## Standard Plots

The initial standardized plot set should be small:

1. stacked bar plot of cluster fraction by dataset
2. stacked bar plot of cluster fraction by condition
3. heatmap of cluster fraction or delta/log2FC across conditions
4. cluster profile plot for interpretation

Optional later additions:

- top-changing region heatmap
- baseline-vs-induction scatter plots
- preset study-specific figure suites

## Error Handling

The workflow should fail fast for structural issues:

- incompatible motifs across datasets
- inconsistent extracted vector lengths
- missing matched-region metadata in `region_anchored`
- zero usable rows after filtering
- duplicate `sample_id`

The workflow should warn but continue when possible for soft issues:

- a dataset contributes fewer rows than requested to the training pool
- region mapping drops a substantial number of reads
- some clusters are absent from some conditions after assignment

## Testing Strategy

### Unit Tests

Add new tests for:

- balanced pooled sampling
- fit-once / assign-many workflow behavior
- condition distribution summaries
- region summary generation
- reference-condition delta calculations
- chunked assignment equivalence versus unchunked assignment
- artifact resolution under `prefer_cached`, `require_cached`, and `rebuild`
- cohort-driven execution that reuses dataset-level artifacts
- `shape_only`, `level_only`, and `shape_plus_level` feature-path behavior
- preservation of global shifts under `signal_normalization="none"`
- expected correction behavior for `per_sample_global` and `control_regions`
- regression coverage showing that pre-existing pileup/extract-driven workflows still work after the new shared clustering layer is added

### Integration Tests

Add small synthetic multi-sample fixtures that verify:

- shared centroids are reused across samples
- `read_global` assignments can be projected back onto region summaries
- `region_anchored` produces stable per-condition cluster distributions
- cohort workflows can mix cache hits and raw-input fallbacks without changing result schema

### Regression Tests

Retain and extend the clustering tests already on the branch rather than reimplementing feature extraction from scratch.

## Documentation Strategy

Add:

- one API-oriented doc for the workflow
- one example notebook for `read_global`
- one example notebook for `region_anchored`

The collaborator notebook should be preserved as exploratory history but should not remain the primary interface for this analysis.

## Implementation Plan Shape

The implementation should be split into small steps:

1. add models
2. add artifact definitions and resolution helpers
3. add batch and cohort orchestration models
4. add distribution summaries
5. add region-aware summarization helpers
6. add workflow orchestration
7. add standardized plots
8. add tests
9. update docs/examples

## Implementation Defaults

These defaults are fixed for the first implementation:

- default clusterer: `MiniBatchKMeans`
- default cluster-boundary source: balanced pooled subset across all datasets
- default result mode: return structured tables and figures, not raw matrices
- default artifact policy: `prefer_cached`
- default signal normalization: `none`
- default feature scaling: `robust_zscore`
- default cluster basis: `shape_plus_level`
- existing low-level parsing entry points remain the continuity-preserving preprocessing surface

No experiment-specific cluster template system is required for the first version of this workflow.
