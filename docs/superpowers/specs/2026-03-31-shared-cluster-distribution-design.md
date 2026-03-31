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

## Non-Goals

- Reproduce every ad hoc plot or helper from the collaborator notebook as package API.
- Ship experiment-specific hard-coded condition schemes as core logic.
- Recluster separately for each condition or each region.
- Store large raw matrices on result objects by default.

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

#### `dimelo/cluster.py`

Retains low-level primitives already present on the branch:

- `extract_read_windows`
- `read_window_feature_matrix`
- `cluster_read_windows`
- generic cluster profile plotting
- generic region-cluster plotting and summaries

The new workflow should call these lower-level functions rather than duplicate them.

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

## Detailed Workflow Behavior

### Shared Validation

Before any feature extraction:

- Require at least two datasets.
- Require `sample_id`, `condition`, and `extract_h5` for every sample.
- Require consistent `motifs` and compatible extracted vector lengths across datasets.
- In `region_anchored`, require resolvable matched-region inputs for every dataset.
- Record dataset sizes and rows remaining after filtering.

### `read_global` Flow

1. For each sample, load read windows from extract output using existing cluster helpers.
2. Generate read-level feature rows.
3. Build a balanced pooled training subset across samples.
4. Fit preprocessing on the pooled subset:
   - optional valid-site filtering
   - optional scaling
   - PCA or other dimensionality reduction
5. Fit `MiniBatchKMeans` on the pooled subset.
6. Transform and assign the full dataset for each sample in chunks.
7. Build a long-form assignments table with:
   - `sample_id`
   - `condition`
   - `replicate`
   - `cluster`
   - read metadata columns
   - optional region keys if present
8. Summarize:
   - global cluster occupancy per sample
   - cluster occupancy per condition
   - read-cluster occupancy by region when matched-region keys are available

### `region_anchored` Flow

1. For each sample, project reads onto matched regions.
2. Aggregate region-level metrics per matched region.
3. Build one aligned feature row per region per dataset.
4. Build a balanced pooled training subset across datasets.
5. Fit preprocessing and one shared clustering model.
6. Assign cluster labels to all region rows.
7. Summarize region-cluster occupancy shifts across conditions.

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

### Integration Tests

Add small synthetic multi-sample fixtures that verify:

- shared centroids are reused across samples
- `read_global` assignments can be projected back onto region summaries
- `region_anchored` produces stable per-condition cluster distributions

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
2. add distribution summaries
3. add region-aware summarization helpers
4. add workflow orchestration
5. add standardized plots
6. add tests
7. update docs/examples

## Implementation Defaults

These defaults are fixed for the first implementation:

- default clusterer: `MiniBatchKMeans`
- default cluster-boundary source: balanced pooled subset across all datasets
- default result mode: return structured tables and figures, not raw matrices

No experiment-specific cluster template system is required for the first version of this workflow.
