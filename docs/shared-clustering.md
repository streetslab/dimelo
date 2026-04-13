# Shared Clustering

`dimelo.workflows.shared_cluster_distribution()` provides a higher-level shared-boundary workflow on top of the existing parsing layer. It does not replace `parse_bam.extract()` or `parse_bam.pileup()`; it consumes those outputs and returns canonical tables that you can plot with Matplotlib, seaborn, Plotly, Altair, or your own stack.

## Preprocessing Choices

- Run `parse_bam.extract()` when you want read-level clustering or single-read pattern analysis.
- Run `parse_bam.pileup()` when you want matched-region clustering from pileup-derived summaries.
- Run both when you want region-level summaries plus read-level follow-up on the same samples.

## Workflow Entry Point

### Read-global clustering

```python
from dimelo import workflows
from dimelo.models import SampleSpec

result = workflows.shared_cluster_distribution(
    samples=[
        SampleSpec(sample_id="s1", condition="NS", extract_h5="output/s1.extract.h5"),
        SampleSpec(sample_id="s2", condition="15min", extract_h5="output/s2.extract.h5"),
    ],
    mode="read_global",
    motifs=["A,0"],
    n_clusters=8,
)
```

### Region-anchored clustering

```python
from dimelo import workflows
from dimelo.models import SampleSpec

result = workflows.shared_cluster_distribution(
    samples=[
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="output/s1.extract.h5",
            regions_bed="controls/s1_controls.bed",
            metadata={"pileup_path": "output/s1.pileup.sorted.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="output/s2.extract.h5",
            regions_bed="controls/s2_controls.bed",
            metadata={"pileup_path": "output/s2.pileup.sorted.bed.gz"},
        ),
    ],
    mode="region_anchored",
    motifs=["A,0"],
    matched_regions="matched_regions.bed",
    n_clusters=6,
)
```

`shared_cluster_distribution()` expects `matched_regions` to already be a serializable region spec. In the discovery-driven flow, `discovery_cluster_workflow()` is the caller that derives that spec from BED-style discovery rows and forwards it into this workflow. That keeps discovery output human-readable while still feeding clustering the compact region representation it needs.

`discovery_cluster_contrast_workflow()` extends the same handoff one step further: it derives the serializable region spec from the BED-style discovery rows, passes that into clustering, normalizes clustering-side `region_id` values to the same `chr:start-end,strand` format used by contrasts, and uses the same selected loci as the default follow-up region set for contrasts unless `contrasts["regions"]` overrides it.

For defined-region occupancy follow-up, pass `result.region_summaries` into `region_contrasts.score_regions(..., analysis_unit="cluster_occupancy", signal_source="cluster_occupancy", occupancy_table=...)`. That keeps the clustering output as the reusable source of truth for per-region read-state mixtures.

## Canonical Outputs

The workflow returns a `SharedClusterResult` with canonical tables:

- `result.assignments`
- `result.cluster_distribution`
- `result.condition_distribution`
- `result.cluster_profiles`
- `result.region_summaries` for `region_anchored`
- `result.plot_data`

`result.region_summaries` is also the standard occupancy handoff into `dimelo.region_contrasts` when you want per-region `cluster_fraction`, `dominant_cluster`, or `cluster_entropy` follow-up instead of only global cluster distributions.

Shared clustering now follows the same layered plotting structure as the other newer analysis families:

- canonical tables on the result object
- renderer-neutral prep helpers in `dimelo.plotting`
- optional built-in Matplotlib renderers in `dimelo.plotting_matplotlib`

## Plotting Prep Helpers

`SharedClusterResult` now has three renderer-neutral plotting helpers in `dimelo.plotting`:

- `prepare_shared_cluster_distribution_data(result=...)`
- `prepare_shared_cluster_profile_data(result=...)`
- `prepare_shared_cluster_region_data(result=...)`

These helpers sit on top of the canonical result tables:

- `result.cluster_distribution`
- `result.condition_distribution`
- `result.distribution_change`
- `result.cluster_profiles`
- `result.region_summaries`

The older lightweight payloads in `result.plot_data["cluster_distribution_bar"]` and
`result.plot_data["cluster_distribution_heatmap"]` remain supported for backward familiarity.

The plotting contract stays compatibility-safe with older package versions:

- legacy Matplotlib wrappers still work
- tables in `result.plot_data` are the stable plotting substrate
- `regions_5to3prime` now maps onto the shared `orientation="region_5to3"` axis concept for region-aligned plots
- fixed-window plotting is available for both aggregate and single-read views
- scaled or segmented metaregion axes are aggregate-only and are not used for single-read plots

## Built-In Matplotlib Renderers

`dimelo.plotting_matplotlib` also provides an optional built-in renderer layer on
top of these prepared payloads. This is a convenience layer, not a replacement
for the stable tables in `result.plot_data` or the renderer-neutral prep helpers.

Cluster distribution example:

```python
from dimelo import plotting, plotting_matplotlib

distribution_payload = plotting.prepare_shared_cluster_distribution_data(result=result)
fig, ax = plotting_matplotlib.plot_shared_cluster_distribution_matplotlib(
    distribution_payload,
    level="sample",
)
plotting_matplotlib.save_figure(fig, "shared-cluster-distribution.png")
```

Condition-level cluster-change heatmap:

```python
fig, ax = plotting_matplotlib.plot_shared_cluster_change_matplotlib(
    distribution_payload,
)
```

Cluster profile heatmap example:

```python
profile_payload = plotting.prepare_shared_cluster_profile_data(result=result)
fig, ax = plotting_matplotlib.plot_shared_cluster_profile_heatmap_matplotlib(profile_payload)
```

Cluster profile series example:

```python
fig, ax = plotting_matplotlib.plot_shared_cluster_profile_series_matplotlib(profile_payload)
```

Region occupancy example:

```python
region_payload = plotting.prepare_shared_cluster_region_data(result=result)
fig, ax = plotting_matplotlib.plot_shared_cluster_region_matplotlib(
    region_payload,
    level="condition",
)
```

The default region occupancy view is condition-level because it is usually the
main biological summary. Sample-level occupancy remains available when you want
replicate inspection or QC:

```python
fig, ax = plotting_matplotlib.plot_shared_cluster_region_matplotlib(
    region_payload,
    level="sample",
)
```

The older `result.plot_data["cluster_distribution_bar"]` and
`result.plot_data["cluster_distribution_heatmap"]` payloads remain supported for
users updating from earlier package versions.

## Custom Plotting

```python
import seaborn as sns

bar_data = result.plot_data["cluster_distribution_bar"]
sns.barplot(data=bar_data, x="sample_id", y="fraction", hue="cluster")
```

For heatmaps:

```python
heatmap = result.plot_data["cluster_distribution_heatmap"].set_index("condition")
```

## Global Composition Inference

Shared clustering now includes a dedicated inference helper for global
composition changes:

```python
from dimelo import shared_cluster_tests
from dimelo.models import ContrastSpec

contrast_result = shared_cluster_tests.shared_cluster_tests(
    result=result,
    contrast=ContrastSpec(
        mode="pairwise",
        numerator=["treated"],
        denominator=["NS"],
    ),
    test="permutation",
)
```

Matched pairwise example:

```python
contrast_result = shared_cluster_tests.shared_cluster_tests(
    result=result,
    contrast=ContrastSpec(
        mode="matched_pairwise",
        numerator=["treated"],
        denominator=["NS"],
        pairing_key="subject_id",
    ),
    test="permutation",
)
```

Time-course example:

```python
contrast_result = shared_cluster_tests.shared_cluster_tests(
    result=result,
    contrast=ContrastSpec(mode="time_course", time_order=["t0", "t1", "t2"]),
    test="permutation",
    include_pairwise=True,
)
```

Pooled screening example:

```python
screen_result = shared_cluster_tests.shared_cluster_tests(
    result=result,
    contrast=ContrastSpec(mode="pairwise", numerator=["treated"], denominator=["NS"]),
    test="chi_squared",
)
```

Inference boundary:

- `shared_cluster_tests(...)` is for global cluster-composition inference from
  `SharedClusterResult`.
- `region_contrasts.score_regions(..., analysis_unit="cluster_occupancy")`
  remains the path for region-level occupancy inference.
- Pooled `chi_squared` and `g_test` are screening-oriented outputs; use
  replicate-aware permutation paths as the primary inferential mode when
  replicate structure exists.

## Current V1 Scope

- `mode="read_global"` supports shared read clustering from extract outputs.
- `mode="region_anchored"` supports matched-region clustering from pileup-derived region vectors.
- `discovery_cluster_workflow()` composes discovery and `region_anchored` clustering by passing derived matched regions into this workflow.
- `discovery_cluster_contrast_workflow()` composes discovery, `region_anchored` clustering, and defined-region contrasts on the selected loci by default, with an explicit custom contrast-region override.
- `result.region_summaries` is the v1 handoff table for `cluster_occupancy` region contrasts.
- The default supported clusterer is `minibatch_kmeans`.
- Results are data-first: tables are the stable contract, while plotting is optional.
