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

## Canonical Outputs

The workflow returns a `SharedClusterResult` with canonical tables:

- `result.assignments`
- `result.cluster_distribution`
- `result.condition_distribution`
- `result.cluster_profiles`
- `result.region_summaries` for `region_anchored`
- `result.plot_data`

Built-in plotting remains intentionally thin. `result.plot_data["cluster_distribution_bar"]` and `result.plot_data["cluster_distribution_heatmap"]` are plot-ready DataFrames, not renderer-specific figure objects.

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

## Current V1 Scope

- `mode="read_global"` supports shared read clustering from extract outputs.
- `mode="region_anchored"` supports matched-region clustering from pileup-derived region vectors.
- The default supported clusterer is `minibatch_kmeans`.
- Results are data-first: tables are the stable contract, while plotting is optional.
