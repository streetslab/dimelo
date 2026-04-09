# Region Contrasts

`dimelo.region_contrasts` scores known regions from pileup-backed inputs or clustering-derived occupancy tables. It is the defined-region comparison layer for cases where you already know the loci you want to test.

## When To Use It

- Use this when you already have a BED or matched region set.
- Use `parse_bam.pileup()` first when you want locus-level motif abundance testing.
- Use `cluster` first when you want read-state or cluster-occupancy follow-up rather than average motif abundance alone.
- Use `region_discovery` when you do not yet know the loci and want de novo locus finding before follow-up contrasts.

## V1 Supported Paths

- `analysis_unit="ensemble_region"`
- `representation="modified_fraction"` or `"modified_count"`
- `signal_source="pileup_counts"`
- `test="effect_size_only"` or `"beta_binomial"`

- `analysis_unit="cluster_occupancy"`
- `signal_source="cluster_occupancy"`
- `representation="cluster_fraction"` with `test="effect_size_only"` or `"fraction_test"`
- `representation="dominant_cluster"` or `"cluster_entropy"` with `test="effect_size_only"`

Current v1 behavior:

- `effect_size_only` supports pooled `pairwise` and `group_vs_group` comparisons.
- `beta_binomial` supports the same pooled comparison modes and adds `p_value` / `adjusted_p_value`.
- `multiple_testing="fdr_bh"` is the only supported correction mode for `beta_binomial`.
- The current beta-binomial path is intentionally simple: pooled region counts with a first-pass per-region predictive test, not replicate-aware hierarchical modeling.
- `cluster_occupancy` currently supports `pairwise` and `group_vs_group` only.
- `fraction_test` is only implemented for `representation="cluster_fraction"` and uses sample-level region occupancy fractions rather than pooled reads.
- `dominant_cluster` and `cluster_entropy` are descriptive in v1. They return ranked summaries without p-values.
- Zero-read sample-region groups keep `fraction=0.0`, `cluster_entropy=0.0`, and `dominant_cluster=None` rather than inventing a winning cluster.

## Example

```python
from dimelo import region_contrasts
from dimelo.models import ContrastSpec, SampleSpec

result = region_contrasts.score_regions(
    samples=[
        SampleSpec(
            sample_id="s1",
            condition="NS",
            extract_h5="output/s1.extract.h5",
            metadata={"pileup_path": "output/s1.pileup.sorted.bed.gz"},
        ),
        SampleSpec(
            sample_id="s2",
            condition="15min",
            extract_h5="output/s2.extract.h5",
            metadata={"pileup_path": "output/s2.pileup.sorted.bed.gz"},
        ),
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

## Cluster Occupancy Example

Use this path after shared clustering when you want to test whether per-region read-state mixtures shift between conditions.

```python
from dimelo import region_contrasts, workflows
from dimelo.models import ContrastSpec

cluster_result = workflows.shared_cluster_distribution(
    samples=samples,
    mode="region_anchored",
    motifs=["A,0"],
    matched_regions="matched_regions.bed",
    n_clusters=6,
)

occupancy_result = region_contrasts.score_regions(
    samples=[],
    regions=None,
    motifs=[],
    contrast=ContrastSpec(
        mode="pairwise",
        numerator=["15min"],
        denominator=["NS"],
        reference_condition="NS",
    ),
    analysis_unit="cluster_occupancy",
    representation="cluster_fraction",
    signal_source="cluster_occupancy",
    test="fraction_test",
    occupancy_table=cluster_result.region_summaries,
)
```

## Single-Read Contrasts

`region_contrasts` now also supports `analysis_unit="single_read"` for
defined-region comparison on extract-backed read data.

First-version representations:

- `representation="read_mod_fraction"` with `signal_source="extract_reads"`
- `representation="read_window_features"` with `signal_source="extract_features"`

Supported contrast modes in v1:

- `pairwise`
- `group_vs_group`
- `matched_pairwise`

Important: reads remain observational units, but samples remain the inferential
units. The current single-read paths summarize reads within each `region x sample`
before comparing conditions, and `matched_pairwise` keeps that same boundary while
requiring explicit, non-null pairing metadata and complete pairs only.

## Canonical Outputs

The workflow returns a `RegionContrastResult` with canonical tables:

- `result.regions`
- `result.summary`
- `result.plot_data`

For `effect_size_only`, the main effect-size columns are:

- `fraction`
- `reference_fraction`
- `delta_fraction`
- `log2_fc`
- `rank`

For `representation="modified_count"`, count-based fields are also included:

- `count`
- `reference_count`
- `delta_count`
- `log2_fc_count`

For `beta_binomial`, statistical columns are included in both `result.regions` and `result.summary`:

- `p_value`
- `adjusted_p_value`

For `analysis_unit="cluster_occupancy"`, the main fields depend on representation:

- `cluster_fraction`
  - `cluster`
  - `fraction`
  - `reference_fraction`
  - `delta_fraction`
  - `log2_fc`
  - `numerator_replicate_n`
  - `denominator_replicate_n`
  - optional `p_value` / `adjusted_p_value` for `fraction_test`
- `dominant_cluster`
  - `dominant_cluster`
  - `reference_dominant_cluster`
  - `dominant_cluster_changed`
  - `numerator_replicate_n`
  - `denominator_replicate_n`
- `cluster_entropy`
  - `cluster_entropy`
  - `reference_cluster_entropy`
  - `delta_cluster_entropy`
  - `numerator_replicate_n`
  - `denominator_replicate_n`

## Custom Plotting

The results are data-first. You can use the built-in `plot_data` payloads or ignore them and plot the returned tables directly.

When you move from legacy plotters into the newer plotting-axis layer, the intended mapping is:

- `regions_5to3prime=True` -> `orientation="region_5to3"`
- `relative=True` -> fixed-window plotting around a shared anchor
- single-read views remain coordinate-preserving
- normalized segment maps are reserved for aggregate plots rather than stretched single-read canvases

```python
import seaborn as sns

sns.scatterplot(
    data=result.summary,
    x="delta_fraction",
    y="adjusted_p_value",
)
```

For region-contrast-specific plots, use the new plotting helpers with a `RegionContrastResult`, an explicit `position_table`, and the shared axis/aggregation spec:

```python
from dimelo import plotting

profile_payload = plotting.prepare_region_contrast_profile_data(
    result=result,
    position_table=position_table,
    axis=plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=2000,
        downstream_bp=2000,
    ),
    aggregation=plotting.AggregationSpec(),
)

heatmap_payload = plotting.prepare_region_contrast_heatmap_data(
    result=result,
    position_table=position_table,
    axis=plotting.AxisSpec(
        orientation="region_5to3",
        coordinate_mode="fixed_window",
        anchor="center",
        upstream_bp=2000,
        downstream_bp=2000,
    ),
    aggregation=plotting.AggregationSpec(),
)
```

`position_table` must already be an aggregated positional substrate from the parsing/loading layer, not raw per-read events. The current helper contract expects at minimum:

- `region_id`
- `position`
- `anchor`
- `value`
- `region_strand`
- either `condition` or `sample_id`

For `prepare_region_contrast_heatmap_data(...)`, `result.summary` must also provide one unambiguous `rank` value per plotted region.

For de novo locus finding before this stage, use [region discovery](region-discovery.md). For broad whole-sample summaries and normalization factors, use [global analysis](global-analysis.md).

## Preprocessing Reminder

- Run `parse_bam.pileup()` when you care about motif abundance, defined-region contrasts, or later de novo discovery.
- Run `parse_bam.extract()` when you care about single-read analysis or clustering.
- Run both when you want formal region-level abundance testing plus downstream read-level follow-up on the same samples.
- Run clustering first when you want occupancy-backed contrasts; then pass `SharedClusterResult.region_summaries` as `occupancy_table`.

## Discovery-Driven Follow-Up

When discovery, clustering, and contrasts are chained together, the combined workflow uses the selected loci as the default follow-up contract:

```python
from dimelo import workflows

result = workflows.discovery_cluster_contrast_workflow(
    samples=samples,
    motifs=["A,0"],
    genome_sizes={"chr1": 248956422},
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
```

The data contract is:

- `result.selected_regions` is the BED-style follow-up set
- clustering receives a serializable region-spec derived from those rows
- clustering-side `result.clustering.assignments["region_id"]` and `result.clustering.region_summaries["region_id"]` are normalized to the same `chr:start-end,strand` key used by `result.contrasts`
- `result.contrasts` scores the same selected loci by default
- pass `contrasts={"regions": ...}` to override that default and set `result.metadata["contrast_scope"] == "custom"` instead
- `result.metadata["full_scan_windows"]` keeps the full discovery table as context
