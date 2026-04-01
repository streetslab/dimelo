# Region Contrasts

`dimelo.region_contrasts` scores known regions from pileup-backed inputs. It is the defined-region comparison layer for cases where you already know the loci you want to test.

## When To Use It

- Use this when you already have a BED or matched region set.
- Use `parse_bam.pileup()` first when you want locus-level motif abundance testing.
- Use `cluster` first when you want read-state or cluster-occupancy follow-up rather than average motif abundance alone.
- Use `region_discovery` later for de novo locus finding once that module is implemented.

## V1 Supported Path

- `analysis_unit="ensemble_region"`
- `representation="modified_fraction"` or `"modified_count"`
- `signal_source="pileup_counts"`
- `test="effect_size_only"` or `"beta_binomial"`

Current v1 behavior:

- `effect_size_only` supports pooled `pairwise` and `group_vs_group` comparisons.
- `beta_binomial` supports the same pooled comparison modes and adds `p_value` / `adjusted_p_value`.
- `multiple_testing="fdr_bh"` is the only supported correction mode for `beta_binomial`.
- The current beta-binomial path is intentionally simple: pooled region counts with a first-pass per-region predictive test, not replicate-aware hierarchical modeling.

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

## Custom Plotting

The results are data-first. You can use the built-in `plot_data` payloads or ignore them and plot the returned tables directly.

```python
import seaborn as sns

sns.scatterplot(
    data=result.summary,
    x="delta_fraction",
    y="adjusted_p_value",
)
```

## Preprocessing Reminder

- Run `parse_bam.pileup()` when you care about motif abundance, defined-region contrasts, or later de novo discovery.
- Run `parse_bam.extract()` when you care about single-read analysis or clustering.
- Run both when you want formal region-level abundance testing plus downstream read-level follow-up on the same samples.

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
