# Region Analysis Architecture Design

## Summary

Add a coherent downstream analysis architecture that separates:

- preprocessing and artifact materialization in `parse_bam`
- genome-scale summarization in `global_analysis`
- de novo locus finding in `region_discovery`
- defined-region comparison in `region_contrasts`
- deeper structure analysis in `cluster`

The goal is to preserve the current region-targeted analysis logic while making it easier to support both broad global analysis and explicit defined-region comparison workflows.

## Goals

- Keep `parse_bam` as the preprocessing entry point rather than replacing it.
- Preserve the current logic that supports region-targeted downstream analysis.
- Add a global-analysis layer for broad dataset-wide summaries and normalization support.
- Reserve `region_discovery` for de novo locus finding when the user does not yet know which regions matter.
- Add `region_contrasts` for significance testing and effect-size analysis over predefined regions.
- Make it explicit what each comparison is run on:
  which datasets are compared,
  what the analysis unit is,
  and what read or region representation is being contrasted.
- Integrate cleanly with shared clustering workflows and pileup-backed summaries.

## Non-Goals

- Collapse all downstream analysis into one giant generic module.
- Hide biological interpretation behind vague names like "region test" or "single-read comparison".
- Replace existing pileup/extract logic with one monolithic artifact type.
- Implement every possible statistical test in the first version.

## Architectural Split

### `dimelo.parse_bam`

Remains the preprocessing and artifact-materialization layer.

Responsibilities:

- run `modkit`-backed pileup generation
- run read-level extraction
- validate inputs
- emit reusable outputs that downstream workflows can consume

User model:

- preprocessing/materialization lives here
- downstream interpretation does not

### `dimelo.global_analysis`

New module for broad dataset-wide summarization.

Responsibilities:

- global accessibility or modification summaries
- chromosome- or window-level aggregated tracks
- QC summaries
- global shift metrics
- control-region normalization support
- generation of broad candidate tracks for de novo discovery

This module is for asking:

- what is happening globally in this dataset?
- do we see broad compaction or decompaction?
- are there obvious global technical shifts?

### `dimelo.region_discovery`

Reserved for de novo locus finding only.

Responsibilities:

- scan broad genomic inputs without requiring predefined regions
- tile the genome or candidate intervals
- compute discovery scores on windows
- merge or refine nearby hits
- export discovered loci as BED-like outputs

This module is for asking:

- where in the genome are the most interesting regions?
- which loci stand out within one dataset?
- which loci change strongly between matched datasets when regions were not known in advance?

If the user already has a BED or defined region set, they should not use `region_discovery`.

### `dimelo.region_contrasts`

New module for defined-region comparison and testing.

Responsibilities:

- summarize signal over known regions
- compute effect sizes across matched datasets
- run region-level statistical tests such as beta-binomial
- support pairwise, grouped, background-adjusted, and time-course contrasts
- export ranked results for downstream region-targeted workflows

This module is for asking:

- how do these known loci differ across conditions?
- which predefined regions change significantly?
- are observed differences at known loci likely to reflect technical scaling or real biology?

### `dimelo.cluster`

Retains responsibility for structure analysis once reads or regions have been selected.

Responsibilities:

- shared-boundary clustering across datasets
- read-global and region-anchored clustering workflows
- cluster occupancy summaries
- integration with region-level summaries and contrasts

This module is for asking:

- what structure exists among reads or regions?
- do read states or region states reorganize across conditions?

## User-Facing Flow

The intended flow should be:

1. `parse_bam`
   materialize reusable preprocessing outputs
2. `global_analysis`
   summarize genome-wide behavior and broad technical or biological shifts
3. `region_discovery`
   discover candidate loci when regions are not known
4. `region_contrasts`
   score and test known loci
5. `cluster`
   perform deep structural analysis on reads or regions of interest

This flow should be composable rather than mandatory. A user may enter at `region_contrasts` with an existing BED, or enter at `cluster` using a previously selected set of regions.

## Why `region_discovery` And `region_contrasts` Must Be Separate

These are different scientific tasks:

- discovery asks where signal or change exists
- contrasts ask what happens at loci already defined

If they are collapsed together:

- de novo scanning logic gets mixed with BED-based testing
- statistical assumptions become unclear
- it becomes hard to tell whether a result was discovered or tested on predefined loci

The package should therefore reserve:

- `region_discovery` for unknown loci
- `region_contrasts` for known loci

## `region_contrasts` Core API

### Primary Entry Point

```python
result = region_contrasts.score_regions(
    samples=samples,
    regions=regions,
    contrast=contrast_spec,
    analysis_unit="ensemble_region",
    representation="modified_fraction",
    signal_source="pileup_counts",
    test="beta_binomial",
    multiple_testing="fdr_bh",
)
```

The API must make these dimensions explicit:

1. `contrast`
   which datasets or conditions are compared
2. `analysis_unit`
   what the rows represent
3. `representation`
   what property of those rows is being contrasted

## Contrast Specification

### `ContrastSpec`

```python
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
```

### Supported Contrast Modes

- `single_dataset`
  rank regions within one dataset
- `pairwise`
  compare two datasets directly
- `matched_pairwise`
  compare matched samples using an explicit pairing key
- `group_vs_group`
  compare multiple samples or conditions against another group
- `background_adjusted`
  compare signal after accounting for matched background measurements
- `time_course`
  compare ordered condition trajectories

The same `ContrastSpec` structure should support all of these rather than fragmenting into unrelated function families.

## Analysis Unit Model

The package must be explicit about what observational level is being compared.

### `analysis_unit="ensemble_region"`

Rows represent aggregated region-level signal.

Examples:

- region-level modified counts
- region-level valid counts
- region-level modified fraction
- region-level enrichment summaries

This is the default path for significance testing with beta-binomial models.

### `analysis_unit="single_read"`

Rows represent individual reads within predefined regions.

Examples:

- per-read modification fraction
- read-level density profile summaries
- read window features
- read span or length metrics

This is for within-region distributional changes, not simple count-based testing.

### `analysis_unit="cluster_occupancy"`

Rows represent region-level summaries derived from clustered reads or clustered regions.

Examples:

- fraction of reads assigned to a cluster within a region
- dominant cluster fraction
- cluster entropy
- state-mixture proportions

This is the bridge between clustering and region-level contrasts.

## Representation Model

Even within a given analysis unit, users need to know what property is being compared.

### Ensemble-Region Representations

- `modified_fraction`
- `modified_count`
- `coverage`
- `enrichment`
- `depth_profile_summary`

### Single-Read Representations

- `read_mod_fraction`
- `read_density_profile`
- `read_window_features`
- `read_shape`
- `read_length_or_span`
- `motif_specific_fraction`

### Cluster-Occupancy Representations

- `cluster_fraction`
- `dominant_cluster`
- `cluster_entropy`
- `state_mixture`

## Why `analysis_unit` And `representation` Must Both Be Explicit

Without these fields, a user cannot tell whether a result means:

- the average region became more or less accessible
- the distribution of reads within the region changed
- the mixture of read states changed
- the local pattern shape changed while the mean remained constant

These are biologically different conclusions and should never be hidden behind a generic label like "region contrast".

## Signal Sources

The first implementation should support multiple signal substrates with a clear default.

### `signal_source="pileup_counts"`

Primary default.

Uses per-region modified and valid counts.

This is the natural substrate for:

- beta-binomial testing
- region-level effect sizes
- broad defined-region comparison

### `signal_source="region_features"`

Optional later mode.

Uses region-level summaries derived from read-level artifacts or feature engineering.

This is useful when pileup counts are not the right abstraction.

### `signal_source="cluster_occupancy"`

Integration mode with clustering workflows.

Uses region-level cluster proportions or occupancy summaries as input to defined-region contrasts.

## Statistical Engines

### First-Pass Default

Support `test="beta_binomial"` for:

- `analysis_unit="ensemble_region"`
- `signal_source="pileup_counts"`

This should be the first formal significance-testing path because it aligns with the current count-oriented architecture and the requested statistical direction.

### Additional First-Version Mode

Support `test="effect_size_only"` for fast exploratory use without p-values.

### Planned But Not Required For V1

- distribution-shift tests for `single_read`
- occupancy/fraction tests for `cluster_occupancy`
- richer time-course inferential models

These should be designed into the schema now, but not block the first implementation.

## `region_discovery` Core API

### Primary Entry Point

```python
result = region_discovery.scan_genome(
    samples=samples,
    contrast=contrast_spec,
    window_size=2000,
    step_size=500,
    signal_source="pileup_counts",
    score="beta_binomial",
    merge_hits=True,
)
```

`region_discovery` should:

- generate candidate windows or intervals
- score them using a discovery-specific engine
- rank them
- optionally merge or refine neighboring hits
- export discovered intervals into BED-like outputs

These outputs should be directly consumable by `region_contrasts` and `cluster`.

## Shared Evidence Table

Both `region_discovery` and `region_contrasts` should rely on a common internal region evidence table when possible.

Example schema:

```python
RegionEvidence
- region_id
- chrom
- start
- end
- sample_id
- condition
- replicate
- modified_count
- valid_count
- modified_fraction
- coverage
- analysis_unit
- representation
- signal_source
- metadata...
```

This shared evidence table keeps discovery and contrast workflows comparable without conflating their user-facing roles.

## Integration With Clustering

### `region_contrasts -> cluster`

Defined-region contrasts should be able to export ranked or significant BEDs for downstream deep clustering analysis.

Examples:

- take top significant regions from beta-binomial testing
- cluster reads within those selected loci
- compare region-level cluster occupancy afterward

### `cluster -> region_contrasts`

Clustering workflows should be able to emit region-level occupancy artifacts that `region_contrasts` can score.

Examples:

- test whether cluster C3 occupancy differs between conditions at known regions
- compare cluster entropy across matched loci

### `global_analysis -> region_discovery`

Global summaries should provide coarse substrate for discovery:

- broad tracks
- normalization factors
- window-level summaries

### `region_discovery -> region_contrasts`

Discovered loci should be easy to hand off into formal defined-region testing and ranking workflows.

## Module Layout

### `dimelo/global_analysis.py`

Responsibilities:

- broad summaries
- QC and global scaling diagnostics
- control-region normalization support

### `dimelo/region_discovery.py`

Responsibilities:

- de novo scans
- candidate interval generation
- discovery scoring
- hit merging and refinement

### `dimelo/region_contrasts.py`

Responsibilities:

- defined-region evidence building
- contrast specification handling
- beta-binomial and effect-size scoring
- result tables and export helpers

### `dimelo/cluster.py`

Responsibilities:

- structure analysis
- cluster occupancy outputs
- integration with defined-region analyses

## Output Contracts

### `RegionContrastResult`

```python
@dataclass
class RegionContrastResult:
    regions: pd.DataFrame
    summary: pd.DataFrame
    contrast: ContrastSpec
    metadata: dict[str, Any]
    figures: dict[str, Any]
```

The metadata must include:

- contrast mode
- analysis unit
- representation
- signal source
- test
- normalization mode

### `RegionDiscoveryResult`

```python
@dataclass
class RegionDiscoveryResult:
    hits: pd.DataFrame
    windows: pd.DataFrame
    contrast: ContrastSpec | None
    metadata: dict[str, Any]
    figures: dict[str, Any]
```

The outputs should be easy to export as BEDs and easy to pass into follow-on workflows.

## Testing Strategy

### Unit Tests

Add tests for:

- `ContrastSpec` validation across supported modes
- legal and illegal combinations of `analysis_unit`, `representation`, `signal_source`, and `test`
- region evidence table construction
- effect-size-only ranking
- beta-binomial path over predefined region counts
- handoff from discovered regions into defined-region contrasts

### Integration Tests

Add synthetic tests that verify:

- a de novo scan can produce intervals that are then consumed by `region_contrasts`
- a clustering workflow can emit cluster occupancy summaries that `region_contrasts` can score
- metadata clearly records what was contrasted and at what observational level

## Documentation Strategy

Add:

- one doc for global summaries and preprocessing handoff
- one doc for de novo region discovery
- one doc for defined-region contrasts
- one doc showing how discovered or contrasted regions feed into clustering

## Implementation Plan Shape

The implementation should be split into small pieces:

1. add `ContrastSpec` and result models
2. add `global_analysis` scaffolding
3. add `region_contrasts` evidence-building and effect-size workflows
4. add beta-binomial testing for `ensemble_region` plus `pileup_counts`
5. add `region_discovery` scanning and BED export
6. add clustering integration points
7. add tests
8. add docs/examples

## Implementation Defaults

These defaults are fixed for the first implementation:

- `region_discovery` is for de novo loci only
- `region_contrasts` is for known regions only
- default `analysis_unit`: `ensemble_region`
- default `representation`: `modified_fraction`
- default `signal_source`: `pileup_counts`
- default `test`: `beta_binomial` when the combination is valid

This architecture should make it obvious what is being compared and why, while still supporting clean handoff into clustering and other downstream workflows.
