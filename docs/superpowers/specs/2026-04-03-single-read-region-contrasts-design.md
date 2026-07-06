# Single-Read Region Contrasts Design

## Summary

Extend `dimelo.region_contrasts` with a first `analysis_unit="single_read"` family for
defined-region comparison on extract-backed read data.

The first version should support two parallel representations:

- `representation="read_mod_fraction"`
- `representation="read_window_features"`

Both representations should remain inside `region_contrasts.score_regions(...)`,
reuse the existing explicit `contrast + analysis_unit + representation` contract,
and keep the inferential unit at the sample level rather than treating reads as
independent biological replicates.

## Goals

- Add a first `single_read` analysis path to `region_contrasts` without creating a new top-level module.
- Support direct read-level defined-region comparison on top of the existing extract-backed workflows.
- Make it explicit what kind of read quantity is being contrasted.
- Preserve the current `region_contrasts` mental model:
  one entry point, explicit contrast semantics, explicit unit of analysis.
- Support two first-class `single_read` representations in parallel:
  - `read_mod_fraction`
  - `read_window_features`
- Keep the inferential unit sample-aware rather than pooling reads as if they were replicates.
- Support built-in read-window features by default while allowing user-supplied feature tables.
- Support `pairwise`, `group_vs_group`, and `matched_pairwise` in the first version.
- Preserve continuity with the current extract-backed clustering/window feature logic.

## Non-Goals

- Create a separate `read_contrasts` module in v1.
- Add `single_read` support for every possible representation in the first release.
- Treat reads as the formal inferential unit across conditions.
- Add time-course or background-adjusted `single_read` contrasts in the first version.
- Build full single-read positional plotting in this slice.
- Replace or redefine existing `ensemble_region` or `cluster_occupancy` behavior.

## Why `single_read` Belongs Inside `region_contrasts`

`single_read` is still a defined-region contrast problem:

- the loci are already known
- the user still needs explicit contrast semantics
- the user still needs explicit unit and representation semantics

So the right public API is still:

```python
region_contrasts.score_regions(...)
```

with a new:

- `analysis_unit="single_read"`

This keeps the package coherent:

- `ensemble_region` asks whether locus-aggregated abundance changed
- `single_read` asks whether read-level behavior within known loci changed
- `cluster_occupancy` asks whether the mixture of inferred read states changed

## First-Version Representations

### `representation="read_mod_fraction"`

Each read contributes one scalar value per `region x sample`:

- modified motifs / valid motifs on that read within that region

This representation is for asking:

- did the distribution of per-read modification levels shift?
- did read-level heterogeneity change even when the locus average looked similar?

This is the strongest first `single_read` path because it is:

- easy to interpret
- close to the current motif-abundance framing
- a natural bridge between `ensemble_region` and richer read-structure analyses

### `representation="read_window_features"`

Each read contributes a feature vector per `region x sample`.

The default first source of features should be the current built-in feature set already
used in clustering/window workflows. The API should also allow a validated user-supplied
feature table override.

This representation is for asking:

- did structured read-level summaries change?
- did central-versus-flank, density-like, or shape-like read summaries reorganize?

This representation should be more flexible than `read_mod_fraction`, but also more
conservative in what it claims inferentially in v1.

## Supported Contrast Modes In v1

The first `single_read` release should support:

- `pairwise`
- `group_vs_group`
- `matched_pairwise`

It should not yet support:

- `time_course`
- `background_adjusted`
- more complex repeated-measures or mixed-model paths

### Why `matched_pairwise` Is Included Early

`single_read` summaries become scientifically much cleaner once the inferential unit
stays at the sample level. That makes matched before/after and targeting versus
nontargeting comparisons natural extensions of the first sample-aware design.

## Replication And Statistical Unit

The central statistical rule is:

- reads are observational units
- samples remain the inferential units

The implementation should therefore follow this pattern:

1. compute read-level values within each `region x sample`
2. summarize those read-level values into sample-level region statistics
3. compare those sample-level summaries across contrast sides

This avoids the major failure mode of pooled-read inference:

- pretending thousands of reads are thousands of biological replicates

### Consequence

`analysis_unit="single_read"` does **not** mean:

- every read is treated as an independent test unit across conditions

It means:

- read-level observations are used to build sample-aware contrast summaries

## Evidence Model

### Shared Required Columns

Both `single_read` representations should compile into an evidence table that preserves:

- `region_id`
- `sample_id`
- `condition`
- `read_id`
- optional `replicate`
- optional pairing metadata copied from the sample metadata when needed

The evidence table then diverges by representation.

### `read_mod_fraction` Evidence

Required columns:

- `region_id`
- `sample_id`
- `condition`
- `read_id`
- `modified_count`
- `valid_count`
- `read_mod_fraction`

Optional columns:

- `replicate`
- sample metadata fields used by pairing or grouping

### `read_window_features` Evidence

Required columns:

- `region_id`
- `sample_id`
- `condition`
- `read_id`
- one or more feature columns

Optional columns:

- `replicate`
- sample metadata fields used by pairing or grouping

## Feature Source Contract

### Built-In Default

If the user does not supply a feature table, `single_read` feature contrasts should use
the current built-in feature set already used in extract-backed clustering/window workflows.

This preserves continuity with the package’s current feature engineering path.

### User-Supplied Override

If the user supplies `feature_table=...`, the contrast layer should validate and use that
table instead of the built-in features.

Required constraints for user feature tables:

- same `region_id`
- same `sample_id`
- same `condition`
- same `read_id`
- explicit feature columns

The downstream contrast logic should not need to care whether the feature rows came from
built-in extraction or user-supplied features once validation is complete.

## API Shape

The public API should remain:

```python
result = region_contrasts.score_regions(
    samples=samples,
    regions=regions,
    contrast=contrast_spec,
    analysis_unit="single_read",
    representation="read_mod_fraction",
    signal_source="extract_reads",
    test="sample_distribution_shift",
)
```

and:

```python
result = region_contrasts.score_regions(
    samples=samples,
    regions=regions,
    contrast=contrast_spec,
    analysis_unit="single_read",
    representation="read_window_features",
    signal_source="extract_features",
    test="feature_summary_shift",
    feature_table=user_feature_table,   # optional
)
```

### Required Validation Rules

`analysis_unit="single_read"` should require:

- a supported `representation`
- a supported `signal_source`
- a supported `test`
- only supported contrast modes

## v1 Support Matrix

### `analysis_unit="single_read"`

#### `representation="read_mod_fraction"`

Allowed:

- `signal_source="extract_reads"`
- `test="effect_size_only"`
- `test="sample_distribution_shift"`
- `contrast.mode in {"pairwise", "group_vs_group", "matched_pairwise"}`

#### `representation="read_window_features"`

Allowed:

- `signal_source="extract_features"`
- `test="effect_size_only"`
- `test="feature_summary_shift"`
- `contrast.mode in {"pairwise", "group_vs_group", "matched_pairwise"}`
- built-in features by default
- user-supplied `feature_table` override

Disallowed in v1:

- `time_course`
- `background_adjusted`
- `beta_binomial`
- pooled-read inferential tests

## Statistical Meaning

### `read_mod_fraction`

This representation should support the stronger first inferential `single_read` path.

Suggested v1 outputs:

- sample-level mean read fraction
- sample-level median read fraction
- sample-level variance or spread summary
- contrast-side deltas across sample summaries

`test="sample_distribution_shift"` should be defined at the sample-summary level rather
than on pooled reads.

### `read_window_features`

This representation should start more conservatively.

Suggested v1 outputs:

- per-feature sample-level means
- per-feature sample-level medians
- per-feature effect-size summaries across contrast sides

`test="feature_summary_shift"` should primarily be a sample-summary comparison path in v1.
It may be descriptive or lightly inferential, but it should not overstate multivariate
read-feature inference.

## Result Contract

`RegionContrastResult` should remain the top-level return type.

For `analysis_unit="single_read"`, result metadata must explicitly carry:

- `analysis_unit`
- `representation`
- `signal_source`
- `test`
- `contrast.mode`

The result tables should make it obvious that the path is sample-aware:

- read-level evidence table retained where appropriate for traceability
- sample-level summary table
- contrast summary table

The summary outputs should not look interchangeable with `ensemble_region` summaries.

## Continuity Requirements

This extension must remain additive around the current parsing and extract-backed workflow model.

Required constraints:

- do not redefine `parse_bam.extract()`
- do not change existing `ensemble_region` or `cluster_occupancy` scoring semantics
- reuse the current extract-backed read or feature outputs where possible
- treat built-in read-window features as the default path for continuity
- keep the explicit `analysis_unit` and `representation` contract rather than inferring behavior

## Testing Requirements

The implementation plan for this spec must include:

- validation tests for legal and illegal `single_read` combinations
- tests for `read_mod_fraction` evidence building
- tests for `read_window_features` evidence building
- tests for built-in feature-source behavior
- tests for user-supplied feature-table validation
- tests for `pairwise`, `group_vs_group`, and `matched_pairwise`
- tests that prove sample-aware summaries are used rather than pooled-read inference
- docs updates describing what these contrasts mean biologically

## Recommended Implementation Order

1. add request validation for `analysis_unit="single_read"`
2. implement `read_mod_fraction` evidence building
3. implement `read_mod_fraction` sample-aware scoring
4. implement `read_window_features` evidence building using the built-in feature set
5. add user-supplied `feature_table` support
6. implement `matched_pairwise` support for both tracks
7. update docs and examples

## Open Items Explicitly Deferred

These are intentionally not part of the first implementation plan:

- time-course `single_read` contrasts
- background-adjusted `single_read` contrasts
- full multivariate inferential modeling for feature vectors
- single-read positional plotting
- additional `single_read` representations beyond `read_mod_fraction` and `read_window_features`
