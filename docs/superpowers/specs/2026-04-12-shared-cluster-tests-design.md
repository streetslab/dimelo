# Shared Cluster Tests Design

## Summary

Add a new inference layer for shared clustering that stays separate from
`workflows.shared_cluster_distribution(...)` and consumes an existing
`SharedClusterResult`.

This layer should answer two related but distinct questions:

- global cluster-composition change across conditions or timepoints
- per-region occupancy change through the existing `region_contrasts` path

The new work in this spec is the global composition side. Region-level occupancy
inference should remain in `region_contrasts.score_regions(...)`, which already
accepts `analysis_unit="cluster_occupancy"` and consumes
`SharedClusterResult.region_summaries`.

## Goals

- Add a user-facing inference entry point for global shared-cluster composition.
- Reuse the existing `ContrastSpec` model for package-wide consistency.
- Keep replicate-aware testing as the default inferential path.
- Support pooled count tests such as chi-squared only as explicit screening
  options.
- Support both pairwise/group contrasts and ordered time-course contrasts.
- Support both unpaired and paired designs through explicit `pairing_key`
  declarations.
- Return a result object shaped similarly to current contrast-style results.
- Provide renderer-neutral `plot_data` payloads for summary and per-cluster
  follow-up views.

## Non-Goals

- Do not add inferential testing directly into `shared_cluster_distribution(...)`.
- Do not replace or subsume `region_contrasts` region-level occupancy inference.
- Do not add formal v1 testing for `dominant_cluster` or `cluster_entropy`.
- Do not implement a batch multi-contrast API in v1.
- Do not add Matplotlib renderers in the same slice unless a later plan
  explicitly calls for them.
- Do not add model-heavy compositional methods such as Dirichlet-multinomial
  regression in v1.

## Why A Separate Shared-Cluster Test Layer

`SharedClusterResult` is currently descriptive. It stores assignments,
distribution summaries, feature profiles, and optional region occupancy tables.
That is the correct boundary for the clustering workflow itself.

Formal testing should sit downstream for the same reasons the plotting layer now
sits downstream:

- clustering should remain a reusable descriptive artifact
- testing has different assumptions and failure modes than clustering
- users should be able to choose descriptive use only, inferential use only, or
  both

The design should therefore add a new analysis helper rather than overloading
the workflow that constructs `SharedClusterResult`.

## Recommended API Shape

Add a new public entry point:

```python
shared_cluster_tests(
    *,
    result: SharedClusterResult,
    contrast: ContrastSpec,
    test: str = "permutation",
    multiple_testing: str = "fdr_bh",
    n_permutations: int = 1000,
    random_state: int | None = 42,
    include_pairwise: bool = False,
) -> SharedClusterContrastResult
```

Exact keyword names can vary slightly, but the public contract should stay:

- one `SharedClusterResult` per call
- one `ContrastSpec` per call
- explicit test selection
- optional paired structure declared explicitly inside `ContrastSpec`
- one result object with canonical tables plus `plot_data`

## Why Reuse `ContrastSpec`

The package already has a contrast vocabulary in `dimelo.models.ContrastSpec`:

- `pairwise`
- `matched_pairwise`
- `group_vs_group`
- `time_course`
- `pairing_key`

Reusing `ContrastSpec` is the least confusing option because users do not need a
second contrast schema for shared-clustering inference.

Shared-cluster testing should still validate which modes are supported. In v1:

- supported:
  - `pairwise`
  - `matched_pairwise`
  - `group_vs_group`
  - `time_course`
- unsupported:
  - `background_adjusted`
  - `single_dataset`

`ContrastSpec.pairing_key` should be the only pairing declaration surface for
shared-cluster testing. The inference entry point should not add a second
parallel `pairing_key` argument.

## Analysis Boundaries

### Global Composition Testing

This new layer should test condition- or time-associated changes in overall
cluster composition using:

- `result.cluster_distribution` for sample-level cluster fractions
- `result.condition_distribution` for descriptive summaries

The biological question is:

- do samples differ consistently in their cluster mixture?

This is distinct from per-region occupancy testing.

### Region Occupancy Testing

Region occupancy inference should remain in `region_contrasts.score_regions(...)`
with:

- `analysis_unit="cluster_occupancy"`
- `signal_source="cluster_occupancy"`
- `occupancy_table=result.region_summaries`

For v1 robustness:

- formal inference remains focused on `representation="cluster_fraction"`
- `dominant_cluster` and `cluster_entropy` remain descriptive

This preserves existing architecture and avoids a confusing “one function does
everything” surface.

## Statistical Policy

### Default Inference Policy

Formal inference should be replicate-aware by default.

Reasoning:

- pooled read counts are not biological replicates
- pooled tests can become misleadingly significant with deep read counts
- sample-level cluster fractions better reflect biological variability

Therefore:

- default v1 global test path should operate at the sample level
- pooled count tests are allowed only as explicit screening options

### Pooled Screening Tests

Allowed only when the user explicitly requests them:

- `chi_squared`
- `g_test`

These tests should be clearly labeled in result metadata as pooled count
screening analyses, not the default primary inference result when replicates are
available.

### Replicate-Aware Primary Tests

The recommended v1 primary test family is permutation-based inference on
sample-level cluster fractions.

Reasoning:

- more robust than asymptotic pooled-count tests
- more flexible for awkward cluster-fraction distributions
- easier to extend to paired and time-course designs
- easier to explain than model-heavy compositional approaches

V1 should support:

- omnibus permutation testing for global composition change
- per-cluster permutation follow-up tests

### Deferred Model-Based Methods

The following should be deferred to a later version:

- Dirichlet-multinomial regression
- multinomial logistic regression
- more complex mixed models for repeated measures

Those methods may become useful later, but they add implementation and
interpretation complexity that is unnecessary for the first version.

## Supported Contrast Families

### Pairwise And Group-Vs-Group

Support:

- unpaired `pairwise`
- paired `matched_pairwise`
- `group_vs_group`

The default result should include:

- omnibus composition-level inference
- per-cluster follow-up inference

### Time-Course

Support ordered time-course testing in v1.

This should include two distinct outputs:

- ordered trend test
- omnibus “any timepoint differs” test

Optional pairwise follow-up may be included when explicitly requested, but it
should not be the primary headline result.

### Paired Longitudinal Designs

Paired time-course support should be allowed when a valid explicit pairing
structure is supplied.

The implementation does not need to solve every longitudinal model in v1; it
does need to:

- validate matched structure explicitly
- apply a paired-compatible permutation or resampling scheme
- fail clearly when the requested design is unsupported by the observed sample
  layout

## Pairing Declaration

Pairing should be declared explicitly through a `pairing_key`, not inferred from
sample IDs and not hidden inside condition group definitions.

The least confusing and most robust policy is:

- unpaired analyses: `pairing_key is None`
- paired analyses: `pairing_key` names a metadata field or canonical sample
  column identifying matched subjects/blocks

The implementation should validate:

- all requested paired samples have non-null pairing keys
- each matched pair/block has the required conditions or timepoints
- no ambiguous duplicate mapping exists within a pair/block for the requested
  test

## Result Object Shape

Add a new result model:

```python
SharedClusterContrastResult
```

It should follow the same broad style as current contrast-style outputs:

- `summary`
- `details`
- optional `pairwise`
- `plot_data`
- `metadata`

### `summary`

This is the headline table.

For pairwise/group contrasts it should usually have one row per contrast.

Recommended fields:

- `contrast_id`
- `mode`
- `test`
- `composition_effect_size`
- `effect_size_metric`
- `omnibus_p_value`
- `omnibus_adjusted_p_value` when relevant
- `trend_p_value` when relevant
- `top_cluster`
- `top_cluster_delta_fraction`
- `numerator_replicate_n`
- `denominator_replicate_n`
- `rank`

### `details`

This is the main per-cluster follow-up table.

Recommended fields for pairwise/group contrasts:

- `cluster`
- `fraction`
- `reference_fraction`
- `delta_fraction`
- `log2_fc`
- `effect_size`
- `p_value`
- `adjusted_p_value`
- `numerator_replicate_n`
- `denominator_replicate_n`
- `rank`

For time-course, `details` may additionally include:

- `timepoint`
- `mean_fraction`
- `trend_statistic`
- `trend_p_value`

### `pairwise`

Optional table produced only when requested for time-course or multi-group
follow-up.

This should stay secondary. Users should not have to read `pairwise` to answer
the primary question.

### `metadata`

This should capture:

- test family
- permutation count
- random seed when used
- multiple-testing policy
- whether the analysis is pooled or replicate-aware
- whether the design is paired
- warnings about missing replicates, degenerate clusters, or unsupported
  follow-up paths

## Effect Size Policy

The most interpretable top-level effect size for global composition is a single
overall composition-distance metric.

V1 should use a simple, explainable measure such as:

- total variation distance / L1-derived composition distance

This is preferable to leading with only a p-value or only the top changed
cluster.

The top-level summary should therefore combine:

- one overall composition effect size
- one overall inferential result
- one concise cluster-level highlight such as `top_cluster`

## Plot Data

The result should include renderer-neutral `plot_data`.

V1 should provide payloads for both:

- global summary views
- per-cluster follow-up views

Recommended payload families:

- `summary_table`
  - headline contrast-level rows for report-card style summaries
- `cluster_effect_table`
  - per-cluster fractions, deltas, effect sizes, p-values, adjusted p-values
- `time_course_table`
  - per-timepoint cluster fractions when `contrast.mode == "time_course"`
- optional `pairwise_table`
  - only when pairwise follow-up is explicitly requested

These payloads should remain renderer-neutral so later Matplotlib or external
plotting code can consume them directly.

## Validation Rules

The implementation should fail early and clearly for:

- unsupported `ContrastSpec.mode`
- missing required tables on `SharedClusterResult`
- missing requested conditions or timepoints
- paired analyses without valid `pairing_key`
- paired analyses with incomplete or ambiguous matching
- replicate-aware tests requested with too few replicates
- unsupported `multiple_testing` policy
- unsupported `test` choice for the requested contrast/design

Where possible, the function should return descriptive effect sizes even when a
formal p-value path is unavailable, but only if that behavior is explicit in the
API and metadata.

## Documentation Expectations

User-facing docs should make three boundaries explicit:

- `shared_cluster_distribution(...)` remains descriptive
- `shared_cluster_tests(...)` handles global composition inference
- `region_contrasts.score_regions(..., analysis_unit="cluster_occupancy")`
  remains the route for per-region occupancy inference

The docs should also be explicit that:

- pooled chi-squared or G-test is screening-oriented
- replicate-aware tests are the default inferential path
- region-level formal testing in v1 focuses on `cluster_fraction`

## Implementation Order

Recommended order:

1. add the new result model
2. add failing tests for global shared-cluster inference entry point
3. implement validation and sample-level summary helpers
4. implement pairwise/group replicate-aware permutation tests
5. implement pooled chi-squared / G-test screening option
6. implement time-course omnibus and trend paths
7. add `plot_data`
8. add docs

## Recommended Approach

Three approaches were considered:

1. add a dedicated `shared_cluster_tests(...)` layer
2. extend `shared_cluster_distribution(...)` directly with testing
3. fold global composition testing into `region_contrasts`

Approach 1 is recommended because it:

- preserves the descriptive boundary of clustering
- matches the existing result-centric architecture
- keeps global composition and per-region occupancy conceptually separate
- allows reuse of `ContrastSpec` without overloading unrelated region-level
  semantics

Approaches 2 and 3 would be more confusing because they would mix descriptive
workflow construction with inference, or mix global and region-level biology
under one API surface.
