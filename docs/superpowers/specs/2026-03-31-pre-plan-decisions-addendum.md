# Pre-Plan Decisions Addendum

## Purpose

This addendum resolves the main remaining policy gaps in the shared clustering and region-analysis specs before writing the implementation plan.

It should be read alongside:

- `2026-03-31-shared-cluster-distribution-design.md`
- `2026-03-31-region-analysis-architecture-design.md`

## Decisions

### 1. Replicate Handling For V1

Formal inferential `region_contrasts` in v1 should be limited to:

- `analysis_unit="ensemble_region"`
- `signal_source="pileup_counts"`
- `test="beta_binomial"`

Replicate policy in v1:

- `single_dataset`
  no replicate modeling needed
- `pairwise`
  pool counts within each side after recording replicate-level summaries
- `group_vs_group`
  pool counts within each side after recording replicate-level summaries
- `matched_pairwise`
  allowed in schema, but v1 should expose effect-size summaries only unless replicate pairing can be handled with a clearly defined paired-count procedure
- `background_adjusted`
  allowed in schema, but inferential testing should initially operate on pooled background-adjusted summaries only if the correction step is explicitly defined and audited
- `time_course`
  effect-size-only in v1, not formal significance testing

Rationale:

- this preserves a clean, defensible inferential path for v1
- it avoids pretending that pooled beta-binomial inference is valid for all replicate structures
- it leaves room for later GLM or hierarchical modeling without blocking the initial release

### 2. Authoritative V1 Support Matrix

The implementation plan should treat the following combinations as officially supported in v1:

#### Fully Supported

- `region_contrasts`
  - `analysis_unit="ensemble_region"`
  - `representation in {"modified_fraction", "modified_count"}`
  - `signal_source="pileup_counts"`
  - `test in {"beta_binomial", "effect_size_only"}`

- `shared_cluster_distribution`
  - `mode in {"read_global", "region_anchored"}`
  - `clusterer="minibatch_kmeans"`
  - `signal_normalization in {"none", "per_sample_global", "control_regions"}`
  - `feature_scaling in {"none", "robust_zscore"}`
  - `cluster_basis in {"shape_only", "shape_plus_level", "level_only"}`

#### Schema-Supported But Exploratory-Only

- `region_contrasts`
  - `analysis_unit="single_read"`
  - `analysis_unit="cluster_occupancy"`
  - `signal_source in {"region_features", "cluster_occupancy"}`
  - `time_course` with inferential testing

These combinations may return tables or effect sizes in v1, but should not promise formal significance testing unless explicitly implemented.

### 3. Artifact Compatibility Rules

Every cached artifact must carry a compatibility fingerprint.

Required metadata:

- artifact schema version
- package version
- source file path(s)
- source file digest(s) or timestamp-plus-size fingerprint
- relevant parameter hash
- upstream artifact lineage if the artifact was derived from another artifact

Workflow-level compatibility checks must include:

- normalization mode
- feature scaling mode
- cluster basis
- motifs
- window size
- region source or region digest when applicable

Rule:

- if any required compatibility field mismatches, `prefer_cached` must rebuild rather than reuse

### 4. `region_discovery` Search-Space Policy For V1

V1 discovery should use a simple, explicit tiled-window model.

Defaults:

- canonical contigs only by default
- user-configurable include/exclude contig filters
- fixed `window_size`
- fixed `step_size`
- minimum coverage threshold before scoring
- FDR computed on raw windows before merging
- adjacent significant windows merged after scoring using a fixed merge-distance rule

Rationale:

- easy to document
- easy to test
- avoids hiding discovery behavior in heuristic candidate generation

This can be extended later with adaptive candidate windows, but the first implementation should stay deterministic.

### 5. Normalization Ownership

Normalization factors should have one canonical home:

- `global_analysis` owns computation of reusable normalization diagnostics and factors
- downstream workflows may compute a factor on demand only if no compatible normalization artifact exists

Workflow rule:

- prefer canonical normalization artifacts when present
- if absent, compute on demand and record that the factor was rebuilt locally

This keeps normalization behavior consistent across workflows while preserving the fallback model.

### 6. User Guidance For `pileup` Versus `extract`

The docs and examples should include one explicit recommendation table:

- run `pileup` when you care about:
  - motif abundance
  - locus-level defined-region contrasts
  - de novo region discovery
  - broad global summaries

- run `extract` when you care about:
  - single-read pattern analysis
  - shared read clustering
  - read-state occupancy summaries

- run both when you want:
  - formal locus-level abundance testing plus downstream read-level structural follow-up
  - discovery or contrast first, then clustering on selected loci

V1 should not introduce a new required preprocessing abstraction; it should instead document the decision clearly and allow later wrappers to remain additive.

## Recommended Implementation Consequences

The implementation plan should:

1. treat pooled beta-binomial testing as the only fully inferential defined-region engine in v1
2. include an explicit support-matrix section rather than leaving legal combinations implicit
3. build artifact compatibility checks before broad cache reuse
4. keep `region_discovery` deterministic and tiled in v1
5. make normalization artifacts reusable across workflows
6. include user-facing guidance for when to run `pileup`, `extract`, or both

## Summary

These decisions keep v1 narrow where statistical claims matter, while still letting the architecture expose broader future-facing schema for single-read and occupancy-based contrasts. They also reduce the risk of cache misuse, inconsistent normalization, and user confusion about how preprocessing connects to downstream analysis.
