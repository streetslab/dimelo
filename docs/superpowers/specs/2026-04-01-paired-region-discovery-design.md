# Paired Region Discovery Design

## Goal

Extend `dimelo.region_discovery` with paired discovery modes that support:

- two-condition matched discovery
- ordered paired time-course discovery

The paired framework should preserve the current pileup-backed, deterministic tiled-window architecture while making pairing explicit and auditable.

## Why This Exists

The current `region_discovery.scan_genome(...)` implementation only supports pooled `pairwise` and `group_vs_group` contrasts. That is adequate for simple discovery, but it is not correct for common experimental designs where samples are matched across conditions:

- non-targeting versus targeting with matched replicates
- before/after on the same biological sample
- paired time-course trajectories across ordered conditions

If these are analyzed as pooled groups, the pairing structure is discarded and the result can be misleading. The paired discovery layer should answer a different biological question:

> Which windows change consistently within matched units?

## Scope

This design covers:

- new paired discovery modes in `region_discovery`
- pairing resolution rules
- paired ranking metrics for two-condition and ordered time-course discovery
- result metadata and output-table requirements

This design does not cover:

- mixed-effects or repeated-measures inferential models
- single-read paired discovery
- cluster-occupancy paired discovery
- adaptive candidate interval generation

## Architectural Rule

Build one shared paired-discovery framework and expose multiple contrast modes on top of it.

That means:

- `matched_pairwise` is a two-condition special case
- paired ordered time-course is another consumer of the same pairing engine
- pairing resolution happens once
- per-window paired summaries are the canonical substrate

This avoids creating separate, incompatible systems for paired comparisons and paired trajectories.

## Core User Model

The user still calls `region_discovery.scan_genome(...)`, but with a paired `ContrastSpec`.

Examples:

```python
ContrastSpec(
    mode="matched_pairwise",
    numerator=["targeting"],
    denominator=["nontargeting"],
    pairing_key="replicate_id",
)
```

```python
ContrastSpec(
    mode="time_course",
    time_order=["0min", "15min", "30min", "45min"],
    pairing_key="sample_id",
)
```

The key principle is:

- pooled modes compare condition groups
- paired modes compare within matched units first, then summarize across units

## Pairing Resolution

### Required Input

Paired discovery requires an explicit `pairing_key`.

The pairing key should identify the matched biological or technical unit across conditions, for example:

- replicate identity
- donor identity
- sample identity
- before/after subject identity

If the requested discovery mode is paired and `pairing_key` is missing, the workflow must fail fast.

### Default Pairing Policy

Default policy:

- `pairing_policy="complete_pairs_only"`

Behavior:

- use only matched units that have all required conditions
- drop unmatched units from scoring
- record dropped units in metadata

This is the most robust default because it preserves paired semantics without making imperfect datasets unusable.

### Strict Option

Support:

- `pairing_policy="error_on_missing"`

Behavior:

- fail if any expected matched unit is incomplete

This is useful for highly controlled analyses where partial data should be treated as an error.

### Not In V1

Do not support a silent fallback from paired discovery to pooled discovery. That would blur the biological interpretation too much.

## Internal Data Model

Add a canonical internal paired table concept:

```text
PairedWindowTable
- window_id
- chromosome
- start
- end
- strand
- pair_id
- condition
- modified_count
- valid_count
- window_fraction
```

For paired discovery:

1. build the ordinary window summary table from `global_analysis.build_window_summary(...)`
2. resolve pair membership from sample metadata using `pairing_key`
3. keep only complete matched units under the active pairing policy
4. compute paired summaries from the per-pair window table

This table is internal, but it should be conceptually stable because it will support both matched-pairwise and paired time-course modes.

## Supported Paired Modes

### `matched_pairwise`

Use when there are exactly two condition sets to compare within matched pairs.

Examples:

- non-targeting versus targeting
- before versus after

Required `ContrastSpec` fields:

- `mode="matched_pairwise"`
- `numerator`
- `denominator`
- `pairing_key`

### `time_course`

Use when the same matched unit is measured across an ordered series of conditions.

Examples:

- `0min -> 15min -> 30min -> 45min`

Required `ContrastSpec` fields:

- `mode="time_course"`
- `time_order`
- `pairing_key`

In the paired-discovery context, `time_course` means ordered paired time-course, not pooled trajectory analysis.

## Ranking Families

Two ranking families should be supported.

### 1. Magnitude-First Ranking

This emphasizes effect size.

#### Matched Pairwise

For each pair and each window:

- `delta_i = numerator_fraction_i - denominator_fraction_i`

Canonical summary columns:

- `mean_delta`
- `mean_abs_delta`
- `median_abs_delta`
- `delta_sd`
- `sign_agreement`
- `n_pairs_used`

Recommended default ranking:

- `rank_by="mean_abs_delta"`

This is the best first-pass default because it is simple, interpretable, and stable for small numbers of pairs.

#### Paired Time Course

For each pair and each window:

- compute trajectory amplitude within that pair:
  `trajectory_amplitude_i = max(fractions_i_over_time) - min(fractions_i_over_time)`

Canonical summary columns:

- `trajectory_amplitude_mean`
- `trajectory_amplitude_median`
- `trajectory_amplitude_sd`
- `n_pairs_used`

Recommended default ranking:

- `rank_by="trajectory_amplitude_mean"`

This is the right first ranking mode because it is general and does not prematurely force windows into pattern classes.

### 2. Consistency-Weighted Ranking

This emphasizes reproducibility across matched units.

#### Matched Pairwise

Candidate summary metrics:

- `consistency_weighted_delta = mean_abs_delta / (delta_sd + epsilon)`
- `mean_abs_delta * sign_agreement`

#### Paired Time Course

Candidate summary metrics:

- `trajectory_amplitude_mean / (trajectory_amplitude_sd + epsilon)`
- trajectory amplitude weighted by pairwise direction agreement between endpoints

These should be available as alternate ranking modes, but not the default.

## Recommended V1 Defaults

### Matched Pairwise

- default `rank_by="mean_abs_delta"`
- always return:
  - `mean_delta`
  - `mean_abs_delta`
  - `delta_sd`
  - `sign_agreement`
  - `n_pairs_used`

### Paired Time Course

- default `rank_by="trajectory_amplitude_mean"`
- always return:
  - `trajectory_amplitude_mean`
  - `trajectory_amplitude_sd`
  - `n_pairs_used`

These defaults balance interpretability and robustness.

## Not In First Paired Implementation

Do not add these yet:

- pattern labels like `monotonic_gain`, `transient_peak`, `late_loss`
- paired inferential models
- automatic fallback to pooled modes

Those can come later once the paired summary tables are stable.

## Result Contract Additions

`RegionDiscoveryResult` can remain the same shape, but paired discovery should add metadata fields and columns.

### Metadata

Required additions when using paired modes:

- `pairing_key`
- `pairing_policy`
- `n_pairs_used`
- `n_pairs_dropped`
- `paired_mode`
- `rank_by`

For time-course:

- `time_order`

### Hit/Window Columns

For matched pairwise windows and hits:

- `mean_delta`
- `mean_abs_delta`
- `delta_sd`
- `sign_agreement`
- `n_pairs_used`

For paired time-course windows and hits:

- `trajectory_amplitude_mean`
- `trajectory_amplitude_sd`
- `n_pairs_used`

## Error Handling

The paired workflow must fail fast when:

- `pairing_key` is missing for paired modes
- required conditions are absent
- `time_order` references conditions not present
- no complete pairs remain after filtering

It should not silently fall back to pooled analysis.

If incomplete pairs are dropped under `complete_pairs_only`, the result metadata should clearly report that.

## Testing Requirements

Add tests for:

- paired pair-resolution with complete pairs only
- strict error on missing pairs
- matched-pairwise ranking by `mean_abs_delta`
- paired time-course ranking by trajectory amplitude
- dropped-pair accounting in metadata
- failure when `pairing_key` is omitted
- failure when `time_order` is incomplete or invalid

Also add one integration test that confirms paired discovery still feeds into:

- `hits_to_bed(...)`
- `region_contrasts.score_regions(...)`

## Implementation Shape

I would implement this in two stages:

1. paired pair-resolution plus `matched_pairwise`
2. paired ordered `time_course`

Both stages should reuse the same internal pairing resolver and paired-window table builder.

## Recommendation

Implement `matched_pairwise` first, but design the internal table and ranking interfaces so `time_course` can land immediately after with minimal churn.

That gives you:

- robust paired discovery for your most common `nontargeting vs targeting` case
- a clean path to ordered paired time-course discovery
- no need to redesign the module later
