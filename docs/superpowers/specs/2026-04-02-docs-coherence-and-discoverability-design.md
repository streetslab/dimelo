# Docs Coherence And Discoverability Design

## Goal

Make the current branch documentation coherent and consistent across the internal planning docs, user-facing guides, and the small package surface users are expected to import from.

This is a cleanup slice, not a feature slice. The goal is to make the repository easier to understand after the large amount of recent implementation work, without changing analysis behavior.

## Scope

This design covers:

- internal status-map accuracy in [docs/superpowers/README.md](../README.md)
- consistency across user-facing guides:
  - [README.md](../../../README.md)
  - [docs/global-analysis.md](../../../docs/global-analysis.md)
  - [docs/region-discovery.md](../../../docs/region-discovery.md)
  - [docs/region-contrasts.md](../../../docs/region-contrasts.md)
  - [docs/shared-clustering.md](../../../docs/shared-clustering.md)
- small package-surface cleanup in [dimelo/__init__.py](../../../dimelo/__init__.py) only if docs and the current public import story are clearly misaligned

This design does not cover:

- new analysis features
- new plotting features
- API redesign
- broad refactors
- rewriting the historical specs/plans themselves beyond status-map updates

## Why This Pass Is Needed

The branch now has:

- multiple new analysis layers
- multiple new plotting-prep helper families
- discovery-to-clustering and discovery-to-clustering-to-contrast workflows
- internal planning docs that were cleaned up once already, but are now behind again

The implementation is ahead of the documentation in a few places. The repo needs one coherence pass so a user or future contributor can understand:

- what is implemented now
- which docs correspond to which module/workflow
- how the modules fit together
- which plotting helpers exist for which analysis families

## Definition Of “Coherent And Consistent”

This pass should satisfy four tests.

### 1. Status-map accuracy

[docs/superpowers/README.md](../README.md) should:

- include the newest plotting specs/plans
- stop omitting implemented slices
- label implemented vs partially implemented work in a way that matches repo reality

The status map should remain lightweight. It is an index, not a changelog.

### 2. User-flow consistency

The user-facing docs should tell one consistent story:

- `parse_bam` remains the preprocessing/materialization foundation
- `global_analysis`, `region_discovery`, `region_contrasts`, and `shared_clustering` are the higher-level analysis layers
- plotting is data-first
- helper coverage differs by analysis family, and the docs should say so honestly

### 3. Cross-link and naming consistency

The same helpers and workflows should be named the same way everywhere they appear. In practice, this means consistently using names such as:

- `shared_cluster_distribution(...)`
- `discovery_cluster_workflow(...)`
- `discovery_cluster_contrast_workflow(...)`
- `prepare_region_contrast_profile_data(...)`
- `prepare_region_contrast_heatmap_data(...)`
- `prepare_region_discovery_scan_data(...)`
- `prepare_region_discovery_hit_context_data(...)`
- `prepare_global_analysis_summary_data(...)`
- `prepare_global_analysis_window_data(...)`

### 4. Package discoverability consistency

If the docs imply that users should import a module or helper from top-level `dimelo`, the top-level package surface should support that expectation, but only where the fix is small and obvious.

This does not mean re-exporting every helper automatically. It means removing clear documentation/import mismatches if they exist.

## Recommended Cleanup Structure

The cleanup should proceed in three layers.

### Layer 1: Internal status/docs map

Update [docs/superpowers/README.md](../README.md) to:

- add the missing plotting specs/plans
- refresh the status labels for the slices that have landed
- keep the wording concise

This is the source of truth for planning status and should reflect the current branch state.

### Layer 2: User-facing guide consistency

Sweep the main user-facing docs:

- [README.md](../../../README.md)
- [docs/global-analysis.md](../../../docs/global-analysis.md)
- [docs/region-discovery.md](../../../docs/region-discovery.md)
- [docs/region-contrasts.md](../../../docs/region-contrasts.md)
- [docs/shared-clustering.md](../../../docs/shared-clustering.md)

Focus on:

- stale “future” wording for already implemented work
- inconsistent helper/workflow naming
- inconsistent descriptions of module roles
- inconsistent plotting-helper coverage claims
- missing short cross-links where users would reasonably look first

This should be a targeted consistency sweep, not a rewrite.

### Layer 3: Small package-surface cleanup

Inspect [dimelo/__init__.py](../../../dimelo/__init__.py) after the docs sweep.

Allowed changes:

- exposing a clearly user-facing module that the docs already treat as a top-level import
- fixing an obvious export omission that makes the docs misleading

Disallowed changes:

- broad export churn
- re-exporting many helpers just for convenience
- changing public import style for aesthetic reasons

## Expected Outcomes

After the pass:

- the internal docs index reflects the actual implemented state
- the top-level README and guide docs tell the same architecture story
- plotting coverage is described consistently:
  - `region_contrasts` plotting prep exists
  - `region_discovery` plotting prep exists
  - `global_analysis` plotting prep exists
  - `shared_clustering` still has lighter plotting coverage by comparison
- if a small `__init__.py` fix is needed, the docs and package surface agree

## Non-Goals

This pass should not try to solve every remaining gap in the branch. In particular, it should not turn into:

- a renderer implementation slice
- a `shared_clustering` plotting expansion
- a `single_read` contrast design
- a general docs rewrite

Those are follow-on projects and should remain separate.

## Validation

The cleanup should validate:

- the internal index links still resolve
- the repo docs do not contain stale contradictory statements about which modules are implemented
- any changed top-level exports still allow `import dimelo` cleanly

## Recommendation

Do this as a docs-first coherence sweep with three ordered layers:

1. internal status-map update
2. user-facing guide consistency sweep
3. minimal package-surface cleanup only if the docs clearly demand it

That keeps the pass focused and lets us build the next implementation plan from a clean, accurate documentation state.
