# Superpowers Docs Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Track the existing internal plan docs in git and add a central `docs/superpowers` index that maps specs and plans to their current implementation status.

**Architecture:** Keep the existing `docs/superpowers/specs` and `docs/superpowers/plans` directories unchanged, then add one central [README.md](../README.md) that explains the layout, groups the existing work by theme, and labels status centrally. This is a docs-only cleanup slice; it should not rewrite historical spec/plan content or touch runtime code.

**Tech Stack:** Markdown docs, git, existing `docs/superpowers/specs` and `docs/superpowers/plans` content

---

## File Structure

- Create: `docs/superpowers/README.md`
  - Central index for specs, plans, themes, and status labels.
- Stage existing untracked files under `docs/superpowers/plans/`
  - No content edits required unless a broken filename/path is discovered.
- Verify: `docs/superpowers/specs/`
- Verify: `docs/superpowers/plans/`

### Task 1: Create the central `docs/superpowers` index

**Files:**
- Create: `docs/superpowers/README.md`
- Reference: `docs/superpowers/specs/2026-04-02-superpowers-docs-cleanup-design.md`
- Reference: `docs/superpowers/specs/`
- Reference: `docs/superpowers/plans/`

- [ ] **Step 1: Write the new index file**

Create `docs/superpowers/README.md` with this exact starter structure:

```md
# Superpowers Docs Index

This directory tracks internal design and implementation work for the current branch.

- `specs/` contains design documents.
- `plans/` contains executable implementation breakdowns.
- Status is tracked here centrally so historical docs can remain stable snapshots.

## Status Labels

- `implemented`: the planned slice has landed in the repo.
- `partially implemented`: some related work has landed, but the broader feature family still has follow-on work.
- `planned`: the design exists, but the slice is not yet implemented.
- `not started`: a plan exists, but execution has not started.

## Specs

## Plans

## Current Themes
```

- [ ] **Step 2: Fill the `Specs` section with grouped entries**

Populate `## Specs` with grouped bullets for:

- shared clustering
- region analysis
- plotting
- cleanup/docs

Use this style:

```md
### Shared Clustering

- [2026-03-31-shared-cluster-distribution-design.md](../specs/2026-03-31-shared-cluster-distribution-design.md) - shared-boundary clustering architecture and artifact model. Status: `partially implemented`
```

Include at least these specs:

- `2026-03-31-shared-cluster-distribution-design.md`
- `2026-03-31-region-analysis-architecture-design.md`
- `2026-03-31-pre-plan-decisions-addendum.md`
- `2026-04-01-paired-region-discovery-design.md`
- `2026-04-01-plotting-axis-architecture-design.md`
- `2026-04-01-region-contrasts-plotting-design.md`
- `2026-04-02-superpowers-docs-cleanup-design.md`

- [ ] **Step 3: Fill the `Plans` section with grouped entries**

Populate `## Plans` with grouped bullets for all current plan files in `docs/superpowers/plans/`, using this style:

```md
### Plotting

- [2026-04-01-plotting-axis-architecture.md](2026-04-01-plotting-axis-architecture.md) - shared plotting-axis implementation slice. Status: `implemented`
```

List all current plan files:

- `2026-03-31-global-analysis-foundations.md`
- `2026-03-31-region-contrasts-foundations.md`
- `2026-03-31-region-discovery-foundations.md`
- `2026-03-31-shared-clustering-foundations.md`
- `2026-04-01-cluster-occupancy-region-contrasts.md`
- `2026-04-01-discovery-cluster-contrast-workflow.md`
- `2026-04-01-paired-region-discovery.md`
- `2026-04-01-plotting-axis-architecture.md`
- `2026-04-01-region-contrasts-plotting.md`
- `2026-04-01-region-discovery-cluster-workflow.md`
- `2026-04-02-superpowers-docs-cleanup.md`

- [ ] **Step 4: Fill the `Current Themes` section**

Add a short theme map that points readers to the most relevant spec/plan chains, for example:

```md
## Current Themes

- Shared clustering: start with the shared clustering design spec, then the foundations plan.
- Region contrasts: start with the region analysis architecture spec, then the contrasts foundations plan and the region-contrasts plotting plan.
- Region discovery: use the region analysis architecture spec, the paired discovery spec, and the discovery workflow plans together.
- Plotting: start with the plotting-axis spec, then the plotting-axis implementation plan and the region-contrasts plotting plan.
```

- [ ] **Step 5: Commit the new index**

```bash
git add docs/superpowers/README.md
git commit -m "docs: add superpowers docs index"
```

### Task 2: Stage and track the existing untracked plan files

**Files:**
- Stage: `docs/superpowers/plans/2026-03-31-global-analysis-foundations.md`
- Stage: `docs/superpowers/plans/2026-03-31-region-contrasts-foundations.md`
- Stage: `docs/superpowers/plans/2026-03-31-region-discovery-foundations.md`
- Stage: `docs/superpowers/plans/2026-04-01-cluster-occupancy-region-contrasts.md`
- Stage: `docs/superpowers/plans/2026-04-01-discovery-cluster-contrast-workflow.md`
- Stage: `docs/superpowers/plans/2026-04-01-paired-region-discovery.md`
- Stage: `docs/superpowers/plans/2026-04-01-plotting-axis-architecture.md`
- Stage: `docs/superpowers/plans/2026-04-01-region-contrasts-plotting.md`
- Stage: `docs/superpowers/plans/2026-04-01-region-discovery-cluster-workflow.md`

- [ ] **Step 1: Confirm the currently untracked plan set**

Run:

```bash
git status --short docs/superpowers/plans
```

Expected: the plan files above appear as untracked (`??`).

- [ ] **Step 2: Stage the plan files without editing them**

Run:

```bash
git add \
  docs/superpowers/plans/2026-03-31-global-analysis-foundations.md \
  docs/superpowers/plans/2026-03-31-region-contrasts-foundations.md \
  docs/superpowers/plans/2026-03-31-region-discovery-foundations.md \
  docs/superpowers/plans/2026-04-01-cluster-occupancy-region-contrasts.md \
  docs/superpowers/plans/2026-04-01-discovery-cluster-contrast-workflow.md \
  docs/superpowers/plans/2026-04-01-paired-region-discovery.md \
  docs/superpowers/plans/2026-04-01-plotting-axis-architecture.md \
  docs/superpowers/plans/2026-04-01-region-contrasts-plotting.md \
  docs/superpowers/plans/2026-04-01-region-discovery-cluster-workflow.md \
  docs/superpowers/plans/2026-04-02-superpowers-docs-cleanup.md
```

- [ ] **Step 3: Commit the tracked plan artifacts**

```bash
git add docs/superpowers/plans
git commit -m "docs: track superpowers implementation plans"
```

### Task 3: Verify the index and status map against the repo

**Files:**
- Verify: `docs/superpowers/README.md`
- Verify: `docs/superpowers/specs/`
- Verify: `docs/superpowers/plans/`

- [ ] **Step 1: Check that every linked file in the index exists**

Run:

```bash
rg -o 'docs/superpowers/[^)]+' docs/superpowers/README.md
```

Then verify those paths exist with:

```bash
test -f docs/superpowers/README.md
```

Expected: every file referenced by the index exists on disk.

- [ ] **Step 2: Check that the plans are no longer untracked**

Run:

```bash
git status --short docs/superpowers
```

Expected: no `??` entries remain for the tracked plan docs.

- [ ] **Step 3: Review status labels for obvious misstatements**

Manually confirm the index does not mark unfinished families as fully complete. In particular:

- plotting should reflect that the axis core and region-contrast plotting slices are implemented
- region discovery should remain `partially implemented` if follow-on contrast modes or plotting slices still remain
- global analysis should remain `partially implemented` unless the family is truly complete

- [ ] **Step 4: Commit any index corrections from verification**

```bash
git add docs/superpowers/README.md
git commit -m "docs: verify superpowers docs status index"
```

If no corrections were needed, skip this commit.

## Self-Review

- Spec coverage:
  - central index: Task 1
  - tracking untracked plans: Task 2
  - central status map and verification: Task 3
- Placeholder scan:
  - no `TODO`, `TBD`, or deferred implementation text remains in the task steps
- Type consistency:
  - all referenced paths use the existing `docs/superpowers/specs` and `docs/superpowers/plans` layout

