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

### Shared Clustering

- [2026-03-31-shared-cluster-distribution-design.md](specs/2026-03-31-shared-cluster-distribution-design.md) - shared-boundary clustering architecture and artifact model. Status: `partially implemented`

### Region Analysis

- [2026-03-31-region-analysis-architecture-design.md](specs/2026-03-31-region-analysis-architecture-design.md) - common architecture for region discovery, global analysis, and contrasts. Status: `partially implemented`
- [2026-03-31-pre-plan-decisions-addendum.md](specs/2026-03-31-pre-plan-decisions-addendum.md) - cross-cutting policy decisions for shared clustering and region analysis follow-on work. Status: `partially implemented`
- [2026-04-01-paired-region-discovery-design.md](specs/2026-04-01-paired-region-discovery-design.md) - paired discovery support for matched and ordered comparisons. Status: `implemented`

### Plotting

- [2026-04-01-plotting-axis-architecture-design.md](specs/2026-04-01-plotting-axis-architecture-design.md) - shared plotting-axis and renderer-neutral payload architecture. Status: `implemented`
- [2026-04-01-region-contrasts-plotting-design.md](specs/2026-04-01-region-contrasts-plotting-design.md) - plotting prep for region contrast profiles and heatmaps. Status: `implemented`

### Cleanup And Docs

- [2026-04-02-superpowers-docs-cleanup-design.md](specs/2026-04-02-superpowers-docs-cleanup-design.md) - central index and status-map cleanup for internal planning docs. Status: `implemented`

## Plans

### Shared Clustering

- [2026-03-31-shared-clustering-foundations.md](plans/2026-03-31-shared-clustering-foundations.md) - shared clustering models, artifacts, plotting payloads, and workflows. Status: `implemented`

### Region Analysis

- [2026-03-31-global-analysis-foundations.md](plans/2026-03-31-global-analysis-foundations.md) - global analysis models, workflows, and artifact plumbing. Status: `partially implemented`
- [2026-03-31-region-contrasts-foundations.md](plans/2026-03-31-region-contrasts-foundations.md) - defined-region contrast foundations and scoring paths. Status: `implemented`
- [2026-04-01-cluster-occupancy-region-contrasts.md](plans/2026-04-01-cluster-occupancy-region-contrasts.md) - cluster-occupancy evidence and scoring for region contrasts. Status: `implemented`

### Region Discovery

- [2026-03-31-region-discovery-foundations.md](plans/2026-03-31-region-discovery-foundations.md) - genome scan foundations, result models, and plotting payloads. Status: `implemented`
- [2026-04-01-paired-region-discovery.md](plans/2026-04-01-paired-region-discovery.md) - paired discovery support for matched pairwise and paired time-course scans. Status: `implemented`
- [2026-04-01-region-discovery-cluster-workflow.md](plans/2026-04-01-region-discovery-cluster-workflow.md) - workflow that chains discovery into shared clustering. Status: `implemented`
- [2026-04-01-discovery-cluster-contrast-workflow.md](plans/2026-04-01-discovery-cluster-contrast-workflow.md) - workflow that chains discovery, clustering, and contrasts. Status: `implemented`

### Plotting

- [2026-04-01-plotting-axis-architecture.md](plans/2026-04-01-plotting-axis-architecture.md) - shared plotting-axis implementation slice. Status: `implemented`
- [2026-04-01-region-contrasts-plotting.md](plans/2026-04-01-region-contrasts-plotting.md) - region contrast plotting prep for profiles and heatmaps. Status: `implemented`

### Cleanup And Docs

- [2026-04-02-superpowers-docs-cleanup.md](plans/2026-04-02-superpowers-docs-cleanup.md) - central docs index and tracked-plan cleanup. Status: `implemented`

## Current Themes

- Shared clustering: start with the shared clustering design spec, then the foundations plan, then the discovery-to-clustering workflows.
- Region contrasts: start with the region analysis architecture spec, then the contrasts foundations plan, cluster-occupancy follow-on work, and the region-contrasts plotting plan.
- Region discovery: use the region analysis architecture spec, the paired discovery spec, and the discovery workflow plans together.
- Global analysis: use the region analysis architecture spec first, then the global analysis foundations plan.
- Plotting: start with the plotting-axis spec, then the plotting-axis implementation plan and the region-contrasts plotting plan.
- Docs cleanup: use the docs cleanup design and plan for the central status map and tracked historical plans.
