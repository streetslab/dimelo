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
- [2026-04-03-single-read-region-contrasts-design.md](specs/2026-04-03-single-read-region-contrasts-design.md) - first `single_read` contrast family for defined-region read-level summaries. Status: `implemented`

### Plotting

- [2026-04-01-plotting-axis-architecture-design.md](specs/2026-04-01-plotting-axis-architecture-design.md) - shared plotting-axis and renderer-neutral payload architecture. Status: `implemented`
- [2026-04-01-region-contrasts-plotting-design.md](specs/2026-04-01-region-contrasts-plotting-design.md) - plotting prep for region contrast profiles and heatmaps. Status: `implemented`
- [2026-04-02-region-discovery-plotting-design.md](specs/2026-04-02-region-discovery-plotting-design.md) - plotting prep for discovery scan-overview and hit-context payloads. Status: `implemented`
- [2026-04-02-global-analysis-plotting-design.md](specs/2026-04-02-global-analysis-plotting-design.md) - plotting prep for global summary and broad-window payloads. Status: `implemented`

### Cleanup And Docs

- [2026-04-02-superpowers-docs-cleanup-design.md](specs/2026-04-02-superpowers-docs-cleanup-design.md) - central index and status-map cleanup for internal planning docs. Status: `implemented`
- [2026-04-02-docs-and-integration-cleanup-design.md](specs/2026-04-02-docs-and-integration-cleanup-design.md) - path portability and discoverability cleanup for tracked superpowers docs and related integration surfaces. Status: `implemented`
- [2026-04-02-docs-coherence-and-discoverability-design.md](specs/2026-04-02-docs-coherence-and-discoverability-design.md) - coherence pass across internal and user-facing docs plus narrow discoverability cleanup. Status: `partially implemented`

## Plans

### Shared Clustering

- [2026-03-31-shared-clustering-foundations.md](plans/2026-03-31-shared-clustering-foundations.md) - shared clustering models, artifacts, plotting payloads, and workflows. Status: `implemented`

### Region Analysis

- [2026-03-31-global-analysis-foundations.md](plans/2026-03-31-global-analysis-foundations.md) - global analysis models, workflows, and artifact plumbing. Status: `implemented`
- [2026-03-31-region-contrasts-foundations.md](plans/2026-03-31-region-contrasts-foundations.md) - defined-region contrast foundations and scoring paths. Status: `implemented`
- [2026-04-01-cluster-occupancy-region-contrasts.md](plans/2026-04-01-cluster-occupancy-region-contrasts.md) - cluster-occupancy evidence and scoring for region contrasts. Status: `implemented`
- [2026-04-03-single-read-region-contrasts.md](plans/2026-04-03-single-read-region-contrasts.md) - `single_read` evidence builders and sample-aware scoring for read-level region contrasts. Status: `implemented`

### Region Discovery

- [2026-03-31-region-discovery-foundations.md](plans/2026-03-31-region-discovery-foundations.md) - genome scan foundations, result models, and plotting payloads. Status: `implemented`
- [2026-04-01-paired-region-discovery.md](plans/2026-04-01-paired-region-discovery.md) - paired discovery support for matched pairwise and paired time-course scans. Status: `implemented`
- [2026-04-01-region-discovery-cluster-workflow.md](plans/2026-04-01-region-discovery-cluster-workflow.md) - workflow that chains discovery into shared clustering. Status: `implemented`
- [2026-04-01-discovery-cluster-contrast-workflow.md](plans/2026-04-01-discovery-cluster-contrast-workflow.md) - workflow that chains discovery, clustering, and contrasts. Status: `implemented`

### Plotting

- [2026-04-01-plotting-axis-architecture.md](plans/2026-04-01-plotting-axis-architecture.md) - shared plotting-axis implementation slice. Status: `implemented`
- [2026-04-01-region-contrasts-plotting.md](plans/2026-04-01-region-contrasts-plotting.md) - region contrast plotting prep for profiles and heatmaps. Status: `implemented`
- [2026-04-02-region-discovery-plotting.md](plans/2026-04-02-region-discovery-plotting.md) - region discovery plotting prep for per-contig scans and local hit-context views. Status: `implemented`
- [2026-04-02-global-analysis-plotting.md](plans/2026-04-02-global-analysis-plotting.md) - global analysis plotting prep for summary and broad-window payloads. Status: `implemented`

### Cleanup And Docs

- [2026-04-02-superpowers-docs-cleanup.md](plans/2026-04-02-superpowers-docs-cleanup.md) - central docs index and tracked-plan cleanup. Status: `implemented`
- [2026-04-02-docs-and-integration-cleanup.md](plans/2026-04-02-docs-and-integration-cleanup.md) - historical path normalization plus narrow discoverability cleanup around completed analysis and plotting work. Status: `implemented`
- [2026-04-02-docs-coherence-and-discoverability.md](plans/2026-04-02-docs-coherence-and-discoverability.md) - coherence pass across internal and user-facing docs plus narrow package-surface cleanup if needed. Status: `partially implemented`

## Current Themes

- Shared clustering: start with the shared clustering design spec, then the foundations plan, then the discovery-to-clustering workflows.
- Region contrasts: start with the region analysis architecture spec, then the contrasts foundations plan, cluster-occupancy follow-on work, and the region-contrasts plotting plan.
- Single-read contrasts: use the single-read region-contrasts design and plan after the region-contrasts foundations work when you need extract-backed defined-region comparison.
- Region discovery: use the region analysis architecture spec, the paired discovery spec, and the discovery workflow plans together.
- Global analysis: use the region analysis architecture spec first, then the global analysis foundations plan.
- Plotting: start with the plotting-axis spec, then the plotting-axis implementation plan, then the region-contrasts, region-discovery, and global-analysis plotting plans for helper-level coverage.
- Docs cleanup: use the docs cleanup design and plan for the central status map and tracked historical plans, then the docs-and-integration and docs-coherence follow-on docs passes for portability and discoverability cleanup.
