# Unified Remaining Work Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans. This plan consolidates outstanding work across historical `docs/superpowers/plans/*.md`.

**Goal:** Maintain one canonical backlog of what is still open after auditing all historical plan docs, and execute remaining work in a clean, commit-friendly sequence.

**Audit Date:** 2026-04-14

## Audited Plan Inventory

The table below captures what remains from every tracked plan document.

| Plan | Status | Remaining Work |
|---|---|---|
| `2026-03-31-global-analysis-foundations.md` | implemented | none |
| `2026-03-31-region-contrasts-foundations.md` | implemented | none |
| `2026-03-31-region-discovery-foundations.md` | implemented | none |
| `2026-03-31-shared-clustering-foundations.md` | implemented | none |
| `2026-04-01-cluster-occupancy-region-contrasts.md` | implemented | none |
| `2026-04-01-discovery-cluster-contrast-workflow.md` | implemented | none |
| `2026-04-01-paired-region-discovery.md` | implemented | none |
| `2026-04-01-plotting-axis-architecture.md` | implemented | none |
| `2026-04-01-region-contrasts-plotting.md` | implemented | none |
| `2026-04-01-region-discovery-cluster-workflow.md` | implemented | none |
| `2026-04-02-docs-and-integration-cleanup.md` | implemented | none |
| `2026-04-02-docs-coherence-and-discoverability.md` | partially implemented | keep status index synchronized with newly added plans and run periodic coherence checks |
| `2026-04-02-global-analysis-plotting.md` | implemented | none |
| `2026-04-02-region-discovery-plotting.md` | implemented | none |
| `2026-04-02-shared-clustering-plotting.md` | implemented | none |
| `2026-04-02-superpowers-docs-cleanup.md` | implemented | none |
| `2026-04-03-single-read-region-contrasts.md` | implemented | none |
| `2026-04-08-matplotlib-renderers.md` | implemented | none |
| `2026-04-10-shared-clustering-matplotlib-renderers.md` | implemented | none |
| `2026-04-12-shared-cluster-tests.md` | implemented | none |
| `2026-04-13-legacy-cleanup-and-shared-cluster-tests-execution.md` | partially implemented | finish legacy cleanup sprints and close commit-hygiene gates |

## Canonical Remaining Backlog

Only two workstreams remain active after the full audit:

1. Legacy cleanup execution (from `2026-04-13-*`)
2. Ongoing docs coherence/index maintenance (from `2026-04-02-docs-coherence-*`)

## Workstream A: Legacy Cleanup Completion

**Scope Source:** `2026-04-13-legacy-cleanup-and-shared-cluster-tests-execution.md`

### A1. Finish Sprint 1 leftovers (parser/runtime debt)

- [ ] `run_modkit.py`: resolve remaining typing and API clarity debt.
- [ ] `parse_bam.py`: close remaining high-value TODO clusters that are still behavior-adjacent, not cosmetic.
- [ ] Verification:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_distribution.py tests/test_cluster.py -q
```

### A2. Finish Sprint 4 (legacy plotting normalization)

- [ ] Normalize dispatch/API patterns in:
  - `dimelo/plot_reads.py`
  - `dimelo/plot_enrichment.py`
  - `dimelo/plot_enrichment_profile.py`
  - `dimelo/plot_depth_profile.py`
  - `dimelo/plot_depth_histogram.py`
- [ ] Remove redefinition/mypy-TODO hotspots without output behavior drift.
- [ ] Verification:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

### A3. Finish Sprint 5 (utils + legacy fixture hygiene)

- [ ] `utils.py`: close reproducibility/type TODOs that affect public behavior or determinism.
- [ ] `dimelo/test_data.py`, `dimelo/test/__init__.py`, `dimelo/test/dimelo_test.py`:
  - reduce external fragility (network/download assumptions),
  - improve restricted-environment compatibility (process/semaphore constraints),
  - keep legacy integration coverage meaningful.
- [ ] Verification:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_distribution.py tests/test_cluster.py -q
```

## Workstream B: Docs Coherence Maintenance

**Scope Source:** `2026-04-02-docs-coherence-and-discoverability.md`

- [ ] Keep `docs/superpowers/README.md` in sync with all files in `docs/superpowers/plans/`.
- [ ] Keep status labels aligned with branch reality for the active cycle.
- [ ] Run link-existence verification after index edits:

```bash
for f in $(rg -o "\\((specs|plans)/[^)]+" docs/superpowers/README.md | tr -d "("); do
  test -f "docs/superpowers/$f" || echo "MISSING $f"
done
```

## Execution Order Going Forward

1. Finish Workstream A1, commit.
2. Finish Workstream A2, commit.
3. Finish Workstream A3, commit.
4. Apply Workstream B coherence pass, commit.
5. Run final verification suite for the branch:

```bash
pytest tests
```

## Commit Hygiene Rules (Canonical)

- [ ] Keep one sprint/theme per commit.
- [ ] Keep `git status --short` clean between sprint commits.
- [ ] Do not introduce new `TODO/FIXME/TBD` markers in touched files.
- [ ] Record verification command output in PR notes before merge.

