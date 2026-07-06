# Shared Cluster Tests Completion And Legacy Cleanup Execution Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finish the active shared-cluster inference slice, then execute a prioritized legacy TODO cleanup program without derailing current delivery.

**Architecture:** Treat shared-cluster inference completion as Sprint 0 (feature lane), then run legacy cleanup in Sprints 1-5 (debt lane). Keep commits narrow and test-gated. Prefer no cross-sprint mixed commits.

**Tech Stack:** Python 3.11, pytest, pandas, numpy, scipy, dimelo workflows/models/plotting stack

---

## Status Note (2026-04-14)

This document is retained as the historical sprint breakdown. The canonical
"what remains" backlog is now tracked in
`docs/superpowers/plans/2026-04-14-unified-remaining-work-plan.md`.

## How Current Track Fits

- Sprint 0 is the current active branch work:
  - `dimelo/shared_cluster_tests.py`
  - `tests/test_shared_cluster_tests.py`
  - `README.md`
  - `docs/shared-clustering.md`
  - `docs/superpowers/README.md`
- Sprints 1-5 are legacy debt cleanup outside the new analysis track.
- Recommended flow: complete Sprint 0 first, then run Sprints 1-5 in order.

## File Map

- Modify: `dimelo/shared_cluster_tests.py`
- Modify: `tests/test_shared_cluster_tests.py`
- Modify: `docs/shared-clustering.md`
- Modify: `README.md`
- Modify: `docs/superpowers/README.md`
- Modify: `dimelo/parse_bam.py`
- Modify: `dimelo/run_modkit.py`
- Modify: `dimelo/load_processed.py`
- Modify: `dimelo/plot_read_browser.py`
- Modify: `dimelo/plot_reads.py`
- Modify: `dimelo/plot_enrichment.py`
- Modify: `dimelo/plot_enrichment_profile.py`
- Modify: `dimelo/plot_depth_profile.py`
- Modify: `dimelo/plot_depth_histogram.py`
- Modify: `dimelo/utils.py`
- Modify: `dimelo/test_data.py`
- Modify: `dimelo/test/__init__.py`
- Modify: `dimelo/test/dimelo_test.py`

## Sprint 0: Complete Shared Cluster Tests (Current Track)

**Files:**
- Modify: `tests/test_shared_cluster_tests.py`
- Modify: `dimelo/shared_cluster_tests.py`
- Modify: `docs/shared-clustering.md`
- Modify: `README.md`
- Modify: `docs/superpowers/README.md`

- [ ] **Step 1: Add failing tests for missing scope (time-course + pooled tests + include_pairwise)**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "time_course or chi_squared or g_test or include_pairwise" -q
```

Expected: FAIL on currently missing behavior.

- [ ] **Step 2: Implement pooled screening paths (`chi_squared`, `g_test`)**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "chi_squared or g_test" -q
```

Expected: PASS for pooled tests.

- [ ] **Step 3: Implement `time_course` and `include_pairwise` follow-up**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py -k "time_course or include_pairwise" -q
```

Expected: PASS for time-course behavior.

- [ ] **Step 4: Update docs and status map**

Update:
- `docs/shared-clustering.md` (global shared-cluster tests usage and boundaries)
- `README.md` (analysis guide note for `shared_cluster_tests`)
- `docs/superpowers/README.md` (status line updates)

- [ ] **Step 5: Run shared-cluster regression subset**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_shared_cluster_tests.py tests/test_models.py tests/test_workflows.py -k "shared_cluster" -q
```

Expected: PASS.

- [ ] **Step 6: Commit Sprint 0**

```bash
git add dimelo/shared_cluster_tests.py tests/test_shared_cluster_tests.py docs/shared-clustering.md README.md docs/superpowers/README.md
git commit -m "feat: complete shared cluster inference v1 scope"
```

## Sprint 1: Parser Correctness And Duplication

**Files:**
- Modify: `dimelo/parse_bam.py`
- Modify: `dimelo/run_modkit.py`
- Test: `tests/test_distribution.py`
- Test: `tests/test_cluster.py`

- [ ] **Step 1: Identify parser TODO targets and lock scope**

Run:

```bash
rg -n "TODO|FIXME|TBD" dimelo/parse_bam.py dimelo/run_modkit.py
```

Expected: Enumerated TODO set for this sprint only.

- [ ] **Step 2: Refactor duplicate pileup/extract setup paths**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_distribution.py -q
```

Expected: PASS with no behavior regression.

- [ ] **Step 3: Harden coordinate/chunk/compression paths**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_cluster.py tests/test_distribution.py -q
```

Expected: PASS.

- [ ] **Step 4: Clean `run_modkit.py` typing TODO and run focused regression**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_distribution.py tests/test_cluster.py -q
```

Expected: PASS.

- [ ] **Step 5: Commit Sprint 1**

```bash
git add dimelo/parse_bam.py dimelo/run_modkit.py tests/test_distribution.py tests/test_cluster.py
git commit -m "refactor: harden parse_bam core paths and reduce duplication"
```

## Sprint 2: Processed Loader Reliability

**Files:**
- Modify: `dimelo/load_processed.py`
- Test: `tests/test_distribution.py`
- Test: `tests/test_workflows.py`

- [ ] **Step 1: Triage loader TODOs and group by behavior**

Run:

```bash
rg -n "TODO|FIXME|TBD" dimelo/load_processed.py
```

Expected: TODO set grouped into naming, progress/parallelization, subsetting semantics.

- [ ] **Step 2: Clean naming/API clarity and keep compatibility**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_distribution.py -q
```

Expected: PASS.

- [ ] **Step 3: Add progress/parallel behavior parity and subsetting fixes**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_distribution.py tests/test_workflows.py -k "load or processed" -q
```

Expected: PASS.

- [ ] **Step 4: Commit Sprint 2**

```bash
git add dimelo/load_processed.py tests/test_distribution.py tests/test_workflows.py
git commit -m "fix: improve load_processed reliability and API clarity"
```

## Sprint 3: Read Browser Behavior And Performance

**Files:**
- Modify: `dimelo/plot_read_browser.py`
- Test: `tests/test_plotting.py`

- [ ] **Step 1: Resolve collapse/sort semantics TODOs**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "read_browser" -q
```

Expected: PASS.

- [ ] **Step 2: Stabilize duplicate-read handling and hover contract**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "read_browser or plot_reads" -q
```

Expected: PASS.

- [ ] **Step 3: Reduce hot-path iteration overhead where safe**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

Expected: PASS.

- [ ] **Step 4: Commit Sprint 3**

```bash
git add dimelo/plot_read_browser.py tests/test_plotting.py
git commit -m "fix: harden read browser behavior and clean plotting debt"
```

## Sprint 4: Legacy Plotting API Normalization

**Files:**
- Modify: `dimelo/plot_reads.py`
- Modify: `dimelo/plot_enrichment.py`
- Modify: `dimelo/plot_enrichment_profile.py`
- Modify: `dimelo/plot_depth_profile.py`
- Modify: `dimelo/plot_depth_histogram.py`
- Test: `tests/test_plotting.py`

- [ ] **Step 1: Normalize dispatch/API patterns across legacy plotting modules**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -k "plot_reads or enrichment or depth" -q
```

Expected: PASS.

- [ ] **Step 2: Remove redefinition/mypy TODO clusters without changing output behavior**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_plotting.py -q
```

Expected: PASS.

- [ ] **Step 3: Commit Sprint 4**

```bash
git add dimelo/plot_reads.py dimelo/plot_enrichment.py dimelo/plot_enrichment_profile.py dimelo/plot_depth_profile.py dimelo/plot_depth_histogram.py tests/test_plotting.py
git commit -m "refactor: normalize legacy plotting APIs and remove type-debt hotspots"
```

## Sprint 5: Utils And Fixture Hygiene

**Files:**
- Modify: `dimelo/utils.py`
- Modify: `dimelo/test_data.py`
- Modify: `dimelo/test/__init__.py`
- Modify: `dimelo/test/dimelo_test.py`

- [ ] **Step 1: Resolve utility reproducibility/type TODOs**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest tests/test_distribution.py tests/test_cluster.py -q
```

Expected: PASS.

- [ ] **Step 2: Improve fixture realism and remove internal test TODOs**

Run:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest -q
```

Expected: PASS.

- [ ] **Step 3: Commit Sprint 5**

```bash
git add dimelo/utils.py dimelo/test_data.py dimelo/test/__init__.py dimelo/test/dimelo_test.py
git commit -m "chore: improve utility typing/reproducibility and test fixtures"
```

## Cross-Sprint Gates

- [ ] `git status --short` is clean before each sprint commit.
- [ ] No mixed-sprint commits (one sprint theme per commit).
- [ ] No new `TODO/FIXME/TBD` markers introduced.
- [ ] Full regression at sprint boundaries:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.11 -m pytest -q
```

## Success Metrics

- Baseline markers: `96` TODO/FIXME/TBD in `dimelo/*.py`.
- Post Sprint 0 target: current feature lane complete and documented.
- Post Sprint 3 target: TODO count < `35`.
- Post Sprint 5 target: TODO count < `10`, with any remainder explicitly deferred in docs.
