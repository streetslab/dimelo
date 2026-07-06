# Region Discovery Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a first `region_discovery` layer for deterministic de novo locus finding from pileup-backed inputs, with clean handoff into downstream known-region workflows.

**Architecture:** Build `dimelo.region_discovery` on top of existing parsing + loader behavior (`parse_bam.pileup`, `load_processed.pileup_counts_from_bedmethyl`, `load_processed.regions_to_list`) and the new `global_analysis` windowing utilities. Keep v1 deterministic and tiled.

**Tech Stack:** Python 3.11, pandas, numpy, pytest

---

## File Map

- `dimelo/models.py`
  Add `RegionDiscoveryResult` as the canonical data-first result container.
- `dimelo/region_discovery.py`
  New module for scanning, ranking, hit merging, and BED export.
- `dimelo/__init__.py`
  Export `region_discovery`.
- `tests/test_models.py`
  Add model-level validation for `RegionDiscoveryResult`.
- `tests/test_region_discovery.py`
  Add focused tests for scan, filtering, ranking, merge, and export helpers.
- `docs/region-discovery.md`
  Add user-facing guide for de novo discovery and handoff.
- `README.md`
  Link the new guide from the analysis guides section.

Scope note:

- This plan covers `region_discovery` foundations only.
- It does not implement `region_contrasts` inferential expansions.
- It does not modify existing parse/load semantics.
- It keeps outputs data-first and plotting-optional.

### Task 1: Add RegionDiscoveryResult Model

**Files:**
- Modify: `dimelo/models.py`
- Modify: `tests/test_models.py`

- [ ] **Step 1: Write failing tests**
  - Add a test that constructs `RegionDiscoveryResult` with `hits`, `windows`, `plot_data`, and metadata.
  - Add a test that rejects `None` for required fields (`hits`, `windows`, `plot_data`).

- [ ] **Step 2: Verify tests fail**
  - Run:
    - `pytest tests/test_models.py::test_region_discovery_result_accepts_tables_and_plot_data tests/test_models.py::test_region_discovery_result_rejects_missing_required_fields -q`

- [ ] **Step 3: Implement minimal model**
  - Add dataclass:
    - `hits: pd.DataFrame`
    - `windows: pd.DataFrame`
    - `contrast: ContrastSpec | None`
    - `plot_data: dict[str, pd.DataFrame | dict[str, Any]]`
    - `metadata: dict[str, Any] = field(default_factory=dict)`
    - `figures: dict[str, Any] = field(default_factory=dict)`
  - Add `__post_init__` required field checks.

- [ ] **Step 4: Verify tests pass**
  - Re-run the focused tests above.

- [ ] **Step 5: Commit**
  - `git add dimelo/models.py tests/test_models.py`
  - `git commit -m "feat: add region discovery result model"`

### Task 2: Add Deterministic Genome Scan And Window Scoring

**Files:**
- Create: `dimelo/region_discovery.py`
- Modify: `dimelo/__init__.py`
- Create: `tests/test_region_discovery.py`

- [ ] **Step 1: Write failing tests**
  - Add tests for:
    - `scan_genome(...)` basic behavior with a mocked window summary table.
    - `scan_genome(...)` minimum coverage filtering.
    - `scan_genome(...)` contig include/exclude handoff to tiled window builder.
    - significance ranking fields (`score_value`, `p_value`, `adjusted_p_value`, `rank`).

- [ ] **Step 2: Verify tests fail**
  - Run the focused test subset from `tests/test_region_discovery.py`.

- [ ] **Step 3: Implement minimal discovery scan**
  - Add in `dimelo/region_discovery.py`:
    - `scan_genome(...)`
      - uses `global_analysis.build_window_summary(...)`
      - requires exactly one motif in v1
      - supports `window_size`, `step_size`, `include_contigs`, `exclude_contigs`
      - applies `min_coverage` filter (`valid_count >= min_coverage`)
      - computes a simple v1 discovery score:
        - default `score="effect_size_only"`: rank by window fraction across contrast groups
        - support `score="beta_binomial"` pooled approximation for matched group comparisons
      - computes BH-adjusted p-values when p-values are present
      - returns a `RegionDiscoveryResult` with:
        - `windows`: per-window evidence table
        - `hits`: ranked significant windows
        - `plot_data`: at least `window_score_table`, `top_hits_table`
        - metadata with scan policy and thresholds
  - Keep deterministic ordering with stable sort keys.

- [ ] **Step 4: Verify tests pass**
  - Run:
    - `pytest tests/test_region_discovery.py -q`

- [ ] **Step 5: Commit**
  - `git add dimelo/region_discovery.py dimelo/__init__.py tests/test_region_discovery.py`
  - `git commit -m "feat: add deterministic region discovery scan"`

### Task 3: Add Hit Merging And BED Export Helpers

**Files:**
- Modify: `dimelo/region_discovery.py`
- Modify: `tests/test_region_discovery.py`

- [ ] **Step 1: Write failing tests**
  - Add tests for:
    - `merge_adjacent_hits(...)` merges neighboring significant windows by chromosome and merge distance.
    - `merge_adjacent_hits(...)` preserves score/rank summary fields deterministically.
    - `hits_to_bed(...)` exports required BED columns and ordering.

- [ ] **Step 2: Verify tests fail**
  - Run focused tests for merge/export helpers.

- [ ] **Step 3: Implement minimal helpers**
  - Add:
    - `merge_adjacent_hits(hits: pd.DataFrame, merge_distance: int) -> pd.DataFrame`
    - `hits_to_bed(hits: pd.DataFrame) -> pd.DataFrame`
  - Ensure merge happens after per-window scoring/FDR as per spec policy.

- [ ] **Step 4: Verify tests pass**
  - Run:
    - `pytest tests/test_region_discovery.py -q`

- [ ] **Step 5: Commit**
  - `git add dimelo/region_discovery.py tests/test_region_discovery.py`
  - `git commit -m "feat: add region discovery merge and bed export helpers"`

### Task 4: Add User Guide And End-to-End Result Wrapper Polish

**Files:**
- Modify: `dimelo/region_discovery.py`
- Modify: `tests/test_region_discovery.py`
- Create: `docs/region-discovery.md`
- Modify: `README.md`

- [ ] **Step 1: Write failing tests**
  - Add/extend a test to verify `scan_genome(...)` returns `RegionDiscoveryResult` with expected `plot_data` keys and metadata policy fields.

- [ ] **Step 2: Verify tests fail**
  - Run the focused test(s) first.

- [ ] **Step 3: Implement wrapper polish and docs**
  - Ensure `scan_genome(...)` includes metadata fields for:
    - `analysis_unit`
    - `representation`
    - `signal_source`
    - `score`
    - `window_size`
    - `step_size`
    - `min_coverage`
    - `merge_hits`
    - `merge_distance`
  - Add `docs/region-discovery.md` with:
    - when to use discovery vs contrasts
    - minimal usage example
    - handoff to `region_contrasts`/clustering
  - Link the guide from README analysis guides.

- [ ] **Step 4: Verify tests pass**
  - Run:
    - `pytest tests/test_region_discovery.py tests/test_models.py -q`

- [ ] **Step 5: Commit**
  - `git add dimelo/region_discovery.py tests/test_region_discovery.py docs/region-discovery.md README.md`
  - `git commit -m "docs: add region discovery workflow guide"`

---

## Verification And Review Gates

- After each task:
  - Run the task-focused pytest subset.
  - Run spec-compliance review.
  - Run code-quality review.
- After Task 4:
  - Run:
    - `pytest tests/test_region_discovery.py tests/test_models.py tests/test_global_analysis.py tests/test_region_contrasts.py tests/test_workflows.py -q`

## V1 Defaults To Preserve

- deterministic tiled-window discovery
- canonical contigs by default with include/exclude override
- minimum coverage threshold before scoring
- FDR on raw windows before merge
- merge significant windows afterward with fixed merge distance
- data-first result tables with optional Matplotlib-agnostic plot payloads
- additive architecture around existing parsing/load modules
