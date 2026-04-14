# Docs Coherence And Discoverability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the internal and user-facing docs coherent with the current branch state, and apply only minimal package-surface cleanup where the docs clearly demand it.

**Architecture:** Execute this as a docs-first cleanup in three layers: refresh the internal status map, sweep the main user-facing guides for naming and flow consistency, then inspect the top-level package surface for any small discoverability mismatch exposed by the docs pass. Keep the slice additive and non-behavioral.

**Tech Stack:** Markdown docs, existing Python package surface, lightweight repo inspection

---

## Status Note (2026-04-14)

This document remains the detailed execution playbook for docs coherence work.
The active cross-plan backlog now lives in
`docs/superpowers/plans/2026-04-14-unified-remaining-work-plan.md`.

## File Map

- Modify: `docs/superpowers/README.md`
  - refresh status-map entries and add missing spec/plan links
- Modify as needed:
  - `README.md`
  - `docs/global-analysis.md`
  - `docs/region-discovery.md`
  - `docs/region-contrasts.md`
  - `docs/shared-clustering.md`
  - only for wording, cross-link, and role-consistency cleanup
- Modify only if clearly needed:
  - `dimelo/__init__.py`
  - only for a small discoverability/export mismatch exposed by the docs

## Task 1: Refresh The Internal Status Map

**Files:**
- Modify: `docs/superpowers/README.md`
- Verify: `docs/superpowers/specs/`
- Verify: `docs/superpowers/plans/`

- [ ] **Step 1: Inventory missing or stale entries**

Inspect the current index and repo contents:

```bash
sed -n '1,260p' docs/superpowers/README.md
ls docs/superpowers/specs docs/superpowers/plans
```

Expected findings:

- the internal index is missing the newest plotting specs/plans
- status labels may understate work that is already implemented

- [ ] **Step 2: Update the status map minimally**

Edit `docs/superpowers/README.md` to:

- add:
  - `specs/2026-04-02-region-discovery-plotting-design.md`
  - `specs/2026-04-02-global-analysis-plotting-design.md`
  - `specs/2026-04-02-docs-coherence-and-discoverability-design.md`
  - `plans/2026-04-02-region-discovery-plotting.md`
  - `plans/2026-04-02-global-analysis-plotting.md`
- refresh status labels where current branch reality clearly differs
- update the `Current Themes` section only as needed to reflect the newer plotting helper coverage

Keep the file concise. Do not turn it into a changelog.

- [ ] **Step 3: Verify index links resolve**

Run:

```bash
for f in $(rg -o "\\((specs|plans)/[^)]+" docs/superpowers/README.md | tr -d "("); do test -f "docs/superpowers/$f" || echo "MISSING $f"; done
```

Expected: no output.

- [ ] **Step 4: Commit the status-map refresh**

```bash
git add docs/superpowers/README.md
git commit -m "docs: refresh superpowers status map"
```

## Task 2: Sweep User-Facing Guides For Consistency

**Files:**
- Modify as needed:
  - `README.md`
  - `docs/global-analysis.md`
  - `docs/region-discovery.md`
  - `docs/region-contrasts.md`
  - `docs/shared-clustering.md`

- [ ] **Step 1: Inventory stale wording and naming drift**

Run focused searches:

```bash
rg -n "once implemented|later for|not yet|planned|prepare_region_contrast_|prepare_region_discovery_|prepare_global_analysis_|shared_cluster_distribution|discovery_cluster_workflow|discovery_cluster_contrast_workflow" README.md docs/*.md
```

Inspect the matching sections and note only concrete consistency issues such as:

- stale “future” wording
- inconsistent helper names
- inconsistent module-role descriptions
- missing short cross-links where users would naturally look

- [ ] **Step 2: Apply a narrow consistency sweep**

Edit only the minimum subset of the user-facing guides needed so they:

- tell one consistent story about the analysis stack
- use the same helper/workflow names everywhere
- describe plotting-helper coverage honestly
- avoid implying that already-shipped helpers or workflows are still future work

Allowed changes:

- wording cleanup
- short cross-link additions
- short clarifications of current helper coverage

Disallowed changes:

- major rewrites
- new feature documentation beyond what is already implemented
- API renames for style

- [ ] **Step 3: Verify the guide sweep**

Run:

```bash
rg -n "prepare_region_contrast_|prepare_region_discovery_|prepare_global_analysis_|shared_cluster_distribution|discovery_cluster_workflow|discovery_cluster_contrast_workflow" README.md docs/*.md
```

Expected: matches should reflect the intended current helper/workflow names consistently.

- [ ] **Step 4: Commit the docs consistency sweep**

```bash
git add README.md docs/global-analysis.md docs/region-discovery.md docs/region-contrasts.md docs/shared-clustering.md
git commit -m "docs: align analysis guide naming"
```

Only include files actually changed.

## Task 3: Apply Small Package-Surface Cleanup If Needed

**Files:**
- Inspect: `dimelo/__init__.py`
- Modify only if clearly needed: `dimelo/__init__.py`

- [ ] **Step 1: Check whether the docs expose an import mismatch**

Inspect the top-level package surface:

```bash
sed -n '1,220p' dimelo/__init__.py
```

Compare that with the docs updated in Task 2.

Only proceed if there is a concrete mismatch such as:

- docs telling users to import a top-level module that is not actually exposed
- docs implying a clearly user-facing module path that the package does not support

- [ ] **Step 2: Make the minimum export fix if needed**

If and only if a concrete mismatch exists, update `dimelo/__init__.py` minimally.

Allowed:

- add one missing module export already treated as top-level in the docs

Disallowed:

- broad re-export churn
- exporting helper functions just for convenience
- aesthetic reshaping of `__all__`

- [ ] **Step 3: Verify imports if `__init__.py` changed**

Run:

```bash
python3.11 - <<'PY'
import dimelo
print("dimelo import ok")
PY
```

Expected: `dimelo import ok`

- [ ] **Step 4: Commit the package-surface cleanup if needed**

```bash
git add dimelo/__init__.py
git commit -m "docs: align package surface with guides"
```

If no package-surface fix was needed, skip this commit.

## Final Verification

**Files:**
- Verify: `docs/superpowers/README.md`
- Verify any touched user-facing docs
- Verify `dimelo/__init__.py` if changed

- [ ] **Step 1: Confirm the internal index resolves**

Run:

```bash
for f in $(rg -o "\\((specs|plans)/[^)]+" docs/superpowers/README.md | tr -d "("); do test -f "docs/superpowers/$f" || echo "MISSING $f"; done
```

Expected: no output.

- [ ] **Step 2: Confirm the working tree is clean**

Run:

```bash
git status --short
```

Expected: no output.

- [ ] **Step 3: Summarize the coherence outcomes**

Capture in the handoff:

- which status-map entries were added or corrected
- which user-facing docs were adjusted
- whether any package-surface change was needed
- the main remaining next-step categories after docs are coherent

## Self-Review

Spec coverage check:

- internal status-map accuracy: covered in Task 1
- user-facing guide consistency: covered in Task 2
- small package-surface cleanup only if needed: covered in Task 3
- validation and clean-state closeout: covered in Final Verification

Placeholder scan:

- no `TODO`, `TBD`, or undefined steps remain
- each task includes exact files, commands, and narrow allowed scope

Type consistency:

- all referenced helper/workflow names match the current approved naming:
  - `shared_cluster_distribution(...)`
  - `discovery_cluster_workflow(...)`
  - `discovery_cluster_contrast_workflow(...)`
  - `prepare_region_contrast_*`
  - `prepare_region_discovery_*`
  - `prepare_global_analysis_*`
