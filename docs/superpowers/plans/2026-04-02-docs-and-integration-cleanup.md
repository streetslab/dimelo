# Docs And Integration Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Normalize machine-specific paths in tracked historical `docs/superpowers` docs, then perform a narrow discoverability-focused integration cleanup around the already-implemented analysis and plotting work.

**Architecture:** Execute this cleanup in two phases. Phase A normalizes tracked spec/plan docs so they are portable in the repo without changing their meaning. Phase B applies only the smallest public-surface cleanup needed to make the recently completed work easier to find and understand. The central status source of truth remains [docs/superpowers/README.md](../README.md).

**Tech Stack:** Markdown docs, existing repo documentation, lightweight Python/package surface inspection where needed

---

## File Structure

- Modify: `docs/superpowers/specs/`
  - Normalize machine-local path references to repo-relative or plain repo paths.
- Modify: `docs/superpowers/plans/`
  - Normalize machine-local path references to repo-relative or plain repo paths.
- Modify: `docs/superpowers/README.md`
  - Only if link/status corrections are needed after normalization.
- Modify if needed for integration cleanup:
  - `README.md`
  - `docs/region-contrasts.md`
  - `docs/shared-clustering.md`
  - `docs/region-discovery.md`
  - `docs/global-analysis.md`
  - `dimelo/__init__.py`

### Task 1: Normalize machine-specific paths in tracked `docs/superpowers` docs

**Files:**
- Modify: `docs/superpowers/specs/*.md`
- Modify: `docs/superpowers/plans/*.md`

- [ ] **Step 1: Inventory the remaining machine-local path references**

Run:

```bash
rg -n "/Users/[^/]+/Documents/GitHub/dimelo-toolkit" docs/superpowers/specs docs/superpowers/plans
```

Expected: a list of tracked historical docs that still contain machine-local absolute paths.

- [ ] **Step 2: Normalize Markdown links to repo-relative links**

For Markdown links inside tracked superpowers docs, convert forms like:

```md
[plotting.py](../../../dimelo/plotting.py)
```

to:

```md
[plotting.py](../../../dimelo/plotting.py)
```

or the correct relative form from the document location.

- [ ] **Step 3: Normalize plain path mentions and command snippets**

For plain path references or shell commands inside tracked superpowers docs, replace machine-local prefixes with repo-relative or repo-root-style paths. For example:

```text
docs/superpowers/plans/2026-04-01-region-contrasts-plotting.md
```

becomes:

```text
docs/superpowers/plans/2026-04-01-region-contrasts-plotting.md
```

Do this conservatively: preserve the original meaning and command intent, only remove machine-specific leakage.

- [ ] **Step 4: Run a second grep to confirm the path leak is gone**

Run:

```bash
rg -n "/Users/[^/]+/Documents/GitHub/dimelo-toolkit" docs/superpowers/specs docs/superpowers/plans
```

Expected: no output, unless an intentionally preserved example remains and is documented in the final summary.

- [ ] **Step 5: Commit the normalization pass**

```bash
git add docs/superpowers/specs docs/superpowers/plans
git commit -m "docs: normalize superpowers doc paths"
```

### Task 2: Verify and repair `docs/superpowers` navigation after normalization

**Files:**
- Modify if needed: `docs/superpowers/README.md`
- Verify: `docs/superpowers/specs/`
- Verify: `docs/superpowers/plans/`

- [ ] **Step 1: Verify the central index links still resolve**

Run:

```bash
for f in $(rg -o "\\((specs|plans)/[^)]+" docs/superpowers/README.md | tr -d "("); do test -f "docs/superpowers/$f" || echo "MISSING $f"; done
```

Expected: no output.

- [ ] **Step 2: Fix any broken links or stale status wording in the index**

If the verification step reports missing links or if the normalization changed path conventions enough to make the README confusing, update `docs/superpowers/README.md` minimally.

Keep edits narrow:

- link correction
- wording clarification
- status correction if a recent implementation slice is now clearly misclassified

- [ ] **Step 3: Verify the `docs/superpowers` tree is fully tracked**

Run:

```bash
git status --short docs/superpowers
```

Expected: no `??` entries remain under `docs/superpowers`.

- [ ] **Step 4: Commit the index repair pass if needed**

```bash
git add docs/superpowers/README.md
git commit -m "docs: refresh superpowers index links"
```

If no README changes were needed, skip this commit.

### Task 3: Apply narrow discoverability-focused integration cleanup

**Files:**
- Inspect: `README.md`
- Inspect: `docs/region-contrasts.md`
- Inspect: `docs/shared-clustering.md`
- Inspect: `docs/region-discovery.md`
- Inspect: `docs/global-analysis.md`
- Inspect: `dimelo/__init__.py`
- Modify only if needed: the minimum subset of the above

- [ ] **Step 1: Inventory obvious discoverability gaps around completed work**

Check for:

- recently added helpers/workflows not mentioned where users would look
- stale references to “later” work that is now implemented
- missing package exports for helpers users reasonably expect to import directly

Run these focused inspections:

```bash
rg -n "prepare_region_contrast_profile_data|prepare_region_contrast_heatmap_data|global_analysis|region_discovery|shared_cluster_distribution" README.md docs/*.md dimelo/__init__.py
```

- [ ] **Step 2: Apply only the minimum integration cleanup needed**

Examples of acceptable changes:

- add one missing cross-link in `README.md`
- fix stale wording in an analysis guide
- expose an already-implemented helper in `dimelo/__init__.py` if that is an obvious user-facing expectation

Examples of disallowed changes:

- broad module refactors
- new plotting features
- renaming stable APIs for style alone

- [ ] **Step 3: Verify touched docs or exports**

If only docs changed, run:

```bash
rg -n "prepare_region_contrast_profile_data|prepare_region_contrast_heatmap_data|shared_cluster_distribution|global_analysis|region_discovery" README.md docs/*.md
```

If `dimelo/__init__.py` changed, also run:

```bash
python3.11 - <<'PY'
import dimelo
print("dimelo import ok")
PY
```

Expected: docs references resolve sensibly, and `dimelo` still imports if exports changed.

- [ ] **Step 4: Commit the integration cleanup**

```bash
git add README.md docs dimelo/__init__.py
git commit -m "docs: tighten analysis discoverability"
```

Only include files actually changed.

### Task 4: Final verification and cleanup summary

**Files:**
- Verify: `docs/superpowers/`
- Verify any files touched in Task 3

- [ ] **Step 1: Confirm no machine-local path leaks remain in tracked superpowers docs**

Run:

```bash
rg -n "/Users/[^/]+/Documents/GitHub/dimelo-toolkit" docs/superpowers/specs docs/superpowers/plans docs/superpowers/README.md
```

Expected: no output unless an intentional exception is documented.

- [ ] **Step 2: Confirm the working tree is clean**

Run:

```bash
git status --short
```

Expected: no unexpected modified/untracked files remain from this cleanup slice.

- [ ] **Step 3: Commit any final small verification-driven fix**

```bash
git add -A
git commit -m "docs: finalize cleanup verification"
```

Only do this if verification required one last small correction. Otherwise skip.

## Self-Review

- Spec coverage:
  - historical path portability: Task 1
  - index/navigation verification: Task 2
  - public-surface discoverability cleanup: Task 3
  - final verification: Task 4
- Placeholder scan:
  - no `TODO`, `TBD`, or “fix later” steps remain
- Type consistency:
  - all path references use the current `docs/superpowers` layout and the current repo root
