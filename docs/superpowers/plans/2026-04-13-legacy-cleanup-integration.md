# Legacy Cleanup Integration

**Date:** 2026-04-13  
**Branch:** `codex/legacy-cleanup`  
**Range for this cleanup track:** `8baefd2..HEAD`

## Goal

Package the parser/loader legacy-cleanup sprint into reviewable PR slices with clear risk boundaries and a repeatable verification matrix.

## Scope Summary

Track-level diff (`8baefd2..HEAD`):

- `13 files changed`
- `1755 insertions(+), 697 deletions(-)`
- Core files touched:
  - `dimelo/parse_bam.py`
  - `dimelo/load_processed.py`
  - `dimelo/run_modkit.py`
  - `dimelo/utils.py`
  - `dimelo/plot_*.py` (debt-only cleanup)
  - parser/loader-focused tests

## What Landed

### Parser lane

- Clarified and stabilized threshold semantics in extract conversion.
- Fixed extract coordinate reconstruction and added explicit bounds validation.
- Consolidated chunk/dataset write flow and CSV row parsing in `read_by_base_txt_to_hdf5`.
- Simplified parser helper contracts (`prep_output_directory`, path coercion, thread-command handling).
- Removed stale TODO debt and repeated cleanup code paths.

### Loader lane

- Hardened subset contract validation (`dict` only, must include `n`/`frac`, reject `array`).
- Fixed empty-region behavior and empty-subset safety.
- Disambiguated duplicate read names across overlapping regions.
- Added optional parallel paths for region selection and readwise processing.
- Added loader progress wiring and public alias helpers:
  - `counts_from_pileup`
  - `vectors_from_pileup`

### Plotting/util lane

- Cleanup-only removal of stale TODO/type-reassignment debt in plotting helpers.
- No intentional behavior changes in plotting outputs.

### Hardening tests lane

- Added parser guardrail tests:
  - invalid strand row raises
  - non-positive read length raises
  - out-of-bounds `pos_in_read` raises
- Added loader contract tests:
  - `subset_parameters` non-dict rejection
  - tuple-form `sort_by` path normalization/ordering

## Verification Matrix

All commands below passed on `codex/legacy-cleanup` after lane completion:

1. `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /Users/ngamarra/.pyenv/versions/3.11.3/bin/python3.11 -m pytest tests/test_parse_bam.py tests/test_parse_bam_vectors.py tests/test_utils_regions.py -q`  
   `25 passed`
2. `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /Users/ngamarra/.pyenv/versions/3.11.3/bin/python3.11 -m pytest tests/test_distribution.py tests/test_cluster.py -q`  
   `20 passed, 2 warnings` (existing `tight_layout` and joblib core-detection warnings)
3. `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /Users/ngamarra/.pyenv/versions/3.11.3/bin/python3.11 -m pytest tests/test_region_contrasts.py -k "read_mod_fraction" -q`  
   `19 passed, 65 deselected`
4. `PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /Users/ngamarra/.pyenv/versions/3.11.3/bin/python3.11 -m pytest tests/test_workflows.py -k "read_mod_fraction or shared_cluster" -q`  
   `14 passed, 18 deselected, 1 warning` (existing joblib core-detection warning)

## PR Split Recommendation

Split this cleanup track into **3 thematic PRs** (recommended) to keep review focused:

1. **PR A: Parser correctness + parser refactor**
   - `parse_bam.py`, `run_modkit.py`, parser tests, `utils.py` parser-adjacent fixes.
   - Include coordinate and threshold semantics changes with their tests.
2. **PR B: Loader contract + parallelization + aliases**
   - `load_processed.py`, loader tests, `test_load_processed_aliases.py`.
   - Keep API alias additions with loader behavior fixes and parallel parity tests.
3. **PR C: Plotting debt cleanup**
   - `plot_depth_histogram.py`, `plot_depth_profile.py`, `plot_enrichment.py`, `plot_enrichment_profile.py`, `plot_reads.py`.
   - Pure cleanup/debt removal; no functional change intent.

Hardening tests can remain with PR A/PR B where they logically belong (preferred), rather than a separate PR.

## Integration Decision

Proceed with the 3-PR split above.  
If reviewer throughput is constrained, fallback to a 2-PR split by combining PR C into PR B.
