# Sprint 1 Parser Kickoff

**Goal:** Start the legacy cleanup sprint by inventorying parser debt in `dimelo/parse_bam.py` and `dimelo/run_modkit.py`, then move the highest-risk parser work into small, testable fixes.

**Scope:** This kickoff only covers the legacy parser lane in the `codex/legacy-cleanup` worktree. It does not touch the main workspace checkout.

---

## TODO Inventory

I found `31` `TODO` markers in the two files in scope.

### Correctness-Risk

- `dimelo/parse_bam.py:802` - thresholding behavior is still ambiguous; the comment asks whether the method should threshold without binarization.
- `dimelo/parse_bam.py:979` - vector length is derived from `max(valid_coordinates_list) + 1`; the comment notes this should switch to `read_end - read_start`, which would change output padding semantics and requires reference regeneration.
- `dimelo/parse_bam.py:1049` - the comment explicitly asks to verify that read position is in the correct reference coordinate system.
- `dimelo/parse_bam.py:1055` - read start/end logic is still a placeholder for future true modkit-provided coordinates.
- `dimelo/parse_bam.py:1063` - same placeholder logic for read end calculation.

### API/Design Cleanup

- `dimelo/parse_bam.py:83` - `regions` is still a TODO in the public docstring.
- `dimelo/parse_bam.py:111` - typed-path handling has a cluster of mypy issues that need a coherent fix.
- `dimelo/parse_bam.py:157` - cleanup naming and flow are confusing; the comment calls out different usage between cleanup and extract.
- `dimelo/parse_bam.py:191` - region-specifier construction should be a method or merged into prep logic.
- `dimelo/parse_bam.py:210` - extract and pileup diverge mainly on printing vs. raising, and the comment suggests unifying that behavior.
- `dimelo/parse_bam.py:245` - question about whether this method should store and use its output.
- `dimelo/parse_bam.py:332` - `regions` is still a TODO in the extract docstring.
- `dimelo/parse_bam.py:360` - second cluster of path/type issues in the extract path.
- `dimelo/parse_bam.py:386` - intermediate mod-specific `.txt` files are still an open API/design question.
- `dimelo/parse_bam.py:480` - same question about whether the progress-run output needs to be retained.
- `dimelo/parse_bam.py:703` - the region-command builder may need to be split into conversion plus command construction.
- `dimelo/parse_bam.py:745` - the HDF5 file layout still needs a clearer documented key/value map.
- `dimelo/parse_bam.py:769` - second cluster of type/path issues in the HDF5 conversion path.
- `dimelo/parse_bam.py:1146` - `input_file` may be the wrong abstraction for this helper’s defaulting behavior.
- `dimelo/run_modkit.py:128` - the type annotation on `finding_progress_dict` is still called out as uncertain.

### Code-Structure Cleanup

- `dimelo/parse_bam.py:273` - cleanup could be centralized instead of duplicated.
- `dimelo/parse_bam.py:782` - the function calls around the HDF5 write path can likely be consolidated.
- `dimelo/parse_bam.py:783` - both files may be opened together instead of separately.
- `dimelo/parse_bam.py:817` - repeated dataset handling could loop over a dict.
- `dimelo/parse_bam.py:903` - repeated dataset handling could loop over a dict here too.
- `dimelo/parse_bam.py:935` - read-name initialization could move out of the hot loop.
- `dimelo/parse_bam.py:952` - chunk accounting should probably use `read_counter % chunk_size`.
- `dimelo/parse_bam.py:969` - CSV parsing should use the `csv` module.
- `dimelo/parse_bam.py:999` - gzip compression for vectors should be factored into a shared helper.
- `dimelo/parse_bam.py:1085` - there is another consolidation opportunity in the final write path.
- `dimelo/parse_bam.py:1106` - gzip compression should be shared across the write code path.

## Sprint 1 Execution Order

1. Lock down the correctness-risk items first with targeted tests around coordinate handling, vector length, and threshold semantics in `dimelo/parse_bam.py`.
2. Triage the public API/design debt next by clarifying the `regions` docstrings, path typing, and helper responsibilities in the parser and `run_modkit.py`.
3. Consolidate the repeated HDF5 write-path structure only after the behavior is pinned down by tests.
4. Re-run the parser-focused regression subset after each slice so the cleanup stays incremental instead of turning into a rewrite.

## Baseline Test

- Command:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 /Users/ngamarra/.pyenv/versions/3.11.3/bin/python3.11 -m pytest tests/test_distribution.py tests/test_cluster.py -q
```

- Result:

```text
20 passed, 2 warnings in 10.79s
```

The warnings were non-fatal:

- `tests/test_cluster.py::test_plot_region_cluster_profiles_motif_index` emitted a `tight_layout` warning from `dimelo/cluster.py:1865`.
- `tests/test_cluster.py::test_cluster_read_windows_kmeans` emitted the usual `joblib` physical-core detection warning.
