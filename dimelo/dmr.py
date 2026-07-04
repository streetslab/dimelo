from __future__ import annotations

import subprocess
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from . import run_modkit, utils
from .models import ModkitDMRMultiResult, ModkitDMRPairResult


def _coerce_path(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def _coerce_sample_items(
    samples: dict[str, str | Path] | Iterable[tuple[str, str | Path]],
) -> list[tuple[str, Path]]:
    if isinstance(samples, dict):
        items = [(str(name), _coerce_path(path)) for name, path in samples.items()]
    else:
        items = [(str(name), _coerce_path(path)) for name, path in samples]
    if len(items) < 2:
        raise ValueError("modkit dmr multi requires at least two samples.")
    seen_names: set[str] = set()
    for name, path in items:
        if len(name) == 0:
            raise ValueError("Sample names for modkit dmr multi cannot be empty.")
        if name in seen_names:
            # modkit supports duplicate names (combining), but keep workflow strict by default.
            raise ValueError(f"Duplicate sample name for modkit dmr multi: {name!r}")
        seen_names.add(name)
        if not path.exists():
            raise FileNotFoundError(f"Sample bedMethyl file does not exist: {path}")
    return items


def _bed_table_from_path(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Expected modkit output file not found: {path}")
    if path.stat().st_size == 0:
        return pd.DataFrame()

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        first_line = handle.readline().strip()

    if first_line.startswith("#"):
        columns = first_line.lstrip("#").split("\t")
        return pd.read_csv(path, sep="\t", comment="#", header=None, names=columns)
    return pd.read_csv(path, sep="\t", header=None)


def _select_high_confidence_sites(
    sites_table: pd.DataFrame,
    *,
    pvalue_max: float = 0.01,
    abs_effect_size_min: float = 0.1,
    min_total_coverage: int | None = None,
) -> pd.DataFrame:
    if sites_table.empty:
        return sites_table.copy()

    out = sites_table.copy()
    mask = pd.Series(True, index=out.index)

    pvalue_col = "map_pvalue" if "map_pvalue" in out.columns else None
    if pvalue_col is not None:
        pvals = pd.to_numeric(out[pvalue_col], errors="coerce")
        mask &= pvals <= float(pvalue_max)

    effect_col = "effect_size" if "effect_size" in out.columns else None
    if effect_col is not None:
        effects = pd.to_numeric(out[effect_col], errors="coerce")
        mask &= effects.abs() >= float(abs_effect_size_min)
    elif "score" in out.columns:
        # Review fix #7: previously, when 'effect_size' was absent we silently
        # applied abs_effect_size_min to the 'score' column. That is wrong: modkit
        # effect sizes are on a small fractional scale, while the BED 'score'
        # column is a different quantity (often 0-1000), so the threshold is
        # meaningless there. Refuse to threshold the unrelated 'score' column as
        # if it were an effect size; the caller must supply a table with a real
        # 'effect_size' column to filter by effect magnitude.
        raise ValueError(
            "_select_high_confidence_sites was asked to filter by "
            "abs_effect_size_min but the table has no 'effect_size' column. The "
            "'score' column is a different, non-comparable scale and must not be "
            "thresholded as an effect size; supply a table with an 'effect_size' "
            "column (e.g. run modkit dmr with --header) to filter by effect size."
        )

    if min_total_coverage is not None:
        coverage_col = None
        for candidate in ("a_total", "b_total", "num_sites"):
            if candidate in out.columns:
                coverage_col = candidate
                break
        if coverage_col is not None:
            coverage = pd.to_numeric(out[coverage_col], errors="coerce")
            mask &= coverage >= int(min_total_coverage)

    return out.loc[mask].reset_index(drop=True)


def _append_if_value(command: list[str], flag: str, value: Any) -> None:
    if value is None:
        return
    command.extend([flag, str(value)])


def run_dmr_pair(
    *,
    control_bed_methyl: str | Path,
    experiment_bed_methyl: str | Path,
    ref_genome: str | Path,
    out_path: str | Path,
    regions_bed: str | Path | None = None,
    segment_path: str | Path | None = None,
    bases: Sequence[str] = ("A",),
    assign_codes: Sequence[str] | None = None,
    min_valid_coverage: int = 0,
    dmr_prior: float | None = None,
    diff_stay: float | None = None,
    significance_factor: float | None = None,
    decay_distance: int | None = None,
    max_gap_size: int | None = None,
    log_transition_decay: bool = False,
    fine_grained: bool = False,
    prior_alpha: float | None = None,
    prior_beta: float | None = None,
    delta: float | None = None,
    n_sample_records: int | None = None,
    max_coverages: tuple[int, int] | None = None,
    cap_coverages: bool = False,
    missing: str | None = None,
    threads: int | None = None,
    io_threads: int | None = None,
    batch_size: int | None = None,
    interval_size: int | None = None,
    header: bool = True,
    force: bool = True,
    suppress_progress: bool = True,
    log_filepath: str | Path | None = None,
    modkit_executable: str | Path | None = None,
    pvalue_max: float = 0.01,
    abs_effect_size_min: float = 0.1,
    min_total_coverage: int | None = None,
) -> ModkitDMRPairResult:
    """
    Run modkit's native pairwise DMR pipeline (including optional HMM segmentation).

    Notes:
    - `regions_bed` and `segment_path` are mutually exclusive in modkit.
    - When `segment_path` is set, the output includes single-site statistics plus
      segmented regions from the HMM model.
    """
    if regions_bed is not None and segment_path is not None:
        raise ValueError(
            "modkit dmr pair does not allow --regions-bed with --segment in one command. "
            "Run region-level and segmented analyses as separate runs."
        )
    if len(bases) == 0:
        raise ValueError("run_dmr_pair requires at least one primary base in `bases`.")

    capabilities = run_modkit._ensure_modkit_available(
        quiet=True,
        executable=modkit_executable,
    )
    executable = capabilities.executable
    control_path = _coerce_path(control_bed_methyl)
    experiment_path = _coerce_path(experiment_bed_methyl)
    reference_path = _coerce_path(ref_genome)
    output_path = _coerce_path(out_path)
    region_path = None if regions_bed is None else _coerce_path(regions_bed)
    segment_output_path = None if segment_path is None else _coerce_path(segment_path)
    if not control_path.exists():
        raise FileNotFoundError(
            f"Control bedMethyl file does not exist: {control_path}"
        )
    if not experiment_path.exists():
        raise FileNotFoundError(
            f"Experiment bedMethyl file does not exist: {experiment_path}"
        )
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference_path}")
    if region_path is not None and not region_path.exists():
        raise FileNotFoundError(f"regions_bed does not exist: {region_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if segment_output_path is not None:
        segment_output_path.parent.mkdir(parents=True, exist_ok=True)

    command: list[str] = [
        executable,
        "dmr",
        "pair",
        "-a",
        str(control_path),
        "-b",
        str(experiment_path),
        "--ref",
        str(reference_path),
        "--out-path",
        str(output_path),
    ]
    if region_path is not None:
        command.extend(["--regions-bed", str(region_path)])
    if segment_output_path is not None:
        command.extend(["--segment", str(segment_output_path)])

    for base in bases:
        command.extend(["--base", str(base)])
    for code in assign_codes or ():
        command.extend(["--assign-code", str(code)])

    _append_if_value(command, "--min-valid-coverage", min_valid_coverage)
    _append_if_value(command, "--dmr-prior", dmr_prior)
    _append_if_value(command, "--diff-stay", diff_stay)
    _append_if_value(command, "--significance-factor", significance_factor)
    _append_if_value(command, "--decay-distance", decay_distance)
    _append_if_value(command, "--max-gap-size", max_gap_size)
    _append_if_value(command, "--delta", delta)
    _append_if_value(command, "--n-sample-records", n_sample_records)
    _append_if_value(command, "--missing", missing)
    _append_if_value(command, "--batch-size", batch_size)
    _append_if_value(command, "--interval-size", interval_size)
    if prior_alpha is not None or prior_beta is not None:
        if prior_alpha is None or prior_beta is None:
            raise ValueError(
                "prior_alpha and prior_beta must both be provided together."
            )
        command.extend(["--prior", str(prior_alpha), str(prior_beta)])
    if max_coverages is not None:
        if len(max_coverages) != 2:
            raise ValueError(
                "max_coverages must be a 2-tuple: (control_max, experiment_max)."
            )
        command.extend(
            ["--max-coverages", str(max_coverages[0]), str(max_coverages[1])]
        )

    resolved_threads = utils.cores_to_run(threads)
    command.extend(["--threads", str(resolved_threads)])
    if io_threads is not None:
        command.extend(["--io-threads", str(utils.cores_to_run(io_threads))])
    if log_transition_decay:
        command.append("--log-transition-decay")
    if fine_grained:
        command.append("--fine-grained")
    if cap_coverages:
        command.append("--cap-coverages")
    if header:
        command.append("--header")
    if force:
        command.append("--force")
    if suppress_progress:
        command.append("--suppress-progress")
    if log_filepath is not None:
        command.extend(["--log-filepath", str(_coerce_path(log_filepath))])

    subprocess.run(command, check=True)

    sites_table = _bed_table_from_path(output_path)
    segments_table = (
        _bed_table_from_path(segment_output_path)
        if segment_output_path is not None
        else None
    )
    high_confidence_sites = _select_high_confidence_sites(
        sites_table,
        pvalue_max=pvalue_max,
        abs_effect_size_min=abs_effect_size_min,
        min_total_coverage=min_total_coverage,
    )

    return ModkitDMRPairResult(
        output_path=output_path,
        segment_path=segment_output_path,
        command=command,
        sites=sites_table,
        segments=segments_table,
        high_confidence_sites=high_confidence_sites,
        metadata={
            "bases": list(bases),
            "regions_bed": None if region_path is None else str(region_path),
            "pvalue_max": float(pvalue_max),
            "abs_effect_size_min": float(abs_effect_size_min),
            "min_total_coverage": min_total_coverage,
            "resolved_threads": resolved_threads,
            "modkit_version": capabilities.version,
            "modkit_executable": executable,
        },
    )


def run_dmr_multi(
    *,
    samples: dict[str, str | Path] | Iterable[tuple[str, str | Path]],
    regions_bed: str | Path,
    ref_genome: str | Path,
    out_dir: str | Path,
    bases: Sequence[str] = ("A",),
    assign_codes: Sequence[str] | None = None,
    min_valid_coverage: int = 0,
    missing: str | None = None,
    threads: int | None = None,
    io_threads: int | None = None,
    prefix: str | None = None,
    header: bool = True,
    force: bool = True,
    suppress_progress: bool = True,
    log_filepath: str | Path | None = None,
    modkit_executable: str | Path | None = None,
) -> ModkitDMRMultiResult:
    """
    Run modkit's native multi-sample DMR workflow.
    """
    if len(bases) == 0:
        raise ValueError("run_dmr_multi requires at least one primary base in `bases`.")
    sample_items = _coerce_sample_items(samples)

    capabilities = run_modkit._ensure_modkit_available(
        quiet=True,
        executable=modkit_executable,
    )
    executable = capabilities.executable
    region_path = _coerce_path(regions_bed)
    reference_path = _coerce_path(ref_genome)
    output_dir = _coerce_path(out_dir)
    if not region_path.exists():
        raise FileNotFoundError(f"regions_bed does not exist: {region_path}")
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference_path}")
    output_dir.mkdir(parents=True, exist_ok=True)

    command: list[str] = [
        executable,
        "dmr",
        "multi",
        "--regions-bed",
        str(region_path),
        "--ref",
        str(reference_path),
        "--out-dir",
        str(output_dir),
    ]
    for name, path in sample_items:
        command.extend(["--sample", str(path), str(name)])
    for base in bases:
        command.extend(["--base", str(base)])
    for code in assign_codes or ():
        command.extend(["--assign-code", str(code)])

    _append_if_value(command, "--min-valid-coverage", min_valid_coverage)
    _append_if_value(command, "--missing", missing)
    if prefix is not None:
        command.extend(["--prefix", str(prefix)])

    resolved_threads = utils.cores_to_run(threads)
    command.extend(["--threads", str(resolved_threads)])
    if io_threads is not None:
        command.extend(["--io-threads", str(utils.cores_to_run(io_threads))])
    if header:
        command.append("--header")
    if force:
        command.append("--force")
    if suppress_progress:
        command.append("--suppress-progress")
    if log_filepath is not None:
        command.extend(["--log-filepath", str(_coerce_path(log_filepath))])

    subprocess.run(command, check=True)

    # Review fix #6: scope the result manifest to this run's prefix instead of
    # globbing every *.bed in out_dir. The unscoped glob picked up stale outputs
    # from previous runs (and ignored `prefix` entirely), so the returned
    # manifest could list files that this invocation never produced. When a
    # prefix is set, modkit names its outputs "<prefix>*.bed"; restrict to that.
    glob_pattern = f"{prefix}*.bed" if prefix is not None else "*.bed"
    pair_paths = sorted(output_dir.glob(glob_pattern))
    pair_rows: list[dict[str, Any]] = []
    for bed_path in pair_paths:
        pair_rows.append(
            {
                "pair_file": bed_path,
                "pair_name": bed_path.stem,
            }
        )
    pair_files = pd.DataFrame(pair_rows)

    return ModkitDMRMultiResult(
        out_dir=output_dir,
        command=command,
        pair_files=pair_files,
        metadata={
            "bases": list(bases),
            "regions_bed": str(region_path),
            "resolved_threads": resolved_threads,
            "n_samples": len(sample_items),
            "modkit_version": capabilities.version,
            "modkit_executable": executable,
        },
    )
