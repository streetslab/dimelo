import csv
import gzip
import itertools
import json
import multiprocessing
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pysam
from tqdm.auto import tqdm

from . import run_modkit, utils

"""
This module contains code to convert .bam files into both human-readable and 
indexed random-access pileup and read-wise processed outputs.
"""

"""
Global variables
"""

# Specifies how many reads to check for the base modifications of interest.
NUM_READS_TO_CHECK = 100
ThresholdInput = int | float | dict[str, int | float] | None


def _should_render_live_progress() -> bool:
    """
    Decide whether tqdm live bars should be rendered in this process context.

    "auto" mode avoids notebook/non-TTY contexts where carriage-return bars are
    frequently rendered as mangled output.
    """
    mode = os.environ.get("DIMELO_PROGRESS_MODE", "auto").strip().lower()
    if mode in {"off", "none", "false", "0"}:
        return False
    if mode in {"on", "force", "true", "1"}:
        return True
    in_notebook = "JPY_PARENT_PID" in os.environ or "IPYKERNEL_PARENT_PID" in os.environ
    return sys.stderr.isatty() and not in_notebook


def _unique_preserve_order(items: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        ordered.append(item)
    return ordered


def _build_pileup_targeting_command_list(
    motifs: list[str],
    capabilities: run_modkit.ModkitCapabilities,
) -> list[str]:
    # modkit 0.6.x requires --modified-bases. In some environments capability
    # detection can be stale during in-session binary upgrades, so treat version
    # 0.6+ as requiring modified-bases even if the help-derived flag is false.
    requires_modified_bases = capabilities.supports_modified_bases or (
        capabilities.version_tuple is not None
        and len(capabilities.version_tuple) >= 2
        and (capabilities.version_tuple[0], capabilities.version_tuple[1]) >= (0, 6)
    )

    motif_command_list: list[str] = []
    modified_bases: list[str] = []
    include_cpg = False
    motif_pairs_seen: set[tuple[str, int]] = set()

    for motif in motifs:
        parsed_motif = utils.ParsedMotif(motif)
        if requires_modified_bases:
            modified_bases.extend(
                f"{parsed_motif.modified_base}:{mod_code}"
                for mod_code in parsed_motif.mod_codes
            )

        if parsed_motif.motif_seq == "CG" and parsed_motif.modified_pos == 0:
            include_cpg = True
            continue

        if (
            requires_modified_bases
            and len(parsed_motif.motif_seq) == 1
            and parsed_motif.modified_pos == 0
        ):
            continue

        key = (parsed_motif.motif_seq, parsed_motif.modified_pos)
        if key in motif_pairs_seen:
            continue
        motif_pairs_seen.add(key)
        motif_command_list.extend(
            ["--motif", parsed_motif.motif_seq, str(parsed_motif.modified_pos)]
        )

    modified_bases_command_list: list[str] = []
    if requires_modified_bases:
        for modified_base in _unique_preserve_order(modified_bases):
            modified_bases_command_list.extend(["--modified-bases", modified_base])
        if len(modified_bases_command_list) == 0:
            raise ValueError(
                "No mod codes were resolved from motifs; cannot build required "
                "--modified-bases arguments for modkit pileup. Pass motifs with "
                "explicit mod codes (for example 'A,0,a' or 'CG,0,m')."
            )

    cpg_command_list = ["--cpg"] if include_cpg else []
    return modified_bases_command_list + cpg_command_list + motif_command_list


def _build_extract_targeting_command_list(
    parsed_motif: utils.ParsedMotif,
) -> list[str]:
    if parsed_motif.motif_seq == "CG" and parsed_motif.modified_pos == 0:
        return ["--cpg"]
    return ["--motif", parsed_motif.motif_seq, str(parsed_motif.modified_pos)]


def _build_mod_threshold_command_list(
    motifs: list[str],
    motif_thresholds: dict[str, float],
    capabilities: run_modkit.ModkitCapabilities,
) -> list[str]:
    if capabilities.supports_mod_threshold:
        threshold_flag = "--mod-threshold"
    elif capabilities.supports_mod_thresholds:
        threshold_flag = "--mod-thresholds"
    else:
        threshold_flag = "--mod-threshold"

    mod_code_thresholds: dict[str, float] = {}
    for motif in motifs:
        parsed_motif = utils.ParsedMotif(motif)
        threshold_for_motif = motif_thresholds[motif]
        for mod_code in parsed_motif.mod_codes:
            existing_threshold = mod_code_thresholds.get(mod_code)
            if (
                existing_threshold is not None
                and abs(existing_threshold - threshold_for_motif) > 1e-12
            ):
                raise ValueError(
                    "Cannot apply different thresholds to motifs that share mod code "
                    f"{mod_code!r}. Received both {existing_threshold} and {threshold_for_motif}."
                )
            mod_code_thresholds[mod_code] = threshold_for_motif

    command: list[str] = []
    for mod_code in _unique_preserve_order(list(mod_code_thresholds.keys())):
        command.extend([threshold_flag, f"{mod_code}:{mod_code_thresholds[mod_code]}"])
    return command


def _build_extract_command_prefix(
    input_file: Path,
    output_txt: Path,
    capabilities: run_modkit.ModkitCapabilities,
) -> list[str | Path]:
    if _extract_requires_subcommands(capabilities):
        return [capabilities.executable, "extract", "full", input_file, output_txt]
    return [capabilities.executable, "extract", input_file, output_txt]


def _build_extract_reference_command_list(
    ref_genome: Path,
    capabilities: run_modkit.ModkitCapabilities,
) -> list[str | Path]:
    if _extract_requires_subcommands(capabilities):
        if capabilities.extract_supports_reference_long:
            return ["--reference", ref_genome]
        # modkit 0.6+ extract subcommands require --reference for motif filtering.
        return ["--reference", ref_genome]
    if capabilities.extract_supports_reference_short:
        return ["--ref", ref_genome]
    if capabilities.extract_supports_reference_long:
        return ["--reference", ref_genome]
    return ["--ref", ref_genome]


def _extract_requires_subcommands(
    capabilities: run_modkit.ModkitCapabilities,
) -> bool:
    if capabilities.supports_extract_subcommands:
        return True
    return bool(
        capabilities.version_tuple is not None
        and len(capabilities.version_tuple) >= 2
        and (capabilities.version_tuple[0], capabilities.version_tuple[1]) >= (0, 6)
    )


def _build_implicit_tag_command_list(
    capabilities: run_modkit.ModkitCapabilities,
) -> list[str]:
    if capabilities.supports_force_allow_implicit:
        return ["--force-allow-implicit"]
    return []


def _modkit_requires_multi_motif_pileup_split(
    capabilities: run_modkit.ModkitCapabilities,
) -> bool:
    if capabilities.version_tuple is None or len(capabilities.version_tuple) < 2:
        return False
    return (capabilities.version_tuple[0], capabilities.version_tuple[1]) >= (0, 6)


def _canonical_motif_key(motif: str) -> str:
    parsed = utils.ParsedMotif(motif)
    return f"{parsed.motif_seq},{parsed.modified_pos}"


def _group_motifs_for_pileup(
    motifs: list[str],
    capabilities: run_modkit.ModkitCapabilities,
) -> list[tuple[str, list[str]]]:
    if (
        not _modkit_requires_multi_motif_pileup_split(capabilities)
        or len({_canonical_motif_key(motif) for motif in motifs}) <= 1
    ):
        return [("combined", motifs)]

    grouped: dict[str, list[str]] = {}
    ordered_keys: list[str] = []
    for motif in motifs:
        canonical_key = _canonical_motif_key(motif)
        if canonical_key not in grouped:
            grouped[canonical_key] = []
            ordered_keys.append(canonical_key)
        grouped[canonical_key].append(motif)
    return [(key, grouped[key]) for key in ordered_keys]


def _resolve_motif_thresholds(
    motifs: list[str],
    thresh: ThresholdInput,
    quiet: bool,
) -> dict[str, float] | None:
    if thresh is None:
        return None

    if isinstance(thresh, (int, float)):
        adjusted = utils.adjust_threshold(float(thresh), quiet=quiet)
        return {motif: adjusted for motif in motifs}

    if not isinstance(thresh, dict):
        raise TypeError(
            "thresh must be None, a number, or a dict mapping motif keys to thresholds."
        )

    resolved: dict[str, float] = {}
    for motif in motifs:
        canonical_key = _canonical_motif_key(motif)
        threshold_value = None
        for key in (motif, canonical_key, "default", "*"):
            if key in thresh:
                threshold_value = thresh[key]
                break
        if threshold_value is None:
            raise ValueError(
                f"Missing threshold for motif {motif!r}. Provide an exact key, canonical key "
                f"{canonical_key!r}, or a default key ('default' or '*')."
            )
        resolved[motif] = utils.adjust_threshold(float(threshold_value), quiet=quiet)
    return resolved


def _threads_command_list(cores: int | None, quiet: bool) -> list[str]:
    cores_avail = multiprocessing.cpu_count()
    if cores is None:
        if not quiet:
            print(
                f"No specified number of cores requested. {cores_avail} available on machine, allocating all."
            )
        return ["--threads", str(cores_avail)]
    if cores > cores_avail:
        if not quiet:
            print(
                f"Warning: {cores} cores request, {cores_avail} available. Allocating {cores_avail}"
            )
        return ["--threads", str(cores_avail)]
    if not quiet:
        print(f"Allocating requested {cores} cores.")
    return ["--threads", str(cores)]


def _compress_uint8_vector(
    vector: np.ndarray,
    compress_level: int,
) -> np.ndarray:
    return np.frombuffer(
        gzip.compress(vector.tobytes(), compresslevel=compress_level),
        dtype=np.uint8,
    )


def _unlink_existing(*paths: Path) -> None:
    for path in paths:
        path.unlink(missing_ok=True)


def _text_file_has_any_rows(path: Path) -> bool:
    try:
        with open(path) as handle:
            for line in handle:
                if line.strip():
                    return True
    except FileNotFoundError:
        return False
    return False


def _sanitize_motif_group_label(label: str) -> str:
    sanitized = "".join(char if char.isalnum() else "_" for char in label).strip("_")
    return sanitized or "motif"


def _reference_oriented_read_offset(
    *,
    pos_in_read: int,
    read_length: int,
    ref_strand: str,
) -> int:
    """
    Convert modkit's read-oriented position into a left-to-right reference offset.
    """
    if read_length <= 0:
        raise ValueError(f"read length must be positive; got {read_length}.")
    if pos_in_read < 0 or pos_in_read >= read_length:
        raise ValueError(
            f"pos_in_read {pos_in_read} is out of bounds for read length {read_length}."
        )
    if ref_strand == "+":
        return pos_in_read
    if ref_strand == "-":
        return read_length - pos_in_read - 1
    raise ValueError(f"Unexpected strand '{ref_strand}' in modkit extract row.")


"""
User-facing parse operations: pileup and extract
"""


def pileup(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path | None = None,
    regions: str | Path | list[str | Path] | None = None,
    motifs: list = ["A,0", "CG,0"],
    thresh: ThresholdInput = None,
    window_size: int | None = None,
    cores: int | None = None,
    log: bool = False,
    cleanup: bool = True,
    overwrite: bool = True,
    quiet: bool = False,
    override_checks: bool = False,
    modkit_executable: str | Path | None = None,
) -> tuple[Path, Path]:
    """
    Takes a bam file containing long read sequencing data aligned
    to a reference genome with modification calls for one or more base/context
    and creates a pileup. A pileup contains genome-position-wise sums of both reads with
    bases that could have the modification in question and of reads that are in
    fact modified.

    The current implementation of this method uses modkit, a tool built by
    Nanopore Technologies, along with htslib tools compress and index the output
    bedmethyl file. The modkit command for this function is `modkit pileup`.

    https://github.com/nanoporetech/modkit/

    The intermediate output file is a standard bedmethyl file containing all
    specified motifs and mod codes. The compressed and indexed output file is
    a .bed.gz file with an accompanying .bed.gz.tbi index.

    Args:
        output_file: a string or Path object pointing to the location of a .bam file.
            The file should follow at least v1.6 of the .bam file specifications,
            found here: https://samtools.github.io/hts-specs/
            https://samtools.github.io/hts-specs/SAMv1.pdf

            The file needs to have modifications stored in the standard format,
            with MM and ML tags (NOT mm and ml) and mod names m for 5mC and a
            for 6mA.

            Furthermore, the file must have a .bam.bai index file with the same name.
            You can create an index if needed using samtools index.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        ref_genome: a string of Path objecting pointing to the .fasta file
            for the reference genome to which the .bam file is aligned.
        output_directory: optional str or Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        regions: optional region selector passed through to modkit via an include-bed
            file. This may be a BED file path, a single region file, or a list of
            region files/paths that define the loci to process. When paired with
            ``window_size``, the provided regions are expanded around their centers
            before being passed to modkit.
        motifs: a list of strings specifying which base modifications to look for.
            The basemods are each specified as {sequence_motif},{position_of_modification}.
            For example, a methylated adenine is specified as 'A,0' and CpG methylation
            is specified as 'CG,0'.
        thresh: base modification threshold specification. Accepted forms are:
            * None: modkit infers thresholds from data
            * float/int: one threshold for all motifs
            * dict: motif-specific thresholds keyed by exact motif string
              (e.g. ``{"A,0": 0.7, "CG,0": 0.8}``), canonical motif key, or
              ``default``/``*`` fallback.
        window_size: an integer specifying a window around the center of each bed_file
            region. If set to None, the bed_file is used unmodified. If set to a non-zero
            positive integer, the bed_file regions are replaced by new regions with that
            window size in either direction of the center of the original bed_file regions.
            This is used for e.g. extracting information from around known motifs or peaks.
        cores: an integer specifying how many parallel cores modkit gets to use.
            By default modkit will use all of the available cores on the machine.
        log: a boolean specifying whether to output logs into the output folder.
        cleanup: a boolean specifying whether to clean up to keep intermediate
            outputs. The final processed files are not human-readable, whereas the intermediate
            outputs are. However, intermediate outputs can also be quite large.
        overwrite: when True (default), existing outputs for `output_name` are
            replaced. When False, existing completed outputs are reused and parsing
            is skipped. If partial/conflicting outputs exist, raises FileExistsError.
        override_checks: convert errors from input checking into warnings if True
        modkit_executable: optional executable name or path to a specific modkit
            binary. If None, dimelo resolves modkit from PATH (or
            DIMELO_MODKIT_EXECUTABLE when set).

    Returns:
        Path object pointing to the compressed and indexed .bed.gz bedmethyl file, ready
        for plotting functions.
        Path object pointing to 'regions.processed.bed', the `--include-bed` file used for `modkit pileup`

    """
    ## Verify and prepare inputs and outputs

    input_file, ref_genome, output_directory = utils.sanitize_path_args(
        input_file, ref_genome, output_directory
    )
    input_file = Path(input_file)
    ref_genome = Path(ref_genome)
    output_directory = None if output_directory is None else Path(output_directory)

    output_path, (output_bedmethyl, output_bedmethyl_sorted, output_pileup_path, _) = (
        prep_output_directory(
            output_directory=output_directory,
            output_name=output_name,
            input_file=input_file.parent,
            output_file_names=[
                "pileup.bed",
                "pileup.sorted.bed",
                "pileup.sorted.bed.gz",
                "pileup.sorted.bed.gz.tbi",
            ],
            overwrite=False,
        )
    )
    output_pileup_tbi = Path(f"{output_pileup_path}.tbi")
    processed_regions_candidate = output_path / "regions.processed.bed"
    processed_regions_existing = (
        processed_regions_candidate if processed_regions_candidate.exists() else None
    )

    if not overwrite:
        final_outputs_exist = output_pileup_path.exists() and output_pileup_tbi.exists()
        regions_ready = regions is None or processed_regions_existing is not None
        if final_outputs_exist and regions_ready:
            if not quiet:
                print(
                    f"overwrite=False and existing outputs found. Reusing {output_pileup_path}."
                )
            return output_pileup_path, processed_regions_existing

        existing_paths = [
            output_bedmethyl,
            output_bedmethyl_sorted,
            output_pileup_path,
            output_pileup_tbi,
        ]
        if processed_regions_existing is not None:
            existing_paths.append(processed_regions_existing)
        conflicting_paths = [path for path in existing_paths if path.exists()]
        if len(conflicting_paths) > 0:
            conflict_names = ", ".join(path.name for path in conflicting_paths)
            raise FileExistsError(
                "overwrite=False but output directory contains existing parse artifacts "
                f"for '{output_name}': {conflict_names}. "
                "Set overwrite=True to regenerate, or remove stale files."
            )
    else:
        _unlink_existing(
            output_bedmethyl,
            output_bedmethyl_sorted,
            output_pileup_path,
            output_pileup_tbi,
            processed_regions_candidate,
        )

    capabilities = run_modkit._ensure_modkit_available(
        quiet=quiet,
        executable=modkit_executable,
    )

    try:
        verify_inputs(input_file, motifs, ref_genome)
    except Exception as e:
        if override_checks:
            if not quiet:
                print(f"WARNING: {e}")
        else:
            raise Exception(
                f'{e}\nIf you are confident that your inputs are ok, pass "override_checks=True" to convert to warning and proceed with processing.'
            ) from e

    ## Build up the command list to be sent to modkit, then run modkit

    region_command_list, processed_regions_path = create_region_command_list(
        output_path,
        regions,
        window_size,
    )

    if len(motifs) == 0:
        raise ValueError("Error: no motifs specified. Nothing to process.")
    cores_command_list = _threads_command_list(cores=cores, quiet=quiet)

    motif_groups = _group_motifs_for_pileup(motifs=motifs, capabilities=capabilities)
    split_pileup_runs = len(motif_groups) > 1
    if split_pileup_runs and not quiet:
        print(
            "Detected multiple motif contexts with modkit 0.6.x; running per-motif pileups "
            "and merging outputs to avoid mixed-motif empty output behavior."
        )

    motif_thresholds = _resolve_motif_thresholds(
        motifs=motifs,
        thresh=thresh,
        quiet=quiet,
    )
    if motif_thresholds is None:
        if not quiet:
            print(
                "No base modification threshold provided. Using adaptive threshold selection via modkit."
            )
    else:
        if any(value < 0.5 for value in motif_thresholds.values()) and not quiet:
            print(
                f"WARNING: thresh {thresh} is very low and may lead to unexpected behavior. Typical thresholds are at least 0.5 or 128."
            )
    ref_genome_command_list = ["--ref", ref_genome]
    filter_command_list = ["--filter-threshold", "0"]
    implicit_tag_command_list = _build_implicit_tag_command_list(capabilities)

    intermediate_pileup_files: list[Path] = []
    for index, (group_key, group_motifs) in enumerate(motif_groups, start=1):
        group_suffix = _sanitize_motif_group_label(group_key)
        group_output_bedmethyl = (
            output_bedmethyl
            if not split_pileup_runs
            else output_path / f"pileup.part{index}.{group_suffix}.bed"
        )
        _unlink_existing(group_output_bedmethyl)
        intermediate_pileup_files.append(group_output_bedmethyl)

        group_motif_command_list = _build_pileup_targeting_command_list(
            motifs=group_motifs,
            capabilities=capabilities,
        )
        group_mod_thresh_command_list: list[str] = []
        if motif_thresholds is not None:
            group_motif_thresholds = {
                motif: motif_thresholds[motif] for motif in group_motifs
            }
            group_mod_thresh_command_list = _build_mod_threshold_command_list(
                motifs=group_motifs,
                motif_thresholds=group_motif_thresholds,
                capabilities=capabilities,
            )

        if log:
            log_path = (
                Path(output_path) / "pileup-log"
                if not split_pileup_runs
                else Path(output_path) / f"pileup-log.{index}.{group_suffix}"
            )
            log_command_list = ["--log-filepath", log_path]
            if not quiet and not split_pileup_runs:
                print("Logging to ", log_path)
        else:
            log_command_list = []

        base_pileup_command = (
            [capabilities.executable, "pileup", input_file, group_output_bedmethyl]
            + region_command_list
            + group_motif_command_list
            + ref_genome_command_list
            + filter_command_list
            + implicit_tag_command_list
            + cores_command_list
            + log_command_list
        )

        def _run_pileup_command(
            extra_args: list[str],
            *,
            base_command: list[str] = base_pileup_command,
            motifs_for_progress: list[tuple[str, int, str]] = group_motifs,
        ) -> None:
            run_modkit.run_with_progress_bars(
                command_list=base_command + extra_args,
                input_file=input_file,
                ref_genome=ref_genome,
                motifs=motifs_for_progress,
                load_fasta_regex=r"\s+\[.*?\]\s+(\d+)\s+Reading",
                find_motifs_regex=r"\s+(\d+)/(\d+)\s+finding\s+([A-Za-z0-9,]+)\s+motifs",
                contigs_progress_regex=r"\s+(\d+)/(\d+)\s+contigs",
                single_contig_regex=r"\s+(\d+)/(\d+)\s+processing\s+([\w]+)[^\w]",
                buffer_size=50,
                progress_granularity=25,
                done_str="Done",
                err_str="Error",
                expect_done=True,
                quiet=quiet,
            )

        _run_pileup_command(group_mod_thresh_command_list)

    if split_pileup_runs:
        with open(output_bedmethyl, "w") as merged_file:
            for intermediate_file in intermediate_pileup_files:
                if not intermediate_file.exists():
                    continue
                with open(intermediate_file) as part_file:
                    for line in part_file:
                        merged_file.write(line)
        for intermediate_file in intermediate_pileup_files:
            _unlink_existing(intermediate_file)

    ## Sort, compress, and index the output bedmethyl file

    with open(output_bedmethyl_sorted, "w") as sorted_file:
        subprocess.run(
            ["sort", "-k1,1", "-k2,2n", output_bedmethyl], stdout=sorted_file
        )
    pysam.tabix_compress(output_bedmethyl_sorted, output_pileup_path, force=True)
    pysam.tabix_index(str(output_pileup_path), preset="bed", force=True)

    if cleanup:
        _unlink_existing(output_bedmethyl, output_bedmethyl_sorted)

    return output_pileup_path, processed_regions_path


def extract(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path | None = None,
    regions: str | Path | list[str | Path] | None = None,
    motifs: list = ["A,0", "CG,0", "GCH,1"],
    thresh: ThresholdInput = None,
    window_size: int | None = None,
    cores: int | None = None,
    log: bool = False,
    cleanup: bool = True,
    overwrite: bool = True,
    quiet: bool = False,
    override_checks: bool = False,
    modkit_executable: str | Path | None = None,
) -> tuple[Path, Path]:
    """
    Takes a bam file containing long read sequencing data aligned
    to a reference genome with modification calls for one or more bases/contexts
    and pulls out data from each individual read.

    The current implementation of this method uses modkit, a tool built by
    Nanopore Technologies, along with h5py to build the final output file. The
    modkit command in this function is `modkit extract`.

    https://github.com/nanoporetech/modkit/

    The intermediate outputs are plain text files containing a list of all base modifications,
    with a file for each motif. The compressed and indexed output contains vectors of valid
    and modified positions within each read.

    Args:
        output_file: a string or Path object pointing to the location of a .bam file.
            The file should follow at least v1.6 of the .bam file specifications,
            found here: https://samtools.github.io/hts-specs/
            https://samtools.github.io/hts-specs/SAMv1.pdf

            The file needs to have modifications stored in the standard format,
            with MM and ML tags (NOT mm and ml) and mod names m for 5mC and a
            for 6mA.

            Furthermore, the file must have a .bam.bai index file with the same name.
            You can create an index if needed using samtools index.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        ref_genome: a string of Path objecting pointing to the .fasta file
            for the reference genome to which the .bam file is aligned.
        output_directory: optional str or Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        regions: optional region selector passed through to modkit via an include-bed
            file. This may be a BED file path, a single region file, or a list of
            region files/paths that define the loci to process. When paired with
            ``window_size``, the provided regions are expanded around their centers
            before being passed to modkit.
        motifs: a list of strings specifying which base modifications to look for.
            The basemods are each specified as {sequence_motif},{position_of_modification}.
            For example, a methylated adenine is specified as 'A,0' and CpG methylation
            is specified as 'CG,0'.
        thresh: base modification threshold specification. Accepted forms are:
            * None: keep raw probabilities in output vectors
            * float/int: one threshold for all motifs
            * dict: motif-specific thresholds keyed by exact motif string
              (e.g. ``{"A,0": 0.7, "CG,0": 0.8}``), canonical motif key, or
              ``default``/``*`` fallback.
        window_size: an integer specifying a window around the center of each bed_file
            region. If set to None, the bed_file is used unmodified. If set to a non-zero
            positive integer, the bed_file regions are replaced by new regions with that
            window size in either direction of the center of the original bed_file regions.
            This is used for e.g. extracting information from around known motifs or peaks.
        cores: an integer specifying how many parallel cores modkit gets to use.
            By default modkit will use all of the available cores on the machine.
        log: a boolean specifying whether to output logs into the output folder.
        cleanup: a boolean specifying whether to clean up to keep intermediate
            outputs. The final processed files are not human-readable, whereas the intermediate
            outputs are. However, intermediate outputs can also be quite large.
        overwrite: when True (default), existing outputs for `output_name` are
            replaced. When False, existing completed outputs are reused and parsing
            is skipped. If partial/conflicting outputs exist, raises FileExistsError.
        override_checks: convert errors from input checking into warnings if True
        modkit_executable: optional executable name or path to a specific modkit
            binary. If None, dimelo resolves modkit from PATH (or
            DIMELO_MODKIT_EXECUTABLE when set).

    Returns:
        Path object pointing to the compressed and indexed output .h5 file, ready for
        plotting functions.
        Path object pointing to 'regions.processed.bed', the `--include-bed` file used for `modkit extract`

    """
    ## Verify and prepare inputs and outputs
    input_file, ref_genome, output_directory = utils.sanitize_path_args(
        input_file, ref_genome, output_directory
    )
    input_file = Path(input_file)
    ref_genome = Path(ref_genome)
    output_directory = None if output_directory is None else Path(output_directory)

    output_path, (output_reads_path,) = prep_output_directory(
        output_directory=output_directory,
        output_name=output_name,
        input_file=input_file.parent,
        output_file_names=["reads.combined_basemods.h5"],
        overwrite=False,
    )
    processed_regions_candidate = output_path / "regions.processed.bed"
    processed_regions_existing = (
        processed_regions_candidate if processed_regions_candidate.exists() else None
    )

    if not overwrite:
        final_output_exists = output_reads_path.exists()
        regions_ready = regions is None or processed_regions_existing is not None
        if final_output_exists and regions_ready:
            if not quiet:
                print(
                    f"overwrite=False and existing outputs found. Reusing {output_reads_path}."
                )
            return output_reads_path, processed_regions_existing

        existing_paths = [output_reads_path]
        if processed_regions_existing is not None:
            existing_paths.append(processed_regions_existing)
        conflicting_paths = [path for path in existing_paths if path.exists()]
        if len(conflicting_paths) > 0:
            conflict_names = ", ".join(path.name for path in conflicting_paths)
            raise FileExistsError(
                "overwrite=False but output directory contains existing parse artifacts "
                f"for '{output_name}': {conflict_names}. "
                "Set overwrite=True to regenerate, or remove stale files."
            )
    else:
        _unlink_existing(output_reads_path, processed_regions_candidate)

    capabilities = run_modkit._ensure_modkit_available(
        quiet=quiet,
        executable=modkit_executable,
    )

    try:
        verify_inputs(input_file, motifs, ref_genome)
    except Exception as e:
        if override_checks:
            if not quiet:
                print(f"WARNING: {e}")
        else:
            raise Exception(
                f'{e}\nIf you are confident that your inputs are ok, pass "override_checks=True" to convert to warning and proceed with processing.'
            ) from e

    ## Build up the command lists shared across motifs to be sent to modkit

    region_command_list, processed_regions_path = create_region_command_list(
        output_path,
        regions,
        window_size,
    )

    cores_command_list = _threads_command_list(cores=cores, quiet=quiet)

    motif_thresholds = _resolve_motif_thresholds(
        motifs=motifs,
        thresh=thresh,
        quiet=quiet,
    )
    if motif_thresholds is None:
        if not quiet:
            print(
                "No valid base modification threshold provided. Raw probabilities will be saved."
            )
    else:
        if any(value < 0.5 for value in motif_thresholds.values()) and not quiet:
            print(
                f"WARNING: thresh {thresh} is very low and may lead to unexpected behavior. Typical thresholds are at least 0.5 or 128."
            )
        if not quiet:
            print(
                "Threshold provided. The extract text will stay probability-valued and read_by_base_txt_to_hdf5 will binarize at write time."
            )

    if log:
        if not quiet:
            print("logging to ", Path(output_path) / "extract-log")
        log_command_list = ["--log-filepath", Path(output_path) / "extract-log"]
    else:
        log_command_list = []

    ref_genome_command_list = _build_extract_reference_command_list(
        ref_genome=ref_genome,
        capabilities=capabilities,
    )
    filter_command_list = (
        []
        if _extract_requires_subcommands(capabilities)
        else ["--filter-threshold", "0"]
    )
    implicit_tag_command_list = _build_implicit_tag_command_list(capabilities)

    # Run modkit once for each motif, because the output .txt can be ambiguous otherwise
    # There is no column currently to specify the motif (e.g. CG,0 vs GCH,1), only canonical
    # base (e.g. C) and mod code (e.g. m)
    # There is a 5mer context so we could technically manually motif check if we want to.
    # Our current design paradigm is to leave all such operations to modkit, hence the loop below.
    for motif in motifs:
        # Here we prepare the motif-specific commands and delete any old .txt file because
        # modkit will crash otherwise
        parsed_motif = utils.ParsedMotif(motif)
        motif_command_list = _build_extract_targeting_command_list(parsed_motif)

        output_txt = Path(output_path) / (f"reads.{motif}.txt")

        output_txt.unlink(missing_ok=True)

        extract_command_list = (
            _build_extract_command_prefix(
                input_file=input_file,
                output_txt=output_txt,
                capabilities=capabilities,
            )
            + region_command_list
            + motif_command_list
            + cores_command_list
            + log_command_list
            + ref_genome_command_list
            + filter_command_list
            + implicit_tag_command_list
        )

        run_modkit.run_with_progress_bars(
            command_list=extract_command_list,
            input_file=input_file,
            ref_genome=ref_genome,
            motifs=[motif],
            load_fasta_regex=r"\s+\[.*?\]\s+(\d+)\s+parsing FASTA",
            find_motifs_regex=r"\s+(\d+)/(\d+)\s+([\w]+)\s+searched",
            contigs_progress_regex=r"\s+(\d+)/(\d+)\s+contigs\s+[^s]",
            single_contig_regex=r"\s+(\d+)/(\d+)\s+processing\s+([\w]+)[^\w]",
            buffer_size=100,
            progress_granularity=50,
            done_str="Done",
            err_str="Error",
            expect_done=False,
            quiet=quiet,
        )

        # Create the compressed and indexed output
        read_by_base_txt_to_hdf5(
            output_txt,
            output_reads_path,
            motif,
            None if motif_thresholds is None else motif_thresholds[motif],
            quiet=quiet,
        )
        # Delete intermediate file
        if cleanup:
            output_txt.unlink()

    return output_reads_path, processed_regions_path


"""
Helper functions to facilitate bam parse operations
"""


def verify_inputs(
    input_file,
    motifs,
    ref_genome,
):
    """
    Checks .bam format and alignment quality (to verify that you are using the right reference genome)

    The correct-bases-called fraction, if under 35%, means the user almost definitely passed the wrong reference genome.
    """
    check_bam_format(input_file, motifs)
    correct_bases, total_bases = get_alignment_quality(input_file, ref_genome)
    if total_bases == 0:
        raise ValueError(
            f"First {NUM_READS_TO_CHECK} reads are empty. Please verify your {input_file.name} contents."
        )
    elif correct_bases / total_bases < 0.35:
        raise ValueError(
            f"First {NUM_READS_TO_CHECK} reads have anomalously low alignment quality: only {100 * correct_bases / total_bases}% of bases align.\nPlease verify that {input_file.name} is actually aligned to {ref_genome.name}."
        )
    return


def check_bam_format(
    bam_file: str | Path,
    motifs: list = ["A,0", "CG,0"],
):
    """
    Check whether a .bam file is formatted appropriately for modkit.
    * bam file has a .bai index
    * modification tags named MM/ML (NOT Mm/Ml)
    * tags contain ambiguity specification (? vs.)
    * bam file contains the expected modifications (motif, mod code)

    Args:
        bam_file: a formatted .bam file with a .bai index
        motifs: a list of base modification motifs

    Returns:
        None. If the function returns, you are ok.

    """
    basemods_found_dict = {}
    mod_codes_dict = {}
    mod_codes_found_dict = defaultdict(set)
    for motif in motifs:
        parsed_motif = utils.ParsedMotif(motif)
        mod_codes_dict[parsed_motif.modified_base] = parsed_motif.mod_codes
        basemods_found_dict[parsed_motif.modified_base] = False

    input_bam = pysam.AlignmentFile(bam_file)

    try:
        reads_checked = 0
        for counter, read in enumerate(
            itertools.islice(input_bam.fetch(), NUM_READS_TO_CHECK)
        ):
            reads_checked = counter + 1
            read_dict = read.to_dict()
            for tag_string in read_dict["tags"]:
                tag = tag_string.split(",")[0].split(":")[0]
                if tag == "Mm" or tag == "Ml":
                    raise ValueError(
                        f'Base modification tags are out of spec (Mm and Ml instead of MM and ML). \n\nConsider using "modkit update-tags {str(bam_file)} new_file.bam" in the command line with your conda environment active and then trying with the new file. For megalodon basecalling/modcalling, you may also need to pass "--mode ambiguous.\nBe sure to index the resulting .bam file."'
                    )
                elif tag == "MM":
                    for tag_substring in tag_string.split(";"):
                        tag_fields = tag_substring.split(",")[0].split(":")
                        if len(tag_fields) >= 3:
                            tag_value = tag_fields[2]
                        else:
                            tag_value = tag_fields[0]
                        if (
                            len(tag_value) > 0
                            and tag_value[-1] != "?"
                            and tag_value[-1] != "."
                        ):
                            raise ValueError(
                                f'Base modification tags are out of spec. Need ? or . in TAG:TYPE:VALUE for MM tag, else modified probability is considered to be implicit. \n\nConsider using "modkit update-tags {str(bam_file)} new_file.bam --mode ambiguous" in the command line with your conda environment active and then trying with the new file.'
                            )
                        else:
                            if (
                                len(tag_value) > 0
                                and tag_value[0] in basemods_found_dict
                            ):
                                correct_mod_codes = mod_codes_dict[tag_value[0]]
                                # valid_mod_codes = mod_codes_dict[tag_value[0]].union(
                                #     utils.BASEMOD_NAMES_DICT[tag_value[0]]
                                # )
                                if tag_value[2] in correct_mod_codes:
                                    basemods_found_dict[tag_value[0]] = True
                                else:
                                    mod_codes_found_dict[tag_value[0]].add(tag_value[2])
                                # With the mode-code-aware motifs, it no longer makes sense to throw this error
                                # This is because the warning the user gets if their mod code isn't found, or (if none is specified)
                                # the default mod codes aren't found, can tell them what mod codes *were* found and they can add them
                                # to their motif or use adjust_mods according to what makes sense. Thus, unexpected codes are not a
                                # problem (in part this is because parse_bam will now set thresholds for the motif-specified OR default mod codes)
                                # elif tag_value[2] not in valid_mod_codes:
                                #     raise ValueError(
                                #         f'Base modification name unexpected: {tag_value[2]} to modify {tag_value[0]}, should be in set {valid_mod_codes}. \n\nIf you know what your mod names correspond to in terms of the latest .bam standard, consider using "modkit adjust-mods {str(bam_file)} new_file.bam --convert 5mC_name m --convert N6mA_name a --convert other_basemod_name correct_label" and then trying with the new file. Note: currently supported mod names are {utils.BASEMOD_NAMES_DICT}'
                                #     )
        if reads_checked == NUM_READS_TO_CHECK:
            missing_bases = []
            for base, found in basemods_found_dict.items():
                if not found:
                    missing_bases.append(base)
            print(
                f"""
WARNING: no modified appropriately-coded values found for {missing_bases} in the first {reads_checked} reads. 
Do you expect this file to contain these modifications? parse_bam is looking for {motifs} but for {missing_bases} found only found {[f"{base}+{mod_codes}" for base, mod_codes in mod_codes_found_dict.items()]}.

Consider passing only the motifs and mod codes (e.g. m,h,a) that you expect to be present in your file. 
You can use modkit adjust-mods --convert <CONVERT> <CONVERT> [OPTIONS] <IN_BAM> <OUT_BAM> to update or consolidate mod codes.
See https://github.com/nanoporetech/modkit/blob/master/book/src/advanced_usage.md
                    """
            )
            return
    except ValueError as e:
        if "fetch called on bamfile without index" in str(e):
            raise ValueError(
                f'{e}. Consider using "samtools index {str(bam_file)}" to create an index if your .bam is already sorted.'
            ) from e
        else:
            raise
    except:
        raise


def get_alignment_quality(
    bam_file,
    ref_genome,
) -> tuple[int, int]:
    """
    Determine fraction of read bases that line up with reference genome in first NUM_READS_TO_CHECK reads in bam file
    """
    ref_genome_index = ref_genome.parent / (ref_genome.name + ".fai")
    if not ref_genome_index.exists():
        print(f"Indexing {ref_genome.name}. This only needs to be done once.")
        pysam.faidx(str(ref_genome))
    input_bam = pysam.AlignmentFile(bam_file, "rb")
    genome_fasta = pysam.FastaFile(str(ref_genome))
    total_bases = 0
    correct_bases = 0
    # For NUM_READS_TO_CHECK=100 this is <1s on most machines
    for index, read in enumerate(input_bam.fetch()):
        if index >= NUM_READS_TO_CHECK:
            return correct_bases, total_bases

        # The query sequence is the entire sequence as stored in the .bam file
        # So it is reverse complemented if it was a reverse read
        # Meaning we can compare it directly against the reference genome
        read_sequence = read.query_sequence

        # print(read.mapping_quality)

        # get_aligned_pairs returns a list of (read_coord,ref_coord) pairs with None values when not aligned
        # So if we just skip Nones and compare the remainder it'll tell us the accuracy
        # (in Dorado-basecalled r10 files, as of July 2024, we observe some fraction of reads
        # with empty read_sequence, despite having intact tags and alignment info. the reason
        # for this isn't currentyl known, but with this None check we avoid errors in this alignment
        # checking stage.)

        if read_sequence is not None:
            for pos_in_read, pos_in_ref in read.get_aligned_pairs():
                if pos_in_read is not None and pos_in_ref is not None:
                    total_bases += 1
                    if read_sequence[pos_in_read] == str(
                        genome_fasta.fetch(
                            read.reference_name, pos_in_ref, pos_in_ref + 1
                        )
                    ):
                        correct_bases += 1

    return correct_bases, total_bases


def create_region_command_list(
    output_path: Path,
    regions: str | Path | list[str | Path] | None,
    window_size: int | None,
) -> tuple[list[str], Path | None]:
    """
    Create modkit `--include-bed` arguments and the processed BED path, if regions are provided.
    """

    bed_filepath_processed = _regions_to_processed_bed(
        output_path=output_path,
        regions=regions,
        window_size=window_size,
    )
    if bed_filepath_processed is None:
        return [], None
    return ["--include-bed", str(bed_filepath_processed)], bed_filepath_processed


def _regions_to_processed_bed(
    output_path: Path,
    regions: str | Path | list[str | Path] | None,
    window_size: int | None,
) -> Path | None:
    """
    Build a normalized BED file from region inputs for downstream modkit calls.
    """
    if regions is None:
        return None
    bed_filepath_processed = output_path / "regions.processed.bed"
    regions_dict = utils.regions_dict_from_input(
        regions,
        window_size,
    )
    utils.bed_from_regions_dict(regions_dict, bed_filepath_processed)
    return bed_filepath_processed


def read_by_base_txt_to_hdf5(
    input_txt: str | Path,
    output_h5: str | Path,
    motif: str,
    thresh: float | None = None,
    quiet: bool = False,
    compress_level: int = 1,
    chunk_size: int = 1000,
) -> None:
    """
    Takes in a txt file generated by modkit extract and appends all the data from a
    specified motif into an hdf5 file. If a thresh is specified, the mod calls are
    binarized at write time and the threshold is stored alongside the vectors.

    If the h5 file does not exist it will be created with top-level datasets named
    ``read_name``, ``chromosome``, ``read_start``, ``read_end``, ``strand``, ``motif``,
    ``mod_vector``, and ``val_vector``. Each dataset stores one value per read, so the arrays
    are parallel and have length ``num_reads``. The optional ``threshold`` dataset stores the
    scalar threshold used for binarization and is not part of the per-read arrays.

    Each read's position data is stored in genomic reference coordinates on the positive strand
    convention. ``read_start`` and ``read_end`` are reconstructed from the aligned modkit extract
    coordinates, not copied from the original BAM CIGAR string or any raw alignment tag. In this
    representation, ``read_start`` is the leftmost aligned reference position for the read and
    ``read_end`` is the rightmost aligned reference position observed while iterating through the
    read. The modification vectors are ordered left to right along genomic coordinates.

    Args:
        input_txt: a string or Path pointing to a modkit extracted base-by-base modifications
            file. This parser supports both legacy and current extract table schemas used by
            modkit 0.2.4 and 0.6.x.
        output_h5: a string or Path pointing to a valid place to save an .h5 file. If this
            file already exists, it will not be cleared and will simply be appended to.
        motif: a string specifying a single base modification. Basemods are specified as
            {sequence_motif},{position_of_modification},{optional mod_code}. For example,
            a methylated adenine is specified as 'A,0' or 'A,0,a' and CpG methylation is
            specified as 'CG,0' or 'CG,0,m'.
        thresh: a floating point threshold for base modification calling, between zero and one.
            If specified as None, raw probabilities will be saved in the .h5 output. If set,
            the stored mod_vector values are binary 0/1 values at motif-valid coordinates.
        quiet: if True, this suppresses outputs
        compress_level: gzip compression level for datasets, specifically for vectors for now
        chunk_size: size of write chunks in reads

    Returns:
        None

    """
    input_txt, output_h5 = utils.sanitize_path_args(input_txt, output_h5)
    input_txt = Path(input_txt)
    output_h5 = Path(output_h5)

    parsed_motif = utils.ParsedMotif(motif)

    read_name = ""
    num_reads = 0
    with input_txt.open(newline="") as txt:
        reader = csv.reader(txt, delimiter="\t")
        header_fields = next(reader, None)
        if header_fields is None:
            raise ValueError(f"modkit extract output is empty: {input_txt}")

        if "read_id" in header_fields:
            first_pass_read_name_idx = header_fields.index("read_id")
        elif "read_name" in header_fields:
            first_pass_read_name_idx = header_fields.index("read_name")
        else:
            first_pass_read_name_idx = 0

        # Check file length
        line_index = 0
        for fields in reader:
            line_index += 1
            read_id_value = fields[first_pass_read_name_idx]
            if read_name != read_id_value:
                read_name = read_id_value
                num_reads += 1
        if line_index == 0:
            raise ValueError(f"modkit extract output has no data rows: {input_txt}")
        num_lines = line_index
        txt.seek(0)

        with h5py.File(output_h5, "a") as h5:
            ## Define hdf5 dataset types for later
            dt_str = h5py.string_dtype(encoding="utf-8")
            # mod and val vectors -> uint8 allows us to just write whatever bytes we want
            # h5py does not appear to otherwise support vlen binary
            dt_vlen = h5py.vlen_dtype(np.dtype("uint8"))

            # Threshold metadata:
            # - `threshold` remains for backwards compatibility and for fast binarized/raw checks.
            # - `threshold_by_motif_json` stores per-motif thresholds so different motifs can use
            #   different binarization cutoffs in one output file.
            threshold_for_motif = None if thresh is None else float(thresh)
            threshold_by_motif: dict[str, float | None] = {}
            if "threshold_by_motif_json" in h5:
                raw_threshold_map = h5["threshold_by_motif_json"][()]
                if isinstance(raw_threshold_map, bytes):
                    raw_threshold_map = raw_threshold_map.decode("utf-8")
                threshold_by_motif = json.loads(str(raw_threshold_map))
            elif "threshold" in h5 and "motif" in h5 and h5["motif"].shape[0] > 0:
                legacy_threshold = h5["threshold"][()]
                legacy_threshold_value = (
                    None if np.isnan(legacy_threshold) else float(legacy_threshold)
                )
                existing_motifs = np.unique(np.array(h5["motif"], dtype=str))
                threshold_by_motif = {
                    existing_motif: legacy_threshold_value
                    for existing_motif in existing_motifs
                }

            threshold_by_motif[motif] = threshold_for_motif

            has_raw = any(value is None for value in threshold_by_motif.values())
            has_binary = any(value is not None for value in threshold_by_motif.values())
            if has_raw and has_binary:
                raise ValueError(
                    "Cannot mix raw-probability and thresholded motifs in one output_h5 file."
                )

            if has_binary:
                unique_thresholds = sorted(
                    {
                        float(value)
                        for value in threshold_by_motif.values()
                        if value is not None
                    }
                )
                if len(unique_thresholds) == 1:
                    threshold_to_store = unique_thresholds[0]
                else:
                    # Keep legacy `threshold` dataset non-NaN to indicate binarized vectors.
                    threshold_to_store = 0.0
            else:
                threshold_to_store = np.nan

            if "threshold" in h5:
                del h5["threshold"]
            h5.create_dataset("threshold", data=threshold_to_store)

            if "threshold_by_motif_json" in h5:
                del h5["threshold_by_motif_json"]
            h5.create_dataset(
                "threshold_by_motif_json",
                data=json.dumps(threshold_by_motif),
            )

            ## Create read metadata datasets
            if "read_name" in h5:
                old_size = h5["read_name"].shape[0]
                h5["read_name"].resize((old_size + num_reads,))
            else:
                old_size = 0
                h5.create_dataset(
                    "read_name",
                    (num_reads,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression="gzip",
                    compression_opts=9,
                )

            def ensure_dataset_capacity(
                *,
                name: str,
                dtype,
                mismatch_message: str,
                compression: str | None = None,
                compression_opts: int | None = None,
            ) -> None:
                if name in h5:
                    if old_size != h5[name].shape[0]:
                        print(mismatch_message)
                    else:
                        h5[name].resize((old_size + num_reads,))
                    return

                create_kwargs: dict[str, object] = {
                    "maxshape": (None,),
                    "dtype": dtype,
                }
                if compression is not None:
                    create_kwargs["compression"] = compression
                if compression_opts is not None:
                    create_kwargs["compression_opts"] = compression_opts

                h5.create_dataset(
                    name,
                    (num_reads,),
                    **create_kwargs,
                )

            for dataset_name, dataset_dtype, mismatch_message in [
                ("chromosome", dt_str, "size mismatch: read_name:chromosome"),
                ("read_start", "i", "size mismatch read_name:read_start"),
                ("read_end", "i", "size mismatch read_name:read_end"),
                ("strand", dt_str, "size mismatch read_name:strand"),
                ("motif", dt_str, "size mismatch read_name:motif"),
            ]:
                ensure_dataset_capacity(
                    name=dataset_name,
                    dtype=dataset_dtype,
                    mismatch_message=mismatch_message,
                    compression="gzip",
                    compression_opts=9,
                )

            ## Create the vector datasets. These will contain raw bytes formatted into a uint8 array
            for dataset_name in ["mod_vector", "val_vector"]:
                ensure_dataset_capacity(
                    name=dataset_name,
                    dtype=dt_vlen,
                    mismatch_message=f"size mismatch read_name:{dataset_name}",
                )

            ## Add data to datasets from txt file
            # Initialize loop vars - these will go into datasets
            read_name: str | None = None
            read_chrom = ""
            ref_strand = ""
            read_start = 0
            read_end = 0
            valid_genomic_positions_list: list[int] = []
            mod_values_list: list[float] = []

            # Count reads for batched write
            read_counter = 0
            # Keys (strings): dataset names, values: lists of dataset values by read; string or ints or arrays
            # Contents reset at the end of each chunk, after writing to h5
            chunk_rows: defaultdict[str, list[str | int | np.ndarray]] = defaultdict(
                list
            )
            chunk_row_count = 0

            def flush_pending_chunk_to_h5() -> None:
                nonlocal chunk_row_count, chunk_rows
                if chunk_row_count <= 0:
                    return
                start_index = old_size + read_counter - chunk_row_count
                end_index = old_size + read_counter
                for dataset, entry in chunk_rows.items():
                    h5[dataset][start_index:end_index] = entry
                chunk_rows = defaultdict(list)
                chunk_row_count = 0

            def flush_current_read() -> None:
                nonlocal read_counter, chunk_row_count, chunk_rows

                if read_name is None:
                    return

                read_len_along_ref = max(read_end - read_start, 1)

                # Populate mod vector array appropriately based on thresh settings.
                # Option 1: when thresh is provided, write binary 0/1 values; otherwise preserve raw probs.
                mod_vector = np.zeros(read_len_along_ref, dtype=np.uint8)
                valid_vector = np.zeros(read_len_along_ref, dtype=np.uint8)

                if len(valid_genomic_positions_list) > 0:
                    valid_coordinates = (
                        np.asarray(valid_genomic_positions_list, dtype=int) - read_start
                    )
                    mod_values = np.asarray(mod_values_list, dtype=float)
                    in_bounds = (valid_coordinates >= 0) & (
                        valid_coordinates < read_len_along_ref
                    )
                    valid_coordinates = valid_coordinates[in_bounds]
                    mod_values = mod_values[in_bounds]

                    if thresh is None:
                        # We subtract 0.25 because in modkit they add 0.5, but our elements are zero when the
                        # base motif isn't present, so to get things to round to the right integers to match the
                        # original .bam file, subtracting 0.25 is good. Anything from 0.001 to 0.4999 would work I think
                        mod_vector[valid_coordinates] = np.rint(
                            mod_values * 256 - 0.25
                        ).astype(np.uint8)
                    else:
                        mod_vector[valid_coordinates] = mod_values.astype(np.uint8)
                    valid_vector[valid_coordinates] = 1

                mod_vector_compressed = _compress_uint8_vector(
                    mod_vector, compress_level=compress_level
                )
                valid_vector_compressed = _compress_uint8_vector(
                    valid_vector, compress_level=compress_level
                )

                chunk_rows["read_name"].append(read_name)
                chunk_rows["chromosome"].append(read_chrom)
                chunk_rows["read_start"].append(read_start)
                chunk_rows["read_end"].append(read_end)
                chunk_rows["strand"].append(ref_strand)
                chunk_rows["motif"].append(motif)
                chunk_rows["mod_vector"].append(mod_vector_compressed)
                chunk_rows["val_vector"].append(valid_vector_compressed)

                read_counter += 1
                chunk_row_count += 1
                if chunk_row_count >= chunk_size:
                    flush_pending_chunk_to_h5()

            # Setting up progress bars if not in quiet mode
            # Skip header
            reader = csv.reader(txt, delimiter="\t")
            header = next(reader)

            def _column_index(*candidates: str, fallback: int) -> int:
                for candidate in candidates:
                    if candidate in header:
                        return header.index(candidate)
                return fallback

            read_name_idx = _column_index("read_id", "read_name", fallback=0)
            pos_in_read_idx = _column_index(
                "forward_read_position", "read_position", fallback=1
            )
            pos_in_genome_idx = _column_index(
                "ref_position", "reference_position", fallback=2
            )
            read_chrom_idx = _column_index("chrom", "chromosome", fallback=3)
            ref_strand_idx = _column_index("ref_strand", fallback=5)
            line_read_len_idx = _column_index("read_length", fallback=9)
            prob_idx = _column_index("mod_qual", fallback=10)
            mod_code_idx = _column_index("mod_code", fallback=11)
            canonical_base_idx = _column_index("canonical_base", fallback=15)

            iterator = reader
            if not quiet and _should_render_live_progress():
                iterator = tqdm(
                    iterator,
                    total=num_lines,
                    desc=f"Transferring {num_reads} from {input_txt.name} into {output_h5.name}, new size {old_size + num_reads}",
                    bar_format="{bar}| {desc} {percentage:3.0f}% | {elapsed}<{remaining}",
                )

            # Loop through txt file
            for fields in iterator:
                pos_in_genome = int(fields[pos_in_genome_idx])
                canonical_base = fields[canonical_base_idx]
                prob = float(fields[prob_idx])
                mod_code = fields[mod_code_idx]
                pos_in_read = int(fields[pos_in_read_idx])
                line_read_len = int(fields[line_read_len_idx])
                line_ref_strand = fields[ref_strand_idx]
                pos_in_read_ref = _reference_oriented_read_offset(
                    pos_in_read=pos_in_read,
                    read_length=line_read_len,
                    ref_strand=line_ref_strand,
                )
                start_candidate = pos_in_genome - pos_in_read_ref
                end_candidate = start_candidate + line_read_len

                if read_name != fields[read_name_idx]:
                    flush_current_read()

                    ## Set up for next read
                    # Metadata
                    read_name = fields[read_name_idx]
                    read_chrom = fields[read_chrom_idx]
                    ref_strand = line_ref_strand
                    read_start = start_candidate
                    read_end = end_candidate
                    # Instantiate lists
                    mod_values_list = []
                    valid_genomic_positions_list = []

                # keep read extents in reference coordinates by using inferred starts/ends from
                # each line, which is robust even when rows are not ordered by genomic position
                read_start = min(read_start, start_candidate)
                read_end = max(read_end, end_candidate)

                # Regardless of whether its a new read or not,
                # add modification to vector if motif type is correct
                # for the motif in question
                if (
                    canonical_base == parsed_motif.modified_base
                    and mod_code in parsed_motif.mod_codes
                ):
                    valid_genomic_positions_list.append(pos_in_genome)
                    if thresh is None:
                        mod_values_list.append(prob)
                    elif prob >= thresh:
                        mod_values_list.append(1)
                    else:
                        mod_values_list.append(0)

            flush_current_read()
            flush_pending_chunk_to_h5()
    return


def prep_output_directory(
    output_directory: Path | None,
    output_name: str,
    input_file: Path,
    output_file_names: list[str],
    overwrite: bool = True,
) -> tuple[Path, list[Path]]:
    """
    As a side effect, if files exist that match the requested outputs, they are deleted.

    Args:
        output_directory: Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        input_file: default output directory when output_directory is None
        output_file_names: list of names of desired output files
        overwrite: when True, existing requested output files are deleted. When
            False, existing files are preserved.

    Returns:
        * Path to top-level output directory
        * List of Paths to requested output files
    """
    if output_directory is None:
        output_directory = input_file
        print(f"No output directory provided, using input directory {output_directory}")

    output_path = output_directory / output_name

    output_files = [output_path / file_name for file_name in output_file_names]

    # Ensure output path exists.
    output_path.mkdir(parents=True, exist_ok=True)
    if overwrite:
        for output_file in output_files:
            output_file.unlink(missing_ok=True)

    return output_path, output_files
