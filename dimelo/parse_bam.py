import gzip
import json
import os
import subprocess
import time
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
VERIFY_CACHE_FILENAME = ".dimelo.verify_cache.json"
VERIFY_CACHE_VERSION = 1
VERIFY_CACHE_MAX_ENTRIES = 64
DEFAULT_MODKIT_MEMORY_PER_THREAD_BYTES = 2 * 1024 * 1024 * 1024
DEFAULT_MODKIT_MEMORY_FRACTION = 0.8
IUPAC_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

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
    thresh: float | None = None,
    window_size: int | None = None,
    cores: int | None = None,
    log: bool = False,
    cleanup: bool = True,
    quiet: bool = False,
    override_checks: bool = False,
    verification_cache: bool = True,
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
        regions: TODO
        motifs: a list of strings specifying which base modifications to look for.
            The basemods are each specified as {sequence_motif},{position_of_modification}.
            For example, a methylated adenine is specified as 'A,0' and CpG methylation
            is specified as 'CG,0'.
        thresh: float point number specifying the base modification probability threshold
            used to delineate modificaton calls as True or False. When set to None, modkit
            will select its own threshold automatically based on the data.
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
        override_checks: convert errors from input checking into warnings if True

    Returns:
        Path object pointing to the compressed and indexed .bed.gz bedmethyl file, ready
        for plotting functions.
        Path object pointing to 'regions.processed.bed', the `--include-bed` file used for `modkit pileup`

    """
    """
    TODO: There are a lot of issues that are all related here:
    dimelo/parse_bam.py:150: error: Incompatible types in assignment (expression has type "Path | None", variable has type "str | Path")  [assignment]
    dimelo/parse_bam.py:169: error: Argument "input_file" to "prep_outputs" has incompatible type "str | Path"; expected "Path"  [arg-type]
    dimelo/parse_bam.py:256: error: Argument "input_file" to "run_with_progress_bars" has incompatible type "str | Path"; expected "Path"  [arg-type]
    dimelo/parse_bam.py:257: error: Argument "ref_genome" to "run_with_progress_bars" has incompatible type "str | Path"; expected "Path"  [arg-type]
    
    I'm not sure of the most elegant way to fix it. Come back and address.
    """

    ## Verify and prepare inputs and outputs

    run_modkit.ensure_modkit_available(quiet=quiet)

    input_file, ref_genome, output_directory = utils.sanitize_path_args(
        input_file, ref_genome, output_directory
    )

    try:
        verify_inputs(
            input_file,
            motifs,
            ref_genome,
            quiet=quiet,
            cache=verification_cache,
        )
    except Exception as e:
        if override_checks:
            if not quiet:
                print(f"WARNING: {e}")
        else:
            raise Exception(
                f'{e}\nIf you are confident that your inputs are ok, pass "override_checks=True" to convert to warning and proceed with processing.'
            ) from e

    output_path, (output_bedmethyl, output_bedmethyl_sorted, output_pileup_path, _) = (
        prep_output_directory(
            output_directory=output_directory,
            output_name=output_name,
            input_file=input_file,
            output_file_names=[
                "pileup.bed",
                "pileup.sorted.bed",
                "pileup.sorted.bed.gz",
                "pileup.sorted.bed.gz.tbi",
            ],
        )
    )

    ## Build up the command list to be sent to modkit, then run modkit

    # TODO: This is mildly confusing. I get what it's doing, but it's hard to follow / names are bad. Also, why is it used in cleanup here, but not in extract?
    region_command_list, processed_regions_path = create_region_command_list(
        output_path,
        regions,
        window_size,
    )

    motif_command_list = []
    if len(motifs) > 0:
        seen_motifs: set[tuple[str, int]] = set()
        for motif in motifs:
            parsed_motif = utils.ParsedMotif(motif)
            motif_key = (parsed_motif.motif_seq, parsed_motif.modified_pos)
            # This motif is already going to be processed; we want to skip adding it a second
            # time because modkit does not like duplicate motifs.
            # It's actually ok if it's a different mod code in the two cases because the pileup
            # operation, under the hood, keeps all mod codes. Filtering is only done when loading.
            if motif_key in seen_motifs:
                continue
            seen_motifs.add(motif_key)
            motif_command_list.extend(
                ["--motif", parsed_motif.motif_seq, str(parsed_motif.modified_pos)]
            )
    else:
        raise ValueError("Error: no motifs specified. Nothing to process.")

    if log:
        if not quiet:
            print("Logging to ", Path(output_path) / "pileup-log")
        log_command_list = ["--log-filepath", Path(output_path) / "pileup-log"]
    else:
        log_command_list = []

    cores_command_list = _resolve_modkit_threads_command(cores=cores, quiet=quiet)

    # TODO: This is SO SO SO similar to extract; just the ValueError vs. printing. I think this can be resolved
    mod_thresh_command_list: list[str] = []
    if thresh is None:
        if not quiet:
            print(
                "No base modification threshold provided. Using adaptive threshold selection via modkit."
            )
    else:
        adjusted_threshold = utils.adjust_threshold(thresh, quiet=quiet)
        if adjusted_threshold < 0.5 and not quiet:
            print(
                f"WARNING: thresh {thresh} is very low and may lead to unexpected behavior. Typical thresholds are at least 0.5 or 128."
            )
        for motif in motifs:
            parsed_motif = utils.ParsedMotif(motif)
            for mod_code in parsed_motif.mod_codes:
                mod_thresh_command_list.extend(
                    ["--mod-thresholds", f"{mod_code}:{adjusted_threshold}"]
                )

    ref_genome_command_list = ["--ref", ref_genome]
    filter_command_list = ["--filter-threshold", "0"]

    pileup_command_list = (
        ["modkit", "pileup", input_file, output_bedmethyl]
        + region_command_list
        + motif_command_list
        + ref_genome_command_list
        + filter_command_list
        + mod_thresh_command_list
        + cores_command_list
        + log_command_list
    )

    # TODO: Do we need to store and use the output from this method? Previously was being printed immediately afterward.
    _ = run_modkit.run_with_progress_bars(
        command_list=pileup_command_list,
        input_file=input_file,
        ref_genome=ref_genome,
        motifs=motifs,
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
    # print(done_string)

    ## Sort, compress, and index the output bedmethyl file

    _sort_bedmethyl_output(
        input_bedmethyl=output_bedmethyl,
        output_bedmethyl_sorted=output_bedmethyl_sorted,
        requested_cores=cores,
        quiet=quiet,
    )
    pysam.tabix_compress(output_bedmethyl_sorted, output_pileup_path, force=True)
    pysam.tabix_index(str(output_pileup_path), preset="bed", force=True)

    # TODO: Can cleanup be consolidated?
    if cleanup:
        if output_bedmethyl.exists():
            output_bedmethyl.unlink()
        if output_bedmethyl_sorted.exists():
            output_bedmethyl_sorted.unlink()

    return output_pileup_path, processed_regions_path


def extract(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path | None = None,
    regions: str | Path | list[str | Path] | None = None,
    motifs: list = ["A,0", "CG,0", "GCH,1"],
    thresh: float | None = None,
    window_size: int | None = None,
    cores: int | None = None,
    log: bool = False,
    cleanup: bool = True,
    quiet: bool = False,
    override_checks: bool = False,
    verification_cache: bool = True,
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
        regions: TODO
        motifs: a list of strings specifying which base modifications to look for.
            The basemods are each specified as {sequence_motif},{position_of_modification}.
            For example, a methylated adenine is specified as 'A,0' and CpG methylation
            is specified as 'CG,0'.
        thresh: float point number specifying the base modification probability threshold
            used to delineate modificaton calls as True or False. When set to None, modkit
            will select its own threshold automatically based on the data.
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
        override_checks: convert errors from input checking into warnings if True

    Returns:
        Path object pointing to the compressed and indexed output .h5 file, ready for
        plotting functions.
        Path object pointing to 'regions.processed.bed', the `--include-bed` file used for `modkit extract`

    """
    """
    TODO: There are a lot of issues that are all related here:
    dimelo/parse_bam.py:374: error: Incompatible types in assignment (expression has type "Path | None", variable has type "str | Path")  [assignment]
    dimelo/parse_bam.py:393: error: Argument "input_file" to "prep_outputs" has incompatible type "str | Path"; expected "Path"  [arg-type]
    dimelo/parse_bam.py:480: error: Argument "input_file" to "run_with_progress_bars" has incompatible type "str | Path"; expected "Path"  [arg-type]
    dimelo/parse_bam.py:481: error: Argument "ref_genome" to "run_with_progress_bars" has incompatible type "str | Path"; expected "Path"  [arg-type]
    
    I'm not sure of the most elegant way to fix it. Come back and address.
    """

    ## Verify and prepare inputs and outputs

    input_file, ref_genome, output_directory = utils.sanitize_path_args(
        input_file, ref_genome, output_directory
    )

    try:
        verify_inputs(
            input_file,
            motifs,
            ref_genome,
            quiet=quiet,
            cache=verification_cache,
        )
    except Exception as e:
        if override_checks:
            if not quiet:
                print(f"WARNING: {e}")
        else:
            raise Exception(
                f'{e}\nIf you are confident that your inputs are ok, pass "override_checks=True" to convert to warning and proceed with processing.'
            ) from e

    # TODO: Add intermediate mod-specific .txt files?
    output_path, (output_reads_path,) = prep_output_directory(
        output_directory=output_directory,
        output_name=output_name,
        input_file=input_file,
        output_file_names=["reads.combined_basemods.h5"],
    )

    ## Build up the command lists shared across motifs to be sent to modkit

    region_command_list, processed_regions_path = create_region_command_list(
        output_path,
        regions,
        window_size,
    )

    cores_command_list = _resolve_modkit_threads_command(cores=cores, quiet=quiet)

    if thresh is None:
        if not quiet:
            print(
                "No valid base modification threshold provided. Raw probs will be saved."
            )
        adjusted_threshold = None
    else:
        adjusted_threshold = utils.adjust_threshold(thresh, quiet=quiet)
        if adjusted_threshold < 0.5 and not quiet:
            print(
                f"WARNING: thresh {thresh} is very low and may lead to unexpected behavior. Typical thresholds are at least 0.5 or 128."
            )

    if log:
        if not quiet:
            print("logging to ", Path(output_path) / "extract-log")
        log_command_list = ["--log-filepath", Path(output_path) / "extract-log"]
    else:
        log_command_list = []

    ref_genome_command_list = ["--ref", ref_genome]
    filter_command_list = ["--filter-threshold", "0"]

    motif_command_list = _build_extract_motif_command_list(motifs)
    output_txt = Path(output_path) / "reads.combined_motifs.txt"
    output_txt.unlink(missing_ok=True)

    extract_command_list = (
        ["modkit", "extract", input_file, output_txt]
        + region_command_list
        + motif_command_list
        + cores_command_list
        + log_command_list
        + ref_genome_command_list
        + filter_command_list
    )

    _ = run_modkit.run_with_progress_bars(
        command_list=extract_command_list,
        input_file=input_file,
        ref_genome=ref_genome,
        motifs=motifs,
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

    for motif in motifs:
        read_by_base_txt_to_hdf5(
            output_txt,
            output_reads_path,
            motif,
            adjusted_threshold,
            quiet=quiet,
            require_motif_context=True,
        )
    if cleanup:
        output_txt.unlink(missing_ok=True)

    return output_reads_path, processed_regions_path


"""
Helper functions to facilitate bam parse operations
"""


def _file_fingerprint(path: Path) -> dict[str, int | str]:
    stat = path.stat()
    return {"path": str(path.resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _verification_cache_key(
    *,
    input_file: Path,
    ref_genome: Path,
    motifs: list[str],
) -> str:
    sorted_motifs = sorted(str(motif) for motif in motifs)
    return "|".join([str(input_file.resolve()), str(ref_genome.resolve()), ",".join(sorted_motifs)])


def _verification_cache_path(input_file: Path) -> Path:
    return input_file.parent / VERIFY_CACHE_FILENAME


def _load_verification_cache(cache_path: Path) -> dict:
    if not cache_path.exists():
        return {"version": VERIFY_CACHE_VERSION, "entries": {}}
    try:
        data = json.loads(cache_path.read_text(encoding="utf-8"))
    except Exception:
        return {"version": VERIFY_CACHE_VERSION, "entries": {}}
    if not isinstance(data, dict):
        return {"version": VERIFY_CACHE_VERSION, "entries": {}}
    if data.get("version") != VERIFY_CACHE_VERSION:
        return {"version": VERIFY_CACHE_VERSION, "entries": {}}
    entries = data.get("entries")
    if not isinstance(entries, dict):
        return {"version": VERIFY_CACHE_VERSION, "entries": {}}
    return {"version": VERIFY_CACHE_VERSION, "entries": entries}


def _save_verification_cache(cache_path: Path, cache_data: dict) -> None:
    temp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
    temp_path.write_text(json.dumps(cache_data), encoding="utf-8")
    temp_path.replace(cache_path)


def _load_cached_verification_result(
    *,
    input_file: Path,
    ref_genome: Path,
    motifs: list[str],
) -> tuple[int, int] | None:
    if os.environ.get("DIMELO_DISABLE_VERIFY_CACHE", "").lower() in {"1", "true", "yes"}:
        return None

    cache_path = _verification_cache_path(input_file)
    cache_data = _load_verification_cache(cache_path)
    key = _verification_cache_key(input_file=input_file, ref_genome=ref_genome, motifs=motifs)
    entry = cache_data["entries"].get(key)
    if not isinstance(entry, dict):
        return None

    try:
        if entry["input"] != _file_fingerprint(input_file):
            return None
        if entry["ref"] != _file_fingerprint(ref_genome):
            return None
        total_bases = int(entry["total_bases"])
        correct_bases = int(entry["correct_bases"])
    except Exception:
        return None

    if total_bases <= 0:
        return None
    if correct_bases / total_bases < 0.35:
        return None
    return correct_bases, total_bases


def _store_cached_verification_result(
    *,
    input_file: Path,
    ref_genome: Path,
    motifs: list[str],
    correct_bases: int,
    total_bases: int,
) -> None:
    if os.environ.get("DIMELO_DISABLE_VERIFY_CACHE", "").lower() in {"1", "true", "yes"}:
        return

    cache_path = _verification_cache_path(input_file)
    cache_data = _load_verification_cache(cache_path)
    entries: dict[str, dict] = cache_data["entries"]
    key = _verification_cache_key(input_file=input_file, ref_genome=ref_genome, motifs=motifs)
    entries[key] = {
        "input": _file_fingerprint(input_file),
        "ref": _file_fingerprint(ref_genome),
        "motifs": sorted(str(motif) for motif in motifs),
        "correct_bases": int(correct_bases),
        "total_bases": int(total_bases),
        "timestamp_ns": time.time_ns(),
    }

    if len(entries) > VERIFY_CACHE_MAX_ENTRIES:
        sorted_keys = sorted(
            entries,
            key=lambda entry_key: int(entries[entry_key].get("timestamp_ns", 0)),
        )
        for stale_key in sorted_keys[: len(entries) - VERIFY_CACHE_MAX_ENTRIES]:
            entries.pop(stale_key, None)

    try:
        _save_verification_cache(cache_path, cache_data)
    except Exception:
        # Best effort cache only.
        return


def verify_inputs(
    input_file: Path,
    motifs: list[str],
    ref_genome: Path,
    *,
    quiet: bool = False,
    cache: bool = True,
):
    """
    Checks .bam format and alignment quality (to verify that you are using the right reference genome)

    The correct-bases-called fraction, if under 35%, means the user almost definitely passed the wrong reference genome.
    """
    if cache:
        cached = _load_cached_verification_result(
            input_file=input_file,
            ref_genome=ref_genome,
            motifs=motifs,
        )
        if cached is not None:
            if not quiet:
                print("Using cached parse_bam input verification results.")
            return

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
    if cache:
        _store_cached_verification_result(
            input_file=input_file,
            ref_genome=ref_genome,
            motifs=motifs,
            correct_bases=correct_bases,
            total_bases=total_bases,
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
        for counter, read in enumerate(input_bam.fetch()):
            if read.has_tag("Mm") or read.has_tag("Ml"):
                raise ValueError(
                    f'Base modification tags are out of spec (Mm and Ml instead of MM and ML). \n\nConsider using "modkit update-tags {str(bam_file)} new_file.bam" in the command line with your conda environment active and then trying with the new file. For megalodon basecalling/modcalling, you may also need to pass "--mode ambiguous.\nBe sure to index the resulting .bam file."'
                )

            if read.has_tag("MM"):
                mm_tag = read.get_tag("MM")
                for tag_substring in mm_tag.split(";"):
                    if not tag_substring:
                        continue
                    tag_value = tag_substring.split(",", 1)[0]
                    if (
                        len(tag_value) > 0
                        and tag_value[-1] != "?"
                        and tag_value[-1] != "."
                    ):
                        raise ValueError(
                            f'Base modification tags are out of spec. Need ? or . in TAG:TYPE:VALUE for MM tag, else modified probability is considered to be implicit. \n\nConsider using "modkit update-tags {str(bam_file)} new_file.bam --mode ambiguous" in the command line with your conda environment active and then trying with the new file.'
                        )
                    if len(tag_value) > 2 and tag_value[0] in basemods_found_dict:
                        correct_mod_codes = mod_codes_dict[tag_value[0]]
                        if tag_value[2] in correct_mod_codes:
                            basemods_found_dict[tag_value[0]] = True
                        else:
                            mod_codes_found_dict[tag_value[0]].add(tag_value[2])

            if all(basemods_found_dict.values()):
                return
            if counter >= NUM_READS_TO_CHECK:
                missing_bases = []
                for base, found in basemods_found_dict.items():
                    if not found:
                        missing_bases.append(base)
                print(
                    f"""
WARNING: no modified appropriately-coded values found for {missing_bases} in the first {counter} reads. 
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
            aligned_pairs = read.get_aligned_pairs(matches_only=True)
            if len(aligned_pairs) == 0:
                continue
            reference_positions = [pos_in_ref for _, pos_in_ref in aligned_pairs]
            ref_start = min(reference_positions)
            ref_end = max(reference_positions) + 1
            ref_window = genome_fasta.fetch(read.reference_name, ref_start, ref_end)

            total_bases += len(aligned_pairs)
            correct_bases += sum(
                1
                for pos_in_read, pos_in_ref in aligned_pairs
                if read_sequence[pos_in_read] == ref_window[pos_in_ref - ref_start]
            )

    return correct_bases, total_bases


def create_region_command_list(
    output_path,
    regions,
    window_size,
):
    """
    Creates commands to pass to modkit for specifying genomic regions.

    TODO: Split into two function? Convert to bed, then construct commands
    """

    if regions is not None:
        bed_filepath_processed = output_path / "regions.processed.bed"
        regions_dict = utils.regions_dict_from_input(
            regions,
            window_size,
        )
        utils.bed_from_regions_dict(regions_dict, bed_filepath_processed)
        region_specifier = ["--include-bed", str(bed_filepath_processed)]

    else:
        bed_filepath_processed = None
        region_specifier = []

    return region_specifier, bed_filepath_processed


def _build_extract_motif_command_list(motifs: list[str]) -> list[str]:
    seen_motif_specs: set[tuple[str, int]] = set()
    motif_command_list: list[str] = []
    for motif in motifs:
        parsed_motif = utils.ParsedMotif(motif)
        motif_spec = (parsed_motif.motif_seq, parsed_motif.modified_pos)
        if motif_spec in seen_motif_specs:
            continue
        seen_motif_specs.add(motif_spec)
        motif_command_list.extend(
            ["--motif", parsed_motif.motif_seq, str(parsed_motif.modified_pos)]
        )
    if not motif_command_list:
        raise ValueError("Error: no motifs specified. Nothing to process.")
    return motif_command_list


def _parse_positive_int(value: str | None) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    digits = []
    for char in text:
        if char.isdigit():
            digits.append(char)
        else:
            break
    if not digits:
        return None
    parsed = int("".join(digits))
    return parsed if parsed > 0 else None


def _slurm_memory_limit_bytes() -> int | None:
    mem_per_node_mb = _parse_positive_int(os.environ.get("SLURM_MEM_PER_NODE"))
    mem_per_cpu_mb = _parse_positive_int(os.environ.get("SLURM_MEM_PER_CPU"))
    cpus_per_task = _parse_positive_int(os.environ.get("SLURM_CPUS_PER_TASK"))
    cpus_on_node = _parse_positive_int(os.environ.get("SLURM_CPUS_ON_NODE"))

    candidates: list[int] = []
    if mem_per_node_mb is not None:
        candidates.append(mem_per_node_mb * 1024 * 1024)
    if mem_per_cpu_mb is not None:
        cpus = cpus_per_task or cpus_on_node or utils.cores_to_run(None)
        candidates.append(mem_per_cpu_mb * cpus * 1024 * 1024)
    if not candidates:
        return None
    return min(candidates)


def _available_memory_bytes_for_modkit() -> int | None:
    override = _parse_positive_int(os.environ.get("DIMELO_AVAILABLE_MEMORY_BYTES"))
    if override is not None:
        return override

    candidates: list[int] = []
    slurm_limit = _slurm_memory_limit_bytes()
    if slurm_limit is not None:
        candidates.append(slurm_limit)

    try:
        import psutil  # type: ignore

        candidates.append(int(psutil.virtual_memory().available))
    except Exception:
        pass

    if not candidates:
        return None
    return min(candidates)


def _memory_limited_modkit_threads(cpu_limited_threads: int) -> int:
    available = _available_memory_bytes_for_modkit()
    if available is None:
        return cpu_limited_threads

    per_thread = _parse_positive_int(os.environ.get("DIMELO_MODKIT_MEMORY_PER_THREAD_BYTES"))
    if per_thread is None:
        per_thread = DEFAULT_MODKIT_MEMORY_PER_THREAD_BYTES

    fraction = DEFAULT_MODKIT_MEMORY_FRACTION
    fraction_env = os.environ.get("DIMELO_MODKIT_MEMORY_FRACTION")
    if fraction_env:
        try:
            parsed_fraction = float(fraction_env)
            if 0 < parsed_fraction <= 1:
                fraction = parsed_fraction
        except ValueError:
            pass

    budget = int(available * fraction)
    if budget <= 0:
        return 1
    memory_cap = max(1, budget // max(1, per_thread))
    return max(1, min(cpu_limited_threads, memory_cap))


def _sort_bedmethyl_output(
    *,
    input_bedmethyl: Path,
    output_bedmethyl_sorted: Path,
    requested_cores: int | None,
    quiet: bool,
) -> None:
    sort_threads = utils.cores_to_run(requested_cores)
    sort_tmpdir = os.environ.get("SLURM_TMPDIR")
    basic_command = ["sort", "-k1,1", "-k2,2n", str(input_bedmethyl)]
    attempts: list[list[str]] = []
    if sort_threads > 1:
        parallel_command = [
            "sort",
            "--parallel",
            str(sort_threads),
            "-k1,1",
            "-k2,2n",
            str(input_bedmethyl),
        ]
        if sort_tmpdir:
            attempts.append(parallel_command[:-1] + ["-T", sort_tmpdir, parallel_command[-1]])
        attempts.append(parallel_command)
    if sort_tmpdir:
        attempts.append(["sort", "-k1,1", "-k2,2n", "-T", sort_tmpdir, str(input_bedmethyl)])
    attempts.append(basic_command)

    last_error: subprocess.CalledProcessError | None = None
    for attempt_index, command in enumerate(attempts):
        with output_bedmethyl_sorted.open("w") as sorted_file:
            try:
                subprocess.run(command, stdout=sorted_file, stderr=subprocess.PIPE, text=True, check=True)
                return
            except subprocess.CalledProcessError as err:
                last_error = err
                if not quiet and attempt_index == 0:
                    print("Falling back to compatible sort options after --parallel sort failure.")
    if last_error is not None:
        raise last_error


def _resolve_modkit_threads_command(*, cores: int | None, quiet: bool) -> list[str]:
    available_cores = utils.cores_to_run(None)
    if cores is None:
        resolved_cores = _memory_limited_modkit_threads(available_cores)
        if not quiet:
            if resolved_cores < available_cores:
                print(
                    f"No specified number of cores requested. CPU allows {available_cores}, memory-aware cap allocates {resolved_cores}."
                )
            else:
                print(
                    f"No specified number of cores requested. {available_cores} available in this runtime, allocating all."
                )
        return ["--threads", str(resolved_cores)]

    requested_cores = int(cores)
    resolved_cores = utils.cores_to_run(requested_cores)
    if requested_cores != resolved_cores and not quiet:
        print(
            f"Warning: requested {requested_cores} cores, allocating {resolved_cores} based on available CPU/cluster constraints."
        )
    elif not quiet:
        print(f"Allocating requested {resolved_cores} cores.")
    return ["--threads", str(resolved_cores)]


def _build_compressed_read_vectors(
    *,
    valid_coordinates: list[int],
    mod_values: list[float],
    read_len: int,
    thresh: float | None,
    compress_level: int,
) -> tuple[np.ndarray, np.ndarray]:
    if len(valid_coordinates) > 0:
        read_len_along_ref = max(valid_coordinates) + 1
    else:
        read_len_along_ref = read_len

    coordinate_array = np.asarray(valid_coordinates, dtype=np.intp)

    mod_vector = np.zeros(read_len_along_ref, dtype=np.uint8)
    if coordinate_array.size > 0:
        if thresh is None:
            # We subtract 0.25 because in modkit they add 0.5, but our elements are zero when the
            # base motif isn't present, so to get things to round to the right integers to match the
            # original .bam file, subtracting 0.25 is good. Anything from 0.001 to 0.4999 would work.
            mod_values_uint8 = np.rint(np.asarray(mod_values) * 256 - 0.25).astype(np.uint8)
        else:
            mod_values_uint8 = np.asarray(mod_values, dtype=np.uint8)
        mod_vector[coordinate_array] = mod_values_uint8

    valid_vector = np.zeros(read_len_along_ref, dtype=np.uint8)
    if coordinate_array.size > 0:
        valid_vector[coordinate_array] = 1

    mod_vector_compressed = np.frombuffer(
        gzip.compress(mod_vector.tobytes(), compresslevel=compress_level),
        dtype=np.uint8,
    )
    valid_vector_compressed = np.frombuffer(
        gzip.compress(valid_vector.tobytes(), compresslevel=compress_level),
        dtype=np.uint8,
    )
    return mod_vector_compressed, valid_vector_compressed


def _iupac_bases_for_char(base_char: str) -> set[str]:
    return IUPAC_BASES.get(base_char.upper(), {base_char.upper()})


def _kmer_matches_parsed_motif(kmer_context: str, parsed_motif: utils.ParsedMotif) -> bool:
    kmer = kmer_context.strip().upper()
    motif_seq = parsed_motif.motif_seq.upper()
    if len(kmer) == 0 or len(kmer) < len(motif_seq):
        return False

    kmer_center = len(kmer) // 2
    for start_index in range(0, len(kmer) - len(motif_seq) + 1):
        if start_index + parsed_motif.modified_pos != kmer_center:
            continue
        matches = True
        for kmer_char, motif_char in zip(
            kmer[start_index : start_index + len(motif_seq)],
            motif_seq,
        ):
            if _iupac_bases_for_char(kmer_char).isdisjoint(
                _iupac_bases_for_char(motif_char)
            ):
                matches = False
                break
        if matches:
            return True
    return False


def read_by_base_txt_to_hdf5(
    input_txt: str | Path,
    output_h5: str | Path,
    motif: str,
    thresh: float | None = None,
    quiet: bool = False,
    compress_level: int = 1,
    chunk_size: int = 1000,
    require_motif_context: bool = False,
) -> None:
    """
    Takes in a txt file generated by modkit extract and appends
    all the data from a specified motif into an hdf5 file. If a thresh is specified, it
    also binarizes the mod calls.

    If the h5 file does not exist it will be created and datasets will be added for read_name,
    chromosome, read_start, read_end, strand, motif, mod_vector, and val_vector.

    All the datasets (exception threshold) are parallel arrays of length num_reads

    Each read's position data is defined in genomic reference coordinates on the positive strand
    (i.e. the read_start is the leftmost aligned position, read_end is the rightmost, vectors
    are left to right along genomic coordinates)

    TODO: Make a nice key:value map of the h5 file structure, make sure start and end are documented
    as reconstructions NOT original cigarstring alignment info. mention pysam

    Args:
        input_txt: a string or Path pointing to a modkit extracted base-by-base modifications
            file. This file is assumed to have been created by modkit v0.2.4, other versions may
            have a different format and may not function normally.
        output_h5: a string or Path pointing to a valid place to save an .h5 file. If this
            file already exists, it will not be cleared and will simply be appended to.
        motif: a string specifying a single base modification. Basemods are specified as
            {sequence_motif},{position_of_modification},{optional mod_code}. For example,
            a methylated adenine is specified as 'A,0' or 'A,0,a' and CpG methylation is
            specified as 'CG,0' or 'CG,0,m'.
        thresh: a floating point threshold for base modification calling, between zero and one.
            If specified as None, raw probabilities will be saved in the .h5 output.
        quiet: if True, this suppresses outputs
        compress_level: gzip compression level for datasets, specifically for vectors for now
        chunk_size: size of write chunks in reads
        require_motif_context: if True, additionally filter rows to motif sequence using the kmer
            context column from modkit extract output (used when one extract pass includes many motifs)

    Returns:
        None

    """
    """
    TODO: There are some issues that are all related here:
    dimelo/parse_bam.py:718: error: Incompatible types in assignment (expression has type "Path | None", variable has type "str | Path")  [assignment]
    dimelo/parse_bam.py:725: error: Item "str" of "str | Path" has no attribute "open"  [union-attr]
    dimelo/parse_bam.py:890: error: Item "str" of "str | Path" has no attribute "name"  [union-attr]
    
    I'm not sure of the most elegant way to fix it. Come back and address.
    """
    input_txt, output_h5 = utils.sanitize_path_args(input_txt, output_h5)

    parsed_motif = utils.ParsedMotif(motif)

    with input_txt.open() as txt:
        with h5py.File(output_h5, "a") as h5:
            ## Define hdf5 dataset types for later
            dt_str = h5py.string_dtype(encoding="utf-8")
            # mod and val vectors -> uint8 allows us to just write whatever bytes we want
            # h5py does not appear to otherwise support vlen binary
            dt_vlen = h5py.vlen_dtype(np.dtype("uint8"))

            ## Format threshold value and create dataset to store whether this data is thresholded (binary) or raw (float16)
            # TODO: should this method thresholding without binarization
            # None becomes NaN
            threshold_to_store = np.nan if thresh is None else thresh
            if "threshold" in h5:
                threshold_from_existing = h5["threshold"][()]
                if threshold_from_existing != threshold_to_store and not (
                    np.isnan(threshold_from_existing) and np.isnan(threshold_to_store)
                ):
                    raise ValueError(
                        "existing threshold in output_h5 does not match provided threshold for read_by_base_txt_to_hdf5."
                    )
            else:
                h5.create_dataset("threshold", data=threshold_to_store)

            if "read_name" in h5:
                old_size = h5["read_name"].shape[0]
            else:
                old_size = 0
                h5.create_dataset(
                    "read_name",
                    (old_size,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression="gzip",
                    compression_opts=9,
                )

            ## Create read metadata datasets
            # TODO: loop through dict instead?
            if "chromosome" in h5:
                if old_size != h5["chromosome"].shape[0]:
                    print("size mismatch: read_name:chromosome")
            else:
                h5.create_dataset(
                    "chromosome",
                    (old_size,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression="gzip",
                    compression_opts=9,
                )
            if "read_start" in h5:
                if old_size != h5["read_start"].shape[0]:
                    print("size mismatch", "read_name", "read_start")
            else:
                h5.create_dataset(
                    "read_start",
                    (old_size,),
                    maxshape=(None,),
                    dtype="i",
                    compression="gzip",
                    compression_opts=9,
                )
            if "read_end" in h5:
                if old_size != h5["read_end"].shape[0]:
                    print("size mismatch", "read_name", "read_end")
            else:
                h5.create_dataset(
                    "read_end",
                    (old_size,),
                    maxshape=(None,),
                    dtype="i",
                    compression="gzip",
                    compression_opts=9,
                )
            if "strand" in h5:
                if old_size != h5["strand"].shape[0]:
                    print("size mismatch", "read_name", "strand")
            else:
                h5.create_dataset(
                    "strand",
                    (old_size,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression="gzip",
                    compression_opts=9,
                )
            if "motif" in h5:
                if old_size != h5["motif"].shape[0]:
                    print("size mismatch", "read_name", "motif")
            else:
                h5.create_dataset(
                    "motif",
                    (old_size,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression="gzip",
                    compression_opts=9,
                )

            ## Create the vector datasets. These will contain raw bytes formatted into a uint8 array
            # TODO: loop through dict instead
            if "mod_vector" in h5:
                if old_size != h5["mod_vector"].shape[0]:
                    print("size mismatch read_name:mod_vector")
            else:
                h5.create_dataset(
                    "mod_vector",
                    (old_size,),
                    maxshape=(None,),
                    dtype=dt_vlen,
                    # compression='gzip', # we are handling compression ourselves because hdf5 is bad at it
                    # compression_opts=9,
                )
            if "val_vector" in h5:
                if old_size != h5["val_vector"].shape[0]:
                    print("size mismatch read_name:val_vector")
            else:
                h5.create_dataset(
                    "val_vector",
                    (old_size,),
                    maxshape=(None,),
                    dtype=dt_vlen,
                    # compression='gzip', # we are handling compression ourselves because hdf5 is bad at it
                    # compression_opts=9,
                )

            read_name_dataset = h5["read_name"]
            chromosome_dataset = h5["chromosome"]
            read_start_dataset = h5["read_start"]
            read_end_dataset = h5["read_end"]
            strand_dataset = h5["strand"]
            motif_dataset = h5["motif"]
            mod_vector_dataset = h5["mod_vector"]
            val_vector_dataset = h5["val_vector"]

            ## Add data to datasets from txt file
            # Initialize loop vars - these will go into datasets
            # TODO: initialize read name to actual first read so we can get rid of the logic in the loop
            read_name = ""
            read_chrom = ""
            read_len = 0
            ref_strand = ""
            reverse_read = False
            read_start = 0
            read_end = 0
            valid_coordinates_list: list[int] = []
            mod_values_list: list[float] = []

            # Count reads for batched write
            write_index = old_size
            safe_chunk_size = max(1, int(chunk_size))
            reads_in_chunk = 0

            chunk_read_names: list[str] = []
            chunk_chromosomes: list[str] = []
            chunk_read_starts: list[int] = []
            chunk_read_ends: list[int] = []
            chunk_strands: list[str] = []
            chunk_motifs: list[str] = []
            chunk_mod_vectors: list[np.ndarray] = []
            chunk_val_vectors: list[np.ndarray] = []

            motif_base = parsed_motif.modified_base
            motif_mod_codes = parsed_motif.mod_codes
            raw_probabilities = thresh is None

            def _write_chunk() -> None:
                nonlocal reads_in_chunk, write_index
                if reads_in_chunk == 0:
                    return
                start_index = write_index
                end_index = write_index + reads_in_chunk
                read_name_dataset.resize((end_index,))
                chromosome_dataset.resize((end_index,))
                read_start_dataset.resize((end_index,))
                read_end_dataset.resize((end_index,))
                strand_dataset.resize((end_index,))
                motif_dataset.resize((end_index,))
                mod_vector_dataset.resize((end_index,))
                val_vector_dataset.resize((end_index,))
                read_name_dataset[start_index:end_index] = chunk_read_names
                chromosome_dataset[start_index:end_index] = chunk_chromosomes
                read_start_dataset[start_index:end_index] = chunk_read_starts
                read_end_dataset[start_index:end_index] = chunk_read_ends
                strand_dataset[start_index:end_index] = chunk_strands
                motif_dataset[start_index:end_index] = chunk_motifs
                mod_vector_dataset[start_index:end_index] = chunk_mod_vectors
                val_vector_dataset[start_index:end_index] = chunk_val_vectors
                write_index = end_index
                chunk_read_names.clear()
                chunk_chromosomes.clear()
                chunk_read_starts.clear()
                chunk_read_ends.clear()
                chunk_strands.clear()
                chunk_motifs.clear()
                chunk_mod_vectors.clear()
                chunk_val_vectors.clear()
                reads_in_chunk = 0

            def _append_current_read() -> None:
                nonlocal reads_in_chunk
                mod_vector_compressed, valid_vector_compressed = (
                    _build_compressed_read_vectors(
                        valid_coordinates=valid_coordinates_list,
                        mod_values=mod_values_list,
                        read_len=read_len,
                        thresh=thresh,
                        compress_level=compress_level,
                    )
                )
                chunk_read_names.append(read_name)
                chunk_chromosomes.append(read_chrom)
                chunk_read_starts.append(read_start)
                chunk_read_ends.append(read_end)
                chunk_strands.append(ref_strand)
                chunk_motifs.append(motif)
                chunk_mod_vectors.append(mod_vector_compressed)
                chunk_val_vectors.append(valid_vector_compressed)
                reads_in_chunk += 1
                if reads_in_chunk >= safe_chunk_size:
                    _write_chunk()

            # Setting up progress bars if not in quiet mode
            # Skip header
            iterator = txt
            header_line = next(iterator, None)
            if require_motif_context:
                if header_line is None:
                    raise ValueError(
                        "Cannot require motif context: extract output is empty (missing header)."
                    )
                header_fields = header_line.rstrip("\n").split("\t")
                if len(header_fields) <= 14:
                    raise ValueError(
                        "Cannot require motif context: extract output is missing the kmer context column."
                    )
            if not quiet:
                iterator = tqdm(
                    iterator,
                    desc=f"Transferring reads from {input_txt.name} into {output_h5.name}",
                )

            # Loop through txt file
            for line in iterator:
                # TODO: use csv module
                fields = line.rstrip("\n").split("\t")
                read_id = fields[0]
                pos_in_read = int(fields[1])
                pos_in_genome = int(fields[2])

                if read_name != read_id:
                    if len(read_name) > 0:
                        _append_current_read()

                    ## Set up for next read
                    # Metadata
                    read_name = read_id
                    read_chrom = fields[3]
                    read_len = int(fields[9])
                    ref_strand = fields[5]
                    reverse_read = ref_strand == "-"
                    if reverse_read:
                        pos_in_read_ref = read_len - pos_in_read - 1
                    else:
                        pos_in_read_ref = pos_in_read
                    # Calculate read start (leftmost position on ref genome)
                    # TODO: logic can be replaced when we switch to true read start/end from modkit
                    read_start = pos_in_genome - pos_in_read_ref
                    read_end = 0
                    # Instantiate lists
                    mod_values_list = []
                    valid_coordinates_list = []

                # TODO: verify that read position is in the right (ref) coordinate system
                if reverse_read:
                    pos_in_read_ref = read_len - pos_in_read - 1
                else:
                    pos_in_read_ref = pos_in_read

                # Adjust the read_end (rightmost position on ref genome) each time there's a new mod
                # This will lead to the most accurate end positions for gapped reads
                # TODO: logic can be replaced when we switch to true read start/end from modkit
                read_end = pos_in_genome + (read_len - pos_in_read_ref)
                # Regardless of whether its a new read or not,
                # add modification to vector if motif type is correct
                # for the motif in question
                canonical_base = fields[15]
                mod_code = fields[11]
                if (
                    canonical_base == motif_base
                    and mod_code in motif_mod_codes
                ):
                    if require_motif_context:
                        if len(fields) <= 14:
                            raise ValueError(
                                "Cannot require motif context: encountered extract row without kmer context column."
                            )
                        if not _kmer_matches_parsed_motif(fields[14], parsed_motif):
                            continue
                    valid_coordinates_list.append(pos_in_genome - read_start)
                    prob = float(fields[10])
                    if raw_probabilities:
                        mod_values_list.append(prob)
                    elif prob >= thresh:
                        mod_values_list.append(1)
                    else:
                        mod_values_list.append(0)

            # Save the last read
            # TODO: try to consolidate
            if len(read_name) > 0:
                _append_current_read()
            _write_chunk()
    return


def prep_output_directory(
    output_directory: Path | None,
    output_name: str,
    input_file: Path,
    output_file_names: list[str],
) -> tuple[Path, list[Path]]:
    """
    As a side effect, if files exist that match the requested outputs, they are deleted.

    TODO: Is it kind of silly that this takes in input_file? Maybe should take in some generic default parameter, or this default should be set outside this method?

    Args:
        output_directory: Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        input_file: Path to input file; used to define default output directory
        output_file_names: list of names of desired output files

    Returns:
        * Path to top-level output directory
        * List of Paths to requested output files
    """
    if output_directory is None:
        output_directory = input_file.parent
        print(f"No output directory provided, using input directory {output_directory}")

    output_path = output_directory / output_name

    output_files = [output_path / file_name for file_name in output_file_names]

    # Ensure output path exists, and that any of the specified output files do not already exist (necessary for some outputs)
    # Delete the files that do already exist
    output_path.mkdir(parents=True, exist_ok=True)
    for output_file in output_files:
        output_file.unlink(missing_ok=True)

    return output_path, output_files
