import csv
import gzip
import itertools
import multiprocessing
import subprocess
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
    ## Verify and prepare inputs and outputs

    run_modkit.ensure_modkit_available(quiet=quiet)

    input_file, ref_genome, output_directory = utils.sanitize_path_args(
        input_file, ref_genome, output_directory
    )
    input_file = Path(input_file)
    ref_genome = Path(ref_genome)
    output_directory = None if output_directory is None else Path(output_directory)

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

    output_path, (output_bedmethyl, output_bedmethyl_sorted, output_pileup_path, _) = (
        prep_output_directory(
            output_directory=output_directory,
            output_name=output_name,
            default_directory=input_file.parent,
            output_file_names=[
                "pileup.bed",
                "pileup.sorted.bed",
                "pileup.sorted.bed.gz",
                "pileup.sorted.bed.gz.tbi",
            ],
        )
    )

    ## Build up the command list to be sent to modkit, then run modkit

    region_command_list, processed_regions_path = create_region_command_list(
        output_path,
        regions,
        window_size,
    )

    motif_command_list = []
    if len(motifs) > 0:
        for motif in motifs:
            parsed_motif = utils.ParsedMotif(motif)
            motif_command_present = False
            for a, b in zip(motif_command_list, motif_command_list[1:]):
                if a == parsed_motif.motif_seq and b == str(parsed_motif.modified_pos):
                    # This motif is already going to be processed; we want to skip adding it a second
                    # time because modkit does not like duplicate motifs.
                    # It's actually ok if it's a different mod code in the two cases because the pileup
                    # operation, under the hood, keeps all mod codes. Filtering is only done when loading.
                    motif_command_present = True
                    break
            if not motif_command_present:
                motif_command_list.append("--motif")
                motif_command_list.append(parsed_motif.motif_seq)
                motif_command_list.append(str(parsed_motif.modified_pos))
    else:
        raise ValueError("Error: no motifs specified. Nothing to process.")

    if log:
        if not quiet:
            print("Logging to ", Path(output_path) / "pileup-log")
        log_command_list = ["--log-filepath", Path(output_path) / "pileup-log"]
    else:
        log_command_list = []

    cores_command_list = _threads_command_list(cores=cores, quiet=quiet)

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
                mod_thresh_command_list = mod_thresh_command_list + [
                    "--mod-thresholds",
                    f"{mod_code}:{adjusted_threshold}",
                ]

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

    run_modkit.run_with_progress_bars(
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
    thresh: float | None = None,
    window_size: int | None = None,
    cores: int | None = None,
    log: bool = False,
    cleanup: bool = True,
    quiet: bool = False,
    override_checks: bool = False,
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
    ## Verify and prepare inputs and outputs

    input_file, ref_genome, output_directory = utils.sanitize_path_args(
        input_file, ref_genome, output_directory
    )
    input_file = Path(input_file)
    ref_genome = Path(ref_genome)
    output_directory = None if output_directory is None else Path(output_directory)

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

    output_path, (output_reads_path,) = prep_output_directory(
        output_directory=output_directory,
        output_name=output_name,
        default_directory=input_file.parent,
        output_file_names=["reads.combined_basemods.h5"],
    )

    ## Build up the command lists shared across motifs to be sent to modkit

    region_command_list, processed_regions_path = create_region_command_list(
        output_path,
        regions,
        window_size,
    )

    cores_command_list = _threads_command_list(cores=cores, quiet=quiet)

    if thresh is None:
        if not quiet:
            print(
                "No valid base modification threshold provided. Raw probabilities will be saved."
            )
        adjusted_threshold = None
    else:
        adjusted_threshold = utils.adjust_threshold(thresh, quiet=quiet)
        if adjusted_threshold < 0.5 and not quiet:
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

    ref_genome_command_list = ["--ref", ref_genome]
    filter_command_list = ["--filter-threshold", "0"]

    # Run modkit once for each motif, because the output .txt can be ambiguous otherwise
    # There is no column currently to specify the motif (e.g. CG,0 vs GCH,1), only canonical
    # base (e.g. C) and mod code (e.g. m)
    # There is a 5mer context so we could technically manually motif check if we want to.
    # Our current design paradigm is to leave all such operations to modkit, hence the loop below.
    for motif in motifs:
        # Here we prepare the motif-specific commands and delete any old .txt file because
        # modkit will crash otherwise
        motif_command_list = []
        parsed_motif = utils.ParsedMotif(motif)
        motif_command_list.append("--motif")
        motif_command_list.append(parsed_motif.motif_seq)
        motif_command_list.append(str(parsed_motif.modified_pos))

        output_txt = Path(output_path) / (f"reads.{motif}.txt")

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
            adjusted_threshold,
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
            file. This file is assumed to have been created by modkit v0.2.4, other versions may
            have a different format and may not function normally.
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
        # Check file length
        line_index = -1
        for line_index, fields in enumerate(reader):
            if line_index > 0 and read_name != fields[0]:
                read_name = fields[0]
                num_reads += 1
        if line_index < 0:
            raise ValueError(f"modkit extract output is empty: {input_txt}")
        num_lines = line_index
        txt.seek(0)

        with h5py.File(output_h5, "a") as h5:
            ## Define hdf5 dataset types for later
            dt_str = h5py.string_dtype(encoding="utf-8")
            # mod and val vectors -> uint8 allows us to just write whatever bytes we want
            # h5py does not appear to otherwise support vlen binary
            dt_vlen = h5py.vlen_dtype(np.dtype("uint8"))

            ## Format threshold value and create dataset to store whether this data is thresholded (binary) or raw.
            # Option 1: thresholded outputs are written as binary 0/1 vectors and the threshold dataset
            # records the cutoff used. Raw-probability outputs are preserved only when thresh is None.
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
            chunk_rows: defaultdict[str, list[str | int | np.ndarray]] = (
                defaultdict(list)
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
            next(reader)
            iterator = reader
            if not quiet:
                iterator = tqdm(
                    iterator,
                    total=num_lines,
                    desc=f"Transferring {num_reads} from {input_txt.name} into {output_h5.name}, new size {old_size + num_reads}",
                    bar_format="{bar}| {desc} {percentage:3.0f}% | {elapsed}<{remaining}",
                )

            # Loop through txt file
            for fields in iterator:
                pos_in_genome = int(fields[2])
                canonical_base = fields[15]
                prob = float(fields[10])
                mod_code = fields[11]
                pos_in_read = int(fields[1])
                line_read_len = int(fields[9])
                line_ref_strand = fields[5]
                # TODO: verify that read position is in the right (ref) coordinate system
                if line_ref_strand == "+":
                    pos_in_read_ref = pos_in_read
                elif line_ref_strand == "-":
                    pos_in_read_ref = line_read_len - pos_in_read - 1
                else:
                    raise ValueError(
                        f"Unexpected strand '{line_ref_strand}' in modkit extract row."
                    )
                start_candidate = pos_in_genome - pos_in_read_ref
                end_candidate = start_candidate + line_read_len

                if read_name != fields[0]:
                    flush_current_read()

                    ## Set up for next read
                    # Metadata
                    read_name = fields[0]
                    read_chrom = fields[3]
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
    default_directory: Path,
    output_file_names: list[str],
) -> tuple[Path, list[Path]]:
    """
    As a side effect, if files exist that match the requested outputs, they are deleted.

    Args:
        output_directory: Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        default_directory: default output directory when output_directory is None
        output_file_names: list of names of desired output files

    Returns:
        * Path to top-level output directory
        * List of Paths to requested output files
    """
    if output_directory is None:
        output_directory = default_directory
        print(f"No output directory provided, using input directory {output_directory}")

    output_path = output_directory / output_name

    output_files = [output_path / file_name for file_name in output_file_names]

    # Ensure output path exists, and that any of the specified output files do not already exist (necessary for some outputs)
    # Delete the files that do already exist
    output_path.mkdir(parents=True, exist_ok=True)
    for output_file in output_files:
        output_file.unlink(missing_ok=True)

    return output_path, output_files
