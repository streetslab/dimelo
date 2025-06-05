import os
from collections import deque
from pathlib import Path

import pyBigWig
import pysam
from tqdm.auto import tqdm

from . import load_processed, utils

"""
This module contains code to export indexed and compressed parse output files to other formats that may be helpful for downstream analysis.
"""


def tail(n, iterable):
    """
    Return an iterator over the last n items.
    Copied from https://docs.python.org/3/library/itertools.html#itertools-recipes
    """

    # tail(3, 'ABCDEFG') → E F G
    return iter(deque(iterable, maxlen=n))


def pileup_to_bigwig(
    bedmethyl_file: str | Path,
    motif: str,
    bigwig_file: str | Path | None = None,
    ref_genome: str | Path | None = None,
    strand: str = ".",
    chunk_size: int = 1000,
):
    """
    Extract a single motif from a pileup and write its mod fractions by position to a bigwig file.

    This function will take the entire contents of the pileup bedmethyl file and create a bigwig header with all of the same contigs, with
    contig lengths in the bigwig header set to the highest motif coordinate for each contig. If strand is specified as + or -, only that
    strand will be written to the output bigwig - this can allow for strand bias analysis in a genome browser. If strand is specified as .,
    as is the default, both strands will be included.

    The operation can be quite slow for large pileups. The current design is that if you want to create a bigwig for a subset of the genome,
    you can specify the regions at parsing time, rather than re-implementing the subset handling logic here.

    Args:
        bedmethyl_file: Path to the input tabix-indexed gzipped bedmethyl file
        motif: type of modification to extract data for
        bigwig_file: Path to the output bigwig destination. If unspecified, a pileup.bw file will be created in the bedmethyl file's directory
        ref_genome: a reference genome to use for constructing the bigwig header, i.e. contig lengths. If None, the bedmethyl file will be used
            to estimate contig lengths, which can take some time.
        strand: the DNA strand to extra, + or - for forward or reverse and . for both
        chunk_size: size for bigwig write chunks, in bedmethyl lines
    """
    bedmethyl_file, bigwig_file, ref_genome = utils.sanitize_path_args(
        bedmethyl_file, bigwig_file, ref_genome
    )

    # Set up output directories if they don't exist; load up objects for bedmethyl tabix file and motif specifier
    output_file_path = (
        bigwig_file
        if bigwig_file is not None
        else bedmethyl_file.parent / "pileup.fractions.bigwig"
    )
    os.makedirs(output_file_path.parent, exist_ok=True)
    tabix = pysam.TabixFile(str(bedmethyl_file))
    parsed_motif = utils.ParsedMotif(motif)

    # Because we need to set up the bigwig header before we start writing data to it, we need to pre-calculate the length of each contig
    # The header essentially needs to contain a list of the contigs/chromosomes to which the data is aligned, and their sizes.
    # There may be a way to adjust this as we write a bigwig file, but my testing with pyBigWig suggests that you must set it upfront

    contig_lengths_tuples = []
    lines_by_contig = {}

    # If we only have a bedmethyl file, we need to go through it to get contig lengths
    if ref_genome is None:
        for contig in tqdm(
            tabix.contigs,
            desc=f"Step 1: Indexing contigs in {bedmethyl_file.name} to set up bigwig header for {output_file_path.name}",
        ):
            # count up the number of rows, for progress tracking, and pull out the last row so as to grab the length of the chromosome
            # note: the tqdm progress bar slows things down by about 33%, which was deemed better at the time of writing this than
            # 90 seconds without any status updates
            rows_count, last_row = list(
                tail(
                    n=1,
                    iterable=enumerate(
                        tqdm(
                            tabix.fetch(contig),
                            mininterval=1.0,
                            desc=f"Indexing {contig}.",
                            leave=False,
                        )
                    ),
                )
            )[0]
            fields = last_row.split("\t")
            max_coord = int(fields[2])
            contig_lengths_tuples.append((contig, max_coord))
            lines_by_contig[contig] = rows_count
    # If we have a fasta file we can just reference that for contig lengths
    else:
        # Open the reference genome fasta file using pysam
        with pysam.FastaFile(ref_genome) as fasta:
            for contig in tqdm(
                tabix.contigs,  # if these are in the wrong order, e.g. the order from the fasta, it is an issue for pyBigWig somehow
                desc=f"Step 1: Indexing contigs in {ref_genome.name} to set up bigwig header for {output_file_path.name}",
            ):
                # Get the length of the contig
                try:
                    contig_length = fasta.get_reference_length(contig)
                except Exception as err:
                    raise ValueError(
                        f"Error loading {contig} length from {ref_genome.name}. Are you certain that {bedmethyl_file.name} is aligned to this reference?"
                    ) from err
                contig_lengths_tuples.append((contig, contig_length))
                # if we used a fasta to calculate contig lengths we actually don't know the lines per contig
                lines_by_contig[contig] = None

    with pyBigWig.open(str(output_file_path), "w") as bw:
        bw.addHeader(contig_lengths_tuples)
        for contig in tqdm(
            tabix.contigs,
            desc=f"Step 2: Writing {bedmethyl_file.name} contents to {output_file_path.name}",
        ):
            contig_list = []
            start_list = []
            end_list = []
            values_list = []

            for row in tqdm(
                tabix.fetch(contig),
                desc=f"Writing {contig}.",
                total=lines_by_contig[contig],
                leave=False,
            ):
                keep_basemod, genomic_coord, modified_in_row, valid_in_row = (
                    load_processed.process_pileup_row(
                        row=row,
                        parsed_motif=parsed_motif,
                        region_strand=strand,
                        single_strand=(strand != "."),
                    )
                )
                if keep_basemod and valid_in_row > 0:
                    contig_list.append(contig)
                    start_list.append(genomic_coord)
                    end_list.append(genomic_coord + 1)
                    values_list.append(modified_in_row / valid_in_row)

                    if len(values_list) > chunk_size:
                        bw.addEntries(
                            contig_list,  # Contig names
                            start_list,  # Start positions
                            ends=end_list,  # End positions
                            values=values_list,  # Corresponding values
                        )
                        contig_list = []
                        start_list = []
                        end_list = []
                        values_list = []
            bw.addEntries(
                contig_list,  # Contig names
                start_list,  # Start positions
                ends=end_list,  # End positions
                values=values_list,  # Corresponding values
            )
