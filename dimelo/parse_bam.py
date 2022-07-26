r"""
=================
parse_bam module
=================
.. currentmodule:: dimelo.parse_bam
.. autosummary::
    parse_bam

parse_bam allows you to summarize modification calls in a sql database

"""


import argparse
import multiprocessing
import os
import sqlite3
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed
from tqdm import tqdm

from dimelo.utils import clear_db, create_sql_table, execute_sql_command

DEFAULT_BASEMOD = "A+CG"
DEFAULT_THRESH_A = 129
DEFAULT_THRESH_C = 129
DEFAULT_WINDOW_SIZE = 1000


class Region(object):
    def __init__(self, region: Union[str, pd.Series]):
        """Represents a region of genetic data.
        Attributes:
                - chromosome: string name of the chromosome to which the region applies
                - begin: integer start position of region
                - end: integer end position of region
                - size: length of region
                - string: string representation of region
                - strand: string specifying forward or reverse strand; either "+" or "-" (default +)
        """
        self.chromosome = None
        self.begin = None
        self.end = None
        self.size = None
        self.string = None
        self.strand = "+"

        if isinstance(region, str):  # ":" in region:
            # String of format "{CHROMOSOME}:{START}-{END}"
            try:
                self.chromosome, interval = region.replace(",", "").split(":")
                try:
                    # see if just integer chromosomes are used
                    self.chromosome = int(self.chromosome)
                except ValueError:
                    pass
                self.begin, self.end = [int(i) for i in interval.split("-")]
            except ValueError:
                raise TypeError(
                    "Invalid region string. Example of accepted format: 'chr5:150200605-150423790'"
                )
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        elif isinstance(region, pd.Series):
            # Ordered sequence containing [CHROMOSOME, START, END] and optionally [STRAND], where STRAND can be either "+" or "-"
            self.chromosome = region[0]
            self.begin = region[1]
            self.end = region[2]
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
            # strand of motif to orient single molecules
            # if not passed just keep as all +
            if len(region) >= 4:
                if (region[3] == "+") or (region[3] == "-"):
                    self.strand = region[3]
                # handle case of bed file with additional field that isn't strand +/-
                else:
                    self.strand = "+"
            else:
                self.strand = "+"
        else:
            raise TypeError(
                "Unknown datatype passed for Region initialization"
            )


def make_db(
    fileName: str,
    sampleName: str,
    outDir: str,
    testMode: bool = False,
    qc: bool = False,
) -> Tuple[str, List[str]]:
    """Sets up the necessary database tables.

    Args:
            :param fileName: name of bam file with Mm and Ml tags
            :param sampleName: name of sample for output SQL table name labelling
            :param outDir: directory where SQL database is stored
            :param testMode: turns on test mode; note that this will clear the database if it exists
            :param qc: turns on qc mode

    Returns:
            - path to the new database
            - list of newly-created table names
    """
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    DATABASE_NAME = (
        outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
    )

    if testMode:
        clear_db(DATABASE_NAME)

    tables = []
    # for qc report
    if qc:
        table_name = "reads_" + sampleName
        cols = [
            "name",
            "chr",
            "start",
            "end",
            "length",
            "strand",
            "mapq",
            "ave_baseq",
            "ave_alignq",
        ]
        dtypes = [
            "TEXT",
            "TEXT",
            "INT",
            "INT",
            "INT",
            "TEXT",
            "INT",
            "INT",
            "INT",
        ]
        create_sql_table(DATABASE_NAME, table_name, cols, dtypes)
        tables.append(table_name)
    # for browser and enrichment plots
    else:
        table_name = "methylationByBase_" + sampleName
        cols = ["id", "read_name", "chr", "pos", "prob", "mod"]
        dtypes = ["TEXT", "TEXT", "TEXT", "INT", "INT", "TEXT"]
        create_sql_table(DATABASE_NAME, table_name, cols, dtypes)
        tables.append(table_name)

        table_name = "methylationAggregate_" + sampleName
        cols = ["id", "pos", "mod", "methylated_bases", "total_bases"]
        dtypes = ["TEXT", "INT", "TEXT", "INT", "INT"]
        create_sql_table(DATABASE_NAME, table_name, cols, dtypes)
        tables.append(table_name)

    return DATABASE_NAME, tables


def parse_bam(
    fileName: str,
    sampleName: str,
    outDir: str,
    bedFile: str = None,
    basemod: str = DEFAULT_BASEMOD,
    center: bool = False,
    windowSize: int = DEFAULT_WINDOW_SIZE,
    region: str = None,
    threshA: int = DEFAULT_THRESH_A,
    threshC: int = DEFAULT_THRESH_C,
    extractAllBases: bool = False,
    cores: int = None,
) -> None:
    """
    fileName
        name of bam file with Mm and Ml tags
    sampleName
        name of sample for output SQL table name labelling. Valid names contain [``a-zA-Z0-9_``].
    outDir
        directory where SQL database is stored
    bedFile
        name of bed file that defines regions of interest over which to extract mod calls. The bed file either defines regions over which to extract mod calls OR defines regions (likely motifs) over which to center positions for mod calls and then parse_bam extracts mod calls over a window flanking that region defined in by ``windowSize``. Optional 4th column in bed file to specify strand of region of interest as ``+`` or ``-``. Default is to consider regions as all ``+``. NB. The ``bedFile`` and ``region`` parameters are mutually exclusive; specify one or the other.
    basemod
        One of the following:

        * ``'A'`` - extract mA only
        * ``'CG'`` - extract mCpG only
        * ``'A+CG'`` - extract mA and mCpG
    center
        One of the following:

        * ``'True'`` - report positions with respect to center of motif window (+/- windowSize); only valid with bed file input
        * ``'False'`` - report positions in original reference space
    windowSize
        window size around center point of feature of interest to plot (+/-); only mods within this window are stored; only used if center=True; still, only reads that span the regions defined in the bed file will be included; default is 1,000 bp
    region
        single region over which to extract base mods, rather than specifying many windows in bedFile; format is chr:start-end. NB. The ``bedFile`` and ``region`` parameters are mutually exclusive; specify one or the other.
    threshA
        threshold above which to call an A base methylated; default is 129
    threshC
        threshold above which to call a C base methylated; default is 129
    extractAllBases
        One of the following:

        * ``'True'`` - Store all base mod calls, regardles of methylation probability threshold. Bases stored are those that can have a modification call (A, CG, or both depending on ``basemod`` parameter) and are sequenced bases, not all bases in the reference.
        * ``'False'`` - Only modifications above specified threshold are stored
    cores
        number of cores over which to parallelize; default is all available

    Valid argument combinations for ``bedFile``, ``center``, and ``windowSize`` are below. Regions of interest generally fall into two categories: small motifs at which to center analysis (use ``center`` = True) or full windows of interest (do not specify ``center`` or ``windowSize``).

        * ``bedFile`` --> extract all modified bases in regions defined in bed file
        * ``bedfile`` + ``center`` --> extract all modified bases in regions defined in bed file, report positions relative to region centers and extract base modificiations within default windowSize of 1kb
        * ``bedfile`` + ``center`` + ``windowSize`` --> extract all modified bases in regions defined in bed file, report positions relative to region centers and extract base modifications within flanking +/- windowSize
        * ``region`` --> extract all modified bases in single region

    **Example**

    For regions defined by ``bedFile``:

    >>> dm.parse_bam("dimelo/test/data/mod_mappings_subset.bam", "test", "dimelo/dimelo_test", bedFile="dimelo/test/data/test.bed", basemod="A+CG", center=True, windowSize=500, threshA=190, threshC=190, extractAllBases=False, cores=8)

    For single region defined with ``region``:

    >>> dm.parse_bam("dimelo/test/data/mod_mappings_subset.bam", "test", "dimelo/dimelo_test", region="chr1:2907273-2909473", basemod="A+CG", threshA=190, threshC=190, cores=8)

    **Return**

    Returns a SQL database in the specified output directory. Database can be converted into pandas dataframe with:

    >>> fileName = "dimelo/test/data/mod_mappings_subset.bam"
    >>> sampleName = "test"
    >>> outDir = "dimelo/dimelo_test"
    >>> all_data = pd.read_sql("SELECT * from methylationByBase_" + sampleName, sqlite3.connect(outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"))
    >>> aggregate_counts = pd.read_sql("SELECT * from methylationAggregate_" + sampleName, sqlite3.connect(outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"))

    Each database contains these two tables with columns listed below:

    1. methylationByBase_sampleName
        * id(read_name:pos)
        * read_name
        * chr
        * pos
        * prob
        * mod
    2. methylationAggregate_sampleName
        * id(pos:mod)
        * pos
        * mod
        * methylated_bases
        * total_bases

    When running parse_bam with a region defined, a summary bed file is also produced to support visualizing aggregate data with any genome browser tool. The columns of this bed file are chr, start, end, methylated_bases, total_bases.

    For example, to take a summary output bed and create a file with fraction of modified bases with a window size of 100 bp for visualization with the WashU browser, you could run the below commands in terminal:

        * ``bedtools makewindows -g ref_genome.chromsizes.txt -w 100 > ref_genome_windows.100.bp.bed``
        * ``bedtools map -a ref_genome_windows.100.bp.bed -b outDir/fileName_sampleName_chr_start_end_A.bed -c 4,5 -o sum,sum -null 0 | awk -v "OFS=\\t" '{if($5>0){print $1,$2,$3,$4/$5}else{print $1,$2,$3,$5}}' > outDir/fileName_sampleName_chr_start_end_A.100.bed``
        * ``bgzip outDir/fileName_sampleName_chr_start_end_A.100.bed``
        * ``tabix -f -p bed outDir/fileName_sampleName_chr_start_end_A.100.bed.gz``

    """
    # Ensure exactly one of bedFile and region are specified
    if sum([arg is None for arg in (bedFile, region)]) != 1:
        raise RuntimeError(
            "Exactly one of the mutually exclusive arguments 'bedFile' or 'region' must be specified."
        )
    # The argument center is incompatible with region
    if region is not None:
        if center:
            raise RuntimeError(
                "Argument 'center' cannot be given alongside 'region'."
            )

    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    make_db(fileName, sampleName, outDir)

    if bedFile is not None:
        # make a region object for each row of bedFile
        bed = pd.read_csv(bedFile, sep="\t", header=None)
        windows = []
        for _, row in bed.iterrows():
            windows.append(Region(row))

    if region is not None:
        windows = [Region(region)]

    # Configure progress reporting
    if len(windows) == 1:
        show_read_progress = True
    else:
        show_read_progress = False
        # Enable top-level progress bar for multi-window processing
        windows = tqdm(windows, desc="Parsing windows", unit="windows")

    # default number of cores is max available
    cores_avail = multiprocessing.cpu_count()
    if cores is None:
        num_cores = cores_avail
    else:
        # if more than available cores is specified, process with available cores
        if cores > cores_avail:
            num_cores = cores_avail
        else:
            num_cores = cores

    batchSize = 100

    Parallel(n_jobs=num_cores)(
        delayed(parse_reads_window)(
            fileName,
            sampleName,
            basemod,
            windowSize,
            window,
            center,
            threshA,
            threshC,
            batchSize,
            outDir,
            extractAllBases,
            showReadProgress=show_read_progress,
        )
        for window in windows
    )

    # create summary bed files
    if region is not None:
        if "A" in basemod:
            make_bed_file_output(fileName, sampleName, outDir, region, "A")
        if "C" in basemod:
            make_bed_file_output(fileName, sampleName, outDir, region, "C")

    # output all files created to std out
    f = fileName.split("/")[-1].replace(".bam", "")
    out_path = f"{outDir}/{f}.db"
    if region is None:
        str_out = f"Outputs\n_______\nDB file: {out_path}"
    else:
        bed_paths = []
        if "A" in basemod:
            bed_path = (
                f"{outDir}/{f}_{sampleName}_{Region(region).string}_A.bed"
            )
            bed_paths.append(bed_path)
        if "C" in basemod:
            bed_path = (
                f"{outDir}/{f}_{sampleName}_{Region(region).string}_CG.bed"
            )
            bed_paths.append(bed_path)
        str_out = (
            f"Outputs\n_______\nDB file: {out_path}\nBED file: {bed_paths}"
        )
    print(str_out)


def make_bed_file_output(fileName, sampleName, outDir, region, mod):
    """
    Make output bed file that can be used to visualize aggregate with other genome browsers
    """
    r = Region(region)
    f = fileName.split("/")[-1].replace(".bam", "")
    out_path = f"{outDir}/{f}.db"
    aggregate_counts_all = pd.read_sql(
        "SELECT * from methylationAggregate_" + sampleName,
        sqlite3.connect(out_path),
    )
    aggregate_counts_mod = aggregate_counts_all[
        aggregate_counts_all["mod"].str.contains(mod)
    ].copy()
    aggregate_counts_mod["chr"] = r.chromosome
    aggregate_counts_mod["end"] = aggregate_counts_mod["pos"] + 1
    # columns are: id(pos:mod), pos, mod, methylated_bases, total_bases
    dictionary_agg = {
        "chr": aggregate_counts_mod["chr"],
        "start": aggregate_counts_mod["pos"],
        "end": aggregate_counts_mod["end"],
        "methylated": aggregate_counts_mod["methylated_bases"],
        "total": aggregate_counts_mod["total_bases"],
    }
    bed_agg = pd.DataFrame(dictionary_agg)
    bed_agg.sort_values(by="start", ascending=True, inplace=True)
    if "A" in mod:
        mod_name = "A"
    if "C" in mod:
        mod_name = "CG"
    bed_agg.to_csv(
        f"{outDir}/{f}_{sampleName}_{r.string}_{mod_name}.bed",
        sep="\t",
        header=False,
        index=False,
    )


def parse_reads_window(
    fileName: str,
    sampleName: str,
    basemod: str,
    windowSize: int,
    window: Region,
    center: bool,
    threshA: int,
    threshC: int,
    batchSize: int,
    outDir: str,
    extractAllBases: bool,
    showReadProgress: bool = False,
) -> None:
    """Parse all reads in window and put data into methylationByBase table.

    Args:
            :param bam: read in bam file with Mm and Ml tags
            :param fileName: name of bam file
            :param sampleName: name of sample for output file name labelling
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param windowSize: window size around center point of feature of interest to plot (+/-); only mods within this window are stored; only applicable for center=True
            :param window: single window
            :param center: report positions with respect to reference center (+/- window size) if True or in original reference space if False
            :param threshA: threshold above which to call an A base methylated
            :param threshC: threshold above which to call a C base methylated
            :param showReadProgress: when true, display progress for read processing
    """
    bam = pysam.AlignmentFile(fileName, "rb")
    data = []
    if showReadProgress:
        total_reads = bam.count(
            reference=window.chromosome, start=window.begin, end=window.end
        )
        bam.reset()
    reads = bam.fetch(
        reference=window.chromosome, start=window.begin, end=window.end
    )
    if showReadProgress:
        reads = tqdm(
            reads, desc="Processing reads", unit="reads", total=total_reads
        )
    for read in reads:
        [
            (mod, positions, probs),
            (mod2, positions2, probs2),
        ] = get_modified_reference_positions(
            read,
            basemod,
            window,
            center,
            threshA,
            threshC,
            windowSize,
            fileName,
            sampleName,
            outDir,
            extractAllBases,
        )
        # Generate rows for methylationByBase database update
        for pos, prob in zip(positions, probs):
            if pos is not None:
                if (center is True and abs(pos) <= windowSize) or (
                    center is False and pos > window.begin and pos < window.end
                ):  # to decrease memory, only store bases within the window
                    d = (
                        read.query_name + ":" + str(pos),
                        read.query_name,
                        window.chromosome,
                        int(pos),
                        int(prob),
                        mod,
                    )
                    data.append(d)
        for pos, prob in zip(positions2, probs2):
            if pos is not None:
                if (center is True and abs(pos) <= windowSize) or (
                    center is False and pos > window.begin and pos < window.end
                ):  # to decrease memory, only store bases within the window
                    d = (
                        read.query_name + ":" + str(pos),
                        read.query_name,
                        window.chromosome,
                        int(pos),
                        int(prob),
                        mod2,
                    )
                    data.append(d)
    if data:
        # data is list of tuples associated with given read
        # or ignore because a read may overlap multiple windows
        DATABASE_NAME = (
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        )
        table_name = "methylationByBase_" + sampleName
        command = (
            """INSERT OR IGNORE INTO """
            + table_name
            + """ VALUES(?,?,?,?,?,?);"""
        )
        connection = sqlite3.connect(DATABASE_NAME, timeout=60.0)
        execute_sql_command(command, DATABASE_NAME, data, connection)
        connection.close()


def get_modified_reference_positions(
    read: pysam.AlignedSegment,
    basemod: str,
    window: Region,
    center: bool,
    threshA: int,
    threshC: int,
    windowSize: int,
    fileName: str,
    sampleName: str,
    outDir: str,
    extractAllBases: bool,
):
    """Extract mA and mC pos & prob information for the read
    Args:
            :param read: single read from bam file
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param window: window from bed file
            :param center: report positions with respect to reference center (+/- window size) if True or in original reference space if False
            :param threshA: threshold above which to call an A base methylated
            :param threshC: threshold above which to call a C base methylated
            :param windowSize: window size around center point of feature of interest to plot (+/-); only mods within this window are stored; only applicable for center=True
    Return:
        For each mod, you get the positions where those mods are and the probabilities for those mods (parallel vectors)
    """
    if (read.has_tag("Mm")) & (";" in read.get_tag("Mm")):
        mod1 = read.get_tag("Mm").split(";")[0].split(",", 1)[0]
        mod2 = read.get_tag("Mm").split(";")[1].split(",", 1)[0]
        base = basemod[0]  # this will be A, C, or A
        if basemod == "A+CG":
            base2 = basemod[2]  # this will be C for A+CG case
        else:  # in the case of a single mod will just be checking that single base
            base2 = base
        if base in mod1 or base2 in mod1:
            mod1_return = get_mod_reference_positions_by_mod(
                read,
                mod1,
                0,
                window,
                center,
                threshA,
                threshC,
                windowSize,
                fileName,
                sampleName,
                outDir,
                extractAllBases,
            )
        else:
            mod1_return = (None, [None], [None])
        if base in mod2 or base2 in mod2:
            mod2_return = get_mod_reference_positions_by_mod(
                read,
                mod2,
                1,
                window,
                center,
                threshA,
                threshC,
                windowSize,
                fileName,
                sampleName,
                outDir,
                extractAllBases,
            )
            return (mod1_return, mod2_return)
        else:
            return (mod1_return, (None, [None], [None]))
    else:
        return ((None, [None], [None]), (None, [None], [None]))


def get_mod_reference_positions_by_mod(
    read: pysam.AlignedSegment,
    basemod: str,
    index: int,
    window: Region,
    center: bool,
    threshA: int,
    threshC: int,
    windowSize: int,
    fileName: str,
    sampleName: str,
    outDir: str,
    extractAllBases: bool,
):
    """Get positions and probabilities of modified bases for a single read
    Args:
            :param read: one read in bam file
            :param mod: which basemod, reported as base+x/y/m
            :param window: window from bed file
            :param center: report positions with respect to reference center (+/- window size) if True or in original reference space if False
            :param threshA: threshold above which to call an A base methylated
            :param threshC: threshold above which to call a C base methylated
            :param windowSize: window size around center point of feature of interest to plot (+/-); only mods within this window are stored; only applicable for center=True
            :param index: 0 or 1
    """
    modsPresent = True
    base, mod = basemod.split("+")
    num_base = len(read.get_tag("Mm").split(";")[index].split(",")) - 1
    # get base_index
    base_index = np.array(
        [
            i
            for i, letter in enumerate(read.get_forward_sequence())
            if letter == base
        ]
    )
    # get reference positons
    refpos = np.array(read.get_reference_positions(full_length=True))
    if read.is_reverse:
        refpos = np.flipud(refpos)
    modified_bases = []
    if num_base == 0:
        modsPresent = False
    if modsPresent:
        deltas = [
            int(i) for i in read.get_tag("Mm").split(";")[index].split(",")[1:]
        ]
        Ml = read.get_tag("Ml")
        if index == 0:
            probabilities = np.array(Ml[0:num_base], dtype=int)
        if index == 1:
            probabilities = np.array(Ml[0 - num_base :], dtype=int)
        # determine locations of the modified bases, where index_adj is the adjustment of the base_index
        # based on the cumulative sum of the deltas
        locations = np.cumsum(deltas)
        # loop through locations and increment index_adj by the difference between the next location and current one + 1
        # if the difference is zero, therefore, the index adjustment just gets incremented by one because no base should be skipped
        index_adj = []
        index_adj.append(locations[0])
        i = 0
        for i in range(len(locations) - 1):
            diff = locations[i + 1] - locations[i]
            index_adj.append(index_adj[i] + diff + 1)
        # get the indices of the modified bases
        modified_bases = base_index[index_adj]

    # extract CpG sites only rather than all mC
    keep = []
    prob_keep = []
    all_bases_index = []
    probs = []
    i = 0
    seq = read.get_forward_sequence()
    # deal with None for refpos from soft clipped / unaligned bases
    if "C" in basemod:
        for b in base_index:
            if (
                b < len(seq) - 1
            ):  # if modified C is not the last base in the read
                if (refpos[b] is not None) & (refpos[b + 1] is not None):
                    if seq[b + 1] == "G":
                        if (
                            abs(refpos[b + 1] - refpos[b]) == 1
                        ):  # ensure there isn't a gap
                            all_bases_index.append(
                                b
                            )  # add to all_bases_index whether or not modified
                            if b in modified_bases:
                                if probabilities[i] >= threshC:
                                    keep.append(b)
                                    prob_keep.append(i)
                            if extractAllBases:
                                if b in modified_bases:
                                    probs.append(probabilities[i])
                                else:
                                    probs.append(0)
            # increment for each instance of modified base
            if b in modified_bases:
                i = i + 1
    else:  # for m6A no need to look at neighboring base; do need to remove refpos that are None
        for b in base_index:
            if refpos[b] is not None:
                all_bases_index.append(
                    b
                )  # add to all_bases_index whether or not modified
                if b in modified_bases:
                    if probabilities[i] >= threshA:
                        keep.append(b)
                        prob_keep.append(i)
                if extractAllBases:
                    if b in modified_bases:
                        probs.append(probabilities[i])
                    else:
                        probs.append(0)
            # increment for each instance of modified base
            if b in modified_bases:
                i = i + 1
    # adjust position to be centered at 0 at the center of the motif; round in case is at 0.5
    # add returning base_index for plotting mod/base_abundance
    if center is True:
        if window.strand == "+":
            refpos_mod_adjusted = np.array(refpos[keep]) - round(
                ((window.end - window.begin) / 2 + window.begin)
            )
            refpos_total_adjusted = np.array(refpos[all_bases_index]) - round(
                ((window.end - window.begin) / 2 + window.begin)
            )
        if window.strand == "-":
            refpos_mod_adjusted = -1 * (
                np.array(refpos[keep])
                - round(((window.end - window.begin) / 2 + window.begin))
            )
            refpos_total_adjusted = -1 * (
                np.array(refpos[all_bases_index])
                - round(((window.end - window.begin) / 2 + window.begin))
            )
        update_methylation_aggregate_db(
            refpos_mod_adjusted,
            refpos_total_adjusted,
            basemod,
            center,
            windowSize,
            window,
            fileName,
            sampleName,
            outDir,
        )
        if extractAllBases:
            return (basemod, refpos_total_adjusted, probs)
        elif not modsPresent:
            return (None, [None], [None])
        else:
            return (basemod, refpos_mod_adjusted, probabilities[prob_keep])
    else:
        update_methylation_aggregate_db(
            refpos[keep],
            refpos[all_bases_index],
            basemod,
            center,
            windowSize,
            window,
            fileName,
            sampleName,
            outDir,
        )
        if extractAllBases:
            return (basemod, np.array(refpos[all_bases_index]), probs)
        elif not modsPresent:
            return (None, [None], [None])
        else:
            return (basemod, np.array(refpos[keep]), probabilities[prob_keep])


def update_methylation_aggregate_db(
    refpos_mod: np.ndarray,
    refpos_total: np.ndarray,
    basemod: str,
    center: bool,
    windowSize: int,
    window: Region,
    fileName: str,
    sampleName: str,
    outDir: str,
) -> None:
    """Updates the aggregate methylation table with all of the methylation information from a single read.
    Args:
            :param refpos_mod: list of modified reference positions
            :param refpos_total: list of all reference positions for the base in question
    df with columns pos:modification, pos, mod, methylated_bases, total_bases
    """
    # store list of entries for a given read
    data = []
    for pos in refpos_total:
        # only store positions within window
        if (center is True and abs(pos) <= windowSize) or (
            center is False and pos > window.begin and pos < window.end
        ):
            # key is pos:mod
            id = str(pos) + ":" + basemod
            if pos in refpos_mod:
                data.append((id, int(pos), basemod, 1, 1))
            else:
                data.append((id, int(pos), basemod, 0, 1))

    if data:  # if data to append is not empty
        DATABASE_NAME = (
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        )
        # set variables for sqlite entry
        table_name = "methylationAggregate_" + sampleName

        # create or ignore if key already exists
        # need to add 0 filler here so later is not incremented during update command
        command = (
            """INSERT OR IGNORE INTO """
            + table_name
            + """ VALUES(?,?,?,?,?);"""
        )

        data_fill = [(x[0], x[1], x[2], 0, 0) for x in data]
        connection = sqlite3.connect(DATABASE_NAME, timeout=60.0)
        execute_sql_command(command, DATABASE_NAME, data_fill, connection)
        connection.close()

        # update table for all entries
        # values: methylated_bases, total_bases, id
        # these are entries 3, 4, 0 in list of tuples
        values_subset = [(x[3], x[4], x[0]) for x in data]
        command = (
            """UPDATE """
            + table_name
            + """ SET methylated_bases = methylated_bases + ?, total_bases = total_bases + ? WHERE id = ?"""
        )
        connection = sqlite3.connect(DATABASE_NAME, timeout=60.0)
        execute_sql_command(command, DATABASE_NAME, values_subset, connection)
        connection.close()


def main():
    parser = argparse.ArgumentParser(
        description="Parse a bam file into DiMeLo database tables"
    )

    # Required arguments
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-f",
        "--fileName",
        required=True,
        help="name of bam file with Mm and Ml tags",
    )
    required_args.add_argument(
        "-s",
        "--sampleName",
        required=True,
        help="name of sample for output SQL table name labelling",
    )
    required_args.add_argument(
        "-o",
        "--outDir",
        required=True,
        help="directory where SQL database is stored",
    )

    # Required, mutually exclusive arguments
    window_group = parser.add_mutually_exclusive_group(required=True)
    window_group.add_argument(
        "-b",
        "--bedFile",
        help="name of bed file that defines regions of interest over which to extract mod calls",
    )
    window_group.add_argument(
        "-r",
        "--region",
        help='single region over which to extract base mods, e.g. "chr1:1-100000"',
    )

    # Optional arguments
    parser.add_argument(
        "-m",
        "--basemod",
        type=str,
        default=DEFAULT_BASEMOD,
        choices=["A", "CG", "A+CG"],
        help="which base modifications to extract",
    )
    parser.add_argument(
        "-A",
        "--threshA",
        type=int,
        default=DEFAULT_THRESH_A,
        help="threshold above which to call an A base methylated",
    )
    parser.add_argument(
        "-C",
        "--threshC",
        type=int,
        default=DEFAULT_THRESH_C,
        help="threshold above which to call a C base methylated",
    )
    parser.add_argument(
        "-e",
        "--extractAllBases",
        action="store_true",
        help="store all base mod calls, regardless of methylation probability threshold",
    )
    parser.add_argument(
        "-p",
        "--cores",
        type=int,
        help="number of cores over which to parallelize",
    )
    parser.add_argument(
        "-c",
        "--center",
        action="store_true",
        help="report positions with respect to center of motif window; only valid with bed file input",
    )
    parser.add_argument(
        "-w",
        "--windowSize",
        type=int,
        default=DEFAULT_WINDOW_SIZE,
        help=f"window size around center point of feature of interest to plot (+/-); only mods within this window are stored (default: {DEFAULT_WINDOW_SIZE} bp)",
    )

    args = parser.parse_args()
    parse_bam(**vars(args))
