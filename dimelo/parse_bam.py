r"""
=================
parse_bam module
=================
.. currentmodule:: dimelo.parse_bam
.. autosummary::
    parse_bam

parse_bam allows you to summarize modification calls in a sql database

"""

import multiprocessing
import os
from typing import List, Tuple

import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed

from dimelo.utils import clear_db, create_sql_table, execute_sql_command


"""
TODO:
    - convert paths over to pathlib where possible, and update relevant type annotations
    - Standardize comment formats
    - convert string concatenations to f-strings, to make things more readable
"""


class Region(object):
    def __init__(self, region):
        """
        TODO:
                - document this
                - it appears that the region parameter can be either a string of some sort or a row of a bed file, as returned from pd.DataFrame.iterrows()
        """
        if isinstance(region, str):  # ":" in region:
            try:
                self.chromosome, interval = region.replace(",", "").split(":")
                try:
                    # see if just integer chromosomes are used
                    self.chromosome = int(self.chromosome)
                except ValueError:
                    pass
                self.begin, self.end = [int(i) for i in interval.split("-")]
            except ValueError:
                pass
                # sys.exit(
                #    "\n\nERROR: Window (-w/--window) inproperly formatted, "
                #    "examples of accepted formats are:\n"
                #    "'chr5:150200605-150423790'\n\n"
                # )
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        else:
            self.chromosome = region[1][0]
            self.begin = region[1][1]
            self.end = region[1][2]
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
            # strand of motif to orient single molecules
            # if not passed just keep as all +
            if len(region[1]) >= 4:
                if (region[1][3] == "+") or (region[1][3] == "-"):
                    self.strand = region[1][3]
                # handle case of bed file with additional field that isn't strand +/-
                else:
                    self.strand = "+"
            else:
                self.strand = "+"


def make_db(
        fileName: str,
        sampleName: str,
        outDir: str,
        testMode: bool=False,
        qc: bool=False,
        joint: bool=False
) -> Tuple[str, List[str]]:
    """Sets up the necessary database tables.

    Args:
            :param fileName: name of bam file with Mm and Ml tags
            :param sampleName: name of sample for output SQL table name labelling
            :param outDir: directory where SQL database is stored
            :param testMode: turns on test mode; note that this will clear the database if it exists
            :param qc: turns on qc mode
            :param joint: turns on joint enrichment mode

    Returns:
            - path to the new database
            - list of newly-created table names
    
    TODO:
            - document testMode, qc, and joint mode more fully; update top-level documentation accordingly
            - make this use pathlib for path operations
            - current uses:
              parse_bam.py:234:    make_db(fileName, sampleName, outDir, testMode, qc, joint)
              plot_joint_enrichment.py:240:    make_db(fileName, sampleName, outDir, joint=True)
              qc.py:96:    DB_NAME, tables = make_db(bamIn, "", outDir, qc=True)
            - ****Right now, joint and qc are basically mutually exclusive. Is this functionality intended? Why should it be so?
            - Modularize and simplify sql table specification. There should be easy-to-reference top-level variables that contain the table names, columns, datatypes, etc. Not sure how they should be stored yet, however, because I do not know if/how they are going to be referenced elsewhere in the codebase.
    """
    # TODO: These should replace the path operations once pathlib conversion is in place
    # filePath = Path(fileName)
    # outPath = Path(outDir)
    
    # outPath.mkdir(parents=True, exist_ok=True)

    # (outPath / filePath.name).withsuffix('.db')
    
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
        table_name = "reads"
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
    # for joint occupancy plots
    elif joint:
        table_name = "methylationByBaseJoint_" + sampleName
        cols = [
            "id",
            "read_windows",
            "read_name",
            "pos",
            "prob",
            "mod",
            "peak_strength",
        ]
        dtypes = ["TEXT", "TEXT", "TEXT", "INT", "INT", "TEXT", "FLOAT"]
        create_sql_table(DATABASE_NAME, table_name, cols, dtypes)
        tables.append(table_name)

        table_name = "methylationAggregate_" + sampleName
        cols = ["id", "pos", "mod", "methylated_bases", "total_bases"]
        dtypes = ["TEXT", "INT", "TEXT", "INT", "INT"]
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
    bedFile: str=None,
    basemod: str="A+CG",
    center: bool=False,
    windowSize: int=None,
    region: str=None,
    threshA: int=129,
    threshC: int=129,
    extractAllBases: bool=False,
    testMode: bool=False,
    qc: bool=False,
    joint: bool=False,
    cores: int=None
) -> None:
    """
    fileName
        name of bam file with Mm and Ml tags
    sampleName
        name of sample for output SQL table name labelling
    outDir
        directory where SQL database is stored
    bedFile
        name of bed file that defines regions of interest over which to extract mod calls
    basemod
        One of the following:

        * ``'A'`` - extract mA only
        * ``'CG'`` - extract mCpG only
        * ``'A+CG'`` - extract mA and mCpG
    center
        One of the following:

        * ``'True'`` - report positions with respect to reference center (+/- window size)
        * ``'False'`` - report positions in original reference space
    windowSize
        window size around center point of feature of interest to plot (+/-); only mods within this window are stored; only specify if center=True
    region
        single region over which to extract base mods, rather than specifying many windows in bedFile; format is chr:start-end
    threshA
        threshold above which to call an A base methylated; default is 129
    threshC
        threshold above which to call a C base methylated; default is 129
    extractAllBases
         One of the following:

        * ``'True'`` - store all base mod calls, regardles of methylation probability threshold
        * ``'False'`` - only modifications above specified threshold are stored
    cores
        number of cores over which to parallelize; default is all available

    **Example**

    >>> dm.parse_bam("dimelo/test/data/mod_mappings_subset.bam", "test", "/dimelo/dimelo_test", bedFile="dimelo/test/data/test.bed", basemod="A+CG", center=True, windowSize=500, threshA=190, threshC=190, extractAllBases=False, cores=8)

    **Return**

    Returns a SQL database in the specified output directory. Database can be converted into pandas dataframe with:

    >>> all_data = pd.read_sql("SELECT * from methylationByBase_" + sampleName, sqlite3.connect(outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"))
    >>> aggregate_counts = pd.read_sql("SELECT * from methylationAggregate_" + sampleName, sqlite3.connect(outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"))

    """
    make_db(fileName, sampleName, outDir, testMode, qc, joint)

    if bedFile is not None:
        # make a region object for each row of bedFile
        bed = pd.read_csv(bedFile, sep="\t", header=None)
        windows = []
        for row in bed.iterrows():
            windows.append(Region(row))

    if region is not None:
        windows = [region]

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

    Parallel(n_jobs=num_cores, verbose=10)(
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
            qc,
        )
        for window in windows
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
    qc: bool,
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

    TODO:
            - Find a way to mention in documentation that this has a side effect of populating the aggregated table as well...
    """
    bam = pysam.AlignmentFile(fileName, "rb")
    data = []
    for read in bam.fetch(
        reference=window.chromosome, start=window.begin, end=window.end
    ):
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
            qc,
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
        execute_sql_command(command, DATABASE_NAME, data)


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
    qc: bool,
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
    TODO:
            - This is referenced in plot_joint_enrichment
            - Oh boy, what does this return? At minimum it should have a type annotation. At maximum, it should have a class.
    """
    if (read.has_tag("Mm")) & (";" in read.get_tag("Mm")):
        mod1 = read.get_tag("Mm").split(";")[0].split(",", 1)[0]
        mod2 = read.get_tag("Mm").split(";")[1].split(",", 1)[0]
        # mod1_list = read.get_tag("Mm").split(";")[0].split(",", 1)
        # mod2_list = read.get_tag("Mm").split(";")[1].split(",", 1)
        base = basemod[0]  # this will be A, C, or A
        if basemod == "A+CG":
            base2 = basemod[2]  # this will be C for A+CG case
        else:  # in the case of a single mod will just be checking that single base
            base2 = base
        # if len(mod1_list) > 1 and (base in mod1 or base2 in mod1):
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
                qc,
            )
        else:
            mod1_return = (None, [None], [None])
        # if len(mod2_list) > 1 and (base in mod2 or base2 in mod2):
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
                qc,
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
    qc: bool,
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

    TODO:
            - What is index? What is its type? What does it represent?
            - Pretty sure this is where I'd actually have to update things for GpC methylation? Maybe elsewhere too though?
            - Oh boy, what does this return? At minimum it should have a type annotation. At maximum, it should have a class.
                -  -> Tuple(str, List[int], List[int])
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
                # TODO: For GpC, I want b-1 to be a G, but don't want b+1 to be a C
                # TODO: need to check first position edge case
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
    elif not qc:
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
    else:
        if not modsPresent:
            return (basemod, len(np.array(refpos[all_bases_index])), [])
        else:
            return (
                basemod,
                len(np.array(refpos[all_bases_index])),
                probabilities[prob_keep],
            )


def update_methylation_aggregate_db(
    refpos_mod: List[int],
    refpos_total: List[int],
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

    TODO:
       - verify the exact types of the refpos arguments
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
        execute_sql_command(command, DATABASE_NAME, data_fill)

        # update table for all entries
        # values: methylated_bases, total_bases, id
        # these are entries 3, 4, 0 in list of tuples
        values_subset = [(x[3], x[4], x[0]) for x in data]
        command = (
            """UPDATE """
            + table_name
            + """ SET methylated_bases = methylated_bases + ?, total_bases = total_bases + ? WHERE id = ?"""
        )
        execute_sql_command(command, DATABASE_NAME, values_subset)
