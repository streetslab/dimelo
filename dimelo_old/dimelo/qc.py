r"""
=================
QC module
=================
.. currentmodule:: dimelo.qc_report
.. autosummary::
    qc_report

qc_report provides a detailed summary report of many important quality control information including read length, mapping quality, etc.

"""
import argparse
import multiprocessing
import os
import sqlite3
import time
from math import log

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from joblib import Parallel, delayed
from tqdm import tqdm

from dimelo.parse_bam import make_db
from dimelo.utils import execute_sql_command

DEFAULT_COLOR_LIST = [
    "#BB4430",
    "#FFBC0A",
    "#053C5E",
    "#A9E5BB",
    "#610345",
    "#2D1E2F",
    "#559CAD",
    "#5E747F",
    "#F343F4",
]


def batch_read_generator(file_bamIn, filename):
    counter = 0
    r_list = []

    lines = pysam.idxstats(filename).splitlines()
    total_reads = sum(
        [
            int(line.split("\t")[2])
            for line in lines
            if not line.startswith("#")
        ]
    )
    batch_size = 0.1 * total_reads

    for read in file_bamIn.fetch(until_eof=True):
        r = [
            read.query_name,
            read.reference_name,
            read.reference_start,
            read.reference_end,
            read.query_length,
            "-" if read.is_reverse else "+",
            read.mapping_quality,
            ave_qual(read.query_qualities),
            ave_qual(read.query_alignment_qualities),
        ]

        r = tuple(r)
        if counter < batch_size:
            r_list.append(r)
            counter += 1
        else:
            yield r_list
            counter = 0
            r_list = [r]
    yield r_list


def logger(statement):
    print(statement)


def prob_bin(bin):
    # probability a base in the window (or across reads or across bases within a read) is methylated by:
    # calculating probability that no base in the window (or across reads) is methylated and then taking the complement
    # treat p=1 as 254/255 for prevent log(0)
    probs = [
        np.log(1 - p) for p in bin if ((p < 1) and (p >= 0.5))
    ]  # only consider probabilities > 0.5 and handle 1 on next line
    probs1 = [np.log(1 - 254 / 255) for p in bin if p == 1]
    probsAll = probs + probs1
    prob = 1 - np.exp(sum(probsAll))
    return prob


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10 ** (q / -10) for q in range(n + 1)]


def ave_qual(quals, qround=False, tab=errs_tab(129)):
    """Calculate average basecall quality of a read.
    Receive the integer quality scores of a read and return the average quality for that read
    First convert Phred scores to probabilities,
    calculate average error probability
    convert average back to Phred scale
    """
    if quals:
        mq = -10 * log(sum([tab[q] for q in quals]) / len(quals), 10)
        if qround:
            return round(mq)
        else:
            return mq
    else:
        return None


def parse_bam_read(bamIn, sampleName, outDir, cores=None):
    file_bamIn = pysam.AlignmentFile(bamIn, "rb")

    DB_NAME, tables = make_db(bamIn, sampleName, outDir, qc=True)
    template_command = (
        """INSERT INTO """ + tables[0] + """ VALUES(?,?,?,?,?,?,?,?,?);"""
    )
    connect = sqlite3.connect(DB_NAME, timeout=60.0, check_same_thread=False)
    cores_avail = multiprocessing.cpu_count()
    if cores is None:
        num_cores = cores_avail
    else:
        # if more than available cores is specified, process with available cores
        if cores > cores_avail:
            num_cores = cores_avail
        else:
            num_cores = cores

    c = connect.cursor()
    c.execute("BEGIN TRANSACTION")

    Parallel(n_jobs=num_cores, backend="threading")(
        delayed(execute_sql_command)(template_command, DB_NAME, i, connect)
        for i in tqdm(
            batch_read_generator(file_bamIn, bamIn),
            total=10,
            desc="Processing reads",
            unit=" batches",
        )
    )

    c.close()
    connect.close()
    return DB_NAME, tables[0]


def get_runtime(f, inp1, inp2, inp3):
    start = time.time()
    re_val = f(inp1, inp2, inp3)
    time.sleep(1)
    end = time.time()
    return f"Runtime of the program is {end - start}", re_val


def qc_plot(x, sampleName, plotType, colors, num, axes):
    an_array = np.array(x)
    if all(v is None for v in an_array):
        return []
    q1 = np.quantile(an_array, 0.25)
    q3 = np.quantile(an_array, 0.75)
    iq = q3 - q1
    outlier = q3 + 3 * iq
    not_outlier = an_array <= outlier
    no_outliers = an_array[not_outlier]

    ptype = ""
    unit = ""
    xlabel = ""
    hasN50 = False
    if plotType == "L":
        ptype = " Read Length"
        xlabel = "Read Length (bp)"
        unit = " bp"
        hasN50 = True
        n50 = calculate_N50(x)
    elif plotType == "M":
        ptype = " Mapping Quality"
        xlabel = "Mapping Quality"
    elif plotType == "B":
        ptype = " Basecall Quality"
        xlabel = "Average Basecall Quality"
    elif plotType == "A":
        ptype = " Alignment Quality"
        xlabel = "Average Alignment Quality"
    plt.hist(no_outliers, bins=200, color=colors[6], density=True)  #
    plt.axvline(
        x.median(),
        color=colors[0],
        linestyle="dashed",
        linewidth=1.3,
        label="median: " + str(round(x.median())) + unit,
    )
    plt.axvline(
        x.mean(),
        color=colors[2],
        linestyle="dashed",
        linewidth=1.3,
        label="mean: " + str(round(x.mean())) + unit,
    )
    if hasN50:
        plt.axvline(
            n50,
            color=colors[1],
            linestyle="dashed",
            linewidth=1,
            label="N50: " + str(round(n50)) + unit,
        )
    plt.title(ptype)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.xlim(
        0,
    )
    plt.plot([], [], " ", label="max: " + str(round(max(x))) + unit)
    plt.legend()

    values = [
        round(min(x)),
        round(q1),
        round(x.median()),
        round(q3),
        round(max(x)),
        round(x.mean()),
    ]
    return values


def calculate_N50(x):
    array_rl = np.array(x)
    N = np.sum(array_rl)
    array_rl[::-1].sort()
    rl_cumsum = np.cumsum(array_rl)
    n50 = array_rl[np.argmax(rl_cumsum > N / 2)]
    return n50


def qc_report(
    fileNames,
    sampleNames,
    outDir,
    colors=DEFAULT_COLOR_LIST,
    cores=None,
):

    """
    fileNames
        list of names of bam files; indexed; or single file name as string
    sampleNames
        list of names of samples for output plot name labelling; or single sample name as string; valid names contain [``a-zA-Z0-9_``].
    outDir
        directory to output QC summary report
    cores
        number of cores over which to parallelize; default is all available
    colors
        color list in hex for overlay plots; default is:
        ["#BB4430","#FFBC0A","#053C5E","#A9E5BB","#610345",
        "#2D1E2F","#559CAD","#5E747F","#F343F4"]

    **Example**

    For single sample:

        >>> dm.qc_report("dimelo/test/data/mod_mappings_subset.bam", "test", "dimelo/dimelo_test")

    For multiple sample files:

        >>> dm.qc_report(["dimelo/test/data/mod_mappings_subset.bam", "dimelo/test/data/winnowmap_guppy_merge_subset.bam"], ["test1", "test2"], "dimelo/dimelo_test")

    **Return**

        * PDF of QC Summary Report which includes:
            * read length histogram
            * mapping quality histogram
            * average alignment quality per read histogram (if basecaller provided information)
            * average basecall quality per read histogram (if basecaller provided information)
            * summary table describing spread of data
            * number of reads, number of basepairs

    Returns a SQL database in the specified output directory. Database can be converted into pandas dataframe with:

    >>> fileName = "dimelo/test/data/mod_mappings_subset.bam"
    >>> sampleName = "test"
    >>> outDir = "dimelo/dimelo_test"
    >>> all_reads = pd.read_sql("SELECT * from reads_" + sampleName, sqlite3.connect(outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"))


    After QC, each database contains this table with columns listed below:

    reads_sampleName
        * name
        * chr
        * start
        * end
        * length
        * strand
        * mapq
        * ave_baseq
        * ave_alignq

    **Example Plots**
    :ref:`sphx_glr_auto_examples_plot_qc_example.py`

    """
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    if type(fileNames) != list:
        fileNames = [fileNames]
        sampleNames = [sampleNames]

    for index in range(len(fileNames)):
        filebamIn = fileNames[index]
        sampleName = sampleNames[index]
        DB_NAME, TABLE_NAME = parse_bam_read(
            filebamIn, sampleName, outDir, cores
        )

        if sampleName is None:
            sampleName = DB_NAME.split("/")[-1][:-3]

        plot_feature_df = pd.read_sql(
            "SELECT * from " + TABLE_NAME, con=sqlite3.connect(DB_NAME)
        )

        fig = plt.figure(figsize=(12, 10))
        grid = plt.GridSpec(3, 2, figure=fig)

        ax_5 = plt.subplot2grid(shape=(3, 2), loc=(0, 0), colspan=2)
        ax_5.axis("off")

        keep_values = [0, 0, 0, 0]
        valRL = []
        valMQ = []
        valBQ = []
        valAQ = []

        # Read Length
        x = plot_feature_df["length"]
        ax_1 = fig.add_subplot(grid[0, 0])
        valRL = qc_plot(x, sampleName, "L", colors, 1, ax_1)
        if valRL:
            keep_values[0] = 1

        # Mapping Quality
        x = plot_feature_df["mapq"]
        ax_2 = fig.add_subplot(grid[0, 1])
        valMQ = qc_plot(x, sampleName, "M", colors, 2, ax_2)
        if valMQ:
            keep_values[1] = 1

        # Basecall Quality
        x = plot_feature_df["ave_baseq"]
        ax_3 = fig.add_subplot(grid[1, 0])
        valBQ = qc_plot(x, sampleName, "B", colors, 3, ax_3)
        if valBQ:
            keep_values[2] = 1

        # Alignment Quality
        x = plot_feature_df["ave_alignq"]
        ax_4 = fig.add_subplot(grid[1, 1])
        valAQ = qc_plot(x, sampleName, "A", colors, 4, ax_4)
        if valAQ:
            keep_values[3] = 1

        val_table = [valRL, valMQ, valBQ, valAQ]
        val_table_new = []
        cols = [
            "Read Length",
            "Mapping Quality",
            "Basecall Quality",
            "Alignment Quality",
        ]
        columns = []
        axes_stored = [ax_1, ax_2, ax_3, ax_4]
        for i in range(len(keep_values)):
            if keep_values[i] == 1:
                val_table_new.append(val_table[i])
                columns.append(cols[i])
            else:
                plt.delaxes(axes_stored[i])

        report_table = np.array(val_table_new).T

        rows = ["Min", "25%", "Median", "75%", "Max", "Mean"]
        print("mean length: ", valRL[5])
        print("num reads: ", len(x))
        print("num bases: ", round(valRL[5] * len(x)))

        if len(columns) <= 2:
            ax_5 = plt.subplot2grid(shape=(3, 2), loc=(1, 0), colspan=2)
        else:
            ax_5 = plt.subplot2grid(shape=(3, 2), loc=(2, 0), colspan=2)
        ax_5.axis("off")
        ax_5.table(
            cellText=report_table,
            rowLabels=rows,
            colLabels=columns,
            loc="center",
        )

        fig.tight_layout(w_pad=2, h_pad=4)

        summary_data = "mean length: " + str(valRL[5]) + " bp"
        summary_data = summary_data + "; num reads: " + str(len(x))
        summary_data = (
            summary_data + "; " + "num bases: " + str(round(valRL[5] * len(x)))
        )
        fig.suptitle(sampleName + " QC Summary Report", y=1.05)

        plt.title(summary_data, y=0.8)

        # saving as PDF
        final_file_name = outDir + "/" + sampleName + "_qc_report"
        plt.savefig(final_file_name + ".pdf", bbox_inches="tight")
        plt.close()

        print("QC report located at: " + final_file_name + ".pdf")
        print("Database located at: " + DB_NAME)


def main():
    parser = argparse.ArgumentParser(description="Generate DiMeLo qc report")

    # Required arguments
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-f", "--fileNames", required=True, nargs="+", help="bam file name(s)"
    )
    required_args.add_argument(
        "-s",
        "--sampleNames",
        required=True,
        nargs="+",
        help="sample name(s) for output labelling",
    )
    required_args.add_argument(
        "-o",
        "--outDir",
        required=True,
        help="directory to output QC summary report",
    )

    # Plotting arguments
    plotting_args = parser.add_argument_group("plotting options")
    plotting_args.add_argument(
        "--colors",
        type=str,
        nargs="+",
        default=DEFAULT_COLOR_LIST,
        help='color list in hex (e.g. "#BB4430") for overlay plots',
    )

    # Optional arguments
    parser.add_argument(
        "-p",
        "--cores",
        type=int,
        help="number of cores over which to parallelize",
    )

    args = parser.parse_args()
    qc_report(**vars(args))
