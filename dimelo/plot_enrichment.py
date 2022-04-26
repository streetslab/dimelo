r"""
=======================
plot_enrichment module
=======================
.. currentmodule:: dimelo.plot_enrichment
.. autosummary::
    plot_enrichment

plot_enrichment plots fraction of bases modified within regions of interest defined by bed file

"""

import multiprocessing
import os
import sqlite3
import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from dimelo.parse_bam import parse_bam


DEFAULT_THRESH_A = 129
DEFAULT_THRESH_C = 129
DEFAULT_COLOR_LIST = ["#2D1E2F", "#A9E5BB", "#610345", "#559CAD", "#5E747F"]


def plot_enrichment(
    fileNames,
    sampleNames,
    bedFiles,
    basemod,
    outDir,
    threshA=DEFAULT_THRESH_A,
    threshC=DEFAULT_THRESH_C,
    colors=DEFAULT_COLOR_LIST,
    cores=None,
):
    """
    fileNames
        name(s) of bam file with Mm and Ml tags
    sampleNames
        name(s) of sample for output file name labelling
    bedFiles
        specified windows for region(s) of interest
    basemod
        One of the following (only valid to look at one type of mod):

        * ``'A'`` - extract mA only
        * ``'CG'`` - extract mCpG only
    outDir
        directory to output plot
    threshA
        threshold for calling mA; default 129
    threshC
        threshold for calling mCG; default 129
    colors
        color list in hex for overlay; default is ["#2D1E2F", "#A9E5BB", "#610345", "#559CAD", "#5E747F"]
    cores
        number of cores over which to parallelize; default is all available

    **Example**

    >>> dm.plot_enrichment(["dimelo/test/data/mod_mappings_subset.bam", "dimelo/test/data/mod_mappings_subset.bam"], ["test1", "test2"], "dimelo/test/data/test.bed", "CG", "dimelo/dimelo_test", threshC=129)
    >>> dm.plot_enrichment("dimelo/test/data/mod_mappings_subset.bam", ["test1", "test2"], ["dimelo/test/data/test.bed", "dimelo/test/data/test.bed"], "CG", "dimelo/dimelo_test", threshC=129)

    **Return**

    Barplot with overall fraction of bases modified within regions of interest specified by bedFile(s)

    """
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

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

    # A+CG is not valid; only valid to look at one type of mod
    if (basemod != "A") and (basemod != "CG"):
        print("valid basemods are A or CG")
        return

    # if  single bam file rather than list is entered, convert to list
    if type(fileNames) != list:
        fileNames = [fileNames]
    # if single sample name rather than list is entered, convert to list
    if type(sampleNames) != list:
        sampleNames = [sampleNames]
    # if single bed file rather than list is entered, convert to list
    if type(bedFiles) != list:
        bedFiles = [bedFiles]

    # extract counts and create barplots
    # get average across all bases for regions defined in the bed file
    columns = ["fileName", "bedFile", "sampleName", "fractionMethylated"]
    data = []
    if len(fileNames) > 1 or len(bedFiles) > 1:
        if len(fileNames) > 1:
            if len(bedFiles) > 1:
                print(
                    "only a single region file can be used when analyzing multiple bam files"
                )
                return
            for f, n in zip(fileNames, sampleNames):
                values = get_counts(
                    f,
                    n,
                    bedFiles[0],
                    basemod,
                    outDir,
                    threshA,
                    threshC,
                    num_cores,
                )
                zipped = zip(columns, values)
                a_dictionary = dict(zipped)
                data.append(a_dictionary)
        if len(bedFiles) > 1:
            if len(fileNames) > 1:
                print(
                    "only a single bam file can be used when analyzing multiple bed file regions"
                )
                return
            for b, n in zip(bedFiles, sampleNames):
                values = get_counts(
                    fileNames[0],
                    n,
                    b,
                    basemod,
                    outDir,
                    threshA,
                    threshC,
                    num_cores,
                )
                zipped = zip(columns, values)
                a_dictionary = dict(zipped)
                data.append(a_dictionary)

    # allow for barplot for single file and region
    if (len(fileNames) == 1) and (len(bedFiles) == 1):
        values = get_counts(
            fileNames[0],
            sampleNames[0],
            bedFiles[0],
            basemod,
            outDir,
            threshA,
            threshC,
            num_cores,
        )
        zipped = zip(columns, values)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)

    df = pd.DataFrame(data)
    # draw from aggregate to calculate modified/total in region of interest
    if len(fileNames) == 1:
        title = "sample_" + fileNames[0].split("/")[-1].replace(".bam", "")
    if len(bedFiles) == 1:
        title = "region_" + bedFiles[0].split("/")[-1].replace(".bed", "")
    if (len(fileNames) == 1) and (len(bedFiles) == 1):
        title = (
            "sample_"
            + fileNames[0].split("/")[-1].replace(".bam", "")
            + "_region_"
            + bedFiles[0].split("/")[-1].replace(".bed", "")
        )
    plot_barchart(df, basemod, outDir, colors, title)

    db_paths = []
    for f in fileNames:
        db = outDir + "/" + f.split("/")[-1].replace(".bam", "") + ".db"
        db_paths.append(db)

    plot_path = f"{outDir}/{title}_{basemod}_enrichment_barplot.png"
    str_out = f"Outputs\n_______\nDB file: {db_paths}\nenrichment barplot: {plot_path}"
    print(str_out)


def get_counts(
    fileName,
    sampleName,
    bedFile,
    basemod,
    outDir,
    threshA,
    threshC,
    num_cores,
):
    # parse_bam for files / regions
    parse_bam(
        fileName,
        sampleName,
        outDir,
        bedFile,
        basemod,
        threshA=threshA,
        threshC=threshC,
        cores=num_cores,
    )

    # get aggregate counts
    aggregate_counts = pd.read_sql(
        "SELECT * from methylationAggregate_" + sampleName,
        sqlite3.connect(
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        ),
    )
    methylated_bases = aggregate_counts["methylated_bases"].sum()
    total_bases = aggregate_counts["total_bases"].sum()
    if total_bases == 0:
        fractionMethylated = 0
    else:
        fractionMethylated = methylated_bases / total_bases

    return [fileName, bedFile, sampleName, fractionMethylated]


def plot_barchart(data, basemod, outDir, colors, title):
    """
    x-axis: sample or region
    y-axis: fraction methylated bases
    """
    fig, ax1 = plt.subplots()
    plt.bar("sampleName", "fractionMethylated", data=data, color=colors)
    print("\nData for barplot")
    print("________________\n")
    print(f"{data.sampleName}")
    print(f"{data.fractionMethylated}")
    print("\n")
    sns.despine(fig)
    plt.ylabel("fraction methylated bases")
    plt.xlabel("")
    plt.savefig(
        outDir + "/" + title + "_" + basemod + "_enrichment_barplot.png",
        dpi=300,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot DiMeLo methylation enrichment"
    )

    # Required arguments
    required_args = parser.add_argument_group("required arguments")
    required_args.add_argument(
        "-f", "--fileNames", required=True,
        nargs="+",
        help="bam file name(s)"
    )
    required_args.add_argument(
        "-s", "--sampleNames", required=True,
        nargs="+",
        help="sample name(s) for output file labelling"
    )
    required_args.add_argument(
        "-b", "--bedFiles", required=True,
        nargs="+",
        help="name of bed file(s) defining region(s) of interest"
    )
    required_args.add_argument(
        "-m", "--basemod", required=True,
        type=str, choices=["A", "CG"],
        help="which base modification to extract"
    )
    required_args.add_argument(
        "-o", "--outDir", required=True,
        help="directory to output plot"
    )

    # Plotting arguments
    plotting_args = parser.add_argument_group("plotting options")
    plotting_args.add_argument(
        "--colors", type=str, nargs="+",
        default=DEFAULT_COLOR_LIST,
        help="color list in hex (e.g. \"#BB4430\") for overlay plots"
    )

    # Optional arguments
    parser.add_argument(
        "-A", "--threshA", type=int,
        default=DEFAULT_THRESH_A,
        help="threshold above which to call an A base methylated"
    )
    parser.add_argument(
        "-C", "--threshC", type=int,
        default=DEFAULT_THRESH_C,
        help="threshold above which to call a C base methylated"
    )
    parser.add_argument(
        "-p", "--cores", type=int,
        help="number of cores over which to parallelize"
    )

    args = parser.parse_args()
    plot_enrichment(**vars(args))
