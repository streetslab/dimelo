r"""
==================================
plot_enrichment_profile module
==================================
.. currentmodule:: dimelo.plot_enrichment_profile
.. autosummary::
    plot_enrichment_profile

plot_enrichment_profile plots single molecules centered at regions of interest defined in bed file and produces aggregate profile

"""

import multiprocessing
import os
import sqlite3
from itertools import cycle
import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from dimelo.parse_bam import parse_bam

COLOR_A = "#053C5E"
COLOR_C = "#BB4430"
COLOR_LIST = ["#2D1E2F", "#A9E5BB", "#610345", "#559CAD", "#5E747F"]


def plot_enrichment_profile(
    fileNames,
    sampleNames,
    bedFiles,
    basemod,
    outDir,
    threshA=129,
    threshC=129,
    windowSize=1000,
    colorA=COLOR_A,
    colorC=COLOR_C,
    colors=COLOR_LIST,
    dotsize=0.5,
    smooth=50,
    min_periods=10,
    cores=None,
):
    """
    fileNames
        name(s) of bam file with Mm and Ml tags
    sampleNames
        name(s) of sample for output file name labelling
    bedFiles
        specified windows for region(s) of interest; optional 4th column in bed file to specify strand of region of interest as ``+`` or ``-``. Default is to consider regions as all ``+``. Reads will be oriented with respect to strand.
    basemod
        One of the following:

        * ``'A'`` - extract mA only
        * ``'CG'`` - extract mCpG only
        * ``'A+CG'`` - extract mA and mCpG
    outDir
        directory to output plot
    threshA
        threshold for calling mA; default 129
    threshC
        threshold for calling mCG; default 129
    windowSize
        window size around center point of feature of interest to plot (+/-); default 1000 bp
    colorA
        color in hex for mA; default #053C5E
    colorC
        color in hex for mCG; default #BB4430
    colors
        color list in hex for overlay plots; default is ["#2D1E2F", "#A9E5BB", "#610345", "#559CAD", "#5E747F"]
    dotsize
        size of points; default is 0.5
    smooth
        window over which to smooth aggregate curve; default of 50 bp
    min_periods
        minimum number of bases to consider for smoothing: default of 10 bp
    cores
        number of cores over which to parallelize; default is all available

    **Example**

    For single file and region:

    >>> dm.plot_enrichment_profile("dimelo/test/data/mod_mappings_subset.bam", "test", "dimelo/test/data/test.bed", "A+CG", "dimelo/dimelo_test", windowSize=500, dotsize=1)

    To overlay multiple regions of interest (can conversely also overlay multiple samples over a single region if a list of files is provided):

    >>> dm.plot_enrichment_profile("dimelo/test/data/mod_mappings_subset.bam", ["test1","test2"], ["dimelo/test/data/test.bed","dimelo/test/data/test.bed"], "A", "dimelo/dimelo_test", windowSize=500, dotsize=1)

    **Return**

        * Aggregate profile of fraction of bases modified centered at features of interest
        * Single molecules centered at features of interest
        * Base abundance centered at features of interest


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

    # if  single bam file rather than list is entered, convert to list
    if type(fileNames) != list:
        fileNames = [fileNames]
    # if single sample name rather than list is entered, convert to list
    if type(sampleNames) != list:
        sampleNames = [sampleNames]
    # if single bed file rather than list is entered, convert to list
    if type(bedFiles) != list:
        bedFiles = [bedFiles]

    db_paths = []
    for f in fileNames:
        db = outDir + "/" + f.split("/")[-1].replace(".bam", "") + ".db"
        db_paths.append(db)

    # overlay condition
    if len(fileNames) > 1 or len(bedFiles) > 1:
        if basemod == "A+CG":
            print(
                "enrichment overlays can only be produced for a single base modification at a time"
            )
            return
        fig = plt.figure()
        if len(fileNames) > 1:
            if len(bedFiles) > 1:
                print(
                    "only a single region file can be used when overlaying multiple bam files"
                )
                return
            for f, n, c in zip(fileNames, sampleNames, cycle(colors)):
                execute_overlay(
                    f,
                    n,
                    c,
                    bedFiles[0],
                    basemod,
                    outDir,
                    threshA,
                    threshC,
                    windowSize,
                    dotsize,
                    smooth,
                    min_periods,
                    num_cores,
                )
        if len(bedFiles) > 1:
            if len(fileNames) > 1:
                print(
                    "only a single bam file can be used when overlaying multiple bed file regions"
                )
                return
            for b, n, c in zip(bedFiles, sampleNames, cycle(colors)):
                execute_overlay(
                    fileNames[0],
                    n,
                    c,
                    b,
                    basemod,
                    outDir,
                    threshA,
                    threshC,
                    windowSize,
                    dotsize,
                    smooth,
                    min_periods,
                    num_cores,
                )
        if len(fileNames) == 1:
            title = "sample_" + fileNames[0].split("/")[-1].replace(".bam", "")
        if len(bedFiles) == 1:
            title = "region_" + bedFiles[0].split("/")[-1].replace(".bed", "")
        plt.title(basemod)
        plt.legend(sampleNames)
        plt.show()
        fig.savefig(
            outDir
            + "/"
            + title
            + "_"
            + basemod
            + "_sm_rolling_avg_overlay.pdf"
        )

        overlay_path = f"{outDir}/{title}_{basemod}_sm_rolling_avg_overlay.pdf"
        str_out = f"Outputs\n_______\nDB file: {db_paths}\noverlay plot: {overlay_path}"
        print(str_out)

    # no overlay condition
    if (len(fileNames) == 1) and (len(bedFiles) == 1):
        execute_single_plot(
            fileNames[0],
            sampleNames[0],
            bedFiles[0],
            basemod,
            outDir,
            threshA,
            threshC,
            windowSize,
            colorA,
            colorC,
            dotsize,
            smooth,
            min_periods,
            num_cores,
        )

        t_paths = []
        if "A" in basemod:
            t_path = (
                outDir + "/" + sampleNames[0] + "_" + "A" + "_base_count.png"
            )
            t_paths.append(t_path)
        if "C" in basemod:
            t_path = (
                outDir + "/" + sampleNames[0] + "_" + "C" + "_base_count.png"
            )
            t_paths.append(t_path)

        enrichment_path = (
            f"{outDir}/{sampleNames[0]}_{basemod}_sm_rolling_avg.pdf"
        )
        sm_path = f"{outDir}/{sampleNames[0]}_{basemod}_sm_scatter.png"
        str_out = f"Outputs\n_______\nDB file: {db_paths}\nenrichment plot: {enrichment_path}\nsingle molecule plot: {sm_path}\nbase count plots: {t_paths}"
        print(str_out)


def execute_overlay(
    fileName,
    sampleName,
    color,
    bedFile,
    basemod,
    outDir,
    threshA,
    threshC,
    windowSize,
    dotsize,
    smooth,
    min_periods,
    num_cores,
):
    parse_bam(
        fileName,
        sampleName,
        outDir,
        bedFile,
        basemod,
        center=True,
        windowSize=windowSize,
        threshA=threshA,
        threshC=threshC,
        cores=num_cores,
    )
    aggregate_counts = pd.read_sql(
        "SELECT * from methylationAggregate_" + sampleName,
        sqlite3.connect(
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        ),
    )
    aggregate_counts["frac"] = (
        aggregate_counts["methylated_bases"] / aggregate_counts["total_bases"]
    )
    if "A" in basemod:
        plot_aggregate_helper(
            aggregate_counts, "A", smooth, min_periods, color
        )

    if "C" in basemod:
        plot_aggregate_helper(
            aggregate_counts, "C", smooth, min_periods, color
        )


def execute_single_plot(
    fileName,
    sampleName,
    bedFile,
    basemod,
    outDir,
    threshA,
    threshC,
    windowSize,
    colorA,
    colorC,
    dotsize,
    smooth,
    min_periods,
    num_cores,
):

    parse_bam(
        fileName,
        sampleName,
        outDir,
        bedFile,
        basemod,
        center=True,
        windowSize=windowSize,
        threshA=threshA,
        threshC=threshC,
        cores=num_cores,
    )

    all_data = pd.read_sql(
        "SELECT * from methylationByBase_" + sampleName,
        sqlite3.connect(
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        ),
    )
    aggregate_counts = pd.read_sql(
        "SELECT * from methylationAggregate_" + sampleName,
        sqlite3.connect(
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        ),
    )

    print(
        "processing "
        + str(len(all_data["read_name"].unique()))
        + " reads with methylation above threshold for "
        + sampleName
        + " for bam: "
        + fileName
    )

    fig, ax = plt.subplots()

    colors = {"A+Y": colorA, "A+a": colorA, "C+Z": colorC, "C+m": colorC}

    sns.scatterplot(
        data=all_data,
        x="pos",
        y="read_name",
        hue="mod",
        palette=colors,
        s=dotsize,
        marker="s",
        linewidth=0,
        legend=None,
    )

    ax.spines[["top", "right", "left"]].set_visible(False)

    plt.yticks([])
    plt.ylabel("")
    plt.xlabel("")
    plt.xlim(-windowSize, windowSize)
    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_sm_scatter.png", dpi=600
    )

    plot_aggregate_me_frac(
        sampleName,
        aggregate_counts,
        smooth,
        min_periods,
        windowSize,
        basemod,
        outDir,
        colorA,
        colorC,
    )


# average profile
def plot_aggregate_me_frac(
    sampleName,
    aggregate_counts,
    smooth,
    min_periods,
    windowSize,
    basemod,
    outDir,
    colorA,
    colorC,
):
    aggregate_counts["frac"] = (
        aggregate_counts["methylated_bases"] / aggregate_counts["total_bases"]
    )

    fig = plt.figure()
    labels = []
    if "A" in basemod:
        plot_aggregate_helper(
            aggregate_counts, "A", smooth, min_periods, colorA
        )
        labels.append("A")
    if "C" in basemod:
        plot_aggregate_helper(
            aggregate_counts, "C", smooth, min_periods, colorC
        )
        labels.append("CG")
    plt.title(basemod)
    plt.legend(labels)
    plt.show()
    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_sm_rolling_avg.pdf"
    )

    if "A" in basemod:
        aggregate_A = aggregate_counts[
            aggregate_counts["mod"].str.contains("A")
        ]
        # need to sort first!
        aggregate_A.sort_values(["pos"], inplace=True)
        plot_base_abundance(
            sampleName,
            aggregate_A,
            "A",
            windowSize,
            outDir,
        )
    if "C" in basemod:
        aggregate_C = aggregate_counts[
            aggregate_counts["mod"].str.contains("C")
        ]
        # need to sort first!
        aggregate_C.sort_values(["pos"], inplace=True)
        plot_base_abundance(
            sampleName,
            aggregate_C,
            "CG",
            windowSize,
            outDir,
        )


# helper function to create smoothed lineplot
def plot_aggregate_helper(aggregate_counts, mod, smooth, min_periods, color):
    aggregate = aggregate_counts[aggregate_counts["mod"].str.contains(mod)]
    # need to sort first!
    aggregate.sort_values(["pos"], inplace=True)
    aggregate_rolling = aggregate.rolling(
        window=smooth, min_periods=min_periods, center=True, on="pos"
    ).mean()
    sns.lineplot(
        x=aggregate_rolling["pos"], y=aggregate_rolling["frac"], color=color
    )


def plot_base_abundance(
    sampleName, aggregate_counts, basemod, windowSize, outDir
):
    cmapPurple = colors.LinearSegmentedColormap.from_list(
        "custom purple", ["white", "#2D1E2F"], N=200
    )
    # plot base abundance
    fig = plt.figure()
    x = aggregate_counts["pos"].to_numpy()
    y = aggregate_counts["total_bases"].to_numpy()  # base_count
    fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
    extent = [x[0] - (x[1] - x[0]) / 2.0, x[-1] + (x[1] - x[0]) / 2.0, 0, 1]
    # extent = [x[0], x[-1], 0, 1]
    im = ax.imshow(
        y[np.newaxis, :], cmap=cmapPurple, aspect="auto", extent=extent
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.25)
    fig.colorbar(im, cax=cax, orientation="horizontal")
    ax.set_yticks([])
    ax.set_xlim(extent[0], extent[1])
    ax2.plot(x, y, "o", ms=0.5, color="#2D1E2F")
    ax2.set_xlim(extent[0], extent[1])
    # ax.spines["top"].set_visible(False)
    # ax.spines["right"].set_visible(False)
    # ax.spines["bottom"].set_visible(False)
    # ax.spines["left"].set_visible(False)
    # ax.get_xaxis().set_ticks([])
    plt.tight_layout()
    plt.show()
    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_base_count.png", dpi=600
    )


def main():
    parser = argparse.ArgumentParser(
        description="Plot DiMeLo enrichment profile"
    )

    # Required arguments
    parser.add_argument(
        "-f", "--fileNames", required=True,
        nargs="+",
        help="bam file name(s)"
    )
    parser.add_argument(
        "-s", "--sampleNames", required=True,
        nargs="+",
        help="sample name(s) for output file labelling"
    )
    parser.add_argument(
        "-b", "--bedFiles", required=True,
        nargs="+",
        help="name of bed file(s) defining region(s) of interest"
    )
    parser.add_argument(
        "-m", "--basemod", required=True,
        type=str, choices=["A", "CG", "A+CG"],
        help="which base modification to extract"
    )
    parser.add_argument(
        "-o", "--outDir", required=True,
        help="directory to output plot"
    )

    # Optional arguments
    parser.add_argument(
        "-A", "--threshA", type=int,
        default=129,
        help="threshold above which to call an A base methylated"
    )
    parser.add_argument(
        "-C", "--threshC", type=int,
        default=129,
        help="threshold above which to call a C base methylated"
    )
    parser.add_argument(
        "-w", "--windowSize", type=int,
        default=1000,
        help="window size around center point of feature of interest to plot (+/-)"
    )
    parser.add_argument(
        "--colorA", type=str,
        default=COLOR_A,
        help="color in hex (e.g. \"#BB4430\") for mA"
    )
    parser.add_argument(
        "--colorC", type=str,
        default=COLOR_C,
        help="color in hex (e.g. \"#BB4430\") for mCG"
    )
    parser.add_argument(
        "-c", "--colors", type=str, nargs="+",
        default=COLOR_LIST,
        help="color list in hex (e.g. \"#BB4430\") for overlay plots"
    )
    parser.add_argument(
        "-d", "--dotsize", type=float,
        default=0.5,
        help="size of points"
    )
    parser.add_argument(
        "-t", "--smooth", type=int,
        default=50,
        help="window over which to smooth aggregate curve"
    )
    parser.add_argument(
        "-n", "--min_periods", type=int,
        default=10,
        help="minimum number of bases to consider for smoothing"
    )
    parser.add_argument(
        "-p", "--cores", type=int,
        help="number of cores over which to parallelize"
    )

    args = parser.parse_args()
    plot_enrichment_profile(**vars(args))
