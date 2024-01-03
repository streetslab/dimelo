r"""
==================================
plot_enrichment_profile module
==================================
.. currentmodule:: dimelo.plot_enrichment_profile
.. autosummary::
    plot_enrichment_profile

plot_enrichment_profile plots single molecules centered at regions of interest defined in bed file and produces aggregate profile

"""

import argparse
import multiprocessing
import os
import sqlite3
from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from dimelo.parse_bam import parse_bam

DEFAULT_THRESH_A = 129
DEFAULT_THRESH_C = 129
DEFAULT_WINDOW_SIZE = 1000
COLOR_A = "#053C5E"
COLOR_C = "#BB4430"
COLOR_LIST = ["#2D1E2F", "#A9E5BB", "#610345", "#559CAD", "#5E747F"]
DEFAULT_DOTSIZE = 0.5
DEFAULT_SMOOTH = 50
DEFAULT_MIN_PERIODS = 10


def plot_enrichment_profile(
    fileNames,
    sampleNames,
    bedFiles,
    basemod,
    outDir,
    threshA=DEFAULT_THRESH_A,
    threshC=DEFAULT_THRESH_C,
    windowSize=DEFAULT_WINDOW_SIZE,
    colorA=COLOR_A,
    colorC=COLOR_C,
    colors=COLOR_LIST,
    dotsize=DEFAULT_DOTSIZE,
    smooth=DEFAULT_SMOOTH,
    min_periods=DEFAULT_MIN_PERIODS,
    cores=None,
):
    """
    fileNames
        name(s) of bam file with Mm and Ml tags
    sampleNames
        name(s) of sample for output file name labelling; valid names contain [``a-zA-Z0-9_``].
    bedFiles
        specified windows for region(s) of interest; optional 4th column in bed file to specify strand of region of interest as ``+`` or ``-``. Default is to consider regions as all ``+``. Reads will be oriented with respect to strand. Only reads overlapping regions defined in bed file will be extracted, regardless of windowSize. Plots are centered at the center of the bed file regions.
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

    **Example Plots**

        * :ref:`sphx_glr_auto_examples_enrichment_profile_single_example.py`
        * :ref:`sphx_glr_auto_examples_enrichment_profile_ma_mc_example.py`
        * :ref:`sphx_glr_auto_examples_enrichment_profile_overlay_example.py`

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
            raise RuntimeError(
                "enrichment overlays can only be produced for a single base modification at a time"
            )
        fig = plt.figure()
        if len(fileNames) > 1:
            if len(bedFiles) > 1:
                raise RuntimeError(
                    "only a single region file can be used when overlaying multiple bam files"
                )
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
                raise RuntimeError(
                    "only a single bam file can be used when overlaying multiple bed file regions"
                )
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
        fig.savefig(
            outDir
            + "/"
            + title
            + "_"
            + basemod
            + "_sm_rolling_avg_overlay.pdf"
        )
        plt.close()

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
                outDir + "/" + sampleNames[0] + "_" + "CG" + "_base_count.png"
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
    plt.close()

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
    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_sm_rolling_avg.pdf"
    )
    plt.close()

    if "A" in basemod:
        aggregate_A = aggregate_counts[
            aggregate_counts["mod"].str.contains("A")
        ].copy()
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
        ].copy()
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
    aggregate = aggregate_counts[
        aggregate_counts["mod"].str.contains(mod)
    ].copy()
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
    aggregate_counts = (
        aggregate_counts.set_index("pos")
        .reindex(
            pd.Index(
                np.arange(
                    aggregate_counts["pos"].min(),
                    aggregate_counts["pos"].max(),
                    1,
                ),
                name="pos",
            )
        )
        .reset_index()
    )
    aggregate_counts = aggregate_counts.fillna(0)
    fig = plt.figure()
    x = aggregate_counts["pos"].to_numpy()
    y = aggregate_counts["total_bases"].to_numpy()  # base_count
    fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True)
    extent = [x[0] - (x[1] - x[0]) / 2.0, x[-1] + (x[1] - x[0]) / 2.0, 0, 1]
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
    plt.tight_layout()
    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_base_count.png", dpi=600
    )
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot DiMeLo enrichment profile"
    )

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
        help="sample name(s) for output file labelling",
    )
    required_args.add_argument(
        "-b",
        "--bedFiles",
        required=True,
        nargs="+",
        help="name of bed file(s) defining region(s) of interest",
    )
    required_args.add_argument(
        "-m",
        "--basemod",
        required=True,
        type=str,
        choices=["A", "CG", "A+CG"],
        help="which base modification to extract",
    )
    required_args.add_argument(
        "-o", "--outDir", required=True, help="directory to output plot"
    )

    # Smoothing options
    smoothing_args = parser.add_argument_group("smoothing options")
    smoothing_args.add_argument(
        "-t",
        "--smooth",
        type=int,
        default=DEFAULT_SMOOTH,
        help="window over which to smooth aggregate curve",
    )
    smoothing_args.add_argument(
        "-n",
        "--min_periods",
        type=int,
        default=DEFAULT_MIN_PERIODS,
        help="minimum number of bases to consider for smoothing",
    )

    # Plotting arguments
    plotting_args = parser.add_argument_group("plotting options")
    plotting_args.add_argument(
        "--colorA",
        type=str,
        default=COLOR_A,
        help='color in hex (e.g. "#BB4430") for mA',
    )
    plotting_args.add_argument(
        "--colorC",
        type=str,
        default=COLOR_C,
        help='color in hex (e.g. "#BB4430") for mCG',
    )
    plotting_args.add_argument(
        "--colors",
        type=str,
        nargs="+",
        default=COLOR_LIST,
        help='color list in hex (e.g. "#BB4430") for overlay plots',
    )
    plotting_args.add_argument(
        "-d",
        "--dotsize",
        type=float,
        default=DEFAULT_DOTSIZE,
        help="size of points",
    )

    # Optional arguments
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
        "-w",
        "--windowSize",
        type=int,
        default=DEFAULT_WINDOW_SIZE,
        help="window size around center point of feature of interest to plot (+/-)",
    )
    parser.add_argument(
        "-p",
        "--cores",
        type=int,
        help="number of cores over which to parallelize",
    )

    args = parser.parse_args()
    plot_enrichment_profile(**vars(args))
