import multiprocessing
import sqlite3

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
    """Create single molecule plots centered at region of interest.
    Overlay EITHER multiple files OR multiple bed regions for a single file
    Args:
            :param fileNames: name(s) of bam file with Mm and Ml tags
            :param sampleNames: name(s) of sample for output file name labelling
            :param bedFiles: specified windows for region(s) of interest
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param outDir: directory to output plot
            :param threshA: threshold for calling mA; default 129
            :param threshC: threshold for calling mCG; default 129
            :param windowSize: window size around center point of feature of interest to plot (+/-); default 1000 bp
            :param colorA: color in hex for mA
            :param colorC: color in hex for mCG
            :param colors: color list in hex for overlay; default is ["#2D1E2F", "#A9E5BB", "#610345", "#559CAD", "#5E747F"]
            :param dotsize: size of points
            :param smooth: window over which to smooth aggregate curve; default of 50 bp
            :param min_periods: minimum number of bases to consider for smoothing: default of 10 bp
    Return:
            plot of single molecules centered at region of interest
    """

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
            for f, n, c in zip(fileNames, sampleNames, colors):
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
            for b, n, c in zip(bedFiles, sampleNames, colors):
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
        plt.title(basemod)
        plt.legend(sampleNames)
        plt.show()
        fig.savefig(outDir + "/" + basemod + "_sm_rolling_avg_overlay.pdf")

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
            outDir + "/" + fileName.split("/")[-1].split(".")[0] + ".db"
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
            outDir + "/" + fileName.split("/")[-1].split(".")[0] + ".db"
        ),
    )
    aggregate_counts = pd.read_sql(
        "SELECT * from methylationAggregate_" + sampleName,
        sqlite3.connect(
            outDir + "/" + fileName.split("/")[-1].split(".")[0] + ".db"
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
        plot_base_abundance(
            sampleName,
            aggregate_counts[
                aggregate_counts["mod"].str.contains("A")
            ].sort_values(["pos"], inplace=True),
            "A",
            windowSize,
            outDir,
        )
    if "C" in basemod:
        plot_base_abundance(
            sampleName,
            aggregate_counts[
                aggregate_counts["mod"].str.contains("C")
            ].sort_values(["pos"], inplace=True),
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
    im = ax.imshow(
        y[np.newaxis, :], cmap=cmapPurple, aspect="auto", extent=extent
    )
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.25)
    fig.colorbar(im, cax=cax, orientation="horizontal")
    ax.set_yticks([])
    ax.set_xlim(extent[0], extent[1])
    ax2.plot(x, y, "o", ms=0.5, color="#2D1E2F")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.get_xaxis().set_ticks([])
    plt.tight_layout()
    plt.show()

    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_base_count.png", dpi=600
    )


def main():
    # TODO add argument parsing
    print("main")
