r"""
=================
plot_joint_enrichment module
=================
.. currentmodule:: dimelo.plot_joint_enrichment
.. autosummary::
    plot_joint_enrichment

plot_joint_enrichment plots single molecules that span two sites of interest

"""

import multiprocessing
import os
import sqlite3

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from joblib import Parallel, delayed
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pybedtools import BedTool
from sklearn.cluster import KMeans

from dimelo.parse_bam import get_modified_reference_positions, make_db
from dimelo.utils import execute_sql_command

COLOR_A = "#053C5E"
COLOR_C = "#BB4430"


class Region(object):
    def __init__(self, region):
        self.chromosome = region[1][0]
        # now being and end correspond to motif begin and end
        self.begin = region[1][1]
        self.end = region[1][2]
        self.size = self.end - self.begin
        self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        # additional info to store about CTCF regions
        self.strand = region[1][3]
        self.peak_strength = region[1][4]


def plot_joint_enrichment(
    fileName,
    sampleName,
    bedFile,
    basemod,
    outDir,
    threshA=129,
    threshC=129,
    windowSize=1000,
    colorA=COLOR_A,
    colorC=COLOR_C,
    dotsize=0.5,
    min_distance=2000,
    max_distance=10000,
    gap=1000,
    cores=None,
    num_clusters=3,
):
    """
    fileName
        name of bam file with Mm and Ml tags
    sampleName
        name of sample for output file name labelling
    bedFile
        specified windows for region(s) of interest
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
    dotsize
        size of points; default is 0.5
    min_distance
        minimum distance between two regions of interest spanned by a single molecule; default is 2000 bp
    max_distance
        maximum distance between two regions of interest spanned by a single molecule; default is 10000 bp
    gap
        gap between neighboring viewpoints; default is 1000 bp
    cores
        number of cores over which to parallelize; default is all available
    num_clusters
        number of clusters for kmeans clustering for visualization; default is 3

    **Example**

    >>> dm.plot_joint_enrichment("dimelo/test/data/mod_mappings_subset.bam", "joint_test", "dimelo/test/data/test.bed", "A",  "dimelo/dimelo_test")

    **Return**

    Single molecules that span two sites of interest with base modifications colored, clustered for visualization with kmeans clustering

    """

    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    peak_left, peak_right = extract_peak_pairs(
        bedFile, min_distance, max_distance, outDir
    )

    peak_left_windows = make_windows(peak_left)
    peak_right_windows = make_windows(peak_right)

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

    parse_bam_paired(
        fileName,
        sampleName,
        outDir,
        peak_left_windows,
        peak_right_windows,
        basemod,
        True,
        windowSize,
        0,  # get all As (for binning); threshold before plotting
        0,  # get all As (for binning); threshold before plotting
        gap,
        num_cores,
    )

    all_data = pd.read_sql(
        "SELECT * from methylationByBaseJoint_" + sampleName,
        sqlite3.connect(
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        ),
    )

    if "A" in basemod:
        all_data_A_binned = bin_probabilities(
            all_data, "A", threshA
        )  # all_data
        print(
            "processing "
            + str(len(all_data_A_binned["read_name"].unique()))
            + " reads for for bam: "
            + fileName
        )
        if all_data_A_binned.shape[0] > 0:
            make_cluster_plot(
                all_data_A_binned,
                sampleName,
                basemod,
                threshA,
                outDir,
                windowSize,
                gap,
                colorA,
                dotsize,
                num_clusters,
            )
    if "C" in basemod:
        all_data_C_binned = bin_probabilities(all_data, "C", threshC)
        print(
            "processing "
            + str(len(all_data_C_binned["read_name"].unique()))
            + " reads for for bam: "
            + fileName
        )
        if all_data_C_binned.shape[0] > 0:
            make_cluster_plot(
                all_data_C_binned,
                sampleName,
                basemod,
                threshC,
                outDir,
                windowSize,
                gap,
                colorC,
                dotsize,
                num_clusters,
            )


def extract_peak_pairs(bed, min_dist, max_dist, outDir):
    b = BedTool(bed)
    b.sort().saveas(outDir + "/tmp.sorted.bed")
    bed = pd.read_csv(outDir + "/tmp.sorted.bed", sep="\t", header=None)
    # find middle
    left_keep = []
    right_keep = []
    bed["middle"] = (bed[1] + bed[2]) / 2
    for i in range(0, len(bed) - 1):
        if (bed.iloc[i + 1].middle - bed.iloc[i].middle >= min_dist) & (
            bed.iloc[i + 1].middle - bed.iloc[i].middle <= max_dist
        ):
            left_keep.append(i)
            right_keep.append(i + 1)
    print("number of peak pairs: " + str(bed.iloc[left_keep].shape[0]))
    return bed.iloc[left_keep], bed.iloc[right_keep]


def make_windows(bed):
    reg = []
    for row in bed.iterrows():
        reg.append(Region(row))
    return reg


def parse_bam_paired(
    fileName,
    sampleName,
    outDir,
    windows1,
    windows2,
    basemod,
    center,
    windowSize,
    threshA,
    threshC,
    gap,
    num_cores,
):

    make_db(fileName, sampleName, outDir, False, False, True)

    Parallel(n_jobs=num_cores, verbose=10)(
        delayed(parse_reads_windows)(
            fileName,
            sampleName,
            basemod,
            windowSize,
            window1,
            window2,
            center,
            threshA,
            threshC,
            outDir,
            True,  # want to extract all bases for binning
            False,
            gap,
        )
        for window1, window2 in zip(windows1, windows2)
    )


def parse_reads_windows(
    fileName,
    sampleName,
    basemod,
    windowSize,
    w1,
    w2,
    center,
    threshA,
    threshC,
    outDir,
    extractAllBases,
    qc,
    gap,
):
    """Parse all reads in batchs
    Args:
            :param fileName: name of bam file
            :param sampleName: name of sample for output file name labelling
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param windowSize: window size around center point of feature of interest to plot (+/-); only mods within this window are stored; only applicable for center=True
            :param windows1: windows from peak_left bed file
            :param windows2: windows from peak_right bed file
            :param center: report positions with respect to reference center (+/- window size) if True or in original reference space if False
            :param threshA: threshold above which to call an A base methylated
            :param threshC: threshold above which to call a C base methylated
            :param gap: gap between neighboring viewpoints
    Return:
            data to put into methylationByBaseJoint table
    """
    data = []

    bam = pysam.AlignmentFile(fileName, "rb")

    reads1 = bam.fetch(reference=w1.chromosome, start=w1.begin, end=w1.end)
    reads2 = bam.fetch(reference=w2.chromosome, start=w2.begin, end=w2.end)
    count = 0
    for read in reads1:
        # only add reads that span both sites
        # includeunique identifier that is w1-w2-read_name
        if read in reads2:
            count = count + 1
            [
                (mod, positionsL, probs),
                (mod2, positions2L, probs2),
            ] = get_modified_reference_positions(
                read,
                basemod,
                w1,
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
            [
                (mod, positionsR, probs),
                (mod2, positions2R, probs2),
            ] = get_modified_reference_positions(
                read,
                basemod,
                w2,
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
            for posL, posR, prob in zip(positionsL, positionsR, probs):
                if posL is not None:
                    if abs(posL) <= windowSize:
                        data.append(
                            (
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string
                                + ":"
                                + str(posL),
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string,
                                read.query_name,
                                int(posL),
                                int(prob),
                                mod,
                                w1.peak_strength,
                            )
                        )
                    if abs(posR) <= windowSize:
                        posR_adj = posR + 2 * windowSize + gap
                        data.append(
                            (
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string
                                + ":"
                                + str(posR_adj),
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string,
                                read.query_name,
                                int(posR_adj),
                                int(prob),
                                mod,
                                w2.peak_strength,
                            )
                        )
            for posL, posR, prob in zip(positions2L, positions2R, probs2):
                if posL is not None:
                    if abs(posL) <= windowSize:
                        data.append(
                            (
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string
                                + ":"
                                + str(posL),
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string,
                                read.query_name,
                                int(posL),
                                int(prob),
                                mod2,
                                w1.peak_strength,
                            )
                        )
                    if abs(posR) <= windowSize:
                        posR_adj = posR + 2 * windowSize + gap
                        data.append(
                            (
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string
                                + ":"
                                + str(posR_adj),
                                read.query_name
                                + ":"
                                + w1.string
                                + "-"
                                + w2.string,
                                read.query_name,
                                int(posR_adj),
                                int(prob),
                                mod2,
                                w2.peak_strength,
                            )
                        )
    if data:
        DATABASE_NAME = (
            outDir + "/" + fileName.split("/")[-1].replace(".bam", "") + ".db"
        )
        table_name = "methylationByBaseJoint_" + sampleName
        command = (
            """INSERT OR IGNORE INTO """
            + table_name
            + """ VALUES(?,?,?,?,?,?,?);"""
        )
        execute_sql_command(command, DATABASE_NAME, data)


def bin_probabilities(all_data, mod, thresh):
    """
    bin data for 10 A's as slide across the read
    calculate the probability at least one base is methylated
    """
    all_data_mod = all_data[all_data["mod"].str.contains(mod)]
    all_data_mod.loc[:, "prob"] = all_data_mod["prob"] / 255
    read_ids = all_data_mod[
        "read_windows"
    ].unique()  # no longer id because id has position
    for r in read_ids:
        subset = all_data_mod[
            all_data_mod["read_windows"] == r
        ]  # no longer id because id has position
        probs = subset["prob"]
        binned_probs = probs.rolling(window=20, center=True).apply(
            lambda b: prob_bin(b, thresh)
        )
        all_data_mod.loc[
            all_data_mod["read_windows"] == r, "prob"
        ] = binned_probs  # no longer id because id has position
    return all_data_mod


def prob_bin(bin, thresh):
    """
    probability a base in the window is methylated by:
    calculating probability that no base in the window is methylated and then taking the complement
    treat p=1 as 254/255 for prevent log(0)
    """
    probs = [
        np.log(1 - p) for p in bin if ((p < 1) and (p >= (thresh / 255)))
    ]  # only consider probabilities > 0.5 and handle 1 later
    probs1 = [np.log(1 - 254 / 255) for p in bin if p == 1]
    probsAll = probs + probs1
    prob = 1 - np.exp(sum(probsAll))
    return prob


def make_cluster_plot(
    all_data,
    sampleName,
    basemod,
    thresh,
    outDir,
    windowSize,
    gap,
    color,
    dotsize,
    num_clusters,
):
    # all_data is already threshold to only contain
    all_data_t = all_data[
        all_data["prob"] > 0.9
    ]  # TODO:stick with 0.9 for now (thresh later)

    # require that quality > thresh be within 100 bp of the peak center on either side
    peak = all_data_t[abs(all_data_t["pos"]) <= 100]
    peak2 = all_data_t[
        (abs(all_data_t["pos"]) > 2 * windowSize + gap - 100)
        & (abs(all_data_t["pos"]) < 2 * windowSize + gap + 100)
    ]
    peak_ids = peak["read_windows"].unique()
    peak_ids2 = peak2["read_windows"].unique()
    boolean_keep_series = all_data_t.read_windows.isin(
        peak_ids
    ) | all_data_t.read_windows.isin(
        peak_ids2
    )  # reads_keep
    all_data_t_p = all_data_t[boolean_keep_series]

    print(
        "number of reads in final plot: "
        + str(len(all_data_t_p.read_name.unique()))
    )

    fig, ax = plt.subplots()
    sns.scatterplot(
        data=all_data_t_p,
        x="pos",
        y="read_name",
        color=color,
        s=dotsize,
        marker="s",
        linewidth=0,
        legend=None,
    )

    # TODO: add box for gap region
    ax.spines[["top", "right", "left"]].set_visible(False)
    plt.yticks([])
    plt.ylabel("")
    plt.xlabel("")
    plt.xlim(-windowSize, 3 * windowSize + gap)

    # Create a Rectangle patch
    rect = patches.Rectangle(
        (windowSize, 0),
        gap,
        len(all_data_t_p.read_windows.unique()),
        linewidth=1,
        edgecolor="grey",
        facecolor="grey",
    )

    # Add the patch to the Axes
    ax.add_patch(rect)

    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_joint_enrichment.png",
        dpi=600,
    )

    # cluster plots
    # TODO: speed this up
    all_data_pivoted = pd.pivot_table(
        all_data_t_p, values="prob", columns="pos", index="read_windows"
    )  # index was read_name #p
    r = range(-windowSize, 4 * windowSize + gap + 1, 1)
    for bp in r:
        if bp not in all_data_pivoted.columns:
            all_data_pivoted[bp] = np.nan
    all_data_pivoted = all_data_pivoted.sort_index(axis=1)
    all_data_pivoted_0 = all_data_pivoted.fillna(0)
    cmapA = colors.LinearSegmentedColormap.from_list(
        "custom", ["white", color], N=255
    )

    all_data_pivoted_mod_0_rolling = pd.DataFrame()
    for i in range(0, all_data_pivoted_0.shape[0]):
        all_data_pivoted_mod_0_rolling = all_data_pivoted_mod_0_rolling.append(
            all_data_pivoted_0.iloc[i, :]
            .rolling(window=5)
            .mean()  # TODO: add on position?
        )  # 33
    all_data_pivoted_mod_0_rolling_0 = all_data_pivoted_mod_0_rolling.fillna(0)
    k = KMeans(n_clusters=num_clusters, random_state=1)
    k.fit(all_data_pivoted_mod_0_rolling_0)

    # sort by left peak signal strength after labels
    all_data_pivoted_0["left_sum"] = 0
    subset_left_peak = all_data_pivoted_0.iloc[
        :, (windowSize - 100) : (windowSize + 100)
    ]
    for idx, row in subset_left_peak.iterrows():
        left_sum = row.sum()
        all_data_pivoted_0.loc[idx, "left_sum"] = left_sum
    all_data_pivoted_0["labels"] = k.labels_
    all_data_pivoted_0 = all_data_pivoted_0.sort_values(
        by=["labels"], axis=0, ascending=False
    )  # sort by labels only
    # all_data_pivoted_0 = all_data_pivoted_0.sort_values(
    #     by=["labels", "left_sum"], axis=0, ascending=False
    # )
    to_plot = all_data_pivoted_0
    to_plot_2 = to_plot.loc[
        :, (to_plot.columns != "labels") & (to_plot.columns != "left_sum")
    ]

    fig, ax = plt.subplots()
    sns.heatmap(
        to_plot_2,
        cmap=cmapA,
        xticklabels=False,
        yticklabels=False,
        cbar_kws=dict(use_gridspec=False, location="top"),
    )

    # Create a Rectangle patch
    rect = patches.Rectangle(
        (2 * windowSize, 0),
        gap,
        len(all_data_t_p.read_windows.unique()),
        linewidth=1,
        edgecolor="lightgrey",
        facecolor="lightgrey",
    )

    # Add the patch to the Axes
    ax.add_patch(rect)

    plt.show()
    fig.savefig(
        outDir
        + "/"
        + sampleName
        + "_"
        + str(thresh)
        + "_cluster_double_peak.4.rolling.png",
        dpi=500,
    )

    # add ChIP-seq signal strength
    cmapPurple = colors.LinearSegmentedColormap.from_list(
        "custom purple", ["white", "#610345"], N=200
    )

    ordered_read_ids = to_plot.index.values
    print(
        "all clusters: "
        + str(len(all_data_pivoted_mod_0_rolling_0.index.values))
    )
    w1_signal_strength = []
    w2_signal_strength = []
    for r in ordered_read_ids:
        read_data = all_data[all_data["read_windows"] == r]
        w1_peak = read_data[read_data["pos"] <= windowSize][
            "peak_strength"
        ].unique()
        w2_peak = read_data[read_data["pos"] >= (windowSize + gap)][
            "peak_strength"
        ].unique()
        w1_signal_strength.append(w1_peak[0])
        w2_signal_strength.append(w2_peak[0])

    fig = plt.figure()
    x = np.linspace(0, len(ordered_read_ids), num=len(ordered_read_ids))
    y = pd.Series(w1_signal_strength)
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
    fig.savefig(outDir + "/" + sampleName + "_peak_signal_w1.png", dpi=300)

    # window 2
    # from positions (windowSize+gap) to (windowSize+gap+2*windowSize)
    fig = plt.figure()
    x = np.linspace(0, len(ordered_read_ids), num=len(ordered_read_ids))
    y = pd.Series(w2_signal_strength)
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
    fig.savefig(outDir + "/" + sampleName + "_peak_signal_w2.png", dpi=300)


def main():
    # TODO add argument parsing for command line
    print("main")
