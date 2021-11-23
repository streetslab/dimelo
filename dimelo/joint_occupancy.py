import multiprocessing
import sqlite3

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from joblib import Parallel, delayed
from matplotlib import colors
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


def joint_occupancy(
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
):
    peak_left, peak_right = extract_peak_pairs(
        bedFile, min_distance, max_distance, outDir
    )

    peak_left_windows = make_windows(peak_left)
    peak_right_windows = make_windows(peak_right)

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
    )

    all_data = pd.read_sql(
        "SELECT * from methylationByBaseJoint_" + sampleName,
        sqlite3.connect(
            outDir + "/" + fileName.split("/")[-1].split(".")[0] + ".db"
        ),
    )

    if "A" in basemod:
        all_data_A_binned = bin_probabilities(all_data, "A")  # all_data
        print(
            "processing "
            + str(len(all_data_A_binned["read_name"].unique()))
            + " reads for for bam: "
            + fileName
        )
        print(all_data_A_binned.shape[0])
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
            )
    if "C" in basemod:
        all_data_C_binned = bin_probabilities(all_data, "C")
        print(
            "processing "
            + str(len(all_data_A_binned["read_name"].unique()))
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
):

    make_db(fileName, sampleName, outDir, False, False, True)

    num_cores = multiprocessing.cpu_count()
    batchSize = 100

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
            batchSize,
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
    batchSize,
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
            :param batchSize: number of reads to submit to db at once
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
    print(w1.string)
    print(w2.string)
    print(count)
    if data:
        DATABASE_NAME = (
            outDir + "/" + fileName.split("/")[-1].split(".")[0] + ".db"
        )
        table_name = "methylationByBaseJoint_" + sampleName
        command = (
            """INSERT OR IGNORE INTO """
            + table_name
            + """ VALUES(?,?,?,?,?,?,?);"""
        )
        execute_sql_command(command, DATABASE_NAME, data)


def bin_probabilities(all_data, mod):
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
            lambda b: prob_bin(b)
        )
        all_data_mod.loc[
            all_data_mod["read_windows"] == r, "prob"
        ] = binned_probs  # no longer id because id has position
    return all_data_mod


def prob_bin(bin):
    """
    probability a base in the window is methylated by:
    calculating probability that no base in the window is methylated and then taking the complement
    treat p=1 as 254/255 for prevent log(0)
    """
    probs = [
        np.log(1 - p) for p in bin if ((p < 1) and (p > 0.5))
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
        outDir + "/" + sampleName + "_" + basemod + "_joint_occupancy.png",
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
    k = KMeans(
        n_clusters=3, random_state=1  # try with 3 clusters
    )  # try different numbers of clusters #2
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
        by=["labels", "left_sum"], axis=0, ascending=False
    )
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


def main():
    # TODO add argument parsing for command line
    print("main")
