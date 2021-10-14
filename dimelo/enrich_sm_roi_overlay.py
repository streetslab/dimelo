import sqlite3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from dimelo.parse_bam import parse_bam


def enrich_sm_roi_overlay(
    fileNames,
    sampleNames,
    bedFile,
    basemod,
    outDir,
    threshA=128,
    threshC=128,
    windowSize=1000,
    colors=["#2D1E2F", "#610345", "#559CAD", "#A9E5BB", "#5E747F"],
    dotsize=0.5,
    smooth=50,
    min_periods=10,
):
    """Create single molecule plots centered at region of interest.
    Args:
            :param fileNames: list names of bam files to overlay with Mm and Ml tags
            :param sampleName: list of names of sample for output file name labelling
            :param bedFile: specified windows for regions of interest
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param outDir: directory to output plot
            :param threshA: threshold for calling mA; default 128
            :param threshC: threshold for calling mCG; default 128
            :param windowSize: window size around center point of feature of interest to plot (+/-); default 1000 bp
            :param colors: list of hex colors for overlay
            :param colorC: color in hex for mCG
            :param dotsize: size of points
            :param smooth: window over which to smooth aggregate curve; default of 50 bp
            :param min_periods: minimum number of bases to consider for smoothing: default of 10 bp
    Return:
            plot of single molecules centered at region of interest
    """
    fig = plt.figure()
    i = 0
    for f in fileNames:
        parse_bam(
            f,
            sampleNames[i],
            bedFile,
            basemod,
            center=True,
            windowSize=windowSize,
            threshA=threshA,
            threshC=threshC,
        )
        aggregate_counts = pd.read_sql(
            "SELECT * from methylationAggregate", sqlite3.connect(f + ".db")
        )
        aggregate_counts["frac"] = (
            aggregate_counts["methylated_bases"]
            / aggregate_counts["total_bases"]
        )
        if "A" in basemod:
            aggregate_A = aggregate_counts[
                aggregate_counts["mod"].str.contains("A")
            ]
            # need to sort first!
            aggregate_A.sort_values(["pos"], inplace=True)
            aggregate_A_rolling = aggregate_A.rolling(
                window=smooth, min_periods=min_periods, center=True, on="pos"
            ).mean()
            sns.lineplot(
                x=aggregate_A_rolling["pos"],
                y=aggregate_A_rolling["frac"],
                color=colors[i],
            )
        if "C" in basemod:
            aggregate_C = aggregate_counts[
                aggregate_counts["mod"].str.contains("C")
            ]
            # need to sort first!
            aggregate_C.sort_values(["pos"], inplace=True)
            aggregate_C_rolling = aggregate_C.rolling(
                window=smooth, min_periods=min_periods, center=True, on="pos"
            ).mean()
            sns.lineplot(
                x=aggregate_C_rolling["pos"],
                y=aggregate_C_rolling["frac"],
                color=colors[i],
            )
        i = i + 1
    plt.title(basemod)
    plt.legend(sampleNames)
    plt.show()
    fig.savefig(outDir + "/" + basemod + "_sm_rolling_avg_overlay.pdf")


def main():
    # TODO add argument parsing
    print("main")
