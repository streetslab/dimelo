import sqlite3

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from dimelo.parse_bam import parse_bam


def enrich_sm_roi_overlay(
    fileNames,
    sampleNames,
    bedFiles,
    basemod,
    outDir,
    threshA=129,
    threshC=129,
    windowSize=1000,
    colors=["#2D1E2F", "#A9E5BB", "#610345", "#559CAD", "#5E747F"],
    dotsize=0.5,
    smooth=50,
    min_periods=10,
):
    """Create single molecule plots centered at region of interest.
    Overlay EITHER multiple files OR multiple bed regions for a single file
    Args:
            :param fileNames: list names of bam files to overlay with Mm and Ml tags
            :param sampleNames: list of names of sample for output file name labelling
            :param bedFiles: list of specified windows for regions of interest to overlay
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param outDir: directory to output plot
            :param threshA: threshold for calling mA; default 129
            :param threshC: threshold for calling mCG; default 129
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
    if len(fileNames) > 1:
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
            )
    if len(bedFiles) > 1:
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
            )
    plt.title(basemod)
    plt.legend(sampleNames)
    plt.show()
    fig.savefig(outDir + "/" + basemod + "_sm_rolling_avg_overlay.pdf")


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
            color=color,
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
            color=color,
        )


def main():
    # TODO add argument parsing
    print("main")
