import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

from dimelo.parse_bam import parse_bam

COLOR_A = "#053C5E"
COLOR_C = "#BB4430"


def enrich_sm_roi(
    fileName,
    sampleName,
    bedFile,
    basemod,
    outDir,
    threshA=128,
    threshC=128,
    windowSize=1000,
    colorA=COLOR_A,
    colorC=COLOR_C,
    dotsize=0.5,
    smooth=50,
    min_periods=10,
):
    """Create single molecule plots centered at region of interest.
    Args:
            :param fileName: name of bam file with Mm and Ml tags
            :param sampleName: name of sample for output file name labelling
            :param bedFile: specified windows for regions of interest
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param outDir: directory to output plot
            :param threshA: threshold for calling mA; default 128
            :param threshC: threshold for calling mCG; default 128
            :param windowSize: window size around center point of feature of interest to plot (+/-); default 1000 bp
            :param colorA: color in hex for mA
            :param colorC: color in hex for mCG
            :param dotsize: size of points
            :param smooth: window over which to smooth aggregate curve; default of 50 bp
            :param min_periods: minimum number of bases to consider for smoothing: default of 10 bp
    Return:
            plot of single molecules centered at region of interest
    """
    all_data, ave_dict = parse_bam(
        fileName,
        sampleName,
        bedFile,
        basemod,
        center=True,
        windowSize=windowSize,
        threshA=threshA,
        threshC=threshC,
    )

    # only keep calls with probability above threshold
    # now only storing above threshold
    # all_data_t = all_data[all_data["prob"] >= thresh]
    # all_data_t = all_data[
    #     ((all_data["prob"] >= threshA) & (all_data["mod"].str.contains("A")))
    #     | ((all_data["prob"] >= threshC) & (all_data["mod"].str.contains("C")))
    # ]

    # print(
    #     "processing "
    #     + str(len(all_data["read_name"].unique()))
    #     + " reads for "
    #     + sampleName
    #     + " for bam: "
    #     + fileName
    # )
    print(
        "processing "
        + str(len(all_data["read_name"].unique()))
        + " reads with methylation above threshold for "
        + sampleName
        + " for bam: "
        + fileName
    )

    fig = plt.figure()

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

    plt.yticks([])
    plt.ylabel("")
    plt.xlabel("")
    plt.xlim(-windowSize, windowSize)
    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_sm_scatter.png", dpi=600
    )

    plot_aggregate_me_frac(
        sampleName, ave_dict, smooth, min_periods, windowSize, basemod, outDir
    )

    return all_data, ave_dict


# average profile
def plot_aggregate_me_frac(
    sampleName, ave_dict, smooth, min_periods, windowSize, basemod, outDir
):
    fig = plt.figure()
    ave_df = pd.DataFrame.from_dict(
        ave_dict, orient="index", columns=["mod_count", "total_count"]
    )  # keys are rows
    ave_df["pos"] = ave_df.index.to_series().str.split(":").str[0].astype(int)
    ave_df["mod"] = ave_df.index.to_series().str.split(":").str[1]
    ave_df["frac"] = ave_df["mod_count"] / ave_df["total_count"]

    r = range(-windowSize, windowSize + 1, 1)
    if "A" in basemod:
        frac_A = ave_df[ave_df["mod"].str.contains("A")]
        for bp in r:
            if bp not in frac_A["pos"].values:
                df2 = {
                    "pos": bp,
                    "mod": "A",
                    "mod_count": 0,
                    "total_count": 0,
                    "frac": 0,
                }
                frac_A = frac_A.append(df2, ignore_index=True)
        frac_A.sort_values(by=["pos"], inplace=True)
        frac_A_rolling = (
            frac_A["frac"]
            .rolling(window=smooth, min_periods=min_periods, center=True)
            .mean()
        )
        sns.lineplot(x=r, y=frac_A_rolling, color="#053C5E")
    if "C" in basemod:
        frac_C = ave_df[ave_df["mod"].str.contains("C")]
        for bp in r:
            if bp not in frac_C["pos"].values:
                df2 = {
                    "pos": bp,
                    "mod": "C",
                    "mod_count": 0,
                    "total_count": 0,
                    "frac": 0,
                }
                frac_C = frac_C.append(df2, ignore_index=True)
        frac_C.sort_values(by=["pos"], inplace=True)
        frac_C_rolling = (
            frac_C["frac"]
            .rolling(window=smooth, min_periods=min_periods, center=True)
            .mean()
        )
        sns.lineplot(x=r, y=frac_C_rolling, color="#BB4430")
    plt.title(basemod)
    plt.show()
    fig.savefig(
        outDir + "/" + sampleName + "_" + basemod + "_sm_rolling_avg.pdf"
    )

    if "A" in basemod:
        plot_base_abundance(
            sampleName,
            frac_A,
            "A",
            windowSize,
            outDir,
        )
    if "C" in basemod:
        plot_base_abundance(
            sampleName,
            frac_C,
            "CG",
            windowSize,
            outDir,
        )


def plot_base_abundance(sampleName, ave_df, basemod, windowSize, outDir):
    cmapPurple = colors.LinearSegmentedColormap.from_list(
        "custom purple", ["white", "#2D1E2F"], N=200
    )
    # plot base abundance
    fig = plt.figure()
    x = np.linspace(-windowSize, windowSize, num=2 * windowSize + 1)
    y = ave_df["total_count"].to_numpy()  # base_count
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
