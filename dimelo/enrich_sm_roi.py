import matplotlib.pyplot as plt
import seaborn as sns

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
    Return:
            plot of single molecules centered at region of interest
    """
    all_data = parse_bam(
        fileName,
        sampleName,
        bedFile,
        basemod,
        center=True,
        windowSize=windowSize,
    )

    # only keep calls with probability above threshold
    # all_data_t = all_data[all_data["prob"] >= thresh]
    all_data_t = all_data[
        (all_data["prob"] >= threshA & all_data["mod"].str.contains("A"))
        | (all_data["prob"] >= threshC & all_data["mod"].str.contains("C"))
    ]

    print(
        "processing "
        + str(len(all_data["read_name"].unique()))
        + " reads for "
        + sampleName
        + " for bam: "
        + fileName
    )
    print(
        "processing "
        + str(len(all_data_t["read_name"].unique()))
        + " reads with methylation above threshold for "
        + sampleName
        + " for bam: "
        + fileName
    )

    fig = plt.figure()

    colors = {"A+Y": colorA, "A+a": colorA, "C+Z": colorC, "C+m": colorC}

    sns.scatterplot(
        data=all_data_t,
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


# def plot_base_abundance:
# def plot_aggregate_me_frac:


def main():
    # TODO add argument parsing
    print("main")
