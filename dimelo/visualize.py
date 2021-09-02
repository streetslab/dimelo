import matplotlib.pyplot as plt
import seaborn as sns

from dimelo.parse_bam import parse_ont_bam

COLOR_A = "#053C5E"
COLOR_C = "#BB4430"


def enrich_sm_roi(
    fileName,
    sampleName,
    bedFile,
    basemod,
    outDir,
    thresh,
    windowSize=1000,
    colorA=COLOR_A,
    colorC=COLOR_C,
    size=0.5,
):
    """Create single molecule plots centered at region of interest.
    Args:
            :param fileName: name of bam file with Mm and Ml tags
            :param sampleName: name of sample for output file name labelling
            :param bedFile: specified windows for regions of interest
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param outDir: directory to output plot
            :param thresh: threshold for coloring a position as methylated
            :param windowSize: window size around center point of feature of interest to plot (+/-); default 1000 bp
            :param colorA: color in hex for mA
            :param colorC: color in hex for mCG
            :param size: size of points
    Return:
            plot of single molecules centered at region of interest
    """
    all_data = parse_ont_bam(
        fileName,
        sampleName,
        bedFile,
        basemod,
        center=True,
        windowSize=windowSize,
    )

    # only keep calls with probability above threshold
    all_data_t = all_data[all_data["prob"] >= thresh]

    fig = plt.figure()

    colors = {"A+Y": colorA, "A+a": colorA, "C+Z": colorC, "C+m": colorC}

    sns.scatterplot(
        data=all_data_t,
        x="pos",
        y="read_name",
        hue="mod",
        palette=colors,
        s=size,
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
# def browser_sm_roi:
