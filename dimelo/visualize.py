import sys

import matplotlib.pyplot as plt
import seaborn as sns

from dimelo.parse_bam import parse_ont_bam
from dimelo.plotting_helper_functions import meth_browser

COLOR_A = "#053C5E"
COLOR_C = "#BB4430"


class Region(object):
    def __init__(self, region, fasta=None):
        if ":" in region:
            try:
                self.chromosome, interval = region.replace(",", "").split(":")
                self.begin, self.end = [int(i) for i in interval.split("-")]
            except ValueError:
                sys.exit(
                    "\n\nERROR: Window (-w/--window) inproperly formatted, "
                    "examples of accepted formats are:\n"
                    "'chr5:150200605-150423790'\n\n"
                )
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"


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
    dotsize=0.5,
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
            :param dotsize: size of points
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
        + " reads with methylation > "
        + thresh
        + " for "
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


def browser_sm_roi(
    fileNames,
    sampleNames,
    window,
    basemod,
    outDir,
    thresh,
    bedFileFeatures=None,
    colorA=COLOR_A,
    colorC=COLOR_C,
    dotsize=2,
    static=False,
):
    """
    Create single molecule plots within a region of interest
    Args:
            :param fileNames: list of names of bam files with Mm and Ml tags
            :param sampleNames: list of names of samples for output plot name labelling
            :param window: formatted as for example: "chr1:1-100000"
            :param basemod: which basemods, currently supported options are 'A', 'CG', 'A+CG'
            :param outDir: directory to output plot
            :param thresh: threshold for coloring a position as methylated
            :param bedFileFeatures: annotation to display in browser (optional); default None
            :param colorA: color in hex for mA
            :param colorC: color in hex for mCG
            :param dotsize: size of points; default 2
            :param static: produce pdf if True, produe html if False; default False
    Return:
            plot of single molecules centered at region of interest
    """

    w = Region(window)

    all_data = [
        parse_ont_bam(f, n, basemod=basemod, region=w)
        for f, n in zip(fileNames, sampleNames)
    ]

    meth_browser(
        meth_data=all_data,
        window=Region(window),
        outDir=outDir,
        bed=bedFileFeatures,
        dotsize=dotsize,
        static=static,
        thresh=thresh,
        colorA=colorA,
        colorC=colorC,
    )
