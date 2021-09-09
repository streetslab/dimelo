import sys

from dimelo.parse_bam import parse_bam
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


def browser_sm_roi(
    fileNames,
    sampleNames,
    window,
    basemod,
    outDir,
    threshA=128,
    threshC=128,
    bedFileFeatures=None,
    colorA=COLOR_A,
    colorC=COLOR_C,
    dotsize=4,
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
            :param threshA: threshold for calling mA; default 128
            :param threshC: threshold for calling mCG; default 128
            :param bedFileFeatures: annotation to display in browser (optional); default None
            :param colorA: color in hex for mA
            :param colorC: color in hex for mCG
            :param dotsize: size of points; default 4
            :param static: produce pdf if True, produe html if False; default False
    Return:
            plot of single molecules centered at region of interest
    """

    w = Region(window)

    all_data = [
        parse_bam(f, n, basemod=basemod, region=w)
        for f, n in zip(fileNames, sampleNames)
    ]

    meth_browser(
        meth_data=all_data,
        window=Region(window),
        outDir=outDir,
        bed=bedFileFeatures,
        dotsize=dotsize,
        static=static,
        threshA=threshA,
        threshC=threshC,
        colorA=colorA,
        colorC=colorC,
    )
