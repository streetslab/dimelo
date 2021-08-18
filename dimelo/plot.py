r"""
========
plot
========
.. currentmodule:: dimelo.plot
.. autosummary::
  :toctree: _generate/
  plot
  main
"""

import argparse
import os
import sys
import traceback

import pysam

from ._version import __version__
# You can import relative imports from your package. To do this, you can specify
# the relative path of the module you would like to import. In this case,
# we will import the module functions, which is a directory up in PlotMAPQ.
from .functions import PrintDictionaryToTab, SaveMAPQHistogram


def plot_mapq(bamIn, sampleName, dirOut, quiet):
    """Plots MAPQ distributions using matplotlib.
    Args:
      :param bamIn: input bam file
      :param sampleName: name of bam sample
      :param dirOut: output directory to save results to
      :param quiet: whether to output progress to stdout
    """

    mapqDict = {}

    file_bamIn = pysam.AlignmentFile(bamIn, "rb")
    i = 0
    for read in file_bamIn:
        iMAPQ = read.mapping_quality
        mapqDict[iMAPQ] = mapqDict.get(iMAPQ, 0) + 1

        # print progress to stderr if not quiet
        if i % 100000 == 0 and not quiet:
            sys.stderr.write("%i reads processed.\n" % i)
        i += 1

    # make output directory if it does not already exist
    os.makedirs(dirOut, exist_ok=True)

    # Print MAPQ results to tab file
    filePath = os.path.join(dirOut, sampleName + "_MAPQScores.tab")
    PrintDictionaryToTab("MAPQ Score", "Num Reads", mapqDict, filePath)

    # Plot MAPQ results and save
    figPath = os.path.join(dirOut, sampleName + "_MAPQhist.pdf")
    SaveMAPQHistogram(mapqDict, figPath, title="MAPQ: %s" % sampleName)


def main():
    """main() function used to run plot_mapq from the command line."""

    parser = argparse.ArgumentParser(
        usage="%(prog)s [options] alignment_bam sample_name out_directory",
        description="Calculate MAPQ distribution from bam file",
        epilog="Written by Alyssa Morrow ",
    )

    parser.add_argument(
        "-v", "--version", action="store_true", help="Show software version and exit"
    )
    args, argv = parser.parse_known_args()

    # Version is the only case where the bam and out are optional
    if args.version:
        print(__version__)
        sys.exit()

    parser.add_argument("bamIn", type=str, help="Enter the input bam file.")
    parser.add_argument("sampleName", type=str, help="Enter the name of your sample.")
    parser.add_argument(
        "dirOut", type=str, help="Output folder for plot and spreadsheet"
    )

    parser.add_argument(
        "-q", "--quiet", action="store_true", dest="quiet", help="Suppress progress"
    )

    args = parser.parse_args()

    try:
        plot_mapq(args.bamIn, args.sampleName, args.dirOut, args.quiet)
    except Exception:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write(
            "  [Exception type: %s, raised in %s:%d]\n"
            % (
                sys.exc_info()[1].__class__.__name__,
                os.path.basename(traceback.extract_tb(sys.exc_info()[2])[-1][0]),
                traceback.extract_tb(sys.exc_info()[2])[-1][1],
            )
        )
        sys.exit(1)
