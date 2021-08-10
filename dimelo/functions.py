r"""
========================
PlotMAPQHelper functions
========================
.. currentmodule:: PlotMAPQ.functions
.. autosummary::
  :toctree: _generate/
  PrintDictionaryToTab
  SaveMAPQHistogram
"""


import matplotlib.pyplot as plt


def PrintDictionaryToTab(keyHeader,
                    valueHeader,
                    dictIn,
                    filePath):
    """ Write out dictIn to a tab-delimited file.
    Args:
        :param keyHeader: title of header for keys
        :param valueHeader: title of header for values
        :param dictIn: dictionary of keys and values
        :param filePath: filePath to save results to
    """

    with open(filePath,'w') as f:
        f.write(keyHeader + "\t" + valueHeader + "\n")
        for key in dictIn.keys():
            f.write(str(key) + "\t" + str(dictIn[key]) + "\n")

def SaveMAPQHistogram(dictIn,
                    filePath,
                    title="MAPQ"):
    """ Plots and saves dictIn to a histogram.
    Args:
        :param dictIn: dictionary of (k,v) where k indicates mapq and v indicates
            frequency.
        :param filePath: filePath to save results to
    """

    # sort MAPQS
    sortedMAPQ = sorted(dictIn.items(), key=lambda x: x[0])

    fig, ax = plt.subplots(1, 1)
    x = [i[0] for i in sortedMAPQ]
    y = [i[1] for i in sortedMAPQ]
    ax.plot(x, y)
    ax.set_title(title)
    fig.savefig(filePath)