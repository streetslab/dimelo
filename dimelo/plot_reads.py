"""
I'm conflicted about how to handle some of this.

There are two different ways of doing single read plotting: "rectangular" and "whole read".
"rectangular" means displaying exactly the requested region.
"whole read" means displaying the entirety of any read overlapping the requested region.
Probably need separate methods for all of this? Is there shared functionality? Do they live in the same file? Etc.

I'm beginning to lose the thread of where we check for regions making sense.
Maybe this is an argument for an internal region class that makes checking easy? I don't know.
"""
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
import matplotlib.pyplot as plt

from . import load_processed


def plot_reads(mod_file_name: str | Path,
               bed_file_name: str | Path,
               mod_names: list[str],
               window_size: int = 0) -> Axes:
    """
    Plots centered single reads as a scatterplot, cut off at the boundaries of the requested regions?

    TODO: Clarify this documentation it's a mess. How do I say this concisely?
    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: This name stinks?
    TODO: So far, this is the only method to do plotting without utility methods. Is this reasonable? Is it that unique?

    Args:
        mod_file_name: path to file containing modification data for single reads
        bed_file_name: path to bed file specifying regions (WHAT DO THESE REPRESENT???)
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files

    Returns:
        Axes object containing the plot
    """ 
    mod_file_name = Path(mod_file_name)
    bed_file_name = Path(bed_file_name)

    match mod_file_name.suffix:
        case _:
            reads, read_names, mods = load_processed.reads_from_hdf5(
                file=mod_file_name,
                bed_file=bed_file_name,
                mod_names=mod_names,
                window_size=window_size,
            )

    # Convert data frame where each row represents a read to a data frame where each row represents a single modified position in a read
    df = pd.DataFrame({
        'read_name': read_names,
        'mod': mods,
        'pos': reads
    }).explode('pos')
    axes = sns.scatterplot(
        data=df,
        x="pos",
        y="read_name",
        hue="mod",
        # palette=colors,
        s=0.5,
        marker="s",
        linewidth=0,
    )
    # Retrieve the existing legend
    legend = axes.legend_

    # Update legend properties
    legend.set_title('Mod')
    for handle in legend.legendHandles:
        handle.set_markersize(10)  # Set a larger marker size for legend
#     handles,labels = axes.get_legend_handles_labels()
#     # Create a new legend with larger marker size
#     new_handles = [plt.Line2D([], [], 
#                               marker='s', 
#                               linestyle='', 
#                               color=handle.get_color(), 
#                               markersize=10) 
#                    for handle in handles] # Skip the first handle as it is the legend title
#     axes.legend(new_handles, labels[1:], title='Mod', loc='upper right')
    
    return axes