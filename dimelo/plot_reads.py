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
from . import utils

def plot_reads(mod_file_name: str | Path,
               regions: str | Path | list[str | Path],
               motifs: list[str],
               window_size: int = None,
               sort_by: str | list[str] = 'shuffle',
               thresh: float = None,
               relative: bool = True,
               **kwargs
              ) -> Axes:
    """
    Plots centered single reads as a scatterplot, cut off at the boundaries of the requested regions?

    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: So far, this is the only method to do plotting without utility methods. Is this reasonable? Is it that unique?

    Args:
        mod_file_name: path to file containing modification data for single reads
        regions: path to bed file specifying regions to extract
        motifs: list of modifications to extract; expected to match mods available in the relevant mod_files

    Returns:
        Axes object containing the plot
    """ 
    mod_file_name = Path(mod_file_name)
    # bed_file_name = Path(bed_file_name)
    size = kwargs.pop('s', 0.5)

    palette = kwargs.pop('palette', {})

    merged_palette = {**utils.DEFAULT_COLORS, **palette}




    match mod_file_name.suffix:
        case '.fake':
            reads,read_names,mods,regions_dict = load_processed.reads_from_fake(
                file = mod_file_name,
                regions = regions,
                motifs = motifs,
            )
        case _:
            reads,read_names,mods,regions_dict = load_processed.readwise_binary_modification_arrays(
            file = mod_file_name,
            regions = regions,
            motifs = motifs,
            window_size = window_size,
            thresh = thresh,
            relative = relative,
            sort_by = sort_by,
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
        s=size,
        marker="s",
        linewidth=0,
        palette = merged_palette,
        **kwargs
    )
    # Retrieve the existing legend
    legend = axes.legend_

    # Update legend properties
    if legend is not None:
        legend.set_title('Mod')
        for handle in legend.legendHandles:
            handle.set_markersize(10)  # Set a larger marker size for legend

    if relative and len(regions_dict)>0:
        region1_start,region1_end,_ = next(iter(regions_dict.values()))[0]
        effective_window_size = (region1_end-region1_start)//2
        axes.set_xlim([-effective_window_size,effective_window_size])
    
    return axes