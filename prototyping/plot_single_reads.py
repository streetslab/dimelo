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

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes

import utils
import test_data


""" TEMPORARY STUB VARS """
STUB_HALFSIZE = 500
STUB_N_READS = 500


def extract_centered_reads_from_UNKNOWN_FILE_TYPE(file: Path,
                                                  bed_file: Path,
                                                  mod_names: list[str]) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
    """
    TODO: What does the bed file represent in this method? This one is breaking my brain a bit.
    TODO: Variable names in this method stink.
    """
    reads = []
    read_names = []
    mods = []
    for mod_name in mod_names:
        match mod_name:
            case 'A':
                mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'peak') for _ in range(STUB_N_READS)]
            case 'C':
                mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'inverse_peak') for _ in range(STUB_N_READS)]
            case _:
                raise ValueError(f'No stub settings for requested mod {mod_name}')
        reads += mod_reads
        read_names.append(np.arange(len(mod_reads)))
        mods.append([mod_name] * len(mod_reads))
    
    read_names = np.concatenate(read_names)
    mods = np.concatenate(mods)
    return reads, read_names, mods

def plot_single_reads_rectangle(mod_file_name: str | Path,
                                bed_file_name: str | Path,
                                mod_names: list[str]) -> Axes:
    """
    Plots centered single reads as a scatterplot, cut off at the boundaries of the requested regions?

    TODO: Clarify this documentation it's a mess. How do I say this concisely?
    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: This name stinks?
    TODO: So far, this is the only method to do plotting without utility methods. Is this reasonable? Is it that unique?
    """
    mod_file_name = Path(mod_file_name)
    bed_file_name = Path(bed_file_name)

    match mod_file_name.suffix:
        case _:
            reads, read_names, mods = extract_centered_reads_from_UNKNOWN_FILE_TYPE(file=mod_file_name,
                                                                                    bed_file=bed_file_name,
                                                                                    mod_names=mod_names)

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
        legend=None,
    )
    return axes