"""
Here's the thought: what if we required that the user pass in a well-formatted bed file?

The idea would be that this method expects regions to be all the same size and set up correctly. This would allow someone to parse the entire genome, then extract only the relevant regions.
Would require there to be separate infrastructure for defining the bed regions correctly, which might be interesting?

Takes the onus of doing centering, etc. away from this method. I think this makes things much cleaner.
"""
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes

from . import utils
from . import load_processed


def plot_enrichment_profile_base(mod_file_names: list[str | Path],
                                 bed_file_names: list[str | Path],
                                 mod_names: list[str],
                                 sample_names: list[str],
                                 window_size: int,
                                 smooth_window: int | None = None) -> Axes:
    """
    Plots enrichment profiles, overlaying the resulting traces on top of each other

    Input files should contain sufficient information to generate modification pileup data.
    All input lists are parallel; each index represents the settings for extracting and plotting one sample.

    This is the most flexile method. Say, for example, you have 3 traces to plot. Each modification file contains different types of modification. Each trace is defined by different regions. You want to name the traces something random. This will let you do all of that at once.

    TODO: Clarify this documentation it's a mess. How do I say this concisely?
    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: This is set up in such a way that you may be required to open the same files multiple times. Depending on the file type and loading operations, this could result in unnecessary slowdown.
    TODO: I think it's reasonable for smoothing min_periods to be always set to 1 for this method, as it's a visualization tool, not quantitative. Is this unreasonable?
    TODO: Should the more restrictive meta versions allow *args, or only **kwargs?
    TODO: It's mildly confusing that there are required args that are only seen as *args or **kwargs in the more restrictive meta versions... But this is so much cleaner...

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: TODO: Documentation for this; I think it's a half-size?
        smooth_window: size of the moving window to use for smoothing. If set to None, no smoothing is performed
    
    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(mod_file_names, bed_file_names, mod_names, sample_names):
        raise ValueError('Unequal number of inputs')
    mod_file_names = [Path(fn) for fn in mod_file_names]
    bed_file_names = [Path(fn) for fn in bed_file_names]

    trace_vectors = []
    for mod_file, bed_file, mod_name in zip(mod_file_names, bed_file_names, mod_names):
        match mod_file.suffix:
            case '.gz':
                trace = load_processed.vector_from_bedmethyl(bedmethyl_file=mod_file,
                                                             bed_file=bed_file,
                                                             mod_name=mod_name,
                                                             window_size=window_size)
            case '.fake':
                trace = load_processed.vector_from_fake(mod_file=mod_file,
                                                        bed_file=bed_file,
                                                        mod_name=mod_name,
                                                        window_size=window_size)
            case _:
                raise ValueError(f'Unsupported file type for {mod_file}')
        if smooth_window is not None:
            trace = utils.smooth_rolling_mean(trace, window=smooth_window)
        trace_vectors.append(trace)
    
    axes = utils.line_plot(x=np.arange(-window_size,window_size),
                           x_label='pos',
                           vectors=trace_vectors,
                           vector_names=sample_names,
                           y_label='fraction modified bases')
    return axes

def plot_enrichment_profile_vary_mod(mod_file_name: str | Path,
                                     bed_file_name: str | Path,
                                     mod_names: list[str],
                                     *args,
                                     **kwargs) -> Axes:
    """
    Plot enrichment profile, holding modification file and regions constant, varying modification types
    """
    n_mods = len(mod_names)
    return plot_enrichment_profile_base(mod_file_names=[mod_file_name] * n_mods,
                                        bed_file_names=[bed_file_name] * n_mods,
                                        mod_names=mod_names,
                                        sample_names=mod_names,
                                        *args,
                                        **kwargs)

def plot_enrichment_profile_vary_regions(mod_file_name: str | Path,
                                         bed_file_names: list[str | Path],
                                         mod_name: str,
                                         sample_names: list[str] = None,
                                         *args,
                                         **kwargs) -> Axes:
    """
    Plot enrichment profile, holding modification file and modification types constant, varying regions

    Sample names default to the names of the bed files (?)
    """
    if sample_names is None:
        sample_names = bed_file_names
    n_beds = len(bed_file_names)
    return plot_enrichment_profile_base(mod_file_names=[mod_file_name] * n_beds,
                                        bed_file_names=bed_file_names,
                                        mod_names=[mod_name] * n_beds,
                                        sample_names=sample_names,
                                        *args,
                                        **kwargs)

def plot_enrichment_profile_vary_experiments(mod_file_names: list[str | Path],
                                             bed_file_name: str | Path,
                                             mod_name: str,
                                             sample_names: list[str] = None,
                                             *args,
                                             **kwargs) -> Axes:
    """
    Plot enrichment profile, holding modification types and regions constant, varying modification files

    Sample names default to the names of the modification files (?)

    TODO: This name stinks
    """
    if sample_names is None:
        sample_names = mod_file_names
    n_mod_files = len(mod_file_names)
    return plot_enrichment_profile_base(mod_file_names=mod_file_names,
                                        bed_file_names=[bed_file_name] * n_mod_files,
                                        mod_names=[mod_name] * n_mod_files,
                                        sample_names=sample_names,
                                        *args,
                                        **kwargs)
