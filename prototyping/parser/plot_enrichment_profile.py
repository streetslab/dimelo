"""
Copied from https://github.com/thekugelmeister/sl_scratch/blob/master/dimelo_2024_plot_interfaces/plot_enrichment_profile.py, Jan 3 2023 12:15pm eastern time

Here's the thought: what if we required that the user pass in a well-formatted bed file?

The idea would be that this method expects regions to be all the same size and set up correctly. This would allow someone to parse the entire genome, then extract only the relevant regions.
Would require there to be separate infrastructure for defining the bed regions correctly, which might be interesting?

Takes the onus of doing centering, etc. away from this method. I think this makes things much cleaner.
"""
from pathlib import Path

import numpy as np
from matplotlib.axes import Axes

import utils
import test_data


""" TEMPORARY STUB VARS """
STUB_HALFSIZE = 100


def extract_vector_from_bedmethyl(bedmethyl_file: Path,
                                  bed_file: Path,
                                  mod_name: str) -> np.ndarray:
    """
    Generate trace for the specified modification aggregated across all regions in the given bed file.

    TODO: How to name this method?
    TODO: I feel like stuff like this should be shared functionality
    TODO: Stub; implement this
    
    Returns:
        vector of values calculated for each position
    """
    return test_data.fake_peak_trace(halfsize=STUB_HALFSIZE)

def get_x_vector_from_bed(bed_file: Path) -> np.ndarray:
    """
    Get a relative coordinate vector from the bed file centered at 0

    TODO: Maybe this can be made generic by providing a centering option? All it does is probably extract the length of the first region; could do so centered or not.
    TODO: Stub; implement this
    """
    return np.arange(STUB_HALFSIZE * 2) - STUB_HALFSIZE


def plot_enrichment_profile_base(mod_file_names: list,
                                 bed_file_names: list,
                                 mod_names: list,
                                 sample_names: list) -> Axes:
    """
    Plots enrichment profiles, overlaying the resulting traces on top of each other

    Input files should contain sufficient information to generate modification pileup data.

    This is the most flexile method. Say, for example, you have 3 traces to plot. Each modification file contains different types of modification. Each trace is defined by different regions. You want to name the traces something random. This will let you do all of that at once.

    TODO: Clarify this documentation it's a mess. How do I say this concisely?
    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: This is set up in such a way that you may be required to open the same files multiple times. Depending on the file type and loading operations, this could result in unnecessary slowdown.
    """
    if not utils.check_len_equal(mod_file_names, bed_file_names, mod_names, sample_names):
        raise ValueError('Unequal number of inputs')
    mod_file_names = [Path(fn) for fn in mod_file_names]
    bed_file_names = [Path(fn) for fn in bed_file_names]

    trace_vectors = []
    for mod_file, bed_file, mod_name in zip(mod_file_names, bed_file_names, mod_names):
        if mod_file.suffix == '.bed':
            trace = extract_vector_from_bedmethyl(bedmethyl_file=mod_file,
                                                      bed_file=bed_file,
                                                      mod_name=mod_name)
        else:
            raise ValueError(f'Unsupported file type for {mod_file}')
        
        trace_vectors.append(trace)
    
    axes = utils.line_plot(x=get_x_vector_from_bed(bed_file=bed_file_names[0]),
                           x_label='pos',
                           vectors=trace_vectors,
                           vector_names=sample_names,
                           y_label='fraction modified bases')
    return axes

def plot_enrichment_profile_vary_mod(mod_file_name: str,
                                     bed_file_name: str,
                                     mod_names: list) -> Axes:
    """
    Plot enrichment profile, holding modification file and regions constant, varying modification types
    """
    n_mods = len(mod_names)
    return plot_enrichment_profile_base(mod_file_names=[mod_file_name] * n_mods,
                                        bed_file_names=[bed_file_name] * n_mods,
                                        mod_names=mod_names,
                                        sample_names=mod_names)

def plot_enrichment_profile_vary_regions(mod_file_name: str,
                                         bed_file_names: list,
                                         mod_name: str,
                                         sample_names: list = None) -> Axes:
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
                                        sample_names=sample_names)

def plot_enrichment_profile_vary_experiments(mod_file_names: list,
                                             bed_file_name: str,
                                             mod_name: str,
                                             sample_names: list = None) -> Axes:
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
                                        sample_names=sample_names)