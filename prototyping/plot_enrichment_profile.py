"""
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

import pysam


""" TEMPORARY STUB VARS """
STUB_HALFSIZE = 100


def extract_vector_from_bedmethyl(bedmethyl_file: Path,
                                  bed_file: Path,
                                  mod_name: str,
                                 window_size: int) -> np.ndarray:
    """
    Generate trace for the specified modification aggregated across all regions in the given bed file.

    TODO: How to name this method?
    TODO: I feel like stuff like this should be shared functionality
    TODO: Stub; implement this
    TODO: I _THINK_ this should return a value for every position, with non-modified positions as zeros

    Args:
        bedmethyl_file: Path to bedmethyl file
        bed_file: Path to bed file specifying centered equal-length regions
        mod_name: type of modification to extract data for
    
    Returns:
        vector of fraction modifiied bases (e.g. mA/A) calculated for each position; float values between 0 and 1
    """
    source_tabix = pysam.TabixFile(bedmethyl_file)
    valid_base_counts = np.zeros(window_size*2)
    modified_base_counts = np.zeros(window_size*2)
    with open(bed_file) as regions_file:
        for line in regions_file:
            fields = line.split('\t')
            center_coord = (int(fields[2])+int(fields[1]))//2
            chromosome = fields[0]
            if chromosome=='chr1':
                if center_coord-window_size>0:
                    for row in source_tabix.fetch(chromosome,center_coord-window_size,center_coord+window_size):
                        tabix_fields = row.split('\t')
#                         print(tabix_fields)
                        pileup_basemod = tabix_fields[3]
                        if mod_name in pileup_basemod:
                            pileup_info = tabix_fields[9].split(' ')
                            pileup_coord_relative = int(tabix_fields[1])-center_coord+window_size
                            valid_base_counts[pileup_coord_relative] += int(pileup_info[0])
                            modified_base_counts[pileup_coord_relative] += int(pileup_info[2])
                        
    
    modified_fractions = np.divide(modified_base_counts,valid_base_counts,where=valid_base_counts!=0)
    return np.nan_to_num(modified_fractions)
#     return test_data.fake_peak_trace(halfsize=STUB_HALFSIZE)

def get_x_vector_from_bed(bed_file: Path) -> np.ndarray[int]:
    """
    Get a relative coordinate vector from the bed file centered at 0

    TODO: Maybe this can be made generic by providing a centering option? All it does is probably extract the length of the first region; could do so centered or not.
    TODO: Stub; implement this

    Args:
        bed_file: Path to bed file specifying centered equal-length regions
    
    Returns:
        Vector of signed integer positions relative to the center of the region
    """
    return np.arange(STUB_HALFSIZE * 2) - STUB_HALFSIZE


def plot_enrichment_profile_base(mod_file_names: list[str | Path],
                                 bed_file_names: list[str | Path],
                                 mod_names: list[str],
                                 sample_names: list[str]) -> Axes:
    """
    Plots enrichment profiles, overlaying the resulting traces on top of each other

    Input files should contain sufficient information to generate modification pileup data.
    All input lists are parallel; each index represents the settings for extracting and plotting one sample.

    This is the most flexile method. Say, for example, you have 3 traces to plot. Each modification file contains different types of modification. Each trace is defined by different regions. You want to name the traces something random. This will let you do all of that at once.

    TODO: Clarify this documentation it's a mess. How do I say this concisely?
    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: This is set up in such a way that you may be required to open the same files multiple times. Depending on the file type and loading operations, this could result in unnecessary slowdown.

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
    
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
            case '.bed':
                trace = extract_vector_from_bedmethyl(bedmethyl_file=mod_file,
                                                      bed_file=bed_file,
                                                      mod_name=mod_name)
            case _:
                raise ValueError(f'Unsupported file type for {mod_file}')
        trace_vectors.append(trace)
    
    axes = utils.line_plot(x=get_x_vector_from_bed(bed_file=bed_file_names[0]),
                           x_label='pos',
                           vectors=trace_vectors,
                           vector_names=sample_names,
                           y_label='fraction modified bases')
    return axes

def plot_enrichment_profile_vary_mod(mod_file_name: str | Path,
                                     bed_file_name: str | Path,
                                     mod_names: list[str]) -> Axes:
    """
    Plot enrichment profile, holding modification file and regions constant, varying modification types
    """
    n_mods = len(mod_names)
    return plot_enrichment_profile_base(mod_file_names=[mod_file_name] * n_mods,
                                        bed_file_names=[bed_file_name] * n_mods,
                                        mod_names=mod_names,
                                        sample_names=mod_names)

def plot_enrichment_profile_vary_regions(mod_file_name: str | Path,
                                         bed_file_names: list[str | Path],
                                         mod_name: str,
                                         sample_names: list[str] = None) -> Axes:
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

def plot_enrichment_profile_vary_experiments(mod_file_names: list[str | Path],
                                             bed_file_name: str | Path,
                                             mod_name: str,
                                             sample_names: list[str] = None) -> Axes:
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
