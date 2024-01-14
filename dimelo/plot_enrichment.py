from pathlib import Path

import numpy as np
from matplotlib.axes import Axes
import pysam

import utils
import test_data


def extract_counts_from_bedmethyl(bedmethyl_file: Path,
                                  bed_file: Path,
                                  mod_name: str) -> tuple[int, int]:
    """
    Extract number of modified bases and total number of bases from the given bedmethyl file

    TODO: How to name this method?
    TODO: I feel like stuff like this should be shared functionality
    TODO: Stub; implement this
    
    Args:
        bedmethyl_file: Path to bedmethyl file
        bed_file: Path to bed file specifying regions
        mod_name: type of modification to extract data for
    
    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """
    # Deleted the windowing stuff from the plot_enrichment_profile

    source_tabix = pysam.TabixFile(str(bedmethyl_file))
    # Don't need vectors, just need counts; also not guaranteed that windows are the same length?
    # valid_base_counts = np.zeros(window_size*2)
    # modified_base_counts = np.zeros(window_size*2)
    valid_base_count = 0
    modified_base_count = 0
    
    mod_motif = mod_name.split(',')[0]
    mod_coord_in_motif = mod_name.split(',')[1]

    with open(bed_file) as regions_file:
        for line in regions_file:
            # pythonic condensed operation for extracting important fields
            # fields = line.split('\t')
            chromosome, start_pos, end_pos, *_ = line.split('\t')
            start_pos = int(start_pos)
            end_pos = int(end_pos)
            # Removing centering
            # center_coord = (int(fields[2])+int(fields[1]))//2
            # chromosome = fields[0]
            if chromosome in source_tabix.contigs:
                # Removing centering
                # if center_coord-window_size>0:
                # for row in source_tabix.fetch(chromosome,center_coord-window_size,center_coord+window_size):
                for row in source_tabix.fetch(chromosome, start_pos, end_pos):
                    tabix_fields = row.split('\t')
                    pileup_basemod = tabix_fields[3]
                    if mod_motif in pileup_basemod and mod_coord_in_motif in pileup_basemod:
                        pileup_info = tabix_fields[9].split(' ')
                        # Removing centering
                        # pileup_coord_relative = int(tabix_fields[1])-center_coord+window_size
                        # But actually don't need this because don't care about positions
                        # pileup_coord_relative = int(tabix_fields[1]) - end_pos
                        # valid_base_counts[pileup_coord_relative] += int(pileup_info[0])
                        # modified_base_counts[pileup_coord_relative] += int(pileup_info[2])
                        valid_base_count += int(pileup_info[0])
                        modified_base_count += int(pileup_info[2])

    return (modified_base_count, valid_base_count)

def extract_counts_fake(mod_file: Path,
                        bed_file: Path,
                        mod_name: str) -> tuple[int, int]:
    """
    Generates a fake set of counts
    """
    window_halfsize = 500
    return test_data.fake_peak_counts(halfsize=window_halfsize, peak_height=0.15)

def plot_enrichment_base(mod_file_names: list[str | Path],
                         bed_file_names: list[str | Path],
                         mod_names: list[str],
                         sample_names: list[str]) -> Axes:
    """
    Plots enrichment comparison barplots using the given list of pre-processed input files.

    Input files should contain sufficient information to generate modification pileup data.
    Each file is expected to entirely represent one relevant test condition. All modification events in the file are used.

    For example, each input file might represent different DiMeLo experiments, evaluated in the same regions. Alternatively, each input file might represent the same DiMeLo experiment, evalueted at different sets of regions.

    TODO: Clarify this documentation it's a mess. How do I say this concisely?
    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    
    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying regions to extract
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling bars in the output; x-axis labels

    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(mod_file_names, bed_file_names, mod_names, sample_names):
        raise ValueError('Unequal number of inputs')
    mod_file_names = [Path(fn) for fn in mod_file_names]
    bed_file_names = [Path(fn) for fn in bed_file_names]

    mod_fractions = []
    for mod_file, bed_file, mod_name in zip(mod_file_names, bed_file_names, mod_names):
        match mod_file.suffix:
            case '.gz':
                n_mod, n_total = extract_counts_from_bedmethyl(bedmethyl_file=mod_file,
                                                               bed_file=bed_file,
                                                               mod_name=mod_name)
            case '.fake':
                n_mod, n_total = extract_counts_fake(mod_file=mod_file,
                                                     bed_file=bed_file,
                                                     mod_name=mod_name)
            case _:
                raise ValueError(f'Unsupported file type for {mod_file}')
        try:
            mod_fractions.append(n_mod / n_total)
        except ZeroDivisionError:
            mod_fractions.append(0)
    
    axes = utils.bar_plot(categories=sample_names, values=mod_fractions, y_label='fraction modified bases')
    return axes


def plot_enrichment_vary_mod(mod_file_name: str | Path,
                             bed_file_name: str | Path,
                             mod_names: list[str],
                             *args,
                             **kwargs) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and regions constant, varying modification types

    TODO: There are no *args or **kwargs currently, but there probably will be for plotting methods?
    """
    n_mods = len(mod_names)
    return plot_enrichment_base(mod_file_names=[mod_file_name] * n_mods,
                                bed_file_names=[bed_file_name] * n_mods,
                                mod_names=mod_names,
                                sample_names=mod_names,
                                *args,
                                **kwargs)

def plot_enrichment_vary_regions(mod_file_name: str | Path,
                                 bed_file_names: list[str | Path],
                                 mod_name: str,
                                 sample_names: list[str] = None,
                                 *args,
                                 **kwargs) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and modification types constant, varying regions

    Sample names default to the names of the bed files (?)
    """
    if sample_names is None:
        sample_names = bed_file_names
    n_beds = len(bed_file_names)
    return plot_enrichment_base(mod_file_names=[mod_file_name] * n_beds,
                                bed_file_names=bed_file_names,
                                mod_names=[mod_name] * n_beds,
                                sample_names=sample_names,
                                *args,
                                **kwargs)

def plot_enrichment_vary_experiments(mod_file_names: list[str | Path],
                                     bed_file_name: str | Path,
                                     mod_name: str,
                                     sample_names: list[str] = None,
                                     *args,
                                     **kwargs) -> Axes:
    """
    Plot enrichment bar plots, holding modification types and regions constant, varying modification files

    Sample names default to the names of the modification files (?)

    TODO: This name stinks
    """
    if sample_names is None:
        sample_names = mod_file_names
    n_mod_files = len(mod_file_names)
    return plot_enrichment_base(mod_file_names=mod_file_names,
                                bed_file_names=[bed_file_name] * n_mod_files,
                                mod_names=[mod_name] * n_mod_files,
                                sample_names=sample_names,
                                *args,
                                **kwargs)
