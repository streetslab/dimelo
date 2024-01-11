from pathlib import Path

from matplotlib.axes import Axes

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
    return (69, 420)

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
                n_mod, n_total = extract_counts_from_bedmethyl(bedmethyl_file=mod_file)
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

# TODO: Wait, this whole file doesn't have modification stuff enabled... It's probably relevant here too?
# def plot_enrichment_vary_mod(mod_file_names: list[str | Path],
#                              sample_names: list[str]) -> Axes:
#     pass