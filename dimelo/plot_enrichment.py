from pathlib import Path

from matplotlib.axes import Axes

from . import utils
from . import load_processed


def plot_enrichment(mod_file_names: list[str | Path],
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
                n_mod, n_total = load_processed.counts_from_bedmethyl(bedmethyl_file=mod_file,
                                                                      bed_file=bed_file,
                                                                      mod_name=mod_name)
            case '.fake':
                n_mod, n_total = load_processed.counts_from_fake(mod_file=mod_file,
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


def by_modification(mod_file_name: str | Path,
                    bed_file_name: str | Path,
                    mod_names: list[str],
                    *args,
                    **kwargs) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and regions constant, varying modification types

    TODO: There are no *args or **kwargs currently, but there probably will be for plotting methods?
    """
    n_mods = len(mod_names)
    return plot_enrichment(mod_file_names=[mod_file_name] * n_mods,
                           bed_file_names=[bed_file_name] * n_mods,
                           mod_names=mod_names,
                           sample_names=mod_names,
                           *args,
                           **kwargs)

def by_regions(mod_file_name: str | Path,
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
    return plot_enrichment(mod_file_names=[mod_file_name] * n_beds,
                           bed_file_names=bed_file_names,
                           mod_names=[mod_name] * n_beds,
                           sample_names=sample_names,
                           *args,
                           **kwargs)

def by_dataset(mod_file_names: list[str | Path],
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
    return plot_enrichment(mod_file_names=mod_file_names,
                           bed_file_names=[bed_file_name] * n_mod_files,
                           mod_names=[mod_name] * n_mod_files,
                           sample_names=sample_names,
                           *args,
                           **kwargs)
