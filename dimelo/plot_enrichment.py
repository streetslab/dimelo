from pathlib import Path

from matplotlib.axes import Axes

from . import utils
from . import load_processed


def plot_enrichment(mod_file_names: list[str | Path],
                    regions_list: list[str | Path | list[str | Path]],
                    motifs: list[str],
                    sample_names: list[str],
                    **kwargs) -> Axes:
    """
    Plot enrichment comparison barplots using the given list of pre-processed input files.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for enrichment plotting. For most use cases, consider
    using one of the plot_enrichment.by_* methods.

    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    
    Args:
        mod_file_names: list of paths to modified base pileup data files
        bed_file_names: list of paths to bed files specifying regions to extract
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling bars in the output; x-axis labels
        kwargs: other keyword parameters passed through to utils.bar_plot

    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(mod_file_names, regions_list, motifs, sample_names):
        raise ValueError('Unequal number of inputs')
    mod_file_names = [Path(fn) for fn in mod_file_names]

    mod_fractions = []
    for mod_file, regions, motif in zip(mod_file_names, regions_list, motifs):
        match mod_file.suffix:
            case '.gz':
                n_mod, n_total = load_processed.pileup_counts_from_bedmethyl(bedmethyl_file=mod_file,
                                                                      regions=regions,
                                                                      motif=motif)
            case '.fake':
                n_mod, n_total = load_processed.counts_from_fake(mod_file=mod_file,
                                                                 regions=regions,
                                                                 motif=motif)
            case _:
                raise ValueError(f'Unsupported file type for {mod_file}')
        try:
            mod_fractions.append(n_mod / n_total)
        except ZeroDivisionError:
            mod_fractions.append(0)
    
    axes = utils.bar_plot(categories=sample_names, values=mod_fractions, y_label='fraction modified bases', **kwargs)
    return axes


def by_modification(mod_file_name: str | Path,
                    regions: str | Path | list[str | Path],
                    motifs: list[str],
                    *args,
                    **kwargs) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and regions constant, varying modification types

    See plot_enrichment for details.
    """
    n_mods = len(motifs)
    return plot_enrichment(mod_file_names=[mod_file_name] * n_mods,
                           regions_list=[regions] * n_mods,
                           motifs=motifs,
                           sample_names=motifs,
                           *args,
                           **kwargs)

def by_regions(mod_file_name: str | Path,
               regions_list: list[str | Path | list[str | Path]],
               motif: str,
               sample_names: list[str] = None,
               *args,
               **kwargs) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and modification types constant, varying regions
    
    Note: Sample names default to the names of the bed files.

    See plot_enrichment for details.
    """
    if sample_names is None:
        sample_names = regions_list
    n_beds = len(regions_list)
    return plot_enrichment(mod_file_names=[mod_file_name] * n_beds,
                           regions_list=regions_list,
                           motifs=[motif] * n_beds,
                           sample_names=sample_names,
                           *args,
                           **kwargs)

def by_dataset(mod_file_names: list[str | Path],
               regions: str | Path | list[str | Path],
               motif: str,
               sample_names: list[str] = None,
               *args,
               **kwargs) -> Axes:
    """
    Plot enrichment bar plots, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_enrichment for details.
    """
    if sample_names is None:
        sample_names = mod_file_names
    n_mod_files = len(mod_file_names)
    return plot_enrichment(mod_file_names=mod_file_names,
                           regions_list=[regions] * n_mod_files,
                           motifs=[motif] * n_mod_files,
                           sample_names=sample_names,
                           *args,
                           **kwargs)
