from pathlib import Path

import numpy as np
from matplotlib.axes import Axes

from . import utils
from . import load_processed


def plot_enrichment_profile(mod_file_names: list[str | Path],
                            bed_file_names: list[str | Path],
                            mod_names: list[str],
                            sample_names: list[str],
                            window_size: int,
                            smooth_window: int | None = None,
                            **kwargs) -> Axes:
    """
    Plot enrichment profiles, overlaying the resulting traces on top of each other.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for enrichment profile plotting. For most use cases, consider
    using one of the plot_enrichment_profile.by_* methods.

    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: I think it's reasonable for smoothing min_periods to be always set to 1 for this method, as it's a visualization tool, not quantitative. Is this unreasonable?
    TODO: Should the more restrictive meta versions allow *args, or only **kwargs?
    TODO: It's mildly confusing that there are required args that are only seen as *args or **kwargs in the more restrictive meta versions... But this is so much cleaner...

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        smooth_window: size of the moving window to use for smoothing. If set to None, no smoothing is performed
        kwargs: other keyword parameters passed through to utils.line_plot
    
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
                modified_base_counts,valid_base_counts = load_processed.pileup_vectors_from_bedmethyl(bedmethyl_file=mod_file,
                                                             bed_file=bed_file,
                                                             mod_name=mod_name,
                                                             window_size=window_size)
                # This probably wants to get changed to NaNs where valid_base_counts==0 but that means fixing smoothing to work with it also
                trace = np.divide(modified_base_counts,valid_base_counts, out=np.zeros_like(modified_base_counts, dtype=float), where=valid_base_counts!=0)
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

    axes = utils.line_plot(indep_vector=np.arange(-len(trace_vectors[0])//2,len(trace_vectors[0])//2+len(trace_vectors[0])%2),
                           indep_name='pos',
                           dep_vectors=trace_vectors,
                           dep_names=sample_names,
                           y_label='fraction modified bases',
                           **kwargs)
    return axes

def by_modification(mod_file_name: str | Path,
                    bed_file_name: str | Path,
                    mod_names: list[str],
                    *args,
                    **kwargs) -> Axes:
    """
    Plot enrichment profile, holding modification file and regions constant, varying modification types

    See plot_enrichment_profile for details.
    """
    n_mods = len(mod_names)
    return plot_enrichment_profile(mod_file_names=[mod_file_name] * n_mods,
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
    Plot enrichment profile, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_enrichment_profile for details.
    """
    if sample_names is None:
        sample_names = bed_file_names
    n_beds = len(bed_file_names)
    return plot_enrichment_profile(mod_file_names=[mod_file_name] * n_beds,
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
    Plot enrichment profile, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_enrichment_profile for details.
    """
    if sample_names is None:
        sample_names = mod_file_names
    n_mod_files = len(mod_file_names)
    return plot_enrichment_profile(mod_file_names=mod_file_names,
                                   bed_file_names=[bed_file_name] * n_mod_files,
                                   mod_names=[mod_name] * n_mod_files,
                                   sample_names=sample_names,
                                   *args,
                                   **kwargs)
