from pathlib import Path

import numpy as np
from matplotlib.axes import Axes

from . import load_processed, utils


def plot_enrichment_profile(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    sample_names: list[str],
    window_size: int,
    relative: bool = True,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    smooth_window: int | None = None,
    quiet: bool = False,
    cores: int | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment profiles, overlaying the resulting traces on top of each other.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for enrichment profile plotting. For most use cases, consider
    using one of the plot_enrichment_profile.by_* methods.

    TODO: I think it's reasonable for smoothing min_periods to be always set to 1 for this method, as it's a visualization tool, not quantitative. Is this unreasonable?
    TODO: Should the more restrictive meta versions allow *args, or only **kwargs?
    No, we want to be able to pass kwargs down to the line plotter, I think. Especially if we swap it out for one that takes more different standard args.
    TODO: It's mildly confusing that there are required args that are only seen as *args or **kwargs in the more restrictive meta versions... But this is so much cleaner...

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        relative: True means x-axis is centered around region centers, False means x-axis is absolute genome positions
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        smooth_window: size of the moving window to use for smoothing. If set to None, no smoothing is performed
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        kwargs: other keyword parameters passed through to utils.line_plot

    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(mod_file_names, regions_list, motifs, sample_names):
        raise ValueError("Unequal number of inputs")

    trace_vectors = get_enrichment_profiles(
        mod_file_names=mod_file_names,
        regions_list=regions_list,
        motifs=motifs,
        window_size=window_size,
        single_strand=single_strand,
        regions_5to3prime=regions_5to3prime,
        smooth_window=smooth_window,
        quiet=quiet,
        cores=cores,
    )

    if relative:
        offset_center = 0
    else:
        regions_dict = utils.regions_dict_from_input(
            regions_list[0],
            window_size,
        )
        if len(regions_dict) == 1 and len(list(regions_dict.values())[0]) == 1:
            region_tuple = list(regions_dict.values())[0][0]
            offset_center = (region_tuple[0] + region_tuple[1]) // 2
        else:
            raise ValueError(
                "relative=False must be used when plotting more than one region."
            )

    axes = make_enrichment_profile_plot(
        trace_vectors=trace_vectors,
        sample_names=sample_names,
        offset_center=offset_center,
        **kwargs,
    )
    return axes


def by_modification(
    mod_file_name: str | Path,
    regions: str | Path,
    motifs: list[str],
    **kwargs,
) -> Axes:
    """
    Plot enrichment profile, holding modification file and regions constant, varying modification types

    See plot_enrichment_profile for details.
    """
    n_mods = len(motifs)
    return plot_enrichment_profile(
        mod_file_names=[mod_file_name] * n_mods,
        regions_list=[regions] * n_mods,
        motifs=motifs,
        sample_names=motifs,
        **kwargs,
    )


"""
TODO: Re-assignment issue:
dimelo/plot_enrichment_profile.py:142: error: Incompatible types in assignment (expression has type "list[str | Path | list[str | Path]]", variable has type "list[str] | None")  [assignment]
dimelo/plot_enrichment_profile.py:148: error: Argument "sample_names" to "plot_enrichment_profile" has incompatible type "list[str] | None"; expected "list[str]"  [arg-type]
dimelo/plot_enrichment_profile.py:168: error: Incompatible types in assignment (expression has type "list[str | Path]", variable has type "list[str] | None")  [assignment]
dimelo/plot_enrichment_profile.py:174: error: Argument "sample_names" to "plot_enrichment_profile" has incompatible type "list[str] | None"; expected "list[str]"  [arg-type]

If sample names is None we assign it non-None values, so it's not clear what the problem is to me. We could make an intermediate dummy variable I guess? If that is the complaint?
"""


def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment profile, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_enrichment_profile for details.
    """
    if sample_names is None:
        sample_names = regions_list
    n_beds = len(regions_list)
    return plot_enrichment_profile(
        mod_file_names=[mod_file_name] * n_beds,
        regions_list=regions_list,
        motifs=[motif] * n_beds,
        sample_names=sample_names,
        **kwargs,
    )


def by_dataset(
    mod_file_names: list[str | Path],
    regions: str | Path | list[str | Path],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment profile, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_enrichment_profile for details.
    """
    if sample_names is None:
        sample_names = mod_file_names
    n_mod_files = len(mod_file_names)
    return plot_enrichment_profile(
        mod_file_names=mod_file_names,
        regions_list=[regions] * n_mod_files,
        motifs=[motif] * n_mod_files,
        sample_names=sample_names,
        **kwargs,
    )


def get_enrichment_profiles(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    window_size: int,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    smooth_window: int | None = None,
    quiet: bool = False,
    cores: int | None = None,
) -> list[np.ndarray]:
    """
    Get the enrichment profile traces, ready for plotting.

    This helper function can be useful during plot prototyping, when repeatedly building plots from the same data.
    Its outputs can be passed as the first argument to make_enrichment_profile_plot().

    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: I think it's reasonable for smoothing min_periods to be always set to 1 for this method, as it's a visualization tool, not quantitative. Is this unreasonable?

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing
        smooth_window: size of the moving window to use for smoothing. If set to None, no smoothing is performed

    Returns:
        List of enrichment profile traces
    """
    if not utils.check_len_equal(mod_file_names, regions_list, motifs):
        raise ValueError("Unequal number of inputs")
    # TODO: redefinition error; still need to figure out how to do this elegantly in a way mypy likes
    # dimelo/plot_enrichment_profile.py:53: error: Item "str" of "str | Path" has no attribute "suffix"  [union-attr]
    mod_file_names = [Path(fn) for fn in mod_file_names]

    trace_vectors = []
    for mod_file, regions, motif in zip(mod_file_names, regions_list, motifs):
        match mod_file.suffix:
            case ".gz":
                modified_base_counts, valid_base_counts = (
                    load_processed.pileup_vectors_from_bedmethyl(
                        bedmethyl_file=mod_file,
                        regions=regions,
                        motif=motif,
                        window_size=window_size,
                        single_strand=single_strand,
                        regions_5to3prime=regions_5to3prime,
                        quiet=quiet,
                        cores=cores,
                    )
                )
                # Default to nan so we can skip over unfilled values when plotting or doing a rolling average
                nans_everywhere = np.full_like(
                    modified_base_counts, np.nan, dtype=float
                )
                trace = np.divide(
                    modified_base_counts,
                    valid_base_counts,
                    out=nans_everywhere,
                    where=valid_base_counts != 0,
                )
            case ".fake":
                trace = load_processed.vector_from_fake(
                    mod_file=mod_file,
                    bed_file=regions,
                    motif=motif,
                    window_size=window_size,
                )
            case _:
                raise ValueError(f"Unsupported file type for {mod_file}")
        if smooth_window is not None:
            trace = utils.smooth_rolling_mean(trace, window=smooth_window)
        trace_vectors.append(trace)
    return trace_vectors


def make_enrichment_profile_plot(
    trace_vectors: list[np.ndarray],
    sample_names: list[str],
    offset_center: int = 0,
    **kwargs,
) -> Axes:
    """
    Plot the given enrichment profile traces.

    This helper function can be useful during plot prototyping, when repeatedly building plots from the same data.
    The first argument should be the output of get_enrichment_profiles().

    Args:
        trace_vectors: list of enrichment profile traces
        sample_names: list of names to use for labeling traces in the output; legend entries
        offset_center: position offset to apply to x-axis (e.g., when plotting absolute genome positions)
        kwargs: other keyword parameters passed through to utils.line_plot

    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(trace_vectors, sample_names):
        raise ValueError("Unequal number of inputs")
    axes = utils.line_plot(
        indep_vector=np.arange(
            offset_center - len(trace_vectors[0]) // 2,
            offset_center + len(trace_vectors[0]) // 2 + len(trace_vectors[0]) % 2,
        ),
        indep_name="pos",
        dep_vectors=trace_vectors,
        dep_names=sample_names,
        y_label="fraction modified bases",
        **kwargs,
    )
    return axes
