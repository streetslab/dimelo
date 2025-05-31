from pathlib import Path

import numpy as np
from matplotlib.axes import Axes

from . import load_processed, utils


def plot_depth_profile(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    sample_names: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    smooth_window: int | None = None,
    quiet: bool = False,
    cores: int | None = None,
    **kwargs,
) -> Axes:
    """
    Plot depth profiles, overlaying the resulting traces on top of each other.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for depth profile plotting. For most use cases, consider
    using one of the plot_depth_profile.by_* methods.

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
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

    trace_vectors = get_depth_profiles(
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

    axes = make_depth_profile_plot(
        trace_vectors=trace_vectors, sample_names=sample_names, **kwargs
    )
    return axes


def by_modification(
    mod_file_name: str | Path,
    regions: str | Path,
    motifs: list[str],
    **kwargs,
) -> Axes:
    """
    Plot depth profile, holding modification file and regions constant, varying modification types

    See plot_depth_profile for details.
    """
    n_mods = len(motifs)
    return plot_depth_profile(
        mod_file_names=[mod_file_name] * n_mods,
        regions_list=[regions] * n_mods,
        motifs=motifs,
        sample_names=[f"{motif} depth" for motif in motifs],
        **kwargs,
    )


def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]],
    motif: str,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot depth profile, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_depth_profile for details.
    """
    if sample_names is None:
        sample_names = regions_list
    n_beds = len(regions_list)
    return plot_depth_profile(
        mod_file_names=[mod_file_name] * n_beds,
        regions_list=regions_list,
        motifs=[motif] * n_beds,
        sample_names=[f"{sample_name} depth" for sample_name in sample_names],
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
    Plot depth profile, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_depth_profile for details.
    """
    if sample_names is None:
        sample_names = mod_file_names
    n_mod_files = len(mod_file_names)
    return plot_depth_profile(
        mod_file_names=mod_file_names,
        regions_list=[regions] * n_mod_files,
        motifs=[motif] * n_mod_files,
        sample_names=[f"{sample_name} depth" for sample_name in sample_names],
        **kwargs,
    )


def get_depth_profiles(
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
    Get the depth profile traces, ready for plotting.

    This helper function can be useful during plot prototyping, when repeatedly building plots from the same data.
    Its outputs can be passed as the first argument to make_depth_profile_plot().

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping
        smooth_window: size of the moving window to use for smoothing. If set to None, no smoothing is performed
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing

    Returns:
        List of depth profile traces
    """
    if not utils.check_len_equal(mod_file_names, regions_list, motifs):
        raise ValueError("Unequal number of inputs")
    # TODO: redefinition error; still need to figure out how to do this elegantly in a way mypy likes
    # dimelo/plot_depth_profile.py:53: error: Item "str" of "str | Path" has no attribute "suffix"  [union-attr]
    mod_file_names = [Path(fn) for fn in mod_file_names]

    trace_vectors = []
    for mod_file, regions, motif in zip(mod_file_names, regions_list, motifs):
        match mod_file.suffix:
            case ".gz":
                _, valid_base_counts = load_processed.pileup_vectors_from_bedmethyl(
                    bedmethyl_file=mod_file,
                    regions=regions,
                    motif=motif,
                    window_size=window_size,
                    single_strand=single_strand,
                    regions_5to3prime=regions_5to3prime,
                    quiet=quiet,
                    cores=cores,
                )
                trace = valid_base_counts.astype(float)
                trace[trace == 0] = np.nan
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


def make_depth_profile_plot(
    trace_vectors: list[np.ndarray],
    sample_names: list[str],
    **kwargs,
) -> Axes:
    """
    Plot the given depth profile traces.

    This helper function can be useful during plot prototyping, when repeatedly building plots from the same data.
    The first argument should be the output of get_depth_profiles().

    Args:
        trace_vectors: list of depth profile traces
        sample_names: list of names to use for labeling traces in the output; legend entries
        kwargs: other keyword parameters passed through to utils.line_plot

    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(trace_vectors, sample_names):
        raise ValueError("Unequal number of inputs")
    axes = utils.line_plot(
        indep_vector=np.arange(
            -len(trace_vectors[0]) // 2,
            len(trace_vectors[0]) // 2 + len(trace_vectors[0]) % 2,
        ),
        indep_name="pos",
        dep_vectors=trace_vectors,
        dep_names=sample_names,
        y_label="per strand reads\nwith motif and mod info",
        **kwargs,
    )
    return axes
