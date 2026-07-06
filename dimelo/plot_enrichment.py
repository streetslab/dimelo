from pathlib import Path

from matplotlib.axes import Axes

from . import load_processed, utils


def plot_enrichment(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    sample_names: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    quiet: bool = False,
    cores: int | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment comparison barplots using the given list of pre-processed input files.

    Each input list is expected to be parallel and the same length. Each index represents one analysis condition across the lists.
    Using the same file for multiple conditions requires adding the same file multiple times, in the appropriate indices.

    This is the most flexible method for enrichment plotting. For most use cases, consider
    using one of the plot_enrichment.by_* methods.

    Args:
        mod_file_names: list of paths to modified base pileup data files
        bed_file_names: list of paths to bed files specifying regions to extract
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling bars in the output; x-axis labels
        window_size: (currently disabled) window around center of region, +-window_size//2
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing. Default to None, which means all available.
        kwargs: other keyword parameters passed through to utils.bar_plot

    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(mod_file_names, regions_list, motifs, sample_names):
        raise ValueError("Unequal number of inputs")

    mod_fractions = get_enrichments(
        mod_file_names=mod_file_names,
        regions_list=regions_list,
        motifs=motifs,
        window_size=window_size,
        single_strand=single_strand,
        quiet=quiet,
        cores=cores,
    )

    axes = make_enrichment_plot(
        mod_fractions=mod_fractions,
        sample_names=sample_names,
        **kwargs,
    )
    return axes


def by_modification(
    mod_file_name: str | Path,
    regions: str | Path | list[str | Path],
    motifs: list[str],
    **kwargs,
) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and regions constant, varying modification types

    See plot_enrichment for details.
    """
    n_mods = len(motifs)
    return plot_enrichment(
        mod_file_names=[mod_file_name] * n_mods,
        regions_list=[regions] * n_mods,
        motifs=motifs,
        sample_names=motifs,
        **kwargs,
    )


def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]] | None = None,
    motif: str | None = None,
    regions: list[str | Path | list[str | Path]] | None = None,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment bar plots, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_enrichment for details.
    """
    if regions is not None:
        if regions_list is not None and regions_list != regions:
            raise ValueError(
                "Pass either regions_list or regions to by_regions, not both with different values."
            )
        regions_list = regions
    if regions_list is None:
        raise ValueError("by_regions requires regions_list (or alias regions).")
    if motif is None:
        raise ValueError("by_regions requires motif.")

    sample_names_for_plot = (
        sample_names
        if sample_names is not None
        else [str(region) for region in regions_list]
    )
    n_beds = len(regions_list)
    return plot_enrichment(
        mod_file_names=[mod_file_name] * n_beds,
        regions_list=regions_list,
        motifs=[motif] * n_beds,
        sample_names=sample_names_for_plot,
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
    Plot enrichment bar plots, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_enrichment for details.
    """
    sample_names_for_plot = (
        sample_names
        if sample_names is not None
        else [str(mod_file) for mod_file in mod_file_names]
    )
    n_mod_files = len(mod_file_names)
    return plot_enrichment(
        mod_file_names=mod_file_names,
        regions_list=[regions] * n_mod_files,
        motifs=[motif] * n_mod_files,
        sample_names=sample_names_for_plot,
        **kwargs,
    )


def get_enrichments(
    mod_file_names: list[str | Path],
    regions_list: list[str | Path | list[str | Path]],
    motifs: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    quiet: bool = False,
    cores: int | None = None,
) -> list[float]:
    """
    Get the enrichment values, ready for plotting.

    This helper function can be useful during plot prototyping, when repeatedly building plots from the same data.
    Its outputs can be passed as the first argument to make_enrichment_plot().

    Args:
        mod_file_names: list of paths to modified base pileup data files
        regions_list: list of paths to bed files specifying regions to extract
        motifs: list of modifications to extract; expected to match mods available in the relevant mod_files
        window_size: (currently disabled) window around center of region, +-window_size//2
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        quiet: disables progress bars
        cores: CPU cores across which to parallelize processing

    Returns:
        List of modified fraction values.
    """
    if not utils.check_len_equal(mod_file_names, regions_list, motifs):
        raise ValueError("Unequal number of inputs")
    mod_file_paths = [Path(fn) for fn in mod_file_names]

    mod_fractions = []
    for mod_file, regions, motif in zip(
        mod_file_paths, regions_list, motifs, strict=False
    ):
        match mod_file.suffix:
            case ".gz":
                n_mod, n_total = load_processed.pileup_counts_from_bedmethyl(
                    bedmethyl_file=mod_file,
                    regions=regions,
                    motif=motif,
                    window_size=window_size,
                    single_strand=single_strand,
                    quiet=quiet,
                    cores=cores,
                )
            case ".fake":
                n_mod, n_total = load_processed.counts_from_fake(
                    mod_file=mod_file, regions=regions, motif=motif
                )
            case _:
                raise ValueError(f"Unsupported file type for {mod_file}")
        try:
            mod_fractions.append(n_mod / n_total)
        except ZeroDivisionError:
            mod_fractions.append(0)

    return mod_fractions


def make_enrichment_plot(
    mod_fractions: list[float],
    sample_names: list[str],
    **kwargs,
) -> Axes:
    """
    Plot the given enrichment values.

    This helper function can be useful during plot prototyping, when repeatedly building plots from the same data.
    The first argument should be the output of get_enrichments().

    Args:
        mod_fractions: list of modified fraction values.
        sample_names: list of names to use for labeling bars in the output; x-axis labels
        kwargs: other keyword parameters passed through to utils.bar_plot

    Returns:
        Axes object containing the plot
    """
    if not utils.check_len_equal(mod_fractions, sample_names):
        raise ValueError("Unequal number of inputs")

    axes = utils.bar_plot(
        categories=sample_names,
        values=mod_fractions,
        y_label="fraction modified bases",
        **kwargs,
    )

    return axes
