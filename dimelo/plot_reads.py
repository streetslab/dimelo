"""
I'm conflicted about how to handle some of this.

There are two different ways of doing single read plotting: "rectangular" and "whole read".
"rectangular" means displaying exactly the requested region.
"whole read" means displaying the entirety of any read overlapping the requested region.
Probably need separate methods for all of this? Is there shared functionality? Do they live in the same file? Etc.

I'm beginning to lose the thread of where we check for regions making sense.
Maybe this is an argument for an internal region class that makes checking easy? I don't know.
"""

from collections.abc import Sequence
import math
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes

from . import load_processed, plotting, utils


def _legacy_single_read_axis_spec(
    *,
    relative: bool,
    regions_5to3prime: bool,
    window_size: int | None,
    reads: Sequence[Sequence[object]],
    regions_dict: dict | None,
) -> plotting.AxisSpec:
    if relative:
        if window_size is not None:
            half_window = window_size
        elif regions_dict is not None and len(regions_dict) > 0:
            region_start, region_end, _ = next(iter(regions_dict.values()))[0]
            half_window = (region_end - region_start) // 2
        else:
            half_window = 0
        return plotting.AxisSpec(
            orientation="region_5to3" if regions_5to3prime else "genomic",
            coordinate_mode="fixed_window",
            anchor="center",
            upstream_bp=half_window,
            downstream_bp=half_window,
        )

    absolute_bound = 0.0
    for read_positions in reads:
        for position in read_positions:
            absolute_bound = max(absolute_bound, abs(float(position)))

    return plotting.AxisSpec(
        orientation="genomic",
        coordinate_mode="fixed_window",
        anchor="absolute",
        upstream_bp=math.ceil(absolute_bound),
        downstream_bp=math.ceil(absolute_bound),
    )


def _legacy_single_read_plot_table(
    reads: Sequence[Sequence[object]],
    read_names: Sequence[object],
    mods: Sequence[object],
    *,
    relative: bool,
    regions_5to3prime: bool,
) -> pd.DataFrame:
    plot_table = pd.DataFrame({"read_name": read_names, "mod": mods, "pos": reads}).explode(
        "pos"
    )
    plot_table = plot_table.reset_index(drop=True)
    plot_table["anchor"] = 0.0
    if relative and regions_5to3prime:
        plot_table["region_strand"] = "+"
    return plot_table


def plot_reads(
    mod_file_name: str | Path,
    regions: str | Path | list[str | Path],
    motifs: list[str],
    window_size: int | None = None,
    single_strand: bool = False,
    regions_5to3prime: bool = False,
    sort_by: str | list[str] = "shuffle",
    thresh: float | None = None,
    relative: bool = True,
    **kwargs,
) -> Axes:
    """
    Plots centered single reads as a scatterplot, cut off at the boundaries of the requested regions?

    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: So far, this is the only method to do plotting without utility methods. Is this reasonable? Is it that unique?

    Args:
        mod_file_name: path to file containing modification data for single reads
        regions: path to bed file specifying regions to extract
        motifs: list of modifications to extract; expected to match mods available in the relevant mod_files
        window_size: we plot +-window_size//2 from the center of the region(s)
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        regions_5to3prime: True means negative strand regions get flipped, False means no flipping. Only works if relative=True
        sort_by: ordered list for hierarchical sort. Currently only smallest to biggest.
        thresh: if no threshold has been applied already, this will threshold the mod calls for plotting (method is only boolean)
        relative: if True, all regions are centered

    Returns:
        Axes object containing the plot
    """
    mod_file_name = Path(mod_file_name)
    # bed_file_name = Path(bed_file_name)
    size = kwargs.pop("s", 0.5)

    palette = kwargs.pop("palette", {})

    merged_palette = {**utils.DEFAULT_COLORS, **palette}

    match mod_file_name.suffix:
        # TODO: Fix how the fake reads options work, and make sure they have the same interface as the real ones.
        # dimelo/plot_reads.py:63: error: Argument "regions" to "reads_from_fake" has incompatible type "str | Path | list[str | Path]"; expected "Path"  [arg-type]
        # Will also fix the following error:
        # dimelo/plot_reads.py:68: error: Incompatible types in assignment (expression has type "dict[Any, Any] | None", variable has type "dict[Any, Any]")  [assignment]
        case ".fake":
            reads, read_names, mods, regions_dict = load_processed.reads_from_fake(
                file=mod_file_name,
                regions=regions,
                motifs=motifs,
            )
        case _:
            reads, read_names, mods, regions_dict = (
                load_processed.readwise_binary_modification_arrays(
                    file=mod_file_name,
                    regions=regions,
                    motifs=motifs,
                    window_size=window_size,
                    single_strand=single_strand,
                    regions_5to3prime=regions_5to3prime,
                    thresh=thresh,
                    relative=relative,
                    sort_by=sort_by,
                )
            )

    plot_table = _legacy_single_read_plot_table(
        reads,
        read_names,
        mods,
        relative=relative,
        regions_5to3prime=regions_5to3prime,
    )
    axis = _legacy_single_read_axis_spec(
        relative=relative,
        regions_5to3prime=regions_5to3prime,
        window_size=window_size,
        reads=reads,
        regions_dict=regions_dict,
    )
    prep_kwargs: dict[str, str] = {
        "position_column": "pos",
        "anchor_column": "anchor",
    }
    if axis.orientation == "region_5to3":
        prep_kwargs["region_strand_column"] = "region_strand"

    prepared = plotting.prepare_single_read_plot_data(
        plot_table,
        plot_family="single_read_raster",
        axis=axis,
        **prep_kwargs,
    )
    axes = sns.scatterplot(
        data=prepared["plot_table"],
        x="plot_x",
        y="read_name",
        hue="mod",
        # palette=colors,
        s=size,
        marker="s",
        linewidth=0,
        palette=merged_palette,
        **kwargs,
    )
    # Retrieve the existing legend
    legend = axes.legend_

    # Retrieve legend handles and labels
    handles, labels = axes.get_legend_handles_labels()

    # Update legend properties
    # TODO: Do we need to do this now and after?
    if legend is not None:
        legend.set_title("Mod")

    # Update marker size for all handles
    for handle in handles:
        if hasattr(handle, "set_markersize"):
            handle.set_markersize(10)  # Set a larger marker size for legend

    # Re-apply the legend with updated handles
    # TODO: Is this step necessary?
    axes.legend(handles, labels, title="Mod")

    # TODO: Technically, regions_dict can be None by this point. In that scenario, it will error out when checking the length.
    # It can be None according to type hints but in the actual logical flow I believe it cannot be None
    # However, we can easily just check whether it is None here as well, in case we change behavior elsewhere.
    # Identified with mypy through the following error:
    # dimelo/plot_reads.py:101: error: Argument 1 to "len" has incompatible type "dict[Any, Any] | None"; expected "Sized"  [arg-type]
    if relative and regions_dict is not None and len(regions_dict) > 0:
        region1_start, region1_end, _ = next(iter(regions_dict.values()))[0]
        effective_window_size = (region1_end - region1_start) // 2
        axes.set_xlim([-effective_window_size, effective_window_size])

    return axes
