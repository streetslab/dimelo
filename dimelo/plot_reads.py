"""Single-read plotting entrypoints and legacy axis routing helpers."""

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
    plot_table = pd.DataFrame(
        {
            "read_name": read_names,
            "read_index": list(range(len(read_names))),
            "mod": mods,
            "pos": reads,
        }
    ).explode("pos")
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
    legend_title = kwargs.pop("legend_title", "Mod, index")
    x_axis_label = kwargs.pop("x_label", "Position (bp)")
    y_axis_mode = kwargs.pop("y_axis", "read_index")
    y_axis_label = kwargs.pop(
        "y_label",
        "Read index" if y_axis_mode == "read_index" else "Read name",
    )

    palette = kwargs.pop("palette", {})

    merged_palette = {**utils.DEFAULT_COLORS, **palette}

    match mod_file_name.suffix:
        # Keep fake-file handling aligned with the non-fake return contract.
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
    if y_axis_mode not in plot_table.columns:
        raise ValueError(f"Unsupported y_axis {y_axis_mode!r}. Use one of: {sorted(plot_table.columns)}")
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
        y=y_axis_mode,
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

    # Update legend properties.
    if legend is not None:
        legend.set_title(legend_title)

    # Update marker size for all handles
    for handle in handles:
        if hasattr(handle, "set_markersize"):
            handle.set_markersize(10)  # Set a larger marker size for legend

    # Re-apply the legend with updated marker size.
    axes.legend(handles, labels, title=legend_title)

    if hasattr(axes, "set_xlabel"):
        axes.set_xlabel(x_axis_label)
    if hasattr(axes, "set_ylabel"):
        axes.set_ylabel(y_axis_label)

    # regions_dict may be absent for some loader paths, so guard before reading bounds.
    if relative and regions_dict is not None and len(regions_dict) > 0:
        region1_start, region1_end, _ = next(iter(regions_dict.values()))[0]
        effective_window_size = (region1_end - region1_start) // 2
        axes.set_xlim([-effective_window_size, effective_window_size])

    return axes
