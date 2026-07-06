from pathlib import Path

import numpy as np
import pandas as pd
import plotly

from . import load_processed, utils


def plot_read_browser(
    mod_file_name: str | Path,
    region: str,
    motifs: list[str],
    thresh: int | float | None = None,
    single_strand: bool = False,
    sort_by: str | list[str] = "shuffle",
    hover: bool = True,
    subset_parameters: dict | None = None,
    span_full_window: bool = False,
    width: int = 1,
    size: int = 4,
    colorscales: dict = utils.DEFAULT_COLORSCALES,
    **kwargs,
) -> plotly.graph_objs.Figure:
    """
    Plot base modifications on single reads in a high-quality, interactive-enabled fashion.

    This method returns a plotly Figure object, which can be used in a number of ways to view and save
    the figure in different formats. To view the figure interactively (in a notebook or python script),
    simply call the show() method of the returned Figure object. See the helper methods below for saving
    figures.

    Additional keyword arguments will be passed down to collapse_rows() if sort_by == "collapse". See that
    method for details.

    Args:
        mod_file_name: path to file containing modification data for single reads
        region: region string specifying the region to plot
        motifs: list of modifications to extract; expected to match mods available in the relevant mod_files
        thresh: A modification calling threshold. While the browser always displays float probabilities, setting
            this to a value will limit display to only modification events over the given threshold. Else, display
            all modifications regardless of probability.
        single_strand: True means we only grab counts from reads from the same strand as
            the region of interest, False means we always grab both strands within the regions
        sort_by: ordered list for hierarchical sort; see load_processed.read_vectors_from_hdf5() for details.
            Can also pass the argument "collapse" to allow multiple reads on single rows of the browser, for a
            more condensed visualization. Note that "collapse" is mutually exclusive with all other sorting options,
            and is only allowed to be passed as a single string option.
        hover: if False, disables display of information on mouse hover
        subset_parameters: Parameters to pass to the utils.random_sample() method, to subset the
            reads to be returned. If not None, at least one of n or frac must be provided.
        span_full_window: if True, only plot reads that fully span the window defined by region_start-region_end
        width: width of the read lines in the browser
        size: size of the modification event markers in the browser
        colorscales: dictionary mapping motif names to plotly colorscale specifications

    Returns:
        plotly Figure object containing the plot

    TODO: Improve color specification? User should be able to set their own colors.
    TODO: Should this let the user set arbitrary thresholds for each motif individually?
    TODO: The way that "collapse" is specified is unintuitive and problematic; what if the user passes "collapse" as
        an element in an array?
    TODO: Is it worth having an option for meta-sorting of collapsed reads? It's here for now to enable testing.
    """
    # If asked to collapse reads, set up the initial read sorting appropriately and prep for later
    collapse = False
    if sort_by == "collapse":
        collapse = True
        sort_by = "read_start"

    read_tuples, entry_labels, _ = load_processed.read_vectors_from_hdf5(
        file=mod_file_name,
        regions=region,
        motifs=motifs,
        single_strand=single_strand,
        sort_by=sort_by,
        calculate_mod_fractions=True,
        subset_parameters=subset_parameters,
        span_full_window=span_full_window,
    )

    mod_vector_index = entry_labels.index("mod_vector")

    if read_tuples[0][mod_vector_index].dtype == np.bool_:
        raise ValueError(
            "A threshold has been applied to this .h5 single read data. plot_read_browser must be used with an .h5 file extracted using thresh=None."
        )

    read_extent_df, mod_event_df = format_browser_data(
        read_tuples=read_tuples, entry_labels=entry_labels
    )

    # Apply threshold to mod_event_df
    if thresh is not None:
        mod_event_df = mod_event_df[mod_event_df.prob > utils.adjust_threshold(thresh)]
    else:
        # Still need to filter out all values that are effectively 0, or the read bars cannot be seen
        # TODO: This seems like the wrong place to be handling this.
        mod_event_df = mod_event_df[mod_event_df.prob > utils.adjust_threshold(2)]

    try:
        chrom, (region_start, region_end, _) = utils.parse_region_string(
            region=region, window_size=None
        )
    except ValueError as err:
        raise ValueError(
            "Invalid region specification: plot_read_browser requires a single genomic locus."
        ) from err

    fig = make_browser_figure(
        read_extent_df=read_extent_df,
        mod_event_df=mod_event_df,
        collapse=collapse,
        chrom=chrom,
        region_start=region_start,
        region_end=region_end,
        hover=hover,
        thresh=thresh,
        width=width,
        size=size,
        colorscales=colorscales,
        **kwargs,
    )

    return fig


def format_browser_data(
    read_tuples: list[tuple],
    entry_labels: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Take data from load_processed.read_vectors_from_hdf5() and format it for browser plotting.

    Argument descriptions taken directly from documentation of load_processed.read_vectors_from_hdf5().

    Args:
        read_tuples: a list of tuples, each tuple containing all datasets corresponding to an individual read that
        was within the specified regions.
        entry_labels: a list of strings, naming the datasets returned.

    Returns:
        * dataframe defining the start, end, name, and desired y-index (sorting) of each read
        * dataframe defining all necessary information to place each modification event on the browser

    TODO: There is an observed issue where there are duplicated reads; duplicates are currently thrown away. This only
        manifests itself when subsetting reads, as different random subsets of the same requested size will show up as
        different numbers of rows.
    """
    # Coerce read tuples to initial dataframe, throwing away unnecessary columns
    to_exclude = ["chromosome", "strand", "region_start", "region_end"]
    read_df = pd.DataFrame.from_records(
        read_tuples, columns=entry_labels, exclude=to_exclude
    )

    # For each row, pull out just the positions of valid bases and the probabilities at those positions.
    # Iterate over NumPy-backed column arrays to avoid iterrows() per-row Series overhead.
    mod_vectors = read_df["mod_vector"].to_numpy(copy=False)
    val_vectors = read_df["val_vector"].to_numpy(copy=False)
    read_starts = read_df["read_start"].to_numpy(copy=False)
    prob_vectors = []
    pos_vectors = []
    for mod_vector, val_vector, read_start in zip(
        mod_vectors, val_vectors, read_starts, strict=False
    ):
        selected_indices = np.flatnonzero(val_vector == 1)
        prob_vectors.append(mod_vector[selected_indices])
        pos_vectors.append(selected_indices + read_start)
    read_df["prob_vector"] = prob_vectors
    read_df["pos_vector"] = pos_vectors

    # assign reads y-axis values based on their unique names, preserving the pre-set order
    # TODO: I'm pretty sure that however I do this, there will possibly be an "overdrawing" problem if the same read shows up multiple times. Do I care? Will it be taken care of by dropping duplicates?
    read_df["y_index"] = read_df.read_name.map(
        {k: v for v, k in enumerate(read_df.read_name.unique())}
    )

    # TODO: Dropping duplicates should hide the "overdrawing" problem when reads are duplicated in the dataset. Is this a problem or the intended behavior?
    # Get two separate dataframes:
    # * represents the read extents, to draw the grey lines

    # Added the read_name subsetting to avoid an index mapping error in make_browser_figure
    # caused by the read metadata for the same read with different motifs having slight
    # start and end offsets in some cases, when using Dorado or PacBio data. These 1-5bp
    # differences are caused by the assumptions made in parse_bam.extract h5 conversion
    # and don't matter for much but do cause the full duplicate check to leave duplicates,
    # which then cause non-unique indices for mapping in collapse mode
    # The sorting below is just to make sure that when we drop duplicates, we keep
    # the longest read extent for each read_name. It doesn't effect final figure sort order.
    read_extent_df = (
        read_df[["read_start", "read_end", "y_index", "read_name"]]
        .assign(read_length=lambda df: df.read_end - df.read_start)
        .sort_values(
            "read_length", ascending=False
        )  # TODO: logic can be removed once we reference true read lengths from modkit
        .drop_duplicates(subset=["read_name"])
        .drop(columns=["read_length"])
    )
    # * represents the methylation events
    mod_event_df = (
        read_df[["y_index", "read_name", "motif", "pos_vector", "prob_vector"]]
        .explode(["pos_vector", "prob_vector"])
        .rename(columns={"pos_vector": "pos", "prob_vector": "prob"})
    )

    return read_extent_df, mod_event_df


def collapse_rows(
    read_extent_df: pd.DataFrame,
    minimum_gap: int = 500,
    meta_sort: str | None = "full_extent",
) -> pd.Series:
    """
    Takes a sorted dataframe of read extents and collapses reads onto a smaller number of rows.

    The input dataframe is expected to be sorted in a sensible fashion and pre-indexed with a set of
    unique index values. This method has been tailored for and verified using read_start pre-sorting.
    Behavior using other starting sorts may be undefined.

    Optionally, performs a "meta-sort" of the resulting rows. Current options work as follows:
    * full_extent: sort by full extent of the covered reads in the row (max end - min start);
        rows covering a larger region at the bottom
    * covered_bases: sort by number of bases covered by reads in the row;
        rows covering more bases at the bottom
    * None: no meta-sorting

    Returns a series that maps the original indices on to the final collapsed and meta-sorted indices.
    This series can be applied to the original data by using the pd.Series.map() method.

    Args:
        read_extent_df: read extent dataframe from format_browser_data()
        minimum_gap: minimum number of bases allowed between the end of one read and the beginning of
            the next for the two reads to be placed on the same row
        meta_sort: type of meta sorting to do; one of ["full_extent", "covered_bases", None]

    Returns:
        Series mapping original indices to meta indices

    TODO: This could be improved by checking for overlaps on both ends of the seed read for each row.
        This might allow other types of pre-sorting to work more effectively.
    """
    read_count = len(read_extent_df)
    # Sentinel value for un-indexed reads is -1
    collapsed_indices = -np.ones(read_count, dtype=int)
    if read_count == 0:
        return pd.Series(collapsed_indices, index=read_extent_df["y_index"])

    starts = read_extent_df["read_start"].to_numpy(copy=False)
    ends = read_extent_df["read_end"].to_numpy(copy=False)

    # parent[i] points to the next unassigned index at or after i; parent[n] is a sentinel.
    parent = np.arange(read_count + 1, dtype=int)

    def _find_next_unassigned(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def _mark_assigned(index: int) -> None:
        parent[index] = _find_next_unassigned(index + 1)

    curr_y_idx = 0
    seed_search_start = 0
    while True:
        seed_read_idx = _find_next_unassigned(seed_search_start)
        if seed_read_idx >= read_count:
            break

        collapsed_indices[seed_read_idx] = curr_y_idx
        _mark_assigned(seed_read_idx)

        curr_row_end = int(ends[seed_read_idx])
        scan_cursor = seed_read_idx
        while True:
            # Preserve the original one-pass scan semantics while skipping impossible ranges.
            start_threshold = int(
                np.searchsorted(starts, curr_row_end + minimum_gap, side="right")
            )
            candidate_start = max(scan_cursor + 1, start_threshold)
            candidate_idx = _find_next_unassigned(candidate_start)
            if candidate_idx >= read_count:
                break

            collapsed_indices[candidate_idx] = curr_y_idx
            _mark_assigned(candidate_idx)
            curr_row_end = int(ends[candidate_idx])
            scan_cursor = candidate_idx

        curr_y_idx += 1
        seed_search_start = seed_read_idx

    # Series mapping original indices to collapsed indices
    idx_map_orig2collapse = pd.Series(
        collapsed_indices, index=read_extent_df["y_index"]
    )

    if meta_sort is not None:
        # Perform meta-sorting
        match meta_sort:
            case "full_extent":
                # sort by full extent of the covered reads in the row (max end - min start)
                grouped = (
                    read_extent_df.assign(_collapsed=collapsed_indices)
                    .groupby("_collapsed", sort=False)
                    .agg(
                        start_min=("read_start", "min"),
                        end_max=("read_end", "max"),
                    )
                )
                idx_map_meta2collapse = (
                    (grouped["end_max"] - grouped["start_min"])
                    .sort_values(ascending=False)
                    .index
                )
            case "covered_bases":
                # sort by number of bases covered by reads in the row
                read_lengths = read_extent_df["read_end"] - read_extent_df["read_start"]
                idx_map_meta2collapse = (
                    read_lengths.groupby(collapsed_indices)
                    .sum()
                    .sort_values(ascending=False)
                    .index
                )
            case _:
                raise ValueError(f"Invalid meta sorting option: {meta_sort}")

        # Series mapping collapsed indices to meta indices
        idx_map_collapse2meta = pd.Series(
            np.arange(len(idx_map_meta2collapse), dtype=int),
            index=idx_map_meta2collapse,
        )

        # Return Series mapping original indices to meta indices
        return idx_map_orig2collapse.map(idx_map_collapse2meta)
    else:
        # Return Series mapping original indices to collapsed indices
        return idx_map_orig2collapse


def make_browser_figure(
    read_extent_df: pd.DataFrame,
    mod_event_df: pd.DataFrame,
    collapse: bool,
    chrom: str,
    region_start: int,
    region_end: int,
    hover: bool = True,
    thresh: int | float | None = None,
    width: int = 1,
    size: int = 4,
    colorscales: dict = utils.DEFAULT_COLORSCALES,
    **kwargs,
) -> plotly.graph_objs.Figure:
    """
    Make a browser figure, using the provided pre-processed data

    Additional keyword arguments will be passed down to collapse_rows() if collapse == True. See that
    method for details.

    Args:
        read_extent_df: read extent dataframe from format_browser_data()
        mod_event_df: mod event dataframe from format_browser_data()
        collapse: if True, allows multiple reads on single rows of the browser for a more condensed
            visualization.
        chrom: chromosome of the region being browsed
        region_start: start position of the region being browsed
        region_end: end position of the region being browsed
        hover: if False, disables display of information on mouse hover
        thresh: pass down threshold for color scaling of modification probabilities
        width: width of the read lines in the browser
        size: size of the modification event markers in the browser
        colorscales: dictionary mapping motif names to plotly colorscale specifications

    TODO: Think about how this interfaces with different types of initial sorting...
    TODO: Should this method do the collapsing, or should this method require collapsing outside?
    """
    if collapse:
        index_map = collapse_rows(read_extent_df, **kwargs)
        # NOTE: `DataFrame.assign` performs a deep copy, preventing the modification of the array that is passed in. For very large datasets, this duplication may cause memory issues.
        read_extent_df = read_extent_df.assign(
            y_index=read_extent_df["y_index"].map(index_map)
        )
        mod_event_df = mod_event_df.assign(
            y_index=mod_event_df["y_index"].map(index_map)
        )

    # Build final figure
    # TODO: Enable setting some relevant parameters

    # TODO: Understand all of the options here; are they all as desired?
    layout = plotly.graph_objs.Layout(
        barmode="overlay",
        title=chrom,
        hovermode="closest",
        plot_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(range=[region_start, region_end]),
    )
    fig = plotly.graph_objects.Figure(layout=layout)
    if not read_extent_df.empty:
        starts = read_extent_df["read_start"].to_numpy(copy=False)
        ends = read_extent_df["read_end"].to_numpy(copy=False)
        y_values = read_extent_df["y_index"].to_numpy(copy=False)
        read_names = read_extent_df["read_name"].astype(str).to_numpy(copy=False)

        n_reads = len(read_extent_df)
        x_coords = np.empty(n_reads * 3, dtype=object)
        y_coords = np.empty(n_reads * 3, dtype=float)
        hover_text = np.empty(n_reads * 3, dtype=object)

        x_coords[0::3] = starts
        x_coords[1::3] = ends
        x_coords[2::3] = None

        y_coords[0::3] = y_values
        y_coords[1::3] = y_values
        y_coords[2::3] = np.nan

        hover_text[0::3] = read_names
        hover_text[1::3] = read_names
        hover_text[2::3] = ""

        fig.add_trace(
            plotly.graph_objects.Scatter(
                x=x_coords,
                y=y_coords,
                mode="lines",
                line=dict(width=width, color="lightgrey"),
                showlegend=False,
                hoverinfo="text",
                hovertext=hover_text,
                connectgaps=False,
            )
        )

    for motif_idx, (motif, motif_df) in enumerate(mod_event_df.groupby("motif")):
        min_overall = motif_df["prob"].min()
        max_overall = motif_df["prob"].max()
        fig.add_trace(
            plotly.graph_objs.Scatter(
                x=motif_df["pos"],
                y=motif_df["y_index"],
                mode="markers",
                showlegend=False,
                customdata=motif_df[["read_name", "prob"]],
                hovertemplate="<br>".join(
                    [
                        "<b>Read</b>: %{customdata[0]}",
                        "<b>Position</b>: %{x:,}",
                        "<b>Probability</b>: %{customdata[1]:.2f}",
                    ]
                ),
                marker=dict(
                    size=size,
                    color=motif_df["prob"],
                    colorscale=colorscales[motif],
                    cmin=thresh,
                    cmax=1.0,
                    colorbar=dict(
                        title=dict(
                            text=f"{motif} probability",
                            side="right",
                        ),
                        tickmode="array",
                        tickvals=[min_overall, max_overall],
                        ticktext=[
                            str(round(min_overall, 2)),
                            str(round(max_overall, 2)),
                        ],
                        ticks="outside",
                        thickness=15,
                        # TODO: Is this positioning system dumb?
                        x=1 + (motif_idx * 0.10),
                    ),
                ),
            )
        )
    if not hover:
        fig.update_layout(hovermode=False)
    return fig


def save_static(
    fig: plotly.graph_objs.Figure,
    output_dir: str | Path,
    output_basename: str,
    format: str | list[str] = "pdf",
    width: int = 1000,
    height: int = 400,
) -> None:
    """
    Helper function for saving static plot browser images.

    Args:
        fig: plotly figure to save
        output_dir: directory in which to save images
        output_basename: descriptive basename of output file (before file extension)
        format: one or more valid output formats for plotly; for valid options, see
            https://plotly.github.io/plotly.py-docs/generated/plotly.io.write_image.html
        width: width of output image in pixels
        height: height of output image in pixels
    """
    # sanitize and prep inputs
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    if isinstance(format, str):
        format = [format]

    # write figures
    for fmt in format:
        fig.write_image(
            output_dir / f"{output_basename}.{fmt}", width=width, height=height
        )


def save_interactive(
    fig: plotly.graph_objs.Figure,
    output_dir: str | Path,
    output_basename: str,
) -> None:
    """
    Helper function for saving interactive plot browsers.

    Args:
        fig: plotly figure to save
        output_dir: directory in which to save images
        output_basename: descriptive basename of output file (before file extension)
    """
    # sanitize inputs
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # write figure
    fig.write_html(output_dir / f"{output_basename}.html", include_plotlyjs="cdn")
