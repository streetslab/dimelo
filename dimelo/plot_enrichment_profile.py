from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.axes import Axes

from . import load_processed, plotting, utils


def _looks_like_motif_spec(label: str) -> bool:
    parts = str(label).split(",")
    if len(parts) < 2:
        return False
    return bool(parts[0]) and parts[1].strip().isdigit()


def _legacy_aggregate_axis_spec(
    *,
    window_size: int,
    relative: bool,
    regions_5to3prime: bool,
) -> plotting.AxisSpec:
    return plotting.AxisSpec(
        orientation="region_5to3" if regions_5to3prime else "genomic",
        coordinate_mode="fixed_window",
        anchor="center" if relative else "absolute",
        upstream_bp=window_size,
        downstream_bp=window_size,
    )


def _legacy_aggregate_profile_table(
    trace_vectors: list[np.ndarray],
    sample_names: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for sample_name, trace_vector in zip(sample_names, trace_vectors, strict=True):
        positions = np.arange(
            -(len(trace_vector) // 2),
            len(trace_vector) // 2 + len(trace_vector) % 2,
        )
        for position, value in zip(positions, trace_vector, strict=True):
            rows.append(
                {
                    "sample_name": sample_name,
                    "position": float(position),
                    "anchor": 0.0,
                    "region_strand": "+",
                    "value": float(value),
                }
            )

    return pd.DataFrame(rows)


def _prepare_legacy_aggregate_profile_vectors(
    prepared_table: pd.DataFrame,
    *,
    sample_names: list[str],
) -> list[np.ndarray]:
    prepared_vectors: list[np.ndarray] = []
    for sample_name in sample_names:
        sample_table = prepared_table.loc[
            prepared_table["sample_name"] == sample_name
        ].copy()
        sample_table = sample_table.sort_values("plot_x")
        prepared_vectors.append(sample_table["value"].to_numpy())
    return prepared_vectors


def _route_legacy_aggregate_axis_through_shared_core(
    *,
    window_size: int,
    relative: bool,
    regions_5to3prime: bool,
    trace_vectors: list[np.ndarray],
    sample_names: list[str],
) -> list[np.ndarray]:
    axis = _legacy_aggregate_axis_spec(
        window_size=window_size,
        relative=relative,
        regions_5to3prime=regions_5to3prime,
    )
    aggregation = plotting.AggregationSpec()
    aggregate_table = _legacy_aggregate_profile_table(
        trace_vectors,
        sample_names,
    )
    prepared_payload = plotting.prepare_aggregate_plot_data(
        aggregate_table,
        plot_family="aggregate_profile",
        axis=axis,
        aggregation=aggregation,
        value_column="value",
        position_column="position",
        anchor_column="anchor",
        region_strand_column="region_strand",
    )
    return _prepare_legacy_aggregate_profile_vectors(
        prepared_payload["plot_table"],
        sample_names=sample_names,
    )


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

    Args:
        mod_file_names: list of paths to modified base data files
        bed_file_names: list of paths to bed files specifying centered equal-length regions
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files
        sample_names: list of names to use for labeling traces in the output; legend entries
        window_size: half-size of the desired window to plot; how far the window stretches on either side of the center point
        relative: True means x-axis is centered around region centers, False means x-axis is absolute genome positions. Must be True when plotting more than one region.
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

    trace_vectors = _route_legacy_aggregate_axis_through_shared_core(
        window_size=window_size,
        relative=relative,
        regions_5to3prime=regions_5to3prime,
        trace_vectors=trace_vectors,
        sample_names=sample_names,
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


def by_regions(
    mod_file_name: str | Path,
    regions_list: list[str | Path | list[str | Path]] | None = None,
    motif: str | None = None,
    regions: list[str | Path | list[str | Path]] | None = None,
    sample_names: list[str] | None = None,
    **kwargs,
) -> Axes:
    """
    Plot enrichment profile, holding modification file and modification types constant, varying regions

    Note: Sample names default to the names of the bed files.

    See plot_enrichment_profile for details.
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
    return plot_enrichment_profile(
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
    Plot enrichment profile, holding modification types and regions constant, varying modification files

    Note: Sample names default to the names of the modification files.

    See plot_enrichment_profile for details.
    """
    sample_names_for_plot = (
        sample_names
        if sample_names is not None
        else [str(mod_file) for mod_file in mod_file_names]
    )
    n_mod_files = len(mod_file_names)
    return plot_enrichment_profile(
        mod_file_names=mod_file_names,
        regions_list=[regions] * n_mod_files,
        motifs=[motif] * n_mod_files,
        sample_names=sample_names_for_plot,
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
    mod_file_paths = [Path(fn) for fn in mod_file_names]

    trace_vectors = []
    for mod_file, regions, motif in zip(
        mod_file_paths, regions_list, motifs, strict=False
    ):
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
    legend_title = kwargs.pop("legend_title", None)
    if (
        legend_title is None
        and sample_names
        and all(_looks_like_motif_spec(name) for name in sample_names)
    ):
        legend_title = "Modifications (motif, mod_index)"
    resolved_legend_title = "variable" if legend_title is None else str(legend_title)
    axes = utils.line_plot(
        indep_vector=np.arange(
            offset_center - len(trace_vectors[0]) // 2,
            offset_center + len(trace_vectors[0]) // 2 + len(trace_vectors[0]) % 2,
        ),
        indep_name="pos",
        dep_vectors=trace_vectors,
        dep_names=sample_names,
        y_label="fraction modified bases",
        legend_title=resolved_legend_title,
        **kwargs,
    )
    return axes
