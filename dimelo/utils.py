import multiprocessing
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes

# This provides the mapping of canonical bases to sets of valid mode names
# It is a defaultdict because any bases without a default should still
# have valid entries, but they should be empty sets
BASEMOD_NAMES_DICT = defaultdict(lambda: set())
BASEMOD_NAMES_DICT.update(
    {
        "A": {"a", "Y"},
        "C": {"m", "Z"},
    }
)

# Default colors for seaborn plots
DEFAULT_COLORS = defaultdict(lambda: "grey")
DEFAULT_COLORS.update(
    {
        "A,0": "blue",
        "A,0,a": "blue",
        "CG,0": "orange",
        "CG,0,m": "yellow",
        "CG,0,h": "red",
        "GCH,1": "purple",
    }
)
# Default colorscales for plotly; based off of DEFAULT_COLORS
DEFAULT_COLORSCALES = defaultdict(lambda: ["white", "grey"])
DEFAULT_COLORSCALES.update([(k, ["white", v]) for k, v in DEFAULT_COLORS.items()])


# Define the source of randomness for a variety of purposes throughout the package
# TODO: how best to initialize seed for random state, to allow for reproducibility?
rng = np.random.default_rng()


def cores_to_run(cores):
    cores_avail = multiprocessing.cpu_count()
    if cores is None or cores > cores_avail:
        return cores_avail
    else:
        return cores


class ParsedMotif:
    def __init__(self, motif_string):
        """
        ParsedMotif takes in a base modification specifier motif, e.g. CG,0 or CG,0,m,
        and parses it for easy use later.

        args:
            -motif_string: a specifier string containing sequence_context,mod_coord(,optional mod_code)
        """
        parts = motif_string.split(",")
        if len(parts) == 2:
            # If a mod code isn't specified, we use the default set
            self.motif_seq = parts[0]
            self.modified_pos = int(parts[1])
            if self.modified_pos >= len(self.motif_seq):
                raise ValueError(f"Motif {motif_string} has an out-of-range mod index.")
            self.modified_base = self.motif_seq[self.modified_pos]
            self.mod_codes = BASEMOD_NAMES_DICT[self.modified_base]
        elif len(parts) == 3:
            # If a mod code is specified, that will be the only one we look for
            self.motif_seq = parts[0]
            self.modified_pos = int(parts[1])
            if self.modified_pos >= len(self.motif_seq):
                raise ValueError(f"Motif {motif_string} has an out-of-range mod index.")
            self.modified_base = self.motif_seq[self.modified_pos]
            self.mod_codes = set(parts[2])
        else:
            # Motifs need both a sequence and an index, separated by a comma
            raise ValueError(
                f"Motif {motif_string} must have 2 or 3 comma-separated elements: sequence, index, and (optionally) mod code."
            )


def adjust_threshold(
    thresh,
    quiet=True,
):
    if thresh > 0:
        if thresh > 1:
            if not quiet:
                print(
                    f"Modification threshold of {thresh} assumed to be for range 0-255. {thresh}/255={thresh/255} will be sent to modkit."
                )
            thresh_scaled = thresh / 255
        else:
            if not quiet:
                print(
                    f"Modification threshold of {thresh} will be treated as coming from range 0-1."
                )
            thresh_scaled = thresh

        return thresh_scaled
    return thresh


def regions_dict_from_input(
    regions: str | Path | list[str | Path] | None = None,
    window_size: int | None = None,
) -> dict:
    """
    Create a dictionary defining every region provided in the regions input
        key: chromosome strings
        value: lists of (start,end,strand)

    TODO: Maybe this dict object should be a more codified type
    """
    # TODO: Why is this declared out here, and not within add_region_to_dict? To my eye, that method should just return the fully-loaded dict.
    # I don't think this approach works because add_region_to_dict can be called many times; the regions parameter can be a single bed path / string OR many in a list
    regions_dict: defaultdict[str, list] = defaultdict(list)

    if window_size is not None and window_size <= 0:
        raise ValueError(
            "Invalid window_size. To disable windowing, set window_size to None or do not pass a value (the default is None)."
        )

    if isinstance(regions, list):
        for region in regions:
            add_region_to_dict(region, window_size, regions_dict)
    else:
        add_region_to_dict(regions, window_size, regions_dict)
    for chrom in regions_dict:
        regions_dict[chrom].sort(key=lambda x: x[0])

    return regions_dict


def add_region_to_dict(
    region: str | Path,
    window_size: int | None,
    regions_dict: dict,
):
    # TODO: The flow of this is very confusing, creating mypy errors, and possibly creates actual errors.
    # mypy error: dimelo/utils.py:110: error: Item "str" of "str | Path" has no attribute "name"  [union-attr]
    # Basically, this method is confusing because the string can be a pathlike or a region string.
    # Find a different way to check whether the string is pathlike or a region string, coerce paths to Path objects, then clean up everything else.

    # Added None as a window_size option, it was already handled below so was always a valid input

    # We check whether the region is a path to a .bed file by seeing if, when coerced into a Path object, it has the suffix 'bed'
    if Path(region).suffix == ".bed":
        with open(region) as bed_regions:
            for line_index, line in enumerate(bed_regions):
                fields = line.split()
                if len(fields) > 2:
                    # Per the bed spec, the 6th column is strand
                    # https://genome.ucsc.edu/FAQ/FAQformat.html
                    if len(fields) > 5:
                        chrom, start, end, strand = (
                            fields[0],
                            int(fields[1]),
                            int(fields[2]),
                            fields[5],
                        )
                    # If strand isn't in the bed file, we set to . (neither/both)
                    else:
                        chrom, start, end, strand = (
                            fields[0],
                            int(fields[1]),
                            int(fields[2]),
                            ".",
                        )
                    if window_size is None:
                        regions_dict[chrom].append((start, end, strand))
                    else:
                        center_coord = (start + end) // 2
                        regions_dict[chrom].append(
                            (
                                center_coord - window_size,
                                center_coord + window_size,
                                strand,
                            )
                        )
                else:
                    raise ValueError(
                        f"Invalid bed format line {line_index} of {Path(region).name}"
                    )
    # If the region is a path but *not* to a bed file, that isn't valid
    elif isinstance(region, Path):
        raise ValueError(
            f"Path object {region} is not pointing to a .bed file. regions must be provided as paths to .bed files or as strings in the format chrX:XXX-XXX,strand."
        )
    # If the region is a string and doesn't convert to a path to a bed file, then it must be a region string else it cannot be parsed
    elif (
        isinstance(region, str)
        and len(region.split(":")) == 2
        and 2 <= len(region.split(":")[1].split("-")) <= 3
    ):
        chrom, (start, end, strand) = parse_region_string(
            region=region, window_size=window_size
        )
        regions_dict[chrom].append((start, end, strand))
    else:
        raise ValueError(
            f"Invalid regions {type(region)}: {region}. Please use the format chrX:XXX-XXX,strand."
        )


def parse_region_string(
    region: str,
    window_size: int | None,
) -> tuple[str, tuple[int, int, str]]:
    """
    Parse a region specification string into its component parts.

    Args:
        region: a region string of the format chrX:XXX-XXX or chrX:XXX-XXX,strand (+/-/.)
        window_size: if present, returns a window of this size around the center of the given region

    Returns:
        chromosome, (start_pos, end_pos, strand)
    """
    try:
        # region strings can be either chrX:XXX-XXX or chrX:XXX-XXX,strand (+/-/.)
        region_coords = region.split(",")
        # The default strand is ., which is neither strand
        strand = region_coords[1] if len(region_coords) > 1 else "."
        chrom, coords = region_coords[0].split(":")
        start, end = map(int, coords.split("-"))
        if window_size is None:
            return chrom, (start, end, strand)
        else:
            center_coord = (start + end) // 2
            return chrom, (
                center_coord - window_size,
                center_coord + window_size,
                strand,
            )
    except (ValueError, AttributeError) as err:
        raise ValueError(
            f"Invalid region string {region}. Region strings can be either chrX:XXX-XXX or chrX:XXX-XXX,strand (+/-/.)."
        ) from err


def bed_from_regions_dict(
    regions_dict: dict,
    save_bed_path: Path,
):
    with open(save_bed_path, "w") as processed_bed:
        for chrom, regions_list in regions_dict.items():
            for start, end, _ in regions_list:
                bed_line = (
                    "\t".join([chrom, str(start), str(end), ".", ".", "."]) + "\n"
                )
                processed_bed.write(bed_line)


def bedmethyl_to_bigwig(input_bedmethyl: str | Path, output_bigwig: str | Path):
    return 0


def sanitize_path_args(*args) -> tuple:
    """
    Coerce all given arguments to Path objects, leaving Nones as Nones.
    """
    return tuple(Path(f) if f is not None else f for f in args)


def check_len_equal(*args: list) -> bool:
    """
    Checks whether all provided lists are the same length.
    """
    return all(len(x) == len(args[0]) for x in args)


def bar_plot(categories: list[str], values: np.ndarray, y_label: str, **kwargs) -> Axes:
    """
    Utility for producing bar plots.

    Args:
        categories: parallel with values; bar labels
        values: parallel with categories: bar heights
        y_label: y-axis label
        kwargs: other keyword parameters passed through to seaborn.barplot

    Returns:
        Axes object containing the plot
    """
    axes = sns.barplot(x=categories, y=values, hue=categories, **kwargs)
    axes.set(ylabel=y_label)
    return axes


def line_plot(
    indep_vector: np.ndarray,
    indep_name: str,
    dep_vectors: list[np.ndarray],
    dep_names: list[str],
    y_label: str,
    **kwargs,
) -> Axes:
    """
    Utility for producing overlayed line plots for data vectors with the same x-axis values.

    Takes in one independent vector and arbitrarily many dependent vectors. Plots all dependent vectors on the same axes against the same dependent vector.
    All vectors must be of equal length.

    TODO: Right now, this always generates a legend with the title "variable". I could add a parameter to specify this (by passing the var_name argument to pd.DataFrame.melt), but then that percolates upwards to other methods. How to do this cleanly?

    Args:
        indep_vector: parallel with each entry in dep_vectors; independent variable values shared across each overlayed line
        indep_name: name of independent variable; set as x axis label
        dep_vectors: outer list parallel with dep_names; each inner vector parallel with indep_vector; dependent variable values for each overlayed line
        dep_names: parallel with dep_vectors; names of each overlayed line; set as legend entries
        y_label: y-axis label
        kwargs: other keyword parameters passed through to seaborn.lineplot

    Returns:
        Axes object containing the plot

    Raises:
        ValueError: raised if any vectors are of unequal length
    """
    # construct dict of {vector_name: vector}, including the x vector using dict union operations
    data_dict = {indep_name: indep_vector} | dict(zip(dep_names, dep_vectors))
    # construct long-form data table for plotting
    try:
        data_table = pd.DataFrame(data_dict).melt(
            id_vars=indep_name, value_name=y_label
        )
    except ValueError as e:
        raise ValueError(
            "All dependent and independent vectors must be the same length"
        ) from e
    # plot lines
    return sns.lineplot(
        data=data_table, x=indep_name, y=y_label, hue="variable", **kwargs
    )


def hist_plot(
    value_vectors: list[np.ndarray],
    value_names: list[str],
    x_label: str,
    y_label: str,
    integer_values: bool = False,
    **kwargs,
) -> Axes:
    """
    Utility for producing overlayed histogram plots for data vectors containing values with some distribution.

    Takes arbitrarily many counts vectors and plots on same histogram.

    Args:
        value_vectors: parallel with value_names; vectors of values to plot histograms of; each vector will be a separate overlayed histogram
        value_names: parallel with value_vectors; names of each overlayed histogram; set as legend entries
        x_label: name of distributed values; set as x axis label
        y_label: y-axis label
        integer_values: True if hist bins are only at integer values, meaning bins shouldn't be auto-determined
        kwargs: other keyword parameters passed through to seaborn.histplot

    Returns:
        Axes object containing the plot

    Raises:
        ValueError: raised if any vectors are of unequal length
    """
    # Flatten the vectors and assign corresponding labels
    data_dict = {
        x_label: np.concatenate(value_vectors),
        y_label: np.repeat(value_names, [len(vec) for vec in value_vectors]),
    }

    # Create DataFrame
    data_table = pd.DataFrame(data_dict)
    if integer_values:
        # Warn user that passed bins are being overwritten
        if "bins" in kwargs:
            print("Warning: bin settings overwritten by defaults")
        kwargs["bins"] = np.arange(
            data_table[x_label].min() - 0.5, data_table[x_label].max() + 1.5, 1
        )

    # plot histogram
    ax = sns.histplot(
        data=data_table,
        x=x_label,
        hue=y_label,
        multiple="dodge",
        **kwargs,
    )

    ax.set_ylabel(y_label)

    return ax


def smooth_rolling_mean(
    vector: np.ndarray[float], window: int, min_periods: int = 1
) -> np.ndarray:
    """
    Smooths the given vector, using rolling centered windows of the given size.
    See pandas rolling documentation for details; documentation for relevant arguments copied here.

    Note: Because this operation is always centered, min_periods only has an effect if it is less than half of window size.

    TODO: Is pandas the most efficient implementation for this?
    TODO: Is it reasonable for min_periods to be default 1? That makes some sense for plotting, but might make analysis misleading in the future, compared to defaulting to window size.

    Args:
        vector: the vector of values to smooth
        window: size of the moving window
        min_periods: minimum number of observations in window to output a value; otherwise, result is np.nan

    Returns:
        Vector of smoothed values
    """
    return (
        pd.Series(vector)
        .rolling(window=window, min_periods=min_periods, center=True)
        .mean()
        .values
    )


def random_sample(
    array: np.ndarray,
    n: int | None = None,
    frac: float | None = None,
    replace: bool = False,
    # shuffle: bool = True,
):
    """
    Utility method for generating a random sample of the elements in an array. Always defaults to sampling along the first axis.
    This means that for 2d arrays this method will return a subsample of the rows.

    Handles n/frac specification like pandas.DataFrame.sample(): https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.sample.html.

    Args:
        array: Array to sample from
        n: Number of elements to return; mutually exclusive with frac
        frac: Fraction of elements to return; mutually exclusive with n
        replace: When True, sample with replacement

    Return: the requested random subset of the given array

    NOTE: Outputs are not guaranteed to be in the same order as the original array. If ordering is required to be preserved, re-sort after.

    TODO: enable non-uniform weights? Seems like too much to me...
    TODO: shuffled outputs originally broke h5 loading, so it was turned off. However, because of the instability of the sampling order,
        this doesn't actually matter. However, it is being left off to prevent confusion for the end user.
    """
    size = pd.core.sample.process_sampling_size(n, frac, replace)
    if size is None:
        assert frac is not None
        size = round(frac * len(array))
    return rng.choice(
        a=array,
        size=size,
        replace=replace,
        p=None,
        axis=0,
        # shuffle=shuffle
        shuffle=False,
    )
