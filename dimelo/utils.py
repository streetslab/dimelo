import numpy as np
import pandas as pd
from matplotlib.axes import Axes
import seaborn as sns
from pathlib import Path

def generate_centered_windows_bed(
    input_bed: str | Path,
    output_bed: str | Path,
    window_size: int,
):
    """
    TODO: Documentation; I think window_size is a half-size?
    TODO: Do we anticipate always working from bed files, or should the centering functionality work on some internal region representation; would just need a bed file wrapper around it, then.
    TODO: I think I like the way this works, where you can pass it a bed file with arbitrary regions and it will find the center --> generate windows of the appropriate size. But then again maybe we should discuss.
    TODO: Right now, this returns even-number-length windows on the half-open interval [center_coord - window_size, center_coord + window_size). I think it should do odd-number-length windows on the closed interval. Thoughts?
    """
    with open(output_bed,'w') as windowed_bed:
        with open(input_bed) as source_bed:
            for line in source_bed:
                fields = line.split('\t')
                windowed_fields = fields
                center_coord = (int(fields[1])+int(fields[2]))//2
                windowed_fields[1] = str(center_coord - window_size)
                windowed_fields[2] = str(center_coord + window_size)
                windowed_line = '\t'.join(windowed_fields)
                windowed_bed.write(windowed_line)

def bedmethyl_to_bigwig(
    input_bedmethyl: str | Path,
    output_bigwig: str | Path
):
    return 0

def check_len_equal(*args: list) -> bool:
    """
    Checks whether all provided lists are the same length.
    """
    return all(len(x) == len(args[0]) for x in args)

def bar_plot(categories: list[str],
             values: np.ndarray,
             y_label: str,
             **kwargs) -> Axes:
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

def line_plot(indep_vector: np.ndarray,
              indep_name: str,
              dep_vectors: list[np.ndarray],
              dep_names: list[str],
              y_label: str,
              **kwargs) -> Axes:
    """
    Utility for producing overlayed line plots for data vectors with the same x-axis values.

    Takes in one independent vector and arbitrarily many dependent vectors. Plots all dependent vectors on the same axes against the same dependent vector.
    All vectors must be of equal length.

    TODO: Right now, this always generates a legend with the title "variable". I could add a parameter to specify this (by passing the var_name argument to pd.DataFrame.melt), but then that percolates upwards to other methods. How to do this cleanly?

    Args:
        indep_vector: parallel with each entry in vectors; independent variable values shared across each overlayed line
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
        data_table = pd.DataFrame(data_dict).melt(id_vars=indep_name, value_name=y_label)
    except ValueError:
        raise ValueError('All dependent and independent vectors must be the same length')
    # plot lines
    return sns.lineplot(data=data_table, x=indep_name, y=y_label, hue='variable', **kwargs)

def smooth_rolling_mean(vector: np.ndarray[float],
                        window: int,
                        min_periods: int = 1) -> np.ndarray:
    """
    Smooths the given vector, using rolling centered windows of the given size.

    See pandas rolling documentation for details; documentation for relevant arguments copied here.

    TODO: It feels a little silly to be delegating smoothing to pandas, but this works. Consider refactoring?
    TODO: Is it reasonable for min_periods to be default 1? That makes some sense for plotting, but might make analysis misleading in the future.

    Args:
        vector: the vector of values to smooth
        window: size of the moving window
        min_periods: minimum number of observations in window to output a value; otherwise, result is np.nan
    """
    return pd.Series(vector).rolling(window=window, min_periods=min_periods, center=True).mean().values
