import numpy as np
import pandas as pd
from matplotlib.axes import Axes
import seaborn as sns

def check_len_equal(*args: list) -> bool:
    """
    Checks whether all provided lists are the same length
    """
    return all(len(x) == len(args[0]) for x in args)

def bar_plot(categories: np.ndarray,
             values: np.ndarray,
             y_label: str) -> Axes:
    """
    Utility for producing bar plots; intended so that we can swap out plotting backends and styles easily.

    TODO: I'm not convinced we should actually be using seaborn here. It's kind of more heavy duty than we need? We're not really relying on pandas or data aggregation. But probably could refactor to match other methods.
    TODO: Set "pallete" appropriately
    """
    axes = sns.barplot(x=categories, y=values, hue=categories)
    axes.set(ylabel=y_label)
    return axes

def line_plot(x: np.ndarray,
              x_label: str,
              vectors: list[np.ndarray],
              vector_names: list[str],
              y_label: str) -> Axes:
    """
    Utility for producing line plots; intended so that we can swap out plotting backends and styles easily.

    Takes in one independent vector and arbitrarily many dependent vectors. Overlays all dependent vectors on the same axes versus the same dependent vector.

    TODO: Don't like the names of the y parameters at all. They're value vectors?
    TODO: Color pallete
    TODO: What happens if this gets vectors of different lengths? What's the intended result?
    TODO: Right now, this always generates a legend with the title "variable". I could add a parameter to specify this (by passing the var_name argument to pd.DataFrame.melt), but then that percolates upwards to other methods and I'm lazy.
    """
    # construct dict of {vector_name: vector}, including the x vector using dict union operations
    data_dict = {x_label: x} | dict(zip(vector_names, vectors))
    # construct long-form data table for plotting
    data_table = pd.DataFrame(data_dict).melt(id_vars=x_label, value_name=y_label)
    # plot lines
    return sns.lineplot(data=data_table, x=x_label, y=y_label, hue='variable')
