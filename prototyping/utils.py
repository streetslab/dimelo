import numpy as np
import pandas as pd
from matplotlib.axes import Axes
import seaborn as sns
from pathlib import Path
import h5py
from collections import defaultdict

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

def read_by_base_txt_to_hdf5(
    input_txt: str | Path,
    output_h5: str | Path,
    basemod: str,
    thresh: float,
):
    motif,modco = tuple(basemod.split(','))
    motif_modified_base = motif[int(modco)]
    print(motif_modified_base)
    with h5py.File(output_h5,'w') as h5, open(input_txt) as txt:
        next(txt)
        read_name = ''
        for index,line in enumerate(txt):
            fields = line.split('\t')
            if read_name!=fields[0]:
                #New read
                #Record name
                read_name = fields[0]
                #Read in relevant values
                pos_in_read = int(fields[1])
                pos_in_genome = int(fields[2])
                read_chrom = fields[3]
                read_len = int(fields[9])
                canonical_base = fields[13]
                prob = float(fields[10])
                #Calculate read info
                read_start = pos_in_genome - pos_in_read
                read_end = read_start + int(fields[9])
                #Build read vectors
                mod_vector = np.zeros(read_end-read_start)
                val_vector = np.zeros(read_end-read_start)
                
                #Add modification to vector if type is correct
                if canonical_base == motif_modified_base:
                    print(pos_in_genome,read_start,len(val_vector))
                    val_vector[pos_in_genome-read_start] = 1
                    if prob>=thresh:
                        mod_vector[pos_in_genome-read_start] = 1
            else:
                pos_in_genome = int(fields[2])
                canonical_base = fields[13]
                prob = float(fields[10])
                if canonical_base == motif_modified_base:
                    print(pos_in_genome,read_start,len(val_vector))
                    val_vector[pos_in_genome-read_start] = 1
                    if prob>=thresh:
                        mod_vector[pos_in_genome-read_start] = 1
    print(read_name,np.sum(val_vector),np.sum(mod_vector))
    return 0

def check_len_equal(*args: list) -> bool:
    """
    Checks whether all provided lists are the same length
    """
    return all(len(x) == len(args[0]) for x in args)

def bar_plot(categories: list[str],
             values: np.ndarray,
             y_label: str) -> Axes:
    """
    Utility for producing bar plots; intended so that we can swap out plotting backends and styles easily.

    TODO: I'm not convinced we should actually be using seaborn here. It's kind of more heavy duty than we need? We're not really relying on pandas or data aggregation. But probably could refactor to match other methods.
    TODO: Set "pallete" appropriately

    Args:
        categories: parallel with values; bar labels
        values: parallel with categories: bar heights
        y_label: y-axis label
    
    Returns:
        Axes object containing the plot
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
    Utility for producing overlayed line plots with the same x-axis values; intended so that we can swap out plotting backends and styles easily.

    Takes in one independent vector and arbitrarily many dependent vectors. Overlays all dependent vectors on the same axes versus the same dependent vector.

    TODO: Don't like the names of the y parameters at all. They're value vectors?
    TODO: Color pallete
    TODO: What happens if this gets vectors of different lengths? What's the intended result?
    TODO: Right now, this always generates a legend with the title "variable". I could add a parameter to specify this (by passing the var_name argument to pd.DataFrame.melt), but then that percolates upwards to other methods and I'm lazy.

    Args:
        x: parallel with each entry in vectors; x values shared across each overlayed line
        x_label: x-axis label
        vectors: whole list parallel with vector_names; each entry parallel with x; y values for each overlayed line
        vector_names: parallel with vectors; names of each overlayed line; legend entries
        y_label: y-axis label

    Returns:
        Axes object containing the plot
    """
    # construct dict of {vector_name: vector}, including the x vector using dict union operations
    data_dict = {x_label: x} | dict(zip(vector_names, vectors))
    # construct long-form data table for plotting
    data_table = pd.DataFrame(data_dict).melt(id_vars=x_label, value_name=y_label)
    # plot lines
    return sns.lineplot(data=data_table, x=x_label, y=y_label, hue='variable')

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
