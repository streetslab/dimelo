"""
I'm conflicted about how to handle some of this.

There are two different ways of doing single read plotting: "rectangular" and "whole read".
"rectangular" means displaying exactly the requested region.
"whole read" means displaying the entirety of any read overlapping the requested region.
Probably need separate methods for all of this? Is there shared functionality? Do they live in the same file? Etc.

I'm beginning to lose the thread of where we check for regions making sense.
Maybe this is an argument for an internal region class that makes checking easy? I don't know.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
from collections import defaultdict
import h5py

import utils
import test_data

""" TEMPORARY STUB VARS """
STUB_HALFSIZE = 500
STUB_N_READS = 500


def binary_search_region(
    regions,
    chrom,
    position,
):
    region_list = regions.get(chrom,[])
    left,right = 0, len(region_list)-1
    while left <= right:
        mid = left + (right-left)//2
        if region_list[mid][0] <= position <= region_list[mid][1]:
            return True
        elif position < region_list[mid][0]:
            right == mid - 1
        else:
            left = mid + 1
    return False

def extract_centered_reads_from_hdf5(
    file: Path,
    bed_file: Path,
    mod_names: list[str],
    window_size:int=0,
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
    """
    TODO: What does the bed file represent in this method? This one is breaking my brain a bit.
    TODO: Variable names in this method stink.
    TODO: Currently assumes mod calling (thresholding probabilities) was already performed elsewhere

    Args:
        file: Path to file containing modification data for single reads
        bed_file: Path to bed file specifying regions (WHAT DO THESE REPRESENT???)
        mod_names: types of modification to extract data for
    
    Returns:
        Returns three parallel arrays, of length (N_READS * len(mod_names)), containing the following for each index:
        * array of positions at which the specified modification was found in a read
        * unique integer ID for the read
        * modification represented by the positions
        For example, if called on a dataset with a single read and two modification types, each array would have two entries. The unique IDs would be the same, as both entries would represent the same single read. The mods and positions would be different, as they would extact different mods.
    """
    if window_size>0:
        bed_filepath = Path(bed_file)
        print(f'Loading regions from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
        bed_filepath_processed = file.parent / (bed_filepath.stem + f'.windowed{window_size}-for-readout' + bed_filepath.suffix)
        print(f'Writing new bed file {bed_filepath_processed.name}')
        utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
    else:
        bed_filepath = Path(bed_file)
        print(f'Using window size defined by bed file {bed_filepath.name}.')
        bed_filepath_processed = bed_filepath
    
    regions_dict = defaultdict(list)
    with open(bed_filepath_processed) as bed_regions:
        for line in bed_regions:
            fields = line.split()
            if len(fields)>2:
                chrom,start,end = fields[0],int(fields[1]),int(fields[2])
                regions_dict[chrom].append((start,end))
            
    for chrom in regions_dict:
        regions_dict[chrom].sort(key=lambda x: x[0])
    
    mod_coords_list = []
    read_ints_list = []
    mod_names_list = []
#     for mod_name in mod_names:
#         match mod_name:
#             case 'A,0':
#                 mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'peak', 0.7) for _ in range(STUB_N_READS)]
#             case 'CG,0':
#                 mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'inverse_peak', 0.4) for _ in range(STUB_N_READS)]
#             case _:
#                 raise ValueError(f'No stub settings for requested mod {mod_name}')
#         reads += mod_reads
#         read_names.append(np.arange(len(mod_reads)))
#         mods.append([mod_name] * len(mod_reads))
    
    in_regions = 0
    out_regions = 0
    with h5py.File(file,'r') as h5:
        read_names = np.array(h5['read_name'],dtype=str)
        unique_read_names, first_indices = np.unique(read_names,return_index=True)
        string_to_int = {read_name: index for index, read_name in zip(first_indices,unique_read_names)}
        read_ints = np.array([string_to_int[read_name] for read_name in read_names])
        read_chromosomes = np.array(h5['chromosome'],dtype=str)
        read_starts = np.array(h5['read_start'])
        read_ends = np.array(h5['read_end'])
        read_motifs = np.array(h5['motif'],dtype=str)   
#         print('starts',read_starts)
#         print('ends',read_ends)
        for chrom,region_list in regions_dict.items():
            for mod_name in mod_names:
                for region_start,region_end in region_list:
#                     print(chrom,region_start,region_end)
                    center_coord = (region_start+region_end)//2
                    relevant_read_indices = np.flatnonzero(
                        (read_ends > region_start) & 
                        (read_starts < region_end) & 
                        (read_motifs == mod_name) & 
                        (read_chromosomes == chrom)
                    )
#                     print(relevant_read_indices)
                    for read_index in relevant_read_indices:
                        mod_vector = np.array(h5['mod_vector'][read_index])
                        read_start = read_starts[read_index]
                        read_int = read_ints[read_index]
                        mod_indices = np.flatnonzero(mod_vector)
                        for mod_index in mod_indices:
                            mod_coord_absolute = mod_index+read_start
                            if region_start < mod_coord_absolute < region_end:
                                mod_coords_list.append(mod_coord_absolute-center_coord)
                                read_ints_list.append(read_int)
                                mod_names_list.append(mod_name)
                
#     print(in_regions,out_regions)
#     for mod_name in mod_names:
#         match mod_name:
#             case 'A':
#                 mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'peak') for _ in range(STUB_N_READS)]
#             case 'C':
#                 mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'inverse_peak') for _ in range(STUB_N_READS)]
#             case _:
#                 raise ValueError(f'No stub settings for requested mod {mod_name}')
#         reads += mod_reads
#         read_names.append(np.arange(len(mod_reads)))
#         mods.append([mod_name] * len(mod_reads))
    
#     read_names = np.concatenate(read_names)
#     mods = np.concatenate(mods)
    return (np.array(mod_coords_list),np.array(read_ints_list),np.array(mod_names_list))

def extract_centered_reads_from_UNKNOWN_FILE_TYPE(file: Path,
                                                  bed_file: Path,
                                                  mod_names: list[str]) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
    """
    TODO: What does the bed file represent in this method? This one is breaking my brain a bit.
    TODO: Variable names in this method stink.
    TODO: Currently assumes mod calling (thresholding probabilities) was already performed elsewhere

    Args:
        file: Path to file containing modification data for single reads
        bed_file: Path to bed file specifying regions (WHAT DO THESE REPRESENT???)
        mod_names: types of modification to extract data for
    
    Returns:
        Returns three parallel arrays, of length (N_READS * len(mod_names)), containing the following for each index:
        * array of positions at which the specified modification was found in a read
        * unique integer ID for the read
        * modification represented by the positions
        For example, if called on a dataset with a single read and two modification types, each array would have two entries. The unique IDs would be the same, as both entries would represent the same single read. The mods and positions would be different, as they would extact different mods.
    """
    reads = []
    read_names = []
    mods = []
    for mod_name in mod_names:
        match mod_name:
            case 'A':
                mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'peak', 0.7) for _ in range(STUB_N_READS)]
            case 'C':
                mod_reads = [test_data.fake_read_mod_positions(STUB_HALFSIZE, 'inverse_peak', 0.4) for _ in range(STUB_N_READS)]
            case _:
                raise ValueError(f'No stub settings for requested mod {mod_name}')
        reads += mod_reads
        read_names.append(np.arange(len(mod_reads)))
        mods.append([mod_name] * len(mod_reads))
    
    read_names = np.concatenate(read_names)
    mods = np.concatenate(mods)
    return reads, read_names, mods

def plot_single_reads_rectangle(mod_file_name: str | Path,
                                bed_file_name: str | Path,
                                mod_names: list[str],
                                window_size: int = 0) -> Axes:
    """
    Plots centered single reads as a scatterplot, cut off at the boundaries of the requested regions?

    TODO: Clarify this documentation it's a mess. How do I say this concisely?
    TODO: I feel like this should be able to take in data directly as vectors/other datatypes, not just read from files.
    TODO: Style-wise, is it cleaner to have it be a match statement or calling a method from a global dict? Cleaner here with a dict, cleaner overall with the match statements?
    TODO: This name stinks?
    TODO: So far, this is the only method to do plotting without utility methods. Is this reasonable? Is it that unique?

    Args:
        mod_file_name: path to file containing modification data for single reads
        bed_file_name: path to bed file specifying regions (WHAT DO THESE REPRESENT???)
        mod_names: list of modifications to extract; expected to match mods available in the relevant mod_files

    Returns:
        Axes object containing the plot
    """ 
    mod_file_name = Path(mod_file_name)
    bed_file_name = Path(bed_file_name)

    match mod_file_name.suffix:
        case _:
            reads, read_names, mods = extract_centered_reads_from_hdf5(
                file=mod_file_name,
                bed_file=bed_file_name,
                mod_names=mod_names,
                window_size=window_size,
            )

    # Convert data frame where each row represents a read to a data frame where each row represents a single modified position in a read
    df = pd.DataFrame({
        'read_name': read_names,
        'mod': mods,
        'pos': reads
    }).explode('pos')
    axes = sns.scatterplot(
        data=df,
        x="pos",
        y="read_name",
        hue="mod",
        # palette=colors,
        s=0.5,
        marker="s",
        linewidth=0,
    )
    # Retrieve the existing legend
    legend = axes.legend_

    # Update legend properties
    legend.set_title('Mod')
    for handle in legend.legendHandles:
        handle.set_markersize(10)  # Set a larger marker size for legend
#     handles,labels = axes.get_legend_handles_labels()
#     # Create a new legend with larger marker size
#     new_handles = [plt.Line2D([], [], 
#                               marker='s', 
#                               linestyle='', 
#                               color=handle.get_color(), 
#                               markersize=10) 
#                    for handle in handles] # Skip the first handle as it is the legend title
#     axes.legend(new_handles, labels[1:], title='Mod', loc='upper right')
    
    return axes