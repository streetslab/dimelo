from pathlib import Path
from collections import defaultdict

import pysam
import numpy as np
import h5py

from . import test_data
from . import utils

def pileup_counts_from_bedmethyl(bedmethyl_file: Path,
                          mod_name: str,
                          bed_file: Path = None,
                          region_str: str = None,
                          ) -> tuple[int, int]:
    """
    Extract number of modified bases and total number of bases from the given bedmethyl file

    TODO: Merge bed_file / region_str handling into a unified function somewhere
    TODO: How to name this method?
    
    Args:
        bedmethyl_file: Path to bedmethyl file
        bed_file: Path to bed file specifying regions
        mod_name: type of modification to extract data for
    
    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """
    # Deleted the windowing stuff from the plot_enrichment_profile

    source_tabix = pysam.TabixFile(str(bedmethyl_file))
    # Don't need vectors, just need counts; also not guaranteed that windows are the same length?
    # valid_base_counts = np.zeros(window_size*2)
    # modified_base_counts = np.zeros(window_size*2)
    valid_base_count = 0
    modified_base_count = 0
    
    mod_motif = mod_name.split(',')[0]
    mod_coord_in_motif = mod_name.split(',')[1]

    if bed_file is not None and region_str is None:
        with open(bed_file) as regions_file:
            for line in regions_file:
                # pythonic condensed operation for extracting important fields
                # fields = line.split('\t')
                chromosome, start_pos, end_pos, *_ = line.split('\t')
                start_pos = int(start_pos)
                end_pos = int(end_pos)
                # Removing centering
                # center_coord = (int(fields[2])+int(fields[1]))//2
                # chromosome = fields[0]
                if chromosome in source_tabix.contigs:
                    # Removing centering
                    # if center_coord-window_size>0:
                    # for row in source_tabix.fetch(chromosome,center_coord-window_size,center_coord+window_size):
                    for row in source_tabix.fetch(chromosome, start_pos, end_pos):
                        tabix_fields = row.split('\t')
                        pileup_basemod = tabix_fields[3]
                        if mod_motif in pileup_basemod and mod_coord_in_motif in pileup_basemod:
                            pileup_info = tabix_fields[9].split(' ')
                            # Removing centering
                            # pileup_coord_relative = int(tabix_fields[1])-center_coord+window_size
                            # But actually don't need this because don't care about positions
                            # pileup_coord_relative = int(tabix_fields[1]) - end_pos
                            # valid_base_counts[pileup_coord_relative] += int(pileup_info[0])
                            # modified_base_counts[pileup_coord_relative] += int(pileup_info[2])
                            valid_base_count += int(pileup_info[0])
                            modified_base_count += int(pileup_info[2])
    elif region_str is not None and bed_file is None:
        chromosome, coords = region_str.split(':')
        start_coord, end_coord = map(int, coords.split('-'))
        if chromosome in source_tabix.contigs:
            for row in source_tabix.fetch(chromosome,start_coord,end_coord):
                tabix_fields = row.split('\t')
                pileup_basemod = tabix_fields[3]
                if mod_motif in pileup_basemod and mod_coord_in_motif in pileup_basemod:
                    pileup_info = tabix_fields[9].split(' ')
                    # Removing centering
                    # pileup_coord_relative = int(tabix_fields[1])-center_coord+window_size
                    # But actually don't need this because don't care about positions
                    # pileup_coord_relative = int(tabix_fields[1]) - end_pos
                    # valid_base_counts[pileup_coord_relative] += int(pileup_info[0])
                    # modified_base_counts[pileup_coord_relative] += int(pileup_info[2])
                    valid_base_count += int(pileup_info[0])
                    modified_base_count += int(pileup_info[2])  
    else:
        raise ValueError("Cannot process both a region string and a regions bed file.")

    return (modified_base_count, valid_base_count)

def counts_from_fake(*args,
                     **kwargs) -> tuple[int, int]:
    """
    Generates a fake set of enrichment counts. Ignores all arguments.

    Returns:
        tuple containing counts of (modified_bases, total_bases)
    """
    window_halfsize = 500
    return test_data.fake_peak_enrichment(halfsize=window_halfsize, peak_height=0.15)

def pileup_vectors_from_bedmethyl(bedmethyl_file: str | Path,
                          mod_name: str,
                          bed_file: str | Path = None,
                          region_str: str = None,
                          window_size: int = None) -> np.ndarray:
    """
    Generate trace for the specified modification aggregated across all regions in the given bed file. Called by profile plotters, can also be used by a user directly.

    TODO: Merge bed_file / region_str / window_size handling into a unified function somewhere
    TODO: How to name this method?
    TODO: I feel like stuff like this should be shared functionality
    TODO: I _THINK_ this should return a value for every position, with non-modified positions as zeros
    TODO: Currently, this redundantly generates a centered bed file and then never uses it; still doing centering in the actual parsing loop.
        Obvious solution is to just remove the second centering operation and reference the windowed file. But is there a smarter way to structure the code?

    Args:
        bedmethyl_file: Path to bedmethyl file
        bed_file: Path to bed file specifying centered equal-length regions
        mod_name: type of modification to extract data for
        window_size: the extent in either direction for windows around the center of regions.
    
    Returns:
        vector of fraction modifiied bases (e.g. mA/A) calculated for each position; float values between 0 and 1
    """
    
    source_tabix = pysam.TabixFile(str(bedmethyl_file))
    
    mod_motif = mod_name.split(',')[0]
    mod_coord_in_motif = mod_name.split(',')[1]
    
    if bed_file is not None and region_str is None:
        if window_size is None:
            bed_filepath_processed = Path(bed_file)
        elif window_size > 0:
            bed_filepath = Path(bed_file)
            print(f'Loading regions from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
            bed_filepath_processed = bedmethyl_file.parent / (bed_filepath.stem + f'.windowed{window_size}-for-readout' + bed_filepath.suffix)
            print(f'Writing new bed file {bed_filepath_processed.name}')
            utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
        else:
            raise ValueError(f'Invalid window size {window_size}.')
        with open(bed_filepath_processed) as regions_file:
            # Allocate the pileup arrays
            for line in regions_file:
                fields = line.split('\t')
                if window_size is None:
                    region_len = abs(int(fields[2])-int(fields[1])) 
                    valid_base_counts = np.zeros(region_len, dtype=int)
                    modified_base_counts = np.zeros(region_len, dtype=int)    
                else:
                    valid_base_counts = np.zeros(window_size*2, dtype=int)
                    modified_base_counts = np.zeros(window_size*2, dtype=int)   
                break
            regions_file.seek(0)
            for line_index,line in enumerate(regions_file):
                fields = line.split('\t')
                start_coord = min(int(fields[1]),int(fields[2]))
                end_coord = max(int(fields[1]),int(fields[2]))
                center_coord = (start_coord+end_coord)//2
                chromosome = fields[0]
                if chromosome in source_tabix.contigs:
                    # TODO: Does this mean that we throw out any windows that are too short on specifically the left side?
                    if window_size is None or center_coord-window_size>0:
                        if window_size is None:
                            search_start = start_coord
                            search_end = end_coord
                        else:
                            search_start = center_coord-window_size
                            search_end = center_coord+window_size
                        for row in source_tabix.fetch(chromosome,search_start,search_end):
                            tabix_fields = row.split('\t')
    #                         print(tabix_fields)
                            pileup_basemod = tabix_fields[3]
                            if mod_motif in pileup_basemod and mod_coord_in_motif in pileup_basemod:
                                pileup_info = tabix_fields[9].split(' ')
                                pileup_coord_relative = int(tabix_fields[1])-search_start
                                if window_size is None and pileup_coord_relative>region_len:
                                    print(f'WARNING: {bed_filepath_processed} line {line_index+1} (1-based) specifies a region that is longer than the first region; the end of the region will be skipped. To make a profile plot with differently-sized region, consider using the window_size parameter to make a profile across centered windows.')
                                    break
                                else:
                                    valid_base_counts[pileup_coord_relative] += int(pileup_info[0])
                                    modified_base_counts[pileup_coord_relative] += int(pileup_info[2])
    elif region_str is not None and bed_file is None:
        if window_size is not None:
            print('Window size is ignored if there is no region_bed')
        chromosome, coords = region_str.split(':')
        start_coord, end_coord = map(int, coords.split('-'))
        valid_base_counts = np.zeros(end_coord-start_coord,dtype=int)
        modified_base_counts = np.zeros(end_coord-start_coord,dtype=int)
        if chromosome in source_tabix.contigs:
            for row in source_tabix.fetch(chromosome,start_coord,end_coord):
                tabix_fields = row.split('\t')
                pileup_basemod = tabix_fields[3]
                if mod_motif in pileup_basemod and mod_coord_in_motif in pileup_basemod:
                    pileup_info = tabix_fields[9].split(' ')       
                    pileup_coord_relative = int(tabix_fields[1])-start_coord
                    valid_base_counts[pileup_coord_relative] += int(pileup_info[0])
                    modified_base_counts[pileup_coord_relative] += int(pileup_info[2])
    else:
        raise ValueError("Cannot processed both a bed file and a region string")
    
    return modified_base_counts,valid_base_counts

def vector_from_fake(window_size: int,
                     *args,
                     **kwargs) -> np.ndarray:
    """
    Generates a fake peak trace. Ignores all arguments except window_size.

    Args:
        window_size: halfsize of the window; how far the window stretches on either side of the center point
    
    Returns:
        vector of fraction modified bases calculated for each position; float values between 0 and 1
    """
    return test_data.fake_peak_enrichment_profile(halfsize=window_size, peak_height=0.15)

def convert_bytes(item):
    """Convert bytes to string if item is bytes, otherwise return as is."""
    if isinstance(item, bytes):
        return item.decode()
    return item

def convert_tuple_elements(tup):
    """Convert all bytes elements in a tuple to strings."""
    return tuple(convert_bytes(item) for item in tup)

def read_vectors_from_hdf5(
        file: str | Path,
        motifs: list[str],
        regions: str | Path | list[str | Path] = None,
        window_size: int = None,
        sort_by: str | list[str] = ['chromosome','region_start','read_start'],
) -> (list[tuple],list[str],dict):
    regions_dict = utils.regions_dict_from_input(
        regions=regions,
        window_size=window_size,
    )

    with h5py.File(file,'r') as h5:
        datasets = [name for name, obj in h5.items() if isinstance(obj, h5py.Dataset)]
        
        read_chromosomes = np.array(h5['chromosome'],dtype=str)
        read_starts = np.array(h5['read_start'])
        read_ends = np.array(h5['read_end'])
        read_motifs = np.array(h5['motif'],dtype=str)         

        
        if len(regions_dict)>0:
            read_data_list = []
            for chrom,region_list in regions_dict.items():
                for region_start,region_end in region_list:
                        relevant_read_indices = np.flatnonzero(
                            (read_ends > region_start) & 
                            (read_starts < region_end) & 
                            np.isin(read_motifs, motifs) & 
                            (read_chromosomes == chrom)
                        )               
                        read_data_list += list(zip(
                            *(h5[dataset][relevant_read_indices] for dataset in datasets),
                            [region_start for _ in relevant_read_indices],
                            [region_end for _ in relevant_read_indices]
                        ))
        else:
            relevant_read_indices = np.flatnonzero(
                np.isin(read_motifs, motifs)
            )
            read_data_list = list(zip(
                *(h5[dataset][relevant_read_indices] for dataset in datasets),
                [region_start for _ in relevant_read_indices],
                [region_end for _ in relevant_read_indices]
            ))
 
    # We add region information (start and end; chromosome is already present!) so that it is possible to sort by these
    datasets += ['region_start','region_end']
    try:
        sort_by_indices = [datasets.index(sort_item) for sort_item in sort_by]
    except ValueError as e:
        raise ValueError(f"Sorting error. {e}. Datasets include {datasets}")
    
    sorted_read_data = sorted(
        read_data_list, 
        key=lambda x: tuple(x[index] for index in sort_by_indices)
        )

    sorted_read_data_converted = [convert_tuple_elements(tup) for tup in sorted_read_data]

    return sorted_read_data_converted, datasets, regions_dict


def readwise_binary_modification_arrays(
    file: str | Path,
    motifs: list[str],
    regions: str | Path | list[str|Path],
    window_size: int = None,
    sort_by: str | list[str] = ['chromosome','read_start'],
    thresh: float = None,
    relative: bool = True,
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
    file = Path(file)
    if file.suffix=='.h5' or file.suffix=='.hdf5':
        sorted_read_data_converted, datasets, regions_dict = read_vectors_from_hdf5(
            file = file,
            motifs = motifs,
            regions = regions,
            window_size = window_size,
            sort_by = sort_by,
        )
        read_name_index = datasets.index('read_name')
        mod_vector_index = datasets.index('mod_vector')
        motif_index = datasets.index('motif')
        region_start_index = datasets.index('region_start')
        region_end_index = datasets.index('region_end')
        read_start_index = datasets.index('read_start')

        # Check that this .h5 file was created with a threshold, i.e. that the mod calls are binarized
        if thresh is None:
            if not (sorted_read_data_converted[0][mod_vector_index].dtype==np.bool_):
                raise ValueError('No threshold has been applied to this .h5 single read data. You must provide a threshold using the thresh parameter in order to extract binarized modification arrays.')
        else:
            thresh = utils.adjust_threshold(thresh)

        read_ints_list = []
        mod_coords_list = []
        motifs_list = []

        read_names = np.array(
            [read_data[read_name_index] 
            for read_data in sorted_read_data_converted]
            )
        # TODO: handle the case where a read shows up in more than one different region
        unique_read_names, first_indices = np.unique(read_names,return_index=True)
        string_to_int = {read_name: index for index, read_name in zip(first_indices,unique_read_names)}
        read_ints = np.array([string_to_int[read_name] for read_name in read_names])
        for read_int,read_data in zip(read_ints,sorted_read_data_converted):
            if thresh is None:
                mod_pos_in_read = np.flatnonzero(read_data[mod_vector_index])
            else:
                mod_pos_in_read = np.flatnonzero(read_data[mod_vector_index]>thresh)

            if relative:
                mod_pos_record = mod_pos_in_read + read_data[read_start_index] - (read_data[region_start_index]+read_data[region_end_index])//2
            else:
                mod_pos_record = mod_pos_in_read + read_data[read_start_index]

            mod_coords_list += list(mod_pos_record)
            read_ints_list += [read_int]*len(mod_pos_record)
            motifs_list += [read_data[motif_index]]*len(mod_pos_record)

        return (np.array(mod_coords_list),np.array(read_ints_list),np.array(motifs_list),regions_dict)

    else:
        raise ValueError(f'File {file} does not have a recognized extension for single read data.')
def reads_from_hdf5(
    file: str | Path,
    mod_names: list[str],
    bed_file: str | Path = None,
    region_str: str = None,
    window_size: int = None,
) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
    """
    Pulls a list of read data out of an hdf5 file containing processed read vectors.

    Args:
        file: Path to an hdf5 (.h5) file containing modification data for single reads,
            stored in datasets read_name, chromosome, read_start,
            read_end, base modification motif, mod_vector, and val_vector.
        bed_file: Path to bed file specifying regions from which to draw reads. These
            should all be regions for which your original .bam file had reads extracted,
            although by design this method will not raise an error if any region contains
            zero reads, as this may simply be a matter of low read depth.
        mod_names: types of modification to extract data for. Basemods are specified as 
            {sequence_motif},{position_of_modification}. For example, a methylated adenine is specified 
            as 'A,0' and CpG methylation is specified as 'CG,0'.
        window_size: 
    
    Returns:
        Returns three parallel arrays, of length (N_READS * len(mod_names)), containing the following for each index:
        * array of positions at which the specified modification was found in a read
        * unique integer ID for the read
        * modification represented by the positions
        For example, if called on a dataset with a single read and two modification types, each array would have two entries. The unique IDs would be the same, as both entries would represent the same single read. The mods and positions would be different, as they would extact different mods.
    """
    regions_dict = defaultdict(list)
    if bed_file is not None and region_str is None:
        if window_size is None:
            bed_filepath = Path(bed_file)
            print(f'Using window size defined by bed file {bed_filepath.name}.')
            bed_filepath_processed = bed_filepath
        elif window_size>0:
            bed_filepath = Path(bed_file)
            print(f'Loading regions from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
            bed_filepath_processed = file.parent / (bed_filepath.stem + f'.windowed{window_size}-for-readout' + bed_filepath.suffix)
            print(f'Writing new bed file {bed_filepath_processed.name}')
            utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
        else:
            raise ValueError(f'Invalid window_size of {window_size}')
        with open(bed_filepath_processed) as bed_regions:
            for line in bed_regions:
                fields = line.split()
                if len(fields)>2:
                    chrom,start,end = fields[0],int(fields[1]),int(fields[2])
                    regions_dict[chrom].append((start,end))
    elif region_str is not None and bed_file is None:
        if window_size is not None:
            print('WARNING: window_size will be ignored, this parameter only works with bed file regions')
        chromosome, coords = region_str.split(':')
        start_coord, end_coord = map(int, coords.split('-'))
        regions_dict[chromosome].append((start_coord,end_coord))
    else:
        raise ValueError('Cannot use both a bed file and a region string')        
    for chrom in regions_dict:
        regions_dict[chrom].sort(key=lambda x: x[0])
    
    mod_coords_list = []
    read_ints_list = []
    mod_names_list = []
    
    with h5py.File(file,'r') as h5:
        read_names = np.array(h5['read_name'],dtype=str)
        unique_read_names, first_indices = np.unique(read_names,return_index=True)
        string_to_int = {read_name: index for index, read_name in zip(first_indices,unique_read_names)}
        read_ints = np.array([string_to_int[read_name] for read_name in read_names])
        read_chromosomes = np.array(h5['chromosome'],dtype=str)
        read_starts = np.array(h5['read_start'])
        read_ends = np.array(h5['read_end'])
        read_motifs = np.array(h5['motif'],dtype=str)   

        # Go through the regions one by one
        for chrom,region_list in regions_dict.items():
            for mod_name in mod_names:
                for region_start,region_end in region_list:
                    # Finds reads that are within the window you want
                    center_coord = (region_start+region_end)//2
                    relevant_read_indices = np.flatnonzero(
                        (read_ends > region_start) & 
                        (read_starts < region_end) & 
                        (read_motifs == mod_name) & 
                        (read_chromosomes == chrom)
                    )
                    # For each read, adds all of its elements to the lists that will be getting returned
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
                
    return (np.array(mod_coords_list),np.array(read_ints_list),np.array(mod_names_list))

def reads_from_fake(mod_names: list[str],
                    *args,
                    **kwargs) -> tuple[list[np.ndarray], np.ndarray[int], np.ndarray[str]]:
    """
    Generates a set of reads with modifications, spanning a region of some specified size. Ignores all arguments except for mod_names.

    TODO: Is integer index still correct for this?

    Args:
        mod_names: types of modification to extract data for; see below for valid choices
    
    Returns:
        Returns three parallel arrays, of length (N_READS * len(mod_names)), containing the following for each index:
        * array of positions at which the specified modification was found in a read
        * unique integer ID for the read
        * modification represented by the positions
        For example, if called on a dataset with a single read and two modification types, each array would have two entries. The unique IDs would be the same, as both entries would represent the same single read. The mods and positions would be different, as they would extact different mods.
    """
    window_halfsize = 500
    n_reads = 500

    reads = []
    read_names = []
    mods = []
    for mod_name in mod_names:
        match mod_name:
            case 'A':
                mod_reads = [test_data.fake_read_mod_positions(window_halfsize, 'peak', 0.7) for _ in range(n_reads)]
            case 'C':
                mod_reads = [test_data.fake_read_mod_positions(window_halfsize, 'inverse_peak', 0.4) for _ in range(n_reads)]
            case _:
                raise ValueError(f'No stub settings for requested mod {mod_name}')
        reads += mod_reads
        read_names.append(np.arange(len(mod_reads)))
        mods.append([mod_name] * len(mod_reads))
    
    read_names = np.concatenate(read_names)
    mods = np.concatenate(mods)
    return reads, read_names, mods