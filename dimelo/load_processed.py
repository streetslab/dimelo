from pathlib import Path
from collections import defaultdict

import pysam
import numpy as np
import h5py

from . import test_data
from . import utils

def counts_from_bedmethyl(bedmethyl_file: Path,
                          bed_file: Path,
                          mod_name: str) -> tuple[int, int]:
    """
    Extract number of modified bases and total number of bases from the given bedmethyl file
    
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

    return (modified_base_count, valid_base_count)

def counts_from_fake(mod_file: Path,
                     bed_file: Path,
                     mod_name: str) -> tuple[int, int]:
    """
    Generates a fake set of counts
    """
    window_halfsize = 500
    return test_data.fake_peak_counts(halfsize=window_halfsize, peak_height=0.15)

def vector_from_bedmethyl(bedmethyl_file: Path,
                          bed_file: Path,
                          mod_name: str,
                          window_size: int) -> np.ndarray:
    """
    Generate trace for the specified modification aggregated across all regions in the given bed file.

    Args:
        bedmethyl_file: Path to bedmethyl file
        bed_file: Path to bed file specifying centered equal-length regions
        mod_name: type of modification to extract data for
        window_size: the extent in either direction for windows around the center of regions.
    
    Returns:
        vector of fraction modifiied bases (e.g. mA/A) calculated for each position; float values between 0 and 1
    """
    if window_size > 0:
        # TODO: I think this should not need to be explicitly done here; the specification for this method is that it's already a Path object. If this is intended to be user facing, we can do this, but I think it's overkill.
        bed_filepath = Path(bed_file)
        print(f'Loading regions from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
        bed_filepath_processed = bedmethyl_file.parent / (bed_filepath.stem + f'.windowed{window_size}-for-readout' + bed_filepath.suffix)
        print(f'Writing new bed file {bed_filepath_processed.name}')
        utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
    else:
        print(f'Invalid window size {window_size}.')
        return -1
    
    source_tabix = pysam.TabixFile(str(bedmethyl_file))
    valid_base_counts = np.zeros(window_size*2, dtype=int)
    modified_base_counts = np.zeros(window_size*2, dtype=int)
    
    mod_motif = mod_name.split(',')[0]
    mod_coord_in_motif = mod_name.split(',')[1]
    
    with open(bed_file) as regions_file:
        for line in regions_file:
            fields = line.split('\t')
            center_coord = (int(fields[2])+int(fields[1]))//2
            chromosome = fields[0]
            if chromosome in source_tabix.contigs:
                # TODO: Does this mean that we throw out any windows that are too short on specifically the left side?
                if center_coord-window_size>0:
                    for row in source_tabix.fetch(chromosome,center_coord-window_size,center_coord+window_size):
                        tabix_fields = row.split('\t')
#                         print(tabix_fields)
                        pileup_basemod = tabix_fields[3]
                        if mod_motif in pileup_basemod and mod_coord_in_motif in pileup_basemod:
                            pileup_info = tabix_fields[9].split(' ')
                            pileup_coord_relative = int(tabix_fields[1])-center_coord+window_size
                            valid_base_counts[pileup_coord_relative] += int(pileup_info[0])
                            modified_base_counts[pileup_coord_relative] += int(pileup_info[2])
                        
    modified_fractions = np.divide(modified_base_counts,valid_base_counts, out=np.zeros_like(modified_base_counts, dtype=float), where=valid_base_counts!=0)
    return modified_fractions

def vector_from_fake(mod_file: Path,
                     bed_file: Path,
                     mod_name: str,
                     window_size: int) -> np.ndarray:
    """
    Generates a fake peak trace.
    """
    return test_data.fake_peak_trace(halfsize=window_size, peak_height=0.15)

""" TEMPORARY STUB VARS """
STUB_HALFSIZE = 500
STUB_N_READS = 500

def reads_from_hdf5(
    file: Path,
    bed_file: Path,
    mod_names: list[str],
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

def reads_from_fake(file: Path,
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