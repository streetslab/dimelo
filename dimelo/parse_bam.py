import sys
import os
import shutil
import subprocess
import multiprocessing
from pathlib import Path

import numpy as np
import h5py
import pysam

from . import utils

"""
This module contains code to convert .bam files into both human-readable and 
indexed random-access pileup and read-wise processed outputs.
"""

"""
Global variables
"""

# This should be updated in tandem with the environment.yml nanoporetech::modkit version
EXPECTED_MODKIT_VERSION = '0.2.4'

# Specifies how many reads to check for the base modifications of interest.
NUM_READS_TO_CHECK = 500


"""
Import checks
"""
# Add conda env bin folder to path if it is not already present
current_interpreter = sys.executable
env_bin_path = os.path.dirname(current_interpreter)
if env_bin_path not in os.environ['PATH']:
    print(f'PATH does not include the conda environment /bin folder. Adding {env_bin_path}.')
    os.environ['PATH'] = f'{env_bin_path}:{os.environ["PATH"]}'
    print(f'PATH is now {os.environ["PATH"]}')

# Check modkit on first import
try:
    result = subprocess.run(['modkit','--version'], stdout=subprocess.PIPE, text=True)
    modkit_version = result.stdout
    if modkit_version.split()[1] == EXPECTED_MODKIT_VERSION:
        print(f'modkit found with expected version {EXPECTED_MODKIT_VERSION}')
    else:
        print(f'modkit found with unexpected version {modkit_version.split()[1]}. Versions other than {EXPECTED_MODKIT_VERSION} may exhibit unexpected behavior. It is recommended that you use v{EXPECTED_MODKIT_VERSION}')
except:
    print('Executable not found for modkit. Install dimelo using "conda env create -f environment.yml" or install modkit manually to your conda environment using "conda install nanoporetech::modkit==0.2.4". Without modkit you cannot run parse_bam functions.')

"""
User-facing parse operations: pileup and extract
"""
def pileup(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    regions: str | Path | list[str | Path] = None,
    motifs: list = ['A,0','CG,0'],
    thresh: float = None,
    window_size: int = None,
    cores: int = None,
    log: bool = False,
    cleanup: bool = True,) -> Path:

    """
    TODO: Merge bed_file / region_str / window_size handling into a unified function somewhere

    Takes a file containing long read sequencing data aligned 
    to a reference genome with modification calls for one or more base/context 
    and creates a pileup. A pileup is a genome-position-wise sum of both reads with
    bases that could have the modification in question and of reads that are in
    fact modified.

    The current implementation of this method uses modkit, a tool built by 
    Nanopore Technologies, along with htslib tools compress and index the output
    bedmethyl file.

    https://github.com/nanoporetech/modkit/

    Args:
        output_file: a string or Path object pointing to the location of a .bam file.
            The file should follow at least v1.6 of the .bam file specifications, 
            found here: https://samtools.github.io/hts-specs/
            https://samtools.github.io/hts-specs/SAMv1.pdf 

            The file needs to have modifications stored in the standard format,
            with MM and ML tags (NOT mm and ml) and mod names m for 5mC and a 
            for 6mA.

            Furthermore, the file must have a .bam.bai index file with the same name.
            You can create an index if needed using samtools index.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        ref_genome: a string of Path objecting pointing to the .fasta file
            for the reference genome to which the .bam file is aligned.
        output_directory: optional str or Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        regions: TODO
        motifs: a list of strings specifying which base modifications to look for.
            The basemods are each specified as {sequence_motif},{position_of_modification}.
            For example, a methylated adenine is specified as 'A,0' and CpG methylation
            is specified as 'CG,0'.
        thresh: float point number specifying the base modification probability threshold
            used to delineate modificaton calls as True or False. When set to None, modkit
            will select its own threshold automatically based on the data.
        window_size: an integer specifying a window around the center of each bed_file
            region. If set to None, the bed_file is used unmodified. If set to a non-zero
            positive integer, the bed_file regions are replaced by new regions with that
            window size in either direction of the center of the original bed_file regions.
            This is used for e.g. extracting information from around known motifs or peaks.
        cores: an integer specifying how many parallel cores modkit gets to use. 
            By default modkit will use all of the available cores on the machine.
        log: a boolean specifying whether to output logs into the output folder.
        clean_intermediate: a boolean specifying whether to clean up to keep intermediate
            outputs. The final processed files are not human-readable, whereas the intermediate
            outputs are. However, intermediate outputs can also be quite large.

    Returns:
        Path object pointing to the compressed and indexed .bed.gz bedmethyl file, ready 
        for plotting functions.
        
    """

    input_file, ref_genome, output_directory = sanitize_path_args(
        input_file, ref_genome, output_directory
    )
    
    check_bam_format(input_file, motifs)
    
    # TODO: Add .tbi file? Add windowed bed file too, maybe?
    output_path, (output_bedmethyl, output_bedmethyl_sorted, output_bedgz_sorted) = prep_outputs(
        output_directory=output_directory,
        output_name=output_name,
        input_file=input_file,
        output_file_names=['pileup.bed', 'pileup.sorted.bed', 'pileup.sorted.bed.gz']
    )
    
    # TODO: This is mildly confusing. I get what it's doing, but it's hard to follow / names are bad. Also, why is it used in cleanup here, but not in extract?
    region_specifier, bed_filepath_processed = create_region_specifier(
        output_path,
        regions,
        window_size,
    )
    
    motif_command_list = []
    if len(motifs)>0:
        for basemod in motifs:
            # TODO: Can be split out to method; same functionality as in extract
            motif_details = basemod.split(',')
            motif_command_list.append('--motif')
            motif_command_list.append(motif_details[0])
            motif_command_list.append(motif_details[1])
    else:
        raise('Error: no motifs specified. Nothing to process.')
    
    if log:
        print('Logging to ',Path(output_path)/'pileup-log')
        log_command=['--log-filepath',Path(output_path)/'pileup-log']
    else:
        log_command=[]
    
    # TODO: This should be a method, like create_region_specifier, or just combined into a prep method for the start...
    cores_avail = multiprocessing.cpu_count()
    if cores is None:
        print(f'No specified number of cores requested. {cores_avail} available on machine, allocating {cores_avail//2}')
        cores_command_list = ['--threads',str(cores_avail//2)]
    elif cores>cores_avail:
        print(f'Warning: {cores} cores request, {cores_avail} available. Allocating {cores_avail}')
        cores_command_list = ['--threads',str(cores_avail)]
    else:
        print(f'Allocating requested {cores} cores.')
        cores_command_list = ['--threads',str(cores)]

    # TODO: This is SO SO SO similar to extract; just the ValueError vs. printing. I think this can be resolved
    mod_thresh_list = []
    if thresh is None:
        print('No base modification threshold provided. Using adaptive threshold selection via modkit.')
    elif thresh<=0:
        raise ValueError(f'Threshold {thresh} cannot be used for pileup, please pick a positive nonzero value.')
    else:
        adjusted_threshold = utils.adjust_threshold(thresh)
        for modnames_set in utils.BASEMOD_NAMES_DICT.values():
            for modname in modnames_set:
                mod_thresh_list = mod_thresh_list + ['--mod-thresholds',f'{modname}:{adjusted_threshold}']
    
    pileup_command_list = (['modkit', 'pileup', input_file, output_bedmethyl]
                           + region_specifier 
                           + motif_command_list 
                           + ['--ref',ref_genome,'--filter-threshold','0'] 
                           + mod_thresh_list 
                           + cores_command_list
                           + log_command)
    
    subprocess.run(pileup_command_list)
    
    with open(output_bedmethyl_sorted,'w') as sorted_file:
        subprocess.run(['sort', '-k1,1', '-k2,2n', output_bedmethyl], stdout=sorted_file)
    pysam.tabix_compress(output_bedmethyl_sorted,output_bedgz_sorted,force=True)
    pysam.tabix_index(str(output_bedgz_sorted),preset='bed',force=True)
    
    # TODO: Can cleanup be consolidated?
    if cleanup:
        if bed_filepath_processed is not None:
            bed_filepath_processed.unlink()
        os.remove(output_bedmethyl)
        os.remove(output_bedmethyl_sorted)
    
    return output_bedgz_sorted

def extract(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    regions: str | Path | list[str | Path] = None,
    motifs: list = ['A,0','CG,0','GCH,1'],
    thresh: float = None,
    window_size: int = 0,
    cores: int = None,
    log: bool = False,
    cleanup: bool = True,) -> Path:

    """
    TODO: Merge bed_file / region_str / window_size handling into a unified function somewhere

    Takes a file containing long read sequencing data aligned 
    to a reference genome with modification calls for one or more base/context 
    and pulls out data from each individual read. The intermediate outputs contain
    a plain-text list of all base modifications, split out by type. The compressed
    and indexed output contains vectors of valid and modified positions within each
    read.

    The current implementation of this method uses modkit, a tool built by 
    Nanopore Technologies, along with h5py to build the final output file.

    https://github.com/nanoporetech/modkit/

    Args:
        output_file: a string or Path object pointing to the location of a .bam file.
            The file should follow at least v1.6 of the .bam file specifications, 
            found here: https://samtools.github.io/hts-specs/
            https://samtools.github.io/hts-specs/SAMv1.pdf 

            The file needs to have modifications stored in the standard format,
            with MM and ML tags (NOT mm and ml) and mod names m for 5mC and a 
            for 6mA.

            Furthermore, the file must have a .bam.bai index file with the same name.
            You can create an index if needed using samtools index.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        ref_genome: a string of Path objecting pointing to the .fasta file
            for the reference genome to which the .bam file is aligned.
        output_directory: optional str or Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        regions: TODO
        motifs: a list of strings specifying which base modifications to look for.
            The basemods are each specified as {sequence_motif},{position_of_modification}.
            For example, a methylated adenine is specified as 'A,0' and CpG methylation
            is specified as 'CG,0'.
        thresh: float point number specifying the base modification probability threshold
            used to delineate modificaton calls as True or False. When set to None, modkit
            will select its own threshold automatically based on the data.
        window_size: an integer specifying a window around the center of each bed_file
            region. If set to None, the bed_file is used unmodified. If set to a non-zero
            positive integer, the bed_file regions are replaced by new regions with that
            window size in either direction of the center of the original bed_file regions.
            This is used for e.g. extracting information from around known motifs or peaks.
        cores: an integer specifying how many parallel cores modkit gets to use. 
            By default modkit will use all of the available cores on the machine.
        log: a boolean specifying whether to output logs into the output folder.
        clean_intermediate: a boolean specifying whether to clean up to keep intermediate
            outputs. The final processed files are not human-readable, whereas the intermediate
            outputs are. However, intermediate outputs can also be quite large.

    Returns:
        Path object pointing to the compressed and indexed output .h5 file, ready for
        plotting functions.

    """
    input_file, ref_genome, output_directory = sanitize_path_args(
        input_file, ref_genome, output_directory
    )
    
    check_bam_format(input_file, motifs)

    # TODO: Add intermediate mod-specific .txt files?
    output_path, (output_h5,) = prep_outputs(
        output_directory=output_directory,
        output_name=output_name,
        input_file=input_file,
        output_file_names=['reads.combined_basemods.h5']
    )
    
    region_specifier, bed_filepath_processed = create_region_specifier(
        output_path,
        regions,
        window_size,
    )
    
    cores_avail = multiprocessing.cpu_count()
    if cores is None:
        print(f'No specified number of cores requested. {cores_avail} available on machine, allocating {cores_avail//2}')
        cores_command_list = ['--threads',str(cores_avail//2)]
    elif cores>cores_avail:
        print(f'Warning: {cores} cores request, {cores_avail} available. Allocating {cores_avail}')
        cores_command_list = ['--threads',str(cores_avail)]
    else:
        print(f'Allocating requested {cores} cores.')
        cores_command_list = ['--threads',str(cores)]
        
    mod_thresh_list = []
    if thresh is None:
        print('No valid base modification threshold provided. Raw probs will be saved.')
        adjusted_threshold = None
    elif thresh<=0:
        print(f'With a thresh of {thresh}, modkit will simply save all tagged modifications.')
        adjusted_threshold = 0
    else:
        adjusted_threshold = utils.adjust_threshold(thresh)
        for modnames_set in utils.BASEMOD_NAMES_DICT.values():
            for modname in modnames_set:
                mod_thresh_list = mod_thresh_list + ['--mod-thresholds',f'{modname}:{adjusted_threshold}']

    if log:
        print('logging to ',Path(output_path)/'extract-log')
        log_command=['--log-filepath',Path(output_path)/f'extract-log']
    else:
        log_command=[]   

    for basemod in motifs:    
        print(f'Extracting {basemod} sites')
        motif_command_list = []
        motif_details = basemod.split(',')
        motif_command_list.append('--motif')
        motif_command_list.append(motif_details[0])
        motif_command_list.append(motif_details[1])

        output_txt = Path(output_path)/(f'reads.{basemod}.txt')
        
        if os.path.exists(output_txt):
            os.remove(output_txt)
        
        extract_command_list = (['modkit', 'extract', input_file, output_txt]
                                + region_specifier
                                + motif_command_list
                                + cores_command_list
                                + log_command
                                + ['--ref',ref_genome, '--filter-threshold','0',])
        subprocess.run(extract_command_list)
        
        print(f'Adding {basemod} to {output_h5}')
        read_by_base_txt_to_hdf5(
            output_txt,
            output_h5,
            basemod,
            adjusted_threshold,
        )
        if cleanup:
            os.remove(output_txt)
    if cleanup:
        if bed_filepath_processed is not None:
            bed_filepath_processed.unlink()
            
    return output_h5

"""
Helper functions to facilitate bam parse operations

check_bam_format: verify that a bam is formatted correctly to be processed.
create_region_specifier: create a list to append to the modkit call for specifying genomic regions.
adjust_threshold: backwards-compatible threshold adjustment, i.e. taking 0-255 thresholds and turning
    them into 0-1.
read_by_base_txt_to_hdf5: convert modkit extract txt into an .h5 file for rapid read access.
"""

def check_bam_format(
    bam_file: str | Path,
    basemods: list = ['A,0','CG,0'],
) -> bool:
    """
    Check whether a .bam file is formatted appropriately for modkit

    Args:
        bam_file: a formatted .bam file with a .bai index
        basemods: a list of base modification motifs
    
    Returns:
        None. If the function returns, you are ok. 
        
    """
    basemods_found_dict = {}
    for basemod in basemods:
        motif,pos = basemod.split(',')
        base = motif[int(pos)]
        basemods_found_dict[base] = False
    
    input_bam = pysam.AlignmentFile(bam_file)
    
    try:
        for counter,read in enumerate(input_bam.fetch()):
            read_dict = read.to_dict()
            for tag_string in read_dict['tags']:
                tag_fields = tag_string.split(',')[0].split(':')
                tag = tag_fields[0]
                # tag_type = tag_fields[1]
                tag_value = tag_fields[2]
                if tag=='Mm' or tag=='Ml':
                    raise ValueError(f'Base modification tags are out of spec (Mm and Ml instead of MM and ML). \n\nConsider using "modkit update-tags {str(bam_file)} new_file.bam" in the command line with your conda environment active and then trying with the new file. For megalodon basecalling/modcalling, you may also need to pass "--mode ambiguous.\nBe sure to index the resulting .bam file."')
                elif tag=='MM':
                    if len(tag_value)>0 and tag_value[-1]!='?' and tag_value[-1]!='.':
                        raise ValueError(f'Base modification tags are out of spec. Need ? or . in TAG:TYPE:VALUE for MM tag, else modified probability is considered to be implicit. \n\nConsider using "modkit update-tags {str(bam_file)} new_file.bam --mode ambiguous" in the command line with your conda environment active and then trying with the new file.')
                    else:
                        if len(tag_value)>0 and tag_value[0] in basemods_found_dict.keys():
                            if tag_value[2] in utils.BASEMOD_NAMES_DICT[tag_value[0]]:
                                basemods_found_dict[tag_value[0]] = True
                            else:
                                raise ValueError(f'Base modification name unexpected: {tag_value[2]} to modify {tag_value[0]}, should be in set {utils.BASEMOD_NAMES_DICT[tag_value[0]]}. \n\nIf you know what your mod names correspond to in terms of the latest .bam standard, consider using "modkit adjust-mods {str(bam_file)} new_file.bam --convert 5mC_name m --convert N6mA_name a --convert other_basemod_name correct_label" and then trying with the new file. Note: currently supported mod names are {utils.BASEMOD_NAMES_DICT}')
            if all(basemods_found_dict.values()):
                return
            if counter>=NUM_READS_TO_CHECK:
                missing_bases = []
                for base,found in basemods_found_dict.items():
                    if not found:
                        missing_bases.append(base)
                print(f'WARNING: no modified values found for {missing_bases} in the first {counter} reads. Do you expect this file to contain these modifications? parse_bam is looking for {basemods} but only found modifications on {[base for base, found in basemods_found_dict.items() if found]}. \n\nConsider passing only the basemods that you expect to be present in your file.')
                return
    except ValueError as e:
        if 'fetch called on bamfile without index' in str(e):
            raise ValueError(f'{e}. Consider using "samtools index {str(bam_file)}" to create an index if your .bam is already sorted.')
        else:
            raise
    except:
        raise
        
def get_alignment_quality(
    bam_file,
    ref_genome,
) -> float:
    input_bam = pysam.AlignmentFile(bam_file,'rb')
    for read in enumerate(input_bam.fetch()):
        read_sequence = read.get_forward_sequence()
        read_alignment = read.get_forward_positions()
        print(read_alignment)
        
def create_region_specifier(
    output_path,
    regions,
    window_size,
):
    """
    Creates commands to pass to modkit based on bed_file regions.
    """
    
    if regions is not None:
        bed_filepath_processed = output_path / 'regions.processed.bed'
        regions_dict = utils.regions_dict_from_input(
            regions,
            window_size,
        )
        region_specifier = ['--include-bed',str(bed_filepath_processed)]
    else:
        bed_filepath_processed = None
        region_specifier = []

    utils.bed_from_regions_dict(regions_dict,bed_filepath_processed)
    

    return region_specifier, bed_filepath_processed

def read_by_base_txt_to_hdf5(
    input_txt: str | Path,
    output_h5: str | Path,
    basemod: str,
    thresh: float=None,
) -> None:
    """
    Takes in a txt file generated by modkit extract and appends
    all the data from a specified basemod into an hdf5 file. If a thresh is specified, it
    also binarizes the mod calls.

    Args:
        input_txt: a string or Path pointing to a modkit extracted base-by-base modifications
            file. This file is assumed to have been created by modkit v0.2.4, other versions may
            have a different format and may not function normally.
        output_h5: a string or Path pointing to a valid place to save an .h5 file. If this
            file already exists, it will not be cleared and will simply be appended to. If it does
            not exist it will be created and datasets will be added for read_name, chromosome, read_start,
            read_end, base modification motif, mod_vector, and val_vector.
        basemod: a string specifying a single base modification. Basemods are specified as 
            {sequence_motif},{position_of_modification}. For example, a methylated adenine is specified 
            as 'A,0' and CpG methylation is specified as 'CG,0'.
        thresh: a floating point threshold for base modification calling, between zero and one. 
            If specified as None, raw probabilities will be saved in the .h5 output.

    Returns:
        None

    """
    motif,modco = tuple(basemod.split(','))
    motif_modified_base = motif[int(modco)]
    read_name = ''
    num_reads = 0
    # TODO: I think the function calls can be consolidated; lots of repetition
    with open(input_txt) as txt:
        for index,line in enumerate(txt):
            fields = line.split('\t')
            if index>0 and read_name!=fields[0]:
                read_name = fields[0]
                num_reads+=1
        print(f'{num_reads} reads found in {input_txt}')
        txt.seek(0)
        with h5py.File(output_h5,'a') as h5:
            # Create datasets
            dt_str = h5py.string_dtype(encoding='utf-8')
            if thresh==None:
                dt_vlen = h5py.vlen_dtype(np.float16)
            else:
                dt_vlen = h5py.vlen_dtype(bool)
            if 'read_name' in h5:

                old_size = h5['read_name'].shape[0]
                print(f'extending from {old_size} to {old_size+num_reads}')
                h5['read_name'].resize((old_size+num_reads,))
            else:
                old_size = 0
                h5.create_dataset(
                    'read_name',
                    (num_reads,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression='gzip',
                    compression_opts=9,
                )
            if 'chromosome' in h5:
                if old_size != h5['chromosome'].shape[0]:
                    print('size mismatch: read_name:chromosome')
                else:
                    h5['chromosome'].resize((old_size+num_reads,))
            else:
                h5.create_dataset(
                    'chromosome',
                    (num_reads,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression='gzip',
                    compression_opts=9,
                )
            if 'read_start' in h5:
                if old_size != h5['read_start'].shape[0]:
                    print('size mismatch','read_name','read_start')
                else:
                    h5['read_start'].resize((old_size+num_reads,))
            else:
                h5.create_dataset(
                    'read_start',
                    (num_reads,),
                    maxshape=(None,),
                    dtype='i',
                    compression='gzip',
                    compression_opts=9,
                )
            if 'read_end' in h5:
                if old_size != h5['read_end'].shape[0]:
                    print('size mismatch','read_name','read_end')
                else:
                    h5['read_end'].resize((old_size+num_reads,))
            else:
                h5.create_dataset(
                    'read_end',
                    (num_reads,),
                    maxshape=(None,),
                    dtype='i',
                    compression='gzip',
                    compression_opts=9,
                )
            if 'strand' in h5:
                if old_size != h5['strand'].shape[0]:
                    print('size mismatch','read_name','strand')
                else:
                    h5['strand'].resize((old_size+num_reads,))
            else:
                h5.create_dataset(
                    'strand',
                    (num_reads,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression='gzip',
                    compression_opts=9,
                )
            if 'motif' in h5:
                if old_size != h5['motif'].shape[0]:
                    print('size mismatch','read_name','motif')
                else:
                    h5['motif'].resize((old_size+num_reads,))
            else:
                h5.create_dataset(
                    'motif',
                    (num_reads,),
                    maxshape=(None,),
                    dtype=dt_str,
                    compression='gzip',
                    compression_opts=9,
                )
            if 'mod_vector' in h5:
                if old_size != h5['mod_vector'].shape[0]:
                    print('size mismatch read_name:mod_vector')
                else:
                    h5['mod_vector'].resize((old_size+num_reads,))
            else:
                h5.create_dataset(
                    'mod_vector',
                    (num_reads,),
                    maxshape=(None,),
                    dtype=dt_vlen,
                    compression='gzip',
                    compression_opts=9,
                )
            if 'val_vector' in h5:
                if old_size != h5['val_vector'].shape[0]:
                    print('size mismatch read_name:val_vector')
                else:
                    h5['val_vector'].resize((old_size+num_reads,))
            else:
                h5.create_dataset(
                    'val_vector',
                    (num_reads,),
                    maxshape=(None,),
                    dtype=dt_vlen,
                    compression='gzip',
                    compression_opts=9,
                )


    #         next(txt)
            read_name = ''
            read_counter = 0
            readlen_sum = 0
            for index,line in enumerate(txt):
                if index==0:
#                     print(line)
                    continue
                fields = line.split('\t')
                pos_in_genome = int(fields[2])
                canonical_base = fields[15]
                prob = float(fields[10])
                
                if read_name!=fields[0]:
                    # Record the read details unless this is the first read
                    if index>1:
                        if len(val_coordinates_list)>0:
                            read_len_along_ref = max(val_coordinates_list)+1
                        else:
                            read_len_along_ref = read_len
                        mod_vector = np.zeros(read_len_along_ref,dtype=float)
                        mod_vector[val_coordinates_list] = mod_values_list
                        val_vector = np.zeros(read_len_along_ref,dtype=int)
                        val_vector[val_coordinates_list] = 1
                        # Build mod_vector and val_vector from lists
                        h5['read_name'][read_counter+old_size]=read_name
                        h5['chromosome'][read_counter+old_size]=read_chrom
                        h5['read_start'][read_counter+old_size]=read_start
                        h5['read_end'][read_counter+old_size]=read_end
                        h5['strand'][read_counter+old_size]=ref_strand
                        h5['motif'][read_counter+old_size]=basemod
                        h5['mod_vector'][read_counter+old_size]=mod_vector
                        h5['val_vector'][read_counter+old_size]=val_vector
                        read_counter+=1
                    # Set the read name of the next read
                    read_name=fields[0]
                    # Store some relevant read metadata
                    read_chrom = fields[3]
                    read_len = int(fields[9])                    
                    ref_strand = fields[5]
                    if ref_strand == '+':
                        pos_in_read_ref = int(fields[1])
                    elif ref_strand == '-':
                        pos_in_read_ref = read_len - int(fields[1]) - 1
                    #Calculate read info
                    read_start = pos_in_genome - pos_in_read_ref
                    read_end = read_start + read_len
                    # Instantiate lists
                    mod_values_list = []
                    val_coordinates_list = []

                # Regardless of whether its a new read or not, 
                # add modification to vector if motif type is correct
                # for the motif in question
                if canonical_base == motif_modified_base:
                    val_coordinates_list.append(pos_in_genome-read_start)
                    if thresh==None:
                        mod_values_list.append(prob)
                    elif prob>=thresh:
                        mod_values_list.append(1)
                    else:
                        mod_values_list.append(0)


            # Save the last read
            if len(read_name)>0:
                # Build the vectors
                if len(val_coordinates_list)>0:
                    read_len_along_ref = max(val_coordinates_list)+1
                else:
                    read_len_along_ref = read_len
                mod_vector = np.zeros(read_len_along_ref,dtype=float)
                mod_vector[val_coordinates_list] = mod_values_list
                val_vector = np.zeros(read_len_along_ref,dtype=int)
                val_vector[val_coordinates_list] = 1
                h5['read_name'][read_counter+old_size]=read_name
                h5['chromosome'][read_counter+old_size]=read_chrom
                h5['read_start'][read_counter+old_size]=read_start
                h5['read_end'][read_counter+old_size]=read_end
                h5['strand'][read_counter+old_size]=ref_strand
                h5['motif'][read_counter+old_size]=basemod
                h5['mod_vector'][read_counter+old_size]=mod_vector
                h5['val_vector'][read_counter+old_size]=val_vector
                read_counter+=1
    return

def sanitize_path_args(*args) -> tuple:
    """
    Coerce all given arguments to Path objects, leaving Nones as Nones.
    """
    return tuple(Path(f) if f is not None else f for f in args)


def prep_outputs(output_directory: Path | None,
                 output_name: str,
                 input_file: Path,
                 output_file_names: list[str]) -> tuple[Path, list[Path]]:
    """
    As a side effect, if files exist that match the requested outputs, they are deleted.

    TODO: Is it kind of silly that this takes in input_file? Maybe should take in some generic default parameter, or this default should be set outside this method?
    Args:
        output_directory: Path pointing to an output directory.
            If left as None, outputs will be stored in a new folder within the input
            directory.
        output_name: a string that will be used to create an output folder
            containing the intermediate and final outputs, along with any logs.
        input_file: Path to input file; used to define default output directory
        output_file_names: list of names of desired output files

    Returns:
        * Path to top-level output directory
        * List of Paths to requested output files
    """
    if output_directory is None:
        output_directory = input_file.parent
        print(f'No output directory provided, using input directory {output_directory}')
          
    output_path = output_directory / output_name

    output_files = [output_path / file_name for file_name in output_file_names]

    # Ensure output path exists, and that any of the specified output files do not already exist (necessary for some outputs)
    output_path.mkdir(parents=True, exist_ok=True)
    for output_file in output_files:
        output_file.unlink(missing_ok=True)

    return output_path, output_files
