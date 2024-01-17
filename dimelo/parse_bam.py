import sys
import os
import shutil
import subprocess
import multiprocessing
from pathlib import Path

import numpy as np
import h5py
import pysam

from .config import EXE_CONFIG

from . import utils

"""
This module contains code to convert .bam files into both human-readable and 
indexed random-access pileup and read-wise processed outputs.
"""

def parse_bam_pileup(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    region_str: str = None,
    bed_file: str | Path = None,
    basemods: list = ['A,0','CG,0','GCH,1'],
    thresh: float = None,
    window_size: int = 0,
    cores: int = None,
    log: bool = False,
    clean_intermediate: bool = True,) -> Path:

    """
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
        region_str: optional string in format chr{num}:{start}-{end} specifying a single
            region to process. Mutually exclusive with bed_file.
        bed_file: optional string or Path for a .bed file containing one or more regions
            of the genome to process. Mutually exclusive with region_str.
        basemods: a list of strings specifying which base modifications to look for.
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
    
    if output_directory is None:
        output_directory = Path(input_file).parent
        print(f'No output directory provided, using input directory {output_directory}')
    
    output_path = Path(output_directory)/output_name
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path,exist_ok=True)
    
    if bed_file is not None and region_str is None:
        if window_size>0:
            bed_filepath = Path(bed_file)
            print(f'Processing from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
            bed_filepath_processed = output_path / (bed_filepath.stem + f'.windowed{window_size}-for-pileup' + bed_filepath.suffix)
            print(f'Writing new bed file {bed_filepath_processed.name}')
            utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
            region_specifier = ['--include-bed',bed_filepath_processed]
        elif window_size==0:
            bed_filepath_processed = Path(bed_file)
            print(f'Processing from {bed_filepath_processed.name} using unmodified bed regions.')
            region_specifier = ['--include-bed',bed_filepath_processed]
        else:
            raise(f'Error: invalid window size {window_size}bp')
    elif bed_file is None and region_str is not None:
        if window_size==0:
            print(f'Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
        else:
            print(f'Warning: window size {window_size}bp will be ignored. Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
    elif bed_file is None and region_str is None:
        print('No region(s) specified, processing the entire genome.')
        region_specifier = []
        if window_size is not None:
            print('A window_size was specified but will be ignored.')
    else:
        raise ValueError('Error: cannot process both a region and a bed file.')
    
    motif_command_list = []
    if len(basemods)>0:
        for basemod in basemods:
            motif_details = basemod.split(',')
            motif_command_list.append('--motif')
            motif_command_list.append(motif_details[0])
            motif_command_list.append(motif_details[1])
    else:
        raise('Error: no basemods specified. Nothing to process.')
    
    if log:
        print('Logging to ',Path(output_path)/'pileup-log')
        log_command=['--log-filepath',Path(output_path)/'pileup-log']
    else:
        log_command=[]
    
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
        
    
    if thresh<=0:
        print('No valid base modification threshold provided. Using adaptive threshold selection via modkit.')
        mod_thresh_list = []
    elif thresh>1:
        print(f'Modification threshold of {thresh} assumed to be for range 0-255. {thresh}/255={thresh/255} will be sent to modkit.')
        mod_thresh_list = ['--mod-thresholds', f'm:{thresh/255}','--mod-thresholds', f'a:{thresh/255}']
    else:
        print(f'Modification threshold of {thresh} will be treated as coming from range 0-1.')
        mod_thresh_list = ['--mod-thresholds', f'm:{thresh}','--mod-thresholds', f'a:{thresh}']
        
    output_bed = Path(output_path)/('pileup.bed')
    output_bed_sorted = Path(output_path)/('pileup.sorted.bed')
    output_bedgz_sorted = Path(output_path)/('pileup.sorted.bed.gz')
    
    pileup_command_list = ([EXE_CONFIG.modkit_exe,
                            'pileup',
                            input_file,
                            output_bed]
                           + region_specifier 
                           + motif_command_list 
                           + ['--ref',ref_genome,'--filter-threshold','0'] 
                           + mod_thresh_list 
                           + cores_command_list
                           + log_command)
    
    subprocess.run(pileup_command_list)
    with open(output_bed_sorted,'w') as sorted_file:
        subprocess.run(['sort', '-k1,1', '-k2,2n', output_bed], stdout=sorted_file)
    # with open(output_bedgz_sorted,'w') as compressed_file:
    #     subprocess.run([EXE_CONFIG.bgzip_exe,
    #                     '-c',output_bed_sorted],
    #                     stdout=compressed_file)
    # subprocess.run([EXE_CONFIG.tabix_exe,
    #                 '-p','bed',output_bedgz_sorted])
    pysam.tabix_compress(output_bed_sorted,output_bedgz_sorted,force=True)
    pysam.tabix_index(str(output_bedgz_sorted),preset='bed',force=True)

    return output_bedgz_sorted


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
                if read_name!=fields[0]:
                    try:
                        if len(read_name)>0:
                            h5['read_name'][read_counter+old_size]=read_name
                            h5['chromosome'][read_counter+old_size]=read_chrom
                            h5['read_start'][read_counter+old_size]=read_start
                            h5['read_end'][read_counter+old_size]=read_end
                            h5['motif'][read_counter+old_size]=basemod
                            h5['mod_vector'][read_counter+old_size]=mod_vector
                            h5['val_vector'][read_counter+old_size]=val_vector
                            read_counter+=1
                    except:
                        pass
                    #New read
                    #Record name
                    read_name = fields[0]
                    #Read in relevant values
                    pos_in_genome = int(fields[2])
                    read_chrom = fields[3]
                    read_len = int(fields[9])
                    readlen_sum+=read_len
                    canonical_base = fields[15]
                    prob = float(fields[10])
                    ref_strand = fields[5]
                    if ref_strand == '+':
                        pos_in_read_ref = int(fields[1])
                    elif ref_strand == '-':
                        pos_in_read_ref = read_len - int(fields[1]) - 1
                    #Calculate read info
                    read_start = pos_in_genome - pos_in_read_ref
                    read_end = read_start + read_len
                    #Build read vectors
                    mod_vector = np.zeros(read_len,dtype=np.float16)
                    val_vector = np.zeros(read_len,dtype=np.float16)

                    #Add modification to vector if type is correct
                    if canonical_base == motif_modified_base:
                        val_vector[pos_in_genome-read_start] = 1
                        if thresh==None:
                            mod_vector[pos_in_genome-read_start] = prob
                        elif prob>=thresh:
                            mod_vector[pos_in_genome-read_start] = 1
                else:
                    pos_in_genome = int(fields[2])
                    canonical_base = fields[15]
                    prob = float(fields[10])
                    if canonical_base == motif_modified_base:
                        val_vector[pos_in_genome-read_start] = 1
                        if thresh==None:
                            mod_vector[pos_in_genome-read_start] = prob
                        elif prob>=thresh:
                            mod_vector[pos_in_genome-read_start] = 1
            try:
                if len(read_name)>0:
                        h5['read_name'][read_counter+old_size]=read_name
                        h5['chromosome'][read_counter+old_size]=read_chrom
                        h5['read_start'][read_counter+old_size]=read_start
                        h5['read_end'][read_counter+old_size]=read_end
                        h5['motif'][read_counter+old_size]=basemod
                        h5['mod_vector'][read_counter+old_size]=mod_vector
                        h5['val_vector'][read_counter+old_size]=val_vector
                        read_counter+=1
            except:
                pass

#             print(readlen_sum/read_counter)
    return

def parse_bam_extract(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    region_str: str = None,
    bed_file: str | Path = None,
    basemods: list = ['A,0','CG,0','GCH,1'],
    thresh: float = 0,
    window_size: int = 0,
    cores: int = None,
    log: bool = False,
    clean_intermediate: bool = True,) -> Path:

    """
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
        region_str: optional string in format chr{num}:{start}-{end} specifying a single
            region to process. Mutually exclusive with bed_file.
        bed_file: optional string or Path for a .bed file containing one or more regions
            of the genome to process. Mutually exclusive with region_str.
        basemods: a list of strings specifying which base modifications to look for.
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
    
    if output_directory is None:
        output_directory = Path(input_file).parent
        print(f'No output directory provided, using input directory {output_directory}')  
        
    output_path = Path(output_directory)/output_name
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    os.makedirs(output_path,exist_ok=True)
    
    if bed_file is not None and region_str is None:
        if window_size>0:
            bed_filepath = Path(bed_file)
            print(f'Processing from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
            bed_filepath_processed = output_path / (bed_filepath.stem + f'.windowed{window_size}-for-extract' + bed_filepath.suffix)
            print(f'Writing new bed file {bed_filepath_processed.name}')
            utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
            region_specifier = ['--include-bed',bed_filepath_processed]
        elif window_size==0:
            bed_filepath_processed = Path(bed_file)
            print(f'Processing from {bed_filepath_processed.name} using unmodified bed regions.')
            region_specifier = ['--include-bed',bed_filepath_processed]
        else:
            raise(f'Error: invalid window size {window_size}bp')
    elif bed_file is None and region_str is not None:
        if window_size==0:
            print(f'Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
        else:
            print(f'Warning: window size {window_size}bp will be ignored. Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
    elif bed_file is None and region_str is None:
        print('No region(s) specified, processing the entire genome.')
        region_specifier = []
        if window_size is not None:
            print('A window_size was specified but will be ignored.')    
    else:
        raise ValueError('Error: cannot process both a region and a bed file.')
    
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
        
    
    if thresh<=0:
        print('No valid base modification threshold provided. Raw probs will be saved.')
        mod_thresh_list = []
    elif thresh>1:
        print(f'Modification threshold of {thresh} assumed to be for range 0-255. {thresh}/255={thresh/255} will be sent to modkit.')
        mod_thresh_list = ['--mod-thresholds', f'm:{thresh/255}','--mod-thresholds', f'a:{thresh/255}']
    else:
        print(f'Modification threshold of {thresh} will be treated as coming from range 0-1.')
        mod_thresh_list = ['--mod-thresholds', f'm:{thresh}','--mod-thresholds', f'a:{thresh}']
    
    output_h5 = Path(output_path)/(f'reads.combined_basemods.h5')
    with open(output_h5,'w') as f:
        pass
    for basemod in basemods:    
        print(f'Extracting {basemod} sites')
        motif_command_list = []
        motif_details = basemod.split(',')
        motif_command_list.append('--motif')
        motif_command_list.append(motif_details[0])
        motif_command_list.append(motif_details[1])

        if log:
            print('logging to ',Path(output_path)/'extract-log')
            log_command=['--log-filepath',Path(output_path)/f'extract-log']
        else:
            log_command=[]



        output_txt = Path(output_path)/(f'reads.{basemod}.txt')
        subprocess.run(
          [EXE_CONFIG.modkit_exe,
          'extract',
          input_file,
          output_txt]
          +region_specifier
          +motif_command_list
          +log_command
          +['--ref',ref_genome,
          '--filter-threshold','0',]
        )
        print(f'Adding {basemod} to {output_h5}')
        read_by_base_txt_to_hdf5(
            output_txt,
            output_h5,
            basemod,
            thresh,
        )
    return output_h5
