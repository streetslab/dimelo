import sys
import os
import shutil
import subprocess
import multiprocessing
from pathlib import Path

from .config import EXE_CONFIG

# # Get the directory one level up
# parent_dir = os.path.dirname(os.getcwd())
# # Add the parent directory to sys.path
# sys.path.append(parent_dir)
from . import utils

"""
TODO: The name "parser" conflicts with a built-in python library, I think. VSCode definitely complains at me.
"""

def parse_bam_pileup(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    region_str=None,
    bed_file=None,
    basemods = ['A,0','CG,0','GCH,1'],
    thresh = 0,
    window_size=0,
    cores=None,
    log=False,
):
    """
    TODO: Documentation
    TODO: Raise errors, maybe return bools, but don't return ints
    TODO: Where should the windowed file be written to? Right now, seems to write to the location of the reference file; this could be confusing.
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
    else:
        raise('Error: cannot process both a region and a bed file.')
    
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
        print(f'Modification threshold of {thresh} assumed to be for range 0-255. {thresh}/255={thesh/255} will be sent to modkit.')
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
    with open(output_bedgz_sorted,'w') as compressed_file:
        subprocess.run([EXE_CONFIG.bgzip_exe,
                        '-c',output_bed_sorted],
                        stdout=compressed_file)
    subprocess.run([EXE_CONFIG.tabix_exe,
                    '-p','bed',output_bedgz_sorted])

    return output_bedgz_sorted

def parse_bam_extract(
    input_file: str | Path,
    output_name: str,
    ref_genome: str | Path,
    output_directory: str | Path = None,
    region_str=None,
    bed_file=None,
    basemods = ['A,0','CG,0','GCH,1'],
    thresh = 0,
    window_size=-1,
    cores=None,
    log=False
):
    """
    TODO: Documentation
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
    else:
        raise('Error: cannot process both a region and a bed file.')
    
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
        print(f'Modification threshold of {thresh} assumed to be for range 0-255. {thresh}/255={thesh/255} will be sent to modkit.')
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
        utils.read_by_base_txt_to_hdf5(
            output_txt,
            output_h5,
            basemod,
            thresh,
        )
    return output_h5
