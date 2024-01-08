import subprocess
import sys
import os
import multiprocessing
# Get the directory one level up
parent_dir = os.path.dirname(os.getcwd())

# Add the parent directory to sys.path
sys.path.append(parent_dir)
import utils
from pathlib import Path
def parse_bam_modkit_pileup(
    input_file: str | Path,
    output_path: str | Path,
    output_name: str,
    ref_genome: str | Path,
    region_str=None,
    bed_file=None,
    basemods = ['A,0','CG,0','GCH,1'],
    thresh = 0,
    window_size=-1,
    cores=None
):
    if bed_file is not None and region_str is None:
        if window_size>0:
            bed_filepath = Path(bed_file)
            print(f'Processing from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
            bed_filepath_processed = bed_filepath.parent / (bed_filepath.stem + '.windowed' + bed_filepath.suffix)
            print(f'Writing new bed file {bed_filepath_processed.name}')
            utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
            region_specifier = ['--include-bed',bed_filepath_processed]
        elif window_size==-1:
            bed_filepath_processed = Path(bed_file)
            print(f'Processing from {bed_filepath_processed.name} using unmodified bed regions.')
            region_specifier = ['--include-bed',bed_filepath_processed]
        else:
            print(f'Error: invalid window size {window_size}bp')
            return -1
    elif bed_file is None and region_str is not None:
        if window_size==-1:
            print(f'Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
        else:
            print(f'Warning: window size {window_size}bp will be ignored. Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
    else:
        print('Error: cannot process both a region and a bed file.')
        return -1
    
    motif_command_list = []
    if len(basemods)>0:
        for basemod in basemods:
            motif_details = basemod.split(',')
            motif_command_list.append('--motif')
            motif_command_list.append(motif_details[0])
            motif_command_list.append(motif_details[1])
    else:
        print('Error: no basemods specified. Nothing to process.')
        return -1
    
    cores_avail = multiprocessing.cpu_count()
    if cores>cores_avail:
        print(f'Warning: {cores} cores request, {cores_avail} available. Allocating {cores_avail}')
        cores_command_list = ['--threads',str(cores_avail)]
    elif cores is None:
        print(f'No specified number of cores requested. {cores_avail} available on machine, allocating {cores_avail//2}')
        cores_command_list = ['--threads',str(cores_avail//2)]
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
        
    output_bed = Path(output_path)/(output_name+'.bed')
    output_bed_sorted = Path(output_path)/(output_name+'.sorted.bed')
    output_bedgz_sorted = Path(output_path)/(output_name+'.sorted.bed.gz')
    
    pileup_command_list = (['../dependencies/modkit/modkit','pileup',
                             input_file,
                             output_bed]
                           + region_specifier 
                           + motif_command_list 
                           + ['--ref',ref_genome,'--filter-threshold','0'] 
                           + mod_thresh_list 
                           + cores_command_list)
    
    subprocess.run(pileup_command_list)
    with open(output_bed_sorted,'w') as sorted_file:
        subprocess.run(['sort', '-k1,1', '-k2,2n', output_bed], stdout=sorted_file)
    with open(output_bedgz_sorted,'w') as compressed_file:
        subprocess.run(['../dependencies/htslib/bin/bgzip',
                    '-c',output_bed_sorted],
                       stdout=compressed_file)
    subprocess.run(['../dependencies/htslib/bin/tabix',
                    '-p','bed',output_bedgz_sorted])
                             
#     result = subprocess.run(['ls', '-l'], capture_output=True, text=True)
    return 0
def parse_bam_modkit_extract(
    input_file: str | Path,
    output_path: str | Path,
    output_name: str,
    ref_genome: str | Path,
    region_str=None,
    bed_file=None,
    basemods = ['A,0','CG,0','GCH,1'],
    thresh = 0,
    window_size=-1,
    cores=None
):
    if bed_file is not None and region_str is None:
        if window_size>0:
            bed_filepath = Path(bed_file)
            print(f'Processing from {bed_filepath.name} using even {window_size}bp windows in either direction from bed region centers.')
            bed_filepath_processed = bed_filepath.parent / (bed_filepath.stem + '.windowed' + bed_filepath.suffix)
            print(f'Writing new bed file {bed_filepath_processed.name}')
            utils.generate_centered_windows_bed(bed_filepath,bed_filepath_processed,window_size)
            region_specifier = ['--include-bed',bed_filepath_processed]
        elif window_size==-1:
            bed_filepath_processed = Path(bed_file)
            print(f'Processing from {bed_filepath_processed.name} using unmodified bed regions.')
            region_specifier = ['--include-bed',bed_filepath_processed]
        else:
            print(f'Error: invalid window size {window_size}bp')
            return -1
    elif bed_file is None and region_str is not None:
        if window_size==-1:
            print(f'Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
        else:
            print(f'Warning: window size {window_size}bp will be ignored. Processing from region {region_str}.')
            region_specifier = ['--region',region_str]
    else:
        print('Error: cannot process both a region and a bed file.')
        return -1
    
    motif_command_list = []
    if len(basemods)>0:
        for basemod in basemods:
            motif_details = basemod.split(',')
            motif_command_list.append('--motif')
            motif_command_list.append(motif_details[0])
            motif_command_list.append(motif_details[1])
    else:
        print('Error: no basemods specified. Nothing to process.')
        return -1
    
    cores_avail = multiprocessing.cpu_count()
    if cores>cores_avail:
        print(f'Warning: {cores} cores request, {cores_avail} available. Allocating {cores_avail}')
        cores_command_list = ['--threads',str(cores_avail)]
    elif cores is None:
        print(f'No specified number of cores requested. {cores_avail} available on machine, allocating {cores_avail//2}')
        cores_command_list = ['--threads',str(cores_avail//2)]
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
    
    output_txt = Path(output_path)/(output_name+'.txt')
    subprocess.run(['../dependencies/modkit/modkit','extract',
                             input_file,
                             output_txt,
                             '--region',region_str,
                             '--motif' ,'CG','0',
                             '--ref',genome_str,
                             '--filter-threshold','0',
                             '--mod-thresholds', f'm:{thresh}',
                             '--mod-thresholds', f'a:{thresh}'])
    return 0