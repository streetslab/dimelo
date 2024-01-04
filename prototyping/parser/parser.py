import subprocess
from pathlib import Path
def parse_bam_modkit_pileup(
    input_file,
    output_path,
    output_name,
    region_str,
    genome_str,
    thresh,
):
    output_bed = Path(output_path)/(output_name+'.bed')
    output_bed_sorted = Path(output_path)/(output_name+'.sorted.bed')
    output_bedgz_sorted = Path(output_path)/(output_name+'.sorted.bed.gz')
    subprocess.run(['../dependencies/modkit/modkit','pileup',
                             input_file,
                             output_bed,
                             '--region',region_str,
                             '--motif' ,'CG','0',
                             '--motif' ,'A','0',
                             '--motif' ,'GCH','1',
                             '--ref',genome_str,
                             '--filter-threshold','0',
                             '--mod-thresholds', f'm:{thresh}',
                             '--mod-thresholds', f'a:{thresh}'])
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
    input_file,
    output_path,
    output_name,
    region_str,
    genome_str,
    thresh,
):
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